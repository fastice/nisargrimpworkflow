#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ian
"""

import argparse
import glob
import json
import os
import re

import yaml
import utilities as u
from osgeo import ogr, osr
from nisargrimpworkflow.FileNISARProducts import parseFileName

defaultOutputDirName = 'frameInventory'

directionNames = {'A': 'ascending', 'D': 'descending'}

baselineNamePattern = re.compile(r'^baseline\.\d+x\d+\.yaml$')
# Accepts both rBaseline.deltab.yaml (old) and rBaseline.deltabp.yaml (renamed
# for consistency with the legacy rBaseline.deltabp) -- not rBaseline.deltabquad.yaml.
rBaselineNamePattern = re.compile(r'^rBaseline\.deltabp?\.yaml$')
sigmaLinePattern = re.compile(r'sigma\*sqrt\(X2/n\)\s*=\s*([\d.eE+-]+)')
tiepointsLinePattern = re.compile(
    r'Ntiepoints/Ngiven\s+used\s*=\s*(\d+)\s*/\s*(\d+)')

fieldDefs = [
    ('track', ogr.OFTInteger),
    ('orbit', ogr.OFTInteger),
    ('virtualFrameId', ogr.OFTString),
    ('orbit1', ogr.OFTInteger),
    ('orbit2', ogr.OFTInteger),
    ('cycle', ogr.OFTInteger),
    ('direction', ogr.OFTString),
    ('frames', ogr.OFTString),
    ('sigmaBaseline', ogr.OFTReal),
    ('sigmaRBaseline', ogr.OFTReal),
    ('sigmaAz', ogr.OFTReal),
    ('sigmaRBaselineWithoutIon', ogr.OFTReal),
    ('usingIonRBaseline', ogr.OFTInteger),
    ('nTiepointsBaseline', ogr.OFTInteger),
    ('nTiepointsGivenBaseline', ogr.OFTInteger),
    ('nTiepointsUsedRBaseline', ogr.OFTInteger),
    ('nTiepointsGivenRBaseline', ogr.OFTInteger),
    ('nTiepointsUsedAz', ogr.OFTInteger),
    ('nTiepointsGivenAz', ogr.OFTInteger),
]


def loadFramePattern(projectDir):
    '''Read framePattern (glob matching virtual-frame directory suffixes,
    e.g. '00??') from projectDir/project.yaml -- same key/default convention
    as setupNISARTracks.py/makeframetie.py. Defaults to '00??' if project.yaml
    or the key is absent.'''
    projPath = os.path.join(projectDir, 'project.yaml')
    proj = {}
    if os.path.exists(projPath):
        with open(projPath) as fp:
            proj = yaml.safe_load(fp) or {}
    return proj.get('framePattern', '00??')


def buildFrameGpkgArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mBuild per-cycle GeoPackages (ascending + '
        'descending layers) of NISAR virtual frames found under '
        "track-N/*_<framePattern> directories (framePattern from "
        "project.yaml, default '00??' -- canonical *_0000 groups and "
        'straggler-frame/fragment groups alike) \033[0m\n\n'
        'Run buildFrameLayers afterward to turn this output into a QGIS '
        '.qlr layer tree.\n\n')
    parser.add_argument('--projectDir', type=str, default='.',
                        help='Root directory containing track-N '
                        'subdirectories [.] (current directory)')
    parser.add_argument('-o', '--outputDir', type=str, default=None,
                        help='Output directory for per-cycle GeoPackages '
                        f'[<projectDir>/{defaultOutputDirName}]')
    args = parser.parse_args()
    outputDir = args.outputDir or \
        os.path.join(args.projectDir, defaultOutputDirName)
    return args.projectDir, outputDir


def findRequired(virtualFrameDir, relPath, pattern=None):
    '''Find a single required file in virtualFrameDir, optionally by glob
    pattern. Returns (path, None) on success, or (None, reason) if missing
    or ambiguous (reason is also printed immediately).'''
    if pattern is None:
        fullPath = os.path.join(virtualFrameDir, relPath)
        if not os.path.exists(fullPath):
            reason = f'missing {relPath}'
            print(f'{virtualFrameDir}: {reason}')
            return None, reason
        return fullPath, None
    matches = sorted(glob.glob(os.path.join(virtualFrameDir, relPath)))
    matches = [m for m in matches if pattern.match(os.path.basename(m))]
    if len(matches) != 1:
        reason = (f'expected exactly 1 match for {relPath} filtered by '
                  f'pattern, found {len(matches)}')
        print(f'{virtualFrameDir}: {reason}')
        return None, reason
    return matches[0], None


def readFrames(framesFile):
    with open(framesFile) as fp:
        frameNumbers = fp.readline().split()
    return ','.join(frameNumbers)


def readGeodat(geodatFile):
    with open(geodatFile) as fp:
        geojson = json.load(fp)
    imageName = geojson['properties']['ImageName']
    return geojson['geometry'], imageName


def readPairInfo(pairInfoFile):
    with open(pairInfoFile) as fp:
        orbit1, orbit2 = fp.readline().split()[:2]
    return int(orbit1), int(orbit2)


def readYamlSigma(yamlFile, usedKey):
    with open(yamlFile) as fp:
        data = yaml.safe_load(fp)
    if not data:
        raise ValueError(f'{yamlFile}: empty or unparseable YAML')
    nUsed = data.get(usedKey, data.get('nTiepointsUsed'))
    return data['sigma'], nUsed, data['nTiepointsGiven']


def readRBaselineYaml(yamlFile):
    with open(yamlFile) as fp:
        data = yaml.safe_load(fp)
    if not data:
        raise ValueError(f'{yamlFile}: empty or unparseable YAML')
    return (data['sigma'], data['nTiepointsUsed'], data['nTiepointsGiven'],
           data['sigmaWithoutIonCorrection'], data['usingIon'])


def readLegacyAzSigma(azFile):
    with open(azFile) as fp:
        text = fp.read()
    sigmaMatch = sigmaLinePattern.search(text)
    tiepointsMatch = tiepointsLinePattern.search(text)
    if not sigmaMatch or not tiepointsMatch:
        raise ValueError(f'could not find sigma/Ntiepoints lines in '
                          f'{azFile}')
    return (float(sigmaMatch.group(1)), int(tiepointsMatch.group(1)),
            int(tiepointsMatch.group(2)))


def readAzEst(virtualFrameDir):
    '''az.est.const.yaml if present (matches project's yaml-output
    convention), else the legacy az.est.const text file.'''
    yamlFile = os.path.join(virtualFrameDir, 'motion', 'az.est.const.yaml')
    if os.path.exists(yamlFile):
        return readYamlSigma(yamlFile, 'nTiepointsUsed')
    legacyFile = os.path.join(virtualFrameDir, 'motion', 'az.est.const')
    if not os.path.exists(legacyFile):
        raise FileNotFoundError('motion/az.est.const(.yaml)')
    return readLegacyAzSigma(legacyFile)


def readVirtualFrame(virtualFrameDir, track, orbit):
    '''Read all required sidecar files for one virtual frame. Returns
    (record, []) on success, or (None, reasons) -- reasons is a list of
    human-readable strings describing every missing/unusable input -- if any
    required input is missing, unparseable, or internally inconsistent.'''
    framesFile, r1 = findRequired(virtualFrameDir, 'frames.txt')
    geodatFile, r2 = findRequired(virtualFrameDir, 'geodat*.geojson',
                                 pattern=re.compile(r'^geodat(?!.*secondary).*\.geojson$'))
    baselineFile, r3 = findRequired(virtualFrameDir, 'motion/baseline.*.yaml',
                                    pattern=baselineNamePattern)
    rBaselineFile, r4 = findRequired(virtualFrameDir, 'motion/rBaseline.deltab*.yaml',
                                     pattern=rBaselineNamePattern)
    pairInfoFile, r5 = findRequired(virtualFrameDir, '*.pairinfo',
                                    pattern=re.compile(r'^\d+\.\d+\.pairinfo$'))
    reasons = [r for r in (r1, r2, r3, r4, r5) if r]
    if reasons:
        return None, reasons
    try:
        frames = readFrames(framesFile)
        geometry, imageName = readGeodat(geodatFile)
        orbit1, orbit2 = readPairInfo(pairInfoFile)
        sigmaBaseline, nTb, nTgb = readYamlSigma(baselineFile, 'nTiepoints')
        (sigmaRBaseline, nTr, nTgr, sigmaRBaselineWithoutIon,
        usingIonRBaseline) = readRBaselineYaml(rBaselineFile)
        sigmaAz, nTa, nTga = readAzEst(virtualFrameDir)
        pDict, _ = parseFileName(imageName)
        cycle = int(pDict['cycle'])
        direction = directionNames[pDict['direction']]
    except (OSError, ValueError, KeyError, FileNotFoundError, AttributeError) as exc:
        reason = f'unparseable input: {exc}'
        print(f'{virtualFrameDir}: {reason}')
        return None, [reason]
    if orbit1 != orbit:
        reason = (f'pairinfo orbit1 {orbit1} does not match directory '
                  f'orbit {orbit}')
        print(f'{virtualFrameDir}: {reason}')
        return None, [reason]
    return {
        'track': track, 'orbit': orbit,
        'virtualFrameId': os.path.basename(virtualFrameDir),
        'orbit1': orbit1, 'orbit2': orbit2,
        'cycle': cycle, 'direction': direction, 'frames': frames,
        'sigmaBaseline': sigmaBaseline, 'sigmaRBaseline': sigmaRBaseline,
        'sigmaAz': sigmaAz,
        'sigmaRBaselineWithoutIon': sigmaRBaselineWithoutIon,
        'usingIonRBaseline': usingIonRBaseline,
        'nTiepointsBaseline': nTb, 'nTiepointsGivenBaseline': nTgb,
        'nTiepointsUsedRBaseline': nTr, 'nTiepointsGivenRBaseline': nTgr,
        'nTiepointsUsedAz': nTa, 'nTiepointsGivenAz': nTga,
        'geometry': geometry,
    }, []


def collectFrames(projectDir):
    records = []
    skipped = []  # [(virtualFrameDir, [reasons]), ...]
    framePattern = loadFramePattern(projectDir)
    virtualFrameDirs = sorted(
        glob.glob(os.path.join(projectDir, 'track-*', f'*_{framePattern}')))
    for virtualFrameDir in virtualFrameDirs:
        track = int(os.path.basename(os.path.dirname(virtualFrameDir))
                    .split('-')[1])
        orbit = int(os.path.basename(virtualFrameDir).split('_')[0])
        record, reasons = readVirtualFrame(virtualFrameDir, track, orbit)
        if record is None:
            skipped.append((virtualFrameDir, reasons))
            continue
        records.append(record)
    return records, len(virtualFrameDirs), skipped


def groupByCycle(records):
    byCycle = {}
    for record in records:
        byCycle.setdefault(record['cycle'], {'ascending': [],
                                             'descending': []})
        byCycle[record['cycle']][record['direction']].append(record)
    return byCycle


def writeCycleGpkg(gpkgPath, recordsByDirection):
    if os.path.exists(gpkgPath):
        os.remove(gpkgPath)
    driver = ogr.GetDriverByName('GPKG')
    ds = driver.CreateDataSource(gpkgPath)
    #
    srcSRS = osr.SpatialReference()
    srcSRS.ImportFromEPSG(4326)
    srcSRS.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dstSRS = osr.SpatialReference()
    dstSRS.ImportFromEPSG(3413)
    dstSRS.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(srcSRS, dstSRS)
    #
    for direction in ('ascending', 'descending'):
        layer = ds.CreateLayer(direction, srs=dstSRS, geom_type=ogr.wkbPolygon)
        for name, fieldType in fieldDefs:
            layer.CreateField(ogr.FieldDefn(name, fieldType))
        layerDefn = layer.GetLayerDefn()
        for record in recordsByDirection[direction]:
            feat = ogr.Feature(layerDefn)
            geom = ogr.CreateGeometryFromJson(json.dumps(record['geometry']))
            geom.Transform(transform)
            feat.SetGeometry(geom)
            for name, _ in fieldDefs:
                feat.SetField(name, record[name])
            layer.CreateFeature(feat)
    ds = None


def main():
    projectDir, outputDir = buildFrameGpkgArgs()
    if not os.path.isdir(projectDir):
        u.myerror(f'No such project directory: {projectDir}')
    records, nTotal, skipped = collectFrames(projectDir)
    if not records:
        u.myerror(f'No usable virtual frames found under {projectDir}')
    os.makedirs(outputDir, exist_ok=True)
    byCycle = groupByCycle(records)
    print(f'\nOutput dir: {outputDir}  (EPSG:3413)')
    for cycle in sorted(byCycle):
        gpkgPath = os.path.join(outputDir, f'cycle{cycle:02d}.gpkg')
        writeCycleGpkg(gpkgPath, byCycle[cycle])
        nAsc = len(byCycle[cycle]['ascending'])
        nDesc = len(byCycle[cycle]['descending'])
        print(f'  cycle{cycle:02d}.gpkg: {nAsc} ascending, {nDesc} descending')
    print(f'\n{len(records)} frame(s) written, {len(skipped)} of {nTotal} skipped')

    motionIssues = [(d, r) for d, reasons in skipped for r in reasons
                    if 'motion/' in r or 'baseline' in r.lower() or 'az.est' in r]
    if motionIssues:
        print(f'\nMotion directories missing/empty baseline files ({len(motionIssues)}):')
        for d, r in motionIssues:
            print(f'  {d}: {r}')
    print(f'\nYou can now run buildFrameLayers --inventoryDir {outputDir} '
          'to create a QGIS .qlr file for browsing these in QGIS.')


if __name__ == '__main__':
    main()
