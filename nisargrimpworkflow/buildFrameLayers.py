#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ian
"""

import argparse
import glob
import os
import uuid

import utilities as u
from osgeo import ogr

defaultInventoryDirName = 'frameInventory'
defaultOffsetsSigmaThresh = 0.5
highSigmaOutlineColor = '230,30,30,255'

sigmaFields = ['sigmaBaseline', 'sigmaRBaseline', 'sigmaAz']
defaultSigmaField = 'sigmaRBaseline'

# 10 classes in 0.1 increments from 0, top class open-ended (catches > 1)
rBaselineClassBreaks = [(None, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4),
                        (0.4, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8),
                        (0.8, 0.9), (0.9, None)]
# ColorBrewer RdYlGn (10-class), green (low/good) -> red (high/bad)
rBaselineClassColors = ['0,104,55', '26,152,80', '102,189,99', '166,217,106',
                        '217,239,139', '254,224,139', '253,174,97',
                        '244,109,67', '215,48,39', '165,0,38']

# track, orbit, virtualFrameId, orbit1, orbit2, cycle, direction, frames,
# plus the 3 sigma/tiepoint triples -- must match buildFrameGpkg.py fieldDefs
frameFields = ['fid', 'track', 'orbit', 'virtualFrameId', 'orbit1', 'orbit2',
              'cycle', 'direction', 'frames', 'sigmaBaseline',
              'sigmaRBaseline', 'sigmaAz', 'sigmaRBaselineWithoutIon',
              'usingIonRBaseline', 'nTiepointsBaseline',
              'nTiepointsGivenBaseline', 'nTiepointsUsedRBaseline',
              'nTiepointsGivenRBaseline', 'nTiepointsUsedAz',
              'nTiepointsGivenAz']

epsg3413Srs = '''<srs>
        <spatialrefsys nativeFormat="Wkt">
          <wkt>PROJCRS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",BASEGEOGCRS["WGS 84",ENSEMBLE["World Geodetic System 1984 ensemble",MEMBER["World Geodetic System 1984 (Transit)"],MEMBER["World Geodetic System 1984 (G730)"],MEMBER["World Geodetic System 1984 (G873)"],MEMBER["World Geodetic System 1984 (G1150)"],MEMBER["World Geodetic System 1984 (G1674)"],MEMBER["World Geodetic System 1984 (G1762)"],MEMBER["World Geodetic System 1984 (G2139)"],ELLIPSOID["WGS 84",6378137,298.257223563,LENGTHUNIT["metre",1]],ENSEMBLEACCURACY[2.0]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4326]],CONVERSION["US NSIDC Sea Ice polar stereographic north",METHOD["Polar Stereographic (variant B)",ID["EPSG",9829]],PARAMETER["Latitude of standard parallel",70,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8832]],PARAMETER["Longitude of origin",-45,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8833]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["easting (X)",south,MERIDIAN[45,ANGLEUNIT["degree",0.0174532925199433]],ORDER[1],LENGTHUNIT["metre",1]],AXIS["northing (Y)",south,MERIDIAN[135,ANGLEUNIT["degree",0.0174532925199433]],ORDER[2],LENGTHUNIT["metre",1]],USAGE[SCOPE["Polar research."],AREA["Northern hemisphere - north of 60°N onshore and offshore, including Arctic."],BBOX[60,-180,90,180]],ID["EPSG",3413]]</wkt>
          <proj4>+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs</proj4>
          <srsid>1371</srsid>
          <srid>3413</srid>
          <authid>EPSG:3413</authid>
          <description>WGS 84 / NSIDC Sea Ice Polar Stereographic North</description>
          <projectionacronym>stere</projectionacronym>
          <ellipsoidacronym>EPSG:7030</ellipsoidacronym>
          <geographicflag>false</geographicflag>
        </spatialrefsys>
      </srs>'''

fieldConfigTemplate = '''        <field configurationFlags="None" name="{f}">
          <editWidget type="">
            <config>
              <Option/>
            </config>
          </editWidget>
        </field>
'''

simpleFillTemplate = '''          <symbol alpha="1" type="fill" clip_to_extent="1" is_animated="0" frame_rate="10" force_rhr="0" name="{symbolName}">
            <layer enabled="1" pass="0" class="SimpleFill" locked="0" id="{{{symbolLayerKey}}}">
              <Option type="Map">
                <Option type="QString" value="3x:0,0,0,0,0,0" name="border_width_map_unit_scale"/>
                <Option type="QString" value="190,190,190,255" name="color"/>
                <Option type="QString" value="bevel" name="joinstyle"/>
                <Option type="QString" value="0,0" name="offset"/>
                <Option type="QString" value="3x:0,0,0,0,0,0" name="offset_map_unit_scale"/>
                <Option type="QString" value="MM" name="offset_unit"/>
                <Option type="QString" value="{outlineColor}" name="outline_color"/>
                <Option type="QString" value="solid" name="outline_style"/>
                <Option type="QString" value="0.5" name="outline_width"/>
                <Option type="QString" value="MM" name="outline_width_unit"/>
                <Option type="QString" value="no" name="style"/>
              </Option>
{dataDefinedBlock}            </layer>
          </symbol>
'''

dataDefinedOutlineColorTemplate = '''              <data_defined_properties>
                <Option type="Map">
                  <Option type="QString" value="" name="name"/>
                  <Option type="Map" name="properties">
                    <Option type="Map" name="outlineColor">
                      <Option type="bool" value="true" name="active"/>
                      <Option type="QString" value="{colorExpr}" name="expression"/>
                      <Option type="int" value="3" name="type"/>
                    </Option>
                  </Option>
                  <Option type="QString" value="collection" name="type"/>
                </Option>
              </data_defined_properties>
'''

maplayerTemplate = '''    <maplayer geometry="Polygon" autoRefreshTime="0" legendPlaceholderImage="" wkbType="Polygon" autoRefreshEnabled="0" readOnly="0" refreshOnNotifyEnabled="0" minScale="100000000" symbologyReferenceScale="-1" maxScale="0" hasScaleBasedVisibilityFlag="0" styleCategories="AllStyleCategories" simplifyDrawingTol="1" type="vector" simplifyMaxScale="1" simplifyAlgorithm="0" simplifyLocal="1" simplifyDrawingHints="1" refreshOnNotifyMessage="" labelsEnabled="0">
      <id>{layerId}</id>
      <datasource>{datasource}</datasource>
      <keywordList>
        <value></value>
      </keywordList>
      <layername>{displayName}</layername>
      {srs}
      <provider encoding="UTF-8">ogr</provider>
      <vectorjoins/>
      <layerDependencies/>
      <dataDependencies/>
      <expressionfields/>
      <map-layer-style-manager current="default">
        <map-layer-style name="default"/>
      </map-layer-style-manager>
      <auxiliaryLayer/>
      <metadataUrls/>
      <flags>
        <Identifiable>1</Identifiable>
        <Removable>1</Removable>
        <Searchable>1</Searchable>
        <Private>0</Private>
      </flags>
{rendererXml}      <customproperties>
        <Option/>
      </customproperties>
      <blendMode>0</blendMode>
      <featureBlendMode>0</featureBlendMode>
      <layerOpacity>1</layerOpacity>
      <geometryOptions removeDuplicateNodes="0" geometryPrecision="0">
        <activeChecks/>
        <checkConfiguration/>
      </geometryOptions>
      <legend type="default-vector" showLabelLegend="0"/>
      <referencedLayers/>
      <fieldConfiguration>
{fieldConfigs}      </fieldConfiguration>
      <aliases>
{aliases}      </aliases>
      <splitPolicies>
{splitPolicies}      </splitPolicies>
      <defaults>
{defaults}      </defaults>
      <constraints>
{constraints}      </constraints>
      <constraintExpressions>
{constraintExpressions}      </constraintExpressions>
      <expressionfields/>
      <attributeactions/>
      <attributetableconfig sortOrder="0" actionWidgetStyle="dropDown" sortExpression="">
        <columns/>
      </attributetableconfig>
      <conditionalstyles>
        <rowstyles/>
        <fieldstyles/>
      </conditionalstyles>
      <storedexpressions/>
      <editform tolerant="1"></editform>
      <editforminit/>
      <editforminitcodesource>0</editforminitcodesource>
      <editforminitfilepath></editforminitfilepath>
      <editforminitcode><![CDATA[]]></editforminitcode>
      <featformsuppress>0</featformsuppress>
      <editorlayout>generatedlayout</editorlayout>
      <editable/>
      <labelOnTop/>
      <reuseLastValue/>
      <dataDefinedFieldProperties/>
      <widgets/>
      <previewExpression></previewExpression>
      <mapTip></mapTip>
    </maplayer>
'''


def buildFrameLayersArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mBuild a QGIS Layer Definition (.qlr) with '
        'a "Frames" group (sigma_field-switchable), an "rBaseline" group '
        '(fixed 10-class sigmaRBaseline coloring) -- both organized '
        'ascending/descending > Cycle N -- and a flat high-sigma flag group '
        '(ascending/descending only, not split by cycle) \033[0m\n\n'
        'Reads the per-cycle GeoPackages written by buildFrameGpkg; run '
        'that first.\n\n')
    parser.add_argument('--projectDir', type=str, default='.',
                        help='Root directory containing track-N '
                        'subdirectories [.] (current directory)')
    parser.add_argument('--inventoryDir', type=str, default=None,
                        help='Directory containing the per-cycle '
                        f'GeoPackages [<projectDir>/{defaultInventoryDirName}]')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output .qlr path [<inventoryDir>/frames.qlr]')
    parser.add_argument('--offsetsSigmaThresh', type=float,
                        default=defaultOffsetsSigmaThresh,
                        help='sigmaRBaseline threshold for the flat '
                        f'high-sigma flag group [{defaultOffsetsSigmaThresh}]')
    args = parser.parse_args()
    inventoryDir = args.inventoryDir or \
        os.path.join(args.projectDir, defaultInventoryDirName)
    output = args.output or os.path.join(inventoryDir, 'frames.qlr')
    return inventoryDir, output, args.offsetsSigmaThresh


def readCycleGpkgs(inventoryDir):
    gpkgPaths = sorted(glob.glob(os.path.join(inventoryDir, 'cycle*.gpkg')))
    if not gpkgPaths:
        u.myerror(f'No cycle*.gpkg files found in {inventoryDir}')
    cycles = {}
    for gpkgPath in gpkgPaths:
        cycleNum = int(os.path.basename(gpkgPath)[len('cycle'):-len('.gpkg')])
        # Resolve to absolute so the .qlr's embedded datasource paths remain
        # valid regardless of where QGIS's cwd is when the file is loaded.
        cycles[cycleNum] = os.path.abspath(gpkgPath)
    return cycles


def sigmaRangesAcrossCycles(cycleGpkgs):
    sigmaRanges = {f: [None, None] for f in sigmaFields}
    for gpkgPath in cycleGpkgs.values():
        ds = ogr.Open(gpkgPath)
        for direction in ('ascending', 'descending'):
            layer = ds.GetLayer(direction)
            for feat in layer:
                for f in sigmaFields:
                    val = feat.GetField(f)
                    lo, hi = sigmaRanges[f]
                    sigmaRanges[f] = [val if lo is None else min(lo, val),
                                     val if hi is None else max(hi, val)]
    return sigmaRanges


def buildSigmaFieldExpr(sigmaRanges, bound):
    whens = ' '.join(
        f"WHEN @sigma_field = '{f}' THEN {sigmaRanges[f][bound]:.6f}"
        for f in sigmaFields)
    return f'CASE {whens} END'


def xmlAttrEscape(s):
    return (s.replace('&', '&amp;').replace('"', '&quot;')
            .replace('<', '&lt;').replace('>', '&gt;'))


def buildContinuousRendererXml(sigmaRanges):
    '''Single rule, fixed outline color overridden per-feature by an
    expression keyed on the @sigma_field project variable.'''
    loExpr = buildSigmaFieldExpr(sigmaRanges, 0)
    hiExpr = buildSigmaFieldExpr(sigmaRanges, 1)
    colorExpr = (f"ramp_color('RdYlGn', scale_linear(attribute("
                f"$currentfeature, @sigma_field), {loExpr}, {hiExpr}, "
                f"1, 0))").replace('"', '&quot;')
    dataDefinedBlock = dataDefinedOutlineColorTemplate.format(
        colorExpr=colorExpr)
    symbolXml = simpleFillTemplate.format(
        symbolName='0', symbolLayerKey=uuid.uuid4(),
        outlineColor='35,35,35,255', dataDefinedBlock=dataDefinedBlock)
    return (f'      <renderer-v2 type="RuleRenderer" forceraster="0" '
            f'symbollevels="0" enableorderby="0" referencescale="-1">\n'
            f'        <rules key="{{{uuid.uuid4()}}}">\n'
            f'          <rule key="{{{uuid.uuid4()}}}" symbol="0"/>\n'
            f'        </rules>\n        <symbols>\n{symbolXml}'
            f'        </symbols>\n      </renderer-v2>\n')


def buildClassedRendererXml(field, breaks, colors):
    '''One rule per class, each with its own fixed outline color -- no
    data-defined expression, so the legend swatches are accurate.'''
    ruleLines, symbolLines = [], []
    for i, ((lo, hi), color) in enumerate(zip(breaks, colors)):
        if lo is None:
            filterExpr, label = f'"{field}" < {hi}', f'< {hi}'
        elif hi is None:
            filterExpr, label = f'"{field}" >= {lo}', f'>= {lo}'
        else:
            filterExpr = f'"{field}" >= {lo} AND "{field}" < {hi}'
            label = f'{lo} - {hi}'
        ruleLines.append(
            f'          <rule key="{{{uuid.uuid4()}}}" '
            f'filter="{xmlAttrEscape(filterExpr)}" symbol="{i}" '
            f'label="{xmlAttrEscape(label)}"/>\n')
        symbolLines.append(simpleFillTemplate.format(
            symbolName=str(i), symbolLayerKey=uuid.uuid4(),
            outlineColor=f'{color},255', dataDefinedBlock=''))
    return (f'      <renderer-v2 type="RuleRenderer" forceraster="0" '
            f'symbollevels="0" enableorderby="0" referencescale="-1">\n'
            f'        <rules key="{{{uuid.uuid4()}}}">\n'
            + ''.join(ruleLines) + '        </rules>\n        <symbols>\n'
            + ''.join(symbolLines) +
            '        </symbols>\n      </renderer-v2>\n')


def buildFieldBlocks():
    fieldConfigs = ''.join(fieldConfigTemplate.format(f=f)
                           for f in frameFields)
    aliases = ''.join(
        f'        <alias index="{i}" field="{f}" name=""/>\n'
        for i, f in enumerate(frameFields))
    splitPolicies = ''.join(
        f'        <policy policy="Duplicate" field="{f}"/>\n'
        for f in frameFields)
    defaults = ''.join(
        f'        <default applyOnUpdate="0" expression="" field="{f}"/>\n'
        for f in frameFields)
    constraints = ''.join(
        f'        <constraint notnull_strength="0" exp_strength="0" '
        f'constraints="0" field="{f}" unique_strength="0"/>\n'
        for f in frameFields)
    constraintExpressions = ''.join(
        f'        <constraint exp="" desc="" field="{f}"/>\n'
        for f in frameFields)
    return (fieldConfigs, aliases, splitPolicies, defaults, constraints,
           constraintExpressions)


def buildMapLayer(gpkgPath, cycleNum, direction, rendererXml, fieldBlocks):
    layerId = f'{direction}_cycle{cycleNum:02d}_{uuid.uuid4().hex}'
    datasource = f'{gpkgPath}|layername={direction}'
    return layerId, buildMapLayerXml(layerId, datasource,
                                     f'Cycle {cycleNum}', rendererXml,
                                     fieldBlocks)


def buildUnionMapLayer(vrtPath, direction, subsetFilter, rendererXml,
                       fieldBlocks):
    layerId = f'{direction}_highsigma_{uuid.uuid4().hex}'
    datasource = f'{vrtPath}|subset={subsetFilter}'
    return layerId, buildMapLayerXml(layerId, datasource, direction,
                                     rendererXml, fieldBlocks)


def buildMapLayerXml(layerId, datasource, displayName, rendererXml,
                     fieldBlocks):
    (fieldConfigs, aliases, splitPolicies, defaults, constraints,
    constraintExpressions) = fieldBlocks
    return maplayerTemplate.format(
        layerId=layerId, datasource=datasource, displayName=displayName,
        srs=epsg3413Srs, rendererXml=rendererXml, fieldConfigs=fieldConfigs,
        aliases=aliases, splitPolicies=splitPolicies, defaults=defaults,
        constraints=constraints,
        constraintExpressions=constraintExpressions)


def buildUnionVrt(cycleGpkgs, direction, vrtPath):
    '''Write an OGR VRT (vector) Union Layer combining the given direction's
    layer across every per-cycle GeoPackage, so it can be queried/filtered
    as a single flat layer instead of one layer per cycle.'''
    layerLines = ''.join(
        f'    <OGRVRTLayer name="cycle{cycleNum:02d}">\n'
        f'      <SrcDataSource>{gpkgPath}</SrcDataSource>\n'
        f'      <SrcLayer>{direction}</SrcLayer>\n'
        f'    </OGRVRTLayer>\n'
        for cycleNum, gpkgPath in sorted(cycleGpkgs.items()))
    vrt = (f'<OGRVRTDataSource>\n  <OGRVRTUnionLayer name="{direction}">\n'
          f'{layerLines}  </OGRVRTUnionLayer>\n</OGRVRTDataSource>\n')
    with open(vrtPath, 'w') as fp:
        fp.write(vrt)
    return os.path.abspath(vrtPath)


def buildFixedOutlineRendererXml(color):
    '''Single catch-all rule, fixed outline color, no fill -- for layers
    whose membership (e.g. a subset filter) already tells the whole story,
    with no further per-feature coloring needed.'''
    symbolXml = simpleFillTemplate.format(
        symbolName='0', symbolLayerKey=uuid.uuid4(),
        outlineColor=color, dataDefinedBlock='')
    return (f'      <renderer-v2 type="RuleRenderer" forceraster="0" '
            f'symbollevels="0" enableorderby="0" referencescale="-1">\n'
            f'        <rules key="{{{uuid.uuid4()}}}">\n'
            f'          <rule key="{{{uuid.uuid4()}}}" symbol="0"/>\n'
            f'        </rules>\n        <symbols>\n{symbolXml}'
            f'        </symbols>\n      </renderer-v2>\n')


def buildLayerTreeLayerXml(layerId, source, name, indent='          '):
    return (f'{indent}<layer-tree-layer legend_split_behavior="0" '
           f'providerKey="ogr" id="{layerId}" checked="Qt::Checked" '
           f'patch_size="-1,-1" expanded="1" source="{source}" '
           f'legend_exp="" name="{name}">\n'
           f'{indent}  <customproperties>\n{indent}    <Option/>\n'
           f'{indent}  </customproperties>\n'
           f'{indent}</layer-tree-layer>\n')


def buildTopGroup(groupName, cycleGpkgs, rendererXml, fieldBlocks,
                  mapLayerBlocks):
    groupBlocks = []
    for direction in ('ascending', 'descending'):
        layerTreeLines = []
        for cycleNum in sorted(cycleGpkgs):
            gpkgPath = cycleGpkgs[cycleNum]
            layerId, mapLayerXml = buildMapLayer(
                gpkgPath, cycleNum, direction, rendererXml, fieldBlocks)
            mapLayerBlocks.append(mapLayerXml)
            layerTreeLines.append(buildLayerTreeLayerXml(
                layerId, f'{gpkgPath}|layername={direction}',
                f'Cycle {cycleNum}'))
        groupBlocks.append(
            f'        <layer-tree-group checked="Qt::Checked" expanded="1" '
            f'groupLayer="" name="{direction}">\n'
            f'          <customproperties>\n            <Option/>\n'
            f'          </customproperties>\n' + ''.join(layerTreeLines) +
            f'        </layer-tree-group>\n')
    return (f'      <layer-tree-group checked="Qt::Checked" expanded="1" '
            f'groupLayer="" name="{groupName}">\n'
            f'        <customproperties>\n          <Option/>\n'
            f'        </customproperties>\n' + ''.join(groupBlocks) +
            f'      </layer-tree-group>\n')


def buildFlatHighSigmaGroup(groupName, cycleGpkgs, inventoryDir, thresh,
                            fieldBlocks, mapLayerBlocks):
    '''A group with one layer per direction (no per-cycle nesting) showing
    every frame across all cycles where sigmaRBaseline > thresh, via an OGR
    VRT Union Layer + subset filter rather than per-cycle layers.'''
    rendererXml = buildFixedOutlineRendererXml(highSigmaOutlineColor)
    layerTreeLines = []
    for direction in ('ascending', 'descending'):
        vrtPath = buildUnionVrt(
            cycleGpkgs, direction,
            os.path.join(inventoryDir, f'all{direction.capitalize()}.vrt'))
        # Raw (unescaped) for <datasource> element text -- literal '"'/'>'
        # are valid there. The XML-attribute-escaped form is needed
        # separately for <layer-tree-layer source="...">, since that's an
        # attribute value and a literal '"' would terminate it early.
        subsetFilter = f'"sigmaRBaseline" > {thresh}'
        layerId, mapLayerXml = buildUnionMapLayer(
            vrtPath, direction, subsetFilter, rendererXml, fieldBlocks)
        mapLayerBlocks.append(mapLayerXml)
        layerTreeLines.append(buildLayerTreeLayerXml(
            layerId, xmlAttrEscape(f'{vrtPath}|subset={subsetFilter}'),
            direction, indent='        '))
    return (f'      <layer-tree-group checked="Qt::Checked" expanded="1" '
            f'groupLayer="" name="{groupName}">\n'
            f'        <customproperties>\n          <Option/>\n'
            f'        </customproperties>\n' + ''.join(layerTreeLines) +
            f'      </layer-tree-group>\n')


def main():
    inventoryDir, output, offsetsSigmaThresh = buildFrameLayersArgs()
    cycleGpkgs = readCycleGpkgs(inventoryDir)
    sigmaRanges = sigmaRangesAcrossCycles(cycleGpkgs)
    fieldBlocks = buildFieldBlocks()
    #
    mapLayerBlocks = []
    framesGroup = buildTopGroup(
        'Frames', cycleGpkgs, buildContinuousRendererXml(sigmaRanges),
        fieldBlocks, mapLayerBlocks)
    rBaselineGroup = buildTopGroup(
        'rBaseline', cycleGpkgs,
        buildClassedRendererXml('sigmaRBaseline', rBaselineClassBreaks,
                               rBaselineClassColors),
        fieldBlocks, mapLayerBlocks)
    highSigmaGroupName = xmlAttrEscape(f'sigmaRBaseline > {offsetsSigmaThresh}')
    highSigmaGroup = buildFlatHighSigmaGroup(
        highSigmaGroupName, cycleGpkgs, inventoryDir, offsetsSigmaThresh,
        fieldBlocks, mapLayerBlocks)
    #
    qlr = (
        "<!DOCTYPE qgis-layer-definition>\n<qlr>\n"
        '  <layer-tree-group checked="Qt::Checked" expanded="1" '
        'groupLayer="" name="">\n'
        '    <customproperties>\n      <Option/>\n    </customproperties>\n'
        + framesGroup + rBaselineGroup + highSigmaGroup +
        '  </layer-tree-group>\n  <maplayers>\n' + ''.join(mapLayerBlocks) +
        '  </maplayers>\n</qlr>\n')
    with open(output, 'w') as fp:
        fp.write(qlr)
    #
    print(f'Layer definition written: {output}')
    print(f'  {len(cycleGpkgs)} cycle group(s): {sorted(cycleGpkgs)}')
    for f in sigmaFields:
        lo, hi = sigmaRanges[f]
        print(f'  {f}: range [{lo:.4f}, {hi:.4f}]')
    print('\nIn QGIS: Layer > Add Layer > Add Layer Definition File... and '
          'pick this .qlr -- you get three top-level groups:')
    print('  "Frames" and "rBaseline" -- each organized '
          '"ascending"/"descending" > one "Cycle N" layer per cycle.')
    print('  "Frames" -- color driven by the sigma_field project variable '
          f'(default {defaultSigmaField}). Set it once via Project '
          'Properties > Variables -- add a variable named "sigma_field" '
          f'with value one of {sigmaFields} -- then Apply; every layer in '
          'this group recolors together.')
    print('  "rBaseline" -- fixed 10-class coloring of sigmaRBaseline in '
          '0.1 increments (green=low/good to red=high/bad), top class '
          '>= 0.9 (covers values above 1). No variable needed.')
    print(f'  "sigmaRBaseline > {offsetsSigmaThresh}" -- flat '
          '"ascending"/"descending" layers (not split by cycle) flagging '
          'every frame across all cycles above the threshold. Change the '
          'threshold with --offsetsSigmaThresh.')


if __name__ == '__main__':
    main()
