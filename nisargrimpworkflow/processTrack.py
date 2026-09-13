import glob
import threading
import utilities as u
import argparse
import datetime
import subprocess


def _isoDateStr(x):
    '''argparse type: validate YYYY-MM-DD, return the original string
    (passed through verbatim to SetupNISAR).'''
    try:
        datetime.datetime.strptime(x, '%Y-%m-%d')
    except ValueError:
        raise argparse.ArgumentTypeError(f"Date must be YYYY-MM-DD, got '{x}'")
    return x

def parseArgs():
    '''
    Handle command line args
    '''
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mMerge offsets from frames into a single '
        'product \033[0m\n\n',
        epilog='Part of the nisargrimpworkflow package.')
    parser.add_argument('track', type=str, nargs=1,
                       help='Path to input products')
    
    parser.add_argument('--overWrite', type=str, default=None,
                        metavar='all|ORBITS',
                        help='Rebuild existing products. REQUIRES an explicit scope: '
                        'either the word "all" (every orbit in the track) or a '
                        'comma/space separated list of orbits, e.g. "--overWrite 5430" '
                        'or "--overWrite 5430,5431". With a list, ONLY those orbits are '
                        'processed, and each is done in full (every virtual-frame '
                        'section of the orbit). Without --overWrite the default is '
                        '--new: every orbit is visited but only virtual frames that do '
                        'not already have a product VRT are built.')
    parser.add_argument('--RUNWOnly', action="store_true",
                        help='RUNWonly')    
    parser.add_argument('--overWritePhase', action="store_true",
                        help='Overwrite if phase products already exist')
    parser.add_argument('--correlationOnly', action='store_true',
                        help='Process correlation products only (no phase/offsets)')
    parser.add_argument('--debugIono', action='store_true',
                        help='Pass --debugIono to SetupNISAR')
    parser.add_argument('--sepIceRock', action='store_true',
                        help='Force --sepIceRock on for every SetupNISAR call. '
                        'Normally left unset: SetupNISAR reads the sepIceRock '
                        'default from the project.yaml key (on for Greenland, off '
                        'for Antarctica), so this flag is only an override.')
    parser.add_argument('--geodatsOnly', action='store_true',
                        help='Pass --geodatsOnly to SetupNISAR: re-merge '
                        'virtual-frame geodats only, no reprocessing')
    parser.add_argument('--clean', action='store_true',
                        help='Pass --clean to SetupNISAR: remove computed output files '
                        '(everything --overWrite would replace) for every orbit found '
                        'in track, then exit without processing')
    parser.add_argument('--cleanDebug', action='store_true',
                        help='Pass --cleanDebug to SetupNISAR: empty the contents of '
                        'all debug/ directories (leaving the empty directory), then exit')
    parser.add_argument('--bakeOnly', action='store_true',
                        help='Pass --bakeOnly to SetupNISAR: just bake existing '
                        'virtual-frame VRTs into flat GeoTIFFs, skip all reprocessing')
    parser.add_argument('--new', action='store_true',
                        help='DEPRECATED and now the default: only build virtual frames '
                        'whose product VRT does not already exist. Accepted so existing '
                        'command lines keep working; it changes nothing. Use --overWrite '
                        'to rebuild existing products.')
    parser.add_argument('-noPrompt', '--noPrompt', action='store_true',
                        help='Pass --noPrompt to SetupNISAR: skip the confirmation '
                        'prompt for --clean/--cleanDebug')
    parser.add_argument('--firstDate', type=_isoDateStr, default=None,
                        metavar='YYYY-MM-DD',
                        help='Pass --firstDate to SetupNISAR: only process pairs '
                        'whose first (reference) acquisition date is on or after '
                        'this date. Omit for no lower bound.')
    parser.add_argument('--lastDate', type=_isoDateStr, default=None,
                        metavar='YYYY-MM-DD',
                        help='Pass --lastDate to SetupNISAR: only process pairs '
                        'whose first (reference) acquisition date is on or before '
                        'this date. Omit for no upper bound (infinity).')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of orbits to process concurrently '
                        '(separate from -ompThreads, which controls OpenMP '
                        'threads within each C binary call); default 1 '
                        '(orbits run one at a time)')
    args = parser.parse_args()
    #
    # --overWrite is scoped rather than a bare switch: it is either "all" or an explicit
    # orbit list. overWriteAll says whether every orbit is rebuilt; overWriteOrbits, when
    # non-empty, restricts the run to exactly those orbits. Neither set means the default,
    # which is --new (build only virtual frames that have no product VRT yet).
    overWriteAll, overWriteOrbits = False, []
    if args.overWrite is not None:
        if args.overWrite.strip().lower() == 'all':
            overWriteAll = True
        else:
            overWriteOrbits = [o for o in args.overWrite.replace(',', ' ').split() if o]
            if not overWriteOrbits:
                parser.error('--overWrite needs "all" or a list of orbits, '
                             f"got '{args.overWrite}'")
    return args.track[0], overWriteAll, overWriteOrbits, args.overWritePhase, \
        args.RUNWOnly, args.correlationOnly, args.debugIono, args.sepIceRock, \
        args.geodatsOnly, args.clean, args.cleanDebug, args.noPrompt, args.bakeOnly, \
        args.threads, args.firstDate, args.lastDate


def main():
    '''
    Organize a directory full of test products into orbit_frame products in
    GrIMP format.

    Returns
    -------
    None.

    '''
    # Get args
    track, overWriteAll, overWriteOrbits, overWritePhase, RUNWOnly, correlationOnly, \
        debugIono, sepIceRock, geodatsOnly, clean, cleanDebug, noPrompt, bakeOnly, \
        threads, firstDate, lastDate = parseArgs()
    orbitDirs =  glob.glob(f'{track}/*_*')
    orbits = sorted(list(set([x.split('/')[-1].split('_')[0] for x in orbitDirs])))
    orbits = [o for o in orbits if 'tie' not in o]
    #
    # An explicit orbit list restricts the run to those orbits, so --overWrite cannot
    # reach healthy orbits by accident. Every section of a listed orbit is rebuilt,
    # because SetupNISAR works a whole orbit at a time.
    if overWriteOrbits:
        unknown = [o for o in overWriteOrbits if o not in orbits]
        if unknown:
            print(f'processTrack: WARNING orbit(s) not found in {track}: '
                  f'{", ".join(unknown)}')
        orbits = [o for o in orbits if o in overWriteOrbits]
        if not orbits:
            print(f'processTrack: none of the requested orbits exist in {track}, '
                  'nothing to do')
            return
    if overWriteAll:
        mode = 'overWrite ALL orbits'
    elif overWriteOrbits:
        mode = f'overWrite orbits {",".join(overWriteOrbits)}'
    else:
        mode = 'new only (default -- existing products are left alone)'
    print(f'processTrack: {track}  orbits {",".join(orbits)}  mode: {mode}')
    if firstDate or lastDate:
        print(f'processTrack: date window [{firstDate or "-inf"}, {lastDate or "+inf"}]'
              ' -- orbits outside it are skipped by SetupNISAR')
    orbitThreads = []
    for orbit in orbits:
        command =['SetupNISAR', f'{orbit}']
        if RUNWOnly:
            command += ['--RUNWOnly']
        if overWriteAll or orbit in overWriteOrbits:
            command += ['--overWrite']
        else:
            # default: only build virtual frames that have no product VRT yet
            command += ['--new']
        if overWritePhase:
            command += ['--overWritePhase']
        if correlationOnly:
            command += ['--correlationOnly']
        if debugIono:
            command += ['--debugIono']
        if sepIceRock:
            command += ['--sepIceRock']
        if geodatsOnly:
            command += ['--geodatsOnly']
        if clean:
            command += ['--clean']
        if cleanDebug:
            command += ['--cleanDebug']
        if noPrompt:
            command += ['--noPrompt']
        if bakeOnly:
            command += ['--bakeOnly']
        if firstDate:
            command += ['--firstDate', firstDate]
        if lastDate:
            command += ['--lastDate', lastDate]
        #command += ['--verbose']
        print(command)
        orbitThreads.append(threading.Thread(target=subprocess.run,
                                             args=[command],
                                             kwargs={'cwd': track}))
    u.runMyThreads(orbitThreads, threads, 'processTrack')

if __name__ == '__main__':
    main()
