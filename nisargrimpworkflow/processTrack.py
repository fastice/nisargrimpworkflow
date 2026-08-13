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
    
    parser.add_argument('--overWrite', action="store_true",
                        help='Overwrite if products already exist')
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
                        help='Pass --new to SetupNISAR: skip any virtual frame whose '
                        'product VRT already exists; only build new virtual frames')
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
    return args.track[0], args.overWrite, args.overWritePhase, args.RUNWOnly, \
        args.correlationOnly, args.debugIono, args.sepIceRock, args.geodatsOnly, \
        args.clean, args.cleanDebug, args.noPrompt, args.bakeOnly, args.new, \
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
    track, overWrite, overWritePhase, RUNWOnly, correlationOnly, debugIono, \
        sepIceRock, geodatsOnly, clean, cleanDebug, noPrompt, bakeOnly, new, \
        threads, firstDate, lastDate = parseArgs()
    orbitDirs =  glob.glob(f'{track}/*_*')
    print(orbitDirs)
    orbits = sorted(list(set([x.split('/')[-1].split('_')[0] for x in orbitDirs])))
    print(orbits)
    orbitThreads = []
    for orbit in orbits:
        if 'tie' in orbit:
            continue
        print(track)
        command =['SetupNISAR', f'{orbit}']
        if RUNWOnly:
            command += ['--RUNWOnly']
        if overWrite:
            command += ['--overWrite']
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
        if new:
            command += ['--new']
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
