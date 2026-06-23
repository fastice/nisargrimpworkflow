import glob
import utilities as u
import argparse
import subprocess

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
                        help='Pass --sepIceRock to SetupNISAR')
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
    parser.add_argument('-noPrompt', '--noPrompt', action='store_true',
                        help='Pass --noPrompt to SetupNISAR: skip the confirmation '
                        'prompt for --clean/--cleanDebug')
    args = parser.parse_args()
    #
    return args.track[0], args.overWrite, args.overWritePhase, args.RUNWOnly, \
        args.correlationOnly, args.debugIono, args.sepIceRock, args.geodatsOnly, \
        args.clean, args.cleanDebug, args.noPrompt


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
        sepIceRock, geodatsOnly, clean, cleanDebug, noPrompt = parseArgs()
    orbitDirs =  glob.glob(f'{track}/*_*')
    print(orbitDirs)
    orbits = sorted(list(set([x.split('/')[-1].split('_')[0] for x in orbitDirs])))
    print(orbits)
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
        #command += ['--verbose']
        print(command)
        subprocess.run(command, cwd=track)

if __name__ == '__main__':
    main()
