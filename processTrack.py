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
    args = parser.parse_args()
    #
    return args.track[0], args.overWrite, args.overWritePhase, args.RUNWOnly


def main():
    '''
    Organize a directory full of test products into orbit_frame products in
    GrIMP format.

    Returns
    -------
    None.

    '''
    # Get args
    track, overWrite, overWritePhase, RUNWOnly  = parseArgs()
    orbitDirs =  glob.glob(f'{track}/*_*')
    print(orbitDirs)
    orbits = sorted(list(set([x.split('/')[-1].split('_')[0] for x in orbitDirs])))
    print(orbits)
    for orbit in orbits:
        if 'tie' in orbit:
            continue
        print(track)
        command =['SetupNISAR.py', f'{orbit}']
        if RUNWOnly:
            command += ['--RUNWOnly']
        if overWrite:
            command += ['--overWrite']
        if overWritePhase:
            command += ['--overWritePhase']
        #command += ['--verbose']                        
        print(command)
        subprocess.run(command, cwd=track)

main()
