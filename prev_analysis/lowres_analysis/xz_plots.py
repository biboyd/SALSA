print("try importing radialplots")
import radialplots
print("successfully imported radialplots")
import argparse
import logging

def main():
    print("Entered main in xz_plots.py")
    parser = argparse.ArgumentParser(description='Make xz plots for all datadumps.')
    parser.add_argument('filenames', metavar='FN_IN', type=str,
                        help='pattern matching for data dumps')
    parser.add_argument('num_procs', metavar='N_PROCS', type=int, nargs='?', default=1,
                        help='number of processes to parallelize over')
    parser.add_argument('--ions', '-i', metavar='IONS', type=lambda s: [ion for ion in s.split(',')],
                        help='list of ion species (form: \"ion1,ion2,ion3\")')
    parser.add_argument('--c_coord', '-cc', metavar='C', type=list,
                        help='coordinates of profile\'s center')
    parser.add_argument('--ylim1', '-y1', metavar='YLIM1', type=float,
                        help='lower y lim for plot')
    parser.add_argument('--ylim2', '-y2', metavar='YLIM2', type=float,
                        help='upper y lim for plot')
    parser.add_argument('output_filename_head', metavar='FN_OUT_HEAD', type=str,
                        help='filename header for each output image')

    args = parser.parse_args()
    
    logging.basicConfig(filename='%scombine_xz_plots.log' % args.output_filename_head, level=logging.DEBUG, 
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger=logging.getLogger(__name__)
    
    

    # example command:
    # python makeradialplots.py "../../../sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD????/DD????" 2 '../../img_dir/test_directory/radprofile' --ions 'O IV,H I' --r_units 'kpc' --r_len 300
    try:
        print("Trying to run function")
        radialplots.combine_xz_plots(filenames = args.filenames,
                                         num_procs = args.num_procs,
                                         ion_species = args.ions, 
                                         center = args.c_coord, 
                                         ylim1x = args.ylim1,
                                         ylim2x = args.ylim2,
                                         output_filename_head=args.output_filename_head)
    except Exception as e:
        print("Error running function")
        logger.error(e)

if __name__ == '__main__':
    main()
