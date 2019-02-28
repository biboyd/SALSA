import radialplots
import argparse
import logging

def main():
    parser = argparse.ArgumentParser(description='Make radial plots of column density vs. impact parameter using uniformly random sightlines.')
    parser.add_argument('filenames', metavar='FN_IN', type=str,
                        help='pattern matching for data dumps')
    parser.add_argument('num_procs', metavar='N_PROCS', type=int, nargs='?', default=1,
                        help='number of processes to parallelize over')
    parser.add_argument('--ions', '-i', metavar='IONS', type=lambda s: [ion for ion in s.split(',')],
                        help='list of ion species (form: \"ion1,ion2,ion3\")')
    parser.add_argument('--c_coord', '-cc', metavar='C', type=list,
                        help='coordinates of profile\'s center')
    parser.add_argument('--c_units', '-cu', metavar='C_UNITS', type=str,
                        help='units of center coordinates')
    parser.add_argument('--r_len', '-rl', metavar='R', type=float,
                        help='length of radius')
    parser.add_argument('--r_units', '-ru', metavar='R_UNITS', type=str,
                        help='units of radius')
    parser.add_argument('--s_len', '-sl', metavar='S', type=float,
                        help='length of sightline')
    parser.add_argument('--s_units', '-su', metavar='S_UNITS', type=str,
                        help='units of sightline')
    parser.add_argument('--num_lines', '-cnl', metavar='CNL', type=int,
                        help='number of lines within central shell')
    parser.add_argument('--num_bins', '-b', metavar='BINS', type=int,
                        help='number of radial bins')
    parser.add_argument('--ylim', '-Y', metavar='YLIM', type=list,
                        help='y limits for cdip plot')
    parser.add_argument('output_filename_head', metavar='FN_OUT_HEAD', type=str,
                        help='filename header for each output image')

    args = parser.parse_args()
    
    logging.basicConfig(filename='%scdip_uniform.log' % args.output_filename_head, level=logging.DEBUG, 
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger=logging.getLogger(__name__)
    
    # example command:
    # python makeradialplots.py "../../../sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD????/DD????" 2 '../../img_dir/test_directory/radprofile' --ions 'O IV,H I' --r_units 'kpc' --r_len 300
    try:
        print("pit stop 0")
        radialplots.create_cdip_uniform(filenames = args.filenames,
                                         num_procs = args.num_procs,
                                         ion_species = args.ions, 
                                         center_coord = args.c_coord, 
                                         center_units = args.c_units, 
                                         radius_length = args.r_len, 
                                         radius_units = args.r_units, 
                                         sightline_length = args.s_len, 
                                         sightline_units = args.s_units, 
                                         center_num_lines = args.num_lines, 
                                         num_bins = args.num_bins,
                                         y_lim1 = args.ylim,
                                         output_filename_head=args.output_filename_head)
        print("pit stop 1")
    except Exception as e:
        print("pit stop something")
        logger.error(e)

if __name__ == '__main__':
    main()
