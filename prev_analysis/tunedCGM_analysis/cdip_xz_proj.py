import radialplots
import argparse
import logging

def main():
    parser = argparse.ArgumentParser(description='Make radial plots of column density vs. impact parameter.')
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
    parser.add_argument('--px_num', '-pn', metavar='PX_NUM', type=int,
                        help='number of pixels per edge (so that image is PX_NUM x PX_NUM)')
    parser.add_argument('--num_bins', '-b', metavar='BINS', type=int,
                        help='number of radial bins')
    parser.add_argument('--ylimx', '-Yx', metavar='YLIM_x', type=list,
                        help='y limits for x projection plot')
    parser.add_argument('--ylimz', '-Yz', metavar='YLIM_z', type=list,
                        help='y limits for z projection plot')
    parser.add_argument('--ylimxz', '-Yxz', metavar='YLIM_xz', type=list,
                        help='y limits for plot of x- and z- projection medians')
    parser.add_argument('output_filename_head', metavar='FN_OUT_HEAD', type=str,
                        help='filename header for each output image')

    args = parser.parse_args()
    
    logging.basicConfig(filename='%scdip_xz_proj.log' % args.output_filename_head, level=logging.DEBUG, 
                        format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logger=logging.getLogger(__name__)
    
    

    # example command:
    # python makeradialplots.py "../../../sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD????/DD????" 2 '../../img_dir/test_directory/radprofile' --ions 'O IV,H I' --r_units 'kpc' --r_len 300
    try:
        radialplots.create_cdip_xz_proj(filenames = args.filenames,
                                         num_procs = args.num_procs,
                                         ion_species = args.ions, 
                                         center_coord = args.c_coord, 
                                         center_units = args.c_units, 
                                         radius_length = args.r_len, 
                                         radius_units = args.r_units, 
                                         pix_num = args.px_num, 
                                         num_bins = args.num_bins,
                                         y_lim1 = args.ylimx, 
                                         y_lim2 = args.ylimz, 
                                         y_lim3 = args.ylimxz, 
                                         output_filename_head=args.output_filename_head)
    except Exception as e:
        logger.error(e)

if __name__ == '__main__':
    main()
