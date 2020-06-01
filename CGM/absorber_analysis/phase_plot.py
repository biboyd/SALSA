import yt
import argparse

from CGM.general_utils.filter_definitions import hist_range_dict
homedir="/mnt/home/boydbre1/"
scratch="/mnt/gs18/scratch/users/boydbre1/"

def main(args):
    dsName=args.ds
    xfld = args.x
    yfld = args.y
    zfld = args.z

    ion = args.ion
    ref = args.refinement
    max_impact = args.m
    cut = args.cut

    yt.enable_parallelism()
    ion_u = "_".join(ion.split())
    cut_u = "_".join(cut.split())

    #setup file paths
    main_dir=f"{homedir}/data/absorber_data/{ref}_refinement/max_impact200/ion_{ion_u}/"
    out_dir = f"{main_dir}/plots/{cut_u}/phase_plots"
    out_file = f"{out_dir}/{xfld}_x_{yfld}_x_{zfld}.png"

    ds_file=f"{scratch}/cosmological/{ref}_refinement/{dsName}/{dsName}"
    points_file=f"{main_dir}/{cut_u}/{dsName}_positions.npy"

    #create data object of cells
    cell_abs = get_absorber_cells(points_file, ds_file)

    #create and save phase plot 
    range_dict=hist_range_dict[ion]
    make_phase_plot(cell_abs, xfld, yfld, zfld,range_dict=range_dict, outf=out_file)

def get_absorber_cells(coordinates_file, ds_file):
    ds = yt.load(ds_file)
    coords = np.load(coordinates_file)
    coords = ds.arr(coords, 'code_length')

    points=[]
    for co in coords:
        points.append(ds.point(co))

    absorber_cells = ds.union(points, ds=ds)
    return absorber_cells

def make_phase_plot(absorber_cells, xvar, yvar, zvar, range_dict={},cmap='viridis',outf="phase_plot.png"):
    phase = yt.PhasePlot(absorber_cells,xvar, yvar, zvar)

    #limits and stuff
    if xvar in range_dict.keys():
        xlim1, xlim2=range_dict[xvar]
        phase.set_xlim(xlim1, xlim2)

    if yvar in range_dict.keys():
        ylim1, ylim2=range_dict[yvar]
        phase.set_ylim(ylim1, ylim2)

    if zvar in range_dict.keys():
        zlim1, zlim2=range_dict[zvar]
        phase.set_zlim(zlim1, zlim2)

    phase.save(outf)
    return phase

if __name__ == '__main__':

    #create parser
    parser = argparse.ArgumentParser(description='Process what fields to plot and stuff')
    parser.add_argument("--ds", type=str, help="The dataset name (ie RD0020)",
                        required=True)
    parser.add_argument("-x", type=str, help="The x-field to plot ",
                        required=True)
    parser.add_argument("-y", type=str, help="The y-field to plot",
                        required=True)
    parser.add_argument("-z", type=str, help="The z-field to plot",
                        required=True)
    parser.add_argument("-i", "--ion", type=str,
                        help='The ion to look at (ie "H I", "O VI")',
                        default="O VI")
    parser.add_argument("-r", "--refinement", type=str,
                        help="Refinement scheme, cool or nat",
                        choices=['cool', 'nat'], default='cool')
    parser.add_argument("-m", "--max-impact", type=int,
                        help="Max impact parameter sampled", default=200)
    parser.add_argument("-c", "--cut", type=str, default='cgm')

    args = parser.parse_args()

    main(args)
