import yt
import argparse
import numpy as np
from os import makedirs

from CGM.general_utils.filter_definitions import hist_range_dict, default_units_dict
homedir="/mnt/home/boydbre1/"
scratch="/mnt/gs18/scratch/users/boydbre1/"

def main(args):
    dsName=args.ds
    xfld = args.x
    yfld = args.y
    zfld = args.z

    ion = args.ion
    ref = args.refinement
    max_impact = args.max_impact
    cut = args.cut
    cmap = args.cmap

    yt.enable_parallelism()
    ion_u = "_".join(ion.split())
    cut_u = "_".join(cut.split())

    main_dir=f"{homedir}/data/absorber_data/{ref}_refinement/max_impact200/ion_{ion_u}/"

    ds_file=f"{scratch}/cosmological/{ref}_refinement/{dsName}/{dsName}"
    points_file=f"{main_dir}/{cut_u}/{dsName}_positions.npy"
    ds = yt.load(ds_file)

    #setup file paths
    out_dir = f"{main_dir}/plots/{cut_u}/phase_plots"
    out_file = f"{out_dir}/{xfld}_x_{yfld}_x_{zfld}_{ds.current_redshift:.2f}.png"

    makedirs(out_dir, exist_ok=True)
    #create data object of cells
    cell_abs = get_absorber_cells(points_file, ds_file)

    #create and save phase plot 
    range_dict=hist_range_dict[ion]

    make_phase_plot(cell_abs, xfld, yfld, zfld,range_dict=range_dict, cmap=cmap, outf=out_file)

def get_absorber_cells(coordinates_file, ds_file, center=None, bv=None):
    ds = yt.load(ds_file)
    coords = np.load(coordinates_file)
    coords = ds.arr(coords, 'code_length')

    points=[]
    for co in coords:
        points.append(ds.point(co))

    absorber_cells = ds.union(points, ds=ds)
    
    if center is not None:
        absorber_cells.set_field_parameter('radius', center)

    if bv is not None:
        absorber_cells.set_field_parameter('bulk velocity', bv)

    return absorber_cells

def make_phase_plot(absorber_cells, xvar, yvar, zvar, range_dict={}, cmap='viridis',outf="phase_plot.png"):
    
    linear=['radius', 'radial_velocity']

    if xvar in linear:
        xlog=False
    else:
        xlog=True

    if yvar in linear:
        ylog=False
    else:
        ylog=True

    lims={}
    if xvar in range_dict.keys():
        if xvar in default_units_dict:
            lims[xvar] = ((range_dict[xvar][0], default_units_dict[xvar]), (range_dict[xvar][1], default_units_dict[xvar]))
        else:
            lims[xvar] = range_dict[xvar]
    
    if yvar in range_dict.keys():
        if yvar in default_units_dict:
            lims[yvar] = (range_dict[yvar], default_units_dict[yvar])
            lims[yvar] = ((range_dict[yvar][0], default_units_dict[yvar]), (range_dict[yvar][1], default_units_dict[yvar]))
        else:
            lims[yvar] = range_dict[yvar]

    if zvar in range_dict.keys():
        if zvar in default_units_dict:
            lims[zvar] = ((range_dict[zvar][0], default_units_dict[zvar]), (range_dict[zvar][1], default_units_dict[zvar]))

        else:
            lims[zvar] = range_dict[zvar]
    if lims == {}:
        lims=None

    unit_dict={}
    #setting units 
    """if xvar in default_units_dict:
        unit_dict[xvar] = default_units_dict[xvar]

    if yvar in default_units_dict:
        unit_dict[yvar] = default_units_dict[yvar]

    if zvar in default_units_dict:
        unit_dict[zvar] = default_units_dict[zvar]
    """
    if unit_dict == {}:
        unit_dict=None
    
    print(xlog, ylog)
    profile = yt.create_profile(absorber_cells, [xvar, yvar], zvar, logs={xvar:xlog, yvar:ylog}, units=unit_dict, extrema=lims)
    phase = yt.PhasePlot.from_profile(profile)


    phase.set_cmap(zvar, cmap)
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
    parser.add_argument("--cmap", type=str, default='viridis')

    args = parser.parse_args()

    main(args)
