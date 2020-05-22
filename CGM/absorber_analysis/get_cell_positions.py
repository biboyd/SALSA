import numpy as np
import yt
from astropy.table import Table

from CGM.general_utils.filter_definitions import parse_cut_filter

# assume we have table, will handle loading above

def main(table_file,raydir, cut_str, outdir):

    table = Table.read(table_file)

    positions =[]
    for index, start, end in table[['absorber_index', 'interval_start', 'interval_end']]:
        ray_num=index[:-1]

        ray = yt.load(f"{raydir}/ray{ray_num}.h5")

        cut_filters = parse_cut_filter(cut_str)

        data = ray.all_data()
        for fil in cut_filters:
            data = curr_data.cut_region(fil)

        coords, =np.dstack([data['x'][s:e].in_units('code_length'),
                            data['y'][s:e].in_units('code_length'),
                            data['z'][s:e].in_units('code_length')])

        positions.append(coords)

        ray.close()

    all_positions = np.vstack(positions)

    return all_positions


if __name__ == '__main__':
    ds=sys.argv[1]

    data_dir="/mnt/home/boydbre1/data/absorber_data/cool_refinement/ion_O_VI/"
    table_file = f"{data_dir}/cgm/{ds}_absorbers.h5"
    raydir=f"/mnt/gs18/scratch/users/boydbre1/analysis/cool_refinement/{ds}/max_impact200/rays/"

    pos = main(table_file, raydir, "cgm", data_dir)
    np.save(f"{data_dir}/{ds}_positions.npy", pos)
