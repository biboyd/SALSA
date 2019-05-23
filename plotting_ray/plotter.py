import yt
import trident
import numpy as np
from sys import argv
import h5py
from os import remove, listdir, makedirs
import errno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy.linalg import norm

from scipy.constants import centi, kilo, parsec

class multi_plot():
    """
    Plots three images side by side for easy analysis. The images are:
    A slice of dataset along lightray's path.
    The number density of the given ion along that path.
    The artificial absorption spectra that would be observed.
    """
    #define conversion factor
    cm_to_kpc = centi/(parsec * kilo)

    def __init__(self,
                ds_filename,
                ray_filename,
                ion_name='H I',
                slice_field='None',
                absorber_fields=[],
                wavelength_center='None',
                wavelength_width = 300,
                resolution = 0.1,
                open=True,
                figure='None'):
        """
        init file names and ion name

        Parameters:
        ds_filename : Path/name of the enzo dataset to be loaded
        ray_filename : Path/name of the hdf5 ray file to be loaded

        ion_name : Name of the ion to plot in number density plot
        slice_field : Field to plot in slice plot. defaults to ion_name's number density
        absorber_fields : Additional ions to include in plots/Spectra, enter as list
        wavelength_center : Wavelength to center spectrum plot on. defaults to
                            a known spectral line of ion_name. in units of Angstrom
        wavelength_width : sets the wavelenth range of the spectrum plot. defaults
                            to 300 Angstroms
        resolution : width of wavelenth bins in spectrum plot. default 0.1 Angstrom
        figure : matplotlib figure where the multiplot will be plotted. creates one if
                none is specified.
        open : option on whether to immediately open dataset and ray files. defaults to True.
        ###NOTE### ion names should be in notaion:
              Element symbol *space* roman numeral of ion level (i.e. "H I", "O VI")
        """
        #set file names and ion name
        self.ds_filename = ds_filename
        self.ray_filename = ray_filename
        self.ion_name = ion_name

        #add ion name to list of all ions to be plotted
        absorber_fields.append(ion_name)
        self.ion_list = absorber_fields

        #set a value for slice
        self.slice = 'None'

        #set slice field to ion name if no field is specified
        if (slice_field == 'None'):
            self.slice_field = self.ion_p_name() + "_number_density"
        else:
            self.slice_field = slice_field

        self.wavelength_width = wavelength_width
        self.resolution = resolution
        #default set the wavelenth center to one of the known spectral lines
        #for ion name. Use tridents line database to search for correct wavelenth
        if (wavelength_center == "None"):
            #open up tridents default line database
            lbd = trident.LineDatabase("lines.txt")
            #find all lines that match ion
            lines = lbd.parse_subset(subsets= [self.ion_name])
            #take first one and return its wavelenth
            self.wavelenth_center = lines[0].wavelength

        #open up a figure if none specified
        if (figure == 'None'):
            self.fig = plt.figure(figsize=(10, 10))

        #open up the dataset and ray files or set their values to None
        if (open):
            self.ds, self.ray, self.ray_h5 = self.open_files()
        else:
            self.ds = 'None'
            self.ray = 'None'
            self.ray_h5 = 'None'

        #optionally set min/max value for number density plot
        self.num_dense_min = 'None'
        self.num_dense_max = 'None'


    def open_files(self, open_ds=True, open_ray=True):
        """
        Opens dataset and ray files int yt and h5py

        Parameters:
            open_ds : bool operator. whether or not the dataset should be opened
            open_ray :bool operator. whether or not the ray should be opened

        Returns:
        ds : dataset object loaded into yt or None if open_ds = False
        ray : ray object loaded into yt or None if open_ray = False
        ray_h5file : ray object loaded into h5py or None if open_ray = False
        """
        #open ds if open_ds = True
        ds = yt.load(self.ds_filename) if open_ds else 'None'

        #open ray file in its two forms if open_ray = True
        if (open_ray):
            ray = yt.load(self.ray_filename)
            ray_h5file = h5py.File(self.ray_filename)
        else:
            ray = 'None'
            ray_h5file = 'None'

        return ds, ray, ray_h5file

    def ion_p_name(self):
        """
        convert ion species name to yt style field convention
        ('H I' -> 'H_p0')
        Returns:
        outname : Name of the ion in the yt form
        """

        ######### Deprecated, no longer needed #########3
        # 'H I' is an exception just return H

        #if self.ion_name == 'H I':
        #    return 'H'

        #split up the words in ion species name
        ion_split = self.ion_name.split()

        #convert num from roman numeral. subtract one to follow convention
        num = trident.from_roman(ion_split[1])-1

        #combine all the names
        outname = ion_split[0] + '_p' + str(num)
        return outname

    def create_slice(self, cmap="magma"):
        """
        Create a slice in the Dataset along the path of the ray.
        Choose to keep the Z direction maintained.

        Parameters:
        field: The yt field to plot for the slice
        cmap='magma' : the colormap to use for the slice

        Returns:
        slice : yt SlicePlot with ray annotated
        """

        #add ion fields to dataset if not already there
        trident.add_ion_fields(self.ds, ions=self.ion_list, ftype='gas')

        # get beginning and end of ray
        x = self.ray_h5['grid']['x']
        y = self.ray_h5['grid']['y']
        z = self.ray_h5['grid']['z']

        ray_begin = np.array([ x[0], y[0], z[0] ])
        ray_end = np.array([ x[-1], y[-1], z[-1] ])

        #construct vector pointing in ray's direction
        ray_vec = ray_end - ray_begin

        #construct vec orthogonal to ray
        norm_vector = [ray_vec[1], -1*ray_vec[0], 0]

        #get center of ray keep z at center
        center = self.ds.domain_center

        ray_center = (ray_begin + ray_end)/2
        ray_center[2] = center[2]

        #handle case where it is an on axis slice in the y plane
        #yt will ignore north_vector and place z-axis on horizontal axis
        if (norm_vector[0] == 0):
            # change yt coordinates so that z-axis is vertical
            self.ds.coordinates.x_axis[1] = 0
            self.ds.coordinates.x_axis['y'] = 0

            self.ds.coordinates.y_axis[1] = 2
            self.ds.coordinates.y_axis['y'] = 2

        #Create slice along ray. keep slice pointed in z-direction
        slice = yt.SlicePlot(self.ds,
                          norm_vector,
                          self.slice_field,
                          north_vector = [0, 0, 1],
                          center = center,
                          width=(norm(ray_vec), "cm"),
                          axes_unit="kpc")


        # add ray to slice
        slice.annotate_ray(self.ray, arrow=True)

        # set y label to Z
        slice.set_ylabel("Z (kpc)")

        # set color map
        slice.set_cmap(field=self.slice_field, cmap = cmap)

        #assign slice
        self.slice = slice
        return slice

    def plot_spect(self, ax, fname=".temp.h5"):
        """
        Use trident to plot the absorption spectrum of the ray.
        currently defaults to using COS wavelength binning and range.
        Parameters:
            ax : a matplotlib axis in which to draw the plot
            fname=".temp.h5" : filename to temporary record spectrum data
        Returns:
            none
        """
        #set max and min wavelength and resolution
        wave_min = self.wavelenth_center - self.wavelength_width/2
        wave_max = self.wavelenth_center + self.wavelength_width/2
        #generate spectrum defined by inputs
        spect_gen = trident.SpectrumGenerator(lambda_min=wave_min, lambda_max=wave_max, dlambda = self.resolution)
        spect_gen.make_spectrum(self.ray, lines=self.ion_list, output_file = fname)

        #save the spectrum to hdf5
        spect_gen.save_spectrum(fname)

        #open saved file to get wavelength and flux
        spect_h5 = h5py.File(fname)

        wavelength = np.array(spect_h5['wavelength'])
        flux = np.array(spect_h5['flux'])

        #close and delete temp file
        spect_h5.close()
        remove(fname)

        #plot values
        ax.plot(wavelength, flux)
        ax.set_ylim(0, 1.05)
        ax.set_title("Spectrum".format(self.ion_name))
        ax.set_xlabel("Wavelength (Angstrom)")
        ax.set_ylabel("Flux")

    def plot_num_density(self, ax):
        """
        Plots the number density at different lengths along the ray

        Parameters:
            ax : a matplotlib axis in which to draw the plot
        Returns:
            none
        """
        #get list of num density and corresponding lengths
        num_density = np.array(self.ray_h5['grid'][self.ion_p_name()+'_number_density'])
        dl_list = np.array(self.ray_h5['grid']['dl'])

        #convert list of dl's to list of lengths from begin of ray
        num_dls = len(dl_list)
        for i in range(1, num_dls):
            dl_list[i] += dl_list[i-1]

        # convert to kpc
        dl_list = dl_list*self.cm_to_kpc

        #make num density plots
        ax.plot(dl_list, num_density)
        ax.set_title("Number Density of {} Along Ray".format(self.ion_name))
        ax.set_xlabel("Length From Start of Ray $(kpc)$")
        ax.set_ylabel("Number Density $(cm^{-3})$")
        ax.set_yscale('log')

        #chech if min/max num dense was called
        if (self.num_dense_min == 'None' and self.num_dense_max == 'None'):
            pass
        else:
            ax.set_ylim(self.num_dense_min, self.num_dense_max)

    def create_multiplot(self, outfname='None', cmap="magma"):
        """
        combines the slice plot, number density plot, and spectrum plot into
        one image.

        Parameters:
            outfname='None' : the file name/path in which to save the file defaults
                              to being unsaved

            cmap='magma' :     the color map to use for the slice plot

        Returns:
            none
        """
        if (self.slice == 'None'):
            #create the slicePlot using the field of the ion density
            self.create_slice(cmap = cmap)

        grid = AxesGrid(self.fig, (0.075,0.075,0.85,0.85),
                        nrows_ncols = (1, 1),
                        axes_pad = 1.0,
                        label_mode = "L",
                        share_all = True,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_size="3%",
                        cbar_pad="0%")


        plot = self.slice.plots[self.slice_field]
        plot.figure = self.fig
        plot.axes = grid[0].axes
        plot.cax = grid.cbar_axes[0]

        self.slice._setup_plots()

        ax2 = self.fig.add_subplot(312)
        self.plot_num_density(ax2)

        ax3 = self.fig.add_subplot(313)
        self.plot_spect(ax3)

        ax2.set_position([1.1, 0.52, 1, 0.42])
        ax3.set_position([1.1, 0, 1, 0.42])

        if (outfname != 'None'):
            self.fig.savefig(outfname, bbox_inches='tight')

    def zoom(self, factor):
        """
        Zoom into the slice by specified factor

        Parameters:
            factor : factor by which to zoom in using yt's zoom mehtod

        Returns:
            none
        """

        self.slice.zoom(factor)

    def close(self):
        """
        close all opened files
        """

        self.ds.close()
        self.ray_h5.close()

    def compute_col_density(self):
        """
        computes the column density along the given ray for a given ion species. This is done
        by summing the product of the number density for a given length by that length.

        Parameters:

        """
        if (self.ray_h5 == 'None'):
            self.ds, self.ray, self.ray_h5 = self.open_files()

        #get list of num density and corresponding length
        num_density = np.array(self.ray_h5['grid'][self.ion_p_name()+'_number_density'])
        dl_array = np.array(self.ray_h5['grid']['dl'])

        #multiply num density by its dl and sum up to get column density
        col_density = np.sum( num_density*dl_array )

        return col_density

class movie_multi_plot(multi_plot):

    def __init__(self,
            ds_filename,
            ray_dir,
            ion_name='H I',
            slice_field=None,
            absorber_fields=[],
            wavelength_center=None,
            wavelength_width = 300,
            resolution = 0.1,
            out_dir="./images"):

        self.ds_filename = ds_filename

        #create directory if doesn't exist
        try:
            makedirs(out_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        #collect only ray files

        #find all ray files in ray_dir. set first one to ray_filename
        self.ray_files = sorted(listdir(ray_dir))
        self.ray_dir = ray_dir

        self.ion_name = ion_name
        #add ion name to list of all ions to be plotted
        absorber_fields.append(ion_name)
        self.ion_list = absorber_fields

        #set a value for slice
        self.slice =None

        #set slice field to ion name if no field is specified
        if (slice_field == None):
            self.slice_field = self.ion_p_name() + "_number_density"
        else:
            self.slice_field = slice_field

        self.wavelength_width = wavelength_width
        self.resolution = resolution
        #default set the wavelenth center to one of the known spectral lines
        #for ion name. Use tridents line database to search for correct wavelenth
        if (wavelength_center == None):
            #open up tridents default line database
            lbd = trident.LineDatabase("lines.txt")
            #find all lines that match ion
            lines = lbd.parse_subset(subsets= [self.ion_name])
            #take first one and return its wavelenth
            self.wavelenth_center = lines[0].wavelength

        #calculate the proper yscale for number density plot
        tot_median=0.
        #for rfile in self.ray_files:
        #    #open ray file and get number density
        #    ray_h5 = h5py.File("{}/{}".format(self.ray_dir, rfile))
        #    num_density = np.array(ray_h5['grid'][self.ion_p_name()+'_number_density'])
        #
        #    #add median num_density to sum
        #    tot_median += np.median(num_density)
        #
        #    #close ray file
        #    ray_h5.close()

        #get middle ray to represent scale
        middle_ray = self.ray_files[ int(len(self.ray_files)/2) ]
        ray_h5 = h5py.File("{}/{}".format(self.ray_dir, middle_ray))

        #get median num density
        num_density = np.array(ray_h5['grid'][self.ion_p_name()+'_number_density'])
        med = np.median(num_density)

        #estimate min max values to plot
        self.num_dense_min = 0.01*med
        self.num_dense_max = 1000*med

        self.out_dir = out_dir

    def create_movie(self, cmap="magma"):
        """
        creates a movie by combining all the plots made from the ray in ray_dir

        Parameters:
            cmap="magma" : the colormap with which to use for the slice plot
        """
        #open up dataset
        self.ds = yt.load(self.ds_filename)
        #create fig to plot on
        print("loaded yt. tryinn fig??")
        self.fig = plt.figure(figsize=(10, 10))
        print("well fig worked")

        num_rays = len(self.ray_files)
        for i in range(num_rays):
            #check that its an .h5 file
            if (self.ray_files[i][-3:] != ".h5"):
                continue

            #assign the current ray filename
            self.ray_filename = "{}/{}".format(self.ray_dir,self.ray_files[i])
            #open the current ray file
            junk_var, self.ray, self.ray_h5 = self.open_files(open_ds=False)

            if self.slice != None:
                #annotate slice
                self.slice.annotate_clear()
                self.slice.annotate_ray(self.ray, arrow=True)

            #create multiplot using slice and current ray plots
            self.create_multiplot(outfname = "{}/mp{:04d}".format(self.out_dir, i), cmap=cmap)

            #close ray files and clear figure
            self.ray_h5.close()

            self.fig.clear()

def construct_rays( dataset,
                    line_list,
                    length=200,
                    n_rays=100,
                    direction='z',
                    angle=0,
                    dist_from_center=200,
                    out_dir='./rays',
                    DEBUG=False):
    """
    Constructs a number of light rays to "scan" a galactic data set using trident

    Parameters:
        dataset: enzo dataset on which to construct the rays
        line_list: list of ions to include in rays
        length: the length of the rays in kpc
        n_rays: the number of rays to construct
        direction: the coordinate direction in which to "scan" the galaxy
        angle: The azimuthal angle around the direction. in degrees
        dist_from_center: range to construct rays. in kpc from center of galaxy
        out_dir: directory in which to save the rays

    Returns:
        none
    """
    if DEBUG:
        print(trident.__version__)
    ds = yt.load(dataset)
    #add ion fields to dataset if not already there
    trident.add_ion_fields(ds, ions=line_list, ftype='gas')

    #create directory if doesn't exist
    try:
        makedirs(out_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    #convert lengths to code_length and angle to radians
    length = ds.quan(length, 'kpc').in_units('code_length')
    dist_from_center = ds.quan(dist_from_center, 'kpc').in_units('code_length')
    angle = np.deg2rad(angle)

    #get right indexing for direction
    if (direction == 'x'):
        dir_index=0
        coord1_index=1
        coord2_index=2
    elif (direction == 'y'):
        dir_index=1
        coord1_index=2
        coord2_index=0
    elif (direction == 'z'):
        dir_index=2
        coord1_index=0
        coord2_index=1
    else:
        raise RuntimeError("direction must be 'x', 'y', or 'z'")

    #calculate the changing coordinate
    begin = ds.domain_center[dir_index] + dist_from_center
    end = ds.domain_center[dir_index] - dist_from_center
    direct_coord = np.linspace(begin, end, n_rays)

    #compute ray length in the coordinate directions
    len_coord1 = length* np.cos(angle)
    len_coord2 = length* np.sin(angle)

    #calc beginning and ending of ray for constant coordinates
    coord1_begin = ds.domain_center[coord1_index] - len_coord1
    coord1_end = ds.domain_center[coord1_index] + len_coord1

    coord2_begin = ds.domain_center[coord2_index] - len_coord2
    coord2_end = ds.domain_center[coord2_index] + len_coord2

    #empty ray coordinates to be filled
    ray_begin = ds.arr(np.zeros(3), "code_length")
    ray_end = ds.arr(np.zeros(3), "code_length")
    for i in range(n_rays):
        #set beginning ray
        ray_begin[dir_index] = direct_coord[i]
        ray_begin[coord1_index] = coord1_begin
        ray_begin[coord2_index] = coord2_begin

        #set ending ray
        ray_end[dir_index] = direct_coord[i]
        ray_end[coord1_index] = coord1_end
        ray_end[coord2_index] = coord2_end

        #construct ray
        trident.make_simple_ray(ds,
                                ray_begin,
                                ray_end,
                                lines=line_list,
                                data_filename="{:s}/ray{:04d}.h5".format(out_dir, i))

if __name__ == '__main__':
    data_set_fname = argv[1]
    ray_fname = argv[2]
    ion = 'C IV'
    absorbers = ['H I', 'O VI']

    mp = multi_plot(data_set_fname, ray_fname, ion_name=ion, absorber_fields=absorbers, wavelength_width = 100)

    outfile = "multiplot_" + ion[0] +".png"
    mp.create_multiplot(outfname=outfile)
