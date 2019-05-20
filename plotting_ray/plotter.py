import yt
import trident
import numpy as np
from sys import argv
import h5py
from os import remove
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


    def open_files(self):
        """
        Opens dataset and ray files int yt and h5py

        Parameters:
            none

        Returns:
        ds : dataset object loaded into yt
        ray : ray object loaded into yt
        ray_h5file : ray object loaded into h5py
        """
        ds = yt.load(self.ds_filename)
        ray = yt.load(self.ray_filename)
        ray_h5file = h5py.File(self.ray_filename)

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

    def create_slice(self, cmap="BLUE"):
        """
        Create a slice in the Dataset along the path of the ray.
        Choose to keep the Z direction maintained.

        Parameters:
        field: The yt field to plot for the slice
        cmap='BLUE' : the colormap to use for the slice

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

        #shift to set center at zero
        dl_list -= dl_list[-1]/2

        #make num density plots
        ax.plot(dl_list, num_density)
        ax.set_title("Number Density of {} Along Ray".format(self.ion_name))
        ax.set_xlabel("Length From Start of Ray $(kpc)$")
        ax.set_ylabel("Number Density $(cm^{-3})$")
        ax.set_yscale('log')

    def create_multiplot(self, outfname='None', cmap="BLUE"):
        """
        combines the slice plot, number density plot, and spectrum plot into
        one image.

        Parameters:
            outfname='None' : the file name/path in which to save the file defaults
                              to being unsaved

            cmap='BLUE' :     the color map to use for the slice plot

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

if __name__ == '__main__':
    data_set_fname = argv[1]
    ray_fname = argv[2]
    ion = 'H I'
    mp = multi_plot(data_set_fname, ray_fname)

    outfile = "class_combined_" + ion[0] +".png"
    mp.create_multiplot(outfile)
