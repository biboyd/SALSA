import yt
import trident as tri
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

    def __init__(self, ds_filename, ray_filename, figure='None', ion='H I', open=True):
        """
        init file names and ion name

        Parameters:
        ds_filename : Path/name of the enzo dataset to be loaded
        ray_filename : Path/name of the hdf5 ray file to be loaded
        ion : Name of the ion in notaion
              Element symbol *space* roman numeral of ion level (i.e. "H I", "O VI")
        """
        self.ds_filename = ds_filename
        self.ray_filename = ray_filename
        self.ion_name = ion

        self.slice = 'None'
        if (figure == 'None'):
            self.fig = plt.figure(figsize=(10, 10))
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
        convert ion species name from trident style to one that
        can be used with h5 files

        Returns:
        outname : Name of the ion in the form found in hdf5 files
        """

        ######### Deprecated, no longer needed #########3
        # 'H I' is an exception just return H

        if self.ion_name == 'H I':
            return 'H'

        #split up the words in ion species name
        ion_split = self.ion_name.split()

        #convert num from roman numeral. subtract run b/c h5
        num = tri.from_roman(ion_split[1])-1

        #combine all the names
        outname = ion_split[0] + '_p' + str(num)
        return outname

    def create_slice(self, field='density', cmap="BLUE"):
        """
        Create a slice in the Dataset along the path of the ray.
        Choose to keep the Z direction maintained.

        Parameters:
        field='density' : The yt field to plot for the slice
        cmap='BLUE' : the colormap to use for the slice

        Returns:
        slice : yt SlicePlot with ray annotated
        """

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
                          field,
                          north_vector = [0, 0, 1],
                          center = center,
                          width=(norm(ray_vec), "cm"),
                          axes_unit="kpc")
        # add ray to slice
        slice.annotate_ray(self.ray)

        # set y label to Z
        slice.set_ylabel("Z (kpc)")

        # set color map
        slice.set_cmap(field=field, cmap = cmap)

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
        line_list = [self.ion_name]
        #generate spectrum defined by inputs
        spect_gen = tri.SpectrumGenerator('COS')
        spect_gen.make_spectrum(self.ray, lines=line_list, output_file = fname)

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
        ax.set_title("Spectra {} ".format(self.ion_name))
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

    def create_multiplot(self, outfname='None', cmap="BLUE",field='density'):

        if (self.slice == 'None'):
            #create the slicePlot
            self.create_slice(field = field, cmap = cmap)

        grid = AxesGrid(self.fig, (0.075,0.075,0.85,0.85),
                        nrows_ncols = (1, 1),
                        axes_pad = 1.0,
                        label_mode = "L",
                        share_all = True,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_size="3%",
                        cbar_pad="0%")


        plot = self.slice.plots[field]
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
