import matplotlib as mpl
mpl.use('Agg')
import yt
import trident
import numpy as np

from sys import argv, path
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import AxesGrid
from astropy.units import angstrom

from salsa.utils.functions import ion_p_num
from salsa.utils.defaults import default_units_dict, default_limits_dict
from salsa.absorber_extractor import SPICEAbsorberExtractor

class AbsorberPlotter():
    """
    Create's plot to easily see where absorbers are found along the light ray
    both directly and through resulting spectra. Uses absorber extractor as base
    class.

    Parameters
    ----------
    abs_ext : AbsorberExtractor object
        Absorber extractor used to extract the data we want to plot.
        Can be any of the child classes; e.g. SPICEAbsorberExtractor

    ion_name: string, optional
        Name of the ion to extract absorbers of and plot
        Default: "H I"

    ftype : string, optional
        Field type for ion fields. Allows ion fields to be added
        if not already present.
        Default: "gas"

    absorber_fields : list of strings, optional
        Additional ions to include in plots/Spectra, enter as list.

    markers :bool, optional
        whether to include markers on light ray and number density plot

    mark_plot_args : dict, optional
        set the property of markers if they are to be plotted.
        optional settings are:

        `marker_spacing`: determines how far apart markers are in kpc.

        `marker_shape`: shape of marker see matplotlib for notation

        `marker_cmap`: colormap used to differentiate markers any other property that can be passer to matplotlib scatter

    figure :matplotlib figure, optional
        where the multi_plot will be plotted. creates one if none is specified.


    Notes
    ------
    ion names should be in following notaion: neutral Hydrogen-->"H I",
    5-times ionized Oxygen --> "O VI"
    """

    def __init__(self,
                 abs_ext,
                 ion_name='H I',
                 ftype='gas',
                 absorber_fields=[],
                 markers=True,
                 mark_plot_args=None,
                 figure=None):

        self.abs_ext = abs_ext
        self.ds = abs_ext.ds

        self.ion_name = ion_name
        #add ion name to list of all ions to be plotted
        self.ion_list = [ion_name] + absorber_fields

        #add ion fields to dataset if not already there
        trident.add_ion_fields(self.ds, ions=self.ion_list, ftype=ftype)

        #set a value for slice
        self.slice = None

        #open up a figure if none specified
        if figure is None:
            self.fig = plt.figure(figsize=(10, 10))
        else:
            self.fig = figure

        #set marker plot properties
        self.markers = markers
        if markers:
            self.mark_kwargs = {'alpha' : 0.45,
                                's' : 100,
                                'edgecolors' : 'black',
                                'linewidth' : 3,
                                'spacing' :50,
                                'marker_cmap' : 'viridis',
                                'marker_shape' :'s'}
            if mark_plot_args != None:
                self.mark_kwargs.update(mark_plot_args)

            self.marker_spacing = self.mark_kwargs.pop('spacing')
            self.marker_cmap = self.mark_kwargs.pop('marker_cmap')
            self.marker_shape = self.mark_kwargs.pop('marker_shape')
        self.mark_dist_arr = None

    def _add_annotations(self, plot_ray=True):
        """
        Adds ray and marker annotations to slice plot

        Parameters
        -----------
        plot_ray : bool
            Whether to annotate ray & markers to the slice plot or to just
            calculate the marker positions for placing on los_velocity plot
        """

        if plot_ray:
            #annotate ray
            self.slice.annotate_ray(self.abs_ext.ray,
                                    arrow=True,
                                    alpha=0.5,
                                    color='white',
                                    linewidth=2)

        if self.markers:
            #get ray positional properties
            ray_begin, ray_end, ray_length, ray_direction = self.abs_ext.ray_position_prop(units='kpc')
            #make marker every x kpc. skip start
            mark_dist = self.marker_spacing #kpc
            mark_dist_arr = np.arange(mark_dist, ray_length.value, mark_dist)
            self.mark_dist_arr = self.ds.arr(mark_dist_arr, 'kpc')

            #define colormap and scale
            mrk_cmap = plt.get_cmap(self.marker_cmap)
            self.colorscale = np.linspace(0, 1, mark_dist_arr.size)

            if plot_ray:
                #construct unit vec from ray
                for i in range(mark_dist_arr.size):
                    #calculate the position
                    mrk_pos = ray_begin + ray_direction * self.mark_dist_arr[i]

                    #choose correct color from cmap
                    mrk_kwargs = self.mark_kwargs.copy()
                    mrk_kwargs['color'] = mrk_cmap(self.colorscale[i])

                    self.slice.annotate_marker(mrk_pos,
                                               marker=self.marker_shape,
                                               **mrk_kwargs)

    def create_slice(self,
                     center = None,
                     slice_field = None,
                     north_vector = [0, 0, 1],
                     plot_ray = True,
                     cmap="magma",
                     height=None,
                     width=None):
        """
        Create a slice in the Dataset along the path of the ray.
        Choose to keep the Z direction maintained.

        Parameters
        ----------

        center : array type, optional
            Center of the slice in code_length.
            If None, then defaults to domain_center

        slice_field : string, optional
            Field to plot. Defaults to the number density for the ion
            associated with the abs_ext attribute.

        north_vector : array type, optional
            vector used to fix the orientation of the slice plot.
            Defaults to z-axis

        plot_ray : bool
            Whether to annotate ray/markers to the slice plot or to just
            calculate the marker positions for placing on los_velocity plot

        cmap: str
            the colormap to use for the slice. Default: 'magma'

        height, width: float, optional
            dimension for the figure

        Returns
        -------
        slice : yt SlicePlot
            Slice with ray annotated. This SlicePlot is also saved to the `.slice` attribute.
        """
        #set slice field to ion name if no field is specified
        if (slice_field is None):
            slice_field = ion_p_num(self.ion_name)
        else:
            slice_field = slice_field

        ray_begin, ray_end, ray_length, ray_unit = self.abs_ext.ray_position_prop(units='kpc')

        #construct vec orthogonal to ray/plane
        north_vector = self.ds.arr(north_vector, 'dimensionless')
        norm_vector = np.cross(ray_unit, north_vector)
        norm_vector = norm_vector/np.linalg.norm(norm_vector)

        #set center to domain_center unless specified
        if center is None:
            center = self.ds.domain_center
        else:
            center = self.ds.arr(center, 'code_length')

        #adjust center so that it is in the plane of ray and north_vector
        ray_center = (ray_begin + ray_end)/2
        ray_center = ray_center.in_units('code_length')
        center_dif = ray_center - center
        scale = self.ds.quan(np.dot(center_dif, norm_vector), 'code_length')
        center = scale*norm_vector + center

        #set width/height
        if width ==None and height==None:
            #default to length of ray
            wid_hght=ray_length
        elif width ==None and height != None:
            #width still defaults to length of ray
            wid_hght=((ray_length.to_value(), 'kpc'), (height, 'kpc'))
        elif width != None and height ==None:
            wid_hght= ((width, 'kpc'), (ray_length.to_value(), 'kpc'))
        else:
            wid_hght= ((width, 'kpc'), (height, 'kpc'))

        #Create slice along ray. keep slice pointed in north_Vec direction
        self.slice = yt.OffAxisSlicePlot(self.ds,
                          norm_vector,
                          slice_field,
                          center = center,
                          north_vector = north_vector,
                          width = wid_hght)

        #set axes to kpc
        self.slice.set_axes_unit('kpc')

        #annotate plot
        self._add_annotations(plot_ray=plot_ray)

        # set color map
        self.slice.set_cmap(field=slice_field, cmap = cmap)

        # set background to bottom of color map
        self.slice.set_background_color(slice_field)

        return self.slice

    def plot_vel_space(self,
                       velocity_width = 3000,
                       ax_vel=None,
                       **kwargs):
        """
        Use trident to plot the absorption spectrum of the ray in velocity
        space.
        Uses the wavelength centers and velocity resolution of the AbsorberExtractor
        to generate the spectrum.

        Parameters
        -----------
        velocity_width : float, optional
            sets the velocity range of spectrum plot in units of km/s

        ax_vel : matplotlib axis object, optional
            an axis in which to draw the velocity plot. If None, no plot is
            not drawn.
            Default: None

        Returns
        ----------
            velocity: YT array
                Array of velocity values of the generated spectra (in km/s)

            flux: YT array
                Array of the normalized flux values of the generated spectra

        """

        # add wav center first so it is set as zero point velocity by trident
        wav = int( np.round(self.abs_ext.wavelength_center) )
        line = f"{self.ion_name} {wav}"
        ion_list = [line] + self.ion_list

        #set up spectra
        vel_min = -velocity_width/2
        vel_max = velocity_width/2
        spect_gen = trident.SpectrumGenerator(lambda_min=vel_min,
                                              lambda_max=vel_max,
                                              dlambda = self.abs_ext.velocity_res,
                                              bin_space="velocity")

        #generate spectra and return fields
        spect_gen.make_spectrum(self.abs_ext.ray_data, lines=ion_list)
        flux = spect_gen.flux_field
        velocity = spect_gen.lambda_field.in_units('km/s')

        if ax_vel is not None:
            #plot values for velocity plot
            ax_vel.plot(velocity[:-1], flux[:-1])
            ax_vel.set_ylim(0, 1.05)
            ax_vel.set_xlim(vel_min, vel_max)
            ax_vel.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax_vel.set_title(f"Rel. to line {self.abs_ext.wavelength_center:.1f}"+r"$\AA$", loc='right')
            ax_vel.set_xlabel("Delta_v (km/s)")
            ax_vel.set_ylabel("Flux")
            ax_vel.grid(zorder=0, which='both')

        return velocity, flux

    def plot_lambda_space(self, 
                          wavelength_width = 30,
                          wavelength_res = 0.1,
                          ax=None,
                          **kwargs):
        """
        Use trident to plot the absorption spectrum of the ray. Plot in
        wavelegnth (lambda) space.

        Parameters
        -----------    
        wavelength_width :float, optional
            sets the wavelength range of the spectrum plot. defaults to 300 Angstroms

        wavelength_res :float, optional
            width of wavelength bins in spectrum plot. default 0.1 Angstrom

        ax : matplotlib axis
            axis in which to draw the spectra plot

        Returns
        --------
            wavelength: YT array
                Array of wavelength values of the generated spectra (in km/s)

            flux: YT array
                Array of the normalized flux values of the generated spectra
        """
        #set which ions to add to spectra
        ion_list = self.ion_list

        #adjust wavelegnth_center for redshift
        rest_wavelength = self.abs_ext.wavelength_center
        wave_min = rest_wavelength - wavelength_width/2
        wave_max = rest_wavelength + wavelength_width/2

        #use wavelength_width to set the range
        spect_gen = trident.SpectrumGenerator(lambda_min=wave_min, lambda_max=wave_max, dlambda = wavelength_res)
        spect_gen.make_spectrum(self.abs_ext.ray_data, lines=ion_list, observing_redshift=self.ds.current_redshift)


        #get fields from spectra and give correct units
        #also convert from astropy to unyt for consistency
        rest_wavelength = rest_wavelength * angstrom
        wavelength = spect_gen.lambda_field * angstrom
        flux = spect_gen.flux_field

        if ax is not None:
            #plot values for spectra
            ax.plot(wavelength[:-1], flux[:-1])
            ax.set_ylim(0, 1.05)
            ax.set_xlim(wave_min, wave_max)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.set_title(f"Spectrum {self.ion_name}", loc='right')
            ax.set_xlabel(r"Wavelength $\AA$")
            ax.set_ylabel("Flux")
            ax.grid(zorder=0, which='both')

        return wavelength, flux

    def plot_num_density(self, 
                         ax_num_dense,
                         num_dense_min=None,
                         num_dense_max=None,
                         ax_prop2=None,
                         prop2_name='velocity_los', 
                         prop2_units=None,
                         plot_spice_intervals=False,
                         plot_kwargs={},
                         **kwargs):
        """
        Plots the number density at different lengths along the ray

        Parameters
        ----------
        ax_num_dense : matplotlib axis
            axis in which to draw the number density plot

        num_dense_min: float, optional
            Sets the lower limit for the number density plot. If None, defaults to
            0.01 times the median number density

        num_dense_max: float, optional
            Sets the upper limit for the number density plot. If None, defaults to
            100 times the median number density

        ax_prop2: matplotlib axis, optional
            axis in which to draw the second field. if None, no plot is made
            Default: None

        prop2_name: str, optional
            Field to plot in second plot.
            Default:'velocity_los'

        prop2_units: str, optional
            The units to use in the second field.
            Defaults: None

        plot_spice_intervals: bool, optional
            If using a SPICEAbsorberExtractor, plot the SPICE intervals.
            Has no effect for any other AbsorberExtractor class.

        plot_kwargs: dict, optional
            A Dictionary of plot kwargs to be passed to the pyplot.plot function.
            Default: {}
        """
        #get list of num density  los velocity and corresponding lengths
        num_density = self.abs_ext.ray_data[ion_p_num(self.ion_name)]
        prop2 = self.abs_ext.ray_data[prop2_name]

        prop2_lb = None
        prop2_ub = None

        #load defaults for prop2
        if prop2_units is None:

            #check if defaults for field
            if prop2_name in default_units_dict.keys():
                prop2_units = default_units_dict[prop2_name]
                prop2 = prop2.in_units(prop2_units)

                # set default bounds
                if prop2_name in default_limits_dict.keys():
                    prop2_lb, prop2_ub = default_limits_dict[prop2_name]

            else:
                # use default of trident (genarlly in cgs units)
                prop2_units = str(prop2.units)
        else:
            # convert to specified units
            prop2 = prop2.in_units(prop2_units)

        #get length data and define x limits
        l_list = self.abs_ext.ray_data['l'].in_units('kpc')
        pad = 0.1*l_list[-1]
        xlimits = [-pad, l_list[-1] + pad]

        # check if l_list is non-empty cuz something went wrong then.
        if l_list.size == 0:
            err_file = open("error_file.txt", 'a')
            err_file.write(f"{self.abs_ext.ray_filename} had an l_list that was of size zero")
            err_file.close()
            return 1

        if ax_num_dense is not None:
            #make num density plots
            ax_num_dense.plot(l_list, num_density, **plot_kwargs)
            ax_num_dense.set_title(f"Number Density of {self.ion_name} Along Ray", loc='right')
            ax_num_dense.set_xlabel("Length From Start of Ray $(kpc)$")
            ax_num_dense.set_ylabel("Number Density \n$(cm^{-3})$")
            ax_num_dense.set_yscale('log')
            ax_num_dense.grid(zorder=0)

            #chech if min/max num dense was called
            if (num_dense_min is None and num_dense_max is None):
                med = np.median(num_density)
                num_dense_min = med*0.01
                num_dense_max = med*1000

            #set axes limits
            ax_num_dense.set_ylim(num_dense_min, num_dense_max)
            ax_num_dense.set_xlim(xlimits[0], xlimits[1])
            #add minor tick marks
            ax_num_dense.xaxis.set_minor_locator(AutoMinorLocator(2))

            #check if should plot spice intervals
            if plot_spice_intervals and isinstance(self.abs_ext, SPICEAbsorberExtractor):
                
                if self.abs_ext.num_feat is None:
                    raise RuntimeError("Please run the AbsorberExtractor's get_current_absorbers() method before plotting.")

                #plot spice intervals
                if self.abs_ext.num_feat > 0:
                    data = self.abs_ext.df['col_dens']
                    for i in range(self.abs_ext.num_feat):
                        b, e = self.abs_ext.features[i]

                        #plot interval
                        ax_num_dense.axvspan(l_list[b], l_list[e], alpha=0.5, edgecolor='black',facecolor='tab:grey')

                        #plot on 2nd prop if axis exists
                        if ax_prop2 is not None:
                            ax_prop2.axvspan(l_list[b], l_list[e], alpha=0.5, edgecolor='black',facecolor='tab:grey')

                    #plot number of intervals found
                    box_props = dict(boxstyle='square', facecolor='white')
                    ax_num_dense.text(0.9, 0.85, f"{self.abs_ext.num_feat} feat.", transform=ax_num_dense.transAxes, bbox = box_props)

                    #take three largest absorbers and sort by position
                    max_indices = data.argsort()
                    max_indices = max_indices[-3:]
                    max_indices.sort()

                    #plot markers from left to right
                    colors=['black', 'magenta', 'yellow']
                    for i,c in zip(max_indices, colors):
                        b, e = self.abs_ext.features[i]
                        lcd = np.log10(data[i].value)
                        mid_point = (l_list[b]+l_list[e])/2

                        ax_num_dense.scatter(mid_point, 0.75*num_dense_max,
                                             marker='v',color=c, edgecolors='black',
                                             label=f"logN={lcd:.1f}", zorder=3)

                    ax_num_dense.legend()

        #make second plot
        if ax_prop2 is not None:
            #plot zero mark if line of sight velocity
            if ax_prop2 == 'los_velocity':
                ax_prop2.hlines(0, l_list[0], l_list[-1], linestyles='dashed',alpha=0.25, zorder=1)

            ax_prop2.plot(l_list, prop2, **plot_kwargs)

            #set title and axis labels
            ax_prop2.set_title(f"{prop2_name} Along Ray", loc='right')
            ax_prop2.set_xlabel("Length From Start of Ray $(kpc)$")
            ax_prop2.set_ylabel(f"{prop2_name} $({prop2_units})$")

            #set limits
            ax_prop2.set_ylim(prop2_lb, prop2_ub)
            ax_prop2.set_xlim(xlimits[0], xlimits[1])

            #set minor ticks and grid lines
            ax_prop2.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax_prop2.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax_prop2.grid(zorder=0, which='major')
            ax_prop2.grid(zorder=0, which='minor', axis='y')

        #add appropriate markers to the plot
        if self.markers and ax_prop2 is not None:
            #check if marker distances have been defined
            if self.mark_dist_arr is None:
                self._add_annotations(plot=False)

            Vys = np.zeros_like(self.mark_dist_arr) - 500
            plot_markers = {}
            plot_markers.update(self.mark_kwargs)
            ax_prop2.scatter(self.mark_dist_arr.value, Vys,zorder=3, c=self.colorscale, marker=self.marker_shape, cmap=self.marker_cmap, **plot_markers)

    def plot_multiplot(self,
                       outfname=None,
                       center = None,
                       slice_field = None,
                       north_vector = [0, 0, 1],
                       cmap="magma",
                       make_spectra = True,
                       plot_spice_intervals=True,
                       annotate_column_density=True, 
                       velocity_width = 3000,
                       **kwargs):
        """
        Combines the slice plot, number density plot, and (optionally) spectrum plot into
        one image and saves it to disk.

        Can supply any optional keyword arguments supported by plot_num_density
        and plot_vel_space.

        Parameters
        -----------
        outfname: str , optional
            the file name/path in which to save the file defaults to being unsaved.
            Default: None

        gal_center : array type, optional
            Center of the slice in code_length.
            If None, then defaults to domain_center

        slice_field : string, optional
            Field to plot. Defaults to the number density for the ion
            associated with the abs_ext attribute.

        north_vector : array type, optional
            vector used to fix the orientation of the slice plot.
            Defaults to z-axis

        cmap:  Colormap, optional
            the color map to use for the slice plot.
            Default: magma

        plot_spice_intervals: bool, optional
            whether or not to shade the absorber intervals found by SPICE in the
            number density plot

        make_spectrum: bool, optional
            whether to include the spectrum plot or not. Given that the transition
            from yt 3 to yt 4 has affected Trident's spectra generation
            (https://github.com/trident-project/trident/issues/217), you may
            wish to not include this plot.
            Default: True

        annotate_column_density: bool, optional
            whether or not to annotate the column densities of absorbers in the
            spectrum plot

        velocity_width: float, optional
            velocity range of the spectrum plot in km/s

        Returns
        ---------
        fig : matplotlib figure:
            figure multi_plot is drawn on
        axes : matplotlib axes
            axes the three lower plots are drawn on
        """
        if (slice_field is None):
            slice_field = ion_p_num(self.ion_name)
        else:
            slice_field = slice_field

        if (self.slice == None):
            #create the slicePlot using the field of the ion density
            self.create_slice(center = center,
                              slice_field = slice_field,
                              north_vector = north_vector,
                              cmap = cmap)

        grid = AxesGrid(self.fig, (0.,0.,0.5,0.5),
                        nrows_ncols = (1, 1),
                        axes_pad = 0.5,
                        label_mode = "L",
                        share_all = False,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_size="3%",
                        cbar_pad="0%")

        #redraw slice plot onto figure
        plot = self.slice.plots[slice_field]
        plot.figure = self.fig
        plot.axes = grid[0].axes
        plot.cax = grid.cbar_axes[0]

        self.slice._setup_plots()

        #set up axes and draw other plots to them
        ax1 = self.fig.add_subplot(411)
        ax2 = self.fig.add_subplot(412)
        if make_spectra:
            ax3 = self.fig.add_subplot(413)

        self.plot_num_density(ax1,
                              ax_prop2=ax2,
                              plot_spice_intervals=plot_spice_intervals,
                              **kwargs)
        if make_spectra:
            self.plot_vel_space(ax_vel=ax3,
                                velocity_width=velocity_width,
                                annotate_column_density=annotate_column_density,
                                **kwargs)

        axes= [ax1, ax2]
        if make_spectra:
            axes.append(ax3)

        #setup positioning for the plots underneath
        strt_pos = -0.255
        ax1.set_position( [0.0, strt_pos, 0.5, 0.15] )
        ax2.set_position( [0.0, strt_pos-0.16, 0.5, 0.15] )
        if make_spectra:
            ax3.set_position( [0.0, strt_pos-0.4, 0.5, 0.15] )

        #set num dense and los vel to share axis
        ax1.set_xlabel("")
        ax2.set_title("", loc='right')
        ax1.sharex(ax2)
        ax1.set_xticklabels([])
        if (outfname != None):
            self.fig.savefig(outfname, bbox_inches='tight')

        return self.fig, axes

    def close(self):
        """
        close all opened files
        """

        self.ds.close()
        self.abs_ext.ray.close()
        plt.close(self.fig)
