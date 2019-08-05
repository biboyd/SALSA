import yt
import trident
import numpy as np
from spectacle.fitting import LineFinder1D
from sys import argv, path
from os import remove, listdir, makedirs
import errno
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy.linalg import norm
from center_finder import find_center
import astropy.units  as u

path.insert(0, "/mnt/home/boydbre1/Repo/absorber_generation_post")
from new_generate_contour_absorbers import identify_intervals_char_length

class multi_plot():
    """
    Plots three images side by side for easy analysis. The images are:
    A slice of dataset along lightray's path.
    The number density of the given ion along that path.
    The artificial absorption spectra that would be observed.
    """

    def __init__(self,
                ds_filename,
                ray_filename,
                ion_name='H I',
                slice_field=None,
                absorber_fields=[],
                north_vector=[0, 0, 1],
                center_gal = None,
                wavelength_center=None,
                wavelength_width = 300,
                resolution = 0.1,
                redshift = 0,
                bulk_velocity=None,
                use_spectacle=False,
                plot_spectacle=False,
                spectacle_defaults=None,
                contour = True,
                num_dense_min=None,
                num_dense_max=None,
                markers_nd_pos=None,
                markers=True,
                mark_plot_args=None,
                figure=None):
        """
        init file names and ion name

        Parameters:
        ds_filename : Path/name of the enzo dataset to be loaded
        ray_filename : Path/name of the hdf5 ray file to be loaded
        ion_name :string: Name of the ion to plot in number density plot
        slice_field :string: Field to plot in slice plot. defaults to ion_name's number density
        absorber_fields :list of strings: Additional ions to include in plots/Spectra, enter as list
        north_vector :array type: vector used to fix the orientation of the slice plot defaults to z-axis
        center_gal :array type: center of galaxy in code_length. if None, then defaults to domain_center
        wavelength_center :float: Wavelength to center spectrum plot on. defaults to the stringest
                        known spectral line of ion_name. in units of Angstrom
        wavelength_width :float: sets the wavelength range of the spectrum plot. defaults
                        to 300 Angstroms
        resolution :float: width of wavelength bins in spectrum plot. default 0.1 Angstrom
        redshift :float: redshift of galaxy's motion. adjusts velocity plot calculation.
        bulk_velocity : array type : bulk velocity of the galaxy in km/s
        use_spectacle : bool: Choose whether to use spectacle fit to compute col dense
        markers :bool: whether to include markers on light ray and number density plot
        mark_plot_args : dict : set the property of markers if they are to be plotted.
                        optional settings are:
                        marker_spacing : determines how far apart markers are in kpc
                        marker_shape : shape of marker see matplotlib for notation
                        marker_cmap : colormap used to differentiate markers
                        any other property that can be passer to matplotlib scatter
        ###NOTE### ion names should be in notaion:
              Element symbol *space* roman numeral of ion level (i.e. "H I", "O VI")
        figure :matplotlib figure: where the multi_plot will be plotted. creates one if
                        none is specified.
        """
        #set file names and ion name
        self.ds_filename = ds_filename
        self.ray_filename = ray_filename
        self.ion_name = ion_name

        #open up the dataset and ray files
        self.ds = yt.load(self.ds_filename)
        self.ray = yt.load(self.ray_filename)

        #add ion name to list of all ions to be plotted
        absorber_fields.append(ion_name)
        self.ion_list = absorber_fields

        #set a value for slice
        self.slice = None
        self.north_vector = north_vector
        self.center_gal= center_gal

        #set slice field to ion name if no field is specified
        if (slice_field == None):
            self.slice_field = self.ion_p_name() + "_number_density"
        else:
            self.slice_field = slice_field

        # set bulk velocity
        if bulk_velocity is None:
            self.bulk_velocity = 0
        else:
            ray_b, ray_e, ray_l, ray_u = self.ray_position_prop()
            self.bulk_velocity = np.dot(ray_u, bulk_velocity)
            self.bulk_velocity = self.ds.quan(self.bulk_velocity, 'km/s')

        self.contour = contour
        self.use_spectacle = use_spectacle
        self.plot_spectacle = plot_spectacle
        self.defaults_dict = {
            'bounds' :{
                'column_density' : (12.0, 23)
            },
            'fixed' : {
                'delta_lambda' : True,
                'column_density' : False
            }
        }
        #add user defined defaults
        if spectacle_defaults is not None:
            self.defaults_dict.update(spectacle_defaults)
        self.redshift = redshift
        self.wavelength_width = wavelength_width
        self.resolution = resolution
        #default set the wavelength center to one of the known spectral lines
        #for ion name. Use tridents line database to search for correct wavelength
        if wavelength_center is None:
            #open up tridents default line database
            lbd = trident.LineDatabase("my_lines.txt")
            #find all lines that match ion
            lines = lbd.parse_subset(subsets= [self.ion_name])
            #take one with largest f_value
            f_val = 0
            for line in lines:
                if line.f_value >= f_val:
                    f_val = line.f_value
                    self.wavelength_center = line.wavelength
        else:
            self.wavelength_center = wavelength_center

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


        #optionally set min/max value for number density plot
        self.num_dense_min = num_dense_min
        self.num_dense_max = num_dense_max

        #arrays storing spectra
        self.lambda_array = None
        self.velocity_array = None
        self.flux_array = None
        #optionally set position of markers on number density plot
        self.markers_nd_pos = markers_nd_pos

    def add_annotations(self):
        """
        Adds ray annotation and marker annotations to slice plot
        """

        #annotate ray
        self.slice.annotate_ray(self.ray, arrow=True, plot_args={'alpha':0.5, 'color':'white', 'linewidth':2})

        if self.markers:
            #get ray positional properties
            ray_begin, ray_end, ray_length, ray_direction = self.ray_position_prop(units='kpc')
            #make marker every x kpc. skip start
            mark_dist = self.marker_spacing #kpc
            mark_dist_arr = np.arange(mark_dist, ray_length.value, mark_dist)
            self.mark_dist_arr = self.ds.arr(mark_dist_arr, 'kpc')

            #define colormap and scale
            mrk_cmap = plt.cm.get_cmap(self.marker_cmap)
            self.colorscale = np.linspace(0, 1, mark_dist_arr.size)

            #construct unit vec from ray
            for i in range(mark_dist_arr.size):
                #calculate the position
                mrk_pos = ray_begin + ray_direction * self.mark_dist_arr[i]

                #choose correct color from cmap
                mrk_kwargs = self.mark_kwargs.copy()
                mrk_kwargs['color'] = mrk_cmap(self.colorscale[i])

                self.slice.annotate_marker(mrk_pos, marker=self.marker_shape, plot_args=mrk_kwargs)

    def ray_position_prop(self, units='code_length'):
        """
        returns positional/directional properties of the ray so that it can be used like a vector

        Parameters:
            units : YT defined units to return arrays in. defaults to 'code length'
        Returns:
            ray_begin : the starting coordinates of ray (YT arr)
            ray_end : the ending coordinates of the ray (YT arr)
            ray_length : the length of the ray (YT quan)
            ray_unit : unit vector showing direction of the ray (YT arr)
        """
        #get start and end points of ray. convert to defined units
        ray_begin = self.ray.light_ray_solution[0]['start']
        ray_end = self.ray.light_ray_solution[0]['end']

        ray_begin = ray_begin.in_units(units)
        ray_end = ray_end.in_units(units)

        #construct vector pointing in ray's direction
        ray_vec = ray_end - ray_begin
        ray_length = self.ds.quan(norm(ray_vec.value), units)

        #normalize vector to unit length
        ray_unit = ray_vec/ray_length

        return ray_begin, ray_end, ray_length, ray_unit

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

    def create_slice(self, cmap="magma", height=None, width=None):
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


        ray_begin, ray_end, ray_length, ray_unit = self.ray_position_prop(units='kpc')

        #construct vec orthogonal to ray/plane
        self.north_vector = self.ds.arr(self.north_vector, 'dimensionless')
        norm_vector = np.cross(ray_unit, self.north_vector)
        norm_vector = norm_vector/np.linalg.norm(norm_vector)


        #handle case where it is an on axis slice in the y plane
        #yt will ignore north_vector and place z-axis on horizontal axis
        """if norm_vector[0] == 0:
            # change yt coordinates so that z-axis is vertical
            self.ds.coordinates.x_axis[1] = 0
            self.ds.coordinates.x_axis['y'] = 0

            self.ds.coordinates.y_axis[1] = 2
            self.ds.coordinates.y_axis['y'] = 2"""

        #set center to domain_center unless specified
        if self.center_gal is None:
            self.center_gal = self.ds.domain_center
        else:
            self.center_gal = self.ds.arr(self.center_gal, 'code_length')

        #adjust center so that it is in the plane of ray and north_vector
        ray_center = (ray_begin + ray_end)/2
        ray_center = ray_center.in_units('code_length')
        center_dif = ray_center - self.center_gal
        scale = self.ds.quan(np.dot(center_dif, norm_vector), 'code_length')
        center = scale*norm_vector + self.center_gal

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
                          self.slice_field,
                          center=center,
                          north_vector = self.north_vector,
                          width = wid_hght)



        #set axes to kpc
        self.slice.set_axes_unit('kpc')

        #annotate plot
        self.add_annotations()

        # set color map
        self.slice.set_cmap(field=self.slice_field, cmap = cmap)

        # set background to bottom of color map
        self.slice.set_background_color(self.slice_field)

        return self.slice

    def plot_spect_vel(self, ax_spect=None, ax_vel=None, single_line=None):
        """
        Use trident to plot the absorption spectrum of the ray. Then
        convert wavelength to line of sight velocity and plot versus flux.
        Uses wavelength_center to calculate the velocity.

        Parameters:
            ax_spect : a matplotlib axis in which to draw the spectra plot
            ax_vel : a matplotlib axis in which to draw the velocity plot
            single_line : str : option to plot a single line ie "HI1216" for
                    fitting purposes.
        Returns:
            wavelength, velocity and flux arrays created by spectrum generator
                units are angstrom, km/s
        """
        if single_line is None:
            ion_list = self.ion_list
        else:
            ion_list = [single_line]
        # calc doppler redshift due to bulk motion
        c = yt.units.c
        beta = self.bulk_velocity/c
        z_dopp = (1 - beta)/np.sqrt(1 +beta**2) -1
        z_dopp = z_dopp.value
        #adjust wavelegnth_center for redshift
        rest_wavelength = self.wavelength_center*(1+self.redshift)*(1+z_dopp)
        #generate spectrum defined by inputs
        max_iter=5
        for _ in range(max_iter):
        #set wavelength limits for plots
            if single_line is None:
                wave_min = rest_wavelength - self.wavelength_width/2
                wave_max = rest_wavelength + self.wavelength_width/2
            else:
                fit_range = 150
                wave_min = rest_wavelength - fit_range
                wave_max = rest_wavelength + fit_range

            spect_gen = trident.SpectrumGenerator(lambda_min=wave_min, lambda_max=wave_max, dlambda = self.resolution, line_database='my_lines.txt')
            spect_gen.make_spectrum(self.ray, lines=ion_list)

            #increase width when doing a fit
            if single_line is None:
                break
            else:
                #check if wings are at flux of 1
                if spect_gen.flux_field[0] >= 0.99:
                    break
                #increase fit range
                else:
                    fit_range +=50

        #get necessary and give correct units
        rest_wavelength = rest_wavelength*u.Unit('angstrom')
        wavelength = spect_gen.lambda_field * u.Unit('angstrom')
        flux = spect_gen.flux_field

        wave_min = wave_min*u.Unit('angstrom')
        wave_max = wave_max*u.Unit('angstrom')


        #calc velocity using relativistic doppler equation
        doppler_equiv = u.equivalencies.doppler_relativistic(rest_wavelength)
        velocity = wavelength.to('km/s', equivalencies=doppler_equiv)
        vel_min = wave_min.to('km/s', equivalencies=doppler_equiv)
        vel_max = wave_max.to('km/s', equivalencies=doppler_equiv)

        #save arrays so they can be acessed
        if single_line is None:
            self.lambda_array = wavelength
            self.velocity_array = velocity
            self.flux_array = flux

        if ax_spect is not None:
            #plot values for spectra
            ax_spect.plot(wavelength[:-1], flux[:-1])
            ax_spect.set_ylim(0, 1.05)
            ax_spect.set_xlim(wave_min.value, wave_max.value)
            ax_spect.set_title(f"Spectrum {self.ion_name}", loc='right')
            ax_spect.set_xlabel("Wavelength $\AA$")
            ax_spect.set_ylabel("Flux")
            ax_spect.grid(zorder=0)
        if ax_vel is not None:
            #plot values for velocity plot
            ax_vel.plot(velocity[:-1], flux[:-1])
            ax_vel.set_ylim(0, 1.05)
            ax_vel.set_xlim(vel_min.value, vel_max.value)
            ax_vel.set_title(f"Rel. to line {self.wavelength_center:.1f} $\AA$", loc='right')
            ax_vel.set_xlabel("Delta_v (km/s)")
            ax_vel.set_ylabel("Flux")
            ax_vel.grid(zorder=0)


            #annotate plot with column densities
            sums_list, line_txt, line_models = self.compute_col_density()
            box_props = dict(boxstyle='square', facecolor='white')
            ax_vel.text(0.8, 0.05, line_txt, transform=ax_vel.transAxes, bbox = box_props)

            colors = ['tab:purple', 'tab:orange', 'tab:green']
            vel= np.linspace(vel_min.value, vel_max.value, 1000) *u.Unit('km/s')
            #plot individual column lines
            if line_models is not None:
                for mod, color in zip(line_models, colors):
                    #plot centroids of largest lines
                    dv = mod.lines[0].delta_v.value
                    cd = mod.lines[0].column_density.value
                    ax_vel.scatter(dv, 1, c=color, marker='v',zorder=1, label="logN={:04.1f}".format(cd))
                    #plott the largest lines
                    if self.plot_spectacle:
                        ax_vel.step(vel, mod(vel), linestyle='--', color=color, alpha=0.75)
                ax_vel.legend(loc='lower left')


        return wavelength, velocity, flux

    def plot_num_dense_los_vel(self, ax_num_dense=None, ax_los_velocity=None):
        """
        Plots the number density at different lengths along the ray

        Parameters:
            ax : a matplotlib axis in which to draw the plot
        Returns:
            none
        """
        #get list of num density  los velocity and corresponding lengths
        num_density = self.ray.data[self.ion_p_name()+'_number_density']
        los_vel = self.ray.data['velocity_los']
        los_vel = los_vel.in_units('km/s') + self.bulk_velocity
        dl_list = self.ray.data['dl']

        #convert dl's to array of lengths from begin of ray
        l_list = dl_list.copy()
        l_list = l_list.in_units('kpc')
        num_dls = l_list.size
        for i in range(1, num_dls):
            l_list[i] += l_list[i-1]

        if ax_num_dense is not None:
            #make num density plots
            ax_num_dense.plot(l_list, num_density)
            ax_num_dense.set_title(f"Number Density of {self.ion_name} Along Ray", loc='right')
            ax_num_dense.set_xlabel("Length From Start of Ray $(kpc)$")
            ax_num_dense.set_ylabel("Number Density \n$(cm^{-3})$")
            ax_num_dense.set_yscale('log')
            ax_num_dense.grid(zorder=0)

            #chech if min/max num dense was called
            if (self.num_dense_min is None and self.num_dense_max is None):
                med = np.median(num_density)
                self.num_dense_min = med*0.01
                self.num_dense_max = med*1000

            ax_num_dense.set_ylim(self.num_dense_min, self.num_dense_max)

            #check if should plot contour intervals
            if self.contour:
                intervals = self.get_contour_intervals()
                lcd_list =[]
                my_intervals=[]
                #check interval has high enough column density
                lim = self.defaults_dict['bounds']['column_density'][0]
                tot_lcd=0
                for i in range(len(intervals)):
                    b, e = intervals[i]
                    #compute log col density
                    curr_lcd = np.log10( np.sum(dl_list[b:e]*num_density[b:e]) )
                    if curr_lcd > lim:
                        lcd_list.append(curr_lcd)
                        my_intervals.append((b, e))
                        #plot interval
                        ax_num_dense.axvspan(l_list[b], l_list[e], alpha=0.5, color='tab:grey')
                        tot_lcd += 10**curr_lcd
                        #plot on los vel if axis exists
                        if ax_los_velocity is not None:
                            ax_los_velocity.axvspan(l_list[b], l_list[e], alpha=0.5, color='tab:grey')

                #add the total column density onto plot
                tot_lcd = np.log10(tot_lcd)
                box_props = dict(boxstyle='square', facecolor='white')
                ax_num_dense.text(0.8, 0.05, f"Total: {tot_lcd:.1f}", transform=ax_num_dense.transAxes, bbox = box_props)

                #take three largest absorbers
                max_indices = np.argsort(lcd_list)
                max_indices = max_indices[-3:]
                max_indices.sort()
                print(max_indices, lcd_list)

                #plot from left to right
                colors=['black', 'magenta', 'yellow']
                for i,c in zip(max_indices, colors):
                    b, e = my_intervals[i]
                    lcd = lcd_list[i]
                    mid_point = (l_list[b]+l_list[e])/2
                    ax_num_dense.scatter(mid_point, 0.75*self.num_dense_max,
                                         marker='v',color=c, edgecolors='black',
                                         label=f"logN={lcd:.1f}", zorder=3)

                ax_num_dense.legend(loc='lower left', bbox_to_anchor=(-0.015, 0.95))

        if ax_los_velocity is not None:
            #make line of sight velocity plots
            ax_los_velocity.hlines(0, l_list[0], l_list[-1], linestyles='dashed',alpha=0.25, zorder=1)
            ax_los_velocity.plot(l_list, los_vel)
            ax_los_velocity.set_title("LOS Velocity Along Ray", loc='right')
            ax_los_velocity.set_xlabel("Length From Start of Ray $(kpc)$")
            ax_los_velocity.set_ylabel("Line of Sight \nVelocity $(km/s)$")
            ax_los_velocity.grid(zorder=0)
            ax_los_velocity.set_ylim(-600, 600)

        #add appropriate markers to the plot
        if self.markers:
            ys = np.zeros_like(self.mark_dist_arr)
            if ax_num_dense is not None:
                if self.markers_nd_pos == None:
                    ys += 50*self.num_dense_min
                else:
                    ys += self.markers_nd_pos

                #ax_num_dense.scatter(self.mark_dist_arr.value, ys,zorder=1, c=self.colorscale, marker=self.marker_shape, cmap=self.marker_cmap, **self.mark_kwargs)
            if ax_los_velocity is not None:
                Vys = np.zeros_like(self.mark_dist_arr) - 500
                plot_markers = {}
                plot_markers.update(self.mark_kwargs)
                plot_markers.update({'alpha':1})
                ax_los_velocity.scatter(self.mark_dist_arr.value, Vys,zorder=3, c=self.colorscale, marker=self.marker_shape, cmap=self.marker_cmap, **plot_markers)

    def create_multi_plot(self, outfname=None, markers=True, cmap="magma"):
        """
        combines the slice plot, number density plot, and spectrum plot into
        one image.

        Parameters:
            outfname=None : the file name/path in which to save the file defaults
                              to being unsaved

            markers=True : boolean. adds markers to slice plot and number density
                            to aid analysis between those plots.

            cmap='magma' :     the color map to use for the slice plot

        Returns:
            fig : matplotlib figure: figure multi_plot is drawn on
            axes : matplotlib axes: axes the three lower plots are drawn on
        """
        if (self.slice == None):
            #create the slicePlot using the field of the ion density
            self.create_slice(cmap = cmap)

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
        plot = self.slice.plots[self.slice_field]
        plot.figure = self.fig
        plot.axes = grid[0].axes
        plot.cax = grid.cbar_axes[0]

        self.slice._setup_plots()

        #set up axes and draw other plots to them
        ax1 = self.fig.add_subplot(411)
        ax2 = self.fig.add_subplot(412)
        ax3 = self.fig.add_subplot(413)
        #ax4 = self.fig.add_subplot(414)
        self.plot_num_dense_los_vel(ax_num_dense=ax1, ax_los_velocity=ax2)
        self.plot_spect_vel(ax_vel=ax3)

        axes= [ax1, ax2, ax3]
        #setup positioning for the plots underneath
        strt_pos = -0.255
        ax1.set_position( [0.0, strt_pos, 0.5, 0.15] )
        ax2.set_position( [0.0, strt_pos-0.16, 0.5, 0.15] )
        ax3.set_position( [0.0, strt_pos-0.4, 0.5, 0.15] )
        #for i in range(len(axes)):
        #    axes[i].set_position( [0.0, strt_pos - i*0.225, 0.5, 0.15] )
        #set num dense and los vel to share axis
        ax1.set_xlabel("")
        ax2.set_title("", loc='right')
        ax1.get_shared_x_axes().join(ax1, ax2)
        ax1.set_xticklabels([])
        if (outfname != None):
            self.fig.savefig(outfname, bbox_inches='tight')

        return self.fig, axes


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
        self.ray.close()

    def compute_col_density(self):
        """
        computes the column density along the given ray for a given ion species.
        This is done by using spectacle if use_spectacle is True. as well as
        by summing the product of the number density for a given length by that length.

        Parameters:
            none

        Returns:
            sums_list : list: list of the sums by the summing line method and the
                        summing all along the ray method
            line_array : array : array of col density and delta_v for lines and
                        sum by summing lines as well as sum by summing number
                        density along ray. col dense in log(N) N in cm^2
                        delta_v in km/s
            line_text : string : a string of properly formatted col dense to
                        be added on to multi_plot
        """
        if self.ray is None:
            self.ray = yt.load(self.ray_filename)

        #compute col density from summing along ray
        num_density = np.array(self.ray.data[self.ion_p_name()+'_number_density'])
        dl_array = np.array(self.ray.data['dl'])

        #multiply num density by its dl and sum up to get column density of ray
        tot_ray_cd= np.sum( num_density*dl_array )
        tot_ray_cd= tot_ray_cd

        line_sum_cd=0
        line_models = None
        #check if should try and fit with spect.
        #don't try to fit lines with over 20 logN. Very oversaturated
        if self.use_spectacle:
            #create spectra for a single line to fit
            wav = int( np.round(self.wavelength_center) )
            line = f"{self.ion_name} {wav}"
            true_wav, true_vel, true_flux=self.plot_spect_vel(single_line=line)
            #format ion correctly to fit
            ion_wav= "".join(line.split())

            #constrain possible column density values
            #create line model
            line_finder = LineFinder1D(ions=[ion_wav], continuum=1, z=0,
                                       defaults=self.defaults_dict,fitter_args={'maxiter':2000},
                                       threshold=0.01, output='flux', auto_fit=True)

            #fit data
            try:
                fit_spec_mod = line_finder(true_vel, true_flux)
            except RuntimeError:
                print('fit failed(prolly hit max iterations)', self.ray)
                fit_spec_mod = None
            except IndexError:
                print('INDEX ERROR on', self.ray)
                fit_spec_mod = None
            #check that fit was succesful/at least one line
            if fit_spec_mod is None:
                print('line could not be fit on ray ', self.ray)
                tot_num_lines=0

            else:
                tot_num_lines = len(fit_spec_mod.lines)
                vel_array = np.linspace(-1500, 1500, 1000)*u.Unit('km/s')
                line_stats = fit_spec_mod.line_stats(vel_array)

                #compute total column density
                line_sum_cd = 0
                for cd in line_stats['col_dens']:
                    line_sum_cd+= 10**cd

                # get biggest lines (max of 3)
                num_lines = line_stats['col_dens'].size
                if num_lines > 3: num_lines = 3

                line_models = []
                indx_max = line_stats.argsort('col_dens')
                for indx in indx_max[-num_lines:]:
                    line = fit_spec_mod.lines[indx]
                    line_models.append( fit_spec_mod.with_line(line, reset=True))

                #sort lines based on delta v
                line_models.sort(key=lambda mod: mod.lines[0].delta_v.value)

        log_line_sum = np.log10(line_sum_cd)
        log_tot_ray = np.log10(tot_ray_cd)
        sums = [log_line_sum, log_tot_ray]

        if self.use_spectacle:
            if line_sum_cd != 0:
                #compute percent difference
                diff = np.abs(line_sum_cd - tot_ray_cd)/tot_ray_cd*100
                #create text for text box
                line_text = "num lines: {}\n".format(tot_num_lines)+\
                           "line sum:  {:04.1f}\n".format(log_line_sum)+\
                           "total sum: {:04.1f}\n".format(log_tot_ray)+\
                           "diff: {:10.2f}%".format(diff)
            else:
                line_text = "num lines: {}\n".format(tot_num_lines)+\
                           "line sum:  {: >4s}\n".format('--')+\
                           "total sum: {:04.1f}\n".format(log_tot_ray)+\
                           "diff: {: >13s}%".format('--')
        else:
            line_text = 'logN={:04.1f}'.format(log_tot_ray)
        return sums, line_text, line_models


    def get_contour_intervals(self, char_density_frac = 0.5):
        """
        get the intervals along the ray where an absorption feature is found.
        currently write hdf5 file and read it b/c quickest implementation. Likely
        want to smoothline this process.

        Parameters:
            char_density_frac : float < 1: Fraction used to define when to end
                an absorption feature. Relative to maximum
        Returns:
            intervals : list of tuples : list of intervals containing absorption
                feature.
        """
        #define intial region to check
        cutoffs = {'H I':1e-11, 'C IV':1e-14, 'O VI':1e-11, 'Si III':1e-11, 'Si II':1e-11}
        init_cutoff = cutoffs[self.ion_name]

        num_density = self.ray.data[self.ion_p_name()+'_number_density']
        intervals = identify_intervals_char_length(num_density, init_cutoff, char_density_frac)

        return np.array(intervals)


if __name__ == '__main__':
    data_set_fname = argv[1]
    ray_fname = argv[2]
    ion = argv[3]
    num=int(argv[4])
    absorbers = [ion] #['H I', 'O VI']
    center, nvec, rshift, bv = find_center(data_set_fname)
    mp = multi_plot(data_set_fname, ray_fname, ion_name=ion, absorber_fields=absorbers,
                    center_gal=center, north_vector=nvec, bulk_velocity=bv,
                    redshift=rshift, wavelength_width = 30)
    makedirs("multi_plot_images", exist_ok=True)
    outfile = f"multi_plot_images/multi_plot_{ion[0]}_{num:02d}.png"
    mp.create_multi_plot(outfname=outfile)
