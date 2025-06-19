import yt
import trident
import numpy as np
import pandas as pd
import warnings

try:
    from spectacle.fitting import LineFinder1D
except ImportError:
    warnings.warn("spectacle not installed.")

from numpy.linalg import norm
import astropy.units  as u

from yt.data_objects.static_output import \
    Dataset

from salsa.utils.functions import ion_p_num, requires_spectacle
from salsa.utils.defaults import default_cloud_dict
from salsa.utils.collect_files import get_ray_num


class AbsorberExtractor():
    """
    Base class for extracting absorbers from a Trident lightray for a given ion species. 

    Only setup is performed; actual extraction is not performed by this class.

    Parameters
    --------------

    ds_filename: str or YT dataset
        Either Path/name of the dataset to be loaded or the dataset itself

    ray_filename: str or Trident ray
        Path/name of the hdf5 ray file to be loaded or the ray already loaded

    ion_name: string, optional
        Name of the ion to extract absorbers of
        Default: "H I"

    wavelegnth_center: float, optional
        The specific absorption line to look at (in unit Angstrom). None
        defaults to strongest absorption line for specified ion
        (using trident's ion table).
        Default: None

    velocity_res: float, optional
        Set minimum resolution (in km/s)  that spectacle will 
        attempt to fit lines to.
        Default: 10

    absorber_min: float, optional
        Minimum Log Column Density that will be used to define an absorber.
        If None, defaults to either default for specific ion or 13
        Default: None
    """

    def __init__(self, ds_filename,
                ion_name='H I', # Not sure this need to be here?
                wavelength_center=None,
                velocity_res=10,
                absorber_min=None):

        #set dataset filename
        if isinstance(ds_filename, str):
            self.ds = yt.load(ds_filename)
        elif isinstance(ds_filename, Dataset):
            self.ds = ds_filename

        self.ion_name = ion_name
        self.ion_list = [ion_name]  # Add ion name to list of all ions to be plotted

        # These will be set by later methods
        self.ray_filename = None
        self.ray = None
        self.features = None
        self.num_feat = None
        self.df = None
        self.data = None

        if absorber_min is None:
            if self.ion_name in default_cloud_dict.keys():
                self.absorber_min = default_cloud_dict[self.ion_name]
            else:
                self.absorber_min=13
        else:
            self.absorber_min = absorber_min

        self.defaults_dict = {
            'bounds' :{
                'column_density' : (self.absorber_min-0.5, 23)
            },
            'fixed' : {
                'delta_lambda' : True,
                'column_density' : False
            }
        }

        self.velocity_res = velocity_res

        #default set the wavelength center to one of the known spectral lines
        #for ion name. Use tridents line database to search for correct wavelength
        if wavelength_center is None:
            #open up tridents default line database
            lbd = trident.LineDatabase('lines.txt')
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

    def load_ray(self, new_ray):
        """
        loads a new ray into the multi_plot class. (same dataset)

        Parameters
        -----------
        new_ray :str or yt.ray
            either filename to rayfile or a trident ray that's already opened

        """
        # reset absorber extraction variables
        self.features=None

        #store number of features found
        self.num_feat = None

        # to store absorber feature table
        self.df=None

        #check if str else assume is ray
        if isinstance(new_ray, str):
            self.ray_filename=new_ray
            self.ray = yt.load(new_ray)
        else:
            self.ray = new_ray
            self.ray_filename=new_ray.filename_template

        # store data
        self.data = self.ray.all_data()

        # Check if ray is empty due to cuts
        if self.data['l'].size == 0:
            print(f'light ray {self.ray} is empty')

    def ray_position_prop(self, units='code_length'):
        """
        returns positional/directional properties of the ray so that it can be used like a vector

        Parameters
        -----------

        units : str, optional
            YT defined units to return arrays in. defaults to 'code length'.
            Default: 'code_length'
        Returns
        -------
        ray_begin : yt.arr
            the starting coordinates of ray

        ray_end : yt.arr
            the ending coordinates of the ray

        ray_length : yt.arr
            the length of the ray

        ray_unit : yt.arr
            unit vector showing direction of the ray
        """
        #get start and end points of ray. convert to defined units
        #print(self.ray)
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

    def close(self):
        """
        close all opened files, dataset, ray
        """

        self.ds.close()
        self.ray.close()

    def get_current_absorbers(self, *args, **kwargs):
        """
        Stub to be implemented by child classes.
        """
        raise NotImplementedError(
            "Method 'get_current_absorbers' is only defined in the child classes "
            "SPICEAbsorberExtractor and SpectacleAbsorberExtractor.")

    def get_all_absorbers(self, ray_list, *args, **kwargs):
        """
        Create catalog of the given rays using absorber extractor

        Parameters
        ----------
        ray_list: list of str or trident.ray objects
            List of ray objects or list of trident rays whose absorbers will be
            extracted

        Returns
        -------
        full_df: pandas.DataFrame
            Catalog of absorber properties in a pandas dataframe.
        """
        df_list=[]

        for ray in ray_list:
            #load new ray and extract absorbers
            self.load_ray(ray)
            df = self.get_current_absorbers(*args, **kwargs)

            if df is not None:
                # add ray index
                ray_num = get_ray_num(ray)
                for i in range(self.num_feat):
                    df.loc[i,'lightray_index'] = ray_num
                df_list.append(df)

        if len(df_list) > 0:
            full_df = pd.concat(df_list, ignore_index=True)
        else:
            full_df = None
        return full_df



@requires_spectacle
class SpectacleAbsorberExtractor(AbsorberExtractor):
    """
    Uses Spectacle to extract absorbers from a Trident lightray for a given ion species.

    Parameters
    --------------

    ds_filename: str or YT dataset
        Either Path/name of the dataset to be loaded or the dataset itself

    ray_filename: str or Trident ray
        Path/name of the hdf5 ray file to be loaded or the ray already loaded

    ion_name: string, optional
        Name of the ion to extract absorbers of
        Default: "H I"

    wavelegnth_center: float, optional
        The specific absorption line to look at (in unit Angstrom). None
        defaults to strongest absorption line for specified ion
        (using trident's ion table).
        Default: None

    velocity_res: float, optional
        Set minimum resolution (in km/s)  that spectacle will 
        attempt to fit lines to.
        Default: 10

    absorber_min: float, optional
        Minimum Log Column Density that will be used to define an absorber.
        If None, defaults to either default for specific ion or 13
        Default: None

    spectacle_defaults: dict, optional
        Dictionary passed to spectacle defining default parameters/ranges
        when fitting absorption lines
        Deafult: None

    """
    def __init__(self,
                 ds_filename,
                 ion_name='H I',
                 wavelength_center=None,
                 velocity_res=10,
                 absorber_min=None,
                 spectacle_defaults=None):
        super().__init__(ds_filename,
                         ion_name, wavelength_center,
                         velocity_res, absorber_min)

        #add user defined defaults
        if spectacle_defaults is not None:
            self.defaults_dict.update(spectacle_defaults)

    def get_current_absorbers(self):
        """
        Extract absorbers from the currently loaded ray using spectacle.

        Spectacle will fit a Trident-made spectra of the specified ion.
        Absorbers correspond to individual line fits.

        Returns
        ----------
        line_stats : pandas.DataFrame
            Table including all line statistics found from spectacle's fit of
            the spectra.
        """
        if self.ray is None:
            raise RuntimeError("You must first load a ray with `load_ray`!")

        #create spectra for a single line to fit
        wav = int( np.round(self.wavelength_center) )
        line = f"{self.ion_name} {wav}"
        #format ion correctly to fit
        ion_wav= "".join(line.split())

        vel_array, flux_array=self._create_spectra()

        #constrain possible column density values
        #create line model
        line_finder = LineFinder1D(ions=[ion_wav], continuum=1, z=0,
                                   defaults=self.defaults_dict,
                                   fitter_args={'maxiter':2000},
                                   threshold=0.01, output='flux',
                                   min_distance=self.velocity_res, 
                                   auto_fit=True)
        #fit data
        try:
            spec_model = line_finder(vel_array*u.Unit('km/s'), flux_array)
        except RuntimeError:
            print('fit failed(prolly hit max iterations)', self.ray)
            spec_model = None
        except IndexError:
            print('INDEX ERROR on', self.ray)
            spec_model = None

        #check if fit found any lines
        if spec_model is None:
            print('line could not be fit on ray ', self.ray)
            self.features = None
            line_stats = None
            self.num_feat = 0

        else:
            init_stats = spec_model.line_stats(vel_array*u.Unit('km/s'))

            # include only lines greater than absorber_min
            line_indxs, = np.where( init_stats['col_dens'] >= self.absorber_min)
            if line_indxs.size == 0:
                print('line could not be fit on ray ', self.ray)
                self.features = None
                line_stats = None
                self.num_feat = 0

            else:
                # retrieve lines that pass col dense threshold
                good_lines=[]
                for i in line_indxs:
                    good_lines.append(spec_model.lines[i])

                #create and save new model with lines desired
                self.features = spec_model.with_lines(good_lines, reset=True)
                self.num_feat = len(good_lines)
                line_stats=self.features.line_stats(vel_array*u.Unit('km/s'))

                #add redshift
                line_stats['redshift'] = self.ds.current_redshift
                line_stats = line_stats.to_pandas()
        self.df=line_stats
        return self.df

    def _create_spectra(self):
        """
        Use trident to create the absorption spectrum of the ray in velocity
        space for use in fitting.

        Returns
        --------
        velocity: YT array
            Array of velocity values of the generated spectra (in km/s)

        flux: YT array
            Array of the normalized flux values of the generated spectra

        """
        #set which ions to add to spectra
        wav = int( np.round(self.wavelength_center) )
        line = f"{self.ion_name} {wav}"
        ion_list = [line]

        #use auto feature to capture full line
        spect_gen = trident.SpectrumGenerator(lambda_min="auto", lambda_max="auto", dlambda = self.velocity_res, bin_space="velocity")
        spect_gen.make_spectrum(self.data, lines=ion_list)

        #get fields from spectra and give correct units
        flux = spect_gen.flux_field
        velocity = spect_gen.lambda_field

        return velocity, flux

class SPICEAbsorberExtractor(AbsorberExtractor):
    """
    Uses SPICE to extract absorbers from a Trident lightray for a given ion species.

    Parameters
    --------------

    ds_filename: str or YT dataset
        Either Path/name of the dataset to be loaded or the dataset itself

    ray_filename: str or Trident ray
        Path/name of the hdf5 ray file to be loaded or the ray already loaded

    ion_name: string, optional
        Name of the ion to extract absorbers of
        Default: "H I"

    wavelegnth_center: float, optional
        The specific absorption line to look at (in unit Angstrom). None
        defaults to strongest absorption line for specified ion
        (using trident's ion table).
        Default: None

    velocity_res: float, optional
        Set minimum resolution (in km/s)  that spectacle will 
        attempt to fit lines to.
        Default: 10

    absorber_min: float, optional
        Minimum Log Column Density that will be used to define an absorber.
        If None, defaults to either default for specific ion or 13
        Default: None

    frac: float, optional
        Parameter defining what fraction of the number density is being
        accounted for in each iteration of the SPICE method. Must be a number
        between 0 and 1.
        Default: 0.8

    """
    def __init__(self,
                 ds_filename, 
                 ion_name='H I', 
                 wavelength_center=None,
                 velocity_res=10,
                 absorber_min=None,
                 frac=0.8):
        super().__init__(ds_filename,
                         ion_name, wavelength_center,
                         velocity_res, absorber_min)

        self.frac = frac

    def get_current_absorbers(self, fields=[], units_dict={}):
        """
        Extract absorbers from the currently loaded ray using SPICE.

        Gas features of the absorbers are also extracted.
        The features include the column density and central velocity of the
        absorption line (delta_v) as well as requested `fields`.

        Parameters
        -----------
        fields : list, optional
            list of yt fields to extract averages of for the absorbers.
            Defalut: []

        units_dict : dict, optional
            dictionary of fields and corresponding units to use for each field.
            Default: {}

        Returns
        ---------
        absorber_info : pandas.DataFrame
            Dataframe of all the absorbers and their corresponding features.

        """
        if self.ray is None:
            raise RuntimeError("You must first load a ray with `load_ray`!")

        # get absorber locations
        self.features = self._run_spice()
        self.num_feat = len(self.features)

        # line information for absorbers
        name_type = [('name', str),
                     ('wave', np.float64),
                     ('redshift', np.float64),
                     ('col_dens', np.float64),
                     ('delta_v', np.float64),
                     ('vel_dispersion', np.float64),
                     ('interval_start', np.int32),
                     ('interval_end', np.int32)]

        # get name of columns and type of data for each
        for f in fields:
            name_type.append( (f, np.float64) )

        #check if any absorbrs were found
        n_abs = len(self.features)
        if n_abs == 0:
            print("No absorbers in ray: ", self.ray)
            return None

        #initialize empty table
        stats_table = pd.DataFrame(np.empty(n_abs , dtype=name_type))

        #add ion name and wavelength
        stats_table['name']= self.ion_name
        stats_table['wave'] = self.wavelength_center
        stats_table['redshift'] = self.ds.current_redshift

        # fill table with absorber features
        for i in range(n_abs):
            #load data for calculating properties
            start, end = self.features[i]
            stats_table.loc[i, 'interval_start'] = start
            stats_table.loc[i, 'interval_end'] = end

            dl = self.data['dl'][start:end].in_units('cm')
            density = self.data[('gas', 'density')][start:end].in_units('g/cm**3')
            tot_density = np.sum(dl*density)

            #calculate column density
            ion_field = ion_p_num(self.ion_name)
            ion_density=self.data[ion_field][start:end].in_units('cm**-3')
            col_density = np.sum(dl*ion_density)

            stats_table.loc[i, 'col_dens'] = np.log10(col_density)

            #calculate delta_v of absorber. ion col dense weighted
            vel_los_dat = self.data['velocity_los'][start:end]
            central_vel = np.sum(dl*ion_density*vel_los_dat)/col_density
            if "delta_v" in units_dict:
                central_vel.convert_to_units(units_dict["delta_v"])
            else:
                central_vel.convert_to_units('km/s')

            stats_table.loc[i, 'delta_v'] = central_vel

            #calculate velocity dispersion

            # set single cell absorber to zero velocity variance
            if end-start == 1:
                vel_variance=np.nan
            else:
                #weighted sample variance
                vel_variance=col_density \
                    * np.sum(dl*ion_density * \
                        (vel_los_dat - central_vel)**2 ) \
                    / (col_density**2 - np.sum( (dl*ion_density)**2 ))
                vel_variance = np.sqrt(vel_variance)

                if "vel_dispersion" in units_dict:
                    vel_variance.convert_to_units(units_dict["vel_dispersion"])
                else:
                    vel_variance.convert_to_units("km/s")

            stats_table.loc[i, 'vel_dispersion'] = vel_variance

            #calculate other field averages. gas col density weighted
            for fld in fields:
                fld_data = self.data[fld][start:end]
                avg_fld = np.sum(dl*density*fld_data)/tot_density

                if fld in units_dict.keys():
                    stats_table.loc[i, fld] = avg_fld.in_units( units_dict[fld] )
                else:
                    stats_table.loc[i, fld] = avg_fld

        self.df = stats_table
        return self.df

    def _run_spice(self):
        """
        iteratively run the cloud method to extract all the absorbers in the
        lightray.

        Returns
        --------
        :final_intervals: list of tuples of int
            List of the indices that indicate the start and end of each absorber.
        """
        num_density = self.data[ion_p_num(self.ion_name)].in_units("cm**(-3)")
        dl_list = self.data['dl'].in_units('cm')

        all_intervals=[]
        curr_num_density = num_density.copy()
        curr_col_density = np.sum(num_density*dl_list)
        min_col_density = 10**self.absorber_min

        while curr_col_density > min_col_density:
            #calc threshold to get fraction from current num density
            curr_thresh = self._cloud_method(curr_num_density, coldens_fraction=self.frac)

            #extract intervals this would cover
            curr_intervals = self._identify_intervals(curr_thresh)
            all_intervals = self._sensible_combination(all_intervals, curr_intervals)

            #mask density array above threshold and apply mask to dl
            curr_num_density = np.ma.masked_greater_equal(num_density, curr_thresh)
            curr_dl =  np.ma.masked_array(dl_list, curr_num_density.mask)

            #calc leftover column density
            curr_col_density = np.sum(curr_num_density*curr_dl)

        #make sure intervals have high enough col density
        final_intervals=[]
        for b, e in all_intervals:
            lcd = np.log10(np.sum(dl_list[b:e]*num_density[b:e]))
            if lcd > self.absorber_min:
                final_intervals.append((b, e))
        return final_intervals

    def _cloud_method(self, num_density_arr, coldens_fraction):
        "run the cloud method"
        cut = 0.999
        total = np.sum(num_density_arr)
        ratio = 0.001
        while ratio < coldens_fraction:
            part = np.sum(num_density_arr[num_density_arr > cut * np.max(num_density_arr)])
            ratio = part / total
            cut = cut - 0.001

        threshold = cut * np.max(num_density_arr)

        return threshold

    def _sensible_combination(self, prev_intervals, curr_intervals):
        """
        adds new intervals by taking into account the velocities when combining them

        Parameters
        -----------
        prev_intervals : list
            the intervals already calculated

        curr_intervals : list
            the intervals that need to be added/combined

        Returns
        --------
        new_intervals : list
            a final list of intervals where prev and curr are properly combined.
        """
        dl_array = self.data['dl'].in_units('cm')
        velocity_array = self.data['velocity_los'].in_units('km/s')
        density_array = self.data['density']

        #check if there are any previous intervals to combine with
        if prev_intervals == []:
            return curr_intervals

        new_intervals=prev_intervals.copy()
        del_v = self.ds.quan(self.velocity_res, 'km/s')

        #loop through current intervals
        for curr_b, curr_e in curr_intervals:

            overlap_intervals=[]
            #loop through all previous intervals
            for b,e in prev_intervals:
                #check if previous interval is nested in curr interval
                if curr_b <= b and curr_e >= b:
                    #print(f"interval ({curr_b}, {curr_e}) overlap with ", b, e)
                    if curr_b <= e and curr_e >= e:
                        overlap_intervals.append((b, e))

                    #check if just beginning point enclose
                    else:
                        err_file = open("error_file.txt", 'a')
                        err_file.write(f"{self.ray_filename} had an intersection that wasn't complete :/")
                        err_file.close()
                #check if just endpoint enclosed
                elif curr_b <= e and curr_e >= e:
                    err_file = open("error_file.txt", 'a')
                    err_file.write(f"{self.ray_filename} had an intersection that wasn't complete :/")
                    err_file.close()

            #
            #
            #This is such a mess below but it works
            #hopefully I'll think of a much cleaner way to do this
            #but for now this is it
            #
            #

            if overlap_intervals == []:
                new_intervals.append((curr_b, curr_e))
            else:

                #collect overlap points into list
                points = [curr_b]
                for b, e in overlap_intervals:
                    #print(f"curr {curr_b, curr_e} ovelaps {b, e}")
                    new_intervals.remove((b, e))
                    points.append(b)
                    points.append(e)
                points.append(curr_e)

                avg_v=[]
                for i in range(len(points)-1):
                    pnt1, pnt2 = points[i], points[i+1]
                    #find weighted avg velocity
                    curr_dense = density_array[pnt1:pnt2]
                    curr_dl = dl_array[pnt1:pnt2]
                    vel = np.sum(curr_dense*curr_dl*velocity_array[pnt1:pnt2]) \
                          /np.sum(curr_dense*curr_dl)
                    avg_v.append((vel, pnt1, pnt2))

                start_b = curr_b
                for i in range(len(avg_v)-1):
                    #if velocity difference is greater than threshold
                    if abs(avg_v[i][0] - avg_v[i+1][0]) > del_v:
                        #create new interval
                        new_intervals.append((start_b, avg_v[i][2]))
                        #change start of next interval
                        start_b = avg_v[i+1][1]
                    #check if this is the last two intervals to check
                    if i == len(avg_v) -2:
                        new_intervals.append((start_b, curr_e))


        return new_intervals

    def _identify_intervals(self, cutoff):
        """
        Find the intervals for absorbers using some cutoff on the number density
        field along lightray.

        Parameters
        -----------
        cutoff : double
            threshold defining where absorbers are.

        Returns
        -----------
        intervals : list of tuples
            list of the intervals defining the absorbers in this ray.
        """
        num_density = self.data[ion_p_num(self.ion_name)].in_units("cm**(-3)")
        in_absorber = False
        intervals = []

        #Iterate over values in field
        for i,value in enumerate(num_density):
            #Check if started an absorber and if above cutoff
            if in_absorber and value < cutoff:
                in_absorber = False
                #add interval to list
                intervals.append((start,i))
            # check if just entered an absorber
            elif not in_absorber and value >= cutoff:
                in_absorber = True
                start = i
            else:
                continue
        #check if was still in absorber when hitting end of ray
        if in_absorber and start != i:
            intervals.append((start, i))
        return intervals

    def get_all_absorbers(self, ray_list, fields=[], units_dict={}):
        """
        Create catalog of the given rays using absorber extractor

        Parameters
        ----------
        ray_list: list of str or trident.ray objects
            List of ray objects or list of trident rays whose absorbers will be
            extracted

        fields : list, optional
            list of yt fields to extract averages of for the absorbers.
            Defalut: []

        units_dict : dict, optional
            dictionary of fields and corresponding units to use for each field.
            Default: {}

        Returns
        -------
        full_df: pandas.DataFrame
            Catalog of absorber properties in a pandas dataframe.
        """
        return super().get_all_absorbers(ray_list, fields, units_dict)
