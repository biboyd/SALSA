import yt
import trident
import numpy as np


cgm_in_radius=10 #kpc
cgm_out_radius=200 #kpc
cgm_density=2e-26 # g/cm^3
cgm_temp=1.5e4 #K
#function to create field in yt
def radius_function(field, data):
    if data.has_field_parameter("center"):
        c = data.get_field_parameter("center")
    else:
        c = data.ds.domain_center

    x = data[('gas', 'x')] - c[0]
    y = data[('gas', 'y')] - c[1]
    z = data[('gas', 'z')] - c[2]
    return np.sqrt(x*x + y*y + z*z)

def parse_cut_filter(cuts):
    """
    Parses a string defining what cuts to apply to analysis. Then returns the proper
    YT cut_region filter
    """
    #define cut dictionary
    cgm_rad = f"(obj[('gas', 'radius')].in_units('kpc') > {cgm_in_radius}.) & (obj[('gas', 'radius')].in_units('kpc') < {cgm_out_radius})"
    cgm_d_t = f"(obj[('gas', 'temperature')].in_units('K') > {cgm_temp}) | (obj[('gas', 'density')].in_units('g/cm**3') < {cgm_density})"
    cgm_filter = f"{cgm_rad} & {cgm_d_t}"

    ism_in = f"(obj[('gas', 'radius')].in_units('kpc') < {cgm_in_radius}.)"
    ism_cold_dense = f"(obj[('gas', 'radius')].in_units('kpc') < {cgm_out_radius}.) & (obj[('gas', 'temperature')].in_units('K') < {cgm_temp}) & (obj[('gas', 'density')].in_units('g/cm**3') > {cgm_density})"

    ism_filter = f"({ism_in}) | ({ism_cold_dense})"# either r<10 kpc or r<200kpc and d/t criteria

    hot = "(obj[('gas', 'temperature')].in_units('K') > 1.e5)"
    cold = "(obj[('gas', 'temperature')].in_units('K') <= 1.e5)"

    inflow = "(obj[('gas', 'radial_velocity')] <= 0.)"
    outflow = "(obj[('gas', 'radial_velocity')] > 0.)"

    high_OVI = "(obj[('gas', 'O_p5_ion_fraction')] >= 0.1)"
    low_OVI = "(obj[('gas', 'O_p5_ion_fraction')] < 0.1)"

    filter_dict=dict(cgm=cgm_filter, ism=ism_filter,
                     hot=hot, cold=cold,
                     inflow=inflow, outflow=outflow,
                     high_OVI=high_OVI, low_OVI=low_OVI)

    #split filter call
    filter_names = cuts.split(' ')
    cut_filters = [ filter_dict[name] for name in filter_names ]

    return cut_filters

cut_alias_dict = dict(cgm="CGM", ism="ISM",
                 cgm_hot="Hot T>1e5", cgm_cold="Cold T<1e5",
                 cgm_inflow="Inflow", cgm_outflow="Outflow",
                 cgm_high_OVI="f_ovi > 0.1", cgm_low_OVI="f_ovi < 0.1",
                 cgm_hot_inflow="Hot Inflow T>1e5", cgm_hot_outflow="Hot Outflow T>1e5",
                 cgm_cold_inflow="Cold Inflow T>1e5", cgm_cold_outflow="Cold Outflow T<1e5")
# labels to use in final plot instead

axis_labels_dict={'log_density':"Log( Density ) ($g/cm^3$)",
                  'log_metallicity':"Log( Metallicity ) ($Z_{\odot}$)",
                  'radius':'Radial Distance ($kpc$)',
                  'log_temperature':"Log( Temperature ) ($K$)",
                  'col_dens':"Log Column Density",
                  'radial_velocity': "Radial Velocity ($km/s$)"}

# Histogram limits dictionary
ovi_hist_dict = dict(col_dens=(13., 16.),
                      log_metallicity=(-2., 0.05),
                      log_temperature=(4., 7.),
                      log_density=(-28.75, -25))
hist_range_dict = {"O VI":ovi_hist_dict}

def ion_p(ion_name):
    """
    convert ion species name from trident style to one that
    matches YT's (ie 'H I' --> 'H_p0')
    """
    #split up the words in ion species name
    ion_split = ion_name.split()
    #convert num from roman numeral. subtract run b/c YT
    num = trident.from_roman(ion_split[1])-1

    #combine all the names
    return f"{ion_split[0]}_p{num}"

def ion_p_num(ion_name):
    """
    convert ion species name from trident style to one that
    matches YT's (ie 'H I' --> 'H_p0_number_density')
    """
    ip = ion_p(ion_name)
    #combine all the names
    outname = f"{ip}_number_density"
    return outname

#ice extraction defaults
default_ice_fields = ['x', 'y', 'z','radius', 'density', 'metallicity', 'temperature', 'radial_velocity']
default_units_dict=dict(velocity_los='km/s',
                       x='kpc',
                       y='kpc',
                       z='kpc',
                       radius='kpc',
                       density='g/cm**3',
                       metallicity='Zsun',
                       temperature='K',
                       radial_velocity='km/s')

default_limits_dict = dict(velocity_los=[-600, 600],
                           metallicity=[0, 1],
                           temperature=[1e4, 1e9],
                           density=[1e-30, 1e-26])
