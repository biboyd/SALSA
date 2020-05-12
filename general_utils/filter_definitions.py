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

    hot = "(obj[('gas', 'temperature')].in_units('K') > 1.e5)"
    cold = "(obj[('gas', 'temperature')].in_units('K') <= 1.e5)"

    inflow = "(obj[('gas', 'radial_velocity')] <= 0.)"
    outflow = "(obj[('gas', 'radial_velocity')] > 0.)"

    filter_dict=dict(cgm=cgm_filter, hot=hot, cold=cold, inflow=inflow, outflow=outflow)

    #split filter call
    filter_names = cuts.split(' ')
    cut_filters = [ filter_dict[name] for name in filter_names ]
    #cut_filter = "&".join(filters)

    return cut_filters

def ion_p_num(ion_name):
    """
    convert ion species name from trident style to one that
    matches YT's (ie 'H I' --> 'H_p0_number_density')
    """
    #split up the words in ion species name
    ion_split = ion_name.split()
    #convert num from roman numeral. subtract run b/c YT
    num = trident.from_roman(ion_split[1])-1

    #combine all the names
    outname = f"{ion_split[0]}_p{num}_number_density"
    return outname

#ice extraction defaults
default_ice_fields = ['x', 'y', 'z','radius', 'density', 'metallicity', 'temperature', 'radial_velocity']
default_unit_dict=dict(velocity_los='km/s',
                       x='kpc',
                       y='kpc',
                       z='kpc',
                       radius='kpc',
                       density='g/cm**3',
                       metallicity='Zsun',
                       temperature='K',
                       radial_velocity='km/s')
