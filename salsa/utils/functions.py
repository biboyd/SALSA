from trident import from_roman


cgm_in_radius=10 #kpc
cgm_out_radius=200 #kpc
cgm_density=2e-26 # g/cm^3
cgm_temp=1.5e4 #K

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

def ion_p(ion_name):
    """
    convert ion species name from trident style to one that
    matches YT's (ie 'H I' --> 'H_p0')
    """
    #split up the words in ion species name
    ion_split = ion_name.split()
    #convert num from roman numeral. subtract run b/c YT
    num = from_roman(ion_split[1])-1

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
