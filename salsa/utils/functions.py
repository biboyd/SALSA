from trident import from_roman


cgm_in_radius=10 #kpc
cgm_out_radius=200 #kpc
cgm_density=2e-26 # g/cm^3
cgm_temp=1.5e4 #K

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

def requires_spectacle(function):
    def wrapped(*args, **kwargs):
        try:
            import spectacle as _
        except ImportError:
            raise ImportError("The spectacle package is required for this feature. Please install salsa[spectacle] for the necessary dependencies.")
        else:
            result = function(*args, **kwargs)
            
        return result
    return wrapped