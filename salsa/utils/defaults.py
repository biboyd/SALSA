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

default_cloud_dict = {'H I': 12.5, 'C IV':13, 'O VI':12.8}
