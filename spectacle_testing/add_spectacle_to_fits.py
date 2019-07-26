'''

Reads in a MISTY fits file that may or may not have spectacle run on it yet;
runs spectacle and creates new fits file with new header information

'''

from __future__ import print_function

import glob
import os

import MISTY
from spectacle.analysis import Resample
from spectacle.modeling.lsfs import GaussianLSFModel

os.sys.path.insert(0, '/Users/molly/Dropbox/foggie/foggie')
from plot_misty_spectra import plot_misty_spectra

import numpy as np

from astropy.io import fits
import astropy.units as u
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.table import Table

from scipy.signal import argrelextrema

def add_spectacle_to_fits(old_fits_name, new_fits_name, **kwargs):
    threshold = kwargs.get('threshold', 0.01)
    plotname = kwargs.get('plotname', 'temp.png')
    line_list = kwargs.get('line_list', ['H I 1216', 'H I 919', \
             'Si II 1260', 'Si IV 1394', 'C IV 1548', 'O VI 1032'])
    plot = kwargs.get('plot', False)

    print('only doing lines in line_list! line_list is: ', line_list)

    orig_hdu = fits.open(old_fits_name)
    new_hdu = fits.HDUList([orig_hdu[0]])
    new_hdu.append(orig_hdu[1])

    keys_to_copy = ('LINENAME',
                    'RESTWAVE',
                    'F_VALUE',
                    'GAMMA',
                    'SIM_TAU_HDENS',
                    'SIM_TAU_TEMP',
                    'SIM_TAU_METAL',
                    'TOT_COLUMN',
                    'EXTNAME',
                    'XTENSION',     ## hack-y since don't want to not delete these
                    'BITPIX',
                    'NAXIS',
                    'NAXIS1',
                    'NAXIS2',
                    'PCOUNT',
                    'GCOUNT',
                    'TFIELDS',
                    'TTYPE1',
                    'TFORM1',
                    'TTYPE2',
                    'TFORM2',
                    'TUNIT2',
                    'TTYPE3',
                    'TFORM3',
                    'TTYPE4',
                    'TFORM4',
                    'TTYPE5',
                    'TFORM5',
                    'TTYPE6',
                    'TFORM6',
                    'TTYPE7',
                    'TFORM7',
                    'TTYPE8',
                    'TFORM8',
                    'TTYPE9',
                    'TFORM9',
                    'TTYPE10',
                    'TFORM10')

    ## now for the individual lines
    nlines = np.int(orig_hdu[0].header['NLINES'])
    for line_num in np.arange(nlines):
        key = 'LINE_'+str(line_num+1)
        line_name = orig_hdu[0].header[key]
        if line_name in line_list:
            print('~~~~> trying',line_name,'~~~~~>>>>')

            if any([x.name.upper() == line_name.upper() for x in orig_hdu]):
                new_ext = orig_hdu[line_name]
                for k in orig_hdu[line_name].header:
                    if k not in keys_to_copy:
                        print("deleting ", k)
                        del new_ext.header[k]

                lambda_0 = orig_hdu[line_name].header['RESTWAVE']
                try:
                    disp = orig_hdu[line_name].data['wavelength']
                    flux = orig_hdu[line_name].data['flux']
                    tau = orig_hdu[line_name].data['tau']
                    redshift = orig_hdu[line_name].data['redshift']
                except:
                    disp = orig_hdu[line_name].data['disp_obs']
                    flux = orig_hdu[line_name].data['flux_obs']
                    tau = orig_hdu[line_name].data['tau_obs']
                    redshift = orig_hdu[line_name].data['redshift_obs']

                zsnap = np.median(redshift)

                ## we want Nmin for a range of thresholds
                Nmin = np.size(np.where(flux[argrelextrema(flux, np.less)[0]] < (1-threshold)))
                new_ext.header['Nmin'] = Nmin
                print('found ', Nmin, ' minima')
                Nmin = np.size(np.where(flux[argrelextrema(flux, np.less)[0]] < (1-0.01)))
                new_ext.header['Nmin001'] = Nmin
                Nmin = np.size(np.where(flux[argrelextrema(flux, np.less)[0]] < (1-0.02)))
                new_ext.header['Nmin002'] = Nmin
                Nmin = np.size(np.where(flux[argrelextrema(flux, np.less)[0]] < (1-0.05)))
                new_ext.header['Nmin005'] = Nmin
                Nmin = np.size(np.where(flux[argrelextrema(flux, np.less)[0]] < (1-0.1)))
                new_ext.header['Nmin010'] = Nmin

                #### lots of assumptions here
                fwhm = 7. # km/s
                dv = 2.
                sigma= (fwhm/dv) / (2*np.sqrt(2.*np.log(2.)))
                lsf = GaussianLSFModel(stddev = sigma)
                print("""
                I AM ASSUMING YOU HAVE AN LSF AND SO IF YOU DON'T YOU SHOULD FIX THIS!!!!!!!!!!!
                """)

                print("~~~~> now trying to run spectacle on line ",line_name, "~~~~~~>")
                lines_properties = MISTY.get_line_info(disp, flux, \
                                                tau=tau, \
                                                redshift=zsnap, \
                                                lambda_0=lambda_0, \
                                                f_value=orig_hdu[line_name].header['F_VALUE'], \
                                                gamma=orig_hdu[line_name].header['GAMMA'], \
                                                ion_name=line_name, \
                                                threshold = threshold,
                                                lsf=lsf)
                                                ##### pass in the LSF here !!!
                print(lines_properties)


                for line_key in lines_properties:
                    if isinstance(lines_properties[line_key], tuple):
                        if np.isnan(lines_properties[line_key][0]):
                            lines_properties[line_key] = -99.
                    new_ext.header[line_key] = lines_properties[line_key]


                new_hdu.append(new_ext)
                print('~~~~> all done with',line_name,'~~~~~<<<')
        else:
            print('<<<<<~~~~ ',line_name,' not found or not in line_list :-(  ~~~~~<<<<<<')

    print("writing out to .... " + new_fits_name)
    new_hdu.writeto(new_fits_name, overwrite=True, output_verify='fix')

    #print('plotting to... ' + plotname)
    #plot_misty_spectra(new_hdu, overplot=True, outname=plotname)



if __name__ == "__main__":

    long_dataset_list = glob.glob(os.path.join(".", 'hlsp*rd0018*axx*v6_lsf.fits.gz'))
    ##long_dataset_list = ['./hlsp_misty_foggie_halo008508_nref11n_nref10f_rd0020_axy_dx042.9_dz082.4_v6_lsf.fits.gz']
    dataset_list = long_dataset_list

    for filename in dataset_list:
        new_filename = '.' + filename.strip('lsf.fits.gz') + 'lsf.fits.gz'
        plotname = '.' + new_filename.strip('.lsf.fits.gz') + 'lsf.png'
        print('adding spectacle to ', filename, ' and saving as ', new_filename)
        add_spectacle_to_fits(filename, new_filename, plot=True, plotname=plotname, threshold=0.05)
