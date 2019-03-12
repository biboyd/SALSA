import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import yt
yt.enable_parallelism()
yt.funcs.mylog.setLevel(50)
from yt.units import *
import trident

import math
import subprocess
from subprocess import call
import os, sys, glob, errno
import logging

import h5py
import PIL
try:
    import Image
except ImportError:
    from PIL import Image

# return formatted string for ion input string
def ion_p_name(ion):
    return ion.split()[0]+'_p'+str(trident.from_roman(ion.split()[1])-1)

# remove space from yt quan; accepts h5file
def yt_remove_spaces(ytquan):
    return str(ytquan).replace(" ", "")


# creates 3-panel plots of column density vs. impact parameter
# compares median, mean, standard deviation in projection along x- and z- axes
# creates cdip plots naively, only looking along two projection axes
#
# takes arguments containing:
# pattern matching for data dumps, number of processes to parallelize over, ion species names, center coordinate,
# center units, radius length, radius units, pixel number (width of image), number of radial bins,
# y limits for each plot, and the header for the output filenames
# 'output_filename_head' can be used to place files inside specific directories (end with a '/')
def create_cdip_xz_proj(filenames,
                        num_procs,
                        ion_species = None,
                        center_coord = None,
                        center_units = None,
                        radius_length = None,
                        radius_units = None,
                        pix_num = None,
                        num_bins = None,
                        y_lim1 = None,
                        y_lim2 = None,
                        y_lim3 = None,
                        output_filename_head = './'):
    if ion_species is None:
        ion_species = ['O VI', 'H I', 'C IV']
    if center_coord is None:
        center_coord = [0.5, 0.5, 0.5]
    if center_units is None:
        center_units = "code_length"
    if radius_length is None:
        radius_length = 300
    if radius_units is None:
        radius_units = "kpc"
    if pix_num is None:
        pix_num = 1000
    if num_bins is None:
        num_bins = 10
    if y_lim1 is None:
        y_lim1 = [0,0]
    if y_lim2 is None:
        y_lim2 = [0,0]
    if y_lim3 is None:
        y_lim3 = [0,0]
    if output_filename_head is None:
        output_filename_head ='./'

    # print("----")
    # print(filenames)
    # print(num_procs)
    # print(ion_species)
    # print(center_coord,center_units)
    # print(radius_length,radius_units)
    # print(pix_num, num_bins)
    # print(y_lim1,y_lim2,y_lim3)
    # print(output_filename_head)
    # print("----")

    # set up logging script parameters
    logging.info('Filenames: %s' % filenames)
    logging.info('Number of processes: %s' % num_procs)
    logging.info('Ion species: %s' % ion_species)
    logging.info('Center: %s %s' % (center_coord,center_units))
    logging.info('Radius: %s %s' % (radius_length,radius_units))
    logging.info('Pixel length/width: %s' % pix_num)
    logging.info('Number of radial bins: %s' % num_bins)
    logging.info('Y-limits, if any set: %s,%s,%s' % (y_lim1,y_lim2,y_lim3))
    logging.info('Output filename head: %s' % output_filename_head)


    ts = yt.DatasetSeries(filenames, parallel = num_procs)
    logging.info('Stop 1')

    for ds in ts.piter():
        logging.info('Stop 2')
        trident.add_ion_fields(ds, ions=ion_species)

        width = ds.quan(radius_length*2., radius_units) # width

        # make cube to limit projection region
        center = ds.arr(center_coord,center_units)
        offset = ds.arr([radius_length,  radius_length,  radius_length],radius_units)
        left_edge = center.in_units(radius_units)-offset
        right_edge = center.in_units(radius_units)+offset
        cube = ds.region(center,left_edge,right_edge)
        d = radius_length*2.

        # iterate over each ion
        for ion in ion_species:
            cdip_filename = (output_filename_head+'cdip_xz_'+ion_p_name(ion)+'_%s_'+str(int(d))+radius_units+'_box.png') % ds
            if not os.path.exists(cdip_filename):

                ion_number_density_name = ion_p_name(ion)+'_number_density'
                proj_ion_0 = yt.ProjectionPlot(ds,0,ion_number_density_name,center,data_source = cube)
                proj_ion_2 = yt.ProjectionPlot(ds,2,ion_number_density_name,center,data_source = cube)


                # pix_num = 1000
                center_pix = (pix_num - 1) / 2.

                # num_bins = 10
                dr = pix_num / (2 * num_bins)
                #num_bins = pix_num / dr

                x = np.arange(0,pix_num)
                r = np.arange(0,num_bins) #* dr
                y_med_ion_0 = []
                y_med_ion_2 = []
                y_mean_ion_0 = []
                y_mean_ion_2 = []
                y_std_ion_0 = []
                y_std_ion_2 = []

                #make list of lists by bin
                val_by_r_ion_0 = [[] for _ in r]
                val_by_r_ion_2 = [[] for _ in r]

                res = [pix_num, pix_num] # create an image with 1000x1000 pixels

                frb_ion_0 = proj_ion_0.data_source.to_frb(width, res)
                proj_array_ion_0 = np.array(frb_ion_0[ion_number_density_name])
                frb_ion_2 = proj_ion_2.data_source.to_frb(width, res)
                proj_array_ion_2 = np.array(frb_ion_2[ion_number_density_name])

                #fill in val_by_r with entire image
                for a in x:
                    for b in x:
                        r_bin = np.floor(np.sqrt((a-center_pix)**2+(b-center_pix)**2) / dr)
                        if r_bin <= num_bins - 1:
                            val_by_r_ion_0[int(r_bin)].append(proj_array_ion_0[a][b])
                            val_by_r_ion_2[int(r_bin)].append(proj_array_ion_2[a][b])

                #reduce val_by_r to median
                for b in r:
                    y_med_ion_0.append(np.median(val_by_r_ion_0[b]))
                    y_mean_ion_0.append(np.mean(val_by_r_ion_0[b]))
                    y_std_ion_0.append(np.std(val_by_r_ion_0[b]))

                    y_med_ion_2.append(np.median(val_by_r_ion_2[b]))
                    y_mean_ion_2.append(np.mean(val_by_r_ion_2[b]))
                    y_std_ion_2.append(np.std(val_by_r_ion_2[b]))

                ###
                # FORMAT FINAL IMAGE
                ###



                plt.close('all')

                fig = plt.figure()
                fig.set_size_inches(6, 12)
                ax = fig.add_subplot(111)
                # Turn off axis lines and ticks of the big subplot
                ax.spines['top'].set_color('none')
                ax.spines['bottom'].set_color('none')
                ax.spines['left'].set_color('none')
                ax.spines['right'].set_color('none')
                ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

                ax1 = fig.add_subplot(311)
                ax2 = fig.add_subplot(312)
                ax3 = fig.add_subplot(313)

                #fig, axarr = plt.subplots(2, 2)
                #fig.tight_layout()
                fig.tight_layout(pad=4, h_pad=3, w_pad=2)
                ax1.set_yscale('log')
                ax2.set_yscale('log')
                ax3.set_yscale('log')

                bar_width = dr * d / pix_num

                ax1.bar((r+0.5) * bar_width, y_med_ion_0, bar_width, color='0.4')
                ax1.errorbar((r+0.5) * bar_width, y_mean_ion_0, yerr=y_std_ion_0, fmt="o", color='0.0', capsize=bar_width/4)
                ax1.set_title(ion+' Projected in x')
                ax1.set_xlim([0,d/2])
                if y_lim1 != [0,0]:
                    ax1.set_ylim(np.float_(y_lim1))

                ax2.bar((r+0.5) * bar_width, y_med_ion_2, bar_width, color='0.8')
                ax2.errorbar((r+0.5) * bar_width, y_mean_ion_2, yerr=y_std_ion_2, fmt="o", color='0.0', capsize=bar_width/4)
                ax2.set_title(ion+' Projected in z')
                ax2.set_xlim([0,d/2])
                if y_lim2 != [0,0]:
                    ax2.set_ylim(np.float_(y_lim2))

                ax3.plot((r+0.5) * bar_width, y_med_ion_0, "-o", color='0.4')
                ax3.plot((r+0.5) * bar_width, y_med_ion_2, "-o", color='0.8')
                ax3.set_title(ion+' for each projection')
                ax3.set_xlim([0,d/2])
                if y_lim3 != [0,0]:
                    ax3.set_ylim(np.float_(y_lim3))

                if ion == 'O VI':
                    ax1.set_ylim(np.float_([1e13,1e17]))
                    ax2.set_ylim(np.float_([1e13,1e17]))
                    ax3.set_ylim(np.float_([1e13,1e17]))
                if ion == 'H I':
                    ax1.set_ylim(np.float_([1e13,1e22]))
                    ax2.set_ylim(np.float_([1e13,1e22]))
                    ax3.set_ylim(np.float_([1e13,1e22]))
                if ion == 'C IV':
                    ax1.set_ylim(np.float_([1e11,1e17]))
                    ax2.set_ylim(np.float_([1e11,1e17]))
                    ax3.set_ylim(np.float_([1e11,1e17]))

                ax.set_title(ion+': Time Step %s' % ds, y=1.03)

                ax.set_xlabel('Impact Parameter (kpc)', labelpad=20)
                ax.set_ylabel(r'Column Density (cm$^{-2}$)', labelpad=20)
                logging.info('Stop 3')
                plt.savefig(cdip_filename)
                logging.info('Stop 4')
                logging.debug('Data dump %s, ion %s complete' % (ds,ion))
                logging.info('Stop 5')

                if yt.is_root():
                    print("test",ds)
                    file = open(output_filename_head+"filenames_cdip_xz_proj.fn",'a+')
                    file.write(filenames[:-13]+"%s/%s -- %s\n" % (ds,ds,ion))
                    file.close()

# Arguments: number of points, sightline length, inner radius, outer radius, yt_bool
# Set yt_bool to true if you're generating sightlines on a loaded dataset ds
# Outputs: list of pairs of endpoints (list of list of lists)
def generate_sightlines(ds,
                        num_points,
                        sightline,
                        center,
                        in_radius,
                        out_radius,
                        yt_bool):
    sightline_array = []
    delta = 5.551115123125783e-17

    if yt_bool == True:
        length = sightline.in_units("code_length")
        c = center.in_units("code_length")
        ri = in_radius.in_units("code_length")
        ro = out_radius.in_units("code_length")
    else:
        length = sightline
        c = center
        ri = in_radius
        ro = out_radius

    #length = sightline_length.in_units("unitary")
    #c = center.in_units("unitary")
    #ri = in_radius.in_units("unitary")
    #ro = out_radius.in_units("unitary")
    for _ in np.arange(num_points):
        # find tangent point within shell
        R = np.power(np.random.uniform(ri**3.,ro**3.),1./3.)
        theta = np.random.uniform(0., 2.*np.pi)
        Z = np.random.uniform(-1.,1.) * R
        X = np.sqrt(1. - (Z/R)**2.) * np.cos(theta) * R
        Y = np.sqrt(1. - (Z/R)**2.) * np.sin(theta) * R

        # establish radial vector r / perpendicular vectors r1 & r2
        if yt_bool == True:
            r = ds.arr([X, Y, Z],"code_length")
        else:
            r = [X,Y,Z]
        #print(r)
        r1_normalization = np.sqrt(X**2. + Y**2.)
        r1 = np.divide([Y, -X, 0], r1_normalization)
        #print([Y, -X, 0],r1_normalization,r1)
        r2temp = np.cross(r, r1)
        r2 = np.divide(r2temp,np.linalg.norm(r2temp))

        # establish direction vector d along random angle alpha from 0 to 2pi
        # d is perpendicular to radial r vector
        alpha = np.random.uniform(0., 2.*np.pi)
        d = r1 * np.cos(alpha) + r2 * np.sin(alpha)
        dnorm = np.linalg.norm(d)
        #d = ds.arr(d, "code_length")

        # calculate endpoints, constrained by [0,1]x[0,1]x[0,1] box
        p_start=[sum(i) for i in zip(c,r,-np.multiply(np.divide(d,dnorm),length/2.))]
        p_end=[sum(i) for i in zip(c,r,np.multiply(np.divide(d,dnorm),length/2.))]

        if yt_bool == True:
            zero = ds.quan(0.,"code_length")
            one = ds.quan(1.,"code_length")
        else:
            zero = 0.
            one = 1.
        # check whether or not the calculated endpoints are within the box
        t_list = [(-c[0]-r[0]+zero)/d[0],(-c[1]-r[1]+zero)/d[1],(-c[2]-r[2]+zero)/d[2],
                  (-c[0]-r[0]+one)/d[0],(-c[1]-r[1]+one)/d[1],(-c[2]-r[2]+one)/d[2]]
        t_list.sort()
        for element_start in p_start:
            if element_start<zero or element_start>one:
                p_start = [sum(i) for i in zip(c,r,np.multiply(t_list[2],d))]
        for element_end in p_end:
            if element_end<zero or element_end>one:
                p_end = [sum(i) for i in zip(c,r,np.multiply(t_list[3],d))]
        # handle floating point error to ensure to doesn't go outside valid domain
        for element in (p_start+p_end):
            if (np.abs(element)<=delta):
                element = zero
            if (np.abs(element - one)<=delta):
                element = one

    # add to final list
        sightline_array.append([p_start,p_end,ds.quan(R,"code_length")]) #last part of array is the start end and R
    return sightline_array

###### put it in terms of ds
# line light, N, sightline length, plot center, inner radius, outer radius
def generate_lightrays(ds,
                       num_points,
                       ion_species,
                       sightline,
                       center,
                       in_radius,
                       out_radius,
                        output_filename_head = './'):

    try:
        os.makedirs(output_filename_head+"rays")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    allpoints = generate_sightlines(ds,
                                    num_points,
                                    sightline,
                                    center,
                                    in_radius,
                                    out_radius,
                                    True)

    #ray_list = []
    index_list = np.arange(0,num_points)
    num_procs = 0
    for index in yt.parallel_objects(index_list,num_procs):
        ray_filename = (output_filename_head+"rays/ray_r%s-%s_%s_%s.h5") % (yt_remove_spaces(in_radius),
                                            yt_remove_spaces(out_radius),ds,index)
        # create initial simple ray
        if not os.path.exists(ray_filename):
            ray = trident.make_simple_ray(ds,
                                          start_position=allpoints[index][0],
                                          end_position=allpoints[index][1],
                                          data_filename=ray_filename,
                                          lines=ion_species,
                                          ftype='gas')
        # create impact_parameter dataset in ray file
        f = h5py.File(ray_filename,'a')
        if not ('/grid/impact_parameter' in f):
            R = allpoints[index][2]
            R.write_hdf5(ray_filename, dataset_name='impact_parameter',group_name='grid')

            #ray_list.append(ray)
            print("%s,%s,%s,%s ray done" % (in_radius,out_radius,ds,index))
        f.close()

    return allpoints[:][2]

# in and out radius in terms of ds
# (accepts 1 ion)
def create_spectrum(ds,
                      index,
                      ion,
                      in_radius,
                      out_radius,
                      spec_lambda_min = 1000,
                      spec_lambda_max = 1050,
                      spec_dlambda = 0.1,
                        output_filename_head = './'):
    try:
        os.makedirs(output_filename_head+"spectra")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    ray = yt.load(output_filename_head+"rays/ray_r%s-%s_%s_%s.h5" % (yt_remove_spaces(in_radius),
        yt_remove_spaces(out_radius),ds,index))
    # create spectrum from sightline
    sg = trident.SpectrumGenerator(spec_lambda_min,spec_lambda_max,spec_dlambda)
    sg.make_spectrum(ray, lines=[ion])
    sg.plot_spectrum(output_filename_head+'spectra/spec_%s_r%s-%s_%s_%s.png' % (ion_p_name(ion),yt_remove_spaces(in_radius),
        yt_remove_spaces(out_radius),ds,index))

# in and out radius in terms of ds
# (accepts 1 ion)
def create_ndpl_plot(ds,
                      index,
                      ion,
                      in_radius,
                      out_radius,
                      ndpl_y_lim = (1e-12,1e-6),
                      ndpl_log = False,
                        output_filename_head = './'):
    try:
        os.makedirs(output_filename_head+"ndpl")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    f = h5py.File(output_filename_head+"rays/ray_r%s-%s_%s_%s.h5" % (yt_remove_spaces(in_radius),
        yt_remove_spaces(out_radius),ds,index), "r")
    y = list(f['grid'][ion_p_name(ion)+'_number_density'])
    dl_list = list(f['grid']['dl'])
    x = [dl_list[0]]
    for i in np.arange(1,len(dl_list)):
        temp = x[i-1]+dl_list[i]
        x.append(temp)
    #print(len(x),len(dl_list),len(y))
    plt.close('all')

    fig = plt.figure()
    fig.set_size_inches(12, 6)
    if ndpl_log:
        plt.yscale('log')
    #plt.xlim(0,1)
    if ndpl_y_lim != [0,0]:
        plt.ylim(np.float_(y_lim))
    #plt.ylim(1e-12,1e-6)
    plt.plot(x,y,'-')
    f.close()

    plt.savefig(output_filename_head+'ndpl/ndpl_%s_r%s-%s_%s_%s.png' % (ion_p_name(ion),yt_remove_spaces(in_radius),
        yt_remove_spaces(out_radius),ds,index))

# produces a column density from ndpl along a sightline
# (accepts 1 ion)
def integrate_ndpl(ds,
                    index,
                    ion,
                    in_radius,
                    out_radius,
                    output_filename_head = './'):
    f = h5py.File(output_filename_head+"rays/ray_r%s-%s_%s_%s.h5" % (yt_remove_spaces(in_radius),
        yt_remove_spaces(out_radius),ds,index), "r")
    # print(f,"test 1")

    # ad hoc
    ion_prefix = ion_p_name(ion)
    if ion_p_name(ion)=="H_p0":
        ion_prefix = "H"
    y = list(f['grid'][ion_prefix+'_number_density'])
    # print(y,'test 2')
    dl_list = list(f['grid']['dl'])
    # print(dl_list,'test 3')
    cd_temp = y[0] * dl_list[0]
    for i in np.arange(1,len(dl_list)):
        cd_temp += y[i] * dl_list[i]
    f.close()
    # print('pls')
    return cd_temp

# creates projection plot along x axis (accepts 1 ion)
# incomplete
def create_annotated_proj(ds,
                            index,
                            ion,
                            axis = 0,
                            center = [0.5,0.5,0.5],
                            width=(600,'kpc'),
                            output_filename_head = './'):
    try:
        os.makedirs(output_filename_head+"annotatedplots")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    p = yt.ProjectionPlot(ds, axis, ion_p_name(ion)+'_number_density', center=center, width=width)
    p.annotate_ray(ray_list[index], arrow=True)
    p.save(output_filename_head+'annotatedplots/annotatedplot_%s_%s_%s.png' % (ion_p_name(ion),ds,index))

# combines all three plots into a 3x1 grid; all line-of-sight (los) plots
def spectrum_los_plot(ds,
                      num_points,
                      ion_species,
                      sightline,
                      center,
                      in_radius,
                      out_radius,
                      spec_lambda_min = 1000,
                      spec_lambda_max = 1050,
                      spec_dlambda = 0.1,
                      ndpl_y_lim = (1e-12,1e-6),
                      ndpl_log = False,
                        output_filename_head = './'):
    try:
        os.makedirs(output_filename_head+"los")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # check that this is correct
    generate_lightrays(ds,
                       num_points,
                       ion_species,
                       sightline,
                       center,
                       in_radius,
                       out_radius,
                       output_filename_head)
    for index in np.arange(0,num_points):
        for ion in ion_species:
            create_spectrum(ds,index,ion,in_radius,out_radius,spec_lambda_min,spec_lambda_max,spec_dlambda,output_filename_head)
            create_ndpl_plot(ds,index,ion,in_radius,out_radius,ndpl_y_lim,ndpl_log,output_filename_head)
            create_annotated_proj(ds,index,ion,0,(600,'kpc'),output_filename_head)

            # string everything together in a single image
            list_im = [output_filename_head+'annotatedplots/annotatedplot_%s_%s_%s.png' % (ion_p_name(ion),ds,index),
                       output_filename_head+'spectra/spec_%s_r%s-%s_%s_%s.png' % (ion_p_name(ion),yt_remove_spaces(in_radius),
                            yt_remove_spaces(out_radius),ds,index),
                       output_filename_head+'ndpl/ndpl_%s_r%s-%s_%s_%s.png' % (ion_p_name(ion),yt_remove_spaces(in_radius),
                            yt_remove_spaces(out_radius),ds,index)]
            imgs    = [ PIL.Image.open(i) for i in list_im ]
            max_shape = sorted( [(i.size[0], i.size ) for i in imgs])[2][1]
            #print([(i.size[0], i.size ) for i in imgs])
            #print(max_shape)

            # for a vertical stacking it is simple: use vstack
            #print((np.asarray( i.resize(max_shape[0],i.size[1]*max_shape[0]/i.size[0]) ) for i in imgs ) )
            temp_list = []
            another_list = []
            for i in imgs:
                x = (int(round(i.size[1]*max_shape[0]/i.size[0])))
                temperino = i.resize((max_shape[0],x), PIL.Image.ANTIALIAS)
                another_list.append(np.asarray(temperino))

            #print(temp_list,another_list)
            #[int(max_shape[0]/i.size[0] * s) for s in i.size] ), PIL.Image.ANTIALIAS) for i in imgs ]
            imgs_comb = np.vstack(another_list)
            #print(ims_comb)
            imgs_comb = PIL.Image.fromarray( imgs_comb)
            imgs_comb.save( output_filename_head+'los/los_%s_r%s-%s_%s_%s.png' % (ion_p_name(ion),yt_remove_spaces(in_radius),
                yt_remove_spaces(out_radius),ds,index) )
            print("%s,%s,%s,%s,%s ray done" % (ion,in_radius,out_radius,ds,index))

# creates projection plots along x and z axes
def create_axes_proj(ds,
                    ion,
                    axis = 0,
                    center = [0.5,0.5,0.5],
                    width=(600,'kpc'),
                    ylim1=1e10,
                    ylim2=1e24,
                    output_filename_head = './'):
    try:
        os.makedirs(output_filename_head+"xz_plots")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    proj_filename = (output_filename_head+'xz_plots/xz_plot_%s_%s_%s.png') % (ion_p_name(ion),ds,axis)
    if not os.path.exists(proj_filename):
        trident.add_ion_fields(ds, ions=[ion])
        p = yt.ProjectionPlot(ds, axis, ion_p_name(ion)+'_number_density', center=center, width=width)
        p.set_cmap(ion_p_name(ion)+'_number_density','STD GAMMA-II')
        p.set_zlim(ion_p_name(ion)+'_number_density',ylim1,ylim2)
        p.annotate_timestamp(redshift=False,draw_inset_box=True)
        p.save(proj_filename)

# combines xz plots into one image, iterate over datadumps
def combine_xz_plots(filenames,
                    num_procs,
                    ion_species = None,
                    center = None,
                    ylim1x = None,
                    ylim2x = None,
                    output_filename_head = './'):
    if ion_species is None:
        ion_species = ['O VI', 'H I', 'C IV']
    if center is None:
        center = [0.5, 0.5, 0.5]

    try:
        os.makedirs(output_filename_head+"combined_xz_plots")
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    ts = yt.DatasetSeries(filenames, parallel = num_procs)

    for ds in ts.piter():
    # for fn in yt.parallel_objects(['/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0008/DD0008', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0006/DD0006', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0007/DD0007', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0010/DD0010', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0001/DD0001', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0005/DD0005', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0004/DD0004', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0011/DD0011', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0002/DD0002', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0000/DD0000', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0012/DD0012', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0009/DD0009', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0013/DD0013', '/mnt/research/galaxies-REU/sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_tunedCGM/DD0003/DD0003'], num_procs):
    #    ds = yt.load(fn)
        for ion in ion_species:
            # suggestions
            if ylim1x is None and ylim2x is None:
                if ion is 'O VI':
                    ylim1 = 1e13
                    ylim2 = 1e18
                elif ion is 'H I':
                    ylim1 = 1e13
                    ylim2 = 1e24
                elif ion is 'C IV':
                    ylim1 = 1e10
                    ylim2 = 1e18
                elif ion is 'Mg II':
                    ylim1 = 1e10
                    ylim2 = 1e20
                else:
                    ylim1 = ylim1x
                    ylim2 = ylim2x

            xz_plot_filename = (output_filename_head+'combined_xz_plots/combined_xz_plot_%s_%s.png') % (ion_p_name(ion),ds)
            if not os.path.exists(xz_plot_filename):
                create_axes_proj(ds,ion,0,center,(600,'kpc'),ylim1,ylim2,output_filename_head)
                create_axes_proj(ds,ion,2,center,(600,'kpc'),ylim1,ylim2,output_filename_head)
                # print("%s,%s individual xz plots created\n" % (ion,ds))

                # string everything together in a single image
                list_im = [output_filename_head+'xz_plots/xz_plot_%s_%s_0.png' % (ion_p_name(ion),ds),
                           output_filename_head+'xz_plots/xz_plot_%s_%s_2.png' % (ion_p_name(ion),ds)]
                imgs    = [ PIL.Image.open(i) for i in list_im ]
                max_shape = sorted( [(i.size[0], i.size ) for i in imgs])[1][1]
                #print([(i.size[0], i.size ) for i in imgs])
                # print("Max shape is",max_shape,"\n")

                # for a vertical stacking it is simple: use vstack
                #print((np.asarray( i.resize(max_shape[0],i.size[1]*max_shape[0]/i.size[0]) ) for i in imgs ) )
                temp_list = []
                another_list = []
                for i in imgs:
                    x = (int(round(i.size[1]*max_shape[0]/i.size[0])))
                    temperino = i.resize((max_shape[0],x), PIL.Image.ANTIALIAS)
                    another_list.append(np.asarray(temperino))
                # print("append to list %s,%s\n" % (ion,ds))
                #print(temp_list,another_list)
                #[int(max_shape[0]/i.size[0] * s) for s in i.size] ), PIL.Image.ANTIALIAS) for i in imgs ]
                imgs_comb = np.vstack(another_list)
                #print(ims_comb)
                imgs_comb = PIL.Image.fromarray( imgs_comb)
                # print("did I make it this far %s,%s\n" % (ion,ds))
                imgs_comb.save(xz_plot_filename)
                # print("%s,%s xz plots done" % (ion,ds))

                # print("i doubt i made it *this* far %s\n" % (ds))
                if yt.is_root():
                    print("Written:",ds)
                    file = open(output_filename_head+"filenames_xz_plots.fn",'a+')
                    file.write(filenames[:-13]+"%s/%s -- %s\n" % (ds,ds,ion))
                    file.close()

# creates a more general version of cdip using uniformly random sightlines
# 'output_filename_head' can be used to place files inside specific directories (end with a '/')
# center_num_lines is the number of sightlines created within the central region of radius dr
def create_cdip_uniform(filenames,
                        num_procs,
                        ion_species = None,
                        center_coord = None,
                        center_units = None,
                        radius_length = None,
                        radius_units = None,
                        sightline_length = None,
                        sightline_units = None,
                        center_num_lines = None,
                        num_bins = None,
                        y_lim1 = None,
                        output_filename_head = './'):
    if ion_species is None:
        ion_species = ['O VI', 'H I', 'C IV']
    if center_coord is None:
        center_coord = [0.5, 0.5, 0.5]
    if center_units is None:
        center_units = "code_length"
    if radius_length is None:
        radius_length = 300
    if radius_units is None:
        radius_units = "kpc"
    if sightline_length is None:
        sightline_length = 600
    if sightline_units is None:
        sightline_units = "kpc"
    if center_num_lines is None:
        center_num_lines = 1
    if num_bins is None:
        num_bins = 10
    if y_lim1 is None:
        y_lim1 = [0,0]
    if output_filename_head is None:
        output_filename_head ='./'

    # print("----")
    # print(filenames)
    # print(num_procs)
    # print(ion_species)
    # print(center_coord,center_units)
    # print(radius_length,radius_units)
    # print(pix_num, num_bins)
    # print(y_lim1,y_lim2,y_lim3)
    # print(output_filename_head)
    # print("----")

    # set up logging script parameters
    logging.info('Filenames: %s' % filenames)
    logging.info('Number of processes: %s' % num_procs)
    logging.info('Ion species: %s' % ion_species)
    logging.info('Center: %s %s' % (center_coord,center_units))
    logging.info('Sightline Length: %s %s' % (sightline_length,sightline_units))
    logging.info('Number of lines in center: %s' % center_num_lines)
    logging.info('Number of radial bins: %s' % num_bins)
    logging.info('Y-limits, if any set: %s' % (y_lim1))
    logging.info('Output filename head: %s' % output_filename_head)

    print("pit stop 1")
    #filename_tracker = open(output_filename_head+"filenames_cdip_uniform.fn",'a+')
    #filename_tracker.seek(0)

    #filelist = glob.glob(filenames)
    #filelist.sort()
    #tracked_filelist = [line.rstrip('\n') for line in filename_tracker]
    #filename_tracker.close()

    #untracked_filelist = list(set(filelist)-set(tracked_filelist))
    ts = yt.DatasetSeries(filenames, parallel = num_procs)

    for ds in ts.piter():
        trident.add_ion_fields(ds, ions=ion_species)

        center = ds.arr(center_coord,center_units)
        dr = ds.quan(radius_length,radius_units) / num_bins
        sightline = ds.quan(sightline_length,sightline_units)
        d = radius_length*2.
        bins = np.arange(0,num_bins)
        impact_parameters=[]
        in_out=[]
        for bin in bins:
            in_radius = bin * dr
            out_radius = (bin + 1) * dr
            in_out.append([in_radius, out_radius])
            # num_lines = center_num_lines * (3 * bin**2 + 3 * bin + 1)
            num_lines = center_num_lines * ((bin + 1)**2 - bin**2)

            if num_lines < 10:
                num_lines = 10

            impact_param = generate_lightrays(ds,
                        num_lines,
                        ion_species,
                        sightline,
                        center,
                        in_radius,
                        out_radius,
                        output_filename_head)

            impact_parameters.append(impact_param)

            print(output_filename_head,ds,"bin",bin,"finished")

        impact_parameters = np.array(impact_parameters)
        in_out = np.array(in_out)
        np.save("in_out_radii", in_out)
        np.save("impact_param", impact_parameters)
        # iterate over each ion
        for ion in ion_species:
            cdip_filename = (output_filename_head+'cdip_uniform_'+ion_p_name(ion)+'_%s_'+str(int(d))+radius_units+'_box.png') % ds
            if not os.path.exists(cdip_filename):

                r = np.arange(0,num_bins) #* x axis in bins
                y_med_ion = []
                y_mean_ion = []
                y_std_ion = []
                print("am i here yet")
                #make list of lists by bin
                val_by_r_ion = [[] for _ in r]

                # iterate over each bin

                for bin in bins:
                    in_radius = bin * dr
                    out_radius = (bin + 1) * dr
                    # num_lines = center_num_lines * (3 * bin**2 + 3 * bin + 1)
                    num_lines = center_num_lines * ((bin + 1)**2 - bin**2)

                    if num_lines < 10:
                        num_lines = 10

                    for index in np.arange(0,num_lines):
                        val_by_r_ion[bin].append(integrate_ndpl(ds,
                                                                index,
                                                                ion,
                                                                in_radius,
                                                                out_radius,
                                                                output_filename_head))
                    print("another point")
                    #reduce val_by_r to median
                    y_med_ion.append(np.median(val_by_r_ion[bin]))
                    y_mean_ion.append(np.mean(val_by_r_ion[bin]))
                    y_std_ion.append(np.std(val_by_r_ion[bin]))

                ###
                # FORMAT FINAL IMAGE
                ###

                plt.close('all')

                fig = plt.figure()
                # fig.set_size_inches(6, 12)
                # ax = fig.add_subplot(111)
                plt.yscale('log')

                bar_width = radius_length / num_bins
                plt.bar((r+0.5) * bar_width, y_med_ion, bar_width, color='0.4')
                plt.errorbar((r+0.5) * bar_width, y_mean_ion, yerr=y_std_ion, fmt="o", color='0.0', capsize=bar_width/4)
                plt.xlim([0,radius_length])
                if y_lim1 != [0,0]:
                    plt.ylim(np.float_(y_lim1))

                if ion == 'O VI':
                    plt.ylim(np.float_([1e13,1e17]))
                if ion == 'H I':
                    plt.ylim(np.float_([1e13,1e22]))
                if ion == 'C IV':
                    plt.ylim(np.float_([1e11,1e17]))

                plt.title(ion+': Time Step %s' % ds, y=1.03)
                plt.xlabel('Impact Parameter (%s)' % (radius_units))
                plt.ylabel(r'Column Density (cm$^{-2}$)')
                plt.savefig(cdip_filename)
                print("would be nice to get this far")
                logging.debug('Data dump %s, ion %s complete' % (ds,ion))
                if yt.is_root():
                    file = open(output_filename_head+"filenames_cdip_uniform.fn",'a+')
                    file.write(filenames[:-13]+"%s/%s -- %s\n" % (ds,ds,ion))
                    file.close()


if __name__ == "__main__":
    # test parameters
    # not sure this one has correct syntax though
    create_cdip_xz_proj("../../sims/isolated-galaxies/MW_1638kpcBox_800pcCGM_200pcDisk_lowres/DD????/DD????",
                         2, ['O VI', 'H I', 'C IV'], [0.5411343,  0.4509816,  0.5134753], "code_length",
                         300, "kpc", pix_num = 1000, num_bins = 10, y_lim1 = [0,0], y_lim2 = [0,0], y_lim3 = [0,0],
                         output_filename_head='../img_dir/test_directory/.//')
