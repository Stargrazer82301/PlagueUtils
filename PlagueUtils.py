#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri 24th June 2022
@author: Chris Clark

This code is designed to generate a very rough, *very* non-scientific plot of
COVID-19 (SARS-CoV-2) viral load versus time, using cropped photographs of
lateral flow tests (AKA rapid antigen tests).

The code is run using a single function, which has only one required argument,
called in the following manner:

PlagueUtils.Run('/some/path/to/photo/folder/')

The function should be given the path to a directory that contains photographs
of lateral flow tests, cropped to only show the test strip itself, oriented
such tha that the control line is the uppermost of the two lines (with the
line itself therefore being horizontal). See the photos provided in this repo
for illustration. It also requires that the lines not fall in the top or
bottom 10% of the image.

The individual photos should be .jpg files (not .jpeg files), and have file
names of the format NAME_YYYY-MM-DD_HH-MM.jpg, where NAME is the same for
all the files.

There are optional arguments to change which colour channel is used in the
input jpgs (by default, green only is used), and to run in debug mode (that
outputs a some intermediate plots)

An example directory called CJRC contains a set of example images, chronicling
my own covid experience. You can therefore generate an example of the output
plots by running the command:

PlagueUtils.Run('CJRC/', green_only=True,  debug=True)

The code uses the C (control) line as a "calibration source"; ie, it measures
the strength of the T (test) line relative to the control line to estimate the
viral load measured by a given test. This should help to account for
differences between tests. It also uses crude implementations of various
techniques used in analysis of astronomical spectra to, eg, account for the
background "continuum" of the test strip.

I wrote this code as an excuse to learn specutils, an astropy-afolicated
Python package for analysing astronomical spectroscopy data. Any horrific
code found herein should be ascribed to the fact I wrote it whilst I had the
plague.
"""




# Import smorgasbord
import os
import copy
import numpy as np
import astropy.units as u
import astropy.modeling.models
import astropy.modeling.fitting
import scipy.ndimage
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import specutils
import specutils.fitting
import time
import datetime
import pytz
import PIL
import warnings
warnings.filterwarnings("ignore")

# Disable specutils continuum warnings
specutils.conf.do_continuum_function_check = False

# Set timezone
timezone = pytz.utc

# Enable seaborn for easy, attractive plots
plt.close('all')
plt.ioff()
sns.set(context='talk') # Possible context settings are 'notebook' (default), 'paper', 'talk', and 'poster'
sns.set_style('darkgrid', {'font.sans-serif':'DejaVu Sans'})
matplotlib.rcParams['font.sans-serif'] = 'DejaVu Sans'
matplotlib.rcParams['mathtext.default'] = 'regular'





def Run(test_dir,
        green_only=True,
        debug=False):
        """
        Function that generates a very rough, *very* non-scientific plot of COVID-19 (SARS-CoV-2) viral load versus
        time, using cropped photographs of lateral flow tests (AKA rapid antigen tests).

        Arguments:
            test_dir:
                    A string, giving the path to a directory that contains photographs of lateral flow tests,
                    oen photo per test, cropped to only show the test strip itself, oriented such tha that the
                    control line is the uppermost of the two lines (with the line itself therefore being
                    horizontal).

        Keyword arguments:
            green_only:
                    A boolean, stating whether to only use the green channel. Default is true. The gold
                    nanoparticles in C and T lines absorb mainly green light, some blue, and little red. So the line
                    has the most contrast with the white background in the green channel; adding other channels
                    usually worsens the signal-to-noise ratio.
            debug:
                    A boolean, stating whether to output intermediate plots. Default is false. This is mainly for
                    diagnosing if something has gone wrong.
        """

        # Identify test images in data directory
        test_files = [test_file for test_file in os.listdir(test_dir) if '.jpg' in test_file]
        test_files.sort()
        n_tests = len(test_files)
        if n_tests == 0:
            raise Exception('All input image filenames must be of format \'NAME_YYYY-MM-DD_HH-MM.jpg\'')

        # Extract name from file names
        names = [test_file.split('_')[0] for test_file in test_files]
        if len(list(set(names))) > 1:
            raise Exception('All input image filenames must be of format \'NAME_YYYY-MM-DD_HH-MM.jpg\'')
        else:
            name = str(names[0])

        # Create lists to hold data
        test_spec_list = []
        raw_spec_list = []
        pro_spec_list = []
        pro_calib_list = []
        test_datetime_list = []
        c_flux_peak_list = []
        c_flux_peak_unc_list = []
        c_flux_int_list = []
        c_flux_int_unc_list = []
        t_flux_peak_list = []
        t_flux_peak_unc_list = []
        t_flux_int_list = []
        t_flux_int_unc_list = []



        # Loop over test images
        for i in range(n_tests):

            # Extract date and time from file name
            try:
                test_year = int(test_files[i].split('_')[1].split('-')[0])
                test_month = int(test_files[i].split('_')[1].split('-')[1])
                test_day = int(test_files[i].split('_')[1].split('-')[2])
                test_hour = int(test_files[i].split('_')[2].split('-')[0])
                test_min = int(test_files[i].split('_')[2].split('-')[1].replace('.jpg',''))
                test_datetime = datetime.datetime(test_year,
                                                  test_month,
                                                  test_day,
                                                  test_hour,
                                                  test_min,
                                                  0,
                                                  tzinfo=timezone)
                test_datetime_list.append(test_datetime)
            except:
                raise Exception('All input image filenames must be of format \'NAME_YYYY-MM-DD_HH-MM.jpg\'')



            # Load current image, and convert to numpy array
            if green_only:
                test_img = PIL.Image.open(os.path.join(test_dir,test_files[i])).split()[1]
            else:
                test_img = PIL.Image.open(os.path.join(test_dir,test_files[i])).convert('L')
            test_img = np.asarray(test_img)

            # Take mean to collapse image down to one dimension
            test_img = np.mean(test_img, axis=1)

            # Invert image, to turn our "emission" feature into an "absorption" feature
            test_img = -1.0 * test_img

            # Use sigma-clipping to set the baseline to zero
            test_clip = SigmaClip(test_img, median=True)
            test_img -= test_clip[1]



            # Create a specutils object to store rapid test "spectrum", with "wavelength" spanning 100 parsecs (YOLO)
            test_wav = np.linspace(0, 100, num=len(test_img))
            test_spec = specutils.Spectrum1D(flux=test_img*u.pix,
                                             spectral_axis=test_wav*u.pc)

            # Plot this tests raw spectrum
            if debug:
                fig, ax = plt.subplots(1, 1, figsize=(6,4))
                ax.step(test_spec.spectral_axis, test_spec.flux)
                ax.set_xlabel('Position', fontname='sans')
                ax.set_ylabel('Strength', fontname='sans')
                ax.set_xlim(0, 100)
                plot_raw_filename = test_files[i].replace('.jpg','_Plot_Raw.png')
                fig.savefig(os.path.join(test_dir, plot_raw_filename),
                            dpi=175, bbox_inches='tight')



            # Identify all regions that are not over 1-sigma; we'll assume these trace reliable baseline
            base_spec = copy.deepcopy(test_spec)
            base_thresh = 0.50 * test_clip[0] * u.pix
            base_good = np.where((np.abs(base_spec.flux) < base_thresh) \
                                 | (test_spec.spectral_axis < (15*u.pc)) \
                                 | (test_spec.spectral_axis > (85*u.pc)))

            # Compute spline over non-baseline pixels
            base_spline = scipy.interpolate.UnivariateSpline(np.array(base_spec.spectral_axis.data)[base_good],
                                                             np.array(base_spec.flux.data)[base_good],
                                                             k=5)
            base_spec.flux.data[:] = base_spline(np.array(base_spec.spectral_axis.data)) * u.pix

            # Plot this tests raw spectrum
            if debug:
                fig, ax = plt.subplots(1, 1, figsize=(6,4))
                ax.scatter(base_spec.spectral_axis, base_spec.flux, s=1)
                ax.set_xlabel('Position', fontname='sans')
                ax.set_ylabel('Strength', fontname='sans')
                ax.set_xlim(0, 100)
                plot_raw_filename = test_files[i].replace('.jpg','_Plot_Baseline.png')
                fig.savefig(os.path.join(test_dir, plot_raw_filename),
                            dpi=175, bbox_inches='tight')



            # Baseline-subtract our actual test spectra (saving copy of raw spectrum)
            raw_spec = copy.deepcopy(test_spec)
            test_spec = test_spec - base_spec

            # Calculate baseline noise level
            test_unc = np.std(test_spec.flux[base_good])

            # Plot this test's processed spectrum
            if debug:
                fig, ax = plt.subplots(1, 1, figsize=(6,4))
                ax.step(test_spec.spectral_axis, test_spec.flux)
                ax.set_xlabel('Position', fontname='sans')
                ax.set_ylabel('Strength', fontname='sans')
                ax.set_xlim(0, 100)
                plot_raw_filename = test_files[i].replace('.jpg','_Plot_Pro.png')
                fig.savefig(os.path.join(test_dir, plot_raw_filename),
                            dpi=175, bbox_inches='tight')



            # Create a smoothed spectra to use for line extraction
            filt_spec = specutils.manipulation.gaussian_smooth(test_spec, test_img.shape[0]*0.01)
            filt_reg = specutils.SpectralRegion(20*u.pc, 85*u.pc)
            filt_spec = specutils.manipulation.extract_region(filt_spec, filt_reg)

            # Extract lines from filtered spectrum (computing fresh sigma clip on baseline-subtracted data)
            test_clip = SigmaClip(test_img, median=True, sigma_thresh=2.0)
            filt_thresh = 0.05 * test_clip[0]
            filt_lines = specutils.fitting.find_lines_derivative(filt_spec, flux_threshold=filt_thresh)

            # Only keep 2 strongest "emission" lines
            filt_lines = filt_lines[filt_lines['line_type'] == 'emission']
            filt_lines_peaks = [filt_spec[filt_lines['line_center_index'][j]].flux.value[0] for j in range(len(filt_lines))]
            filt_lines_peaks_sorted = copy.deepcopy(filt_lines_peaks)
            filt_lines_peaks_sorted.sort()
            filt_lines = filt_lines[np.where(filt_lines_peaks >= filt_lines_peaks_sorted[-2])[0]]



            # Extract region containing C line, using midpoint between C and T as bisectior
            bisect_wav = np.min(filt_lines['line_center']) + (0.5 * (np.max(filt_lines['line_center']) - np.min(filt_lines['line_center'])))
            c_proreg = specutils.SpectralRegion(0*u.pc, bisect_wav)
            c_prospec = specutils.manipulation.extract_region(test_spec, c_proreg)

            # Measure peak "flux" in C line region, and estimate uncertinty
            c_flux_peak = c_prospec.max()
            c_flux_peak_list.append(c_flux_peak.value)
            c_flux_peak_unc = 10.0 * test_unc
            c_flux_peak_unc_list.append(c_flux_peak_unc.value)

            # Fit Gaussian to C line, to roughly estimate line width,
            c_model = astropy.modeling.models.Lorentz1D(amplitude=c_flux_peak,
                                                        x_0=np.min(filt_lines['line_center']),
                                                        fwhm=2*u.pc)
            c_line = specutils.fitting.fit_lines(c_prospec, c_model)

            # Now define more specific region around C line, based on line width
            c_proreg = specutils.SpectralRegion(c_line.x_0-(2.0*c_line.fwhm),
                                                c_line.x_0+(2.0*c_line.fwhm))

            # Measure integrated flux of C line, and estimate uncertinty
            c_flux_int = specutils.analysis.line_flux(test_spec, regions=c_proreg)
            c_flux_int_list.append(c_flux_int.value)
            c_flux_int_unc = test_unc * (len(c_prospec.flux.data))**0.5
            c_flux_int_unc_list.append(c_flux_int_unc.value)



            # Extract region containing T line, moving bisector if needed so it's not inside the C line
            try:
                bisect_wav = max(c_line.x_0+(5.0*c_line.fwhm), bisect_wav)
                t_proreg = specutils.SpectralRegion(bisect_wav, 100*u.pc)
                t_prospec = specutils.manipulation.extract_region(test_spec, t_proreg)

                # Measure peak "flux" in T line region, and estimate uncertinty
                t_flux_peak = t_prospec.max()
                t_flux_peak_list.append(t_flux_peak.value)
                t_flux_peak_unc = 10.0 * test_unc
                t_flux_peak_unc_list.append(t_flux_peak_unc.value)

                # Fit Gaussian to T line, to roughly estimate line width,
                t_model = astropy.modeling.models.Lorentz1D(amplitude=t_flux_peak,
                                                            x_0=np.max(filt_lines['line_center']),
                                                            fwhm=2*u.pc)
                t_line = specutils.fitting.fit_lines(t_prospec, t_model)

                # Now define more specific region around T line, based on line width
                t_line_with_min = max(2*u.pc, 2.0*t_line.fwhm)
                t_proreg = specutils.SpectralRegion(t_line.x_0-t_line_with_min,
                                                    t_line.x_0+t_line_with_min)

                # Measure integrated flux of T line, and estimate uncertinty
                t_flux_int = specutils.analysis.line_flux(test_spec, regions=t_proreg)
                t_flux_int_list.append(t_flux_int.value)
                t_flux_int_unc = test_unc * (len(t_prospec.flux.data))**0.5
                t_flux_int_unc_list.append(t_flux_int_unc.value)



                # Calibrate spectrum
                pro_spec = copy.deepcopy(test_spec)
                pro_calib = c_flux_peak.value
                pro_spec /= pro_calib
                pro_calib_list.append(pro_calib)

                # Put spectra in list
                test_spec_list.append(test_spec)
                raw_spec_list.append(raw_spec)
                pro_spec_list.append(pro_spec)
            except:
                breakpoint()



        # Plot raw spectra
        fig_raw_dims = (int(np.floor(np.sqrt(n_tests)+0.5)), int(np.ceil(np.sqrt(n_tests))))
        fig_raw, axes_raw = plt.subplots(fig_raw_dims[0], fig_raw_dims[1],
                                         figsize=(fig_raw_dims[0]*5.0, fig_raw_dims[0]*3.5),
                                         constrained_layout=True)
        raw_ypeak = np.concatenate((t_flux_peak_list, c_flux_peak_list)).max()
        for i in range(n_tests):
            raw_spec = raw_spec_list[i]
            axes_raw.flatten()[i].step(raw_spec.spectral_axis, raw_spec.flux)
            if i >= (0.5 * n_tests):
                axes_raw.flatten()[i].set_xlabel('Position', fontname='sans')
            if i % fig_raw_dims[1] == 0:
                axes_raw.flatten()[i].set_ylabel('Strength', fontname='sans')
            axes_raw.flatten()[i].set_xlim(0, 100)
            axes_raw.flatten()[i].set_ylim(-0.15 * raw_ypeak, 1.2 * raw_ypeak)
            axes_raw.flatten()[i].text(0.025,
                             0.9,
                             time.strftime('%d %b %Y, %H:%M', test_datetime_list[i].timetuple()),
                             transform=axes_raw.flatten()[i].transAxes)
        for i in range(len(axes_raw.flatten())):
            if i > (n_tests - 1):
                axes_raw.flatten()[i].remove()
        fig_raw.savefig(os.path.join(test_dir, name+'_Spectra_Raw.png'),
                        dpi=175, bbox_inches='tight')



        # Plot processed spectra
        fig_pro_dims = (int(np.floor(np.sqrt(n_tests)+0.5)), int(np.ceil(np.sqrt(n_tests))))
        fig_pro, axes_pro = plt.subplots(fig_pro_dims[0], fig_pro_dims[1],
                                         figsize=(fig_pro_dims[0]*5.0, fig_pro_dims[0]*3.5),
                                         constrained_layout=True)
        pro_ypeak = np.max(np.concatenate((t_flux_peak_list, c_flux_peak_list)) / \
                           np.concatenate((pro_calib_list, pro_calib_list)))
        for i in range(n_tests):
            pro_spec = pro_spec_list[i]
            axes_pro.flatten()[i].step(pro_spec.spectral_axis, pro_spec.flux)
            if i >= (0.5 * n_tests):
                axes_pro.flatten()[i].set_xlabel('Position', fontname='sans')
            if i % fig_pro_dims[1] == 0:
                axes_pro.flatten()[i].set_ylabel('C-Normalised Strength', fontname='sans')
            axes_pro.flatten()[i].set_xlim(0, 100)
            axes_pro.flatten()[i].set_ylim(-0.15 * pro_ypeak, 1.2 * pro_ypeak)
            axes_pro.flatten()[i].text(0.025,
                             0.9,
                             time.strftime('%d %b %Y, %H:%M', test_datetime_list[i].timetuple()),
                             transform=axes_pro.flatten()[i].transAxes)
        for i in range(len(axes_raw.flatten())):
            if i > (n_tests - 1):
                axes_pro.flatten()[i].remove()
        fig_pro.savefig(os.path.join(test_dir, name+'_Spectra_Pro.png'),
                        dpi=175, bbox_inches='tight')



        # Turn flux lists into arrays
        t_fluxes_int = np.array(t_flux_int_list)
        t_fluxes_int_unc = np.array(t_flux_int_unc_list)
        c_fluxes_int = np.array(c_flux_int_list)
        t_fluxes_peak = np.array(t_flux_peak_list)
        t_fluxes_peak_unc = np.array(t_flux_peak_unc_list)
        c_fluxes_peak = np.array(c_flux_peak_list)

        # Create array of plottable times
        datetime_array = np.array([test_datetime_list[i].astimezone(timezone) for i in range(len(test_datetime_list))])



        # Calibrate peak T fluxes and uncertainties using C values
        t_fluxes_peak_calib = t_fluxes_peak / c_fluxes_peak
        t_fluxes_peak_unc_calib = t_fluxes_peak_unc / c_fluxes_peak

        # Plot peak strength versus time
        fig, ax = plt.subplots(1, 1, figsize=(8,6), constrained_layout=True)
        ax.plot(datetime_array, t_fluxes_peak/c_fluxes_peak, lw=4)
        ax.fill_between(datetime_array,
                        t_fluxes_peak_calib - t_fluxes_peak_unc_calib,
                        t_fluxes_peak_calib + t_fluxes_peak_unc_calib,
                        alpha=0.15)
        ax.set_xlabel('Time & Date', fontname='sans')
        ax.set_ylabel('C-Calibrated Peak Plague Level', fontname='sans')

        # Format date & time tick labels
        xlocator = matplotlib.dates.AutoDateLocator(minticks=6, maxticks=20)
        xformatter = matplotlib.dates.ConciseDateFormatter(xlocator, tz=timezone)
        ax.xaxis.set_major_locator(xlocator)
        ax.xaxis.set_major_formatter(xformatter)

        # Save plot
        fig.savefig(os.path.join(test_dir, name+'_Peak_vs_Time.png'),
                        dpi=175, bbox_inches='tight')



        # Calibrate integrated T fluxes and uncertainties using C values
        t_fluxes_int_calib = t_fluxes_int / c_fluxes_int
        t_fluxes_int_unc_calib = t_fluxes_int_unc / c_fluxes_int
        t_fluxes_int_unc_calib = (t_fluxes_int_unc_calib**2.0 + 0.05**2.0)**0.5

        # Plot integrated strength versus time
        fig, ax = plt.subplots(1, 1, figsize=(8,6), constrained_layout=True)
        ax.plot(datetime_array, t_fluxes_int/c_fluxes_int, lw=4)
        ax.fill_between(datetime_array,
                        t_fluxes_int_calib - t_fluxes_int_unc_calib,
                        t_fluxes_int_calib + t_fluxes_int_unc_calib,
                        alpha=0.15)
        ax.set_xlabel('Time & Date', fontname='sans')
        ax.set_ylabel('C-Calibrated Integrated Plague Level', fontname='sans')

        # Format date & time tick labels for CO2 axis
        xlocator = matplotlib.dates.AutoDateLocator(minticks=6, maxticks=20)
        xformatter = matplotlib.dates.ConciseDateFormatter(xlocator, tz=timezone)
        ax.xaxis.set_major_locator(xlocator)
        ax.xaxis.set_major_formatter(xformatter)

        # Save plot
        fig.savefig(os.path.join(test_dir, name+'_Integrated_vs_Time.png'),
                        dpi=175, bbox_inches='tight')

        # Report completion and return
        print('Processing complete - get well soon!')
        breakpoint()
        return





# Define sigma-clipping function
def SigmaClip(values, tolerance=0.001, median=False, sigma_thresh=3.0,):

    # Remove NaNs from input values
    values = np.array(values)
    values = values[np.where(np.isnan(values)==False)]

    # Continue loop until result converges
    diff = 10E10
    while diff>tolerance:

        # Assess current input iteration
        if median == False:
            average = np.mean(values)
        elif median == True:
            average = np.median(values)
        sigma_old = np.std(values)

        # Mask those pixels that lie more than 3 stdev away from mean
        check = np.zeros([len(values)])
        check[np.where(values > (average + (sigma_thresh * sigma_old)))] = 1
        check[np.where(values < (average - (sigma_thresh * sigma_old)))] = 1
        values = values[ np.where(check<1) ]

        # Re-measure sigma and test for convergence
        sigma_new = np.std(values)
        diff = abs(sigma_old-sigma_new) / sigma_old

    # Perform final mask
    check = np.zeros([len(values)])
    check[np.where(values > (average + (sigma_thresh * sigma_old)))] = 1
    check[np.where(values < (average - (sigma_thresh * sigma_old)))] = 1
    values = values[np.where(check<1)]

    # Return results
    return [sigma_new, average, values]





# Example use
Run('CJRC/', green_only=True,  debug=False)