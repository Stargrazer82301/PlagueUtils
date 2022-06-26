# PlagueUtils
Python code designed to generate a very rough, *very* non-scientific plot of COVID-19 (SARS-CoV-2) viral load versus time, using cropped photographs of lateral flow tests (AKA rapid antigen tests).  The code is run using a single function, which has only one required argument, called in the following manner:  

`PlagueUtils.Run('/some/path/to/photo/folder/')`

The function should be given the path to a directory that contains photographs of lateral flow tests, cropped to only show the test strip itself, oriented such tha that the control line is the uppermost of the two lines (with the line itself therefore being horizontal). See the photos provided in this repo for illustration.  

The individual photos should be `.jpg` files (not `.jpeg` files), and have file names of the format `NAME_YYYY-MM-DD_HH-MM.jpg`, where `NAME` is the same for all the files.  

There are optional arguments to change which colour channel is used in the input jpgs (by default, green only is used), and to run in debug mode (that outputs a some intermediate plots)  

An example directory called `CJRC` contains a set of example images, chronicling my own covid experience. You can therefore generate an example of the output plots by running the command:  

`PlagueUtils.Run('CJRC/')  `
 
The code uses the C (control) line as a "calibration source"; ie, it measures the strength of the T (test) line relative to the control line to estimate the viral load measured by a given test. This should help to account for differences between tests. It also uses crude implementations of various techniques used in analysis of astronomical spectra to, eg, account for the background "continuum" of the test strip.  

I wrote this code as an excuse to learn specutils, an astropy-afolicated Python package for analysing astronomical spectroscopy data. Any horrific code found herein should be ascribed to the fact I wrote it whilst I had the plague.
