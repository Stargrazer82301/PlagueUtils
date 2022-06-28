# PlagueUtils

Python code designed to generate a very rough, *very* non-scientific plot of COVID-19 (SARS-CoV-2) viral load versus time, using cropped photographs of lateral flow tests (AKA rapid antigen tests).  

![Example of a PlagueUtils output plot](https://raw.githubusercontent.com/Stargrazer82301/PlagueUtils/main/Readme_Images/CJRC_Integrated_vs_Time.png)

## Usage

The code is run using a single function, which has only one required argument, called in the following manner:  

`PlagueUtils.Run('/some/path/to/photo/folder/')`

The function should be given the path to a directory that contains photographs of lateral flow tests, cropped to only show the test strip itself, oriented such tha that the control line is the uppermost of the two lines (with the line itself therefore being horizontal). See the photos provided in this repo for illustration.  

The individual photos should be `.jpg` files (not `.jpeg` files), and have file names of the format `NAME_YYYY-MM-DD_HH-MM.jpg`, where `NAME` is the same for all the files.  

There are optional arguments to change which colour channel is used in the input jpgs (by default, green only is used), and to run in debug mode (that outputs a some intermediate plots)  

An example directory called `CJRC` contains a set of example images, chronicling my own covid experience. You can therefore generate an example of the output plots by running the command:  

`PlagueUtils.Run('CJRC/')  `
 
The code uses the C (control) line as a "calibration source"; ie, it measures the strength of the T (test) line relative to the control line to estimate the viral load measured by a given test. This should help to account for differences between tests. It also uses crude implementations of various techniques used in analysis of astronomical spectra to, eg, account for the background "continuum" of the test strip. 

## Example Outputs

PlagueUtils outputs 4 plots (into the working directory from which you are running PlagueUtils). The `Spectra_Raw` plot shows the "spectrum" of all of tests provided:

![Example of Plague
Utils raw spectra plot](https://raw.githubusercontent.com/Stargrazer82301/PlagueUtils/main/Readme_Images/CJRC_Spectra_Raw.png)

The `Spectra_Pro` plot shows the processed versions of the spectra, with thebaseline "continuum" subtracted:

![Example of PlagueUtils processed spectra plot](https://raw.githubusercontent.com/Stargrazer82301/PlagueUtils/main/Readme_Images/CJRC_Spectra_Pro.png)

Two plots of plagueiness ("viral load") against time out produced. The first one shows the *peak* of each T line, relative to the peak of the corresponding C lines, versus time:  

![Example of PlagueUtils plot of peak plagueiness versus time](https://raw.githubusercontent.com/Stargrazer82301/PlagueUtils/main/Readme_Images/CJRC_Peak_vs_Time.png)

The second time series plot shows the *integrated* strength of each T line, again calibrated relative to the integrated strengths of the corresponding C lines:

![Example of PlagueUtils plot of peak plagueiness versus time](https://raw.githubusercontent.com/Stargrazer82301/PlagueUtils/main/Readme_Images/CJRC_Integrated_vs_Time.png)

I consider the integrated plot to be the primary output. Although I'm pleasantly surprised how well mine seems to agree with the peak plot. Let me know if yours don't agree!

## Final Words

I wrote this code as an excuse to learn specutils, an astropy-affiliated Python package for analysing astronomical spectroscopy data. Any horrific code found herein should be ascribed to the fact I wrote it whilst I had the plague.

Thank you to the Hagens, for not kicking me out when I tested positive 36 hours after arriving at their house, and taking such good care of me. And thank you to Em for still coming out to California, to hang out out intermediate distance for a frustrating long amount of time.
