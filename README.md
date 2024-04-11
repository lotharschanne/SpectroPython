# SpectroPython
Simple scripts written in Python for processing and evaluating single or series of optical astronomical 1d spectra in fits or ascii format.

Features

For series of 1d spectra (fits or ascii format):
  Clipping to one wavelength range
  Plotting of all spectra with an offset with selectable wavelength range and selectable observation time range
  Color-coded dynamic plot based on the observation date
  Plotting a spectrum and marking lines of an selectable element (ion)
  Header display and modification
  Conversion of fits 1d spectra into ascii files (csv, tab)
  Printout of all observation times of a series
  Smoothing of a spectra series
  Barycentric correction of a spectrum series
  Determination of the signal-to-noise ratio
  Radial velocity measurements by determining a line minimum (maximum) using various methods (regression, spline, Gaussian fit, RBF, interactively)
  Correction of a spectrum series by a radial velocity, which is determined by comparison with a template or measurement of terrestrial lines or a radial velocity to be entered manually
  Graphical representation of a line in velocity space
  New binning to a selectable step size and calculation of an average spectrum as well as calculation of difference spectra using a template
  1 point Normalization of a spectrum series
  Automatic normalization of spectrum series to the continuum
  Line fitting according to various models (Gauss, double gauss, voigt, double voigt)
  Determination of the equivalent width, FWHM and line depth of a line
  Determination of the equivalent width graphically interactive, for different wavelength ranges of a line
  Cross correlation of a spectrum series with a template. Complete spectrum or a selectable wavelength section

Theoretical spectra
  Cropping, rebinning and convolving of a spectrum in ascii format to a certain wavelength range, a certain step size and a certain resolution

Binaries
  Determination of the period from a list of observation times (JD) and radial velocities
  Calculation of the phases for a list of observation times (JD) and radial velocities (RV) with specified To and period
  Calculation of the phases for a list of spectra with specified To and period
  Calculation of a phase plot of a list of observation times (JD) and radial velocities
  Calculation of radial velocities from the orbital parameters at given observation times
  Calculation of orbital elements from a list of observation times (JD) and radial velocities
  Calculation of theoretical sum spectra of multiple systems and formation of the differences to measured spectra to determine effects (e.g. emissions) that are not contained in the theoretical spectra.
