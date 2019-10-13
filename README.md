# SingExoplanets

=====
OptimalApertureSize determines the highest average photon count among the rows, which indicates where the horizontal aperture begins.

Ppm_vs_aperture_fit plots the aperture range against the ppm for that range.
SubpixelCode uses a set aperture position and range and varies subpixels.
=====
MaxPhotonCountandPosition does ??

OptimalApertureCurve creates a best fit line for the aperture for each image by determining the max photon count per column and position of that photon count, then applying a fit to that data.
=====


First, use OptimalApertureSize to find the starting row for the horizontal aperture. Then use that value and run transit_fit_wasp39.py to measure the photon count for the chosen aperture and aperture range. Thirdly, use Ppm_vs_aperture_fit to visualize the ppm for corresponding aperture ranges. After determining the lowest ppms, select the three corresponding apertures and run SubpixelCode to look for even lower ppms.


The improved version of this process should find the max photon count for each column, then create a best fit for these points. Then, several aperture ranges should be measured. The lowest three will then be measured with several subpixel aperture sizes. Finally, these points could be visualized on a plot of aperture +- vs ppm.
