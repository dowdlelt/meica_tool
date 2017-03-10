# meica_tool

I've created a handy matlab script that works with meica.py (https://bitbucket.org/prantikk/me-ica) v3 - from the experimental branch. 

It creates a series of figures that are useful for visualizing the output in a quick manner, including component timeseries as determined but by meica.py, color coded on whether they were:

* BOLD-like - green 
* Non-BOLD - Red
* r2 weighted - pink
* Ignored - black. 

Each plot includes brain slices of the component beta values (from TED/betas_OC.nii)

* motion parameters and framewise displacement
* kappa vs rho scatter plot, where size is proportaional to variance, colors as above
* kappa vs rho line plot
* Bar plot of variance explained

It then creates a bar plot showing the relative variance of each of those categories. 

Its (still) ugly code, but effective...for now. 
