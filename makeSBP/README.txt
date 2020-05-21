OK, all of the useful scripts are living in this directory: /projects/b1095/kylekremer/python_code/CMC_Grid_March2019

Here is summary of the important things and what does what:

make_SBP_and_vel_profile_one.py:  creates 2D projections, cluster params files, SBPs, velocity dispersion profiles for all snapshots for a given path. By adjusting the Delta parameter, you can adjust how many snapshots you create these for. Currently I have Delta=-5000, which means I only create for the last snapshot.

plot_SBP_sigma.py: plots the SBP and velocity dispersion profile for a given model and snapshot compared to a given observed cluster of interest. So for example, to plot the N=3.5e6 model we were looking at compared to 47 Tuc, run:

>>>python plot_SBP_sigma.py /projects/b1095/syr904/cmc/cmc-mpi-09/rundir/3.5e6rv1fb10kick1.0Z0.0002/ 0656 ngc104

find_best_fit_SBP_sigma.py: finds the best fit model for a given observed cluster from the complete CMC catalog list. To run use the command:

>>>python find_best_fit_SBP_sigma.py ngc104 1

this will identify best fit for 47 Tuc. The 1 at the end tells it to find the single best fit. If you enter 2 instead, it will plot the two best fit models, 3 will plot three models, etc. This script also outputs the best fit model mass, rc, rh, etc for 47 Tuc compared to observations from Harris.


RUN WITH PYTHON2!!!
