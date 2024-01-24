# ElasticBlister
Python code for modelling the behaviour of a 2D elastic blister with a viscous interior placed on a sloping base.

This is a cleaned up version of the code used in Tobin, S., & Neufeld, J. 2023 Evolution of an elastic blister in the presence of sloping topography. *J. Fluid Mech.* **967**, A5.

## How to run

The numerical simulation is contained in `blister.py`. This produces output files containing the location of key points of the blister at pre-defined times, as well as less frequent full profiles of the blister. The files `frontplot.py` and `plotprofiles.py` are examples of code used to plot this output. In all files, the depth of the pre-wetted layer is taken to be $\delta = 10^{-3}$ by default, and output is saved in the current working directory. To run the code without changing this, first run `blister.py`, and then either of the plotting programmes.

> [!NOTE]
> Currently, `blister.py` does not terminate until $t = 1,000,000$. You will probably want to manually terminate the code before this happens.

## Notes

1. The code uses a finer grid spacing at the downslope edge of the blister, but not at the upslope edge. This is fine when you are modelling the movement of the blister down a slope, but means it is not suited to horizontally spreading blisters where you would like well resolved upslope peeling.
2. It is slow. The original version of this code was written in FORTRAN, but I switched to Python when I realised that I would have to regrid regularly and it would be nice to be able to just add extra elements to arrays. The trade off between time spent writing the code and time spent running the code worked for me, but if you are starting from scratch you may be best using this as a reference and rewriting it to be more efficient.
