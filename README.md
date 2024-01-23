# ElasticBlister
Python code for modelling the behaviour of a 2D elastic blister with a viscous interior placed on a sloping base.

This is a cleaned up version of the code used in Tobin, S., & Neufeld, J. 2023 Evolution of an elastic blister in the presence of sloping topography. *J. Fluid Mech.* **967**, A5.

## Notes

1. I have tried to remove 'mess', including things like functions used in previous versions but which are now never called, in order to make it easier to understand. I have not tested it since doing so.
2. The code uses a finer grid spacing at the downslope edge of the blister, but not at the upslope edge. This is fine when you are modelling the movement of the blister down a slope, but means it is not suited to horizontally spreading blisters where you would like well resolved upslope peeling.
3. It is slow. The original version of this code was written in FORTRAN, but I switched to Python when I realised that I would have to regrid regularly and it would be nice to be able to just add extra elements to arrays. The trade off between time spent writing the code and time spent running the code worked for me, but if you are starting from scratch you may be best using this as a reference and rewriting it to be more efficient.
