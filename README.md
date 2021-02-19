# Umbrella sampling simulation uses RMSD as the reaction coordinate. 
# The RMSD at each timestep is measured by the Singular Value Decomposition method. Protein have a transition between the two states of Outward Facing Open (OF) and Outwar Faceing Occluded state by RMSD difference at each soordinate. The reaction coordinate is rescaled to have a RMSD between 0 and 1. RC=O is OC and  RC=1 is OF states.
# The code uses the Eigen library of C++
# Force_ext.sh needs for compiling Frc_ext.c
