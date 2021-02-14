This project serves as an extension of LAMMPS to implement a pair style called 'ops', which can compute the pairwise interactions in an oriented particle system. 

Steps to use:
1. Download pair_ops.cpp, pair_ops.h and save them within the folder /lammps/src.
2. Build LAMMPS.
3. Use it like any other pair_style.

# pair_style ops

### Syntax

`pair_style style args`

* style = *ops*
* args = list of arguments for a particular style
```
    ops args = cutoff
    cutoff = global cutoff for Morse interactions (distance units)
```
### Examples
```
   pair_style ops 2.5
   pair_coeff * * 100.0 2.0 1.5 1 0 1 1 1 1 1 1
``` 
### Description

Style *ops* computes pairwise interactions with the formula:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?E%3D%5Calpha_m%5Cphi_m%28r_%7Bij%7D%29&plus;%5Calpha_p%5Cphi_p%28n_i%2Cr_%7Bij%7D%29&plus;%5Calpha_n%5Cphi_n%28n_i%2Cn_j%2Cr_%7Bij%7D%29&plus;%5Calpha_c%5Cphi_c%28n_i%2Cn_j%2Cr_%7Bij%7D%29%5Cquad%20r%3Cr_c" align=middle /></p>

where

<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cphi_m%28r_%7Bij%7D%29%3DD_0%5Cleft%5Be%5E%7B-2%5Calpha%28r-r_0%29%7D-2e%5E%7B-%5Calpha%28r-r_0%29%7D%5Cright%5D" align=middle /></p>
<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cphi_p%28n_i%2Cr_%7Bij%7D%29%3D%28n_i.r_%7Bij%7D%29%5E2%5Cpsi%28r_%7Bij%7D%29" align=middle /></p>
<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cphi_n%28n_i%2Cn_j%2Cr_%7Bij%7D%29%3D%7Cn_i-n_j%7C%5E2%5Cpsi%28r_%7Bij%7D%29" align=middle /></p>
<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cphi_c%28n_i%2Cn_j%2Cr_%7Bij%7D%29%3D%28%28n_i&plus;n_j%29.r_%7Bij%7D%29%5E2%5Cpsi%28r_%7Bij%7D%29" align=middle /></p>
<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cpsi%28r_%7Bij%7D%29%3DKe%5E%7B%28-%5Cfrac%7Bx%5E2%7D%7B2a%5E2%7D-%5Cfrac%7By%5E2%7D%7B2b%5E2%7D-%5Cfrac%7Bz%5E2%7D%7B2c%5E2%7D%29%7D" align=middle /></p>
   
- <img src="https://latex.codecogs.com/gif.latex?r_c"/> is the cutoff.
- Do, K are constants.
- <img src="https://latex.codecogs.com/gif.latex?%5Calpha_m" />, <img src="https://latex.codecogs.com/gif.latex?%5Calpha_p" />, <img src="https://latex.codecogs.com/gif.latex?%5Calpha_n" />, <img src="https://latex.codecogs.com/gif.latex?%5Calpha_c" /> are weight constants for Morse potential, Co-planarity potential, Co-normality potential and Co-Circularity potential, respectively.
- a, b, c are constants for weighting function <img src="https://latex.codecogs.com/gif.latex?%5Cpsi%28r_%7Bij%7D%29" />. The weighting function provides either circular symmetricity <img src="https://latex.codecogs.com/gif.latex?%5Cpsi%28r_%7Bij%7D%29" />, or weights particles more in certain planes.

The 4 potentials used are Morse Potential, Co-Planarity Potential, Co-Normality Potential and Co-Corcularity Potential. Together, their weighted sum defines the potential energy of oriented particles.

The first three coefficients and the cutoff can be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands (from LAMMPS doc). The remaining coefficients remain same for all atom pairs, hence needs to be specified in the main input file alongside `pair_coeff * *`:

* <img src="https://latex.codecogs.com/gif.latex?D_0"/> (energy units)
* <img src="https://latex.codecogs.com/gif.latex?%5Calpha"/> (1/distance units)
* <img src="https://latex.codecogs.com/gif.latex?r_0"/> (distance units)
* <img src="https://latex.codecogs.com/gif.latex?%5Calpha_m"/> (energy units) 
* <img src="https://latex.codecogs.com/gif.latex?%5Calpha_p"/> (energy units)
* <img src="https://latex.codecogs.com/gif.latex?%5Calpha_n"/> (energy units)
* <img src="https://latex.codecogs.com/gif.latex?%5Calpha_c"/> (energy units)
* K
* a (distance units)
* b (distance units)
* c (distance units)
* cutoff <img src="https://latex.codecogs.com/gif.latex?r_c"/> (distance units)

The last coefficient is optional. If not specified, the global morse cutoff is used.

### Restrictions

These pairs only work with `atom_style sphere`. In the input file, a fix needs to be defined to input the normals of each of the atoms. An example is given below:

   `fix Vectors all property/atom d_normalX d_normalY d_normalZ`

With this example, the `read_data` command will look like this:

   `read_data <data.filename> fix Vectors NULL Normals`

'Normals' will be the header in the data file, below which the atom ID along with three normal components (along x, y, and z axes) given for each atom. 

This example is implemented in the test folder given with this project.

### References

[1] Szeliski, R. & Tonnesen, D. (1992). Surface Modeling with Oriented Particle Systems. Computer Graphics, 26(2):185-194.

[2] Tonnesen, D. (1998). Dynamically Coupled Particle Systems for Geometric Modeling, Reconstruction, and Animation. University of Toronto. Retrieved from http://www.dgp.toronto.edu/~davet/phd/index.html

[3] Singh, A.R. (2018). Study of Zero and Finite Temperature Response of Discrete Deformable Surfaces. UCLA. ProQuest ID: Singh_ucla_0031D_17348. Merritt ID: ark:/13030/m5sz1b5x. Retrieved from https://escholarship.org/uc/item/6cg965pq


