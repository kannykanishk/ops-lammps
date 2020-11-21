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

![f1](https://latex.codecogs.com/svg.latex?\small&space;E&space;=&space;\alpha_m&space;\phi_m(r_{ij})&space;&plus;&space;\alpha_p&space;\phi_p(n_i,r_{ij})&space;&plus;&space;\alpha_n&space;\phi_n(n_i,n_j,r_{ij})&space;&plus;&space;\alpha_c&space;\phi_c(n_i,,n_j,r_{ij})&space;\qquad&space;r&space;<&space;r_c)

where

   ![f2](https://latex.codecogs.com/svg.latex?\small&space;\phi_m(r_{ij})&space;=&space;D_0&space;\left[&space;e^{-&space;2&space;\alpha&space;(r&space;-&space;r_0)}&space;-&space;2&space;e^{-&space;\alpha&space;(r&space;-&space;r_0)}&space;\right])
   
   ![f3](https://latex.codecogs.com/svg.latex?\small&space;\phi_p(n_i,r_{ij})&space;=&space;(n_i&space;\cdot&space;r_{ij})^2\psi(r_{ij}))
   
   ![f4](https://latex.codecogs.com/svg.latex?\small&space;\phi_n(n_i,n_j,r_{ij})&space;=&space;|n_i&space;-&space;n_j|^2\psi(r_{ij}))
   
   ![f5](https://latex.codecogs.com/svg.latex?\small&space;\phi_c(n_i,n_j,r_{ij})&space;=&space;((n_i&space;&plus;&space;n_j)&space;\cdot&space;r_{ij})^2\psi(r_{ij}))
   
   ![f6](https://latex.codecogs.com/svg.latex?\small&space;\psi(r_{ij})=Ke^{(-\frac{x^2}{2a^2}-\frac{y^2}{2b^2}-\frac{z^2}{2c^2})})
   
- Rc is the cutoff.
- D₀, K are constants.
- ![f7](https://latex.codecogs.com/svg.latex?\small&space;\alpha_m), ![f8](https://latex.codecogs.com/svg.latex?\small&space;\alpha_p), ![f9](https://latex.codecogs.com/svg.latex?\small&space;\alpha_n), ![f10](https://latex.codecogs.com/svg.latex?\small&space;\alpha_c) are weight constants for Morse potential, Co-planarity potential, Co-normality potential and Co-Circularity potential, respectively.
- a, b, c are constants for weighting function ![f11](https://latex.codecogs.com/svg.latex?\small&space;\psi(r_{ij})). The weighting function provides either circular symmetricity (a=b=c), or weights particles more in certain planes.

The 4 potentials used are Morse Potential, Co-Planarity Potential, Co-Normality Potential and Co-Corcularit Potential. Together, their weighted sum defines the potential energy of oriented particles.

The first three coefficients and the cutoff can be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands (from LAMMPS doc). The remaining coefficients remain same for all atom pairs, hence needs to be specified in the main input file alongside `pair_coeff * *`:

* D₀ (energy units)
* ![f12](https://latex.codecogs.com/svg.latex?\small&space;\alpha) (1/distance units)
* r₀ (distance units)
* ![f7](https://latex.codecogs.com/svg.latex?\small&space;\alpha_m) (energy units) 
* ![f8](https://latex.codecogs.com/svg.latex?\small&space;\alpha_p) (energy units)
* ![f9](https://latex.codecogs.com/svg.latex?\small&space;\alpha_n) (energy units)
* ![f10](https://latex.codecogs.com/svg.latex?\small&space;\alpha_c) (energy units)
* K
* a (1/distance units)
* b (1/distance units)
* c (1/distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global morse cutoff is used.

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

