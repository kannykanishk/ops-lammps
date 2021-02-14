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

<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/9e752f3e479f4aaeb86c0db3bace4e65.svg?invert_in_darkmode" align=middle width=1410.9420731999999pt height=17.031940199999998pt/></p>

where

<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/566dea83f7d54e27fc2a0218f8b9f575.svg?invert_in_darkmode" align=middle width=1013.26864155pt height=124.01050529999999pt/></p>
   
- Rc is the cutoff.
- D₀, K are constants.
- <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/1d697f543cf2bf4c0de3dd3d432b37e8.svg?invert_in_darkmode" align=middle width=68.91952154999998pt height=14.15524440000002pt/>, <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/728f77c88885f33fa56823dbf6119c1f.svg?invert_in_darkmode" align=middle width=64.03115069999998pt height=14.15524440000002pt/>, <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/015924a4daa52ba03d8e9f17bc18c245.svg?invert_in_darkmode" align=middle width=65.3806956pt height=14.15524440000002pt/>, <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/316834b0d33ea5e719cc38f24a9d6a7d.svg?invert_in_darkmode" align=middle width=63.12932504999999pt height=14.15524440000002pt/> are weight constants for Morse potential, Co-planarity potential, Co-normality potential and Co-Circularity potential, respectively.
- a, b, c are constants for weighting function <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/79f249a9adbe4e092b9224d6f33ed554.svg?invert_in_darkmode" align=middle width=89.81564954999999pt height=24.65753399999998pt/>. The weighting function provides either circular symmetricity (a=b=c), or weights particles more in certain planes.

The 4 potentials used are Morse Potential, Co-Planarity Potential, Co-Normality Potential and Co-Corcularity Potential. Together, their weighted sum defines the potential energy of oriented particles.

The first three coefficients and the cutoff can be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands (from LAMMPS doc). The remaining coefficients remain same for all atom pairs, hence needs to be specified in the main input file alongside `pair_coeff * *`:

* D₀ (energy units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/d0d718f9fc184f8a992bae19b12785dc.svg?invert_in_darkmode" align=middle width=57.31553024999999pt height=14.15524440000002pt/> (1/distance units)
* r₀ (distance units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/1d697f543cf2bf4c0de3dd3d432b37e8.svg?invert_in_darkmode" align=middle width=68.91952154999998pt height=14.15524440000002pt/> (energy units) 
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/728f77c88885f33fa56823dbf6119c1f.svg?invert_in_darkmode" align=middle width=64.03115069999998pt height=14.15524440000002pt/> (energy units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/015924a4daa52ba03d8e9f17bc18c245.svg?invert_in_darkmode" align=middle width=65.3806956pt height=14.15524440000002pt/> (energy units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/316834b0d33ea5e719cc38f24a9d6a7d.svg?invert_in_darkmode" align=middle width=63.12932504999999pt height=14.15524440000002pt/> (energy units)
* K
* a (distance units)
* b (distance units)
* c (distance units)
* cutoff (distance units)

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

