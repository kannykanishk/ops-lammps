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

<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/439980efd4213d466cbdae06aef61cfb.svg?invert_in_darkmode" align=middle width=558.050394pt height=17.031940199999998pt/></p>

where

<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/61791d2792d07f249c60858eba21f8e1.svg?invert_in_darkmode" align=middle width=283.73085345pt height=29.58934275pt/></p>
<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/6d3d32472a734e3c93ef5e99efab2688.svg?invert_in_darkmode" align=middle width=195.87115845pt height=18.905967299999997pt/></p>
<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/7497228fe82b44cf6aadbce7c121a0ed.svg?invert_in_darkmode" align=middle width=230.99154375pt height=18.905967299999997pt/></p>
<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/128054bb9bfd59e13821ef64f1400fab.svg?invert_in_darkmode" align=middle width=268.73841884999996pt height=18.905967299999997pt/></p>
<p align="center"><img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/a85311bb572d87be46bd80d72aa228d6.svg?invert_in_darkmode" align=middle width=193.65341819999998pt height=26.077949699999998pt/></p>
   
- <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/dc118d444558d9150b042da57bed75ea.svg?invert_in_darkmode" align=middle width=13.290972749999991pt height=14.15524440000002pt/> is the cutoff.
- <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/42fbe61bcc38cd8c7cb8a2f4abc9b4c7.svg?invert_in_darkmode" align=middle width=20.16214364999999pt height=22.465723500000017pt/>, <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.13700594999999pt height=22.465723500000017pt/> are constants.
- <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/812375f957f3b78a317a4d7150c3ae73.svg?invert_in_darkmode" align=middle width=22.18049624999999pt height=14.15524440000002pt/>, <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/6c437e6db92ada6d18b2c92382e20721.svg?invert_in_darkmode" align=middle width=17.292125399999993pt height=14.15524440000002pt/>, <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/17275f74aac37adec4b7e565a0f8199a.svg?invert_in_darkmode" align=middle width=18.64167029999999pt height=14.15524440000002pt/>, <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/971651e87c9aebb6ec102860c98ae161.svg?invert_in_darkmode" align=middle width=16.390298099999992pt height=14.15524440000002pt/> are weight constants for Morse potential, Co-planarity potential, Co-normality potential and Co-Circularity potential, respectively.
- <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/c7511ce56cd9c8457f7a29917f39df8d.svg?invert_in_darkmode" align=middle width=37.46952164999999pt height=22.831056599999986pt/> are constants for weighting function <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/a5c4119e5685d6ff54503ba7c63bbcdb.svg?invert_in_darkmode" align=middle width=43.07662424999999pt height=24.65753399999998pt/>. The weighting function provides either circular symmetricity (<img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/b11e6eeff71a96de798a5e4addc53da9.svg?invert_in_darkmode" align=middle width=66.6930165pt height=22.831056599999986pt/>), or weights particles more in certain planes.

The 4 potentials used are Morse Potential, Co-Planarity Potential, Co-Normality Potential and Co-Corcularity Potential. Together, their weighted sum defines the potential energy of oriented particles.

The first three coefficients and the cutoff can be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands (from LAMMPS doc). The remaining coefficients remain same for all atom pairs, hence needs to be specified in the main input file alongside `pair_coeff * *`:

* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/084ad41df3595da9790e07d7d2fc9849.svg?invert_in_darkmode" align=middle width=20.09819954999999pt height=22.465723500000017pt/> (energy units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/c745b9b57c145ec5577b82542b2df546.svg?invert_in_darkmode" align=middle width=10.57650494999999pt height=14.15524440000002pt/> (1/distance units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/c1deb563f9f030f805445b5f3197054c.svg?invert_in_darkmode" align=middle width=13.904922899999988pt height=14.15524440000002pt/> (distance units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/812375f957f3b78a317a4d7150c3ae73.svg?invert_in_darkmode" align=middle width=22.18049624999999pt height=14.15524440000002pt/> (energy units) 
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/6c437e6db92ada6d18b2c92382e20721.svg?invert_in_darkmode" align=middle width=17.292125399999993pt height=14.15524440000002pt/> (energy units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/17275f74aac37adec4b7e565a0f8199a.svg?invert_in_darkmode" align=middle width=18.64167029999999pt height=14.15524440000002pt/> (energy units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/971651e87c9aebb6ec102860c98ae161.svg?invert_in_darkmode" align=middle width=16.390298099999992pt height=14.15524440000002pt/> (energy units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/d6328eaebbcd5c358f426dbea4bdbf70.svg?invert_in_darkmode" align=middle width=15.13700594999999pt height=22.465723500000017pt/>
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/44bc9d542a92714cac84e01cbbb7fd61.svg?invert_in_darkmode" align=middle width=8.68915409999999pt height=14.15524440000002pt/> (distance units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/4bdc8d9bcfb35e1c9bfb51fc69687dfc.svg?invert_in_darkmode" align=middle width=7.054796099999991pt height=22.831056599999986pt/> (distance units)
* <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/3e18a4a28fdee1744e5e3f79d13b9ff6.svg?invert_in_darkmode" align=middle width=7.11380504999999pt height=14.15524440000002pt/> (distance units)
* cutoff <img src="https://rawgit.com/kannykanishk/ops-lammps/main/svgs/dc118d444558d9150b042da57bed75ea.svg?invert_in_darkmode" align=middle width=13.290972749999991pt height=14.15524440000002pt/> (distance units)

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


