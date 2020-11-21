This project serves as an extension of LAMMPS to implement a pair style called 'ops', which can compute the pairwise interactions in an oriented particle system. 

Steps to use:
1. Download pair_ops.cpp, pair_ops.h and save them within the folder /lammps/src.
2. Build LAMMPS.
3. Use it like any other pair_style.

pair_style ops
======================================

Syntax
""""""

.. code-block:: LAMMPS
   pair_style style args

* style = *ops*
* args = list of arguments for a particular style

.. parsed-literal::

    *ops* args = cutoff
      cutoff = global cutoff for Morse interactions (distance units)

Examples
""""""""
.. code-block:: LAMMPS

   pair_style ops 2.5
   pair_coeff * * 100.0 2.0 1.5 1 0 1 1 1 1 1 1
   
Description
"""""""""""

Style *ops* computes pairwise interactions with the formula:

.. math::

   E = \alpha_m \phi_m(r_{ij}) + \alpha_p \phi_p(n_i,r_{ij}) + \alpha_n \phi_n(n_i,n_j,r_{ij}) + \alpha_c \phi_c(n_i,,n_j,r_{ij})
       \qquad r < r_c

where

.. math::

   \phi_m(r_{ij}) = D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right]
   \phi_p(n_i,r_{ij}) = (n_i \cdot r_{ij})^2\psi(r_{ij})
   \phi_n(n_i,n_j,r_{ij}) = |n_i - n_j|^2\psi(r_{ij})
   \phi_c(n_i,n_j,r_{ij}) = ((n_i + n_j) \cdot r_{ij})^2\psi(r_{ij})
   \psi(r_{ij}) = K e^{(-\frac{x^2}{2a^2} -\frac{y^2}{2b^2} -\frac{z^2}{2c^2})}
   
Rc is the cutoff.
:math:'D_0, K' are constants.
:math:'\alpha_m, \alpha_p, \alpha_n, \alpha_c' are weight constants for Morse potential, Co-planarity potential, Co-normality potential and Co-Circularity potential, respectively.
:math:'a, b, c' are constants for weighting function :math:'\psi(r_{ij})'. The weighting function provides either circular symmetricity (a=b=c), or weights particles more in certain planes.

The 4 potentials used are Morse Potential, Co-Planarity Potential, Co-Normality Potential and Co-Corcularit Potential. Together, their weighted sum defines the potential energy of oriented particles.

The following coefficients must be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands:

* :math:`D_0` (energy units)
* :math:`\alpha` (1/distance units)
* :math:`r_0` (distance units)
* :math:`\alpha_m` (energy units)
* :math:`\alpha_p` (energy units)
* :math:`\alpha_n` (energy units)
* :math:`\alpha_c` (energy units)
* :math:`K`
* :math:`a` (1/distance units)
* :math:`b` (1/distance units)
* :math:`c` (1/distance units)
* cutoff (distance units)

The last coefficient is optional.  If not specified, the global morse cutoff is used.

These pairs only work with atom_style sphere. In the input file, a fix needs to be defined to input the normals of each of the atoms. An example is given below:

.. code-block:: LAMMPS
   fix Vectors all property/atom d_normalX d_normalY d_normalZ

With this example, the read_data command will look like this:
.. code-block:: LAMMPS
   read_data <data.filename> fix Vectors NULL Normals

'Normals' will be the header in the data file, below which the atom ID along with three normal components (along x, y, and z axes) given for each atom. 

This example is implemented in the test folder given with this project.
