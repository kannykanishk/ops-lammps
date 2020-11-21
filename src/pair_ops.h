/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(ops,PairOPS)

#else

#ifndef LMP_PAIR_OPS_H
#define LMP_PAIR_OPS_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOPS : public Pair {
 public:
  PairOPS(class LAMMPS *);
  virtual ~PairOPS();
  virtual void compute(int, int);

  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);
    
    //Custom functions for calculating
    double dotProduct(double *, double *);
    double *crossProduct(double *, double *);
    void normalize(double *);
    double getModulus(double *);
    double *getNormalVector(double *, double);
    void multiplyMatrices(double **, double **, double **, int, int, int, int);
    double phi_p(double *, double *);
    double phi_n(double *, double *);
    double phi_c(double *, double *, double *);
    double WeightFunction(double *);
    double WeightDerivative(double *, int);
    double *f_p(double *, double *);
    double *t_p(double *, double *);
    double *f_n(double *, double *, double *);
    double *t_n(double *, double *, double *);
    double *f_c(double *, double *, double *);
    double *t_c(double *, double *, double *);

 protected:
  double cut_global, alphaM, alphaC, alphaP, alphaN, Kconstant, a, b, c;
  double **cut;
  double **d0,**alpha,**r0;
  double **morse1;
  double **offset;

  double *normalX, *normalY, *normalZ;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
