/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_ops.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOPS::PairOPS(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  centroidstressflag = 1;
}

/* ---------------------------------------------------------------------- */

PairOPS::~PairOPS()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(r0);
    memory->destroy(morse1);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairOPS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair, fx, fy, fz, tx, ty, tz;
  double rsq,r,dr,dexp,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *fp, *fn, *fc, *tp, *tn, *tc, r_ij[3];
  double u[3], v[3];

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
                                            
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
      
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
        
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      r_ij[0] = delx;
      r_ij[1] = dely;
      r_ij[2] = delz;
        
      if (rsq < cutsq[itype][jtype]) {
        u[0] = normalX[i];
        u[1] = normalY[i];
        u[2] = normalZ[i];

        v[0] = normalX[j];
        v[1] = normalY[j];
        v[2] = normalZ[j];
          
        r = sqrt(rsq);
        dr = r - r0[itype][jtype];
        dexp = exp(-alpha[itype][jtype] * dr);
          
        fpair = factor_lj * morse1[itype][jtype] * (dexp*dexp - dexp) / r;
        fp = f_p(u, r_ij);
        fc = f_c(u, v, r_ij);
        fn = f_n(u, v, r_ij);
        fx = alphaM*delx*fpair + alphaP*fp[0] + alphaC*fc[0] + alphaN*fn[0];
        fy = alphaM*dely*fpair + alphaP*fp[1] + alphaC*fc[1] + alphaN*fn[1];
        fz = alphaM*delz*fpair + alphaP*fp[2] + alphaC*fc[2] + alphaN*fn[2];
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
          
        tp = t_p(u, r_ij);
        tc = t_c(u, v, r_ij);
        tn = t_n(u, v, r_ij);
        tx = alphaC*tc[0] + alphaN*tn[0];
        ty = alphaC*tc[1] + alphaN*tn[1];
        tz = alphaC*tc[2] + alphaN*tn[2];
        torque[i][0] += tx + alphaP*tp[0];
        torque[i][1] += ty + alphaP*tp[0];
        torque[i][2] += tz + alphaP*tp[0];
          
        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= tx;
          torque[j][1] -= ty;
          torque[j][2] -= tz;
        }

        if (eflag) {
          evdwl = (d0[itype][jtype] * (dexp*dexp - 2.0*dexp)) - offset[itype][jtype] + alphaC*phi_c(u, v, r_ij) + alphaN*phi_n(u, v) + alphaP*phi_p(u, r_ij);
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairOPS::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(d0,n+1,n+1,"pair:d0");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(morse1,n+1,n+1,"pair:morse1");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairOPS::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOPS::coeff(int narg, char **arg)
{
  if (narg < 13 || narg > 14)                                       //Extra arguments
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double d0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double alpha_one = utils::numeric(FLERR,arg[3],false,lmp);
  double r0_one = utils::numeric(FLERR,arg[4],false,lmp);
  double alphaM_one = utils::numeric(FLERR,arg[5],false,lmp);         //Weight of Morse argument from input file
  double alphaP_one = utils::numeric(FLERR,arg[6],false,lmp);         //Weight of Co-planarity argument from input file
  double alphaN_one = utils::numeric(FLERR,arg[7],false,lmp);         //Weight of Co-normality argument from input file
  double alphaC_one = utils::numeric(FLERR,arg[8],false,lmp);         //Weight of Co-circularity argument from input file
  double Kconstant_one = utils::numeric(FLERR,arg[9],false,lmp);      //Kconstant argument from input file
  double a_weight = utils::numeric(FLERR,arg[10],false,lmp);          //a Weighting function argument from input file
  double b_weight = utils::numeric(FLERR,arg[11],false,lmp);          //b Weighting function argument from input file
  double c_weight = utils::numeric(FLERR,arg[12],false,lmp);          //c Weighting function argument from input file


  double cut_one = cut_global;
  if (narg == 14) cut_one = utils::numeric(FLERR,arg[13],false,lmp);    //Cutoff distance coefficient (optional)

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      d0[i][j] = d0_one;
      alpha[i][j] = alpha_one;
      r0[i][j] = r0_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  alphaM = alphaM_one;
  alphaP = alphaP_one;
  alphaN = alphaN_one;
  alphaC = alphaC_one;
  Kconstant = Kconstant_one;
  a = a_weight;
  b = b_weight;
  c = c_weight;
    
  int fixFlag = 1;
  int fixIndex;
  fixIndex = atom->find_custom("normalX", fixFlag);
  normalX = atom->dvector[fixIndex];
  fixIndex = atom->find_custom("normalY", fixFlag);
  normalY = atom->dvector[fixIndex];
  fixIndex = atom->find_custom("normalZ", fixFlag);
  normalZ = atom->dvector[fixIndex];
    
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairOPS::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  morse1[i][j] = 2.0*d0[i][j]*alpha[i][j];

  if (offset_flag) {
    double alpha_dr = -alpha[i][j] * (cut[i][j] - r0[i][j]);
    offset[i][j] = d0[i][j] * (exp(2.0*alpha_dr) - 2.0*exp(alpha_dr));
  } else offset[i][j] = 0.0;

  d0[j][i] = d0[i][j];
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  morse1[j][i] = morse1[i][j];
  offset[j][i] = offset[i][j];
  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOPS::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOPS::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&d0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&alpha[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&r0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairOPS::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairOPS::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairOPS::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,d0[i][i],alpha[i][i],r0[i][i]);         //!
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairOPS::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",
              i,j,d0[i][j],alpha[i][j],r0[i][j],cut[i][j]);                 //!
}

/* ---------------------------------------------------------------------- */

double PairOPS::single(int i, int j, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj,
                         double &fforce)
{
  double r,dr,dexp,phi;

  r = sqrt(rsq);
  dr = r - r0[itype][jtype];
  dexp = exp(-alpha[itype][jtype] * dr);
  fforce = factor_lj * morse1[itype][jtype] * (dexp*dexp - dexp) / r;

  phi = d0[itype][jtype] * (dexp*dexp - 2.0*dexp) - offset[itype][jtype];
  return factor_lj*phi;
}

/* ---------------------------------------------------------------------- */

void *PairOPS::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"d0") == 0) return (void *) d0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"alpha") == 0) return (void *) alpha;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

/* Custom functions for vector operations */

double PairOPS::dotProduct(double *v1, double *v2) {
    double ans = 0;
    for (int i=0; i<=2; i++) {
        ans += v1[i] * v2[i];
    }
    return ans;
}

double *PairOPS::crossProduct(double *a, double *b) {
    static double ans[3];
    ans[0] = a[1]*b[2] - a[2]*b[1];
    ans[1] = a[2]*b[0] - a[0]*b[2];
    ans[2] = a[0]*b[1] - a[1]*b[0];
    return ans;
}

double PairOPS::getModulus(double *a) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

void PairOPS::normalize(double *a) {
    double modValue = getModulus(a);
    for (int i=0; i<=2; i++) {
        a[i] = a[i] / modValue;
    }
}

double *PairOPS::getNormalVector(double *u, double ang) {
    static double normalVector[3] = {0,0,0};
    double *p;
    p=normalVector;
    double U = getModulus(u);
    p[0] = 2 * sin(ang) * ( u[1]*U*cos(ang) + u[0]*u[2]*sin(ang) ) / (U*U);
    p[1] = 2 * sin(ang) * (-u[0]*U*cos(ang) + u[1]*u[2]*sin(ang) ) / (U*U);
    p[2] = cos(ang)*cos(ang) + sin(ang)*sin(ang)*(u[2]*u[2] - u[1]*u[1] - u[0]*u[0])/(U*U);
    return p;
}

void PairOPS::multiplyMatrices(double **res, double **a, double **b, int aRows, int aCols, int bRows, int bCols)
{
    if (aCols != bRows)                 //Matrix Multiplication Check
        return;

    for (int i = 0; i < aRows; i++)
    {
        for (int j = 0; j < bCols; j++)
        {
            res[i][j] = 0;
            for (int k = 0; k < aCols; k++)
            {
                res[i][j] += a[i][k]*b[k][j];
            }
        }
    }
}

double PairOPS::phi_p(double *n, double *r) {
    double Phi_p = 0;
    Phi_p = dotProduct(n, r);
    return Phi_p*Phi_p;
}

double PairOPS::phi_n(double *p, double *q) {
    double Phi_n = 0, k;
    for (int i=0; i<3; i++) {
        k = p[i] - q[i];
        Phi_n += k*k;
    }
    return Phi_n;
}

double PairOPS::phi_c(double *p, double *q, double *r) {
    double m[3];
    for (int i=0; i<3; i++) {
        m[i] = p[i] + q[i];
    }
    double Phi_c = dotProduct(m, r);
    return Phi_c*Phi_c;
}

double PairOPS::WeightFunction(double *r) {
    return Kconstant*exp(-(r[0]*r[0])/(2*a*a) - (r[1]*r[1])/(2*b*b) - (r[2]*r[2])/(2*c*c));
}

double PairOPS::WeightDerivative(double *r, int i) {
    return 0;
    // Uncomment the below lines, comment the top line to implement the weight derivative formula
    /*double phi = WeightFunction(r);
    if (i==0) {
        return phi*(-r[0]/(a*a));
    }
    else if (i==1) {
        return phi*(-r[1]/(b*b));
    }
    else {
        return phi*(-r[2]/(c*c));
    }*/
}

double *PairOPS::f_p(double *n, double *r) {
    static double fp[3] = {0,0,0};
    double *p;
    double dot = dotProduct(n, r);
    p=fp;
    for(int i=0; i<2; i++) {
        p[i] = -n[i] * dot * WeightFunction(r) -
                (r[i]/getModulus(r)) * dot * dot * WeightDerivative(r, i);
    }
    return p;
}

double *PairOPS::f_n(double *ni, double *nj, double *r) {
    static double fn[3] = {0,0,0};
    double *p;
    double k[3];
    p=fn;
    for(int i=0; i<2; i++) {
        k[i] = ni[i] - nj[i];
    }
    double modK = getModulus(k);
    for(int i=0; i<2; i++) {
        p[i] = -r[i] * modK * modK * WeightDerivative(r,i);
    }
    return p;
}

double *PairOPS::f_c(double *ni, double *nj, double *r) {
    static double fc[3] = {0,0,0};
    double *p;
    double k[3];
    for (int i=0; i<2; i++) {
        k[i] = ni[i] + nj[i];
    }
    double dot = dotProduct(k, r);
    p=fc;
    for(int i=0; i<2; i++) {
        p[i] = -ni[i] * dot * WeightFunction(r) -
                (r[i]/getModulus(r)) * dot * dot * WeightDerivative(r, i);
    }
    return p;
}

double *PairOPS::t_p(double *ni, double *r) {
    double *p;
    double f[3];
    double dot = dotProduct(ni, r);
    double w = WeightFunction(r);
    for(int i=0; i<2; i++) {
        f[i] = ni[i] * dot * w;
    }
    p = crossProduct(r, f);
    return p;
}

double *PairOPS::t_n(double *ni, double *nj, double *r) {
    double *p;
    p = crossProduct(ni, nj);
    double w = WeightFunction(r);
    for (int i=0; i<2; i++) {
        p[i] *= w;
    }
    return p;
}

double *PairOPS::t_c(double *ni, double *nj, double *r) {
    double *p;
    double f[3], k[3];
    for (int i=0; i<2; i++) {
        k[i] = ni[i] + nj[i];
    }
    double w = WeightFunction(r);
    double dot = dotProduct(k, r);
    for(int i=0; i<2; i++) {
        f[i] = ni[i] * dot * w;
    }
    p = crossProduct(r, f);
    return p;
}
