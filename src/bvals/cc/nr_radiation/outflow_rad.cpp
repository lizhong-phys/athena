//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow_rad.cpp
//  \brief implementation of outflow BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "./bvals_rad.hpp"


// The angular grid can change
// The angular octant is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in radiatin class, n_ang is angles per octant, noct is the number of octant
// radiation relfection means set specific intensitiy corresponding
// the value in the opposite octant


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::OutflowInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief outflow boundary conditions, inner x1 boundary

void RadBoundaryVariable::OutflowInnerX1(Real time, Real dt, int il, int jl, int ju,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  NRRadiation *prad = pmy_block_->pnrrad;
  const int& nang  = prad->nang; // angles per octant
  const int& nfreq = prad->nfreq; // number of frequency bands
  const int& nstok = prad->num_stokes; // number of stokes parameters

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        if (!prad->use_pol_rad) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            for (int n=0; n<nang; ++n) {
              int ang=ifr*nang+n;
              (*var_cc)(k,j,il-i,ang) = (*var_cc)(k,j,il,ang);
            } // endfor n
          } // endfor ifr
        } else { // modifications for polarization
          for (int m=0; m<nstok; ++m) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              for (int n=0; n<nang; ++n) {
                int ang=ifr*nang+n;
                (*var_cc)(k,j,il-i,m,ang) = (*var_cc)(k,j,il,m,ang);
              } // endfor n
            } // endfor ifr
          } // endfor m
        } // endelse !prad->use_pol_rad
      } // endfor i
    } // endfor j
  } // endfor k

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::OutflowOuterX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief Outflow boundary conditions, outer x1 boundary

void RadBoundaryVariable::OutflowOuterX1(Real time, Real dt, int iu, int jl, int ju,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  NRRadiation *prad = pmy_block_->pnrrad;
  const int& nang  = prad->nang; // angles per octant
  const int& nfreq = prad->nfreq; // number of frequency bands
  const int& nstok = prad->num_stokes; // number of stokes parameters

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        if (!prad->use_pol_rad) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            for (int n=0; n<nang; ++n) {
              int ang=ifr*nang+n;
              (*var_cc)(k,j,iu+i,ang) = (*var_cc)(k,j,iu,ang);
            } // endfor n
          } // endfor ifr
        } else { // modifications for polarization
          for (int m=0; m<nstok; ++m) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              for (int n=0; n<nang; ++n) {
                int ang=ifr*nang+n;
                (*var_cc)(k,j,iu+i,m,ang) = (*var_cc)(k,j,iu,m,ang);
              } // endfor n
            } // endfor ifr
          } // endfor m
        } // endelse !prad->use_pol_rad
      } // endfor i
    } // endfor j
  } // endfor k

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::OutflowInnerX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief Outflow boundary conditions, inner x2 boundary

void RadBoundaryVariable::OutflowInnerX2(Real time, Real dt, int il, int iu, int jl,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  NRRadiation *prad = pmy_block_->pnrrad;
  int nang = prad->nang;
  const int& nfreq = prad->nfreq; // number of frequency bands
  const int& nstok = prad->num_stokes; // number of stokes parameters

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        if (!prad->use_pol_rad) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            for (int n=0; n<nang; ++n) {
              int ang=ifr*nang+n;
              (*var_cc)(k,jl-j,i,ang) = (*var_cc)(k,jl,i,ang);
            } // endfor n
          } // endfor ifr
        } else { // modifications for polarization
          for (int m=0; m<nstok; ++m) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              for (int n=0; n<nang; ++n) {
                int ang=ifr*nang+n;
                (*var_cc)(k,jl-j,i,m,ang) = (*var_cc)(k,jl,i,m,ang);
              } // endfor n
            } // endfor ifr
          } // endfor m
        } // endelse !prad->use_pol_rad
      } // endfor i
    } // endfor j
  } // endfor k

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::OutflowOuterX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief Outflow boundary conditions, outer x2 boundary

void RadBoundaryVariable::OutflowOuterX2(Real time, Real dt, int il, int iu, int ju,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  NRRadiation *prad = pmy_block_->pnrrad;
  const int& nang  = prad->nang; // angles per octant
  const int& nfreq = prad->nfreq; // number of frequency bands
  const int& nstok = prad->num_stokes; // number of stokes parameters

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        if (!prad->use_pol_rad) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            for(int n=0; n<nang; ++n) {
              int ang=ifr*nang+n;
              (*var_cc)(k,ju+j,i,ang) = (*var_cc)(k,ju,i,ang);
            } // endfor n
          } // endfor ifr
        } else { // modifications for polarization
          for (int m=0; m<nstok; ++m) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              for(int n=0; n<nang; ++n) {
                int ang=ifr*nang+n;
                (*var_cc)(k,ju+j,i,m,ang) = (*var_cc)(k,ju,i,m,ang);
              } // endfor n
            } // endfor ifr
          } // endfor m
        } // endelse !prad->use_pol_rad
      } // endfor i
    } // endfor j
  } // endfor k

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::OutflowInnerX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)

void RadBoundaryVariable::OutflowInnerX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int kl, int ngh) {
  // copy radiation variables into ghost zones,
  NRRadiation *prad = pmy_block_->pnrrad;
  const int& nang  = prad->nang;
  const int& nfreq = prad->nfreq; // number of frequency bands
  const int& nstok = prad->num_stokes; // number of stokes parameters

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        if (!prad->use_pol_rad) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            for (int n=0; n<nang; ++n) {
              int ang=ifr*nang+n;
              (*var_cc)(kl-k,j,i,ang) = (*var_cc)(kl,j,i,ang);
            } // endfor n
          } // endfor ifr
        } else { // modifications for polarization
          for (int m=0; m<nstok; ++m) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              for (int n=0; n<nang; ++n) {
                int ang=ifr*nang+n;
                (*var_cc)(kl-k,j,i,m,ang) = (*var_cc)(kl,j,i,m,ang);
              } // endfor n
            } // endfor ifr
          } // endfor m
        } // endelse !prad->use_pol_rad
      } // endfor i
    } // endfor j
  } // endfor k

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::OutflowInnerX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief Outflow boundary conditions, outer x3 boundary

void RadBoundaryVariable::OutflowOuterX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  NRRadiation *prad = pmy_block_->pnrrad;
  const int& nang  = prad->nang; // angles per octant
  const int& nfreq = prad->nfreq; // number of frequency bands
  const int& nstok = prad->num_stokes; // number of stokes parameters

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        if (!prad->use_pol_rad) {
          for (int ifr=0; ifr<nfreq; ++ifr) {
            for (int n=0; n<nang; ++n) {
              int ang=ifr*nang+n;
              (*var_cc)(ku+k,j,i,ang) = (*var_cc)(ku,j,i,ang);
            } // endfor n
          } // endfor ifr
        } else { // modifications for polarization
          for (int m=0; m<nstok; ++m) {
            for (int ifr=0; ifr<nfreq; ++ifr) {
              for (int n=0; n<nang; ++n) {
                int ang=ifr*nang+n;
                (*var_cc)(ku+k,j,i,m,ang) = (*var_cc)(ku,j,i,m,ang);
              } // endfor n
            } // endfor ifr
          } // endfor m
        } // endelse !prad->use_pol_rad
      } // endfor i
    } // endfor j
  } // endfor k

  return;
}
