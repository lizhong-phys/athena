//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file absorption_scattering.cpp
//  \brief implementation of absorption-scattering source terms
//======================================================================================

// C headers

// C++ headers
#include <cmath>

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../eos/eos.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../utils/utils.hpp"
#include "../../radiation.hpp"

// this class header
#include "../rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::CalPolAux()
//  \brief

// Note that coefficent matrix M_coeff and vector A are independent of Stokes parameters,
// and therefore can be computed very accurately

void RadIntegrator::CalPolAux(
      AthenaArray<Real> &wmu_cm, AthenaArray<Real> &tran_coef,
      AthenaArray<Real> &nx_cm, AthenaArray<Real> &ny_cm, AthenaArray<Real> &nz_cm,
      AthenaArray<Real> &wmu_cm_lbd, AthenaArray<Real> &tran_coef_lbd,
      AthenaArray<Real> &nx_cm_lbd, AthenaArray<Real> &ny_cm_lbd, AthenaArray<Real> &nz_cm_lbd,
      Real *sigma_a, Real *sigma_p, Real *sigma_pe, Real *sigma_s,
      Real dt, AthenaArray<Real> &ir_cm) {

  // parameters
  bool& refine_pol_coeff_ = pmy_rad->refine_pol_coeff;
  const int& nang_lbd = pmy_rad->nang_lbd;
  const int& nang  = pmy_rad->nang;
  const int& nfreq = pmy_rad->nfreq;
  Real cdt = dt * pmy_rad->crat;
  cdt *= pmy_rad->reduced_c/pmy_rad->crat;

  // reset M_coeff_, polVecA_ and polVecB_
  for (int m=0; m<23; ++m) {
    polVecA_(m) = 0.0;
    polVecB_(m) = 0.0;
    for (int n; n<23; ++n) {
      M_coeff_(m,n) = 0.0;
    }
  }

  // compute coefficents
  Real *wmu = &(wmu_cm(0));
  Real *nx  = &(nx_cm(0));
  Real *ny  = &(ny_cm(0));
  Real *nz  = &(nz_cm(0));
  Real *coef_l = &(tran_coef(0));
  Real *wmu_lbd = &(wmu_cm_lbd(0));
  Real *nx_lbd  = &(nx_cm_lbd(0));
  Real *ny_lbd  = &(ny_cm_lbd(0));
  Real *nz_lbd  = &(nz_cm_lbd(0));
  Real *coef_l_lbd = &(tran_coef_lbd(0));
  for (int ifr=0; ifr<nfreq; ++ifr) {
    // opacity
    Real chi_r = sigma_a[ifr];  // Rosseland mean absorption
    Real chi_s = sigma_s[ifr];  // Scattering
    Real chi_f = chi_s + chi_r;

    // stokes parameters
    Real *si_cm  = &(ir_cm(0,nang*ifr));
    Real *sq_cm  = &(ir_cm(1,nang*ifr));
    Real *su_cm  = &(ir_cm(2,nang*ifr));
    Real *sv_cm  = &(ir_cm(3,nang*ifr));

    // unique moments in vector A and matrix M_coeff
    Real hol_A     = 0.0;
    Real hol_Axx   = 0.0; Real hol_Ayy   = 0.0; Real hol_Azz   = 0.0;
    Real hol_Axy   = 0.0; Real hol_Axz   = 0.0; Real hol_Ayz   = 0.0;
    Real hol_Axxxx = 0.0; Real hol_Ayyyy = 0.0; Real hol_Azzzz = 0.0;
    Real hol_Axxyy = 0.0; Real hol_Axxzz = 0.0; Real hol_Ayyzz = 0.0;
    Real hol_Axxxy = 0.0; Real hol_Axxxz = 0.0;
    Real hol_Axyyy = 0.0; Real hol_Ayyyz = 0.0;
    Real hol_Axzzz = 0.0; Real hol_Ayzzz = 0.0;
    Real hol_Axxyz = 0.0; Real hol_Axyyz = 0.0; Real hol_Axyzz = 0.0;
    Real hol_Ac    = 0.0; Real hol_As    = 0.0;
    Real hol_Acc   = 0.0; Real hol_Ass   = 0.0; Real hol_Acs   = 0.0;
    Real hol_Axxc  = 0.0; Real hol_Ayyc  = 0.0; Real hol_Azzc  = 0.0;
    Real hol_Axxs  = 0.0; Real hol_Ayys  = 0.0; Real hol_Azzs  = 0.0;
    Real hol_Axyc  = 0.0; Real hol_Axzc  = 0.0; Real hol_Ayzc  = 0.0;
    Real hol_Axys  = 0.0; Real hol_Axzs  = 0.0; Real hol_Ayzs  = 0.0;
    Real hol_Azzcc = 0.0; Real hol_Azzcs = 0.0; Real hol_Azzss = 0.0;

    // angle integration
    for (int n=0; n<nang; n++) {
      Real fac = 1. / (1.0 + coef_l[n]*cdt*chi_f);
      Real fac_bi = fac * wmu[n] * si_cm[n];
      Real fac_bq = fac * wmu[n] * sq_cm[n];
      Real fac_bu = fac * wmu[n] * su_cm[n];
      Real fac_bv = fac * wmu[n] * sv_cm[n];
      Real sum_nx2ny2 = SQR(nx[n]) + SQR(ny[n]);
      Real nc = 1.0; Real ns = 0.0; // corresponding to phi=0 along the pole
      if (sum_nx2ny2 > 0) {
        nc  = (SQR(nx[n]) - SQR(ny[n])) / sum_nx2ny2;
        ns  = 2 * nx[n] * ny[n] / sum_nx2ny2;
      }
      Real nc_adhoc = (sum_nx2ny2==0) ? 0.0 : nc; // use nc=0 along the pole to minimize truncation errors in Ac and Azzc

      // Vector B for 23-moment closure:
      // (MJI, MKI11, MKI22, MKI33, MKI12, MKI13, MKI23)
      // (BI,  BI11,  BI22,  BI33,  BI12,  BI13,  BI23 )
      polVecB_(MJI)   += fac_bi;
      polVecB_(MKI11) += fac_bi * nx[n] * nx[n];
      polVecB_(MKI22) += fac_bi * ny[n] * ny[n];
      polVecB_(MKI33) += fac_bi * nz[n] * nz[n];
      polVecB_(MKI12) += fac_bi * nx[n] * ny[n];
      polVecB_(MKI13) += fac_bi * nx[n] * nz[n];
      polVecB_(MKI23) += fac_bi * ny[n] * nz[n];
      // (MJQ, MKQ11, MKQ22, MKQ33, MKQ12, MKQ13, MKQ23, MPQC, MPQS)
      // (BQ,  BQ11,  BQ22,  BQ33,  BQ12,  BQ13,  BQ23,  BQC,  BQS )
      polVecB_(MJQ)   += fac_bq;
      polVecB_(MKQ11) += fac_bq * nx[n] * nx[n];
      polVecB_(MKQ22) += fac_bq * ny[n] * ny[n];
      polVecB_(MKQ33) += fac_bq * nz[n] * nz[n];
      polVecB_(MKQ12) += fac_bq * nx[n] * ny[n];
      polVecB_(MKQ13) += fac_bq * nx[n] * nz[n];
      polVecB_(MKQ23) += fac_bq * ny[n] * nz[n];
      polVecB_(MPQC)  += fac_bq * nc;
      polVecB_(MPQS)  += fac_bq * ns;
      // (MHU1, MHU2, MPUZC, MPUZS)
      // (BU1,  BU2,  BUZC,  BUZS )
      polVecB_(MHU1)  += fac_bu * nx[n];
      polVecB_(MHU2)  += fac_bu * ny[n];
      polVecB_(MPUZC) += fac_bu * nz[n] * nc;
      polVecB_(MPUZS) += fac_bu * nz[n] * ns;
      // (MHV1, MHV2, MHV3)
      // (BV1,  BV2,  BV3 )
      polVecB_(MHV1)  += fac_bv * nx[n];
      polVecB_(MHV2)  += fac_bv * ny[n];
      polVecB_(MHV3)  += fac_bv * nz[n];

      // Compute unique moments in vector A and matrix M_coeff
      if (!refine_pol_coeff_) {
        Real fac_a = fac * coef_l[n]*wmu[n]*cdt;
        hol_A     += fac_a;
        hol_Axx   += fac_a * nx[n] * nx[n];
        hol_Ayy   += fac_a * ny[n] * ny[n];
        hol_Azz   += fac_a * nz[n] * nz[n];
        hol_Axy   += fac_a * nx[n] * ny[n];
        hol_Axz   += fac_a * nx[n] * nz[n];
        hol_Ayz   += fac_a * ny[n] * nz[n];
        hol_Axxxx += fac_a * nx[n] * nx[n] * nx[n] * nx[n];
        hol_Ayyyy += fac_a * ny[n] * ny[n] * ny[n] * ny[n];
        hol_Azzzz += fac_a * nz[n] * nz[n] * nz[n] * nz[n];
        hol_Axxyy += fac_a * nx[n] * nx[n] * ny[n] * ny[n];
        hol_Axxzz += fac_a * nx[n] * nx[n] * nz[n] * nz[n];
        hol_Ayyzz += fac_a * ny[n] * ny[n] * nz[n] * nz[n];
        hol_Axxxy += fac_a * nx[n] * nx[n] * nx[n] * ny[n];
        hol_Axxxz += fac_a * nx[n] * nx[n] * nx[n] * nz[n];
        hol_Axyyy += fac_a * nx[n] * ny[n] * ny[n] * ny[n];
        hol_Ayyyz += fac_a * ny[n] * ny[n] * ny[n] * nz[n];
        hol_Axzzz += fac_a * nx[n] * nz[n] * nz[n] * nz[n];
        hol_Ayzzz += fac_a * ny[n] * nz[n] * nz[n] * nz[n];
        hol_Axxyz += fac_a * nx[n] * nx[n] * ny[n] * nz[n];
        hol_Axyyz += fac_a * nx[n] * ny[n] * ny[n] * nz[n];
        hol_Axyzz += fac_a * nx[n] * ny[n] * nz[n] * nz[n];
        hol_Ac    += fac_a * nc_adhoc;
        hol_As    += fac_a * ns;
        hol_Acc   += fac_a * nc    * nc;
        hol_Ass   += fac_a * ns    * ns;
        hol_Acs   += fac_a * nc    * ns;
        hol_Axxc  += fac_a * nx[n] * nx[n] * nc;
        hol_Ayyc  += fac_a * ny[n] * ny[n] * nc;
        hol_Azzc  += fac_a * nz[n] * nz[n] * nc_adhoc;
        hol_Axxs  += fac_a * nx[n] * nx[n] * ns;
        hol_Ayys  += fac_a * ny[n] * ny[n] * ns;
        hol_Azzs  += fac_a * nz[n] * nz[n] * ns;
        hol_Axyc  += fac_a * nx[n] * ny[n] * nc;
        hol_Axzc  += fac_a * nx[n] * nz[n] * nc;
        hol_Ayzc  += fac_a * ny[n] * nz[n] * nc;
        hol_Axys  += fac_a * nx[n] * ny[n] * ns;
        hol_Axzs  += fac_a * nx[n] * nz[n] * ns;
        hol_Ayzs  += fac_a * ny[n] * nz[n] * ns;
        hol_Azzcc += fac_a * nz[n] * nz[n] * nc    * nc;
        hol_Azzcs += fac_a * nz[n] * nz[n] * nc    * ns;
        hol_Azzss += fac_a * nz[n] * nz[n] * ns    * ns;
      } // endif (!refine_pol_coeff_)
    } // endfor n

    // Compute unique moments in vector A and matrix M_coeff using Lebedev quadrature
    if (refine_pol_coeff_) {
      for (int n=0; n<nang_lbd; n++) {
        Real sum_nx2ny2 = SQR(nx_lbd[n]) + SQR(ny_lbd[n]);
        Real nc_lbd = 1.0; Real ns_lbd = 0.0; // corresponding to phi=0 along the pole
        if (sum_nx2ny2 > 0) {
          nc_lbd = (SQR(nx_lbd[n]) - SQR(ny_lbd[n])) / sum_nx2ny2;
          ns_lbd = 2 * nx_lbd[n] * ny_lbd[n] / sum_nx2ny2;
        }
        Real nc_adhoc_lbd = (sum_nx2ny2==0) ? 0.0 : nc_lbd; // use nc=0 along the pole to minimize truncation errors in Ac and Azzc

        Real fac = 1. / (1.0 + coef_l_lbd[n]*cdt*chi_f);
        Real fac_a = fac * coef_l_lbd[n]*wmu_lbd[n]*cdt;
        hol_A     += fac_a;
        hol_Axx   += fac_a * nx_lbd[n] * nx_lbd[n];
        hol_Ayy   += fac_a * ny_lbd[n] * ny_lbd[n];
        hol_Azz   += fac_a * nz_lbd[n] * nz_lbd[n];
        hol_Axy   += fac_a * nx_lbd[n] * ny_lbd[n];
        hol_Axz   += fac_a * nx_lbd[n] * nz_lbd[n];
        hol_Ayz   += fac_a * ny_lbd[n] * nz_lbd[n];
        hol_Axxxx += fac_a * nx_lbd[n] * nx_lbd[n] * nx_lbd[n] * nx_lbd[n];
        hol_Ayyyy += fac_a * ny_lbd[n] * ny_lbd[n] * ny_lbd[n] * ny_lbd[n];
        hol_Azzzz += fac_a * nz_lbd[n] * nz_lbd[n] * nz_lbd[n] * nz_lbd[n];
        hol_Axxyy += fac_a * nx_lbd[n] * nx_lbd[n] * ny_lbd[n] * ny_lbd[n];
        hol_Axxzz += fac_a * nx_lbd[n] * nx_lbd[n] * nz_lbd[n] * nz_lbd[n];
        hol_Ayyzz += fac_a * ny_lbd[n] * ny_lbd[n] * nz_lbd[n] * nz_lbd[n];
        hol_Axxxy += fac_a * nx_lbd[n] * nx_lbd[n] * nx_lbd[n] * ny_lbd[n];
        hol_Axxxz += fac_a * nx_lbd[n] * nx_lbd[n] * nx_lbd[n] * nz_lbd[n];
        hol_Axyyy += fac_a * nx_lbd[n] * ny_lbd[n] * ny_lbd[n] * ny_lbd[n];
        hol_Ayyyz += fac_a * ny_lbd[n] * ny_lbd[n] * ny_lbd[n] * nz_lbd[n];
        hol_Axzzz += fac_a * nx_lbd[n] * nz_lbd[n] * nz_lbd[n] * nz_lbd[n];
        hol_Ayzzz += fac_a * ny_lbd[n] * nz_lbd[n] * nz_lbd[n] * nz_lbd[n];
        hol_Axxyz += fac_a * nx_lbd[n] * nx_lbd[n] * ny_lbd[n] * nz_lbd[n];
        hol_Axyyz += fac_a * nx_lbd[n] * ny_lbd[n] * ny_lbd[n] * nz_lbd[n];
        hol_Axyzz += fac_a * nx_lbd[n] * ny_lbd[n] * nz_lbd[n] * nz_lbd[n];
        hol_Ac    += fac_a * nc_adhoc_lbd;
        hol_As    += fac_a * ns_lbd;
        hol_Acc   += fac_a * nc_lbd    * nc_lbd;
        hol_Ass   += fac_a * ns_lbd    * ns_lbd;
        hol_Acs   += fac_a * nc_lbd    * ns_lbd;
        hol_Axxc  += fac_a * nx_lbd[n] * nx_lbd[n] * nc_lbd;
        hol_Ayyc  += fac_a * ny_lbd[n] * ny_lbd[n] * nc_lbd;
        hol_Azzc  += fac_a * nz_lbd[n] * nz_lbd[n] * nc_adhoc_lbd;
        hol_Axxs  += fac_a * nx_lbd[n] * nx_lbd[n] * ns_lbd;
        hol_Ayys  += fac_a * ny_lbd[n] * ny_lbd[n] * ns_lbd;
        hol_Azzs  += fac_a * nz_lbd[n] * nz_lbd[n] * ns_lbd;
        hol_Axyc  += fac_a * nx_lbd[n] * ny_lbd[n] * nc_lbd;
        hol_Axzc  += fac_a * nx_lbd[n] * nz_lbd[n] * nc_lbd;
        hol_Ayzc  += fac_a * ny_lbd[n] * nz_lbd[n] * nc_lbd;
        hol_Axys  += fac_a * nx_lbd[n] * ny_lbd[n] * ns_lbd;
        hol_Axzs  += fac_a * nx_lbd[n] * nz_lbd[n] * ns_lbd;
        hol_Ayzs  += fac_a * ny_lbd[n] * nz_lbd[n] * ns_lbd;
        hol_Azzcc += fac_a * nz_lbd[n] * nz_lbd[n] * nc_lbd    * nc_lbd;
        hol_Azzcs += fac_a * nz_lbd[n] * nz_lbd[n] * nc_lbd    * ns_lbd;
        hol_Azzss += fac_a * nz_lbd[n] * nz_lbd[n] * ns_lbd    * ns_lbd;
      } // endfor n
    } // endelse (refine_pol_coeff_)

    // Vector A for 23-moment closure:
    // (MJI, MKI11, MKI22, MKI33, MKI12, MKI13, MKI23)
    // (A,   A_xx,  A_yy,  A_zz,  A_xy,  A_xz,  A_yz)
    polVecA_(MJI)   = hol_A;
    polVecA_(MKI11) = hol_Axx;
    polVecA_(MKI22) = hol_Ayy;
    polVecA_(MKI33) = hol_Azz;
    polVecA_(MKI12) = hol_Axy;
    polVecA_(MKI13) = hol_Axz;
    polVecA_(MKI23) = hol_Ayz;

    // Matrix M_coeff for 23-moment closure:
    // JI moment equation
    M_coeff_(MJI,MJI)   = hol_A;
    M_coeff_(MJI,MKI11) = hol_Axx;
    M_coeff_(MJI,MKI22) = hol_Ayy;
    M_coeff_(MJI,MKI33) = hol_Azz;
    M_coeff_(MJI,MKI12) = hol_Axy * 2;
    M_coeff_(MJI,MKI13) = hol_Axz * 2;
    M_coeff_(MJI,MKI23) = hol_Ayz * 2;
    M_coeff_(MJI,MJQ)   = hol_Azz * (-1);
    M_coeff_(MJI,MKQ11) = hol_Axx;
    M_coeff_(MJI,MKQ22) = hol_Ayy;
    M_coeff_(MJI,MKQ33) = hol_Azz;
    M_coeff_(MJI,MKQ12) = hol_Axy * 2;
    M_coeff_(MJI,MKQ13) = hol_Axz * 2;
    M_coeff_(MJI,MKQ23) = hol_Ayz * 2;
    M_coeff_(MJI,MPQC)  = hol_Ayy - hol_Axx;
    M_coeff_(MJI,MPQS)  = hol_Axy * (-2);
    M_coeff_(MJI,MHU1)  = hol_Ayz * 2;
    M_coeff_(MJI,MHU2)  = hol_Axz * (-2);
    M_coeff_(MJI,MPUZC) = hol_Axy * (-2);
    M_coeff_(MJI,MPUZS) = hol_Axx - hol_Ayy;

    // KI11 moment equation
    M_coeff_(MKI11,MJI)   = hol_Axx;
    M_coeff_(MKI11,MKI11) = hol_Axxxx;
    M_coeff_(MKI11,MKI22) = hol_Axxyy;
    M_coeff_(MKI11,MKI33) = hol_Axxzz;
    M_coeff_(MKI11,MKI12) = hol_Axxxy * 2;
    M_coeff_(MKI11,MKI13) = hol_Axxxz * 2;
    M_coeff_(MKI11,MKI23) = hol_Axxyz * 2;
    M_coeff_(MKI11,MJQ)   = hol_Axxzz * (-1);
    M_coeff_(MKI11,MKQ11) = hol_Axxxx;
    M_coeff_(MKI11,MKQ22) = hol_Axxyy;
    M_coeff_(MKI11,MKQ33) = hol_Axxzz;
    M_coeff_(MKI11,MKQ12) = hol_Axxxy * 2;
    M_coeff_(MKI11,MKQ13) = hol_Axxxz * 2;
    M_coeff_(MKI11,MKQ23) = hol_Axxyz * 2;
    M_coeff_(MKI11,MPQC)  = hol_Axxyy - hol_Axxxx;
    M_coeff_(MKI11,MPQS)  = hol_Axxxy * (-2);
    M_coeff_(MKI11,MHU1)  = hol_Axxyz * 2;
    M_coeff_(MKI11,MHU2)  = hol_Axxxz * (-2);
    M_coeff_(MKI11,MPUZC) = hol_Axxxy * (-2);
    M_coeff_(MKI11,MPUZS) = hol_Axxxx - hol_Axxyy;

    // KI22 moment equation
    M_coeff_(MKI22,MJI)   = hol_Ayy;
    M_coeff_(MKI22,MKI11) = hol_Axxyy;
    M_coeff_(MKI22,MKI22) = hol_Ayyyy;
    M_coeff_(MKI22,MKI33) = hol_Ayyzz;
    M_coeff_(MKI22,MKI12) = hol_Axyyy * 2;
    M_coeff_(MKI22,MKI13) = hol_Axyyz * 2;
    M_coeff_(MKI22,MKI23) = hol_Ayyyz * 2;
    M_coeff_(MKI22,MJQ)   = hol_Ayyzz * (-1);
    M_coeff_(MKI22,MKQ11) = hol_Axxyy;
    M_coeff_(MKI22,MKQ22) = hol_Ayyyy;
    M_coeff_(MKI22,MKQ33) = hol_Ayyzz;
    M_coeff_(MKI22,MKQ12) = hol_Axyyy * 2;
    M_coeff_(MKI22,MKQ13) = hol_Axyyz * 2;
    M_coeff_(MKI22,MKQ23) = hol_Ayyyz * 2;
    M_coeff_(MKI22,MPQC)  = hol_Ayyyy - hol_Axxyy;
    M_coeff_(MKI22,MPQS)  = hol_Axyyy * (-2);
    M_coeff_(MKI22,MHU1)  = hol_Ayyyz * 2;
    M_coeff_(MKI22,MHU2)  = hol_Axyyz * (-2);
    M_coeff_(MKI22,MPUZC) = hol_Axyyy * (-2);
    M_coeff_(MKI22,MPUZS) = hol_Axxyy - hol_Ayyyy;

    // KI33 moment equation
    M_coeff_(MKI33,MJI)   = hol_Azz;
    M_coeff_(MKI33,MKI11) = hol_Axxzz;
    M_coeff_(MKI33,MKI22) = hol_Ayyzz;
    M_coeff_(MKI33,MKI33) = hol_Azzzz;
    M_coeff_(MKI33,MKI12) = hol_Axyzz * 2;
    M_coeff_(MKI33,MKI13) = hol_Axzzz * 2;
    M_coeff_(MKI33,MKI23) = hol_Ayzzz * 2;
    M_coeff_(MKI33,MJQ)   = hol_Azzzz * (-1);
    M_coeff_(MKI33,MKQ11) = hol_Axxzz;
    M_coeff_(MKI33,MKQ22) = hol_Ayyzz;
    M_coeff_(MKI33,MKQ33) = hol_Azzzz;
    M_coeff_(MKI33,MKQ12) = hol_Axyzz * 2;
    M_coeff_(MKI33,MKQ13) = hol_Axzzz * 2;
    M_coeff_(MKI33,MKQ23) = hol_Ayzzz * 2;
    M_coeff_(MKI33,MPQC)  = hol_Ayyzz - hol_Axxzz;
    M_coeff_(MKI33,MPQS)  = hol_Axyzz * (-2);
    M_coeff_(MKI33,MHU1)  = hol_Ayzzz * 2;
    M_coeff_(MKI33,MHU2)  = hol_Axzzz * (-2);
    M_coeff_(MKI33,MPUZC) = hol_Axyzz * (-2);
    M_coeff_(MKI33,MPUZS) = hol_Axxzz - hol_Ayyzz;

    // KI12 moment equation
    M_coeff_(MKI12,MJI)   = hol_Axy;
    M_coeff_(MKI12,MKI11) = hol_Axxxy;
    M_coeff_(MKI12,MKI22) = hol_Axyyy;
    M_coeff_(MKI12,MKI33) = hol_Axyzz;
    M_coeff_(MKI12,MKI12) = hol_Axxyy * 2;
    M_coeff_(MKI12,MKI13) = hol_Axxyz * 2;
    M_coeff_(MKI12,MKI23) = hol_Axyyz * 2;
    M_coeff_(MKI12,MJQ)   = hol_Axyzz * (-1);
    M_coeff_(MKI12,MKQ11) = hol_Axxxy;
    M_coeff_(MKI12,MKQ22) = hol_Axyyy;
    M_coeff_(MKI12,MKQ33) = hol_Axyzz;
    M_coeff_(MKI12,MKQ12) = hol_Axxyy * 2;
    M_coeff_(MKI12,MKQ13) = hol_Axxyz * 2;
    M_coeff_(MKI12,MKQ23) = hol_Axyyz * 2;
    M_coeff_(MKI12,MPQC)  = hol_Axyyy - hol_Axxxy;
    M_coeff_(MKI12,MPQS)  = hol_Axxyy * (-2);
    M_coeff_(MKI12,MHU1)  = hol_Axyyz * 2;
    M_coeff_(MKI12,MHU2)  = hol_Axxyz * (-2);
    M_coeff_(MKI12,MPUZC) = hol_Axxyy * (-2);
    M_coeff_(MKI12,MPUZS) = hol_Axxxy - hol_Axyyy;

    // KI13 moment equation
    M_coeff_(MKI13,MJI)   = hol_Axz;
    M_coeff_(MKI13,MKI11) = hol_Axxxz;
    M_coeff_(MKI13,MKI22) = hol_Axyyz;
    M_coeff_(MKI13,MKI33) = hol_Axzzz;
    M_coeff_(MKI13,MKI12) = hol_Axxyz * 2;
    M_coeff_(MKI13,MKI13) = hol_Axxzz * 2;
    M_coeff_(MKI13,MKI23) = hol_Axyzz * 2;
    M_coeff_(MKI13,MJQ)   = hol_Axzzz * (-1);
    M_coeff_(MKI13,MKQ11) = hol_Axxxz;
    M_coeff_(MKI13,MKQ22) = hol_Axyyz;
    M_coeff_(MKI13,MKQ33) = hol_Axzzz;
    M_coeff_(MKI13,MKQ12) = hol_Axxyz * 2;
    M_coeff_(MKI13,MKQ13) = hol_Axxzz * 2;
    M_coeff_(MKI13,MKQ23) = hol_Axyzz * 2;
    M_coeff_(MKI13,MPQC)  = hol_Axyyz - hol_Axxxz;
    M_coeff_(MKI13,MPQS)  = hol_Axxyz * (-2);
    M_coeff_(MKI13,MHU1)  = hol_Axyzz * 2;
    M_coeff_(MKI13,MHU2)  = hol_Axxzz * (-2);
    M_coeff_(MKI13,MPUZC) = hol_Axxyz * (-2);
    M_coeff_(MKI13,MPUZS) = hol_Axxxz - hol_Axyyz;

    // KI23 moment equation
    M_coeff_(MKI23,MJI)   = hol_Ayz;
    M_coeff_(MKI23,MKI11) = hol_Axxyz;
    M_coeff_(MKI23,MKI22) = hol_Ayyyz;
    M_coeff_(MKI23,MKI33) = hol_Ayzzz;
    M_coeff_(MKI23,MKI12) = hol_Axyyz * 2;
    M_coeff_(MKI23,MKI13) = hol_Axyzz * 2;
    M_coeff_(MKI23,MKI23) = hol_Ayyzz * 2;
    M_coeff_(MKI23,MJQ)   = hol_Ayzzz * (-1);
    M_coeff_(MKI23,MKQ11) = hol_Axxyz;
    M_coeff_(MKI23,MKQ22) = hol_Ayyyz;
    M_coeff_(MKI23,MKQ33) = hol_Ayzzz;
    M_coeff_(MKI23,MKQ12) = hol_Axyyz * 2;
    M_coeff_(MKI23,MKQ13) = hol_Axyzz * 2;
    M_coeff_(MKI23,MKQ23) = hol_Ayyzz * 2;
    M_coeff_(MKI23,MPQC)  = hol_Ayyyz - hol_Axxyz;
    M_coeff_(MKI23,MPQS)  = hol_Axyyz * (-2);
    M_coeff_(MKI23,MHU1)  = hol_Ayyzz * 2;
    M_coeff_(MKI23,MHU2)  = hol_Axyzz * (-2);
    M_coeff_(MKI23,MPUZC) = hol_Axyyz * (-2);
    M_coeff_(MKI23,MPUZS) = hol_Axxyz - hol_Ayyyz;

    // JQ moment equation
    M_coeff_(MJQ,MKI11) =  hol_Axx - hol_Ac;
    M_coeff_(MJQ,MKI22) =  hol_Ayy + hol_Ac;
    M_coeff_(MJQ,MKI33) =  hol_Azz - hol_A;
    M_coeff_(MJQ,MKI12) = (hol_Axy - hol_As) * 2;
    M_coeff_(MJQ,MKI13) =  hol_Axz * 2;
    M_coeff_(MJQ,MKI23) =  hol_Ayz * 2;
    M_coeff_(MJQ,MJQ)   =  hol_A - hol_Azz;
    M_coeff_(MJQ,MKQ11) =  hol_Axx - hol_Ac;
    M_coeff_(MJQ,MKQ22) =  hol_Ayy + hol_Ac;
    M_coeff_(MJQ,MKQ33) =  hol_Azz - hol_A;
    M_coeff_(MJQ,MKQ12) = (hol_Axy - hol_As) * 2;
    M_coeff_(MJQ,MKQ13) =  hol_Axz * 2;
    M_coeff_(MJQ,MKQ23) =  hol_Ayz * 2;
    M_coeff_(MJQ,MPQC)  = (2*hol_Ac - hol_Axx + hol_Ayy);
    M_coeff_(MJQ,MPQS)  = (hol_As - hol_Axy) * 2;
    M_coeff_(MJQ,MHU1)  =  hol_Ayz * 2;
    M_coeff_(MJQ,MHU2)  =  hol_Axz * (-2);
    M_coeff_(MJQ,MPUZC) = (hol_As - hol_Axy) * 2;
    M_coeff_(MJQ,MPUZS) = -(2*hol_Ac - hol_Axx + hol_Ayy);

    // KQ11 moment equation
    M_coeff_(MKQ11,MKI11) =  hol_Axxxx - hol_Axxc;
    M_coeff_(MKQ11,MKI22) =  hol_Axxyy + hol_Axxc;
    M_coeff_(MKQ11,MKI33) =  hol_Axxzz - hol_Axx;
    M_coeff_(MKQ11,MKI12) = (hol_Axxxy - hol_Axxs) * 2;
    M_coeff_(MKQ11,MKI13) =  hol_Axxxz * 2;
    M_coeff_(MKQ11,MKI23) =  hol_Axxyz * 2;
    M_coeff_(MKQ11,MJQ)   =  hol_Axx - hol_Axxzz;
    M_coeff_(MKQ11,MKQ11) =  hol_Axxxx - hol_Axxc;
    M_coeff_(MKQ11,MKQ22) =  hol_Axxyy + hol_Axxc;
    M_coeff_(MKQ11,MKQ33) =  hol_Axxzz - hol_Axx;
    M_coeff_(MKQ11,MKQ12) = (hol_Axxxy - hol_Axxs) * 2;
    M_coeff_(MKQ11,MKQ13) =  hol_Axxxz * 2;
    M_coeff_(MKQ11,MKQ23) =  hol_Axxyz * 2;
    M_coeff_(MKQ11,MPQC)  = (2*hol_Axxc - hol_Axxxx + hol_Axxyy);
    M_coeff_(MKQ11,MPQS)  = (hol_Axxs - hol_Axxxy) * 2;
    M_coeff_(MKQ11,MHU1)  =  hol_Axxyz * 2;
    M_coeff_(MKQ11,MHU2)  =  hol_Axxxz * (-2);
    M_coeff_(MKQ11,MPUZC) = (hol_Axxs - hol_Axxxy) * 2;
    M_coeff_(MKQ11,MPUZS) = -(2*hol_Axxc - hol_Axxxx + hol_Axxyy);

    // KQ22 moment equation
    M_coeff_(MKQ22,MKI11) =  hol_Axxyy - hol_Ayyc;
    M_coeff_(MKQ22,MKI22) =  hol_Ayyyy + hol_Ayyc;
    M_coeff_(MKQ22,MKI33) =  hol_Ayyzz - hol_Ayy;
    M_coeff_(MKQ22,MKI12) = (hol_Axyyy - hol_Ayys) * 2;
    M_coeff_(MKQ22,MKI13) =  hol_Axyyz * 2;
    M_coeff_(MKQ22,MKI23) =  hol_Ayyyz * 2;
    M_coeff_(MKQ22,MJQ)   =  hol_Ayy - hol_Ayyzz;
    M_coeff_(MKQ22,MKQ11) =  hol_Axxyy - hol_Ayyc;
    M_coeff_(MKQ22,MKQ22) =  hol_Ayyyy + hol_Ayyc;
    M_coeff_(MKQ22,MKQ33) =  hol_Ayyzz - hol_Ayy;
    M_coeff_(MKQ22,MKQ12) = (hol_Axyyy - hol_Ayys) * 2;
    M_coeff_(MKQ22,MKQ13) =  hol_Axyyz * 2;
    M_coeff_(MKQ22,MKQ23) =  hol_Ayyyz * 2;
    M_coeff_(MKQ22,MPQC)  = (2*hol_Ayyc - hol_Axxyy + hol_Ayyyy);
    M_coeff_(MKQ22,MPQS)  = (hol_Ayys - hol_Axyyy) * 2;
    M_coeff_(MKQ22,MHU1)  =  hol_Ayyyz * 2;
    M_coeff_(MKQ22,MHU2)  =  hol_Axyyz * (-2);
    M_coeff_(MKQ22,MPUZC) = (hol_Ayys - hol_Axyyy) * 2;
    M_coeff_(MKQ22,MPUZS) = -(2*hol_Ayyc - hol_Axxyy + hol_Ayyyy);

    // KQ33 moment equation
    M_coeff_(MKQ33,MKI11) =  hol_Axxzz - hol_Azzc;
    M_coeff_(MKQ33,MKI22) =  hol_Ayyzz + hol_Azzc;
    M_coeff_(MKQ33,MKI33) =  hol_Azzzz - hol_Azz;
    M_coeff_(MKQ33,MKI12) = (hol_Axyzz - hol_Azzs) * 2;
    M_coeff_(MKQ33,MKI13) =  hol_Axzzz * 2;
    M_coeff_(MKQ33,MKI23) =  hol_Ayzzz * 2;
    M_coeff_(MKQ33,MJQ)   =  hol_Azz - hol_Azzzz;
    M_coeff_(MKQ33,MKQ11) =  hol_Axxzz - hol_Azzc;
    M_coeff_(MKQ33,MKQ22) =  hol_Ayyzz + hol_Azzc;
    M_coeff_(MKQ33,MKQ33) =  hol_Azzzz - hol_Azz;
    M_coeff_(MKQ33,MKQ12) = (hol_Axyzz - hol_Azzs) * 2;
    M_coeff_(MKQ33,MKQ13) =  hol_Axzzz * 2;
    M_coeff_(MKQ33,MKQ23) =  hol_Ayzzz * 2;
    M_coeff_(MKQ33,MPQC)  = (2*hol_Azzc - hol_Axxzz + hol_Ayyzz);
    M_coeff_(MKQ33,MPQS)  = (hol_Azzs - hol_Axyzz) * 2;
    M_coeff_(MKQ33,MHU1)  =  hol_Ayzzz * 2;
    M_coeff_(MKQ33,MHU2)  =  hol_Axzzz * (-2);
    M_coeff_(MKQ33,MPUZC) = (hol_Azzs - hol_Axyzz) * 2;
    M_coeff_(MKQ33,MPUZS) = -(2*hol_Azzc - hol_Axxzz + hol_Ayyzz);

    // KQ12 moment equation
    M_coeff_(MKQ12,MKI11) =  hol_Axxxy - hol_Axyc;
    M_coeff_(MKQ12,MKI22) =  hol_Axyyy + hol_Axyc;
    M_coeff_(MKQ12,MKI33) =  hol_Axyzz - hol_Axy;
    M_coeff_(MKQ12,MKI12) = (hol_Axxyy - hol_Axys) * 2;
    M_coeff_(MKQ12,MKI13) =  hol_Axxyz * 2;
    M_coeff_(MKQ12,MKI23) =  hol_Axyyz * 2;
    M_coeff_(MKQ12,MJQ)   =  hol_Axy - hol_Axyzz;
    M_coeff_(MKQ12,MKQ11) =  hol_Axxxy - hol_Axyc;
    M_coeff_(MKQ12,MKQ22) =  hol_Axyyy + hol_Axyc;
    M_coeff_(MKQ12,MKQ33) =  hol_Axyzz - hol_Axy;
    M_coeff_(MKQ12,MKQ12) = (hol_Axxyy - hol_Axys) * 2;
    M_coeff_(MKQ12,MKQ13) =  hol_Axxyz * 2;
    M_coeff_(MKQ12,MKQ23) =  hol_Axyyz * 2;
    M_coeff_(MKQ12,MPQC)  = (2*hol_Axyc - hol_Axxxy + hol_Axyyy);
    M_coeff_(MKQ12,MPQS)  = (hol_Axys - hol_Axxyy) * 2;
    M_coeff_(MKQ12,MHU1)  =  hol_Axyyz * 2;
    M_coeff_(MKQ12,MHU2)  =  hol_Axxyz * (-2);
    M_coeff_(MKQ12,MPUZC) = (hol_Axys - hol_Axxyy) * 2;
    M_coeff_(MKQ12,MPUZS) = -(2*hol_Axyc - hol_Axxxy + hol_Axyyy);

    // KQ13 moment equation
    M_coeff_(MKQ13,MKI11) =  hol_Axxxz - hol_Axzc;
    M_coeff_(MKQ13,MKI22) =  hol_Axyyz + hol_Axzc;
    M_coeff_(MKQ13,MKI33) =  hol_Axzzz - hol_Axz;
    M_coeff_(MKQ13,MKI12) = (hol_Axxyz - hol_Axzs) * 2;
    M_coeff_(MKQ13,MKI13) =  hol_Axxzz * 2;
    M_coeff_(MKQ13,MKI23) =  hol_Axyzz * 2;
    M_coeff_(MKQ13,MJQ)   =  hol_Axz - hol_Axzzz;
    M_coeff_(MKQ13,MKQ11) =  hol_Axxxz - hol_Axzc;
    M_coeff_(MKQ13,MKQ22) =  hol_Axyyz + hol_Axzc;
    M_coeff_(MKQ13,MKQ33) =  hol_Axzzz - hol_Axz;
    M_coeff_(MKQ13,MKQ12) = (hol_Axxyz - hol_Axzs) * 2;
    M_coeff_(MKQ13,MKQ13) =  hol_Axxzz * 2;
    M_coeff_(MKQ13,MKQ23) =  hol_Axyzz * 2;
    M_coeff_(MKQ13,MPQC)  = (2*hol_Axzc - hol_Axxxz + hol_Axyyz);
    M_coeff_(MKQ13,MPQS)  = (hol_Axzs - hol_Axxyz) * 2;
    M_coeff_(MKQ13,MHU1)  =  hol_Axyzz * 2;
    M_coeff_(MKQ13,MHU2)  =  hol_Axxzz * (-2);
    M_coeff_(MKQ13,MPUZC) = (hol_Axzs - hol_Axxyz) * 2;
    M_coeff_(MKQ13,MPUZS) = -(2*hol_Axzc - hol_Axxxz + hol_Axyyz);

    // KQ23 moment equation
    M_coeff_(MKQ23,MKI11) =  hol_Axxyz - hol_Ayzc;
    M_coeff_(MKQ23,MKI22) =  hol_Ayyyz + hol_Ayzc;
    M_coeff_(MKQ23,MKI33) =  hol_Ayzzz - hol_Ayz;
    M_coeff_(MKQ23,MKI12) = (hol_Axyyz - hol_Ayzs) * 2;
    M_coeff_(MKQ23,MKI13) =  hol_Axyzz * 2;
    M_coeff_(MKQ23,MKI23) =  hol_Ayyzz * 2;
    M_coeff_(MKQ23,MJQ)   =  hol_Ayz - hol_Ayzzz;
    M_coeff_(MKQ23,MKQ11) =  hol_Axxyz - hol_Ayzc;
    M_coeff_(MKQ23,MKQ22) =  hol_Ayyyz + hol_Ayzc;
    M_coeff_(MKQ23,MKQ33) =  hol_Ayzzz - hol_Ayz;
    M_coeff_(MKQ23,MKQ12) = (hol_Axyyz - hol_Ayzs) * 2;
    M_coeff_(MKQ23,MKQ13) =  hol_Axyzz * 2;
    M_coeff_(MKQ23,MKQ23) =  hol_Ayyzz * 2;
    M_coeff_(MKQ23,MPQC)  = (2*hol_Ayzc - hol_Axxyz + hol_Ayyyz);
    M_coeff_(MKQ23,MPQS)  = (hol_Ayzs - hol_Axyyz) * 2;
    M_coeff_(MKQ23,MHU1)  =  hol_Ayyzz * 2;
    M_coeff_(MKQ23,MHU2)  =  hol_Axyzz * (-2);
    M_coeff_(MKQ23,MPUZC) = (hol_Ayzs - hol_Axyyz) * 2;
    M_coeff_(MKQ23,MPUZS) = -(2*hol_Ayzc - hol_Axxyz + hol_Ayyyz);

    // PQC moment equation
    M_coeff_(MPQC,MKI11) =  hol_Axxc - hol_Acc;
    M_coeff_(MPQC,MKI22) =  hol_Ayyc + hol_Acc;
    M_coeff_(MPQC,MKI33) =  hol_Azzc - hol_Ac;
    M_coeff_(MPQC,MKI12) = (hol_Axyc - hol_Acs) * 2;
    M_coeff_(MPQC,MKI13) =  hol_Axzc * 2;
    M_coeff_(MPQC,MKI23) =  hol_Ayzc * 2;
    M_coeff_(MPQC,MJQ)   =  hol_Ac - hol_Azzc;
    M_coeff_(MPQC,MKQ11) =  hol_Axxc - hol_Acc;
    M_coeff_(MPQC,MKQ22) =  hol_Ayyc + hol_Acc;
    M_coeff_(MPQC,MKQ33) =  hol_Azzc - hol_Ac;
    M_coeff_(MPQC,MKQ12) = (hol_Axyc - hol_Acs) * 2;
    M_coeff_(MPQC,MKQ13) =  hol_Axzc * 2;
    M_coeff_(MPQC,MKQ23) =  hol_Ayzc * 2;
    M_coeff_(MPQC,MPQC)  = (2*hol_Acc - hol_Axxc + hol_Ayyc);
    M_coeff_(MPQC,MPQS)  = (hol_Acs - hol_Axyc) * 2;
    M_coeff_(MPQC,MHU1)  =  hol_Ayzc * 2;
    M_coeff_(MPQC,MHU2)  =  hol_Axzc * (-2);
    M_coeff_(MPQC,MPUZC) = (hol_Acs - hol_Axyc) * 2;
    M_coeff_(MPQC,MPUZS) = -(2*hol_Acc - hol_Axxc + hol_Ayyc);

    // PQS moment equation
    M_coeff_(MPQS,MKI11) =  hol_Axxs - hol_Acs;
    M_coeff_(MPQS,MKI22) =  hol_Ayys + hol_Acs;
    M_coeff_(MPQS,MKI33) =  hol_Azzs - hol_As;
    M_coeff_(MPQS,MKI12) = (hol_Axys - hol_Ass) * 2;
    M_coeff_(MPQS,MKI13) =  hol_Axzs * 2;
    M_coeff_(MPQS,MKI23) =  hol_Ayzs * 2;
    M_coeff_(MPQS,MJQ)   =  hol_As - hol_Azzs;
    M_coeff_(MPQS,MKQ11) =  hol_Axxs - hol_Acs;
    M_coeff_(MPQS,MKQ22) =  hol_Ayys + hol_Acs;
    M_coeff_(MPQS,MKQ33) =  hol_Azzs - hol_As;
    M_coeff_(MPQS,MKQ12) = (hol_Axys - hol_Ass) * 2;
    M_coeff_(MPQS,MKQ13) =  hol_Axzs * 2;
    M_coeff_(MPQS,MKQ23) =  hol_Ayzs * 2;
    M_coeff_(MPQS,MPQC)  = (2*hol_Acs - hol_Axxs + hol_Ayys);
    M_coeff_(MPQS,MPQS)  = (hol_Ass - hol_Axys) * 2;
    M_coeff_(MPQS,MHU1)  =  hol_Ayzs * 2;
    M_coeff_(MPQS,MHU2)  =  hol_Axzs * (-2);
    M_coeff_(MPQS,MPUZC) = (hol_Ass - hol_Axys) * 2;
    M_coeff_(MPQS,MPUZS) = -(2*hol_Acs - hol_Axxs + hol_Ayys);

    // HU1 moment equation
    M_coeff_(MHU1,MKI11) = hol_Axzs;
    M_coeff_(MHU1,MKI22) = hol_Axzs * (-1);
    M_coeff_(MHU1,MKI12) = hol_Axzc * (-2);
    M_coeff_(MHU1,MKI13) = hol_Axy  * (-2);
    M_coeff_(MHU1,MKI23) = hol_Axx  * 2;
    M_coeff_(MHU1,MKQ11) = hol_Axzs;
    M_coeff_(MHU1,MKQ22) = hol_Axzs * (-1);
    M_coeff_(MHU1,MKQ12) = hol_Axzc * (-2);
    M_coeff_(MHU1,MKQ13) = hol_Axy  * (-2);
    M_coeff_(MHU1,MKQ23) = hol_Axx  * 2;
    M_coeff_(MHU1,MPQS)  = hol_Axzc * 2;
    M_coeff_(MHU1,MPQC)  = hol_Axzs * (-2);
    M_coeff_(MHU1,MHU1)  = hol_Axx  * 2;
    M_coeff_(MHU1,MHU2)  = hol_Axy  * 2;
    M_coeff_(MHU1,MPUZC) = hol_Axzc * 2;
    M_coeff_(MHU1,MPUZS) = hol_Axzs * 2;

    // HU2 moment equation
    M_coeff_(MHU2,MKI11) = hol_Ayzs;
    M_coeff_(MHU2,MKI22) = hol_Ayzs * (-1);
    M_coeff_(MHU2,MKI12) = hol_Ayzc * (-2);
    M_coeff_(MHU2,MKI13) = hol_Ayy  * (-2);
    M_coeff_(MHU2,MKI23) = hol_Axy  * 2;
    M_coeff_(MHU2,MKQ11) = hol_Ayzs;
    M_coeff_(MHU2,MKQ22) = hol_Ayzs * (-1);
    M_coeff_(MHU2,MKQ12) = hol_Ayzc * (-2);
    M_coeff_(MHU2,MKQ13) = hol_Ayy  * (-2);
    M_coeff_(MHU2,MKQ23) = hol_Axy  * 2;
    M_coeff_(MHU2,MPQS)  = hol_Ayzc * 2;
    M_coeff_(MHU2,MPQC)  = hol_Ayzs * (-2);
    M_coeff_(MHU2,MHU1)  = hol_Axy  * 2;
    M_coeff_(MHU2,MHU2)  = hol_Ayy  * 2;
    M_coeff_(MHU2,MPUZC) = hol_Ayzc * 2;
    M_coeff_(MHU2,MPUZS) = hol_Ayzs * 2;

    // PUZC moment equation
    M_coeff_(MPUZC,MKI11) = hol_Azzcs;
    M_coeff_(MPUZC,MKI22) = hol_Azzcs * (-1);
    M_coeff_(MPUZC,MKI12) = hol_Azzcc * (-2);
    M_coeff_(MPUZC,MKI13) = hol_Ayzc  * (-2);
    M_coeff_(MPUZC,MKI23) = hol_Axzc  * 2;
    M_coeff_(MPUZC,MKQ11) = hol_Azzcs;
    M_coeff_(MPUZC,MKQ22) = hol_Azzcs * (-1);
    M_coeff_(MPUZC,MKQ12) = hol_Azzcc * (-2);
    M_coeff_(MPUZC,MKQ13) = hol_Ayzc  * (-2);
    M_coeff_(MPUZC,MKQ23) = hol_Axzc  * 2;
    M_coeff_(MPUZC,MPQS)  = hol_Azzcc * 2;
    M_coeff_(MPUZC,MPQC)  = hol_Azzcs * (-2);
    M_coeff_(MPUZC,MHU1)  = hol_Axzc  * 2;
    M_coeff_(MPUZC,MHU2)  = hol_Ayzc  * 2;
    M_coeff_(MPUZC,MPUZC) = hol_Azzcc * 2;
    M_coeff_(MPUZC,MPUZS) = hol_Azzcs * 2;

    // PUZS moment equation
    M_coeff_(MPUZS,MKI11) = hol_Azzss;
    M_coeff_(MPUZS,MKI22) = hol_Azzss * (-1);
    M_coeff_(MPUZS,MKI12) = hol_Azzcs * (-2);
    M_coeff_(MPUZS,MKI13) = hol_Ayzs  * (-2);
    M_coeff_(MPUZS,MKI23) = hol_Axzs  * 2;
    M_coeff_(MPUZS,MKQ11) = hol_Azzss;
    M_coeff_(MPUZS,MKQ22) = hol_Azzss * (-1);
    M_coeff_(MPUZS,MKQ12) = hol_Azzcs * (-2);
    M_coeff_(MPUZS,MKQ13) = hol_Ayzs  * (-2);
    M_coeff_(MPUZS,MKQ23) = hol_Axzs  * 2;
    M_coeff_(MPUZS,MPQS)  = hol_Azzcs * 2;
    M_coeff_(MPUZS,MPQC)  = hol_Azzss * (-2);
    M_coeff_(MPUZS,MHU1)  = hol_Axzs  * 2;
    M_coeff_(MPUZS,MHU2)  = hol_Ayzs  * 2;
    M_coeff_(MPUZS,MPUZC) = hol_Azzcs * 2;
    M_coeff_(MPUZS,MPUZS) = hol_Azzss * 2;

    // HV1 moment equation
    M_coeff_(MHV1,MHV1) = hol_Axx * 2;
    M_coeff_(MHV1,MHV2) = hol_Axy * 2;
    M_coeff_(MHV1,MHV3) = hol_Axz * 2;

    // HV2 moment equation
    M_coeff_(MHV2,MHV1) = hol_Axy * 2;
    M_coeff_(MHV2,MHV2) = hol_Ayy * 2;
    M_coeff_(MHV2,MHV3) = hol_Ayz * 2;

    // HV3 moment equation
    M_coeff_(MHV3,MHV1) = hol_Axz * 2;
    M_coeff_(MHV3,MHV2) = hol_Ayz * 2;
    M_coeff_(MHV3,MHV3) = hol_Azz * 2;

    // complete the calculation of M_coeff_
    for (int m=0; m<23; ++m) {
      for (int n=0; n<23; ++n) {
        M_coeff_(m,n) *= -0.75*chi_s;
        if (m == n) M_coeff_(m,n) += 1.0;
      } // endfor n
    } // endfor m
  } // endfor ifr

} // end RadIntegrator::CalPolAux


//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::PolAbsScat()
//  \brief

// wmu_cm is the weight in the co-moving frame
// wmu_cm=wmu * 1/(1-vdotn/Crat)^2 / Lorz^2
// tran_coef is (1-vdotn/Crat)*Lorz, or hollow-L in LZ's notes
// rho is gas density
// tgas is gas temperature
// This function updates normal absorption plus scattering opacity together

Real RadIntegrator::PolAbsScat(
    AthenaArray<Real> &wmu_cm, AthenaArray<Real> &tran_coef,
    AthenaArray<Real> &nx_cm, AthenaArray<Real> &ny_cm, AthenaArray<Real> &nz_cm,
    Real *sigma_a, Real *sigma_p, Real *sigma_pe, Real *sigma_s,
    Real dt, Real lorz, Real rho, Real &tgas, AthenaArray<Real> &ir_cm) {

  /*************** Step 0: Prepare Auxiliary for Computation ***************/
  const Real& prat  = pmy_rad->prat;
  const int&  nang  = pmy_rad->nang;
  const int&  nfreq = pmy_rad->nfreq;

  // velocity reduction
  Real cdt = dt * pmy_rad->crat;
  cdt *= pmy_rad->reduced_c/pmy_rad->crat;

  // gas adiabatic index
  Real gamma = pmy_rad->pmy_block->peos->GetGamma();
  Real gm1 = gamma - 1;

  // angles and coefficents
  Real *wmu = &(wmu_cm(0));
  Real *nx  = &(nx_cm(0));
  Real *ny  = &(ny_cm(0));
  Real *nz  = &(nz_cm(0));
  Real *coef_l = &(tran_coef(0));

  // start iteration on frequeny (only have one for current polarized radiation)
  bool badcell=false;
  Real tgas_new = tgas;
  Real coef[2];
  coef[0] = 0.0; coef[1] = 0.0;
  for (int ifr=0; ifr<nfreq; ++ifr) {
    /*************** Step 1: Update Gas Temperature ***************/
    // opacities
    Real chi_r = sigma_a[ifr];  // Rosseland mean absorption
    Real chi_s = sigma_s[ifr];  // Scattering
    Real chi_p = sigma_p[ifr];  // Planck mean absorption
    Real chi_e = sigma_pe[ifr]; // Energy (Planck) mean absorption

    // coefficients for temperature equation
    Real *si_cm  = &(ir_cm(0,nang*ifr));
    Real *sq_cm  = &(ir_cm(1,nang*ifr));
    Real *su_cm  = &(ir_cm(2,nang*ifr));
    Real *sv_cm  = &(ir_cm(3,nang*ifr));
    Real j0_prev=0.0;
    for (int n=0; n<nang; n++) {
      j0_prev += wmu[n] * si_cm[n];
    }

    // compute coefficents
    M_coeff_tmp_ = M_coeff_; // input matrix is going to be modified in LU decomposition process
    InverseMatrix(23, M_coeff_tmp_, M_inv_);
    Real coeff_a=0; Real coeff_b=0;
    for (int m=0; m<23; ++m) {
      coeff_b += M_inv_(0,m) * polVecB_(m);
      coeff_a += M_inv_(0,m) * polVecA_(m);
    }
    coeff_a *= chi_p;

    // solve the temperature equation
    coef[1] = prat * gm1/rho * coeff_a;
    coef[0] = -tgas + prat * gm1/rho * (coeff_b - j0_prev);
    if (std::abs(coef[1]) > TINY_NUMBER) {
      int flag = FouthPolyRoot(coef[1], coef[0], tgas_new);
      if (flag == -1 || (tgas_new != tgas_new)) {
        badcell = true;
        tgas_new = tgas;
      }
    } else {
      tgas_new = -coef[0];
    } // endelse (std::abs(coef[1]) > TINY_NUMBER)

    /*************** Step 2: Update Moments and Stokes Parameters in Fluid Frame ***************/
    // update necessary moments in fluid frame
    for (int m=0; m<23; ++m) {
      pol_mom_(m) = 0.0;
      for (int n=0; n<23; ++n) {
        pol_mom_(m) += M_inv_(m,n)*polVecB_(n);
        pol_mom_(m) += M_inv_(m,n)*polVecA_(n)*chi_p*SQR(SQR(tgas_new));
      } // endfor n
    } // endfor m

    // update stokes parameters in fluid frame
    for (int n=0; n<nang; ++n) {
      Real fac = 1.0 / (1.0 + coef_l[n]*cdt*(chi_r+chi_s));
      Real fac_p = chi_p*fac*coef_l[n]*cdt;
      Real fac_s = 0.75*chi_s*fac*coef_l[n]*cdt;
      Real nc = (SQR(nx[n]) - SQR(ny[n])) / (SQR(nx[n]) + SQR(ny[n]));
      Real ns = 2 * nx[n] * ny[n] / (SQR(nx[n]) + SQR(ny[n]));
      Real terms_in_bracket_i=0;
      Real terms_in_bracket_q=0;
      Real terms_in_bracket_u=0;
      Real terms_in_bracket_v=0;

      // update I
      terms_in_bracket_i += pol_mom_(MJI);
      terms_in_bracket_i += pol_mom_(MKI11) * nx[n] * nx[n];
      terms_in_bracket_i += pol_mom_(MKI22) * ny[n] * ny[n];
      terms_in_bracket_i += pol_mom_(MKI33) * nz[n] * nz[n];
      terms_in_bracket_i += pol_mom_(MKI12) * nx[n] * ny[n] * 2;
      terms_in_bracket_i += pol_mom_(MKI13) * nx[n] * nz[n] * 2;
      terms_in_bracket_i += pol_mom_(MKI23) * ny[n] * nz[n] * 2;
      terms_in_bracket_i += pol_mom_(MJQ)   * nz[n] * nz[n] * (-1);
      terms_in_bracket_i += pol_mom_(MKQ11) * nx[n] * nx[n];
      terms_in_bracket_i += pol_mom_(MKQ22) * ny[n] * ny[n];
      terms_in_bracket_i += pol_mom_(MKQ33) * nz[n] * nz[n];
      terms_in_bracket_i += pol_mom_(MKQ12) * nx[n] * ny[n] * 2;
      terms_in_bracket_i += pol_mom_(MKQ13) * nx[n] * nz[n] * 2;
      terms_in_bracket_i += pol_mom_(MKQ23) * ny[n] * nz[n] * 2;
      terms_in_bracket_i += pol_mom_(MPQC)  * (SQR(ny[n]) - SQR(nx[n]));
      terms_in_bracket_i += pol_mom_(MPQS)  * nx[n] * ny[n] * (-2);
      terms_in_bracket_i += pol_mom_(MHU1)  * ny[n] * nz[n] * 2;
      terms_in_bracket_i += pol_mom_(MHU2)  * nx[n] * nz[n] * (-2);
      terms_in_bracket_i += pol_mom_(MPUZC) * nx[n] * ny[n] * (-2);
      terms_in_bracket_i += pol_mom_(MPUZS) * (SQR(nx[n]) - SQR(ny[n]));
      si_cm[n] *= fac;
      si_cm[n] += fac_p * SQR(SQR(tgas_new));
      si_cm[n] += fac_s * terms_in_bracket_i;

      // update Q
      terms_in_bracket_q += pol_mom_(MKI11) * (SQR(nx[n]) - nc);
      terms_in_bracket_q += pol_mom_(MKI22) * (SQR(ny[n]) + nc);
      terms_in_bracket_q += pol_mom_(MKI33) * (SQR(nz[n]) - 1);
      terms_in_bracket_q += pol_mom_(MKI12) * (nx[n]*ny[n] - ns) * 2;
      terms_in_bracket_q += pol_mom_(MKI13) * nx[n] * nz[n] * 2;
      terms_in_bracket_q += pol_mom_(MKI23) * ny[n] * nz[n] * 2;
      terms_in_bracket_q += pol_mom_(MJQ)   * (1 - SQR(nz[n]));
      terms_in_bracket_q += pol_mom_(MKQ11) * (SQR(nx[n]) - nc);
      terms_in_bracket_q += pol_mom_(MKQ22) * (SQR(ny[n]) + nc);
      terms_in_bracket_q += pol_mom_(MKQ33) * (SQR(nz[n]) - 1);
      terms_in_bracket_q += pol_mom_(MKQ12) * (nx[n]*ny[n] - ns) * 2;
      terms_in_bracket_q += pol_mom_(MKQ13) * nx[n] * nz[n] * 2;
      terms_in_bracket_q += pol_mom_(MKQ23) * ny[n] * nz[n] * 2;
      terms_in_bracket_q += pol_mom_(MPQC)   * (2*nc - SQR(nx[n]) + SQR(ny[n]));
      terms_in_bracket_q += pol_mom_(MPQS)   * (ns - nx[n]*ny[n]) * 2;
      terms_in_bracket_q += pol_mom_(MHU1)  * ny[n] * nz[n] * 2;
      terms_in_bracket_q += pol_mom_(MHU2)  * nx[n] * nz[n] * (-2);
      terms_in_bracket_q += pol_mom_(MPUZC) * (ns - nx[n]*ny[n]) * 2;
      terms_in_bracket_q += pol_mom_(MPUZS) * (-2*nc + SQR(nx[n]) - SQR(ny[n]));
      sq_cm[n] *= fac;
      sq_cm[n] += fac_s * terms_in_bracket_q;

      // update U
      terms_in_bracket_u += pol_mom_(MKI11) * nz[n] * ns;
      terms_in_bracket_u += pol_mom_(MKI22) * nz[n] * ns * (-1);
      terms_in_bracket_u += pol_mom_(MKI12) * nz[n] * nc * (-2);
      terms_in_bracket_u += pol_mom_(MKI13) * ny[n] * (-2);
      terms_in_bracket_u += pol_mom_(MKI23) * nx[n] * 2;
      terms_in_bracket_u += pol_mom_(MKQ11) * nz[n] * ns;
      terms_in_bracket_u += pol_mom_(MKQ22) * nz[n] * ns * (-1);
      terms_in_bracket_u += pol_mom_(MKQ12) * nz[n] * nc * (-2);
      terms_in_bracket_u += pol_mom_(MKQ13) * ny[n] * (-2);
      terms_in_bracket_u += pol_mom_(MKQ23) * nx[n] * 2;
      terms_in_bracket_u += pol_mom_(MPQS)   * nz[n] * nc * 2;
      terms_in_bracket_u += pol_mom_(MPQC)   * nz[n] * ns * (-2);
      terms_in_bracket_u += pol_mom_(MHU1)  * nx[n] * 2;
      terms_in_bracket_u += pol_mom_(MHU2)  * ny[n] * 2;
      terms_in_bracket_u += pol_mom_(MPUZC) * nz[n] * nc * 2;
      terms_in_bracket_u += pol_mom_(MPUZS) * nz[n] * ns * 2;
      su_cm[n] *= fac;
      su_cm[n] += fac_s * terms_in_bracket_u;

      // update V
      terms_in_bracket_v += pol_mom_(MHV1) * nx[n];
      terms_in_bracket_v += pol_mom_(MHV2) * ny[n];
      terms_in_bracket_v += pol_mom_(MHV3) * nz[n];
      sv_cm[n] *= fac;
      su_cm[n] += fac_s * 2*terms_in_bracket_v;
    } // endfor n

  } // endfor ifr

  return tgas_new;
}


// use Newton-Raphson scheme
// Real RadIntegrator::PolAbsScat(
//     AthenaArray<Real> &wmu_cm, AthenaArray<Real> &tran_coef,
//     AthenaArray<Real> &nx_cm, AthenaArray<Real> &ny_cm, AthenaArray<Real> &nz_cm,
//     Real *sigma_a, Real *sigma_p, Real *sigma_pe, Real *sigma_s,
//     Real dt, Real lorz, Real rho, Real &tgas, AthenaArray<Real> &ir_cm) {
//
//   /*************** Step 0: Prepare Auxiliary for Computation ***************/
//   Real tol = 1e-12;
//   Real num_max_itr = 50;
//   bool& tst_tgas_ini_guess_ = pmy_rad->tst_tgas_ini_guess;
//   const Real& prat  = pmy_rad->prat;
//   const int&  nang  = pmy_rad->nang;
//   const int&  nfreq = pmy_rad->nfreq;
//
//   // velocity reduction
//   Real cdt = dt * pmy_rad->crat;
//   cdt *= pmy_rad->reduced_c/pmy_rad->crat;
//
//   // gas adiabatic index
//   Real gamma = pmy_rad->pmy_block->peos->GetGamma();
//   Real gm1 = gamma - 1;
//
//   // angles and coefficents
//   Real *wmu = &(wmu_cm(0));
//   Real *nx  = &(nx_cm(0));
//   Real *ny  = &(ny_cm(0));
//   Real *nz  = &(nz_cm(0));
//   Real *coef_l = &(tran_coef(0));
//
//   // start iteration on frequeny (only have one for current polarized radiation)
//   bool badcell=false;
//   Real tgas_guess = tgas; Real tgas_new = tgas;
//   Real coef[2];
//   coef[0] = 0.0; coef[1] = 0.0;
//   for (int ifr=0; ifr<nfreq; ++ifr) {
//     /*************** Step 1: Initial Guess of Gas Temperature ***************/
//     // opacities
//     Real chi_r = sigma_a[ifr];  // Rosseland mean absorption
//     Real chi_s = sigma_s[ifr];  // Scattering
//     Real chi_p = sigma_p[ifr];  // Planck mean absorption
//     Real chi_e = sigma_pe[ifr]; // Energy (Planck) mean absorption
//
//     // coefficients for temperature equation
//     Real *si_cm  = &(ir_cm(0,nang*ifr));
//     Real *sq_cm  = &(ir_cm(1,nang*ifr));
//     Real *su_cm  = &(ir_cm(2,nang*ifr));
//     Real *sv_cm  = &(ir_cm(3,nang*ifr));
//     Real c1_sum=0.0, c2_sum=0.0, j0_prev=0.0;
//     for (int n=0; n<nang; n++) {
//       Real fac = 1.0 / (1.0 + coef_l[n]*cdt*(chi_r+chi_s));
//       c1_sum  += fac * coef_l[n] * wmu[n] * cdt;
//       c2_sum  += fac * wmu[n] * si_cm[n];
//       j0_prev += wmu[n] * si_cm[n];
//     }
//     Real c3_sum = 1.0 / (1.0 - c1_sum*(chi_r+chi_s-chi_p));
//
//     // solve the temperature equation assuming non-polarized radiation
//     coef[1] = prat * gm1/rho * c3_sum*c1_sum*chi_p;
//     coef[0] = -tgas + prat * gm1/rho * (c3_sum*c2_sum - j0_prev);
//     if (std::abs(coef[1]) > TINY_NUMBER) {
//       int flag = FouthPolyRoot(coef[1], coef[0], tgas_guess);
//       if (flag == -1 || (tgas_guess != tgas_guess)) {
//         badcell = true;
//         tgas_guess = tgas;
//       }
//     } else {
//       tgas_guess = -coef[0];
//     } // endelse (std::abs(coef[1]) > TINY_NUMBER)
//
//     // if only test temperature initial guess, reset co-moving I
//     if (tst_tgas_ini_guess_ && !badcell) {
//       Real tgas4 = SQR(SQR(tgas_guess));
//       Real j0_new = c3_sum * (c1_sum*chi_p*tgas4 + c2_sum);
//       // update the co-moving frame specific intensity
//       for (int n=0; n<nang; n++) {
//         Real fac1 = 1.0 / (1.0 + coef_l[n]*cdt*(chi_r+chi_s));
//         Real fac2 = fac1*coef_l[n]*cdt;
//         si_cm[n] = fac1*si_cm[n];
//         si_cm[n] += fac2 * (chi_p*tgas4 + (chi_r+chi_s-chi_p)*j0_new);
//         si_cm[n] = std::max(si_cm[n],static_cast<Real>(TINY_NUMBER));
//       } // endfor n
//       tgas_new = tgas_guess;
//     } else {
//       /*************** Step 2: Newton-Raphson Iteration to Update Gas Temperature ***************/
//       // compute coefficents
//      InverseMatrix(23, M_coeff_, M_inv_);
//       Real coeff_a=0; Real coeff_b=0;
//       for (int m=0; m<23; ++m) {
//         coeff_b += M_inv_(0,m) * polVecB_(m);
//         coeff_a += M_inv_(0,m) * polVecA_(m);
//       }
//       coeff_a *= chi_p;
//
//       // update gas temperature
//       Real func, dfunc, tgas_1;
//       Real tgas_0 = tgas_guess;
//       Real j0_update = coeff_b + coeff_a*SQR(SQR(tgas_0));
//       int count; Real l1_err = 1.0;
//       for (count=0; count<num_max_itr; ++count) {
//         func = rho/gm1*(tgas_0-tgas) + prat*(j0_update-j0_prev);
//         dfunc = rho/gm1 + 4*prat*coeff_a*(tgas_0*tgas_0*tgas_0);
//         tgas_1 = tgas_0 - func/dfunc;
//         l1_err = fabs(tgas_1-tgas_0);
//         tgas_0 = tgas_1;
//         j0_update = coeff_b + coeff_a*SQR(SQR(tgas_0));
//         if (tgas_0 < TINY_NUMBER) { // unphysical value
//           badcell = true;
//           break;
//         }
//         if (l1_err < tol) // solution converged
//           break;
//       } // endfor m
//       if (count==num_max_itr) badcell=true; // not converging
//       if (badcell) tgas_0 = tgas_guess; // fix the temperature using the initial guess
//
//       // updated temperature
//       tgas_new = tgas_0;
//
//       /*************** Step 3: Update Moments and Stokes Parameters in Fluid Frame ***************/
//       // update necessary moments in fluid frame
//       for (int m=0; m<23; ++m) {
//         pol_mom_(m) = 0.0;
//         for (int n=0; n<23; ++n) {
//           pol_mom_(m) += M_inv_(m,n)*polVecB_(n);
//           pol_mom_(m) += M_inv_(m,n)*polVecA_(n)*chi_p*SQR(SQR(tgas_new));
//         } // endfor n
//       } // endfor m
//
//       // update stokes parameters in fluid frame
//       for (int n=0; n<nang; ++n) {
//         Real fac = 1.0 / (1.0 + coef_l[n]*cdt*(chi_r+chi_s));
//         Real fac_p = chi_p*fac*coef_l[n]*cdt;
//         Real fac_s = 0.75*chi_s*fac*coef_l[n]*cdt;
//         Real nc = (SQR(nx[n]) - SQR(ny[n])) / (SQR(nx[n]) + SQR(ny[n]));
//         Real ns = 2 * nx[n] * ny[n] / (SQR(nx[n]) + SQR(ny[n]));
//         Real terms_in_bracket_i=0;
//         Real terms_in_bracket_q=0;
//         Real terms_in_bracket_u=0;
//         Real terms_in_bracket_v=0;
//
//         // update I
//         terms_in_bracket_i += pol_mom_(MJI);
//         terms_in_bracket_i += pol_mom_(MKI11) * nx[n] * nx[n];
//         terms_in_bracket_i += pol_mom_(MKI22) * ny[n] * ny[n];
//         terms_in_bracket_i += pol_mom_(MKI33) * nz[n] * nz[n];
//         terms_in_bracket_i += pol_mom_(MKI12) * nx[n] * ny[n] * 2;
//         terms_in_bracket_i += pol_mom_(MKI13) * nx[n] * nz[n] * 2;
//         terms_in_bracket_i += pol_mom_(MKI23) * ny[n] * nz[n] * 2;
//         terms_in_bracket_i += pol_mom_(MJQ)   * nz[n] * nz[n] * (-1);
//         terms_in_bracket_i += pol_mom_(MKQ11) * nx[n] * nx[n];
//         terms_in_bracket_i += pol_mom_(MKQ22) * ny[n] * ny[n];
//         terms_in_bracket_i += pol_mom_(MKQ33) * nz[n] * nz[n];
//         terms_in_bracket_i += pol_mom_(MKQ12) * nx[n] * ny[n] * 2;
//         terms_in_bracket_i += pol_mom_(MKQ13) * nx[n] * nz[n] * 2;
//         terms_in_bracket_i += pol_mom_(MKQ23) * ny[n] * nz[n] * 2;
//         terms_in_bracket_i += pol_mom_(MPQC)  * (SQR(ny[n]) - SQR(nx[n]));
//         terms_in_bracket_i += pol_mom_(MPQS)  * nx[n] * ny[n] * (-2);
//         terms_in_bracket_i += pol_mom_(MHU1)  * ny[n] * nz[n] * 2;
//         terms_in_bracket_i += pol_mom_(MHU2)  * nx[n] * nz[n] * (-2);
//         terms_in_bracket_i += pol_mom_(MPUZC) * nx[n] * ny[n] * (-2);
//         terms_in_bracket_i += pol_mom_(MPUZS) * (SQR(nx[n]) - SQR(ny[n]));
//         si_cm[n] *= fac;
//         si_cm[n] += fac_p * SQR(SQR(tgas_new));
//         si_cm[n] += fac_s * terms_in_bracket_i;
//
//         // update Q
//         terms_in_bracket_q += pol_mom_(MKI11) * (SQR(nx[n]) - nc);
//         terms_in_bracket_q += pol_mom_(MKI22) * (SQR(ny[n]) + nc);
//         terms_in_bracket_q += pol_mom_(MKI33) * (SQR(nz[n]) - 1);
//         terms_in_bracket_q += pol_mom_(MKI12) * (nx[n]*ny[n] - ns) * 2;
//         terms_in_bracket_q += pol_mom_(MKI13) * nx[n] * nz[n] * 2;
//         terms_in_bracket_q += pol_mom_(MKI23) * ny[n] * nz[n] * 2;
//         terms_in_bracket_q += pol_mom_(MJQ)   * (1 - SQR(nz[n]));
//         terms_in_bracket_q += pol_mom_(MKQ11) * (SQR(nx[n]) - nc);
//         terms_in_bracket_q += pol_mom_(MKQ22) * (SQR(ny[n]) + nc);
//         terms_in_bracket_q += pol_mom_(MKQ33) * (SQR(nz[n]) - 1);
//         terms_in_bracket_q += pol_mom_(MKQ12) * (nx[n]*ny[n] - ns) * 2;
//         terms_in_bracket_q += pol_mom_(MKQ13) * nx[n] * nz[n] * 2;
//         terms_in_bracket_q += pol_mom_(MKQ23) * ny[n] * nz[n] * 2;
//         terms_in_bracket_q += pol_mom_(MPQC)   * (2*nc - SQR(nx[n]) + SQR(ny[n]));
//         terms_in_bracket_q += pol_mom_(MPQS)   * (ns - nx[n]*ny[n]) * 2;
//         terms_in_bracket_q += pol_mom_(MHU1)  * ny[n] * nz[n] * 2;
//         terms_in_bracket_q += pol_mom_(MHU2)  * nx[n] * nz[n] * (-2);
//         terms_in_bracket_q += pol_mom_(MPUZC) * (ns - nx[n]*ny[n]) * 2;
//         terms_in_bracket_q += pol_mom_(MPUZS) * (-2*nc + SQR(nx[n]) - SQR(ny[n]));
//         sq_cm[n] *= fac;
//         sq_cm[n] += fac_s * terms_in_bracket_q;
//
//         // update U
//         terms_in_bracket_u += pol_mom_(MKI11) * nz[n] * ns;
//         terms_in_bracket_u += pol_mom_(MKI22) * nz[n] * ns * (-1);
//         terms_in_bracket_u += pol_mom_(MKI12) * nz[n] * nc * (-2);
//         terms_in_bracket_u += pol_mom_(MKI13) * ny[n] * (-2);
//         terms_in_bracket_u += pol_mom_(MKI23) * nx[n] * 2;
//         terms_in_bracket_u += pol_mom_(MKQ11) * nz[n] * ns;
//         terms_in_bracket_u += pol_mom_(MKQ22) * nz[n] * ns * (-1);
//         terms_in_bracket_u += pol_mom_(MKQ12) * nz[n] * nc * (-2);
//         terms_in_bracket_u += pol_mom_(MKQ13) * ny[n] * (-2);
//         terms_in_bracket_u += pol_mom_(MKQ23) * nx[n] * 2;
//         terms_in_bracket_u += pol_mom_(MPQS)   * nz[n] * nc * 2;
//         terms_in_bracket_u += pol_mom_(MPQC)   * nz[n] * ns * (-2);
//         terms_in_bracket_u += pol_mom_(MHU1)  * nx[n] * 2;
//         terms_in_bracket_u += pol_mom_(MHU2)  * ny[n] * 2;
//         terms_in_bracket_u += pol_mom_(MPUZC) * nz[n] * nc * 2;
//         terms_in_bracket_u += pol_mom_(MPUZS) * nz[n] * ns * 2;
//         su_cm[n] *= fac;
//         su_cm[n] += fac_s * terms_in_bracket_u;
//
//         // update V
//         terms_in_bracket_v += pol_mom_(MHV1) * nx[n];
//         terms_in_bracket_v += pol_mom_(MHV2) * ny[n];
//         terms_in_bracket_v += pol_mom_(MHV3) * nz[n];
//         sv_cm[n] *= fac;
//         su_cm[n] += fac_s * 2*terms_in_bracket_v;
//       } // endfor n
//
//     } // endelse (tst_tgas_ini_guess && !badcell)
//   } // endfor ifr
//
//   return tgas_new;
// }
