//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rec_shear.cpp
//! \brief Problem generator for double-layer shearing sheets.
//! HGB: short for Hawley, Gammie & Balbus (see "Local Three-dimensional Magnetohydrodynamic
//!  Simulations of Accretion Disks")
//============================================================================

// C++ headers
#include <cmath>      // sqrt()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp" // ran2()

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

// Explicit source/sink terms
void SrcTerms(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

// Opacity function
void FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

// Hydro and radiation boundary conditions
// void HydroInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//                   Real time, Real dt,
//                   int il, int iu, int jl, int ju, int kl, int ku, int ngh);
//
// void HydroOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//                   Real time, Real dt,
//                   int il, int iu, int jl, int ju, int kl, int ku, int ngh);
//
// void RadInnerX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
//                 const AthenaArray<Real> &w, FaceField &b,
//                 AthenaArray<Real> &ir, Real time, Real dt,
//                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
//
// void RadOuterX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
//                 const AthenaArray<Real> &w, FaceField &b,
//                 AthenaArray<Real> &ir, Real time, Real dt,
//                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
  Real z_mid, z_max, z_min;                // coordinate paramters
  Real z_lower, z_upper, jwidth;           // current sheet paramters
  Real iso_cs, gm1, d0, p0;                // fluid parameters
  Real beta, beta_mri, bx2by, bz2bx;       // magnetic parameters
  Real dfloor, pfloor;                     // floors
  Real Lx, Ly, Lz;                         // root grid size, global to share with output functions
  Real Omega_0, qshear;                    // shearing box parameters
  Real nwx, nwy;                           // wavenumbers
  Real amp;                                // random perturbation amplitude
  Real grav, sinkwidth, tau_sink;          // source and sink parameters
  Real kappa0_sct, kappa0_ffr, kappa0_ffp; // opacity paramters

  Real HistoryBxBy(MeshBlock *pmb, int iout);
  Real HistorydVxVy(MeshBlock *pmb, int iout);
} // end namespace

// ===================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {

  Real float_min = std::numeric_limits<float>::min();
  z_min = pin->GetReal("mesh", "x3min");
  z_max = pin->GetReal("mesh", "x3max");
  z_mid = 0.5*(z_max+z_min);
  z_lower = (z_mid + z_min)/2.0;
  z_upper = (z_mid + z_max)/2.0;
  gm1 = pin->GetOrAddReal("hydro", "gamma", 5./3) - 1;
  grav = pin->GetOrAddReal("problem", "grav", 0.0);
  sinkwidth = pin->GetOrAddReal("problem", "sinkwidth", 0.0);
  tau_sink = pin->GetOrAddReal("problem", "tau_sink", 1.0);
  dfloor = pin->GetOrAddReal("hydro", "dfloor", 1024*(float_min));
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 1024*(float_min));
  kappa0_sct = pin->GetOrAddReal("radiation", "electron_scattering", 0.34);
  kappa0_ffr = pin->GetOrAddReal("radiation", "kappa_ffr", 1024*(float_min));
  kappa0_ffp = pin->GetOrAddReal("radiation", "kappa_ffp", 1024*(float_min));

  // Overwrite periodic boundaries in z-direction
  // if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
  //   EnrollUserBoundaryFunction(BoundaryFace::inner_x3, HydroInnerX3);
  //   EnrollUserBoundaryFunction(BoundaryFace::outer_x3, HydroOuterX3);
  //   EnrollUserRadBoundaryFunction(BoundaryFace::inner_x3, RadInnerX3);
  //   EnrollUserRadBoundaryFunction(BoundaryFace::outer_x3, RadOuterX3);
  // }

  // User-defined outputs
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, HistoryBxBy, "-BxBy");
  EnrollUserHistoryOutput(1, HistorydVxVy, "dVxVy");

  // User-defined physical source/sink terms
  EnrollUserExplicitSourceFunction(SrcTerms);

  // Module check
  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in rec_shear.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in rec_shear.cpp ProblemGenerator" << std::endl
        << "This problem generator works only in 2D or 3D." << std::endl;
    ATHENA_ERROR(msg);
  }

  if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) && (!NON_BAROTROPIC_EOS)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in rec_shear.cpp ProblemGenerator" << std::endl
        << "EoS cannot be barotropic when radiation is included in this problem generator." << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
} // end Mesh::InitUserMeshData

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)
    pnrrad->EnrollOpacityFunction(FreeFreeOpacity);

  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief mhd shearing waves and unstratified disk problem generator for
//  3D problems.
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real SumRvx=0.0, SumRvy=0.0, SumRvz=0.0;

  // gamma, press, sound speed
  Real gamma = 1.0;
  d0     = pin->GetOrAddReal("problem", "d0",     1.0);
  jwidth = pin->GetOrAddReal("problem", "jwidth", 0.1);
  bx2by  = pin->GetOrAddReal("problem", "bx2by",  0.02);
  bz2bx  = pin->GetOrAddReal("problem", "bz2bx",  0.01);
  if (NON_BAROTROPIC_EOS) {
    p0     = pin->GetReal("problem","p0");
    gamma  = peos->GetGamma();
    iso_cs = std::sqrt(gamma*p0/d0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0 = d0*SQR(iso_cs);
  }

  // shearing box parameter
  if (porb->shboxcoord != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in rec_shear.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box in x-y plane." << std::endl
        << "Check <orbital_advection> shboxcoord parameter." << std::endl;
    ATHENA_ERROR(msg);
  }
  Omega_0 = porb->Omega0;
  qshear  = porb->qshear;

  // Read problem parameters for initial conditions
  amp = pin->GetReal("problem","amp");
  beta = pin->GetReal("problem", "beta");
  beta_mri = pin->GetReal("problem", "beta_mri");

  // Compute field strength based on beta.
  Real By0 = std::sqrt(static_cast<Real>(2.0*p0/beta));
  // Real Bz0 = std::sqrt(static_cast<Real>(2.0*p0/beta_mri));
  Real Bz0 = 0.0;

  // Ensure a different initial random seed for each meshblock.
  std::int64_t iseed = - 1 - gid;

  // Initialize boxsize
  Lx = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  Ly = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  Lz = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

  // initialize wavenumbers
  nwx = pin->GetOrAddInteger("problem","nwx",1);
  nwy = pin->GetOrAddInteger("problem","nwy",1);

  Real kx = (TWO_PI/Lx)*(static_cast<Real>(nwx)); // nxw = -ve for leading wave
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));

  // Initialize perturbations
  // hydro
  Real x1, x3;
  Real rd(0.0), rp(0.0), rvx(0.0), rvy(0.0), rvz(0.0);
  Real rval;
  for (int k=ks; k<=ke; k++) {
    x3 = pcoord->x3v(k);
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        x1  = pcoord->x1v(i);
        rval = amp*(ran2(&iseed) - 0.5);
        if (NON_BAROTROPIC_EOS) {
          rp = p0*(1.0 + 2.0*rval);
          rd = d0;
        } else {
          rd  = d0;
          rd += SQR(By0/iso_cs)*(1.0+SQR(bx2by))/2.0 * (1.0 - SQR((std::tanh((x3-z_upper)/jwidth)-std::tanh((x3-z_lower)/jwidth)+1.0)));
        }

        // Following HGB: the perturbations to V/Cs are
        // (1/5)amp/std::sqrt(gamma)
        rval = amp*(ran2(&iseed) - 0.5);
        rvx = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
        SumRvx += rvx;

        rval = amp*(ran2(&iseed) - 0.5);
        rvy = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
        SumRvy += rvy;

        rval = amp*(ran2(&iseed) - 0.5);
        rvz = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
        SumRvz += rvz;

        // Initialize (d, M, P)
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = rd*rvx;
        phydro->u(IM2,k,j,i) = rd*rvy;
        if(!porb->orbital_advection_defined)
          phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
        phydro->u(IM3,k,j,i) = rd*rvz;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = rp/(gamma-1.0)
                                 + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                       +SQR(phydro->u(IM2,k,j,i))
                                       +SQR(phydro->u(IM3,k,j,i)))/rd;
        }
      }  // endfor i
    }  // endfor j
  }  // endfor k

  // For random perturbations as in HGB, ensure net momentum is zero by
  // subtracting off mean of perturbations
  int cell_num = block_size.nx1*block_size.nx2*block_size.nx3;
  SumRvx /= cell_num;
  SumRvy /= cell_num;
  SumRvz /= cell_num;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IM1,k,j,i) -= rd*SumRvx;
        phydro->u(IM2,k,j,i) -= rd*SumRvy;
        phydro->u(IM3,k,j,i) -= rd*SumRvz;
      } // endfor i
    } // endfor j
  } // endfor k

  // Initialize magnetic fields for a double-layer Harris sheet
  for (int k=ks; k<=ke; k++) {
    Real x3 = pcoord->x3v(k);
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x1  = pcoord->x1v(i);
        pfield->b.x1f(k,j,i) = (By0*bx2by) * (std::tanh((x3-z_upper)/jwidth) - std::tanh((x3-z_lower)/jwidth) + 1.0);
        pfield->b.x2f(k,j,i) =  By0        * (std::tanh((x3-z_upper)/jwidth) - std::tanh((x3-z_lower)/jwidth) + 1.0);
        pfield->b.x3f(k,j,i) =  By0*bx2by*bz2bx;
        pfield->b.x3f(k,j,i) += Bz0*(std::sin(static_cast<Real>(kx)*x1)) * (std::exp(-SQR((x3-z_upper)/jwidth)/2.0) + std::exp(-SQR((x3-z_lower)/jwidth)/2.0));
        // let gas pressure balance the lack of B pressure in current sheets
        if (i==ie) pfield->b.x1f(k,j,ie+1) = (By0*bx2by) * (std::tanh((x3-z_upper)/jwidth) - std::tanh((x3-z_lower)/jwidth) + 1.0);
        if (j==je) pfield->b.x2f(k,je+1,i) =  By0        * (std::tanh((x3-z_upper)/jwidth) - std::tanh((x3-z_lower)/jwidth) + 1.0);
        if (k==ke) {
          pfield->b.x3f(ke+1,j,i) =  By0*bx2by*bz2bx;
          pfield->b.x3f(ke+1,j,i) += Bz0*(std::sin(static_cast<Real>(kx)*x1)) * (std::exp(-SQR((x3-z_upper)/jwidth)/2.0) + std::exp(-SQR((x3-z_lower)/jwidth)/2.0));
        }
        //if (k==ke) pfield->b.x3f(ke+1,j,i) = By0*bx2by*bz2bx;
      } // endfor i
    } // endfor j
  } // endfor k

  // Add magnetic energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) +=
              0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
                   SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
                   SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i))));
        } // endfor i
      } // endfor j
    } // endfor k
  } // endif (NON_BAROTROPIC_EOS)

  // Initialize radiation
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for(int ifr=0; ifr<pnrrad->nfreq; ++ifr){
            // initialize opacity
            pnrrad->sigma_s(k,j,i,ifr)  = 0.0;
            pnrrad->sigma_a(k,j,i,ifr)  = 0.0;
            pnrrad->sigma_pe(k,j,i,ifr) = 0.0;
            pnrrad->sigma_p(k,j,i,ifr)  = 0.0;
            // initialize radiation
            for(int n=0; n<pnrrad->nang; ++n){
              pnrrad->ir(k,j,i,ifr*pnrrad->nang+n) = 0.0;
            } // endfor n
          } // endfor ifr
        } // endfor i
      } // endfor j
    } // endfor k
  } // endif radiation_enable

  return;
}


// User-defined physical source/sink terms
void SrcTerms(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar) {

  const Real Omega_0 = pmb->porb->Omega0;
  NRRadiation *prad = (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) ? pmb->pnrrad : NULL;
  Real zsink_l = z_mid - sinkwidth/2;
  Real zsink_u = z_mid + sinkwidth/2;
  Real zbot_l  = z_min;
  Real zbot_u  = z_min + sinkwidth/2;
  Real ztop_l  = z_max - sinkwidth/2;
  Real ztop_u  = z_max;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real x3 = pmb->pcoord->x3v(k);
    Real z0 = (x3 > z_mid) ? z_upper : z_lower;

    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // gravitational source terms
        Real den = prim(IDN,k,j,i);
        cons(IM3,k,j,i) -= grav*dt*den*SQR(Omega_0)*(x3-z0);
        if (NON_BAROTROPIC_EOS) {
          cons(IEN,k,j,i) -= grav*dt*den*SQR(Omega_0)*prim(IVZ,k,j,i)*(x3-z0);
        }

        // gas and radiation sink at the central region
        bool in_mid_sink = (sinkwidth > 0.0)
                           && (x3 >= zsink_l) && (x3 <= zsink_u);
        bool in_bot_sink = (sinkwidth > 0.0)
                           && (x3 >= zbot_l)  && (x3 <= zbot_u);
        bool in_top_sink = (sinkwidth > 0.0)
                           && (x3 >= ztop_l)  && (x3 <= ztop_u);

        if ((NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)
            && (in_mid_sink || in_bot_sink || in_top_sink)) {

          // exponential decay factor that is stable for any dt/tau ratio
          Real window = 0; // smooth window
          if (in_mid_sink) {
            window = SQR(std::cos(M_PI*(x3 - z_mid)/sinkwidth));
          } else if (in_bot_sink) {
            window = SQR(std::cos(M_PI/2.0*(x3 - zbot_l)/sinkwidth));
          } else if (in_top_sink) {
            window = SQR(std::cos(M_PI/2.0*(x3 - ztop_u)/sinkwidth));
          }
          Real fac = std::exp(-window * dt / tau_sink);

          // extract velocity and magnetic energy
          Real vx = cons(IM1,k,j,i)/cons(IDN,k,j,i);
          Real vy = cons(IM2,k,j,i)/cons(IDN,k,j,i);
          Real vz = cons(IM3,k,j,i)/cons(IDN,k,j,i);
          Real emag = (!NON_BAROTROPIC_EOS) ? 0
                      : 0.5*(SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)));

          // apply sink
          Real rho_update = std::max(cons(IDN,k,j,i)*fac, dfloor);
          cons(IDN,k,j,i) = rho_update;
          cons(IM1,k,j,i) = rho_update*vx;
          cons(IM2,k,j,i) = rho_update*vy;
          cons(IM3,k,j,i) = rho_update*vz;
          if (NON_BAROTROPIC_EOS) {
            Real p_update = std::max(prim(IPR,k,j,i)*fac, pfloor);
            cons(IEN,k,j,i) = emag + p_update/gm1
                            + 0.5*(SQR(cons(IM1,k,j,i))
                                  +SQR(cons(IM2,k,j,i))
                                  +SQR(cons(IM3,k,j,i))
                                  )/cons(IDN,k,j,i);
          }

          // zeros radiation intensity in sink
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
            for(int n=0; n<prad->nang; ++n){
              prad->ir(k,j,i,ifr*prad->nang+n) = 0;
            } // endfor n
          } // endfor ifr
        } // endif radiation_enable

      } // endfor i
    } // endfor j
  } // endfor k

  return;
} // end SrcTerms


// Opacity function
void FreeFreeOpacity(MeshBlock *pmb, AthenaArray<Real> &prim) {
  NRRadiation *prad = pmb->pnrrad;
	int il = pmb->is; int jl = pmb->js; int kl = pmb->ks;
	int iu = pmb->ie; int ju = pmb->je; int ku = pmb->ke;

	il -= NGHOST;
	iu += NGHOST;

	if(ju > jl){
		jl -= NGHOST;
		ju += NGHOST;
	}

	if(ku > kl){
		kl -= NGHOST;
		ku += NGHOST;
	}

  for (int k=kl; k<=ku; ++k) {
		for (int j=jl; j<=ju; ++j) {
			for (int i=il; i<=iu; ++i) {
				for (int ifr=0; ifr<prad->nfreq; ++ifr){
          Real rho  = std::max(prim(IDN,k,j,i), dfloor);
          Real pgas = std::max(prim(IPR,k,j,i), pfloor);
          Real tgas_inv  = rho/pgas;
          Real power_law = rho * tgas_inv*SQR(tgas_inv)*sqrt(tgas_inv);
					prad->sigma_s(k,j,i,ifr)  = rho * kappa0_sct;
					prad->sigma_a(k,j,i,ifr)  = rho * kappa0_ffr * power_law;
          prad->sigma_p(k,j,i,ifr)  = rho * kappa0_ffp * power_law;
          prad->sigma_pe(k,j,i,ifr) = rho * kappa0_ffp * power_law;
				} // endfor ifr
			} // endfor i
		} // endfor j
	} // endfor k

} // end FreeFreeOpacity


// Hydro and radiation boundary conditions
// void HydroInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//                   Real time, Real dt,
//                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
//   // only overwrite density and pressure
//   for (int k=1; k<=ngh; ++k) {
//     for (int j=jl; j<=ju; ++j) {
//       for (int i=il; i<=iu; ++i) {
//         prim(IDN,kl-k,j,i) = dfloor;
//         if (NON_BAROTROPIC_EOS)
//           prim(IPR,kl-k,j,i) = pfloor;
//         // everything else already set by periodic exchange
//       } // endfor i
//     } // endfor j
//   } // endfor k
// } // end HydroInnerX3
//
// void HydroOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
//                   Real time, Real dt,
//                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
//   // only overwrite density and pressure
//   for (int k=1; k<=ngh; ++k) {
//     for (int j=jl; j<=ju; ++j) {
//       for (int i=il; i<=iu; ++i) {
//         prim(IDN,ku+k,j,i) = dfloor;
//         if (NON_BAROTROPIC_EOS)
//           prim(IPR,ku+k,j,i) = pfloor;
//         // everything else already set by periodic exchange
//       } // endfor i
//     } // endfor j
//   } // endfor k
// } // end HydroOuterX3
//
// void RadInnerX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
//                 const AthenaArray<Real> &w, FaceField &b,
//                 AthenaArray<Real> &ir, Real time, Real dt,
//                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
//
//   const int& nang = pmb->pnrrad->nang; // angles per octant
//   const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands
//
//   for (int k=1; k<=ngh; ++k) {
//     for (int j=jl; j<=ju; ++j) {
//       for (int i=il; i<=iu; ++i) {
//         for (int ifr=0; ifr<nfreq; ++ifr) {
//           for(int n=0; n<nang; ++n) {
//             int ang = ifr*nang + n;
//             const Real& miuz=pmb->pnrrad->mu(2,kl,j,i,n);
//             if (miuz < 0.0) {
//               ir(kl-k,j,i,ang) = ir(kl,j,i,ang);
//             } else {
//               ir(kl-k,j,i,ang) = 0.0;
//             }
//           } // endfor n
//         } // endfor ifr
//       } // endfor i
//     } // endfor j
//   } // endfor k
//
// } // end RadInnerX3
//
// void RadOuterX3(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
//                 const AthenaArray<Real> &w, FaceField &b,
//                 AthenaArray<Real> &ir, Real time, Real dt,
//                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
//
//   const int& nang = pmb->pnrrad->nang; // angles per octant
//   const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands
//
//   for (int k=1; k<=ngh; ++k) {
//     for (int j=jl; j<=ju; ++j) {
//       for (int i=il; i<=iu; ++i) {
//         for (int ifr=0; ifr<nfreq; ++ifr) {
//           for(int n=0; n<nang; ++n) {
//             int ang = ifr*nang + n;
//             const Real& miuz=pmb->pnrrad->mu(2,ku,j,i,n);
//             if (miuz > 0.0) {
//               ir(ku+k,j,i,ang) = ir(ku,j,i,ang);
//             } else {
//               ir(ku+k,j,i,ang) = 0.0;
//             }
//           } // endfor n
//         } // endfor ifr
//       } // endfor i
//     } // endfor j
//   } // endfor k
//
// } // end RadOuterX3


// User-defined outputs
namespace {
Real HistoryBxBy(MeshBlock *pmb, int iout) {
  Real bxby = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        bxby-=volume(i)*b(IB1,k,j,i)*b(IB2,k,j,i);
      }
    }
  }
  return bxby;
} // end HistoryBxBy

Real HistorydVxVy(MeshBlock *pmb, int iout) {
  Real dvxvy = 0.0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &w = pmb->phydro->w;
  Real vshear = 0.0;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  const Real qshear = pmb->porb->qshear;
  const Real Omega_0 = pmb->porb->Omega0;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        if(!pmb->porb->orbital_advection_defined) {
          vshear = -qshear*Omega_0*pmb->pcoord->x1v(i);
        } else {
          vshear = 0.0;
        }
        dvxvy += volume(i)*w(IDN,k,j,i)*w(IVX,k,j,i)*(w(IVY,k,j,i)+vshear);
      }
    }
  }
  return dvxvy;
} // end HistorydVxVy
} // end namespace
