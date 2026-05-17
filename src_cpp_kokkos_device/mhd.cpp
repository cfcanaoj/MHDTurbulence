/**
 * @file mhd.cpp
 * @brief 
 * @author Tomoya Takiwaki
 * @date 2025-08-21
*/

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <omp.h>

#include "config.hpp"
#include "mhd.hpp"

#include "mpi_config.hpp"

using namespace resolution_mod;

namespace hydflux_mod {
// #pragma omp declare target
  GridArray<double> G;
  FieldArray<double> U; //! U(mconsv,ktot,jtot,itot)
  FieldArray<double> Fx,Fy,Fz;
  FieldArray<double> P; //! P(nprim,ktot,jtot,itot)
  double chg;
  double gam = 1.4;// adiabatic
  //  double csiso; // isothermal
   
// #pragma omp end declare target


void AllocateHydroVariables(GridArray<double>& G,FieldArray<double>& U,FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P){
  
  int dev = omp_get_default_device();
  (void)dev;
  
  G.allocate(ktot,jtot,itot);
  
  for (int i=0; i<itot; i++) {
    G.x1a(i) = 0.0;
    G.x1b(i) = 0.0;
  }
  for (int j=0; j<jtot; j++) {
    G.x2a(j) = 0.0;
    G.x2b(j) = 0.0;
  }
  for (int k=0; k<ktot; k++) {
    G.x3a(k) = 0.0;
    G.x3b(k) = 0.0;
  }
  
  U.allocate(mconsv,ktot,jtot,itot);
  P.allocate(nprim ,ktot,jtot,itot);
      
  Fx.allocate(mconsv,ktot,jtot,itot);
  Fy.allocate(mconsv,ktot,jtot,itot);
  Fz.allocate(mconsv,ktot,jtot,itot);

  for (int m=0; m<mconsv; m++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
  	   U(m,k,j,i) = 0.0;
	  Fx(m,k,j,i) = 0.0;
	  Fy(m,k,j,i) = 0.0;
	  Fz(m,k,j,i) = 0.0;
  }

  for (int n=0; n<nprim; n++)
    for (int k=0; k<ktot; k++)
      for (int j=0; j<jtot; j++)
	for (int i=0; i<itot; i++) {
	  P(n,k,j,i) = 0.0;
  }
  
}

void DeallocateHydroVariables(GridArray<double>& G,FieldArray<double>& U,FieldArray<double>& Fx,FieldArray<double>& Fy,FieldArray<double>& Fz,FieldArray<double>& P){
	G .deallocate();
	U .deallocate();
	Fx.deallocate();
	Fy.deallocate();
	Fz.deallocate();
	P .deallocate();

}

KOKKOS_FUNCTION
void vanLeer(const double& dsvp,const double& dsvm,double& dsv){
    if(dsvp * dsvm > 0.0e0){
      dsv = 2.0e0*dsvp *dsvm/(dsvp + dsvm);
    }else{
      dsv = 0.0e0;
    }
}

KOKKOS_FUNCTION
void HLLD(const double (&leftst)[2*mconsv+madd],const double (&rigtst)[2*mconsv+madd],double (&nflux)[mconsv]){

//=====================================================================
//
// HLLD Scheme
//
// Purpose
// Calculation of Numerical Flux by HLLD method
//
// Reference
//
// Input
// Output
//=====================================================================

//----- U -----
// qql :: left state
// qqr :: right state
      double  rol,vxl,vyl,vzl,ptl,eel;
      double  ror,vxr,vyr,vzr,ptr,eer;
      double  rxl,ryl,rzl;
      double  rxr,ryr,rzr;
      double  bxs,byl,bzl;
      double      byr,bzr;
      double  ptst;
      double scl[ncomp],scr[ncomp];
//----- U* ----
// qqlst ::  left state
// qqrst :: right state
      double  rolst,vxlst,vylst,vzlst,eelst;
      double  rorst,vxrst,vyrst,vzrst,eerst;
      double  rxlst,rylst,rzlst;
      double  rxrst,ryrst,rzrst;
      double        bylst,bzlst;
      double        byrst,bzrst;
      double sclst[ncomp],scrst[ncomp];
      
//----- U** ----
// qqlst ::  left state
// qqrst :: right state
      double  vyldst,vzldst,eeldst;
      double  vyrdst,vzrdst,eerdst;
      double  ryldst,rzldst;
      double  ryrdst,rzrdst;
      double        byldst,bzldst;
      double        byrdst,bzrdst;
      
//----- flux ---
// fqql ::  left physical flux
// fqqr :: right physical flux
      double  frol,frxl,fryl,frzl,feel;
      double            fbyl,fbzl;
      double  fror,frxr,fryr,frzr,feer;
      double            fbyr,fbzr;
      double fscl[ncomp],fscr[ncomp];

//----- wave speed ---
// sl ::  left-going fastest signal velocity
// sr :: right-going fastest signal velocity
// sm :: contact discontinuity velocity
// slst ::  left-going alfven velocity
// srst :: right-going alfven velocity
      double  sm,sl,sr,slst,srst;

// cfl :: left-state Fast wave velocity
// cfr :: right-sate Fast wave velocity
      double  cfl,cfr;

//--------------------
// temporary variables
      double  sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr;
      double  temp;
  
// no if
      double  sign1,maxs1,mins1;
      double  msl,msr,mslst,msrst,temp1,invsumro,sqrtror,sqrtrol,abbx;
      double  bxsq,temp_fst,eps,itf,vdbstl,vdbstr,signbx;

//----- Step 0. ----------------------------------------------------------|
      eps = 1.0e-30;
//---- Left state
        
        rol = leftst[mudn];
        eel = leftst[muet];
        rxl = leftst[muvu];
        ryl = leftst[muvv];
        rzl = leftst[muvw];
        vxl = leftst[muvu]/leftst[mudn];
        vyl = leftst[muvv]/leftst[mudn];
        vzl = leftst[muvw]/leftst[mudn];
        byl = leftst[mubv];
        bzl = leftst[mubw];
        ptl = leftst[mpre];
	for(int n=0;n<ncomp;n++){
	  scl[n] = leftst[mst+n];
	}
//---- Right state
        
        ror = rigtst[mudn];
        eer = rigtst[muet];
        rxr = rigtst[muvu];
        ryr = rigtst[muvv];
        rzr = rigtst[muvw];
        vxr = rigtst[muvu]/rigtst[mudn];
        vyr = rigtst[muvv]/rigtst[mudn];
        vzr = rigtst[muvw]/rigtst[mudn];
        byr = rigtst[mubv];
        bzr = rigtst[mubw];
        ptr = rigtst[mpre];
	for(int n=0;n<ncomp;n++){
	  scr[n] = rigtst[mst+n];
	}
//----- Step 1. ----------------------------------------------------------|
// Compute wave left & right wave speed
//
        cfl = leftst[mcsp];
        cfr = rigtst[mcsp];

        sl = std::fmin(vxl,vxr)-std::fmax(cfl,cfr); // note sl is negative
        sr = std::fmax(vxl,vxr)+std::fmax(cfl,cfr);
//----- Step 2. ----------------------------------------------------------|
// compute L/R fluxs
//
// Left value
        frol = leftst[mfdn];
        feel = leftst[mfet];
        frxl = leftst[mfvu];
        fryl = leftst[mfvv];
        frzl = leftst[mfvw];
        fbyl = leftst[mfbv];
        fbzl = leftst[mfbw];
	for(int n=0;n<ncomp;n++){
	  fscl[n] = leftst[mfst+n];
	}

// Right value
        fror = rigtst[mfdn];
        feer = rigtst[mfet];
        frxr = rigtst[mfvu];
        fryr = rigtst[mfvv]; 
        frzr = rigtst[mfvw];
        fbyr = rigtst[mfbv];
        fbzr = rigtst[mfbw];

	for(int n=0;n<ncomp;n++){
	  fscr[n] = rigtst[mfst+n];
	}

//----- Step 4. ----------------------------------------------------------|
// compute middle and alfven wave
//
        sdl = std::fmin(-1e-20,sl - vxl);
        sdr = std::fmax( 1e-20,sr - vxr);
        rosdl = rol*sdl;
        rosdr = ror*sdr;

        temp = 1.0e0/(rosdr - rosdl);
// Eq. 45
        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp;
           
        sdml = std::fmin(-1e-20,sl - sm); isdml = 1.0e0/sdml;
        sdmr = std::fmax( 1e-20,sr - sm); isdmr = 1.0e0/sdmr;

//----- Step 5. ----------------------------------------------------------|
// compute intermediate states
//
// Eq. 49
        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp;
		
//----- Step 5A. ----------------------------------------------------------|
// compute Ul*
//
           bxs = 0.5e0*(leftst[mubu]+rigtst[mubu]);
           bxsq = bxs*bxs;
           temp_fst = rosdl*sdml - bxsq;
           sign1 = std::copysign(1.0e0,std::abs(temp_fst)-eps);

           maxs1 = std::fmax(0.0e0,sign1);
           mins1 = std::fmin(0.0e0,sign1);

           itf = 1.0e0/(temp_fst+mins1);
           isdml = 1.0e0/sdml;

           temp = bxs*(sdl-sdml)*itf;
           rolst = maxs1*(rosdl*isdml) - mins1*rol;
	   for(int n=0;n<ncomp;n++){
	     sclst[n] = rolst/rol*scl[n];
	   }
           vxlst = maxs1*sm - mins1*vxl;
           rxlst = rolst*vxlst;
           
           vylst = maxs1*(vyl-byl*temp) - mins1*vyl;
           rylst = rolst*vylst;
           vzlst = maxs1*(vzl-bzl*temp) - mins1*vzl;
           rzlst = rolst*vzlst;
           
           temp = (rosdl*sdl-bxsq)*itf;
           bylst = maxs1*(byl*temp) - mins1*byl;
           bzlst = maxs1*(bzl*temp) - mins1*bzl;

           vdbstl = vxlst*bxs+vylst*bylst+vzlst*bzlst;
           eelst = maxs1*(sdl*eel - ptl*vxl + ptst*sm +      
               bxs*(vxl*bxs+vyl*byl+vzl*bzl-vdbstl))*isdml  
               - mins1*eel		;
           
//----- Step 5B. ----------------------------------------------------------|
// compute Ur*
//
           temp_fst = rosdr*sdmr - bxsq;
           sign1 = std::copysign(1.0e0,std::abs(temp_fst)-eps);
           maxs1 = std::fmax(0.0e0,sign1);
           mins1 = std::fmin(0.0e0,sign1);

           itf = 1.0e0/(temp_fst+mins1);
           isdmr = 1.0e0/sdmr;
           
           temp = bxs*(sdr-sdmr)*itf;
           rorst = maxs1*(rosdr*isdmr) - mins1*ror;
	   for(int n=0;n<ncomp;n++){
	     scrst[n] = rorst/ror*scr[n];
	   }
           vxrst = maxs1*sm - mins1*vxr;
           rxrst = rorst*vxrst;
           
           vyrst = maxs1*(vyr-byr*temp) - mins1*vyr;
           ryrst = rorst*vyrst;
           vzrst = maxs1*(vzr-bzr*temp) - mins1*vzr;
           rzrst = rorst*vzrst;
           
           temp = (rosdr*sdr-bxsq)*itf;
           byrst = maxs1*(byr*temp) - mins1*byr;
           bzrst = maxs1*(bzr*temp) - mins1*bzr;
				
           vdbstr = vxrst*bxs+vyrst*byrst+vzrst*bzrst;
           eerst = maxs1*((sdr*eer - ptr*vxr  + ptst*sm)*isdmr +   
               bxs*(vxr*bxs+vyr*byr+vzr*bzr-vdbstr)*isdmr)        
               - mins1*eer;
;
//----- Step 5C. ----------------------------------------------------------|;
// compute Ul** and Ur**;
//;
           sqrtrol = sqrt(rolst);
           sqrtror = sqrt(rorst);

           abbx = std::abs(bxs);
           signbx = std::copysign(1.0e0,bxs);
           sign1  = std::copysign(1.0e0,abbx-eps);

           maxs1 =  std::fmax(0e0,sign1);
           mins1 = -std::fmin(0e0,sign1);
           invsumro = maxs1/(sqrtrol + sqrtror);

           temp = invsumro*(sqrtrol*vylst + sqrtror*vyrst
               + signbx*(byrst-bylst));
           vyldst = vylst*mins1 + temp;
           vyrdst = vyrst*mins1 + temp;
	   (void)vyrdst;
	   [[maybe_unused]] auto vyrdst_ = vyrdst;
           ryldst = rylst*mins1 + rolst * temp;
           ryrdst = ryrst*mins1 + rorst * temp;

           temp = invsumro*(sqrtrol*vzlst + sqrtror*vzrst  
               + signbx*(bzrst-bzlst));
           vzldst = vzlst*mins1 + temp;
           vzrdst = vzrst*mins1 + temp;
	   (void)vzrdst;
	   [[maybe_unused]] auto vzrdst_ = vzrdst;
           rzldst = rzlst*mins1 + rolst * temp;
           rzrdst = rzrst*mins1 + rorst * temp;

           temp = invsumro*(sqrtrol*byrst + sqrtror*bylst
               + signbx*sqrtrol*sqrtror*(vyrst-vylst));
           byldst = bylst*mins1 + temp;
           byrdst = byrst*mins1 + temp;

           temp = invsumro*(sqrtrol*bzrst + sqrtror*bzlst
                + signbx*sqrtrol*sqrtror*(vzrst-vzlst));
           bzldst = bzlst*mins1 + temp;
           bzrdst = bzrst*mins1 + temp;
           
           temp = sm*bxs + vyldst*byldst + vzldst*bzldst;
           eeldst = eelst - sqrtrol*signbx*(vdbstl - temp)*maxs1;
           eerdst = eerst + sqrtror*signbx*(vdbstr - temp)*maxs1;
           
//----- Step 6. ----------------------------------------------------------|;
// compute flux;
           slst = (sm - abbx/sqrtrol)*maxs1;
           srst = (sm + abbx/sqrtror)*maxs1;

           sign1 =  std::copysign(1.0e0,sm);
           maxs1 =  std::fmax(0.0e0,sign1);
           mins1 = -std::fmin(0.0e0,sign1);

           msl = std::fmin(sl,0.0e0);
           msr = std::fmax(sr,0.0e0);
           mslst = std::fmin(slst,0.0e0);
           msrst = std::fmax(srst,0.0e0);

           temp = mslst-msl;
           temp1 = msrst-msr;
	   (void)temp1;
	   [[maybe_unused]] auto temp1_ = temp1;

           nflux[mden] = (frol+(rolst-rol)*msl)*maxs1
                        +(fror+(rorst-ror)*msr)*mins1;
           nflux[meto] = (feel+(eelst-eel)*msl+(eeldst-eelst)*mslst)*maxs1
                        +(feer+(eerst-eer)*msr+(eerdst-eerst)*msrst)*mins1;
           nflux[mrvu] = (frxl+(rxlst-rxl)*msl)*maxs1
                        +(frxr+(rxrst-rxr)*msr)*mins1;
           nflux[mrvv] = (fryl+(rylst-ryl)*msl+(ryldst-rylst)*mslst)*maxs1
                        +(fryr+(ryrst-ryr)*msr+(ryrdst-ryrst)*msrst)*mins1;
           nflux[mrvw] = (frzl+(rzlst-rzl)*msl+(rzldst-rzlst)*mslst)*maxs1
                        +(frzr+(rzrst-rzr)*msr+(rzrdst-rzrst)*msrst)*mins1;
           nflux[mbmu] = 0e0;
           nflux[mbmv] = (fbyl+(bylst-byl)*msl+(byldst-bylst)*mslst)*maxs1
                        +(fbyr+(byrst-byr)*msr+(byrdst-byrst)*msrst)*mins1;
           nflux[mbmw] = (fbzl+(bzlst-bzl)*msl+(bzldst-bzlst)*mslst)*maxs1
                        +(fbzr+(bzrst-bzr)*msr+(bzrdst-bzrst)*msrst)*mins1;
	   for (int n=0;n<ncomp;n++){
	     nflux[mst+n] = (fscl[n]+(sclst[n]-scl[n])*msl)*maxs1
	                   +(fscr[n]+(scrst[n]-scr[n])*msr)*mins1;
	   }

}

void GetNumericalFluxD(
		const int idir, 
		const GridArray<double>&G,
		const FieldArray<double>& P,
		FieldArray<double> &F)
{
	auto chg = hydflux_mod::chg;
	auto Norm2 = KOKKOS_LAMBDA(const double &x, const double &y, const double &z){
		return x*x + y*y + z*z;
	};

	auto Prim2ConsD = KOKKOS_LAMBDA (
			const double (&Prim)[nprim], 
			double (&Cons)[2*mconsv+madd],
			const int muv1,
			const int muv2,
			const int muv3,
			const int mub1,
			const int mub2,
			const int mub3,
			const int mfv1,
			const int mfv2,
			const int mfv3,
			const int nved,
			const int nbmd,
			const int nvef, // f = d + 1 (mod 3)
			const int nbmf,
			const int nveb, // b = d - 1 (mod 3)
			const int nbmb)
	{
		const double rho = Prim[nden];
		const double rhoinv = 1.0 / rho;
		const double vsq = Norm2(Prim[nve1], Prim[nve2], Prim[nve3]);
		const double bsq = Norm2(Prim[nbm1], Prim[nbm2], Prim[nbm3]);
		// direction dependent u,v,w
		Cons[mudn] = rho; // rho
		Cons[muv1] = rho*Prim[nve1]; // rho v_x
		Cons[muv2] = rho*Prim[nve2]; // rho v_y
		Cons[muv3] = rho*Prim[nve3]; // rho v_z
		Cons[muet] = rho*Prim[nene]  // e_i
			+0.5e0*rho*vsq  // + rho v^2/2
			+0.5e0*    bsq; // + B^2/2

		// direction dependent u,v,w
		Cons[mub1] = Prim[nbm1]; // b_x
		Cons[mub2] = Prim[nbm2]; // b_y
		Cons[mub3] = Prim[nbm3]; // b_z
		Cons[mubp] = Prim[nbps]; // psi
		for(int n=0; n<ncomp;n++){
			Cons[must+n] = rho*Prim[nst+n]; // composition
		}
		// total pressure
		double  ptl = Prim[npre] + bsq/2.0e0;

		// direction dependent, nve1 or nbm1
		// direction dependent mfv-u,v,w
		Cons[mfdn] =   rho           *Prim[nved];
		Cons[mfv1] =   rho*Prim[nve1]*Prim[nved] 
		                  -Prim[nbm1]*Prim[nbmd];
		Cons[mfv2] =   rho*Prim[nve2]*Prim[nved]
		                  -Prim[nbm2]*Prim[nbmd];
		Cons[mfv3] =   rho*Prim[nve3]*Prim[nved]
		                  -Prim[nbm3]*Prim[nbmd];
		Cons[mfet] = (Cons[muet]+ptl)*Prim[nved]
		                         -bsq*Prim[nbmd];
		Cons[mfvu] += ptl; // p diagnonal

		// direction dependent 2, 1, 3, 1
		Cons[mfbu] =  0.0e0;
		Cons[mfbv] =  Prim[nbmf]*Prim[nved]
		             -Prim[nvef]*Prim[nbmd];
		Cons[mfbw] =  Prim[nbmb]*Prim[nved]
		             -Prim[nveb]*Prim[nbmd];
		Cons[mfbp] = 0.0e0;  // psi

		for(int n=0; n<ncomp;n++){
			// direction dependent, nve1
			Cons[mfst+n] = rho*Prim[nst+n]*Prim[nved]; // composition
		}
		double css = Prim[ncsp]*Prim[ncsp];
		double cts =  css // c_s^2*c_a^2;
			+ bsq*rhoinv;

		// direction dependent, nve1 or nbm1
		Cons[mcsp] = sqrt((cts +sqrt(cts*cts
						-4.0e0*css*Prim[nbmd]*Prim[nbmd]  //direction dependent
						*rhoinv)   
				  )/2.0e0);
		Cons[mvel] = Prim[nved]; //direction dependent
		Cons[mpre] = ptl;
	};

	auto Prim2Cons1 = KOKKOS_LAMBDA (
			const double (&Prim)[nprim], 
			double (&Cons)[2*mconsv+madd])
	{
		// u,v,w = x,y,z
		Prim2ConsD(Prim, Cons,
				muvu, muvv, muvw,
				mubu, mubv, mubw,
				mfvu, mfvv, mfvw,
				nve1, nbm1,
				nve2, nbm2,
				nve3, nbm3);
	};

	auto Prim2Cons2 = KOKKOS_LAMBDA (
			const double (&Prim)[nprim], 
			double (&Cons)[2*mconsv+madd])
	{
		// u,v,w = y,z,x
		Prim2ConsD(Prim, Cons,
				muvw, muvu, muvv,
				mubw, mubu, mubv,
				mfvw, mfvu, mfvv,
				nve2, nbm2,
				nve3, nbm3,
				nve1, nbm1);
	};

	auto Prim2Cons3 = KOKKOS_LAMBDA (
			const double (&Prim)[nprim], 
			double (&Cons)[2*mconsv+madd])
	{
		// u,v,w = z,x,y
		Prim2ConsD(Prim, Cons,
				muvv, muvw, muvu,
				mubv, mubw, mubu,
				mfvv, mfvw, mfvu,
				nve3, nbm3,
				nve1, nbm1,
				nve2, nbm2);
	};

	auto CalcFlux = KOKKOS_LAMBDA (
			const double  (&Pleftc1)[nprim], 
			const double  (&Pleftc2)[nprim], 
			const double  (&Prigtc1)[nprim], 
			const double  (&Prigtc2)[nprim], 
			double (&numflux)[mconsv], 
			double (&Clefte)[2*mconsv+madd], 
			double (&Crigte)[2*mconsv+madd],
			const int idir)
	{
		// Calculte Left state
		double Plefte [nprim];
		/* | Pleftc1   | Pleftc2 =>| Prigtc1   | Prigtc2   |  */
		/*                     You are here                   */
		for (int n=0; n<nprim; n++){
			double dsvp =  Prigtc1[n]- Pleftc2[n];
			double dsvm =              Pleftc2[n]- Pleftc1[n];
			double dsv;
			vanLeer(dsvp,dsvm,dsv);
			Plefte[n] = Pleftc2[n] + 0.5e0*dsv;
		}
		if(1 == idir )Prim2Cons1(Plefte, Clefte);
		if(2 == idir )Prim2Cons2(Plefte, Clefte);
		if(3 == idir )Prim2Cons3(Plefte, Clefte);

		// Calculte Right state

		double Prigte [nprim];
		/* | Pleftc1   | Pleftc2 |<= Prigtc1   | Prigtc2   |  */
		/*                     You are here                   */
		for (int n=0; n<nprim; n++){
			double dsvp =  Prigtc2[n]- Prigtc1[n];
			double dsvm =              Prigtc1[n]- Pleftc2[n];
			double dsv;
			vanLeer(dsvp,dsvm,dsv);
			Prigte[n] = Prigtc1[n] - 0.5e0*dsv;
		}
		if(1 == idir )Prim2Cons1(Prigte, Crigte);
		if(2 == idir )Prim2Cons2(Prigte, Crigte);
		if(3 == idir )Prim2Cons3(Prigte, Crigte);

		HLLD(Clefte, Crigte, numflux);
	};

	if(1 == idir){
#if 0
#pragma omp target teams distribute parallel for collapse(3)
		for (int k=ks; k<=ke; k++) {
			for (int j=js; j<=je; j++){
				for (int i=is; i<=ie+1; i++) {
					double Pleftc1[nprim];
					double Pleftc2[nprim];
					double Prigtc1[nprim];
					double Prigtc2[nprim];  

					double numflux [mconsv];

					double Clefte [2*mconsv+madd];
					double Crigte [2*mconsv+madd];

					for (int n=0; n<nprim; n++){
						Pleftc1[n] = P(n,k,j,i-2);
						Pleftc2[n] = P(n,k,j,i-1);
						Prigtc1[n] = P(n,k,j,i  );
						Prigtc2[n] = P(n,k,j,i+1);
					}
					CalcFlux(Pleftc1, Pleftc2, Prigtc1, Prigtc2, numflux, Clefte, Crigte, Prim2Cons1);

					F(mden,k,j,i) = numflux[mden];
					F(mrv1,k,j,i) = numflux[mrvu];
					F(mrv2,k,j,i) = numflux[mrvv];
					F(mrv3,k,j,i) = numflux[mrvw];
					F(meto,k,j,i) = numflux[meto];
					//F(mbm1,k,j,i) = numflux[mbmu];
					F(mbm2,k,j,i) = numflux[mbmv];
					F(mbm3,k,j,i) = numflux[mbmw];

					F(mbm1,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
					               -0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
					F(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
					               -0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
					for(int n=0; n<ncomp;n++){
						F(mst+n,k,j,i) = numflux[mst+n]; // composition
					}
				}
			}
		}
#else
		Kokkos::parallel_for("Flux1",
		Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+2, je+1, ke+1}),
		KOKKOS_LAMBDA(const int i, const int j, const int k) {
			double Pleftc1[nprim];
			double Pleftc2[nprim];
			double Prigtc1[nprim];
			double Prigtc2[nprim];  

			double numflux [mconsv];

			double Clefte [2*mconsv+madd];
			double Crigte [2*mconsv+madd];

			for (int n=0; n<nprim; n++){
			Pleftc1[n] = P.dev(n,k,j,i-2);
			Pleftc2[n] = P.dev(n,k,j,i-1);
			Prigtc1[n] = P.dev(n,k,j,i  );
			Prigtc2[n] = P.dev(n,k,j,i+1);
			}
			CalcFlux(Pleftc1, Pleftc2, Prigtc1, Prigtc2, numflux, Clefte, Crigte, 1);

			F.dref(mden,k,j,i) = numflux[mden];
			F.dref(mrv1,k,j,i) = numflux[mrvu];
			F.dref(mrv2,k,j,i) = numflux[mrvv];
			F.dref(mrv3,k,j,i) = numflux[mrvw];
			F.dref(meto,k,j,i) = numflux[meto];
			//F.dref(mbm1,k,j,i) = numflux[mbmu];
			F.dref(mbm2,k,j,i) = numflux[mbmv];
			F.dref(mbm3,k,j,i) = numflux[mbmw];

			F.dref(mbm1,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
				-0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
			F.dref(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
					-0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
			for(int n=0; n<ncomp;n++){
				F.dref(mst+n,k,j,i) = numflux[mst+n]; // composition
			}
		});
#endif
	}

	if(2 == idir){
#if 0
#pragma omp target teams distribute parallel for collapse(3)
		for (int k=ks; k<=ke; k++){
			for (int j=js; j<=je+1; j++){
				for (int i=is; i<=ie; i++) {
					double Pleftc1[nprim];
					double Pleftc2[nprim];
					double Prigtc1[nprim];
					double Prigtc2[nprim];	  
					for (int n=0; n<nprim; n++){
						Pleftc1[n] = P(n,k,j-2,i);
						Pleftc2[n] = P(n,k,j-1,i);
						Prigtc1[n] = P(n,k,j  ,i);
						Prigtc2[n] = P(n,k,j+1,i);
					}

					double numflux [mconsv];
					double Clefte [2*mconsv+madd];
					double Crigte [2*mconsv+madd];
					CalcFlux(Pleftc1, Pleftc2, Prigtc1, Prigtc2, numflux, Clefte, Crigte, Prim2Cons2);

					F(mden,k,j,i) = numflux[mden];
					F(mrv1,k,j,i) = numflux[mrvw];
					F(mrv2,k,j,i) = numflux[mrvu];
					F(mrv3,k,j,i) = numflux[mrvv];
					F(meto,k,j,i) = numflux[meto];
					F(mbm1,k,j,i) = numflux[mbmw];
					//F(mbm2,k,j,i) = numflux[mbmu];
					F(mbm3,k,j,i) = numflux[mbmv];

					F(mbm2,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
						-0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
					F(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
							-0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
					for(int n=0; n<ncomp;n++){
						F(mst+n,k,j,i) = numflux[mst+n]; // composition
					}
				}// j-loop
			}// k,i-loop
		}
#else
		Kokkos::parallel_for("Flux2",
		Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+1, je+2, ke+1}),
		KOKKOS_LAMBDA(const int i, const int j, const int k) {
			double Pleftc1[nprim];
			double Pleftc2[nprim];
			double Prigtc1[nprim];
			double Prigtc2[nprim];	  
			for (int n=0; n<nprim; n++){
			Pleftc1[n] = P.dev(n,k,j-2,i);
			Pleftc2[n] = P.dev(n,k,j-1,i);
			Prigtc1[n] = P.dev(n,k,j  ,i);
			Prigtc2[n] = P.dev(n,k,j+1,i);
			}

			double numflux [mconsv];
			double Clefte [2*mconsv+madd];
			double Crigte [2*mconsv+madd];
			CalcFlux(Pleftc1, Pleftc2, Prigtc1, Prigtc2, numflux, Clefte, Crigte, 2);

			F.dref(mden,k,j,i) = numflux[mden];
			F.dref(mrv1,k,j,i) = numflux[mrvw];
			F.dref(mrv2,k,j,i) = numflux[mrvu];
			F.dref(mrv3,k,j,i) = numflux[mrvv];
			F.dref(meto,k,j,i) = numflux[meto];
			F.dref(mbm1,k,j,i) = numflux[mbmw];
			//F.dref(mbm2,k,j,i) = numflux[mbmu];
			F.dref(mbm3,k,j,i) = numflux[mbmv];

			F.dref(mbm2,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
				-0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
			F.dref(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
					-0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
			for(int n=0; n<ncomp;n++){
				F.dref(mst+n,k,j,i) = numflux[mst+n]; // composition
			}
		});
#endif
	}

	if(3 == idir){
#if 0
#pragma omp target teams distribute parallel for collapse(3)
		for(int j=js; j<=je; ++j){
			for(int k=ks; k<=ke+1; ++k){
				for(int i=is; i<=ie; ++i){
					double Pleftc1[nprim];
					double Pleftc2[nprim];
					double Prigtc1[nprim];
					double Prigtc2[nprim];  

					double numflux [mconsv];

					double Clefte [2*mconsv+madd];
					double Crigte [2*mconsv+madd];

					for (int n=0; n<nprim; n++){
						Pleftc1[n] = P(n,k-2,j,i);
						Pleftc2[n] = P(n,k-1,j,i);
						Prigtc1[n] = P(n,k  ,j,i);
						Prigtc2[n] = P(n,k+1,j,i);
					}
					CalcFlux(Pleftc1, Pleftc2, Prigtc1, Prigtc2, numflux, Clefte, Crigte, Prim2Cons3);

					F(mden,k,j,i) = numflux[mden];
					F(mrv1,k,j,i) = numflux[mrvv];
					F(mrv2,k,j,i) = numflux[mrvw];
					F(mrv3,k,j,i) = numflux[mrvu];
					F(meto,k,j,i) = numflux[meto];
					F(mbm1,k,j,i) = numflux[mbmv];
					F(mbm2,k,j,i) = numflux[mbmw];
					//F(mbm3,k,j,i) = numflux[mbm];

					F(mbm3,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
					                -0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
					F(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
					               -0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
					for(int n=0; n<ncomp;n++){
						F(mst+n,k,j,i) = numflux[mst+n]; // composition
					}
				}
			}
		}
#else
		Kokkos::parallel_for("Flux3",
		Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, ks, js}, {ie+1, ke+2, je+1}),
		KOKKOS_LAMBDA(const int i, const int k, const int j) {
			double Pleftc1[nprim];
			double Pleftc2[nprim];
			double Prigtc1[nprim];
			double Prigtc2[nprim];  

			double numflux [mconsv];

			double Clefte [2*mconsv+madd];
			double Crigte [2*mconsv+madd];

			for (int n=0; n<nprim; n++){
				Pleftc1[n] = P.dev(n,k-2,j,i);
				Pleftc2[n] = P.dev(n,k-1,j,i);
				Prigtc1[n] = P.dev(n,k  ,j,i);
				Prigtc2[n] = P.dev(n,k+1,j,i);
			}
			CalcFlux(Pleftc1, Pleftc2, Prigtc1, Prigtc2, numflux, Clefte, Crigte, 3);

			F.dref(mden,k,j,i) = numflux[mden];
			F.dref(mrv1,k,j,i) = numflux[mrvv];
			F.dref(mrv2,k,j,i) = numflux[mrvw];
			F.dref(mrv3,k,j,i) = numflux[mrvu];
			F.dref(meto,k,j,i) = numflux[meto];
			F.dref(mbm1,k,j,i) = numflux[mbmv];
			F.dref(mbm2,k,j,i) = numflux[mbmw];
			//F.dref(mbm3,k,j,i) = numflux[mbm];

			F.dref(mbm3,k,j,i) = 0.5e0*    (Clefte[mubp]+Crigte[mubp])
					-0.5e0*chg*(Crigte[mubu]-Clefte[mubu]);
			F.dref(mbps,k,j,i) =(0.5e0*    (Clefte[mubu]+Crigte[mubu])
				       -0.5e0/chg*(Crigte[mubp]-Clefte[mubp]))*chg*chg;
			for(int n=0; n<ncomp;n++){
				F.dref(mst+n,k,j,i) = numflux[mst+n]; // composition
			}
		});
#endif
	}
}

void UpdateConservU(const GridArray<double>& G,const FieldArray<double>& Fx,const FieldArray<double>& Fy,const FieldArray<double>& Fz,FieldArray<double>& U){
	auto dt = resolution_mod::dt;

  //printf("pre  U:%e %e %e %e\n", U(mden,ks,js,is), U(mrv3,ks,js,is), U(mbm3,ks,js,is), U(meto,ks,js,is));
  //printf("pre Fx:%e %e %e %e\n",Fx(mden,ks,js,is),Fx(mrv3,ks,js,is),Fx(mbm3,ks,js,is),Fx(meto,ks,js,is));
  //printf("pre Fy:%e %e %e %e\n",Fy(mden,ks,js,is),Fy(mrv3,ks,js,is),Fy(mbm3,ks,js,is),Fy(meto,ks,js,is));
  //printf("pre Fz:%e %e %e %e\n",Fz(mden,ks,js,is),Fz(mrv3,ks,js,is),Fz(mbm3,ks,js,is),Fz(meto,ks,js,is));
  
#if 0
#pragma omp target teams distribute parallel for collapse(4)
  for (int m=0; m<mconsv; m++)
    for (int k=ks; k<=ke; k++)
      for (int j=js; j<=je; j++)
	for (int i=is; i<=ie; i++) {
	  U(m,k,j,i) -= dt * ( (Fx(m,k,j,i+1) - Fx(m,k,j,i)) / (G.x1a(i+1)-G.x1a(i))
			      +(Fy(m,k,j+1,i) - Fy(m,k,j,i)) / (G.x2a(j+1)-G.x2a(j))
			      +(Fz(m,k+1,j,i) - Fz(m,k,j,i)) / (G.x3a(k+1)-G.x3a(k)) );
  }
  //printf("aft U:%e %e %e %e\n",U(mden,ks,js,is),U(mrv1,ks,js,is),U(mbm1,ks,js,is),U(meto,ks,js,is));
#else
	Kokkos::parallel_for("UpdateConservU",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+1, je+1, ke+1}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		const double xinv = 1.0 / (G.dev_x1a(i+1)-G.dev_x1a(i));
		const double yinv = 1.0 / (G.dev_x2a(i+1)-G.dev_x2a(i));
		const double zinv = 1.0 / (G.dev_x3a(i+1)-G.dev_x3a(i));
		for (int m=0; m<mconsv; m++) {
			U.dref(m,k,j,i) -= dt * ( (Fx.dev(m,k,j,i+1) - Fx.dev(m,k,j,i)) * xinv
			                         +(Fy.dev(m,k,j+1,i) - Fy.dev(m,k,j,i)) * yinv
			                         +(Fz.dev(m,k+1,j,i) - Fz.dev(m,k,j,i)) * zinv );
		}
	});
#endif
}


void UpdatePrimitvP(const FieldArray<double>& U,FieldArray<double>& P){
	auto gam = hydflux_mod::gam;
  	auto Norm2 = KOKKOS_LAMBDA(double x, double y, double z) { return x*x + y*y + z*z; };

  //printf("U:et b1 b2 b3=%e %e %e %e\n",U(meto,ks,js,is),U(mbm1,ks,js,is),U(mbm2,ks,js,is),U(mbm3,ks,js,is));
#if 0
#pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) {
         P(nden,k,j,i) = U(mden,k,j,i);
	 P(nve1,k,j,i) = U(mrv1,k,j,i)/U(mden,k,j,i);
         P(nve2,k,j,i) = U(mrv2,k,j,i)/U(mden,k,j,i);
         P(nve3,k,j,i) = U(mrv3,k,j,i)/U(mden,k,j,i);
	 double ekin = 0.5e0*( U(mrv1,k,j,i)*U(mrv1,k,j,i)
			      +U(mrv2,k,j,i)*U(mrv2,k,j,i)
			      +U(mrv3,k,j,i)*U(mrv3,k,j,i))/U(mden,k,j,i);
	 double emag = 0.5e0*( U(mbm1,k,j,i)*U(mbm1,k,j,i)
			      +U(mbm2,k,j,i)*U(mbm2,k,j,i)
			      +U(mbm3,k,j,i)*U(mbm3,k,j,i));
         P(nene,k,j,i) =  (U(meto,k,j,i)-ekin-emag)/U(mden,k,j,i);//specific internal energy
	 //P(npre,k,j,i) =  U(mden,k,j,i) * csiso * csiso;
         P(npre,k,j,i) = P(nene,k,j,i) * P(nden,k,j,i) * (gam-1.0); 
	 //P(ncsp,k,j,i) =  csiso;
	 P(ncsp,k,j,i) = sqrt(P(nene,k,j,i) * gam * (gam-1.0));
	 P(nbm1,k,j,i) =  U(mbm1,k,j,i);
	 P(nbm2,k,j,i) =  U(mbm2,k,j,i);
	 P(nbm3,k,j,i) =  U(mbm3,k,j,i);
	 P(nbps,k,j,i) =  U(mbps,k,j,i);
	 for(int n=0;n<ncomp;n++){
	   P(nst+n,k,j,i) = U(mst+n,k,j,i)/U(mden,k,j,i);
	 }
	 
      }
#else
	Kokkos::parallel_for("UpdatePrimitvP",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+1, je+1, ke+1}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
		const double rho =  U.dev(mden,k,j,i);
		const double rhoinv = 1.0 / rho;
		P.dref(nden,k,j,i) = rho;
		P.dref(nve1,k,j,i) = U.dev(mrv1,k,j,i) * rhoinv;
		P.dref(nve2,k,j,i) = U.dev(mrv2,k,j,i) * rhoinv;
		P.dref(nve3,k,j,i) = U.dev(mrv3,k,j,i) * rhoinv;
		double ekin = 0.5e0*Norm2(U.dev(mrv1,k,j,i), U.dev(mrv2,k,j,i), U.dev(mrv3,k,j,i)) * rhoinv;
		double emag = 0.5e0*Norm2(U.dev(mbm1,k,j,i), U.dev(mbm2,k,j,i), U.dev(mbm3,k,j,i));
		P.dref(nene,k,j,i) =  (U.dev(meto,k,j,i)-ekin-emag)/U.dev(mden,k,j,i);//specific internal energy
		//P.dref(npre,k,j,i) =  U.dev(mden,k,j,i) * csiso * csiso;
		P.dref(npre,k,j,i) = P.dev(nene,k,j,i) * P.dev(nden,k,j,i) * (gam-1.0); 
		//P.dref(ncsp,k,j,i) =  csiso;
		P.dref(ncsp,k,j,i) = sqrt(P.dev(nene,k,j,i) * gam * (gam-1.0));
		P.dref(nbm1,k,j,i) =  U.dev(mbm1,k,j,i);
		P.dref(nbm2,k,j,i) =  U.dev(mbm2,k,j,i);
		P.dref(nbm3,k,j,i) =  U.dev(mbm3,k,j,i);
		P.dref(nbps,k,j,i) =  U.dev(mbps,k,j,i);
		for(int n=0;n<ncomp;n++){
			P.dref(nst+n,k,j,i) = U.dev(mst+n,k,j,i) * rhoinv;
		}
	});
#endif
}

void ControlTimestep(const GridArray<double>& G, const FieldArray<double>& P){
  using namespace mpi_config_mod;
  const double huge = 1.0e90;
  double dtminl = huge;
// #pragma omp target teams distribute parallel for reduction(min:dtminl) collapse(3)
#if 0
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double ctot = sqrt(     P(ncsp,k,j,i)*P(ncsp,k,j,i)
			   + (  P(nbm1,k,j,i)*P(nbm1,k,j,i)
			      + P(nbm2,k,j,i)*P(nbm2,k,j,i)
			      + P(nbm3,k,j,i)*P(nbm3,k,j,i) ) / P(nden,k,j,i)); 
	double dtminloc = std::min({ (G.x1a(i+1)-G.x1a(i))/(std::abs(P(nve1,k,j,i))+ctot)
				    ,(G.x2a(j+1)-G.x2a(j))/(std::abs(P(nve2,k,j,i))+ctot)
				    ,(G.x3a(k+1)-G.x3a(k))/(std::abs(P(nve3,k,j,i))+ctot)});
	dtminl = std::min(dtminloc,dtminl);
	/*  if(dtminloc < dtmin){//for debug
	  dtmin = dtminloc;
	  ip=i;
	  jp=j;
	  kp=k;
	  }*/
	
      }
#else
  	auto Sqr = KOKKOS_LAMBDA(double x) { return x*x; };
  	auto Norm2 = KOKKOS_LAMBDA(double x, double y, double z) { return x*x + y*y + z*z; };
	// auto P = hydflux_mod::P;
	Kokkos::parallel_reduce("ControlTimestep",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+1, je+1, ke+1}),
	KOKKOS_LAMBDA(const int i, const int j, const int k, double& dtminloc) {
	  double ctot = sqrt( 
		      	  Sqr(P.dev(ncsp,k,j,i))
			  + Norm2(P.dev(nbm1,k,j,i), P.dev(nbm2,k,j,i), P.dev(nbm3,k,j,i)) / P.dev(nden,k,j,i)
			  ); 
#if 0
	  double dt= std::min({ 
			  (G.x1a(i+1)-G.x1a(i))/(std::abs(P(nve1,k,j,i))+ctot),
			  (G.x2a(j+1)-G.x2a(j))/(std::abs(P(nve2,k,j,i))+ctot),
			  (G.x3a(k+1)-G.x3a(k))/(std::abs(P(nve3,k,j,i))+ctot)
		  });
	  dtminloc = std::fmin(dtminloc, dt);
#else
	  dtminloc = Kokkos::min({
			  dtminloc, 
			  (G.dev_x1a(i+1)-G.dev_x1a(i))/(std::abs(P.dev(nve1,k,j,i))+ctot),
			  (G.dev_x2a(j+1)-G.dev_x2a(j))/(std::abs(P.dev(nve2,k,j,i))+ctot),
			  (G.dev_x3a(k+1)-G.dev_x3a(k))/(std::abs(P.dev(nve3,k,j,i))+ctot)
		  });
#endif
	},
	Kokkos::Min<double>(dtminl)
	);
#endif
  //  Here dtminl is in host
  double dtming;
  int myid_wg;
  MPIminfind(dtminl,myid_w,dtming,myid_wg);
  dt = 0.05e0*dtming;
}

void EvaluateCh(const FieldArray<double>& P){
  using namespace mpi_config_mod;
  // auto P = hydflux_mod::P;
  double chgloc = 0.0e0;
// #pragma omp target teams distribute parallel for collapse(3) reduction(max:chgloc)
#if 0
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double css = P(ncsp,k,j,i)*P(ncsp,k,j,i);
	double ca1 = P(nbm1,k,j,i)*P(nbm1,k,j,i) / P(nden,k,j,i);
	double ca2 = P(nbm2,k,j,i)*P(nbm2,k,j,i) / P(nden,k,j,i);
	double ca3 = P(nbm3,k,j,i)*P(nbm3,k,j,i) / P(nden,k,j,i);
	double cts = css + ca1 + ca2 + ca3;
        double cm1 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca1))/2.0e0);
        double ch1 = (std::abs(P(nve1,k,j,i))+cm1);
        double cm2 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca2))/2.0e0);
        double ch2 = (std::abs(P(nve2,k,j,i))+cm2);
        double cm3 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca3))/2.0e0);
        double ch3 = (std::abs(P(nve3,k,j,i))+cm3);

        chgloc = std::max({chgloc,ch1,ch2,ch3});
      }
#else
	Kokkos::parallel_reduce("EvaluateCh",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+1, je+1, ke+1}),
	KOKKOS_LAMBDA(const int i, const int j, const int k, double& chgll) {
		double rhoinv = 1.0 / P.dev(nden,k,j,i);
		double css = P.dev(ncsp,k,j,i)*P.dev(ncsp,k,j,i);
		double ca1 = P.dev(nbm1,k,j,i)*P.dev(nbm1,k,j,i) * rhoinv;
		double ca2 = P.dev(nbm2,k,j,i)*P.dev(nbm2,k,j,i) * rhoinv;
		double ca3 = P.dev(nbm3,k,j,i)*P.dev(nbm3,k,j,i) * rhoinv;
		double cts = css + ca1 + ca2 + ca3;
		double cm1 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca1))/2.0e0);
		double ch1 = (std::abs(P.dev(nve1,k,j,i))+cm1);
		double cm2 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca2))/2.0e0);
		double ch2 = (std::abs(P.dev(nve2,k,j,i))+cm2);
		double cm3 = sqrt((cts+sqrt(cts*cts-4.0e0*css*ca3))/2.0e0);
		double ch3 = (std::abs(P.dev(nve3,k,j,i))+cm3);

		chgll = Kokkos::max({chgll, ch1, ch2, ch3});
	},
	Kokkos::Max<double>(chgloc)
	);
#endif
  // Here chgg is in host
  double chgg;
  int  myid_wg;
  MPImaxfind(chgloc,myid_w,chgg,myid_wg);
  chg = chgg;
}

void DampPsi(const GridArray<double>& G,FieldArray<double>& U){
  const double alphabp = 0.1e0;
  auto chg = hydflux_mod::chg;
  auto dt = resolution_mod::dt;
#if 0
// #pragma omp target teams distribute parallel for collapse(3)
  for (int k=ks; k<=ke; k++)
    for (int j=js; j<=je; j++)
      for (int i=is; i<=ie; i++) {
	double dhl = std::min({G.x1a(i+1)-G.x1a(i),G.x2a(j+1)-G.x2a(j),G.x3a(k+1)-G.x3a(k)});
	double taui = alphabp * chg/dhl;
	  U(mbps,k,j,i) = U(mbps,k,j,i) *(1.0e0-dt*taui);
      }
#else
	Kokkos::parallel_for("DampPsi",
	Kokkos::MDRangePolicy<Kokkos::Rank<3>>({is, js, ks}, {ie+1, je+1, ke+1}),
	KOKKOS_LAMBDA(const int i, const int j, const int k) {
#if 0
		double dhl = std::min({G.x1a(i+1)-G.x1a(i), G.x2a(j+1)-G.x2a(j), G.x3a(k+1)-G.x3a(k)});
#endif
		double dx = G.dev_x1a(i+1)-G.dev_x1a(i);
		double dy = G.dev_x2a(j+1)-G.dev_x2a(j);
		double dz = G.dev_x3a(k+1)-G.dev_x3a(k);

		double dhl = Kokkos::min({dx, dy, dz});
		double taui = alphabp * chg/dhl;
		U.dref(mbps,k,j,i) = U.dev(mbps,k,j,i) * (1.0e0-dt*taui);
	});
#endif
}


}; // end of namespace
