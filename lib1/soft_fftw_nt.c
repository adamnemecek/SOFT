/***************************************************************************
  **************************************************************************
  
  SOFT: SO(3) Fourier Transforms
  Version 2.0

  Copyright (c) 2003, 2004, 2007 Peter Kostelec, Dan Rockmore
  
  This file is part of SOFT.

  SOFT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  SOFT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/

/*

  functions that do the FORWARD and INVERSE transforms on the
  full group SO(3)

  will use fftw

  will write over input data!!!

  won't have to transpose - let fftw do via the plan (to save a little memory)

  sample size = (2*bw)^3
  coefficient size = on the order of bw*(4*bw^2-1)/3

  functions in here:

  Forward_SO3_Naive_fftw_nt();
  Inverse_SO3_Naive_fftw_nt();

*/



/* the following is for memset */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fftw3.h"

#include "utils_vec_cx.h"
#include "utils_so3.h"
#include "fft_grids_so3.h"
#include "makeWigner.h"
#include "wignerTransforms_fftw.h"

/*
  ok, the forward transform

  Function arguments are as follows:

  bw = bandwidth of transform
  data: FFTW_COMPLEX array of size (2*bw)^3 containing the input
        signal ->
	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZEOF
	FFTW_COMPLEX!!! (although I will sometimes treat it
	as the plain COMPLEX array it is)
  coeffs: plain COMPLEX array of size (4*bw^3-bw)/3, will contain the
          coefficients of the signal

  workspace_cx: scratch space FFTW_COMPLEX array of size (2*bw)^3
	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
	FFTW_COMPLEX!!! (although I will sometimes treat it
	as the plain COMPLEX array it is)

  workspace_re: double scratch space of size 12*n + n*bw
		where n = 2*bw

  weights: ptr to double array of size 2*bw - this array holds
	   the PRECOMPUTED quadrature weights

  p1: pointer to FFTW plan for correctly ffting DATA array and
      placing the result in DATA (so it's in-place); this is
      an INVERSE FFT you're doing!!!

  flag: = 0 : data is COMPLEX
        = 1 : data is double


*/

void Forward_SO3_Naive_fftw_nt( int bw,
				fftw_complex *data,
				fftw_complex *coeffs,
				fftw_complex *workspace_cx,
				double *workspace_re,
				double *weights,
				fftw_plan *p1,
				int flag )
{
  int j, n ;
  int m1, m2 ;
  int sampHere, coefHere ;
  int coefHere2 ;
  int tmpInt ;
  double *sinPts, *cosPts, *sinPts2, *cosPts2 ;
  double *wigners, *scratch ;
  double fudge ;
  fftw_complex *coeffsPtr ;
  fftw_complex *dataPtr ;
  double dn ;

  n = 2 * bw ;

  sinPts = workspace_re ;
  cosPts = sinPts + n ;
  sinPts2 = cosPts + n ;
  cosPts2 = sinPts2 + n ;
  wigners = cosPts2 + n ;
  scratch = wigners + ( bw * n ) ; /* wigners need at most bw*n space AT
				      ANY given orders m1, m2 */
  /*
    before going further, let's precompute all the sines
    and cosines I'll need. No matter what order transform
    I'm doing, these'll stay the same.
  */
  SinEvalPts( n, sinPts );
  CosEvalPts( n, cosPts );
  SinEvalPts2( n, sinPts2 );
  CosEvalPts2( n, cosPts2 );

  /*
    Stage 1: FFT
  */
  fftw_execute( *p1 ) ;
  /* this is silly, I know, but I can' figure out how to
     do an in-place fft and transpose at the same time */
  memcpy( data, workspace_cx, sizeof(fftw_complex) * n * n * n) ;
   
  /*
    Stage 2: Do the Wigner transforms. This is the tricky bit.

    Since I'm working with two order indeces, m1 and m2, the
    for-loops will be more intricate than in the case of the
    "ordinary" spherical transform on S^2.

    Also, I will be taking advantage of the symmetries of the
    Wigner-d functions. As long as I keep my signs and flips
    right, the Wigner-d's I precompute for an order (m1, m2)
    transform can generally  be used in seven more transforms:
    (m1,-m2), (m2,m1), (m2,-m1), (-m2,m1), (-m2,-m1), (-m1,m2)
    and (-m1,-m2).

    I say "general" because, of course, I'll be transforming
    at orders (m1,m1), (m1,0) and (0,m1), so I won't get such
    a huge reduction. Still, this should save time.

    If assumptions are made regarding the original input signal,
    e.g. it's strictly real, then one may take advantage of
    symmetries of the big D wigners (i.e. function of all 3
    parameters alpha, beta, gamma) and so simplify the for-loops
    some and hence increase the speed of the program. However,
    the for-loops to follow will make no such assumptions.
    Whether the signal is real or complex, these for-loops will
    handle it.

    The for-loops will be "designed" as follows. They will be
    divided into cases according to the orders:

    0) {f_{0,0}}

    1) for 0 <= m1 <= bw-1
    compute the coefficients
    i)   {f_{ m1, m1}}
    ii)  {f_{-m1,-m1}}
    iii) {f_{-m1, m1}}
    iv)  {f_{ m1,-m1}}

    2) for 1 <= m1 <= bw-1
    compute the coefficients
    i)   {f_{ m1,  0}}
    ii)  {f_{-m1,  0}}
    iii) {f_{  0, m1}}
    iv)  {f_{  0,-m1}}

    3) for 1 <= m1 <= bw-1
    for m1+1 <= m2 <= bw-1
    compute the coefficients
    i)    {f_{ m1, m2}}
    ii)   {f_{-m1,-m2}}
    iii)  {f_{ m1,-m2}}
    iv)   {f_{-m1, m2}}
    v)    {f_{ m2, m1}}
    vi)   {f_{-m2,-m1}}
    vii)  {f_{ m2,-m1}}
    viii) {f_{-m2, m1}}


    Fasten your seatbelt, folks. It's going to be a bumpy ride.

  */

  /***************************/
  /*                         */
  /* {f_{0,0}} coefficient   */
  /*                         */
  /***************************/


  /* compute the wigners I'll need */
  genWig_L2( 0, 0, bw,
	     sinPts, cosPts,
	     sinPts2, cosPts2,
	     wigners, scratch ) ;

  /* now, get the locations of where the
     samples I have to transform are, and
     where the coefficients have to go */
  
  sampHere = sampLoc_so3( 0, 0, bw ) ;
  coefHere = coefLoc_so3( 0, 0, bw ) ;

  /* ok, reset sample, coef ptrs */
  coeffsPtr = coeffs ;
  dataPtr = data ;
  
  /* now advance by the computed amounts */
  dataPtr += sampHere ;
  coeffsPtr += coefHere ;
  
  /* now transform the real and imaginary parts
     of the data */
  
  wigNaiveAnalysis_fftw( 0, 0, bw, dataPtr,
			 wigners, weights,
			 coeffsPtr,
			 workspace_cx ) ;
  

  /*** 0 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {

      /* compute the wigners I'll need */
      genWig_L2( m1, m1, bw,
		 sinPts, cosPts,
		 sinPts2, cosPts2,
		 wigners, scratch ) ;

      /***************************/
      /*                         */
      /* {f_{m1,m1}} coefficient */
      /*                         */
      /***************************/

      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */
      
      sampHere = sampLoc_so3( m1, m1, bw ) ;
      coefHere = coefLoc_so3( m1, m1, bw ) ;
      
      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;
      
      /* now transform the real and imaginary parts
	 of the data */
      
      wigNaiveAnalysis_fftw( m1, m1, bw, dataPtr,
			     wigners, weights,
			     coeffsPtr,
			     workspace_cx ) ;
      
      /*****************************/
      /*                           */
      /* {f_{-m1,-m1}} coefficient */
      /*                           */
      /*****************************/

      if ( flag == 0 ) /* if data is complex */
	{
	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( -m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( -m1, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_fftw( -m1, -m1, bw, dataPtr,
				 wigners, weights,
				 coeffsPtr,
				 workspace_cx ) ;

	}
      else  /* data is real, so use symmetry */
	{
	  coefHere = coefLoc_so3( m1, m1, bw ) ;
	  coefHere2 = coefLoc_so3( -m1, -m1, bw ) ;

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      coeffs[coefHere2+j][0] = coeffs[coefHere+j][0];
	      coeffs[coefHere2+j][1] = -coeffs[coefHere+j][1];
	    }

	}

      /*****************************/
      /*                           */
      /* {f_{-m1,m1}} coefficient  */
      /*                           */
      /*****************************/


      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( -m1, m1, bw ) ;
      coefHere = coefLoc_so3( -m1, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveAnalysis_fftwY( -m1, m1, bw, dataPtr,
			      wigners, weights,
			      coeffsPtr,
			      workspace_cx ) ;
      
      /*****************************/
      /*                           */
      /* {f_{m1,-m1}} coefficient  */
      /*                           */
      /*****************************/

      if ( flag == 0 ) /* data is complex */
	{
	  
	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( m1, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_fftwY( m1, -m1, bw, dataPtr,
				  wigners, weights,
				  coeffsPtr,
				  workspace_cx ) ;
	  
	}
      else /* data is real, so use symmetry */
	{
	  coefHere = coefLoc_so3( -m1, m1, bw );
	  coefHere2 = coefLoc_so3( m1, -m1, bw );

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      coeffs[coefHere2+j][0] = coeffs[coefHere+j][0];
	      coeffs[coefHere2+j][1] = -1.*coeffs[coefHere+j][1];
	    }

	}

    }

  /*** for 1 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      
      /* compute the wigners I'll need */
      genWig_L2( m1, 0, bw,
		 sinPts, cosPts,
		 sinPts2, cosPts2,
		 wigners, scratch ) ;

      /***************************/
      /*                         */
      /* {f_{m1,0}} coefficient */
      /*                         */
      /***************************/


      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( m1, 0, bw ) ;
      coefHere = coefLoc_so3( m1, 0, bw ) ;

      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveAnalysis_fftw( m1, 0, bw, dataPtr,
			     wigners, weights,
			     coeffsPtr,
			     workspace_cx ) ;
      
      /***************************/
      /*                         */
      /* {f_{-m1,0}} coefficient */
      /*                         */
      /***************************/


      if ( flag == 0 ) /* data is complex */
	{
	        
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( -m1, 0, bw ) ;
	  coefHere = coefLoc_so3( -m1, 0, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
      
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveAnalysis_fftwX( -m1, 0, bw, dataPtr,
				  wigners, weights,
				  coeffsPtr,
				  workspace_cx ) ;
	}
      else  /* data is real, so use symmetry */
	{
	  coefHere = coefLoc_so3( m1, 0, bw );
	  coefHere2 = coefLoc_so3( -m1, 0, bw );

	  if ( (m1 % 2) == 0 )
	    fudge = 1.0 ;
	  else
	    fudge = -1.0 ;

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      coeffs[coefHere2+j][0] = fudge * coeffs[coefHere+j][0];
	      coeffs[coefHere2+j][1] = -fudge * coeffs[coefHere+j][1];
	    }
	  
	}

      /***************************/
      /*                         */
      /* {f_{0,m1}} coefficient */
      /*                         */
      /***************************/

      
      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( 0, m1, bw ) ;
      coefHere = coefLoc_so3( 0, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveAnalysis_fftwX( 0, m1, bw, dataPtr,
			      wigners, weights,
			      coeffsPtr,
			      workspace_cx ) ;
      

      /***************************/
      /*                         */
      /* {f_{0,-m1}} coefficient */
      /*                         */
      /***************************/


      if ( flag == 0 ) /* data is complex */
	{
      
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( 0, -m1, bw ) ;
	  coefHere = coefLoc_so3( 0, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
      
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveAnalysis_fftw( 0, -m1, bw, dataPtr,
				 wigners, weights,
				 coeffsPtr,
				 workspace_cx ) ;
	}
      else  /* data is real, so use symmetry */
	{
   	  coefHere = coefLoc_so3( 0, m1, bw );
	  coefHere2 = coefLoc_so3( 0, -m1, bw );

	  if ( (m1 % 2) == 0 )
	    fudge = 1.0 ;
	  else
	    fudge = -1.0 ;

	  for ( j = 0 ; j < bw - m1 ; j ++ )
	    {
	      coeffs[coefHere2+j][0] = fudge * coeffs[coefHere+j][0];
	      coeffs[coefHere2+j][1] = -fudge * coeffs[coefHere+j][1];
	    }
	  
	}
    }


  /***
      1 <= m1 <= bw-1
      m1+1 <= m2 <= bw-1
  ***/

  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      for ( m2 = m1 + 1 ; m2 < bw ; m2 ++ )
	{

	  
	  /* compute the wigners I'll need */
	  genWig_L2( m1, m2, bw,
		     sinPts, cosPts,
		     sinPts2, cosPts2,
		     wigners, scratch ) ;


	  /***************************/
	  /*                         */
	  /* {f_{m1,m2}} coefficient */
	  /*                         */
	  /***************************/

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, m2, bw ) ;
	  coefHere = coefLoc_so3( m1, m2, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_fftw( m1, m2, bw, dataPtr,
				 wigners, weights,
				 coeffsPtr,
				 workspace_cx ) ;
	  
	  /*****************************/
	  /*                           */
	  /* {f_{-m1,-m2}} coefficient */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* data is complex */
	    {

	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m1, -m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, -m2, bw ) ;

	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ;
	  
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis_fftwX( -m1, -m2, bw, dataPtr,
				      wigners, weights,
				      coeffsPtr,
				      workspace_cx ) ;
	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m1, m2, bw );
	      coefHere2 = coefLoc_so3( -m1, -m2, bw );
	  
	      
	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;

	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  coeffs[coefHere2+j][0] = fudge * coeffs[coefHere+j][0];
		  coeffs[coefHere2+j][1] = -fudge * coeffs[coefHere+j][1];
		}
	     
	    }


	  /****************************/
	  /*                          */
	  /* {f_{m1,-m2}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, -m2, bw ) ;
	  coefHere = coefLoc_so3( m1, -m2, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_fftwY( m1, -m2, bw, dataPtr,
				  wigners, weights,
				  coeffsPtr,
				  workspace_cx ) ;

	  /*****************************/
	  /*                           */
	  /* {f_{-m1,m2}} coefficient  */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* data is complex */
	    {

	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m1, m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, m2, bw ) ;

	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ;
	  
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis_fftwY( -m1, m2, bw, dataPtr,
				      wigners, weights,
				      coeffsPtr,
				      workspace_cx ) ;

	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m1, -m2, bw );
	      coefHere2 = coefLoc_so3( -m1, m2, bw );

	      
	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;
	      
	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  coeffs[coefHere2+j][0] = fudge * coeffs[coefHere+j][0];
		  coeffs[coefHere2+j][1] = -fudge * coeffs[coefHere+j][1];
		}
	      
	    }


	  
	  /***************************/
	  /*                         */
	  /* {f_{m2,m1}} coefficient */
	  /*                         */
	  /***************************/
	  
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, m1, bw ) ;
	  coefHere = coefLoc_so3( m2, m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_fftwX( m2, m1, bw, dataPtr,
				  wigners, weights,
				  coeffsPtr,
				  workspace_cx ) ;


	  /*****************************/
	  /*                           */
	  /* {f_{-m2,-m1}} coefficient */
	  /*                           */
	  /*****************************/
	  
	  if ( flag == 0 ) /* data is complex */
	    {


	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, -m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, -m1, bw ) ;

	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ;
	  
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis_fftw( -m2, -m1, bw, dataPtr,
				     wigners, weights,
				     coeffsPtr,
				     workspace_cx ) ;

	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m2, m1, bw );
	      coefHere2 = coefLoc_so3( -m2, -m1, bw );

	      
	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;

	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  coeffs[coefHere2+j][0] = fudge * coeffs[coefHere+j][0];
		  coeffs[coefHere2+j][1] = -fudge * coeffs[coefHere+j][1];
		}
	      
	    }


	  /****************************/
	  /*                          */
	  /* {f_{m2,-m1}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, -m1, bw ) ;
	  coefHere = coefLoc_so3( m2, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveAnalysis_fftwY( m1, -m2, bw, dataPtr,
				  wigners, weights,
				  coeffsPtr,
				  workspace_cx ) ;


	  /****************************/
	  /*                          */
	  /* {f_{-m2,m1}} coefficient */
	  /*                          */
	  /****************************/
  

	  if ( flag == 0 ) /* data is complex */
	    {
	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, m1, bw ) ;

	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ;
	  
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveAnalysis_fftwY( -m1, m2, bw, dataPtr,
				      wigners, weights,
				      coeffsPtr,
				      workspace_cx ) ;

	    }
	  else  /* data is real, so use symmetry */
	    {
	      coefHere = coefLoc_so3( m2, -m1, bw );
	      coefHere2 = coefLoc_so3( -m2, m1, bw );
	      

	      if ( ((m2-m1) % 2) == 0 )
		fudge = 1.0 ;
	      else
		fudge = -1.0 ;
	      
	      for ( j = 0 ; j < bw - m2 ; j ++ )
		{
		  coeffs[coefHere2+j][0] = fudge * coeffs[coefHere+j][0];
		  coeffs[coefHere2+j][1] = -fudge * coeffs[coefHere+j][1];
		}

	    }
	}
    }

  	  
  /* reset coef ptrs */
  coeffsPtr = coeffs ;

  /* need to normalize, one last time */
  dn = (M_PI /  ( (double) (bw * n )) );
  tmpInt = totalCoeffs_so3( bw ) ;
  for ( j = 0 ; j < tmpInt ; j ++ )
    {
      coeffsPtr[ j ][0] *= dn ;
      coeffsPtr[ j ][1] *= dn ;
    }

  /*** and we're done ! ***/
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

/*
  ok, the inverse transform

  Function arguments are as follows:

  bw = bandwidth of transform
  coeffs: plain COMPLEX array of size (4*bw^3-bw)/3, will contain the
          coefficients of the signal
  data: FFTW_COMPLEX array of size (2*bw)^3 containing the (output)
        signal ->
	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
	FFTW_COMPLEX!!! (although I will sometimes treat it
	as the plain COMPLEX array it is)

  workspace_cx: scratch space FFTW_COMPLEX array of size (2*bw)^3
	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
	FFTW_COMPLEX!!! (although I will sometimes treat it
	as the plain COMPLEX array it is)
  NOTE: if the fftw_plan is for an IN_PLACE transform, then
        workspace_cx needs to be only of size (2*bw)

  workspace_re: double scratch space of size 12*n + n*bw
		where n = 2*bw

  p1: pointer to FFTW plan for correctly ffting DATA array and
      placing the result in DATA (so it's an in-place FFT);
      this is a FORWARD FFT you're doing!!!

  flag: = 0 : data is COMPLEX
        = 1 : data is double

*/

void Inverse_SO3_Naive_fftw_nt( int bw,
				fftw_complex *coeffs,
				fftw_complex *data,
				fftw_complex *workspace_cx,
				double *workspace_re,
				fftw_plan *p1,
				int flag )
{
  int j, n ;
  int m1, m2 ;
  int sampHere , coefHere ;
  int sampHere2 ;
  fftw_complex *coeffsPtr, *dataPtr ;
  double *sinPts, *cosPts, *sinPts2, *cosPts2 ;
  double *wignersTrans, *scratch ;
  double dn ;

  n = 2 * bw ;

  sinPts = workspace_re ;
  cosPts = sinPts + n ;
  sinPts2 = cosPts + n ;
  cosPts2 = sinPts2 + n ;
  wignersTrans = cosPts2 + n ;
  scratch = wignersTrans + ( bw * n ) ; /* wignersTrans need at most bw*n
					   space AT ANY given orders m1, m2 */

  /*
    before going further, let's precompute all the sines
    and cosines I'll need. No matter what order transform
    I'm doing, these'll stay the same.
  */
  SinEvalPts( n, sinPts );
  CosEvalPts( n, cosPts );
  SinEvalPts2( n, sinPts2 );
  CosEvalPts2( n, cosPts2 );


  /*
    Stage 1: Do the Inverse Wigner transform. The rcoeffs, icoeffs
    arrays are assumed to be in the same "arrangement" as that produced
    by Forward_SO3_Naive().

    Since I'm working with two order indeces, m1 and m2, the
    for-loops will be more intricate than in the case of the
    "ordinary" spherical transform on S^2.

    Also, I will be taking advantage of the symmetries of the
    Wigner-d functions. As long as I keep my signs and flips
    right, the Wigner-d's I precompute for an order (m1, m2)
    transform can generally  be used in seven more transforms:
    (m1,-m2), (m2,m1), (m2,-m1), (-m2,m1), (-m2,-m1), (-m1,m2)
    and (-m1,-m2).


    The for-loops will be "designed" as follows. They will be
    divided into cases according to the orders:

    0) {f_{0,0}} inverse transform

    1) for 0 <= m1 <= bw-1
    compute inverse transform of
    i)   {f_{ m1, m1}}
    ii)  {f_{-m1,-m1}}
    iii) {f_{-m1, m1}}
    iv)  {f_{ m1,-m1}}

    2) for 1 <= m1 <= bw-1
    compute inverse transform of
    i)   {f_{ m1,  0}}
    ii)  {f_{-m1,  0}}
    iii) {f_{  0, m1}}
    iv)  {f_{  0,-m1}}

    3) for 1 <= m1 <= bw-1
    for m1+1 <= m2 <= bw-1
    compute inverse transform 
    i)    {f_{ m1, m2}}
    ii)   {f_{-m1,-m2}}
    iii)  {f_{ m1,-m2}}
    iv)   {f_{-m1, m2}}
    v)    {f_{ m2, m1}}
    vi)   {f_{-m2,-m1}}
    vii)  {f_{ m2,-m1}}
    viii) {f_{-m2, m1}}

    If assumptions are made regarding the original input signal,
    e.g. it's strictly real, then one may take advantage of
    symmetries of the big D wigners (i.e. function of all 3
    parameters alpha, beta, gamma) and so simplify the for-loops
    some and hence increase the speed of the program. However,
    the for-loops to follow will make no such assumptions.
    Whether the signal is real or complex, these for-loops will
    handle it.


    Fasten your seatbelt, folks. It's going to be a bumpy ride.

  */


  /* NOTE that I'm using the rdata, idata arrays as tmp space
     in the early going of the function */

  /***************************/
  /*                         */
  /* {f_{0,0}} coefficient   */
  /*                         */
  /***************************/
      
  /* compute the wigners I'll need */
  genWigTrans_L2( 0, 0, bw,
		  sinPts, cosPts,
		  sinPts2, cosPts2,
		  wignersTrans, scratch ) ;
  
  /* now, get the locations of where the
     samples I have to transform are, and
     where the coefficients have to go */
  
  sampHere = sampLoc_so3( 0, 0, bw ) ;
  coefHere = coefLoc_so3( 0, 0, bw ) ;
  
  /* ok, reset sample, coef ptrs */
  coeffsPtr = coeffs ;
  dataPtr = data ;
  
  /* now advance by the computed amounts */
  dataPtr += sampHere ;
  coeffsPtr += coefHere ;
  
  /* now transform the real and imaginary parts
     of the data */
  
  
  wigNaiveSynthesis_fftw( 0, 0, bw, coeffsPtr,
			  wignersTrans, dataPtr,
			  workspace_cx ) ;


  /*** 0 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      
      /* compute the wigners I'll need */
      genWigTrans_L2( m1, m1, bw,
		      sinPts, cosPts,
		      sinPts2, cosPts2,
		      wignersTrans, scratch ) ;
      
      /***************************/
      /*                         */
      /* {f_{m1,m1}} coefficient */
      /*                         */
      /***************************/


      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( m1, m1, bw ) ;
      coefHere = coefLoc_so3( m1, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ; ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveSynthesis_fftw( m1, m1, bw, coeffsPtr,
			      wignersTrans, dataPtr,
			      workspace_cx ) ;
     
      /*****************************/
      /*                           */
      /* {f_{-m1,-m1}} coefficient */
      /*                           */
      /*****************************/
      
      if ( flag == 0 ) /* if data is complex */
	{
	  
	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( -m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( -m1, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveSynthesis_fftw( -m1, -m1, bw, coeffsPtr,
				  wignersTrans, dataPtr,
				  workspace_cx ) ;
	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( m1, m1, bw );
	  sampHere2 = sampLoc_so3( -m1, -m1, bw );
	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      data[sampHere2+j][0] = data[sampHere+j][0];
	      data[sampHere2+j][1] = -data[sampHere+j][1];
	    }

	}


      /*****************************/
      /*                           */
      /* {f_{-m1,m1}} coefficient  */
      /*                           */
      /*****************************/


      /* now, get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( -m1, m1, bw ) ;
      coefHere = coefLoc_so3( -m1, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ; ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveSynthesis_fftwY( -m1, m1, bw, coeffsPtr,
			       wignersTrans, dataPtr,
			       workspace_cx ) ;

      /*****************************/
      /*                           */
      /* {f_{m1,-m1}} coefficient  */
      /*                           */
      /*****************************/

      if ( flag == 0 )  /* if data is complex */
	{
	  /* now, get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, -m1, bw ) ;
	  coefHere = coefLoc_so3( m1, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveSynthesis_fftwY( m1, -m1, bw, coeffsPtr,
				   wignersTrans, dataPtr,
				   workspace_cx ) ;
	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( -m1, m1, bw );
	  sampHere2 = sampLoc_so3( m1, -m1, bw );
	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      data[sampHere2+j][0] = data[sampHere+j][0];
	      data[sampHere2+j][1] = -data[sampHere+j][1];
	    }
	  
	}

    }

  /*** for 1 <= m1 <= bw-1 ***/
  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      
      /* compute the wigners I'll need */
      genWigTrans_L2( m1, 0, bw,
		      sinPts, cosPts,
		      sinPts2, cosPts2,
		      wignersTrans, scratch ) ;
      

      /***************************/
      /*                         */
      /* {f_{m1,0}} coefficient */
      /*                         */
      /***************************/


      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( m1, 0, bw ) ;
      coefHere = coefLoc_so3( m1, 0, bw ) ;

      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ; ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */

      wigNaiveSynthesis_fftw( m1, 0, bw, coeffsPtr,
			      wignersTrans, dataPtr,
			      workspace_cx ) ;

      /***************************/
      /*                         */
      /* {f_{-m1,0}} coefficient */
      /*                         */
      /***************************/

      if ( flag == 0 ) /* if data is complex */
	{
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */

	  sampHere = sampLoc_so3( -m1, 0, bw ) ;
	  coefHere = coefLoc_so3( -m1, 0, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
      
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;

	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_fftwX( -m1, 0, bw, coeffsPtr,
				   wignersTrans, dataPtr,
				   workspace_cx ) ;

	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( m1, 0, bw );
	  sampHere2 = sampLoc_so3( -m1, 0, bw );

	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      data[sampHere2+j][0] = data[sampHere+j][0];
	      data[sampHere2+j][1] = -data[sampHere+j][1];
	    }

	}


      /***************************/
      /*                         */
      /* {f_{0,m1}} coefficient */
      /*                         */
      /***************************/

      
      /* get the locations of where the
	 samples I have to transform are, and
	 where the coefficients have to go */

      sampHere = sampLoc_so3( 0, m1, bw ) ;
      coefHere = coefLoc_so3( 0, m1, bw ) ;

      /* ok, reset sample, coef ptrs */
      coeffsPtr = coeffs ;
      dataPtr = data ; ;
      
      /* now advance by the computed amounts */
      dataPtr += sampHere ;
      coeffsPtr += coefHere ;

      /* now transform the real and imaginary parts
	 of the data */
      
      wigNaiveSynthesis_fftwX( 0, m1, bw, coeffsPtr,
			       wignersTrans, dataPtr,
			       workspace_cx ) ;

      /***************************/
      /*                         */
      /* {f_{0,-m1}} coefficient */
      /*                         */
      /***************************/

      if ( flag == 0 ) /* if data is complex */
	{
	  
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( 0, -m1, bw ) ;
	  coefHere = coefLoc_so3( 0, -m1, bw ) ;
	  
	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveSynthesis_fftw( 0, -m1, bw, coeffsPtr,
				  wignersTrans, dataPtr,
				  workspace_cx ) ;
	}
      else  /* otherwise, use symmetry */
	{
	  sampHere = sampLoc_so3( 0, m1, bw );
	  sampHere2 = sampLoc_so3( 0, -m1, bw );



	  for ( j = 0 ; j < 2*bw ; j ++ )
	    {
	      data[sampHere2+j][0] = data[sampHere+j][0];
	      data[sampHere2+j][1] = -data[sampHere+j][1];
	    }

	}
    }


  /***
      1 <= m1 <= bw-1
      m1+1 <= m2 <= bw-1
  ***/

  for ( m1 = 1 ; m1 < bw ; m1 ++ )
    {
      for ( m2 = m1 + 1 ; m2 < bw ; m2 ++ )
	{

	  
	  /* compute the wigners I'll need */
	  genWigTrans_L2( m1, m2, bw,
			  sinPts, cosPts,
			  sinPts2, cosPts2,
			  wignersTrans, scratch ) ;


	  /***************************/
	  /*                         */
	  /* {f_{m1,m2}} coefficient */
	  /*                         */
	  /***************************/

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, m2, bw ) ;
	  coefHere = coefLoc_so3( m1, m2, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */
	  
	  wigNaiveSynthesis_fftw( m1, m2, bw, coeffsPtr,
				  wignersTrans, dataPtr,
				  workspace_cx ) ;

	  /*****************************/
	  /*                           */
	  /* {f_{-m1,-m2}} coefficient */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* if data is complex */
	    {
	      
	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	      
	      sampHere = sampLoc_so3( -m1, -m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, -m2, bw ) ;
	      
	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ; ;
	      
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	      
	      /* now transform the real and imaginary parts
		 of the data */
	      
	      wigNaiveSynthesis_fftwX( -m1, -m2, bw, coeffsPtr,
				       wignersTrans, dataPtr,
				       workspace_cx ) ;
	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m1, m2, bw );
	      sampHere2 = sampLoc_so3( -m1, -m2, bw );


	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  data[sampHere2+j][0] = data[sampHere+j][0];
		  data[sampHere2+j][1] = -data[sampHere+j][1];
		}

	    }
	  

	  /****************************/
	  /*                          */
	  /* {f_{m1,-m2}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m1, -m2, bw ) ;
	  coefHere = coefLoc_so3( m1, -m2, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_fftwY( m1, -m2, bw, coeffsPtr,
				   wignersTrans, dataPtr,
				   workspace_cx ) ;

	  /*****************************/
	  /*                           */
	  /* {f_{-m1,m2}} coefficient  */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* if data is complex */
	    {
	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	      
	      sampHere = sampLoc_so3( -m1, m2, bw ) ;
	      coefHere = coefLoc_so3( -m1, m2, bw ) ;
	      
	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ; ;
	      
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	      
	      /* now transform the real and imaginary parts
		 of the data */
	      
	      wigNaiveSynthesis_fftwY( -m1, m2, bw, coeffsPtr,
				       wignersTrans, dataPtr,
				       workspace_cx ) ;
	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m1, -m2, bw );
	      sampHere2 = sampLoc_so3( -m1, m2, bw );

	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  data[sampHere2+j][0] = data[sampHere+j][0];
		  data[sampHere2+j][1] = -data[sampHere+j][1];
		}

	    }


	  /***************************/
	  /*                         */
	  /* {f_{m2,m1}} coefficient */
	  /*                         */
	  /***************************/
	  
	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, m1, bw ) ;
	  coefHere = coefLoc_so3( m2, m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_fftwX( m2, m1, bw, coeffsPtr,
				   wignersTrans, dataPtr,
				   workspace_cx ) ;


	  /*****************************/
	  /*                           */
	  /* {f_{-m2,-m1}} coefficient */
	  /*                           */
	  /*****************************/

	  if ( flag == 0 ) /* if data is complex */
	    {
	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, -m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, -m1, bw ) ;

	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ; ;
	  
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */
	  
	      wigNaiveSynthesis_fftw( -m2, -m1, bw, coeffsPtr,
				      wignersTrans, dataPtr,
				      workspace_cx ) ;
	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m2, m1, bw );
	      sampHere2 = sampLoc_so3( -m2, -m1, bw );

	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  data[sampHere2+j][0] = data[sampHere+j][0];
		  data[sampHere2+j][1] = -data[sampHere+j][1];
		}

	    }

	  
	  /****************************/
	  /*                          */
	  /* {f_{m2,-m1}} coefficient */
	  /*                          */
	  /****************************/
  

	  /* get the locations of where the
	     samples I have to transform are, and
	     where the coefficients have to go */
	  
	  sampHere = sampLoc_so3( m2, -m1, bw ) ;
	  coefHere = coefLoc_so3( m2, -m1, bw ) ;

	  /* ok, reset sample, coef ptrs */
	  coeffsPtr = coeffs ;
	  dataPtr = data ; ;
	  
	  /* now advance by the computed amounts */
	  dataPtr += sampHere ;
	  coeffsPtr += coefHere ;
	  
	  /* now transform the real and imaginary parts
	     of the data */

	  wigNaiveSynthesis_fftwY( m1, -m2, bw, coeffsPtr,
				   wignersTrans, dataPtr,
				   workspace_cx ) ;
	  

	  /****************************/
	  /*                          */
	  /* {f_{-m2,m1}} coefficient */
	  /*                          */
	  /****************************/
  
	  if ( flag == 0 ) /* if data is complex */
	    {
	      /* get the locations of where the
		 samples I have to transform are, and
		 where the coefficients have to go */
	  
	      sampHere = sampLoc_so3( -m2, m1, bw ) ;
	      coefHere = coefLoc_so3( -m2, m1, bw ) ;

	      /* ok, reset sample, coef ptrs */
	      coeffsPtr = coeffs ;
	      dataPtr = data ; ;
	  
	      /* now advance by the computed amounts */
	      dataPtr += sampHere ;
	      coeffsPtr += coefHere ;
	  
	      /* now transform the real and imaginary parts
		 of the data */

	      wigNaiveSynthesis_fftwY( -m1, m2, bw, coeffsPtr,
				       wignersTrans, dataPtr,
				       workspace_cx ) ;
	    }
	  else  /* otherwise, use symmetry */
	    {
	      sampHere = sampLoc_so3( m2, -m1, bw );
	      sampHere2 = sampLoc_so3( -m2, m1, bw );

	      for ( j = 0 ; j < 2*bw ; j ++ )
		{
		  data[sampHere2+j][0] = data[sampHere+j][0];
		  data[sampHere2+j][1] = -data[sampHere+j][1];
		}

	    }
	 
	}
    }


  /* I need to set some zeros in the intermediate t2r, t2i arrays,
     so that I can take the fft correctly */


  /* reset ptrs to correct starting positions */
  dataPtr = data + (n)*(bw) ;

  for ( m1 = 0 ; m1 < bw  ; m1 ++ )
    {
      memset( dataPtr, 0, sizeof(fftw_complex) * n );
      dataPtr += (2*n)*(bw) ;
    }

  dataPtr = data + bw*n*(n);
  memset( dataPtr, 0, sizeof(fftw_complex) * n * n );
  dataPtr += n * n + n*bw;

  for ( m1 = 1 ; m1 < bw  ; m1 ++ )
    {
      memset( dataPtr, 0, sizeof(fftw_complex) * n );
      dataPtr += (2*n)*(bw) ;
    }


  /*
    Stage 2: FFT 
  */

  fftw_execute( *p1 ) ;


  /*
    IF the fftw plan you make for the inverse SO(3) transform
    is IN-PLACE, then you can comment out the following line.
    I'm leaving it in to make it better correspond to the
    forward SO(3) transform defined in this file - I can' seem
    to make an in-place fftw_plan for the forward SO(3) transform.
  */
  memcpy( data, workspace_cx, sizeof(fftw_complex) * n * n * n) ;

  /* normalize the Fourier coefficients (sorry, have to do it) */
  dn = 1./( (double) n ); 
  dn *= ( ((double) bw) / M_PI ) ;
   
  for ( j = 0 ; j < n*n*n; j++ )
    {
      data[ j ][0] *= dn ;
      data[ j ][1] *= dn ;
    }

  /* and that's all, folks */
}
