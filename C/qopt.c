#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <time.h>

#include <omp.h>

#define MAX(a,b) ((a)>(b) ? (a) : (b))

int read_clock(struct timeval *t)
{
    return gettimeofday(t, NULL);
}

int elapsed_time(struct timeval *start, struct timeval *end)
{
    return (int)(end->tv_sec  - start->tv_sec) * 1e6 + 
		   (int)(end->tv_usec - start->tv_usec);
}

double qopt_ref( double bx, double by, double bz,
		         int gx, int gy, int gz,
				 int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
//  double zprd = prd[2];
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;
  int nx_pppm_6 = gx,
	  ny_pppm_6 = gy,
	  nz_pppm_6 = gz;

  int order_6 = P;
  double g_ewald_6 = ewald;

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1.0/inv2ew;
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);      // nz_pppm6 -> grid size in z dir

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + // k^2
          pow(unitkz*mper,2.0);

		// unitk is the vector component 2Pi/L for k's (will be multiplied by "n" for the actual ks)
		// kper is the n?
		// q_ is the actual vector k+2Pi/h m
		// (nx,ny,nz) is the m
		// w_ is tilde(W)^(P)
        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) { // -2 to +2
            qx = unitkx*(kper+nx_pppm_6*nx);
            sx = exp(-qx*qx*inv2ew*inv2ew);
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm_6;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
            for (ny = -nby; ny <= nby; ny++) {
              qy = unitky*(lper+ny_pppm_6*ny);
              sy = exp(-qy*qy*inv2ew*inv2ew);
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm_6;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
              for (nz = -nbz; nz <= nbz; nz++) {
                qz = unitkz*(mper+nz_pppm_6*nz);
                sz = exp(-qz*qz*inv2ew*inv2ew);
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm_6;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew*inv2ew)*sx*sy*sz + // first half
                       2*dot2*rtdot2*inv2ew*inv2ew*inv2ew*rtpi*erfc(rtdot2*inv2ew); // second half or tilde(R)(k)  (p.137 of rolf)
                term *= g_ewald_6*g_ewald_6*g_ewald_6;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += term*term*M_PI*M_PI*M_PI/9.0 * dot2;
                sum2 += -u2*term*M_PI*rtpi/3.0*dot1;
                sum3 += u2;
              }
            }
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
		  /*printf("QOPT: %.15e\n", qopt);*/
        }
      }
    }
  }
  return qopt;
}

// Precomputes erfc
double qopt_v1( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
//  double zprd = prd[2];
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;
  int nx_pppm_6 = gx,
	  ny_pppm_6 = gy,
	  nz_pppm_6 = gz;

  int order_6 = P;
  double g_ewald_6 = ewald;

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1.0/inv2ew;
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int vx, vy, vz;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr= 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  /*unitky*20*unitky*20 +*/
  /*unitkz*64*unitkz*64)*inv2ew));*/
  /*dot2 = qx*qx+qy*qy+qz*qz;*/
  /*rtdot2 = sqrt(dot2);*/
  /*erfc(rtdot2*inv2ew); // second half or tilde(R)(k)  (p.137 of rolf)*/

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);      // nz_pppm6 -> grid size in z dir

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + // k^2
          pow(unitkz*mper,2.0);

		// unitk is the vector component 2Pi/L for k's (will be multiplied by "n" for the actual ks)
		// kper is the n?
		// q_ is the actual vector k+2Pi/h m
		// (nx,ny,nz) is the m
		// w_ is tilde(W)^(P)
        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) { // -2 to +2
			vx = kper+nx_pppm_6*nx;
			vx *= vx;
            qx = unitkx*(kper+nx_pppm_6*nx);
            sx = exp(-qx*qx*inv2ew*inv2ew);
            wx = 1.0;
            argx = 0.5*qx*xprd/nx_pppm_6;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
            for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+ny_pppm_6*ny;
			  vy *= vy;
              qy = unitky*(lper+ny_pppm_6*ny);
              sy = exp(-qy*qy*inv2ew*inv2ew);
              wy = 1.0;
              argy = 0.5*qy*yprd/ny_pppm_6;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
              for (nz = -nbz; nz <= nbz; nz++) {
			    vz = mper+nz_pppm_6*nz;
			    vz *= vz;
                qz = unitkz*(mper+nz_pppm_6*nz);
                sz = exp(-qz*qz*inv2ew*inv2ew);
                wz = 1.0;
                argz = 0.5*qz*zprd_slab/nz_pppm_6;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);

                dot1 = unitkx*kper*qx + unitky*lper*qy + unitkz*mper*qz;
                dot2 = qx*qx+qy*qy+qz*qz;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew*inv2ew)*sx*sy*sz + // first half
                       2*dot2*rtdot2*inv2ew*inv2ew*inv2ew*rtpi*erfc_prec[vx+vy+vz]; // second half or tilde(R)(k)  (p.137 of rolf)
                term *= g_ewald_6*g_ewald_6*g_ewald_6;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += term*term*M_PI*M_PI*M_PI/9.0 * dot2;
                sum2 += -u2*term*M_PI*rtpi/3.0*dot1;
                sum3 += u2;
              }
            }
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }

  free( erfc_prec );

  return qopt;
}

// Precomputes some subexpressions
double qopt_v2( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
//  double zprd = prd[2];
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;
  int nx_pppm_6 = gx,
	  ny_pppm_6 = gy,
	  nz_pppm_6 = gz;

  int order_6 = P;
  double g_ewald_6 = ewald;

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1.0/inv2ew;
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int vx, vy, vz;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr= 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  //
  double inv2ew2 = inv2ew * inv2ew;
  double inv2ew3 = inv2ew * inv2ew2;
  double ewald3 = ewald*ewald*ewald;
  double pi3o9 = M_PI * M_PI * M_PI / 9.0;
  double pi3ho3 = M_PI * rtpi / 3.0;
  double hx = bx/gx, hy = by/gy, hz = bz/gz, hzh = hz*0.5;
  double kx, ky, kz;
  double qx2, qy2, qz2;
  double wxwy, wxyz, sxsy;

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);      // nz_pppm6 -> grid size in z dir
	kz = mper*unitkz;

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);
	  ky = lper*unitky;

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
	    kx = kper*unitkx;
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + // k^2
          pow(unitkz*mper,2.0);

		// unitk is the vector component 2Pi/L for k's (will be multiplied by "n" for the actual ks)
		// kper is the n?
		// q_ is the actual vector k+2Pi/h m
		// (nx,ny,nz) is the m
		// w_ is tilde(W)^(P)
        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) { // -2 to +2
			vx = kper+nx_pppm_6*nx;
			vx *= vx;
            qx = unitkx*(kper+nx_pppm_6*nx);
			qx2 = qx*qx;
            sx = exp(-qx2*inv2ew2);
            wx = 1.0;
            argx = 0.5*qx*hx;
            if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
            for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+ny_pppm_6*ny;
			  vy *= vy;
              qy = unitky*(lper+ny_pppm_6*ny);
			  qy2 = qy*qy;
              sy = exp(-qy2*inv2ew2);
			  sxsy = sx*sy;
              wy = 1.0;
              argy = 0.5*qy*hy;
              if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
			  wxwy = wx*wy;
			  //
              for (nz = -nbz; nz <= nbz; nz++) {
			    vz = mper+nz_pppm_6*nz;
			    vz *= vz;
                qz = unitkz*(mper+nz_pppm_6*nz);
				qz2 = qz*qz;
                sz = exp(-qz2*inv2ew2);
                wz = 1.0;
                argz = qz*hzh;
                if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);

                dot1 = kx*qx + ky*qy + kz*qz;
                dot2 = qx2+qy2+qz2;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew2)*sxsy*sz + // first half
                       2*dot2*rtdot2*inv2ew3*rtpi*erfc_prec[vx+vy+vz]; // second half or tilde(R)(k)  (p.137 of rolf)
                term *= ewald3;
				wxyz = wxwy * wz;
                u2 = wxyz * wxyz;
                sum1 += term * term * pi3o9 * dot2;
                sum2 += -u2 * term * pi3ho3 * dot1;
                sum3 += u2;
              }
            }
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }
  
  free( erfc_prec );

  return qopt;
}

// sin and pow gone
double qopt_v3( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;
  int nx_pppm_6 = gx,
	  ny_pppm_6 = gy,
	  nz_pppm_6 = gz;

  int order_6 = P;
  double g_ewald_6 = ewald;

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1.0/inv2ew;
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int vx, vy, vz;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr = 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  //
  double inv2ew2 = inv2ew * inv2ew;
  double inv2ew3 = inv2ew * inv2ew2;
  double ewald3 = ewald*ewald*ewald;
  double pi3o9 = M_PI * M_PI * M_PI / 9.0;
  double pi3ho3 = M_PI * rtpi / 3.0;
  double hx = bx/gx, hy = by/gy, hz = bz/gz, hzh = hz*0.5;
  double kx, ky, kz;
  double qx2, qy2, qz2;
  double wxwy, wxyz, sxsy;
  // sin + pow
  double *sinx_prec = (double *) malloc ( (5*gx+1) * sizeof(double) );
  double *siny_prec = (double *) malloc ( (5*gy+1) * sizeof(double) );
  double *sinz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  int i, gshift,
	  gx_shift = gx/2+2*gx,
	  gy_shift = gy/2+2*gy,
	  gz_shift = gz/2+2*gz;
  gshift = gx/2+2*gx;
  for ( i = -((gx/2)+2*gx); i < 0; i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinx_prec[gshift] = 1.0;
  for ( i = 1; i < (gx/2+2*gx); i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gy/2+2*gy;
  for ( i = -((gy/2)+2*gy); i < 0; i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  siny_prec[gshift] = 1.0;
  for ( i = 1; i < (gy/2+2*gy); i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gz/2+2*gz;
  for ( i = -((gz/2)+2*gz); i < 0; i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinz_prec[gshift] = 1.0;
  for ( i = 1; i < (gz/2+2*gz); i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);      // nz_pppm6 -> grid size in z dir
	kz = mper*unitkz;

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);
	  ky = lper*unitky;

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
	    kx = kper*unitkx;
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + // k^2
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) { // -2 to +2
			vx = kper+nx_pppm_6*nx;
            qx = unitkx*(kper+nx_pppm_6*nx);
			qx2 = qx*qx;
            sx = exp(-qx2*inv2ew2);
            //wx = 1.0;
            //argx = 0.5*qx*hx;
            //if (argx != 0.0) wx = pow(sin(argx)/argx,order_6);
			wx = sinx_prec[vx+gx_shift];
			//wx = pow( sinx_prec[vx+gx_shift], order_6 );
			vx *= vx;
            for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+ny_pppm_6*ny;
              qy = unitky*(lper+ny_pppm_6*ny);
			  qy2 = qy*qy;
              sy = exp(-qy2*inv2ew2);
			  sxsy = sx*sy;
              //wy = 1.0;
              //argy = 0.5*qy*hy;
              //if (argy != 0.0) wy = pow(sin(argy)/argy,order_6);
			  wy = siny_prec[vy+gy_shift];
			  //wy = pow( siny_prec[vy+gy_shift], order_6 );
			  vy *= vy;
			  wxwy = wx*wy;
			  // 
              for (nz = -nbz; nz <= nbz; nz++) {
			    vz = mper+nz_pppm_6*nz;
                qz = unitkz*(mper+nz_pppm_6*nz);
				qz2 = qz*qz;
                sz = exp(-qz2*inv2ew2);
				//wz = 1.0;
                //argz = qz*hzh;
				//if (argz != 0.0) wz = pow(sin(argz)/argz,order_6);
				wz = sinz_prec[vz+gz_shift];
				//wz = pow( sinz_prec[vz+gz_shift], order_6 );
			    vz *= vz;

                dot1 = kx*qx + ky*qy + kz*qz;
                dot2 = qx2+qy2+qz2;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew2)*sxsy*sz + // first half
                       2*dot2*rtdot2*inv2ew3*rtpi*erfc_prec[vx+vy+vz]; // second half or tilde(R)(k)  (p.137 of rolf)
                term *= ewald3;
				wxyz = wxwy * wz;
                u2 = wxyz * wxyz;
                sum1 += term * term * pi3o9 * dot2;
                sum2 += -u2 * term * pi3ho3 * dot1;
                sum3 += u2;
              }
            }
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }

  free( erfc_prec );
  free( sinx_prec );
  free( siny_prec );
  free( sinz_prec );

  return qopt;
}

// exp gone
double qopt_v4( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
//  double zprd = prd[2];
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;
  int nx_pppm_6 = gx,
	  ny_pppm_6 = gy,
	  nz_pppm_6 = gz;

  double g_ewald_6 = ewald;

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1.0/inv2ew;
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int vx, vy, vz;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr = 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  //
  double inv2ew2 = inv2ew * inv2ew;
  double inv2ew3 = inv2ew * inv2ew2;
  double ewald3 = ewald*ewald*ewald;
  double pi3o9 = M_PI * M_PI * M_PI / 9.0;
  double pi3ho3 = M_PI * rtpi / 3.0;
//  double hx = bx/gx, hy = by/gy, hz = bz/gz, hzh = hz*0.5;
  double kx, ky, kz;
  double qx2, qy2, qz2;
  double wxwy, wxyz, sxsy;
  // sin + pow
  double *sinx_prec = (double *) malloc ( (5*gx+1) * sizeof(double) );
  double *siny_prec = (double *) malloc ( (5*gy+1) * sizeof(double) );
  double *sinz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  int i, gshift,
	  gx_shift = gx/2+2*gx,
	  gy_shift = gy/2+2*gy,
	  gz_shift = gz/2+2*gz;
  gshift = gx/2+2*gx;
  for ( i = -((gx/2)+2*gx); i < 0; i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinx_prec[gshift] = 1.0;
  for ( i = 1; i < (gx/2+2*gx); i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gy/2+2*gy;
  for ( i = -((gy/2)+2*gy); i < 0; i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  siny_prec[gshift] = 1.0;
  for ( i = 1; i < (gy/2+2*gy); i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gz/2+2*gz;
  for ( i = -((gz/2)+2*gz); i < 0; i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinz_prec[gshift] = 1.0;
  for ( i = 1; i < (gz/2+2*gz); i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  // exp
  double *expx_prec = (double *) malloc ( (max_x+1) * sizeof(double) );
  double *expy_prec = (double *) malloc ( (max_y+1) * sizeof(double) );
  double *expz_prec = (double *) malloc ( (max_z+1) * sizeof(double) );
  double temp;
  for ( i = 0; i <= max_x; i++ ) {
	  temp = i*unitkx;
	  arg = -temp*temp*inv2ew2;
	  expx_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_y; i++ ) {
	  temp = i*unitky;
	  arg = -temp*temp*inv2ew2;
	  expy_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_z; i++ ) {
	  temp = i*unitkz;
	  arg = -temp*temp*inv2ew2;
	  expz_prec[i] = exp(arg);
                //sz = exp(-qz2*inv2ew2);
  }

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);      // nz_pppm6 -> grid size in z dir
	kz = mper*unitkz;

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);
	  ky = lper*unitky;

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
	    kx = kper*unitkx;
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + // k^2
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) { // -2 to +2
			vx = kper+nx_pppm_6*nx;
            qx = unitkx*(vx);
			qx2 = qx*qx;
            sx = expx_prec[vx < 0 ? -vx : vx];
			wx = sinx_prec[vx+gx_shift];
			vx *= vx;
            for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+ny_pppm_6*ny;
              qy = unitky*(vy);
			  qy2 = qy*qy;
              sy = expy_prec[vy < 0 ? -vy : vy];
			  sxsy = sx*sy;
			  wy = siny_prec[vy+gy_shift];
			  vy *= vy;
			  wxwy = wx*wy;
			  // 
              for (nz = -nbz; nz <= nbz; nz++) {
			    vz = mper+nz_pppm_6*nz;
                qz = unitkz*(vz);
				qz2 = qz*qz;
                sz = expz_prec[vz < 0 ? -vz : vz];
				wz = sinz_prec[vz+gz_shift];
			    vz *= vz;

                dot1 = kx*qx + ky*qy + kz*qz;
                dot2 = qx2+qy2+qz2;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew2)*sxsy*sz + // first half
                       2*dot2*rtdot2*inv2ew3*rtpi*erfc_prec[vx+vy+vz]; // second half or tilde(R)(k)  (p.137 of rolf)
                term *= ewald3;
				wxyz = wxwy * wz;
                u2 = wxyz * wxyz;
                sum1 += term * term * pi3o9 * dot2;
                sum2 += -u2 * term * pi3ho3 * dot1;
                sum3 += u2;
              }
            }
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }

  free( erfc_prec );
  free( sinx_prec );
  free( siny_prec );
  free( sinz_prec );
  free( expx_prec );
  free( expy_prec );
  free( expz_prec );

  return qopt;
}

// more minor + iterate more naturally
double qopt_v5( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 1/(2*ewald);
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int vx, vy, vz;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr = 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  //
  double inv2ew2 = inv2ew * inv2ew;
  double inv2ew3 = inv2ew * inv2ew2;
  double ewald3 = ewald*ewald*ewald;
  double pi3o9 = M_PI * M_PI * M_PI / 9.0;
  double pi3ho3 = M_PI * rtpi / 3.0;
  double kx, ky, kz;
  double qx2, qy2, qz2;
  double wxwy, wxyz, sxsy;
  // sin + pow
  double *sinx_prec = (double *) malloc ( (5*gx+1) * sizeof(double) );
  double *siny_prec = (double *) malloc ( (5*gy+1) * sizeof(double) );
  double *sinz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  int i, gshift,
	  gx_shift = gx/2+2*gx,
	  gy_shift = gy/2+2*gy,
	  gz_shift = gz/2+2*gz;
  gshift = gx/2+2*gx;
  for ( i = -((gx/2)+2*gx); i < 0; i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinx_prec[gshift] = 1.0;
  for ( i = 1; i < (gx/2+2*gx); i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gy/2+2*gy;
  for ( i = -((gy/2)+2*gy); i < 0; i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  siny_prec[gshift] = 1.0;
  for ( i = 1; i < (gy/2+2*gy); i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gz/2+2*gz;
  for ( i = -((gz/2)+2*gz); i < 0; i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinz_prec[gshift] = 1.0;
  for ( i = 1; i < (gz/2+2*gz); i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  // exp
  double *expx_prec = (double *) malloc ( (max_x+1) * sizeof(double) );
  double *expy_prec = (double *) malloc ( (max_y+1) * sizeof(double) );
  double *expz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  double temp;
  for ( i = 0; i <= max_x; i++ ) {
	  temp = i*unitkx;
	  arg = -temp*temp*inv2ew2;
	  expx_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_y; i++ ) {
	  temp = i*unitky;
	  arg = -temp*temp*inv2ew2;
	  expy_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_z; i++ ) {
	  temp = i*unitkz;
	  arg = -temp*temp*inv2ew2;
	  expz_prec[i] = exp(arg);
	  //expz_prec[i+gz_shift] = exp(arg);
	  //expz_prec[-i+gz_shift] = expz_prec[i+gz_shift];
  }
  // more minor
  double dot1_x, dot1_y;
  double inv2ew3_rtpi_2 = 2 * inv2ew3 * rtpi;
  int v_shift;
  int vz2;

  //for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
  //  mper = m - gz*(2*m/gz);
  for (mper = -(gz/2); mper < gz/2; mper++) {
	kz = mper*unitkz;

    //for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
    //  lper = l - gy*(2*l/gy);
    for (lper = -(gy/2); lper < gy/2; lper++) {
	  ky = lper*unitky;

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - gx*(2*k/gx);
      //for (kper = -(gx/2); kper < gx/2; kper++) {
	    kx = kper*unitkx;
      
        sqk = kx*kx + ky*ky + kz*kz;
        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) { // -2 to +2
			vx = kper+gx*nx;
            qx = unitkx*(vx);
			qx2 = qx*qx;
            sx = expx_prec[vx < 0 ? -vx : vx];
			wx = sinx_prec[vx+gx_shift];
			vx *= vx;
			dot1_x = kx*qx;
            for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+gy*ny;
              qy = unitky*(vy);
			  qy2 = qy*qy;
              sy = expy_prec[vy < 0 ? -vy : vy];
			  sxsy = sx*sy;
			  wy = siny_prec[vy+gy_shift];
			  vy *= vy;
			  wxwy = wx*wy;
			  dot1_y = ky*qy;
			  //
              for (nz = -nbz; nz <= nbz; nz++) {
			    vz = mper+gz*nz;
				v_shift = vz+gz_shift;
                qz = unitkz*(vz);
				qz2 = qz*qz;
                //sz = expz_prec[v_shift];
                sz = expz_prec[vz < 0 ? -vz : vz];
				wz = sinz_prec[v_shift];
			    vz *= vz;

                dot1 = dot1_x + dot1_y + kz*qz;
                dot2 = qx2+qy2+qz2;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew2)*sxsy*sz + // first half
                       dot2*rtdot2*inv2ew3_rtpi_2*erfc_prec[vx+vy+vz]; // second half or tilde(R)(k)  (p.137 of rolf)
                term *= ewald3;
				wxyz = wxwy * wz;
                u2 = wxyz * wxyz;
                sum1 += term * term * pi3o9 * dot2;
                sum2 += -u2 * term * pi3ho3 * dot1;
                sum3 += u2;
              }
            }
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }

  free( erfc_prec );
  free( sinx_prec );
  free( siny_prec );
  free( sinz_prec );
  free( expx_prec );
  free( expy_prec );
  free( expz_prec );

  return qopt;
}

// cache
double qopt_v6( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 1/(2*ewald);
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int vx, vy, vz;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr = 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  //
  double inv2ew2 = inv2ew * inv2ew;
  double inv2ew3 = inv2ew * inv2ew2;
  double ewald3 = ewald*ewald*ewald;
  double pi3o9 = M_PI * M_PI * M_PI / 9.0;
  double pi3ho3 = M_PI * rtpi / 3.0;
  double kx, ky, kz;
  double wxwy, wxyz, sxsy;
  // sin + pow
  double *sinx_prec = (double *) malloc ( (5*gx+1) * sizeof(double) );
  double *siny_prec = (double *) malloc ( (5*gy+1) * sizeof(double) );
  double *sinz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  int i, gshift,
	  gx_shift = gx/2+2*gx,
	  gy_shift = gy/2+2*gy,
	  gz_shift = gz/2+2*gz;
  gshift = gx/2+2*gx;
  for ( i = -((gx/2)+2*gx); i < 0; i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinx_prec[gshift] = 1.0;
  for ( i = 1; i < (gx/2+2*gx); i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gy/2+2*gy;
  for ( i = -((gy/2)+2*gy); i < 0; i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  siny_prec[gshift] = 1.0;
  for ( i = 1; i < (gy/2+2*gy); i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gz/2+2*gz;
  for ( i = -((gz/2)+2*gz); i < 0; i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinz_prec[gshift] = 1.0;
  for ( i = 1; i < (gz/2+2*gz); i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  // exp
  double *expx_prec = (double *) malloc ( (max_x+1) * sizeof(double) );
  double *expy_prec = (double *) malloc ( (max_y+1) * sizeof(double) );
  double *expz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  double temp;
  for ( i = 0; i <= max_x; i++ ) {
	  temp = i*unitkx;
	  arg = -temp*temp*inv2ew2;
	  expx_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_y; i++ ) {
	  temp = i*unitky;
	  arg = -temp*temp*inv2ew2;
	  expy_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_z; i++ ) {
	  temp = i*unitkz;
	  arg = -temp*temp*inv2ew2;
	  expz_prec[i] = exp(arg);
	  //expz_prec[i+gz_shift] = exp(arg);
	  //expz_prec[-i+gz_shift] = expz_prec[i+gz_shift];
  }
  // more minor
  double inv2ew3_rtpi_2 = 2 * inv2ew3 * rtpi;
  // cache
  int vxs[10], vys[10], vzs[10];
  double qxs[25], qys[25], qzs[25];
  int *vxs_p, *vys_p, *vzs_p;
  double *qxs_p, *qys_p, *qzs_p;

  //for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
  //  mper = m - gz*(2*m/gz);
  for (mper = -(gz/2); mper < gz/2; mper++) {
	kz = mper*unitkz;

		  int zv_idx = 0, qz_idx = 0;
          for (nz = -nbz; nz <= nbz; nz++) {
			  vz = mper+gz*nz;
			  vzs[zv_idx++] = vz;    // vz
			  vzs[zv_idx++] = vz*vz; // vz2
			  //
			  qz = unitkz*(vz);
			  qzs[qz_idx++] = qz;    // qz
			  qzs[qz_idx++] = qz*qz; // qz2
			  qzs[qz_idx++] = expz_prec[vz < 0 ? -vz : vz]; // sz
			  qzs[qz_idx++] = sinz_prec[vz+gz_shift];       // wz
			  qzs[qz_idx++] = kz*qz; // dot1_z
		  }
    //for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
    //  lper = l - gy*(2*l/gy);
    for (lper = -(gy/2); lper < gy/2; lper++) {
	  ky = lper*unitky;
		  int yv_idx = 0, qy_idx = 0;
          for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+gy*ny;
			  vys[yv_idx++] = vy;    // vy
			  vys[yv_idx++] = vy*vy; // vy2
			  //
			  qy = unitky*(vy);
			  qys[qy_idx++] = qy;    // qy
			  qys[qy_idx++] = qy*qy; // qy2
			  qys[qy_idx++] = expy_prec[vy < 0 ? -vy : vy]; // sy
			  qys[qy_idx++] = siny_prec[vy+gy_shift];       // wy
			  qys[qy_idx++] = ky*qy; // dot1_y
		  }

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - gx*(2*k/gx);
      //for (kper = -(gx/2); kper < gx/2; kper++) {
	    kx = kper*unitkx;
      
        sqk = kx*kx + ky*ky + kz*kz;
        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
		  //
		  int xv_idx = 0, qx_idx = 0;
          for (nx = -nbx; nx <= nbx; nx++) {
			  vx = kper+gx*nx;
			  vxs[xv_idx++] = vx;    // vx
			  vxs[xv_idx++] = vx*vx; // vx2
			  //
			  qx = unitkx*(vx);
			  qxs[qx_idx++] = qx;    // qx
			  qxs[qx_idx++] = qx*qx; // qx2
			  qxs[qx_idx++] = expx_prec[vx < 0 ? -vx : vx]; // sx
			  qxs[qx_idx++] = sinx_prec[vx+gx_shift];       // wx
			  qxs[qx_idx++] = kx*qx; // dot1_x
		  }

		  vxs_p = &vxs[0];
		  qxs_p = &qxs[0];
          for (nx = 0; nx <= 2*nbx; nx++) { // 0 to 4
		    vys_p = &vys[0];
		    qys_p = &qys[0];
            for (ny = 0; ny <= 2*nby; ny++) {
			  sxsy = qxs_p[2]*qys_p[2];
			  wxwy = qxs_p[3]*qys_p[3];
		      vzs_p = &vzs[0];
		      qzs_p = &qzs[0];
              for (nz = 0; nz <= 2*nbz; nz++) {
                dot1 = qxs_p[4] + qys_p[4] + qzs_p[4];
                dot2 = qxs_p[1] + qys_p[1] + qzs_p[1];
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew2)*sxsy*qzs_p[2] +
                       dot2*rtdot2*inv2ew3_rtpi_2*erfc(rtdot2*inv2ew);
                term *= ewald3;
				wxyz = wxwy * qzs_p[3];
                u2 = wxyz * wxyz;
                sum1 += term * term * pi3o9 * dot2;
                sum2 += -u2 * term * pi3ho3 * dot1;
                sum3 += u2;
				/*lper + (ny-2)*gy,*/
				/*mper + (nz-2)*gz,*/
				/*term * term * pi3o9 * dot2);*/
				vzs_p += 2;
				qzs_p += 5;
              }
			  vys_p += 2;
			  qys_p += 5;
            }
			vxs_p += 2;
			qxs_p += 5;
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }

  free( erfc_prec );
  free( sinx_prec );
  free( siny_prec );
  free( sinz_prec );
  free( expx_prec );
  free( expy_prec );
  free( expz_prec );

  return qopt;
}

// full erfc grid
double qopt_v7( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int i, j, k;

  double unitkx = (2.0*M_PI/bx);
  double unitky = (2.0*M_PI/by);
  double unitkz = (2.0*M_PI/bz);

  int nx,ny,nz,kper,lper,mper;
  double sqk, u2;
  double qx,qy,qz;
  double sum1,sum2, sum3;
  double dot1,dot2, rtdot2, term;
  double inv2ew = 1/(2*ewald);
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int vx, vy, vz;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr = 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  //
  double inv2ew2 = inv2ew * inv2ew;
  double inv2ew3 = inv2ew * inv2ew2;
  double ewald3 = ewald*ewald*ewald;
  double pi3o9 = M_PI * M_PI * M_PI / 9.0;
  double pi3ho3 = M_PI * rtpi / 3.0;
  double kx, ky, kz;
  double wxwy, wxyz, sxsy;
  // sin + pow
  double *sinx_prec = (double *) malloc ( (5*gx+1) * sizeof(double) );
  double *siny_prec = (double *) malloc ( (5*gy+1) * sizeof(double) );
  double *sinz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  int gshift,
	  gx_shift = gx/2+2*gx,
	  gy_shift = gy/2+2*gy,
	  gz_shift = gz/2+2*gz;
  gshift = gx/2+2*gx;
  for ( i = -((gx/2)+2*gx); i < 0; i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinx_prec[gshift] = 1.0;
  for ( i = 1; i < (gx/2+2*gx); i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gy/2+2*gy;
  for ( i = -((gy/2)+2*gy); i < 0; i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  siny_prec[gshift] = 1.0;
  for ( i = 1; i < (gy/2+2*gy); i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gz/2+2*gz;
  for ( i = -((gz/2)+2*gz); i < 0; i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinz_prec[gshift] = 1.0;
  for ( i = 1; i < (gz/2+2*gz); i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  // exp
  double *expx_prec = (double *) malloc ( (max_x+1) * sizeof(double) );
  double *expy_prec = (double *) malloc ( (max_y+1) * sizeof(double) );
  double *expz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  double temp;
  for ( i = 0; i <= max_x; i++ ) {
	  temp = i*unitkx;
	  arg = -temp*temp*inv2ew2;
	  expx_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_y; i++ ) {
	  temp = i*unitky;
	  arg = -temp*temp*inv2ew2;
	  expy_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_z; i++ ) {
	  temp = i*unitkz;
	  arg = -temp*temp*inv2ew2;
	  expz_prec[i] = exp(arg);
	  //expz_prec[i+gz_shift] = exp(arg);
	  //expz_prec[-i+gz_shift] = expz_prec[i+gz_shift];
  }
  // more minor
  double inv2ew3_rtpi_2 = 2 * inv2ew3 * rtpi;
  // cache
  int vxs[10], vys[10], vzs[10];
  double qxs[25], qys[25], qzs[25];
  int *vxs_p, *vys_p, *vzs_p;
  double *qxs_p, *qys_p, *qzs_p;
  // full erfc
  long erfc_grid_size = 125*gx*gy*gz;
  double *erfc_grid = (double *) malloc ( erfc_grid_size * sizeof(double) );
  double *erfc_grid_p = erfc_grid;
  int ii,jj,kk;
  for (k = -(gz/2); k < gz/2; k++) {
    for (j = -(gy/2); j < gy/2; j++) {
      for (i = -(gx/2); i < gx/2; i++) {
		if (i || j || k) {
        for ( ii = -nbx; ii <= nbx; ii++ ) {
          vx = i + ii * gx;
          vx *= vx;
	      for ( jj = -nby; jj <= nby; jj++ ) {
	        vy = j + jj * gy;
	        vy *= vy;
            for ( kk = -nbz; kk <= nbz; kk++ ) {
              vz = k + kk * gz;
              vz *= vz;
              arg_idx = vx + vy + vz;
              (*erfc_grid_p) = erfc_prec[arg_idx];
			  erfc_grid_p++;
			}
		  }
		}
		}
	  }
	}
  }
  erfc_grid_p = erfc_grid;


  //for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
  //  mper = m - gz*(2*m/gz);
  for (mper = -(gz/2); mper < gz/2; mper++) {
	kz = mper*unitkz;

		  int zv_idx = 0, qz_idx = 0;
          for (nz = -nbz; nz <= nbz; nz++) {
			  vz = mper+gz*nz;
			  vzs[zv_idx++] = vz;    // vz
			  vzs[zv_idx++] = vz*vz; // vz2
			  //
			  qz = unitkz*(vz);
			  qzs[qz_idx++] = qz;    // qz
			  qzs[qz_idx++] = qz*qz; // qz2
			  qzs[qz_idx++] = expz_prec[vz < 0 ? -vz : vz]; // sz
			  qzs[qz_idx++] = sinz_prec[vz+gz_shift];       // wz
			  qzs[qz_idx++] = kz*qz; // dot1_z
		  }
    //for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
    //  lper = l - gy*(2*l/gy);
    for (lper = -(gy/2); lper < gy/2; lper++) {
	  ky = lper*unitky;
		  int yv_idx = 0, qy_idx = 0;
          for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+gy*ny;
			  vys[yv_idx++] = vy;    // vy
			  vys[yv_idx++] = vy*vy; // vy2
			  //
			  qy = unitky*(vy);
			  qys[qy_idx++] = qy;    // qy
			  qys[qy_idx++] = qy*qy; // qy2
			  qys[qy_idx++] = expy_prec[vy < 0 ? -vy : vy]; // sy
			  qys[qy_idx++] = siny_prec[vy+gy_shift];       // wy
			  qys[qy_idx++] = ky*qy; // dot1_y
		  }

      //for (k = 0; k <= gx-1; k++) {
      //  kper = k - gx*(2*k/gx);
      for (kper = -(gx/2); kper < gx/2; kper++) {
	    kx = kper*unitkx;
      
        sqk = kx*kx + ky*ky + kz*kz;
        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
		  //
		  int xv_idx = 0, qx_idx = 0;
          for (nx = -nbx; nx <= nbx; nx++) {
			  vx = kper+gx*nx;
			  vxs[xv_idx++] = vx;    // vx
			  vxs[xv_idx++] = vx*vx; // vx2
			  //
			  qx = unitkx*(vx);
			  qxs[qx_idx++] = qx;    // qx
			  qxs[qx_idx++] = qx*qx; // qx2
			  qxs[qx_idx++] = expx_prec[vx < 0 ? -vx : vx]; // sx
			  qxs[qx_idx++] = sinx_prec[vx+gx_shift];       // wx
			  qxs[qx_idx++] = kx*qx; // dot1_x
		  }

		  vxs_p = &vxs[0];
		  qxs_p = &qxs[0];
          for (nx = 0; nx <= 2*nbx; nx++) { // 0 to 4
		    vys_p = &vys[0];
		    qys_p = &qys[0];
            for (ny = 0; ny <= 2*nby; ny++) {
			  sxsy = qxs_p[2]*qys_p[2];
			  wxwy = qxs_p[3]*qys_p[3];
		      vzs_p = &vzs[0];
		      qzs_p = &qzs[0];
              for (nz = 0; nz <= 2*nbz; nz++) {
                dot1 = qxs_p[4] + qys_p[4] + qzs_p[4];
                dot2 = qxs_p[1] + qys_p[1] + qzs_p[1];
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew2)*sxsy*qzs_p[2] +
                       dot2*rtdot2*inv2ew3_rtpi_2*(*erfc_grid_p);
				erfc_grid_p++;
				term *= ewald3;
				wxyz = wxwy * qzs_p[3];
                u2 = wxyz * wxyz;
                sum1 += term * term * pi3o9 * dot2;
                sum2 += -u2 * term * pi3ho3 * dot1;
                sum3 += u2;
				vzs_p += 2;
				qzs_p += 5;
              }
			  vys_p += 2;
			  qys_p += 5;
            }
			vxs_p += 2;
			qxs_p += 5;
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }

  free( erfc_prec );
  free( sinx_prec );
  free( siny_prec );
  free( sinz_prec );
  free( expx_prec );
  free( expy_prec );
  free( expz_prec );
  free( erfc_grid );

  return qopt;
}

// pragma
double qopt_v8( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;

  double xprd = bx;
  double yprd = by;
  double zprd_slab = bz;

  double unitkx = (2.0*M_PI/xprd);
  double unitky = (2.0*M_PI/yprd);
  double unitkz = (2.0*M_PI/zprd_slab);

  int nxlo_fft_6 = 0,
	  nxhi_fft_6 = gx-1,
      nylo_fft_6 = 0,
	  nyhi_fft_6 = gy-1,
      nzlo_fft_6 = 0,
	  nzhi_fft_6 = gz-1;

  int nx,ny,nz,kper,lper,mper;
  double inv2ew = 1/(2*ewald);
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  // Precompute erfc
  int max_x = gx/2 + 2*gx,
      max_y = gy/2 + 2*gy,
      max_z = gz/2 + 2*gz;
  int ctr, 
	  ctr_x, ctr_y, ctr_z;
  int arg_idx;
  double arg;
  double *erfc_prec = (double *) malloc ( (max_x*max_x + max_y*max_y + max_z*max_z + 1) * sizeof(double) );

  ctr = 0;
  for ( ctr_x = 0; ctr_x <= max_x; ctr_x++ )
	  for ( ctr_y = ctr_x; ctr_y <= max_y; ctr_y++ )
		  for ( ctr_z = ctr_y; ctr_z <= max_z; ctr_z++ )
		  {
			  arg_idx = ctr_x*ctr_x + ctr_y*ctr_y + ctr_z*ctr_z;
			  arg = unitkx*unitkx*ctr_x*ctr_x +
				    unitky*unitky*ctr_y*ctr_y +
				    unitkz*unitkz*ctr_z*ctr_z;
			  erfc_prec[arg_idx] = erfc( sqrt(arg)*inv2ew );
			  ctr++;
		  }
  //
  double inv2ew2 = inv2ew * inv2ew;
  double inv2ew3 = inv2ew * inv2ew2;
  double ewald3 = ewald*ewald*ewald;
  double pi3o9 = M_PI * M_PI * M_PI / 9.0;
  double pi3ho3 = M_PI * rtpi / 3.0;
  // sin + pow
  double *sinx_prec = (double *) malloc ( (5*gx+1) * sizeof(double) );
  double *siny_prec = (double *) malloc ( (5*gy+1) * sizeof(double) );
  double *sinz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  int i, gshift,
	  gx_shift = gx/2+2*gx,
	  gy_shift = gy/2+2*gy,
	  gz_shift = gz/2+2*gz;
  gshift = gx/2+2*gx;
  for ( i = -((gx/2)+2*gx); i < 0; i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinx_prec[gshift] = 1.0;
  for ( i = 1; i < (gx/2+2*gx); i++ ) {
	  arg = i * M_PI / gx;
	  //sinx_prec[i+gshift] = sin( arg ) / arg;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gy/2+2*gy;
  for ( i = -((gy/2)+2*gy); i < 0; i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  siny_prec[gshift] = 1.0;
  for ( i = 1; i < (gy/2+2*gy); i++ ) {
	  arg = i * M_PI / gy;
	  //siny_prec[i+gshift] = sin( arg ) / arg;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gz/2+2*gz;
  for ( i = -((gz/2)+2*gz); i < 0; i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinz_prec[gshift] = 1.0;
  for ( i = 1; i < (gz/2+2*gz); i++ ) {
	  arg = i * M_PI / gz;
	  sinz_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  // exp
  double *expx_prec = (double *) malloc ( (max_x+1) * sizeof(double) );
  double *expy_prec = (double *) malloc ( (max_y+1) * sizeof(double) );
  double *expz_prec = (double *) malloc ( (5*gz+1) * sizeof(double) );
  double temp;
  for ( i = 0; i <= max_x; i++ ) {
	  temp = i*unitkx;
	  arg = -temp*temp*inv2ew2;
	  expx_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_y; i++ ) {
	  temp = i*unitky;
	  arg = -temp*temp*inv2ew2;
	  expy_prec[i] = exp(arg);
  }
  for ( i = 0; i <= max_z; i++ ) {
	  temp = i*unitkz;
	  arg = -temp*temp*inv2ew2;
	  expz_prec[i] = exp(arg);
	  //expz_prec[i+gz_shift] = exp(arg);
	  //expz_prec[-i+gz_shift] = expz_prec[i+gz_shift];
  }
  // more minor
  double inv2ew3_rtpi_2 = 2 * inv2ew3 * rtpi;
  // cache
  int *vxs_p, *vys_p, *vzs_p;
  double *qxs_p, *qys_p, *qzs_p;

  //for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) { // nzlo/hi -> local grid points in z dir
  //  mper = m - gz*(2*m/gz);
  omp_set_num_threads(2);
  #pragma omp parallel for private(mper,lper,kper,nz,ny,nx) reduction(+:qopt) schedule(static)
  for (mper = -(gz/2); mper < gz/2; mper++) {
	  printf("%d\n",omp_get_thread_num());
		  int vx, vy, vz;
		  int vxs[10], vys[10], vzs[10];
		  double qxs[25], qys[25], qzs[25];
		  double kx, ky, kz;
		  double wxwy, wxyz, sxsy;
		  double sqk, u2;
		  double qx,qy,qz;
		  double sum1,sum2, sum3;
		  double dot1,dot2, rtdot2, term;
	kz = mper*unitkz;

		  int zv_idx = 0, qz_idx = 0;
          for (nz = -nbz; nz <= nbz; nz++) {
			  vz = mper+gz*nz;
			  vzs[zv_idx++] = vz;    // vz
			  vzs[zv_idx++] = vz*vz; // vz2
			  //
			  qz = unitkz*(vz);
			  qzs[qz_idx++] = qz;    // qz
			  qzs[qz_idx++] = qz*qz; // qz2
			  qzs[qz_idx++] = expz_prec[vz < 0 ? -vz : vz]; // sz
			  qzs[qz_idx++] = sinz_prec[vz+gz_shift];       // wz
			  qzs[qz_idx++] = kz*qz; // dot1_z
		  }
    //for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
    //  lper = l - gy*(2*l/gy);
    for (lper = -(gy/2); lper < gy/2; lper++) {
	  ky = lper*unitky;
		  int yv_idx = 0, qy_idx = 0;
          for (ny = -nby; ny <= nby; ny++) {
			  vy = lper+gy*ny;
			  vys[yv_idx++] = vy;    // vy
			  vys[yv_idx++] = vy*vy; // vy2
			  //
			  qy = unitky*(vy);
			  qys[qy_idx++] = qy;    // qy
			  qys[qy_idx++] = qy*qy; // qy2
			  qys[qy_idx++] = expy_prec[vy < 0 ? -vy : vy]; // sy
			  qys[qy_idx++] = siny_prec[vy+gy_shift];       // wy
			  qys[qy_idx++] = ky*qy; // dot1_y
		  }

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - gx*(2*k/gx);
      //for (kper = -(gx/2); kper < gx/2; kper++) {
	    kx = kper*unitkx;
      
        sqk = kx*kx + ky*ky + kz*kz;
        if (sqk != 0.0) {
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
		  //
		  int xv_idx = 0, qx_idx = 0;
          for (nx = -nbx; nx <= nbx; nx++) {
			  vx = kper+gx*nx;
			  vxs[xv_idx++] = vx;    // vx
			  vxs[xv_idx++] = vx*vx; // vx2
			  //
			  qx = unitkx*(vx);
			  qxs[qx_idx++] = qx;    // qx
			  qxs[qx_idx++] = qx*qx; // qx2
			  qxs[qx_idx++] = expx_prec[vx < 0 ? -vx : vx]; // sx
			  qxs[qx_idx++] = sinx_prec[vx+gx_shift];       // wx
			  qxs[qx_idx++] = kx*qx; // dot1_x
		  }

		  vxs_p = &vxs[0];
		  qxs_p = &qxs[0];
          for (nx = 0; nx <= 2*nbx; nx++) { // 0 to 4
		    vys_p = &vys[0];
		    qys_p = &qys[0];
            for (ny = 0; ny <= 2*nby; ny++) {
			  sxsy = qxs_p[2]*qys_p[2];
			  wxwy = qxs_p[3]*qys_p[3];
		      vzs_p = &vzs[0];
		      qzs_p = &qzs[0];
              for (nz = 0; nz <= 2*nbz; nz++) {
                dot1 = qxs_p[4] + qys_p[4] + qzs_p[4];
                dot2 = qxs_p[1] + qys_p[1] + qzs_p[1];
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew2)*sxsy*qzs_p[2] +
                       dot2*rtdot2*inv2ew3_rtpi_2*erfc_prec[vxs_p[1]+vys_p[1]+vzs_p[1]];
                term *= ewald3;
				wxyz = wxwy * qzs_p[3];
                u2 = wxyz * wxyz;
                sum1 += term * term * pi3o9 * dot2;
                sum2 += -u2 * term * pi3ho3 * dot1;
                sum3 += u2;
				vzs_p += 2;
				qzs_p += 5;
              }
			  vys_p += 2;
			  qys_p += 5;
            }
			vxs_p += 2;
			qxs_p += 5;
          }
          sum2 *= sum2;
          sum3 *= sum3*sqk;
          qopt += sum1 -sum2/sum3;
        }
      }
    }
  }

  free( erfc_prec );
  free( sinx_prec );
  free( siny_prec );
  free( sinz_prec );
  free( expx_prec );
  free( expy_prec );
  free( expz_prec );

  return qopt;
}

double* kspace_error( FILE *fp )
{
	double bx, by, bz;
	int N;
	double accuracy;

	double ewald_min, ewald_max, d_ewald;
	int p_min, p_max;
	int np, newald, ng;

	/*fp = fopen( cfg_file, "r" );*/
	fscanf( fp, "%lf %lf %lf", &bx, &by, &bz );
	fscanf( fp, "%d", &N );
	fscanf( fp, "%lf", &accuracy );
	fscanf( fp, "%lf %lf %lf", &ewald_min, &ewald_max, &d_ewald );
	fscanf( fp, "%d %d", &p_min, &p_max );

	np = p_max - p_min + 1;
	newald = (ewald_max - ewald_min) / d_ewald + 1;
	fscanf( fp, "%d", &ng );

	int *grid = (int *) malloc (3 * ng * sizeof(int));
	int i;
	for ( i = 0; i < 3*ng; i+=3 )
		fscanf( fp, "%d %d %d", &grid[i], &grid[i+1], &grid[i+2] );

	double V = bx*by*bz;
	double C = 2.0*2.0*N;
	double q, err;

	int P, ig, gx,gy,gz;
	double ewald;

	double *err_grid = (double *) malloc (np*newald*ng*sizeof(double));
	int ierr = 0;
	for ( ierr = 0; ierr < np*newald*ng; ierr++ )
		err_grid[ierr] = 1.0;

	P = p_min;
	while (P <= p_max)
	{
		ig = 0;
		int last_grids_ewald_e = 0;
		while (ig < ng*3)
		{
			gx = grid[ig];
			gy = grid[ig+1];
			gz = grid[ig+2];

			int ewald_o = 0, ewald_e = newald;

			q = qopt_v6(bx,by,bz, gx,gy,gz, P, ewald_min);
			err = sqrt(q/N)*C/V;
/*			err_grid[(P-p_min) *   ng*newald +
				     ig/3      *   newald    +
					 ewald_o] = err;
*/			while (ewald_e - ewald_o > 1)
			{
				int mid_idx = (ewald_e + ewald_o) / 2;
				ewald = ewald_min +  mid_idx * d_ewald;
						
				q = qopt_v6(bx,by,bz, gx,gy,gz, P, ewald);
				err = sqrt(q/N)*C/V;
				err_grid[(P-p_min) *   ng*newald +
						 ig/3      *   newald    +
						 mid_idx] = err;

				if (err > accuracy) // not accurate enough, decrease max_ewald (_e)
					ewald_e = mid_idx;
				else
					ewald_o = mid_idx;

			}
			// from ewald_e (first inacc for this grid) until lasts inacc, 
			//   min_grid for these ewalds is current
			int idx;
			printf("Last grids ewald_e: %d\n", last_grids_ewald_e);
			printf("%d - %d\n", ewald_o, ewald_e);
			printf("%f\n", err_grid[(P-p_min) *   ng*newald +
							 ig/3      *   newald    +
							 ewald_o]);
			printf("%f\n", err_grid[(P-p_min) *   ng*newald +
							 ig/3      *   newald    +
							 ewald_e]);

			for ( idx = last_grids_ewald_e; idx < ewald_e; idx++ )
			{
				q = qopt_v6(bx,by,bz, gx,gy,gz, P, ewald_min+idx*d_ewald);
				err = sqrt(q/N)*C/V;
				err_grid[(P-p_min) *   ng*newald +
						 ig/3      *   newald    +
						 idx] 
					= err;
				/*
					err_grid[(P-p_min) *   ng*newald +
							 ig/3      *   newald    +
							 ewald_e];
				*/
			}
			for ( idx = 0; idx < last_grids_ewald_e; idx++ )
			{
				err_grid[(P-p_min) *   ng*newald +
						 ig/3      *   newald    +
						 idx] 
					= 
					err_grid[(P-p_min) *   ng*newald +
							 ig/3      *   newald    +
							 last_grids_ewald_e];
			}
			last_grids_ewald_e = ewald_e;
			ig += 3;
		}
		P++;
	}

	int j, k;
	ierr = 0;
	for ( i = 0; i < np; i++ )
	{
		for ( j = 0; j < ng; j++ )
			for ( k = 0; k < newald; k++ )
			{
				printf("%d (%3d, %3d, %3d) %.2f: %.15e\n", 
						i+p_min, grid[3*j], grid[3*j+1], grid[3*j+2], 
						ewald_min + k*d_ewald, 
						// grid[i][j][k]);
						err_grid[ierr]);
				ierr++;
			}
		fflush(stdout);
	}

	return err_grid;
}

int main(void)
{
	double bx, by, bz;
	int gx, gy, gz;
	int P;
	double ewald;

	int N = 2000;
	double V;
	double C = 2.0*2.0*N;
	double q, err;

	struct timeval tstart, tend;
	float time, time_ref;

	read_clock(&tstart);
	/*double *err_grid = kspace_error( "kspace_error-LargeCube.cfg" );*/
	double *err_grid = kspace_error( stdin );
	read_clock(&tend);
	time_ref = elapsed_time(&tstart, &tend)/1e6;
	printf(" %7.2f  ", elapsed_time(&tstart, &tend)/1e6);

	return 0;
}
