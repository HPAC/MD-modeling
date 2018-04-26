#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "qopt-lib.h"

void parse_config_file( char *cfg_file, config *cfg )
{
	FILE *fp;
	int i;

	fp = fopen( cfg_file, "r" );
	fscanf( fp, "%lf %lf %lf", &cfg->domx, &cfg->domy, &cfg->domz );
	fscanf( fp, "%d", &cfg->npart );
	fscanf( fp, "%lf", &cfg->accuracy );
	fscanf( fp, "%lf %lf %lf", &cfg->ewald_min, &cfg->ewald_max, &cfg->d_ewald );
	fscanf( fp, "%d %d", &cfg->p_min, &cfg->p_max );

	fscanf( fp, "%d", &cfg->ngrids );

	cfg->grids = (int *) malloc (3 * cfg->ngrids * sizeof(int));
	for ( i = 0; i < 3*cfg->ngrids; i+=3 )
		fscanf( fp, "%d %d %d", &cfg->grids[i], &cfg->grids[i+1], &cfg->grids[i+2] );
}

void free_config( config *cfg )
{
	free( cfg->grids );
	cfg->grids = NULL;
}


double qopt_ik_ref( double bx, double by, double bz,
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

  int nbx = 1;
  int nby = 1;
  int nbz = 1;

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
        }
      }
    }
  }
  return qopt;
}

// v6
double qopt_ik( double bx, double by, double bz,
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
  int vx, vy, vz;
  int arg_idx;
  double arg;
  
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
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  sinx_prec[gshift] = 1.0;
  for ( i = 1; i < (gx/2+2*gx); i++ ) {
	  arg = i * M_PI / gx;
	  sinx_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  //
  gshift = gy/2+2*gy;
  for ( i = -((gy/2)+2*gy); i < 0; i++ ) {
	  arg = i * M_PI / gy;
	  siny_prec[i+gshift] = pow(sin( arg ) / arg, P);
  }
  siny_prec[gshift] = 1.0;
  for ( i = 1; i < (gy/2+2*gy); i++ ) {
	  arg = i * M_PI / gy;
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

  free( sinx_prec );
  free( siny_prec );
  free( sinz_prec );
  free( expx_prec );
  free( expy_prec );
  free( expz_prec );

  return qopt;
}

double qopt_ad_ref( double bx, double by, double bz,
		         int gx, int gy, int gz,
				 int P, double ewald )
{
  double qopt = 0.0;
  int k,l,m;
  double *prd;

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
  double argx,argy,argz,wx,wy,wz,sx,sy,sz,qx,qy,qz;
  double u2, sqk;
  double sum1,sum2,sum3,sum4;
  double dot2, rtdot2, term;
  double inv2ew = 2*g_ewald_6;
  inv2ew = 1/inv2ew;
  double rtpi = sqrt(M_PI);

  int nbx = 2;
  int nby = 2;
  int nbz = 2;

  for (m = nzlo_fft_6; m <= nzhi_fft_6; m++) {
    mper = m - nz_pppm_6*(2*m/nz_pppm_6);

    for (l = nylo_fft_6; l <= nyhi_fft_6; l++) {
      lper = l - ny_pppm_6*(2*l/ny_pppm_6);

      for (k = nxlo_fft_6; k <= nxhi_fft_6; k++) {
        kper = k - nx_pppm_6*(2*k/nx_pppm_6);
      
        sqk = pow(unitkx*kper,2.0) + pow(unitky*lper,2.0) + 
          pow(unitkz*mper,2.0);

        if (sqk != 0.0) {
    
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          for (nx = -nbx; nx <= nbx; nx++) {
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

                dot2 = qx*qx+qy*qy+qz*qz;
                rtdot2 = sqrt(dot2);
                term = (1-2*dot2*inv2ew*inv2ew)*sx*sy*sz +
		       2*dot2*rtdot2*inv2ew*inv2ew*inv2ew*rtpi*erfc(rtdot2*inv2ew);
                term *= g_ewald_6*g_ewald_6*g_ewald_6;
                u2 =  pow(wx*wy*wz,2.0);
                sum1 += term*term*M_PI*M_PI*M_PI/9.0 * dot2;
                sum2 += -term*M_PI*rtpi/3.0 * u2 * dot2;
                sum3 += u2;
                sum4 += dot2*u2;
              }
            }
          }
          sum2 *= sum2;
          qopt += sum1 - sum2/(sum3*sum4);
        }
      }
    }
  }
  return qopt;
}

