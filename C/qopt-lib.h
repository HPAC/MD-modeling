#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct 
{
	double domx;
	double domy;
	double domz;
	int npart;
	double accuracy;

	double ewald_min;
	double ewald_max;
	double d_ewald;
	int p_min;
	int p_max;
	int ngrids;
	int *grids;
} config;

// Load and free space configuration
void parse_config_file( char *cfg_file, config *cfg );
void free_config( config *cfg );

// Define as a type a pointer to the functions below
// for readibility
typedef double (*qopt_funcp)(double, double, double, int, int, int, int, double);
// ik error bounds
double qopt_ik_ref( double bx, double by, double bz,
		         int gx, int gy, int gz,
				 int P, double ewald );
// v6
double qopt_ik( double bx, double by, double bz,
		        int gx, int gy, int gz,
			    int P, double ewald );

// ad error bounds
double qopt_ad_ref( double bx, double by, double bz,
		         int gx, int gy, int gz,
				 int P, double ewald );
