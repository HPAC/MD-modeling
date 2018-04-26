#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <time.h>

#include <mpi.h>

#include "qopt-lib.h"

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


double* kspace_error( config *cfg, qopt_funcp qopt )
{
	int Ps = cfg->p_max - cfg->p_min + 1,
		newald = (cfg->ewald_max - cfg->ewald_min) / cfg->d_ewald + 1,
		proc_id, nprocs, max_range, my_range_o, my_range_e;

	int N = cfg->npart;
	double V = cfg->domx * cfg->domy * cfg->domz;
	double C = 2.0 * 2.0 * N;
	double q, err;

	int i, P, ig, gx,gy,gz;
	double ewald;

	double *err_grid = (double *) malloc ( Ps * newald * cfg->ngrids * sizeof(double));
	int ierr = 0;
	for ( ierr = 0; ierr < Ps * newald * cfg->ngrids; ierr++ )
		err_grid[ierr] = 1.0;

	MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
	MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
	max_range = Ps * cfg->ngrids;
	my_range_o =  proc_id    * max_range / nprocs;
	my_range_e = (proc_id+1) * max_range / nprocs - 1;

	int *bounds = (int *) malloc ( 2 * max_range * sizeof(int) );
	int range_it;

	for ( range_it = my_range_o; range_it <= my_range_e; range_it++ )
	{
		P  = range_it / cfg->ngrids + cfg->p_min;
		ig = range_it % cfg->ngrids;

		gx = cfg->grids[3*ig];
		gy = cfg->grids[3*ig+1];
		gz = cfg->grids[3*ig+2];

		int ewald_o = 0, ewald_e = newald-1;

		// calculate error at min and max ewald
		// for a correct steering of the search
		q = qopt( cfg->domx, cfg->domy, cfg->domz, 
					 gx, gy, gz, P, ewald_o );
		err = sqrt(q/N) * C/V;
		err_grid[(P - cfg->p_min) * cfg->ngrids*newald +
				 ig               *             newald +
				 ewald_o] = err;
		q = qopt( cfg->domx, cfg->domy, cfg->domz, 
					 gx, gy, gz, P, ewald_e );
		err = sqrt(q/N) * C/V;
		err_grid[(P - cfg->p_min) * cfg->ngrids*newald +
				 ig               *             newald +
				 ewald_e] = err;
		while (ewald_e - ewald_o > 1)
		{
			int mid_idx = (ewald_e + ewald_o) / 2;
			ewald = cfg->ewald_min + mid_idx * cfg->d_ewald;
					
			q = qopt( cfg->domx, cfg->domy, cfg->domz, 
					     gx, gy, gz, P, ewald );
			err = sqrt(q/N) * C/V;
			err_grid[(P - cfg->p_min) * cfg->ngrids*newald +
					 ig               *             newald +
					 mid_idx] = err;

			if (err > cfg->accuracy) // not accurate enough, decrease max_ewald (_e)
				ewald_e = mid_idx;
			else
				ewald_o = mid_idx;

		}
		bounds[2 * range_it    ] = ewald_o;
		bounds[2 * range_it + 1] = ewald_e;
	}

	// Communicate
	int *recv_counts, *displs;
	int range_o, range_e, displ;
	recv_counts = (int *) malloc (nprocs * sizeof(int));
	displs      = (int *) malloc (nprocs * sizeof(int));
	displ = 0;
	for ( i = 0; i < nprocs; i++ )
	{
		range_o =  i    * max_range / nprocs;
		range_e = (i+1) * max_range / nprocs;
		recv_counts[i] = 2 * (range_e - range_o);
		displs[i] = displ;
		displ += recv_counts[i];
	}
	MPI_Allgatherv( &bounds[2*my_range_o], 2*(my_range_e - my_range_o + 1), MPI_INT,
			bounds, recv_counts, displs, MPI_INT, MPI_COMM_WORLD );
	
	// Compute interface
	for ( range_it = my_range_o; range_it <= my_range_e; range_it++ )
	{
		P  = range_it / cfg->ngrids + cfg->p_min;
		ig = range_it % cfg->ngrids;

		gx = cfg->grids[3*ig];
		gy = cfg->grids[3*ig+1];
		gz = cfg->grids[3*ig+2];

		// from last grid's ewald_e (first inacc) until current ewald_o (last acc),
		// compute the error
		int ewald_e, last_grids_ewald_e;
		if ( ig == 0 )
			last_grids_ewald_e = 0;
		else
			last_grids_ewald_e = bounds[2*(range_it-1)+1];
		ewald_e = bounds[2*range_it+1];
		int idx;
		for ( idx = last_grids_ewald_e; idx < ewald_e; idx++ )
		{
			q = qopt( cfg->domx, cfg->domy, cfg->domz, 
					     gx, gy, gz, 
						 P, cfg->ewald_min + idx * cfg->d_ewald);
			err = sqrt(q/N) * C/V;
			err_grid[(P - cfg->p_min) * cfg->ngrids*newald +
					 ig               *             newald +
					 idx] = err;
		}
		// Propagate to smaller ewalds
		for ( idx = 0; idx < last_grids_ewald_e; idx++ )
		{
			err_grid[(P - cfg->p_min) * cfg->ngrids * newald +
					 ig               *               newald +
					 idx] 
				= 
				err_grid[(P - cfg->p_min) * cfg->ngrids * newald +
						 ig               *               newald +
						 last_grids_ewald_e-1];
		}
	}

	// Gather
	recv_counts = (int *) malloc (nprocs * sizeof(int));
	displs      = (int *) malloc (nprocs * sizeof(int));
	displ = 0;
	for ( i = 0; i < nprocs; i++ )
	{
		range_o =  i    * max_range / nprocs;
		range_e = (i+1) * max_range / nprocs;
		recv_counts[i] = newald * (range_e - range_o);
		displs[i] = displ;
		displ += recv_counts[i];
	}
	MPI_Allgatherv( &err_grid[my_range_o*newald], newald*(my_range_e - my_range_o + 1), MPI_DOUBLE,
			err_grid, recv_counts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

	free(recv_counts);
	free(displs);
	
	return err_grid;
}

void print_error_grid( double *error_grid, config *cfg )
{
	int Ps, newald;
	int i, j, k, ierr;

	Ps = cfg->p_max - cfg->p_min + 1;
	newald = (cfg->ewald_max - cfg->ewald_min) / cfg->d_ewald + 1,

	ierr = 0;
	for ( i = 0; i < Ps; i++ )
		for ( j = 0; j < cfg->ngrids; j++ )
			for ( k = 0; k < newald; k++ )
			{
				printf("%d (%3d, %3d, %3d) %.2f: %.15e\n", 
						i + cfg->p_min, 
						cfg->grids[3*j], cfg->grids[3*j+1], cfg->grids[3*j+2], 
						cfg->ewald_min + k*cfg->d_ewald, 
						error_grid[ierr]);
				ierr++;
			}
}

int main( int argc, char *argv[] )
{
	struct timeval tstart, tend;
	float time, max_time;
	int proc_id;
	config cfg;

	qopt_funcp qopt;

	MPI_Init( &argc, &argv );

	MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );

	if ( argc != 3 && proc_id == 0 )
	{
		fprintf( stderr, "Usage: %s [ik|ad] ConfigFile\n", argv[0] );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	if (strncmp( argv[1], "ik", 255 ) == 0)
	{
		qopt = &qopt_ik;
	}
	else if (strncmp( argv[1], "ad", 255 ) == 0)
	{
		qopt = &qopt_ad_ref;
	}
	else 
	{
		fprintf( stderr, "[ERROR] Unknown diff method. Valid values: ik|ad\n" );
		MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
	}

	// Check file actually exists

	read_clock(&tstart);
	
	// Load config file
	parse_config_file( argv[2], &cfg );
	// Estimate my subset and combine
	double *error_grid = kspace_error( &cfg, qopt );
	// Write grid
	if ( proc_id == 0 )
		print_error_grid( error_grid, &cfg );

	read_clock(&tend);
	time = elapsed_time(&tstart, &tend)/1e6;

	MPI_Reduce( &time, &max_time, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD );
	if ( proc_id == 0 )
	{
		printf( "  %7.2f  ", max_time );
		fflush(stdout);
	}

	free( error_grid );
	error_grid = NULL;
	free_config( &cfg );

	MPI_Finalize();

	return 0;
}
