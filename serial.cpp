#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double davg, dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );  // Figure out how many particle the user wanted

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )  // Loops through NSTEPS of the simulation to move the particles
    {
		  navg = 0;  // A place to count the number of neighboring particles within the interaction distance
        davg = 0.0;  // A place to acucmulate distances to later average them
	     dmin = 1.0;  // A variable to hold the running minimum delta R to the focus particle
        //
        //  compute forces
        //
        for( int i = 0; i < n; i++ )   // Loop through all the particles as the focus particle
        {
            particles[i].ax = particles[i].ay = 0; //Set acceleration of the focus particle to Zero
            for (int j = 0; j < n; j++ )				// Loop though all the neighboring particles
					apply_force( particles[i], particles[j],&dmin,&davg,&navg); // Accumulate the force of each particle on the focus particle
        }
 
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) // Loop throguh al the particles and move them
            move( particles[i] );		

        if( find_option( argc, argv, "-no" ) == -1 )  // If the correctness check is turned on, prepare statistics for the check
        {
          //
          // Computing statistical data
          //
          if (navg)  // If there are particle in the count from the force subroutine 
			 {
            absavg +=  davg/navg; // Determine the absolute value of the average of the delta Rs across all particles
            nabsavg++;				 // 
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
	 
	 // The simulaiton is complete.
	 //
	 // Check for run time
    //
	 simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);


	// Check the statisitics to assess if the run was valid.
    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
		 // 
		 //  -the minimum distance absmin between 2 particles during the run of the simulation
		 //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
		 //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
		 //
		 //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
		 //
		 printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
		 if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
		 if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum ) fclose( fsum );    
    free( particles );
    if( fsave ) fclose( fsave );
    
    return 0;
}
