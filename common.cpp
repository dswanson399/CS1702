#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01  //Interaction Radius
#define min_r   (cutoff/100)  // Limits the distance a particle can be close to anothero particle (keeps sqrt(r) from blowing up) 
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//  Send is the structure containg info on the;
//       primary particle, 
//       neighbor particle, 
//       the minimum distance for interactions
//			Place holders for accumulating delta R (davg - distance for average) and the number of particles assessed (navg - number in the average).
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

	double dx = neighbor.x - particle.x;  // Find the Delta between X postions
	double dy = neighbor.y - particle.y;  // Find the delta between Y postions
	double r2 = dx * dx + dy * dy;        // Find Delta R^2
	if( r2 > cutoff*cutoff )  // If Delta R is bigger than the interaction radius, return to the calling routine
	  return;               //   without increaing the particle count (navg) and adding a delta R (davg)
		  
	if (r2 != 0) // If the focus particle and neighbor particle are not the same
	{
		if (r2/(cutoff*cutoff) < *dmin * (*dmin))  // If the Particle is closer than any other particle
			*dmin = sqrt(r2)/cutoff;					 // Set the minimum distance the this value
			
		(*davg) += sqrt(r2)/cutoff;  // Add the delta R to the average accumulator
		(*navg) ++;                  // Increase the particle count inside the intgeraction radius
// } original location of the (r2 != 0) closure brace
		
		r2 = fmax( r2, min_r*min_r );  // Limits the distance between particles
		double r = sqrt( r2 );         // Detemines the real distance between the particles
	
		//
		//  very simple short-range repulsive force
		//
		double coef = ( 1 - cutoff / r ) / r2 / mass;
		particle.ax += coef * dx;  // Force in the direction of dx
		particle.ay += coef * dy;  // Force in the direction of dy
	}  // Moved this braces to eliminate the work to figure out that a particle 
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;  // If p.x is outside the box refelct is back in
        p.vx = -p.vx;  								// Reflect velocity in x direction
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y; // If p.y is outside the box refelct is back in
        p.vy = -p.vy;  							  // Reflect velocity in y direction 
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
