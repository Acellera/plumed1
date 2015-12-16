#include <stdio.h>
#include <stdlib.h>
#include "tcl.h"
#include "aceplug.h"
#include "stdlib.h"
#include "string.h"
#include <math.h>
#include "metadyn.h"

#define EXCH_TEMP 1
#define EXCH_SYSTEM 2

int do_remd_setup( struct aceplug_sim_t *s, int argc, char **argkey, char**argval ) ;

static unsigned long get_startstep( int argc, char **argkey, char **argval ) {
  int i;
  for( i=0; i<argc; i++ ) {
		if( !strcmp( argkey[i], "startstep" ) ) {
			return atol( argval[i] );
		}
	}
	return 0;		
}


static double get_arg( int argc, char **argkey, char **argval, char *needle ) {
  int i;
  for( i=0; i<argc; i++ ) {
    if( !strcmp( argkey[i], needle ) ) { return atof(argval[i]); }
  }
  return -1.0;
}
static int get_arg_int( int argc, char **argkey, char **argval, char *needle ) {
  int i;
  for( i=0; i<argc; i++ ) {
    if( !strcmp( argkey[i], needle ) ) { return atoi(argval[i]); }
  }
	printf ("# Argument not found : [%s]\n", needle );
	exit(-1);
}



aceplug_err_t  aceplug_init( 
	struct aceplug_sim_t *s,
	int argc,
	char **argkey,
	char **argval
) {
	int i;

	s->privdata = NULL;


	{
	 // Metadynamics init
		real box[3];
		int nrepl = 0;
		int repl = 0;
		double rte0 = 0.;
		double rteio = 0.;
		real *charge = NULL;
		real *mass = NULL;
		int i = 0 ;
		double dt = 0.;
		char *metainp = NULL;

	 for( i=0; i<argc; i++ ) {
    	if( !strcmp( argkey[i], "input" ) ) { metainp = argval[i]; }
  	}

		if( metainp == NULL ) {
			printf( "# PLUMED input file not specified\n" );
			exit(1);
		}


	  box[0]=s->box.x;    /* Box */
  	box[1]=s->box.y;
  	box[2]=s->box.z;

	  charge    = (real *)calloc(s->natoms,sizeof(real));
  	mass      = (real *)calloc(s->natoms,sizeof(real));
	  for(i=0;i<s->natoms;i++) {
  	  mass[i]   =  s->mass[i];
    	charge[i] =  s->charge[i];
  	}

		dt = s->timestep_fs;

	nrepl = 1;
	repl  = 1;
	rte0 = 298.0;
	rteio= rte0;

 	 printf("PLUMED: initializing with control file '%s', initial box (%f,%f,%f), timestep %f fs\n",metainp,box[0],box[1],box[2],dt);
  init_metadyn(s->natoms, charge, mass,
         dt, repl, nrepl, rte0, rteio, metainp, box );

	}

	return 0;
		
}



aceplug_err_t  aceplug_calcforces( struct aceplug_sim_t *s) {
	// Do nothing here
	struct remd_t *r  = (struct remd_t*) s->privdata;

	if ( s->plugin_load_positions() ) { return 1; }
	if ( s->plugin_load_forces() ) { return 1; }
  meta_force_calculation(s, 0, 0 );

	if ( ! s->plugin_update_forces() ) { return 1; }

	return 0;
}

aceplug_err_t aceplug_terminate( struct aceplug_sim_t *s ) {
	return 0;
/*
	struct remd_t *r  = (struct remd_t*) privdata;
	if(privdata==NULL) { return 0 ;}
	
	if( r->flog ) {
		fclose( r->flog );
	}
	free( r->Tall );
	free( privdata );

	return 0;
*/
}

aceplug_err_t aceplug_endstep( struct aceplug_sim_t *s) {
	struct remd_t *r = (struct remd_t*) s->privdata;
	int peer;
	int i;
	double energy_mine[3];
	double energy_peer[3];
	int do_exchange = 0;
	double vel_scale_factor = 0.;


//	if( ACEPLUG_OK !=	s->plugin_get_energy_temp( energy_mine ) ) {
	//	fprintf( stderr, "# REMD : Unable to get energies: not an output step. Set pluginfreq to be a multiple of energyfreq\n" );
//		return 0;
//	}

	return 0;
	
}

                 
