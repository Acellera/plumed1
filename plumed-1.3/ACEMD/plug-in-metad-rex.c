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

static int get_exchange_type( int argc, char **argkey, char **argval ) {
  int i;
  for( i=0; i<argc; i++ ) {
		if( !strcmp( argkey[i], "type" ) ) {
    if( !strcmp( argval[i], "temperature" ) ) { 
			printf( "# REMD exchange type is TEMPERATURE (trajectories continuous in time)\n" );
			return EXCH_TEMP; 
		}
    else if( !strcmp( argval[i], "system" ) ) { 
			printf( "# REMD exchange type is WHOLE SYSTEM (trajectories are discontinuous)\n" );
			return EXCH_SYSTEM; 
		}
		else {
			printf("# Exchange type invalid (temperature|system)\n" );
		}	
		}
  }
	printf("# Exchange type not specified (temperature|system)\n" );
	exit(-1);

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




struct remd_t {
//	double Ti;
	double *Tall;
	unsigned long count;
	int whoami;
	FILE *flog;
	int doing_remd;
	int exchange_type;
	unsigned long startstep;
//	int scale_after_exchange;
	int exchangefreq;
	double last_energy_pe;
};

aceplug_err_t  aceplug_init( 
	struct aceplug_sim_t *s,
	int argc,
	char **argkey,
	char **argval
) {
	int i;

	struct remd_t *p; 
	p = (struct remd_t*) malloc( sizeof( struct remd_t ) );
	s->privdata = p;

	if( s->ensemble_size > 1 ) {
		if ( s->ensemble_size % 2 ) {
			fprintf( stderr, "# REMD : Need an even number of replicas.\n" );
			return -1;
		}
 	 if ( 0 != do_remd_setup( s, argc, argkey, argval )  ) {
    exit(1);
	 }
	}

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
			printf( "# PUMED input file not specified\n" );
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
		rte0 = ((struct remd_t*) (s->privdata))->Tall[ s->ensemble_rank ];
	  rteio = rte0;
		nrepl = s->ensemble_size;
		repl  = s->ensemble_rank;


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
  meta_force_calculation(s, r->last_energy_pe, 1 );

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


	if( ACEPLUG_OK !=	s->plugin_get_energy_temp( energy_mine ) ) {
	//	fprintf( stderr, "# REMD : Unable to get energies: not an output step. Set pluginfreq to be a multiple of energyfreq\n" );
		return 0;
	}
	r->last_energy_pe = energy_mine[ ENERGY_PE ];

	if( s->ensemble_size == 1 ) { return 0; }

	if( s->step < r->startstep ) { return 0; }

	if( s->step % r->exchangefreq ) { return 0; } 


	// decide which neighbour to attempt exchange with
	// and get energies and temp
	// even rank is always receiver
	// and makes the metropolis test
	int receiving = 0;
	int err=0;

	energy_mine[ TEMPERATURE ] = r->Tall[ s->ensemble_rank ];

	if( r->count %2 ) {
		if( s->ensemble_rank % 2 ) {
			// send +
			peer = ( s->ensemble_rank +1 ) % s->ensemble_size;
		}
		else {
			// receive -
			peer = ( s->ensemble_rank -1 ); 
			if( peer < 0 ) { peer += s->ensemble_size; }
		}
	}
	else {
		// echange down
		if( s->ensemble_rank % 2 ) {
			// send -
			peer = ( s->ensemble_rank -1 ); 
			if( peer < 0 ) { peer += s->ensemble_size; }
		}
		else {
			// receive +
			peer = ( s->ensemble_rank +1 ) % s->ensemble_size;
		}
	}

//	receiving = (peer<s->ensemble_rank );
	if( peer > s->ensemble_rank ) { receiving = 1; }// r->Tall[peer] > r->Tall[ s->ensemble_rank ] ) { receiving = 1; }

	if( receiving ) {
		s->plugin_ensemble_send( peer, (void*) energy_mine, 3 * sizeof(double) );
		s->plugin_ensemble_recv( peer, (void*) energy_peer, 3 * sizeof(double) );
	}
	else {
		s->plugin_ensemble_recv( peer, (void*) energy_peer, 3 * sizeof(double) );
		s->plugin_ensemble_send( peer, (void*) energy_mine, 3 * sizeof(double) );
	}


  real *Vbias=NULL,*Vbiasx=NULL;
  snew(Vbias, s->ensemble_size);
  snew(Vbiasx, s->ensemble_size);

  for(i=0; i< s->ensemble_size; i++) {Vbias[i]=0.;Vbiasx[i]=0.;}

  if(logical.remd) {
    // each replica sends its position in CV space to each other, so that everyone knows where everyone is, then
    // each replica calculates its bias (summing its hills) and evaluate it in its own position and in the peer
    // position : 
    // Vbias[me]   == bias felt by me on my hills;    Vbiasx[0]    == bias felt by peer on my hills
    // Vbias[peer] == bias felt by peer on his hills; Vbiasx[peer] == bias felt by me on peer's hills
    // in case I don't exchange (I'm a replica on the edge and it is not my turn) then pass peer=-1

    ptmetad_vbias(peer,Vbias,Vbiasx);
  }

  

	double p = 0.;

//printf("£ STEP %lu exchange with %d\n", s->step, peer );
//printf("£ %f %f : %f %f\n", energy_mine[ENERGY_PE], energy_mine[TEMPERATURE], energy_peer[ENERGY_PE], energy_peer[TEMPERATURE] );

		double cond = .0;
	if( receiving ) {
    // determine whether an exchange should happen
    double delta=0.0;

    // META: calculate delta energy contribution
    real beta_mine, beta_peer, delta_mine, delta_peer;
    delta_mine = Vbias[s->ensemble_rank] - Vbiasx[s->ensemble_rank];
    delta_peer = Vbias[peer] - Vbiasx[peer];
    beta_mine = 1.0 / ( MD_BOLTZMAN * energy_mine[ TEMPERATURE ] );
    beta_peer = 1.0 / ( MD_BOLTZMAN * energy_peer[ TEMPERATURE ] );
    delta = - (beta_mine * delta_mine + beta_peer * delta_peer);

    delta += (beta_mine - beta_peer)*(energy_peer[ ENERGY_PE ]-energy_mine[ ENERGY_PE ]);

    // Metropolis criterion test here
		if ( delta < 0. ) { cond = 1.; }
		else {
      cond  = exp( - delta );
    }

    p    = (float)rand() / RAND_MAX;
      // Metropolis criterion test here
    printf("# Acceptance P=%.2e\n", cond );
    if ( p <= cond ) {
      p=1.0;
    } else {
      p=0.0;
    }

		// p==1.0 if an exchange is to happen.
		// tell the peer in this exchange attempt whether to exchange
		// then, if true, exchange

		s->plugin_ensemble_send( peer, (void*) &p ,1 * sizeof(double) );
	}
	else {
		s->plugin_ensemble_recv( peer, (void*) &p ,1 * sizeof(double) );

	}


		if( p == 1.0 ) { // do coordinate exchange

			r->last_energy_pe = energy_peer[ ENERGY_PE ];

			printf( "# REMD Exchange %d<->%d at step %lu  p=%.3f T%d=%.2f T%d=%.2f\n ", s->ensemble_rank, peer, s->step, cond, s->ensemble_rank,  energy_mine[TEMPERATURE], peer, energy_peer[TEMPERATURE] );
			fflush( stdout );
	
			if( receiving ) {
				double tmp = 0;
				s->plugin_ensemble_send( peer, (void*) &(r->whoami), sizeof(int) );
				s->plugin_ensemble_recv( peer, (void*) &(r->whoami), sizeof(int) );
			} else {
				int tmp = r->whoami;
				s->plugin_ensemble_recv( peer, (void*) &(r->whoami),  sizeof(int));
				s->plugin_ensemble_send( peer, (void*) &(tmp), sizeof(int) );
			}


			switch( r->exchange_type ) {
				case EXCH_TEMP:
					r->Tall[ s->ensemble_rank ] = energy_peer[ TEMPERATURE ];
					s->plugin_set_temperature( r->Tall[ s->ensemble_rank ] );
					vel_scale_factor = sqrt ( r->Tall[ s->ensemble_rank ]/ r->Tall[ peer ] );
				break;
				case EXCH_SYSTEM:
					s->plugin_exchange_system( peer );
					vel_scale_factor = sqrt ( r->Tall[ s->ensemble_rank ]/ r->Tall[ peer ] );
				break;
			}
		
			do_exchange = 1;
		}

	if( do_exchange ) { 
		printf("# REMD Scaling vel to %fK\n", r->Tall[ s->ensemble_rank ] );
  	s->plugin_scale_velocities( vel_scale_factor ); 
	}

	// Now log the trajectories after this exchange
	// so that a continuous traj cpuld be reassembled later
 
	if( s->ensemble_rank == 0 ) {
		fprintf( r->flog, "%lu\t%lu\t", s->step, r->whoami );
		int t;
		for( i=1; i < s->ensemble_size; i++ ) {
			s->plugin_ensemble_recv( i, (void*) &(t) , sizeof(int) );
			fprintf( r->flog, "%lu\t", t );
		}
		fprintf( r->flog, "\n" );
		fflush( r->flog );
	}
	else {
			s->plugin_ensemble_send( 0,  &(r->whoami), sizeof(int) );
	}

	r->count++;

	return 0;
	
}

int do_remd_setup( struct aceplug_sim_t *s, int argc, char **argkey, char**argval ) {

// MJH REMD setup start
  int i;

  struct remd_t *p;
  p = (struct remd_t*) malloc( sizeof( struct remd_t ) );
  p->doing_remd = 0;
	p->exchange_type = EXCH_TEMP;

  s->privdata = p;
  if( s->ensemble_size > 1 ) {
    printf( "# PLUMED Ensemble replica exchange on %d replicas.\n", s->ensemble_size );
    p->doing_remd = 1;
  }
  else {
    printf( "# PLUMED Metadynamics on single simulation.\n" );
    return 0;
  }



  p->count = 0;
//  p->T0 = get_arg( argc , argkey, argval, "T0" );
//  p->k  = get_arg( argc, argkey, argval, "k" );
  p->whoami = s->ensemble_rank;
	p->exchangefreq = get_arg_int( argc, argkey, argval, "exchangefreq" );

  p->flog = NULL;
  if( s->ensemble_rank == 0 ) {
    for(i=0;i<argc;i++ ) {
      if(!strcmp(argkey[i], "log" ) ) {
        p->flog = fopen( argval[i], "w" ); // don't try to append, even if this is a reset (we can't tell)
      }
    }
  }

  if( p->flog == NULL ) {
    if( !s->ensemble_rank ) {
      fprintf( stderr, "# REMD : and log=filename for trace log\n");
      return -1;
    }
  }

	
//	p->exchange_type = get_exchange_type(  argc, argkey, argval );
	p->exchange_type = EXCH_SYSTEM;
//	p->scale_after_exchange = get_scale_after_exchange( argc, argkey, argval );
	p->startstep     = get_startstep( argc, argkey, argval );

  p->Tall = (double*) malloc( sizeof(double) * s->ensemble_size );
//  if( !s->ensemble_rank ) {
    printf( "#REMD Setup:\n# Replica\tTemperature\n" );
    for( i=0; i< s->ensemble_size; i++ ) {
      char tstring[10];
      snprintf( tstring, 9, "T%d", i );
      p->Tall[i] = get_arg( argc, argkey, argval, tstring );
//      if( i == s->ensemble_rank ) { p = p->Tall[i];}
      printf( "# %d\t%f\n", i, p->Tall[i] );
    }
    printf("# \n");
//  }

    for( i=0; i< s->ensemble_size; i++ ) {
			if( p->Tall[i]<=0. || ( i>0 && (p->Tall[i] <= p->Tall[i-1] ) ) ) {
				printf("# Invalid temperature for replica %d (%f)\n", i, p->Tall[i] );
				exit(-1);
			}
		}
  printf("# Rank %d Temperature %f\n", s->ensemble_rank, p->Tall[ s->ensemble_rank ] );

	

  s->plugin_set_temperature( p->Tall[ s->ensemble_rank ] );

	return 0;
// MJH REMD setup done
}
                  
