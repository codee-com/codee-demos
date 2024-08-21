/**
 * ZPIC - em2d
 *
 * Weibel instability
 */

#include <stdlib.h>
#include "../simulation.h"

void sim_init( t_simulation* sim ){

	// Time step
	float dt = 0.07;
	float tmax = 35.0;

	// Simulation box
	int   nx[2]  = { 128, 128 };
	float box[2] = { 12.8, 12.8 };

	// Diagnostic frequency
	int ndump = 10;

    // Initialize particles
	const int n_species = 2;
	t_species* species = (t_species *) malloc( n_species * sizeof( t_species ));

	// Use 4x2 particles per cell
	int ppc[] = {4,2};

	// Initial fluid and thermal velocities
	t_part_data ufl[] = { 0.0, 0.0, 0.6 };
	t_part_data uth[] = { 0.1, 0.1, 0.1 };

	spec_new( &species[0], "electrons", -1.0, ppc, ufl, uth, nx, box, dt, NULL );

	ufl[2] = -ufl[2];
	spec_new( &species[1], "positrons", +1.0, ppc, ufl, uth, nx, box, dt, NULL );

	// Initialize Simulation data
	sim_new( sim, nx, box, dt, tmax, ndump, species, n_species );

}

void sim_report( t_simulation* sim ){

	// All magnetic field components
	emf_report( &sim->emf, BFLD, 0 );
	emf_report( &sim->emf, BFLD, 1 );
	emf_report( &sim->emf, BFLD, 2 );

}
