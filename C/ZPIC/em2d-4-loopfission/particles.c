// clang-format off
/*
 *  particles.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 11/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particles.h"

#include "current.h"
#include "emf.h"
#include "one-malloc.h"
#include "random.h"

#include "zdf.h"
#include "timer.h"

static double _spec_time = 0.0;
static double _spec_npush = 0.0;

void spec_sort( t_species *spec );

/**
 * Returns the total time spent pushing particles (includes boundaries and moving window)
 * @return  Total time in seconds
 */
double spec_time( void )
{
	return _spec_time;
}

/**
 * Returns the performance achieved by the code (push time)
 * @return  Performance in seconds per particle
 */
double spec_perf( void )
{
	return (_spec_npush > 0 )? _spec_time / _spec_npush: 0.0;
}

/*********************************************************************************************
 Initialization
 *********************************************************************************************/

/**
 * Sets the momentum of the range of particles supplieds using a thermal distribution
 * @param spec  Particle species
 * @param start Index of the first particle to set the momentum
 * @param end   Index of the last particle to set the momentum
 */
void spec_set_u( t_species* spec, const int start, const int end )
{
	int i;    

	for (i = start; i <= end; i++) {
		spec->part[i].ux = spec -> ufl[0] + spec -> uth[0] * rand_norm(); 
		spec->part[i].uy = spec -> ufl[1] + spec -> uth[1] * rand_norm(); 
		spec->part[i].uz = spec -> ufl[2] + spec -> uth[2] * rand_norm(); 
	}

}	

void spec_set_x( t_species* spec, int range[2][2] )
{

	int i, j, k, ip;
	
	float* poscell;
	float start, end;
	
	// Calculate particle positions inside the cell
	const int npc = spec->ppc[0]*spec->ppc[1];
	t_part_data const dpcx = 1.0f/spec->ppc[0];
	t_part_data const dpcy = 1.0f/spec->ppc[1];
	
	poscell = malloc( 2 * npc * sizeof( t_part_data ) );
	ip = 0;
	for (j =0; j<spec->ppc[1]; j++) {
		for (i=0; i<spec->ppc[0]; i++) {
			poscell[ip]   = dpcx * ( i + 0.5 );
			poscell[ip+1] = dpcy * ( j + 0.5 );
			ip+=2;
		}
	}

	ip = spec -> np;
	
	// Set position of particles in the specified grid range according to the density profile
	switch ( spec -> density.type ) {
	case STEP: // Step like density profile
		
		// Get edge position normalized to cell size;
		start = spec -> density.start / spec -> dx[0] - spec -> n_move;

		for (j = range[1][0]; j <= range[1][1]; j++) {
			for (i = range[0][0]; i <= range[0][1]; i++) {

				for (k=0; k<npc; k++) {
					if ( i + poscell[2*k] > start ) {
						spec->part[ip].ix = i;
						spec->part[ip].iy = j;
						spec->part[ip].x = poscell[2*k];
						spec->part[ip].y = poscell[2*k+1];
						ip++;
					}
				}
			}
		}
		break;

	case SLAB: // Slab like density profile
		
		// Get edge position normalized to cell size;
		start = spec -> density.start / spec -> dx[0] - spec -> n_move;
		end   = spec -> density.end / spec -> dx[0] - spec -> n_move;

		for (j = range[1][0]; j <= range[1][1]; j++) {
			for (i = range[0][0]; i <= range[0][1]; i++) {

				for (k=0; k<npc; k++) {
					if ( i + poscell[2*k] > start &&  i + poscell[2*k] < end ) {
						spec->part[ip].ix = i;
						spec->part[ip].iy = j;
						spec->part[ip].x = poscell[2*k];
						spec->part[ip].y = poscell[2*k+1];
						ip++;
					}
				}
			}
		}
		break;

	default: // Uniform density
		for (j = range[1][0]; j <= range[1][1]; j++) {
			for (i = range[0][0]; i <= range[0][1]; i++) {

				for (k=0; k<npc; k++) {
					spec->part[ip].ix = i;
					spec->part[ip].iy = j;
					spec->part[ip].x = poscell[2*k];
					spec->part[ip].y = poscell[2*k+1];
					ip++;
				}
			}
		}
	}
	
	spec -> np = ip;
	
	free(poscell);
	
}

void spec_inject_particles( t_species* spec, int range[2][2] )
{
	int start = spec -> np;

	// Get maximum number of particles to inject
	int np_inj = ( range[0][1] - range[0][0] + 1 ) * ( range[1][1] - range[1][0] + 1 ) * 
	             spec -> ppc[0] * spec -> ppc[1];

	// Check if buffer is large enough and if not reallocate
	if ( spec -> np + np_inj > spec -> np_max ) {
        spec -> np_max = (( spec -> np_max + np_inj )/1024 + 1) * 1024;
		spec -> part = realloc( (void*) spec -> part, spec -> np_max * sizeof(t_part) );
	}

	// Set particle positions
	spec_set_x( spec, range );
	
	// Set momentum of injected particles
	spec_set_u( spec, start, spec -> np - 1 );

}

void spec_new( t_species* spec, char name[], const t_part_data m_q, const int ppc[], 
			  const t_part_data *ufl, const t_part_data * uth,
			  const int nx[], t_part_data box[], const float dt, t_density* density )
{

	int i, npc;
	
	// Species name
	strncpy( spec -> name, name, MAX_SPNAME_LEN );
	
	npc = 1;
	// Store species data
	for (i=0; i<2; i++) {
		spec->nx[i] = nx[i];
		spec->ppc[i] = ppc[i];
		npc *= ppc[i];
		
		spec->box[i] = box[i];
		spec->dx[i] = box[i] / nx[i];
	}
	
	spec -> m_q = m_q;
	spec -> q = copysign( 1.0f, m_q ) / npc;

	spec -> dt = dt;
	
	// Initialize particle buffer
	spec->np_max = 0;
	spec->part = NULL;
	
	
	// Initialize density profile
	if ( density ) {
		spec -> density = *density;
		if ( spec -> density.n == 0. ) spec -> density.n = 1.0;
	} else {
		// Default values
		spec -> density = (t_density) { .type = UNIFORM, .n = 1.0 };
	}

	// Initialize temperature profile
	if ( ufl ) {
		for(i=0; i<3; i++) spec -> ufl[i] = ufl[i];
	} else {
		for(i=0; i<3; i++) spec -> ufl[i] = 0;
	}

	// Density multiplier
	spec ->q *= fabsf( spec -> density.n );

	if ( uth ) {
		for(i=0; i<3; i++) spec -> uth[i] = uth[i];
	} else {
		for(i=0; i<3; i++) spec -> uth[i] = 0;
	}

	// Reset iteration number
	spec -> iter = 0;

	// Reset moving window information
	spec -> moving_window = 0;
	spec -> n_move = 0;

    // Inject initial particle distribution
    spec -> np = 0;

    int range[2][2];
    range[0][0] = 0;
    range[0][1] = nx[0]-1;
    range[1][0] = 0;
    range[1][1] = nx[1]-1;

    spec_inject_particles( spec, range );

}

void spec_move_window( t_species *spec ){

	if ((spec->iter * spec->dt ) > (spec->dx[0] * (spec->n_move + 1)))  {
        
        // shift all particles left
        // particles leaving the box will be removed later
        int i;
        for( i = 0; i < spec->np; i++ ) {
        	spec->part[i].ix--;
        }

		// Increase moving window counter
		spec -> n_move++;

        // Inject particles in the right edge of the simulation box
	    int range[2][2];
      range[0][0] = spec->nx[0]-1;
      range[0][1] = spec->nx[0]-1;
      range[1][0] = 0;
      range[1][1] = spec->nx[1]-1;
	    spec_inject_particles( spec, range );

	}

}

void spec_delete( t_species* spec )
{
	free(spec->part);
	spec->np = -1;
}


/*********************************************************************************************
 
 Cuurent deposition
 
 *********************************************************************************************/

void dep_current_esk( int ix0, int iy0, int di, int dj,
							 t_part_data x0, t_part_data y0, t_part_data x1, t_part_data y1, 
							 t_part_data qvx, t_part_data qvy, t_part_data qvz, 
							 t_current *current )
{

	int i, j;
	t_fld S0x[4], S0y[4], S1x[4], S1y[4], DSx[4], DSy[4];
	t_fld Wx[16], Wy[16], Wz[16];
	
	S0x[0] = 0.0f;
	S0x[1] = 1.0f - x0;
	S0x[2] = x0;
	S0x[3] = 0.0f;

	S0y[0] = 0.0f;
	S0y[1] = 1.0f - y0;
	S0y[2] = y0;
	S0y[3] = 0.0f;
	
	for (i=0; i<4; i++) {
		S1x[i] = 0.0f;
		S1y[i] = 0.0f;
	}
	
	S1x[ 1 + di ] = 1.0f - x1;
	S1x[ 2 + di ] = x1;

	S1y[ 1 + dj ] = 1.0f - y1;
	S1y[ 2 + dj ] = y1;

	for (i=0; i<4; i++) {
		DSx[i] = S1x[i] - S0x[i];
		DSy[i] = S1y[i] - S0y[i];
	}
	
	for (j=0; j<4; j++) {
		for (i=0; i<4; i++) {
			Wx[i + 4*j] = DSx[i] * ( S0y[j] + DSy[j]/2.0f );
			Wy[i + 4*j] = DSy[j] * ( S0x[i] + DSx[i]/2.0f );
			Wz[i + 4*j] = S0x[i] * S0y[j] + DSx[i]*S0y[j]/2.0f +
			              S0x[i]*DSy[j]/2.0f + DSx[i]*DSy[j]/3.0f;
		}
	}
		
	// jx
	const int nrow = current -> nrow;
	t_vfld* restrict const J = current -> J;
	
	for (j=0; j<4; j++) {
		t_fld c;
		
		c = -qvx * Wx[4*j];
		J[ ix0 - 1 + (iy0 - 1 + j)*nrow ].x += c;
		for (i=1; i<4; i++) {
			c -= qvx * Wx[i+4*j];
			J[ ix0 + i - 1 + (iy0 -1 + j)*nrow ].x += c;
		}
	}

	// jy
	for (i=0; i<4; i++) {
		t_fld c;

		c = -qvy * Wy[i];
		J[ ix0 + i - 1 + (iy0 - 1)*nrow ].y += c;
		for (j=1; j<4; j++) {
			c -= qvy * Wy[i+4*j];
			J[ ix0 + i - 1 + (iy0 -1 + j)*nrow ].y += c;
		}
	}
	
	// jz
	for (j=0; j<4; j++) {
		for (i=0; i<4; i++) {
			J[ ix0 + i - 1 + (iy0 -1 + j)*nrow ].z += qvz * Wz[ i + 4*j ];
		}
	}
	
	
}

void dep_current_zamb(int ix, int iy, int di, int dj, 
				      float x0, float y0, float dx, float dy,
					  float qnx, float qny, float qvz,
					  int current_nrow,
            		  t_fld * restrict const Jx, t_fld * restrict const Jy, t_fld * restrict const Jz)
{
	// Split the particle trajectory
	
	typedef struct {
		float x0, x1, y0, y1, dx, dy, qvz;
		int ix, iy;
	} t_vp;
	
	t_vp vp[3];
	int vnp = 1;
	
	// split 
	vp[0].x0 = x0;
	vp[0].y0 = y0;

	vp[0].dx = dx;
	vp[0].dy = dy;
	
	vp[0].x1 = x0+dx;
	vp[0].y1 = y0+dy;
	
	vp[0].qvz = qvz/2.0;

	vp[0].ix = ix;
	vp[0].iy = iy;
		
	// x split
	if ( di != 0 ) {
		
		//int ib = ( di+1 )>>1;
		int ib = ( di == 1 );
		
		float delta = (x0+dx-ib)/dx;
		
		// Add new particle
		vp[1].x0 = 1-ib;
		vp[1].x1 = (x0 + dx) - di;
		vp[1].dx = dx*delta;
		vp[1].ix = ix + di;
		
		float ycross = y0 + dy*(1.0f-delta);

		vp[1].y0 = ycross; 
		vp[1].y1 = vp[0].y1;
		vp[1].dy = dy*delta;	
		vp[1].iy = iy;
		
		vp[1].qvz = vp[0].qvz*delta;
		
		// Correct previous particle
		vp[0].x1 = ib;
		vp[0].dx *= (1.0f-delta);
		
		vp[0].dy *= (1.0f-delta);
		vp[0].y1  = ycross;
		
		vp[0].qvz *= (1.0f-delta);
		
        vnp++;		
	}
	
	// ysplit
	if ( dj != 0 ) {
		int isy = 1 - ( vp[0].y1<0.0f || vp[0].y1>=1.0f );
		
		// int jb = ( dj+1 )>>1; 
		int jb = (dj == 1);
		
		// The static analyser gets confused by this but it is correct
		float delta = (vp[isy].y1-jb)/vp[isy].dy;
		
		// Add new particle
		vp[vnp].y0 = 1-jb;
		vp[vnp].y1 = vp[isy].y1 - dj;
		vp[vnp].dy = vp[isy].dy*delta;
		vp[vnp].iy = vp[isy].iy + dj;
		
		float xcross = vp[isy].x0 + vp[isy].dx*(1.0f-delta); 
		
		vp[vnp].x0 = xcross;
		vp[vnp].x1 = vp[isy].x1;
		vp[vnp].dx = vp[isy].dx*delta;
		vp[vnp].ix = vp[isy].ix;

		vp[vnp].qvz = vp[isy].qvz*delta;
		
		// Correct previous particle
		vp[isy].y1  = jb;
		vp[isy].dy *= (1.0f-delta);
		
		vp[isy].dx *= (1.0f-delta);
		vp[isy].x1  = xcross;
		
		vp[isy].qvz *= (1.0f-delta);
		
		// Correct extra vp if needed
 		if ( isy < vnp -1) {
			vp[1].y0 -= dj;
			vp[1].y1 -= dj;
			vp[1].iy += dj;
		} 
		vnp++;		
	}

	// Deposit virtual particle currents
	int k;
	const int nrow = current_nrow;

	for (k = 0; k < vnp; k++) {
		float S0x[2], S1x[2], S0y[2], S1y[2];
		float wl1, wl2;
		float wp1[2],wp2[2];
		
		S0x[0] = 1.0f - vp[k].x0;
		S0x[1] = vp[k].x0;

		S1x[0] = 1.0f - vp[k].x1;
		S1x[1] = vp[k].x1;

		S0y[0] = 1.0f - vp[k].y0;
		S0y[1] = vp[k].y0;

		S1y[0] = 1.0f - vp[k].y1;
		S1y[1] = vp[k].y1;

		wl1 = qnx * vp[k].dx;
		wl2 = qny * vp[k].dy;
		
		wp1[0] = 0.5f*(S0y[0] + S1y[0]);
		wp1[1] = 0.5f*(S0y[1] + S1y[1]);
		
		wp2[0] = 0.5f*(S0x[0] + S1x[0]);
		wp2[1] = 0.5f*(S0x[1] + S1x[1]);
		
		Jx[ vp[k].ix + nrow*vp[k].iy     ] += wl1 * wp1[0];
		Jx[ vp[k].ix + nrow*(vp[k].iy+1) ] += wl1 * wp1[1];

		Jy[ vp[k].ix   + nrow*vp[k].iy ] += wl2 * wp2[0];
		Jy[ vp[k].ix+1 + nrow*vp[k].iy ] += wl2 * wp2[1];

		Jz[ vp[k].ix   + nrow*vp[k].iy    ] += vp[k].qvz * (S0x[0]*S0y[0]+S1x[0]*S1y[0]+(S0x[0]*S1y[0]-S1x[0]*S0y[0])/2.0f);
		Jz[ vp[k].ix+1 + nrow*vp[k].iy    ] += vp[k].qvz * (S0x[1]*S0y[0]+S1x[1]*S1y[0]+(S0x[1]*S1y[0]-S1x[1]*S0y[0])/2.0f);
		Jz[ vp[k].ix   + nrow*(vp[k].iy+1)] += vp[k].qvz * (S0x[0]*S0y[1]+S1x[0]*S1y[1]+(S0x[0]*S1y[1]-S1x[0]*S0y[1])/2.0f);
		Jz[ vp[k].ix+1 + nrow*(vp[k].iy+1)] += vp[k].qvz * (S0x[1]*S0y[1]+S1x[1]*S1y[1]+(S0x[1]*S1y[1]-S1x[1]*S0y[1])/2.0f);
	}

}

/*********************************************************************************************
 
 Sorting
 
 *********************************************************************************************/

void spec_sort( t_species* spec )
{
	int *idx, *npic;

	int ncell = spec->nx[0]*spec->nx[1];
	
	// Allocate index memory
	idx  = malloc(spec->np*sizeof(int));

	// Allocate temp. array with number of particles in cell
	npic = malloc( ncell * sizeof(int));
	memset( npic, 0, ncell * sizeof(int));

	// Generate sorted index
    int i;
	for (i=0; i<spec->np; i++) {
		idx[i] = spec->part[i].ix + spec->part[i].iy * spec->nx[0];
		npic[idx[i]]++;
	}
	
	int isum = 0, j;
	for (i=0; i<ncell; i++) {
		j = npic[i];
		npic[i] = isum;
		isum += j;
	}
	
	for (i=0; i< spec->np; i++) {
		j = idx[i];
		idx[i] = npic[j]++;
	}
	
	// free temp. array
	free(npic);
/*
	// Rearrange particle buffer
	t_part *tmp = malloc( spec->np * sizeof( t_part ) );
	for (i=0; i< spec->np; i++) {
        tmp[idx[i]] = spec->part[i];
	}
	free(spec->part);
	spec->part = tmp;
*/

	// low mem
	for (i=0; i < spec->np; i++) {
		t_part tmp;
		int k;
		
		k = idx[i];
		while ( k > i ) {
			int t;
			
			tmp = spec->part[k];
			spec->part[k] = spec->part[i];
			spec->part[i] = tmp;
			
			t = idx[k];
			idx[k] = -1;
			k = t;
		}
	}

	free(idx);

}


/*********************************************************************************************
 
 Particle advance
 
 *********************************************************************************************/

void interpolate_fld(const t_fld *restrict const Ex, const t_fld *restrict const Ey, const t_fld *restrict const Ez,
                     const t_fld *restrict const Bx, const t_fld *restrict const By, const t_fld *restrict const Bz,
                     const int nrow,
				     
					 // Struct part
					 const int ix, const int iy, const t_fld x, const t_fld y,

                     // Struct to fields of Ep and Bp
                     t_fld *Epx, t_fld *Epy, t_fld *Epz,
                     t_fld *Bpx, t_fld *Bpy, t_fld *Bpz)
{
	register int i, j, ih, jh;
	register t_fld w1, w2, w1h, w2h;
	
	i = ix;
	j = iy;
	
	w1 = x;
	w2 = y;
	
	ih = (w1 <0.5f)? -1 : 0;
	jh = (w2 <0.5f)? -1 : 0;
	
	// w1h = w1 - 0.5f - ih;
	// w2h = w2 - 0.5f - jh;
	w1h = w1 + ((w1 <0.5f)?0.5f:-0.5f);
	w2h = w2 + ((w2 <0.5f)?0.5f:-0.5f);

	ih += i;
	jh += j;
	
  *Epx = ( Ex[ih +     j *nrow] * (1.0f - w1h) + Ex[ih+1 +     j*nrow] * w1h ) * (1.0f -  w2 ) +
        ( Ex[ih + (j +1)*nrow] * (1.0f - w1h) + Ex[ih+1 + (j+1)*nrow] * w1h ) * w2;

  *Epy = ( Ey[i  +     jh*nrow] * (1.0f -  w1) + Ey[i+1  +     jh*nrow] * w1 ) * (1.0f - w2h ) +
        ( Ey[i  + (jh+1)*nrow] * (1.0f -  w1) + Ey[i+1  + (jh+1)*nrow] * w1 ) * w2h;

  *Epz = ( Ez[i  +     j *nrow] * (1.0f - w1) + Ez[i+1 +     j*nrow] * w1 ) * (1.0f - w2 ) +
        ( Ez[i  + (j +1)*nrow] * (1.0f - w1) + Ez[i+1 + (j+1)*nrow] * w1 ) * w2;

  *Bpx = ( Bx[i  +     jh*nrow] * (1.0f - w1) + Bx[i+1 +     jh*nrow] * w1 ) * (1.0f - w2h ) +
        ( Bx[i  + (jh+1)*nrow] * (1.0f - w1) + Bx[i+1 + (jh+1)*nrow] * w1 ) * w2h;

  *Bpy = ( By[ih +     j*nrow] * (1.0f - w1h) + By[ih+1 +     j*nrow] * w1h ) * (1.0f - w2 ) +
        ( By[ih + (j +1)*nrow] * (1.0f - w1h) + By[ih+1 + (j+1)*nrow] * w1h ) * w2;

  *Bpz = ( Bz[ih +     jh*nrow] * (1.0f - w1h) + Bz[ih+1 +     jh*nrow] * w1h ) * (1.0f - w2h ) +
        ( Bz[ih + (jh+1)*nrow] * (1.0f - w1h) + Bz[ih+1 + (jh+1)*nrow] * w1h ) * w2h;

}	

int ltrim( t_part_data x )
{
	return ( x >= 1.0f ) - ( x < 0.0f );
}

// clang-format on
void spec_advance_outlining(
    // Spec data
    const int spec_np, const t_part_data spec_q,
    const t_part_data qnx, const t_part_data qny, const t_part_data tem,
    const t_part_data dt_dx, const t_part_data dt_dy, int *restrict spec_ix,
    int *restrict spec_iy, t_part_data *restrict spec_x,
    t_part_data *restrict spec_y, t_part_data *restrict spec_ux,
    t_part_data *restrict spec_uy, t_part_data *restrict spec_uz,
    // Emf data
    const int emf_nrow, t_fld *restrict Ex, t_fld *restrict Ey,
    t_fld *restrict Ez, t_fld *restrict Bx, t_fld *restrict By,
    t_fld *restrict Bz,
    // Current data
    const int current_nrow, t_fld *restrict current_Jx, t_fld *restrict current_Jy,
    t_fld *restrict current_Jz) {
  t_part_data qvz;
  t_part_data x1, y1;

  int di, dj;
  float dx, dy;

  // Advance particles
  for (int i = 0; i < spec_np; i++) {
    t_fld Epx, Epy, Epz;
    t_fld Bpx, Bpy, Bpz;

    t_part_data utx, uty, utz;
    t_part_data ux, uy, uz, rg;
    t_part_data gtem, otsq;

    // Load particle momenta
    ux = spec_ux[i];
    uy = spec_uy[i];
    uz = spec_uz[i];

    // interpolate fields
    interpolate_fld(Ex, Ey, Ez, Bx, By, Bz, emf_nrow,
					spec_ix[i], spec_iy[i], spec_x[i], spec_y[i],
					&Epx, &Epy, &Epz, &Bpx, &Bpy, &Bpz);

    // advance u using Boris scheme
    Epx *= tem;
    Epy *= tem;
    Epz *= tem;

    utx = ux + Epx;
    uty = uy + Epy;
    utz = uz + Epz;

    // Perform first half of the rotation
    gtem = tem / sqrtf(1.0f + utx * utx + uty * uty + utz * utz);

    Bpx *= gtem;
    Bpy *= gtem;
    Bpz *= gtem;

    otsq = 2.0f / (1.0f + Bpx * Bpx + Bpy * Bpy + Bpz * Bpz);

    ux = utx + uty * Bpz - utz * Bpy;
    uy = uty + utz * Bpx - utx * Bpz;
    uz = utz + utx * Bpy - uty * Bpx;

    // Perform second half of the rotation

    Bpx *= otsq;
    Bpy *= otsq;
    Bpz *= otsq;

    utx += uy * Bpz - uz * Bpy;
    uty += uz * Bpx - ux * Bpz;
    utz += ux * Bpy - uy * Bpx;

    // Perform second half of electric field acceleration
    ux = utx + Epx;
    uy = uty + Epy;
    uz = utz + Epz;

    // Store new momenta
    spec_ux[i] = ux;
    spec_uy[i] = uy;
    spec_uz[i] = uz;
  }

  for (int i = 0; i < spec_np; i++) {
    // Load particle momenta
    t_part_data ux = spec_ux[i];
    t_part_data uy = spec_uy[i];
    t_part_data uz = spec_uz[i];

    // push particle
    t_part_data rg = 1.0f / sqrtf(1.0f + ux * ux + uy * uy + uz * uz);

    dx = dt_dx * rg * ux;
    dy = dt_dy * rg * uy;

    x1 = spec_x[i] + dx;
    y1 = spec_y[i] + dy;

    di = ltrim(x1);
    dj = ltrim(y1);

    x1 -= di;
    y1 -= dj;

    qvz = spec_q * uz * rg;

    // dep_current_zamb(spec_ix[i], spec_iy[i], di, dj, spec_x[i], spec_y[i],
    //                  dx, dy, qnx, qny, qvz, current_nrow, current_Jx,
    //                  current_Jy, current_Jz);

    //////////////////////////////////////
    //   dep_current_zamb inlined START //
    //////////////////////////////////////
    const float ix = spec_ix[i];
    const float iy = spec_iy[i];
    const float x0 = spec_x[i];
    const float y0 = spec_y[i];

    // Split the particle trajectory
    float vp_x0[3], vp_x1[3], vp_y0[3], vp_y1[3], vp_dx[3], vp_dy[3], vp_qvz[3];
    int vp_ix[3], vp_iy[3];

    int vnp = 1;

    // split
    vp_x0[0] = x0;
    vp_y0[0] = y0;

    vp_dx[0] = dx;
    vp_dy[0] = dy;

    vp_x1[0] = x0 + dx;
    vp_y1[0] = y0 + dy;

    vp_qvz[0] = qvz / 2.0;

    vp_ix[0] = ix;
    vp_iy[0] = iy;

    // x split
    if (di != 0) {

      // int ib = ( di+1 )>>1;
      int ib = (di == 1);

      float delta = (x0 + dx - ib) / dx;

      // Add new particle
      vp_x0[1] = 1 - ib;
      vp_x1[1] = (x0 + dx) - di;
      vp_dx[1] = dx * delta;
      vp_ix[1] = ix + di;

      float ycross = y0 + dy * (1.0f - delta);

      vp_y0[1] = ycross;
      vp_y1[1] = vp_y1[0];
      vp_dy[1] = dy * delta;
      vp_iy[1] = iy;

      vp_qvz[1] = vp_qvz[0] * delta;

      // Correct previous particle
      vp_x1[0] = ib;
      vp_dx[0] *= (1.0f - delta);

      vp_dy[0] *= (1.0f - delta);
      vp_y1[0] = ycross;

      vp_qvz[0] *= (1.0f - delta);

      vnp++;
    }

    // ysplit
    if (dj != 0) {
      int isy = 1 - (vp_y1[0] < 0.0f || vp_y1[0] >= 1.0f);

      // int jb = ( dj+1 )>>1;
      int jb = (dj == 1);

      // The static analyser gets confused by this but it is correct
      float delta = (vp_y1[isy] - jb) / vp_dy[isy];

      // Add new particle
      vp_y0[vnp] = 1 - jb;
      vp_y1[vnp] = vp_y1[isy] - dj;
      vp_dy[vnp] = vp_dy[isy] * delta;
      vp_iy[vnp] = vp_iy[isy] + dj;

      float xcross = vp_x0[isy] + vp_dx[isy] * (1.0f - delta);

      vp_x0[vnp] = xcross;
      vp_x1[vnp] = vp_x1[isy];
      vp_dx[vnp] = vp_dx[isy] * delta;
      vp_ix[vnp] = vp_ix[isy];

      vp_qvz[vnp] = vp_qvz[isy] * delta;

      // Correct previous particle
      vp_y1[isy] = jb;
      vp_dy[isy] *= (1.0f - delta);

      vp_dx[isy] *= (1.0f - delta);
      vp_x1[isy] = xcross;

      vp_qvz[isy] *= (1.0f - delta);

      // Correct extra vp if needed
      if (isy < vnp - 1) {
        vp_y0[1] -= dj;
        vp_y1[1] -= dj;
        vp_iy[1] += dj;
      }
      vnp++;
    }

    // Deposit virtual particle currents
    int k;
    const int nrow = current_nrow;

    for (k = 0; k < vnp; k++) {
      float S0x[2], S1x[2], S0y[2], S1y[2];
      float wl1, wl2;
      float wp1[2], wp2[2];

      S0x[0] = 1.0f - vp_x0[k];
      S0x[1] = vp_x0[k];

      S1x[0] = 1.0f - vp_x1[k];
      S1x[1] = vp_x1[k];

      S0y[0] = 1.0f - vp_y0[k];
      S0y[1] = vp_y0[k];

      S1y[0] = 1.0f - vp_y1[k];
      S1y[1] = vp_y1[k];

      wl1 = qnx * vp_dx[k];
      wl2 = qny * vp_dy[k];

      wp1[0] = 0.5f * (S0y[0] + S1y[0]);
      wp1[1] = 0.5f * (S0y[1] + S1y[1]);

      wp2[0] = 0.5f * (S0x[0] + S1x[0]);
      wp2[1] = 0.5f * (S0x[1] + S1x[1]);

      current_Jx[vp_ix[k] + nrow * vp_iy[k]] += wl1 * wp1[0];
      current_Jx[vp_ix[k] + nrow * (vp_iy[k] + 1)] += wl1 * wp1[1];

      current_Jy[vp_ix[k] + nrow * vp_iy[k]] += wl2 * wp2[0];
      current_Jy[vp_ix[k] + 1 + nrow * vp_iy[k]] += wl2 * wp2[1];

      current_Jz[vp_ix[k] + nrow * vp_iy[k]] +=
          vp_qvz[k] * (S0x[0] * S0y[0] + S1x[0] * S1y[0] +
                       (S0x[0] * S1y[0] - S1x[0] * S0y[0]) / 2.0f);
      current_Jz[vp_ix[k] + 1 + nrow * vp_iy[k]] +=
          vp_qvz[k] * (S0x[1] * S0y[0] + S1x[1] * S1y[0] +
                       (S0x[1] * S1y[0] - S1x[1] * S0y[0]) / 2.0f);
      current_Jz[vp_ix[k] + nrow * (vp_iy[k] + 1)] +=
          vp_qvz[k] * (S0x[0] * S0y[1] + S1x[0] * S1y[1] +
                       (S0x[0] * S1y[1] - S1x[0] * S0y[1]) / 2.0f);
      current_Jz[vp_ix[k] + 1 + nrow * (vp_iy[k] + 1)] +=
          vp_qvz[k] * (S0x[1] * S0y[1] + S1x[1] * S1y[1] +
                       (S0x[1] * S1y[1] - S1x[1] * S0y[1]) / 2.0f);
    }

    /////////////////////////////////////
    //   dep_current_zamb inlined  END //
    /////////////////////////////////////

    // Store results
    spec_x[i] = x1;
    spec_y[i] = y1;
    spec_ix[i] += di;
    spec_iy[i] += dj;
  }
}
// clang-format off

void pack_spec(const int size, const t_part *part, int **ix, int **iy,
               t_part_data **x, t_part_data **y, t_part_data **ux,
               t_part_data **uy, t_part_data **uz) {

  *ix = malloc(size * sizeof(int));
  *iy = malloc(size * sizeof(int));

  *x = malloc(size * sizeof(t_part_data));
  *y = malloc(size * sizeof(t_part_data));
  *ux = malloc(size * sizeof(t_part_data));
  *uy = malloc(size * sizeof(t_part_data));
  *uz = malloc(size * sizeof(t_part_data));

  assert(*ix);
  assert(*iy);
  assert(*x);
  assert(*y);
  assert(*ux);
  assert(*uy);
  assert(*uz);

  for (int i = 0; i < size; ++i) {
    (*ix)[i] = part[i].ix;
    (*iy)[i] = part[i].iy;

    (*x)[i] = part[i].x;
    (*y)[i] = part[i].y;
    (*ux)[i] = part[i].ux;
    (*uy)[i] = part[i].uy;
    (*uz)[i] = part[i].uz;
  }
}

void unpack_spec(const int size, t_part *part, int **ix, int **iy,
                 t_part_data **x, t_part_data **y, t_part_data **ux,
                 t_part_data **uy, t_part_data **uz) {

  for (int i = 0; i < size; ++i) {
    part[i].ix = (*ix)[i];
    part[i].iy = (*iy)[i];

    part[i].x = (*x)[i];
    part[i].y = (*y)[i];
    part[i].ux = (*ux)[i];
    part[i].uy = (*uy)[i];
    part[i].uz = (*uz)[i];
  }

  free(*ix);
  free(*iy);
  free(*x);
  free(*y);
  free(*ux);
  free(*uy);
  free(*uz);

  *ix = *iy = NULL;
  *x = *y = *ux = *uy = *uz = NULL;
}

// clang-format on
void spec_advance(t_species *spec, t_emf *emf, t_current *current,
                  t_all_data *data) {
  int i;
  t_part_data qnx, qny, qvz;

  uint64_t t0;
  t0 = timer_ticks();

  const t_part_data tem = 0.5 * spec->dt / spec->m_q;
  const t_part_data dt_dx = spec->dt / spec->dx[0];
  const t_part_data dt_dy = spec->dt / spec->dx[1];

  // Auxiliary values for current deposition
  qnx = spec->q * spec->dx[0] / spec->dt;
  qny = spec->q * spec->dx[1] / spec->dt;

  const int nx0 = spec->nx[0];
  const int nx1 = spec->nx[1];

  /////////////////////
  // Pack spec->part //
  /////////////////////
  int *spec_ix, *spec_iy;
  t_part_data *spec_x, *spec_y, *spec_ux, *spec_uy, *spec_uz;

  pack_spec(spec->np_max, spec->part, &spec_ix, &spec_iy, &spec_x, &spec_y,
            &spec_ux, &spec_uy, &spec_uz);

  copy_in_all_data(data, current, emf);

  spec_advance_outlining( // spec data
      spec->np, spec->q, qnx, qny, tem, dt_dx, dt_dy, spec_ix,
      spec_iy, spec_x, spec_y, spec_ux, spec_uy, spec_uz,
      // emf data
      emf->nrow, data->emf_data.Ex, data->emf_data.Ey, data->emf_data.Ez,
      data->emf_data.Bx, data->emf_data.By, data->emf_data.Bz,
      // current data
      current->nrow, data->current_data.Jx, data->current_data.Jy,
      data->current_data.Jz);

  copy_out_all_data(data, current, emf);

  /////////////////////
  //   Unpack spec   //
  /////////////////////
  unpack_spec(spec->np_max, spec->part, &spec_ix, &spec_iy, &spec_x, &spec_y,
              &spec_ux, &spec_uy, &spec_uz);

  // Advance internal iteration number
  spec->iter += 1;
  _spec_npush += spec->np;

  // Check for particles leaving the box
  if (spec->moving_window) {

    // Move simulation window if needed
    spec_move_window(spec);

    // Use absorbing boundaries along x, periodic along y
    i = 0;
    while (i < spec->np) {
      if ((spec->part[i].ix < 0) || (spec->part[i].ix >= nx0)) {
        spec->part[i] = spec->part[--spec->np];
        continue;
      }
      spec->part[i].iy += ((spec->part[i].iy < 0) ? nx1 : 0) -
                          ((spec->part[i].iy >= nx1) ? nx1 : 0);
      i++;
    }

  } else {
    // Use periodic boundaries in both directions
    for (i = 0; i < spec->np; i++) {
      spec->part[i].ix += ((spec->part[i].ix < 0) ? nx0 : 0) -
                          ((spec->part[i].ix >= nx0) ? nx0 : 0);
      spec->part[i].iy += ((spec->part[i].iy < 0) ? nx1 : 0) -
                          ((spec->part[i].iy >= nx1) ? nx1 : 0);
    }
  }

  // Sort species at every 16 time steps
  if (!(spec->iter % 16))
    spec_sort(spec);

  _spec_time += timer_interval_seconds(t0, timer_ticks());
}
// clang-format off
/*********************************************************************************************
 
 Charge Deposition
 
 *********************************************************************************************/


void spec_deposit_charge( const t_species* spec, t_part_data* charge )
{
	int i,j;
	
	// Charge array is expected to have 1 guard cell at the upper boundary
	int nrow = spec -> nx[0] + 1;
	t_part_data q = spec -> q;
	
	for (i=0; i<spec->np; i++) {
		int idx = spec->part[i].ix + nrow*spec->part[i].iy;
		t_fld w1, w2;
		
		w1 = spec->part[i].x;
		w2 = spec->part[i].y;
		
		charge[ idx            ] += ( 1.0f - w1 ) * ( 1.0f - w2 ) * q;
		charge[ idx + 1        ] += (        w1 ) * ( 1.0f - w2 ) * q;
		charge[ idx     + nrow ] += ( 1.0f - w1 ) * (        w2 ) * q;
		charge[ idx + 1 + nrow ] += (        w1 ) * (        w2 ) * q;
	}

	// Correct boundary values

	// x
	if ( ! spec -> moving_window ){
		for (j = 0; j < spec -> nx[1] + 1; j++) {
			charge[ 0 + j*nrow ] += charge[ spec -> nx[0] + j*nrow ];
		}
	}
	
	// y - Periodic boundaries
	for (i = 0; i < spec->nx[0]+1; i++) {
		charge[ i + 0 ] += charge[ i + spec -> nx[1] * nrow ];
	}

}

/*********************************************************************************************
 
 Diagnostics
 
 *********************************************************************************************/

void spec_rep_particles( const t_species *spec )
{
	
	t_zdf_file part_file;

	int i;
	
	const char * quants[] = {
	    "x1","x2",
	    "u1","u2","u3"
	};

	const char * units[] = {
	    "c/\\omega_p", "c/\\omega_p",
	    "c","c","c"
	};

    t_zdf_iteration iter;
    iter.n = spec->iter;
    iter.t = spec -> iter * spec -> dt;
    iter.time_units = "1/\\omega_p";

	// Allocate buffer for positions
	
	t_zdf_part_info info;

	info.name = (char *) spec -> name;
	info.nquants = 5;
	info.quants = (char **) quants;
	info.units = (char **) units;
	info.np = spec ->np;

	// Create file and add description
	zdf_part_file_open( &part_file, &info, &iter, "PARTICLES" );

	// Add positions and generalized velocities
	size_t size = ( spec -> np ) * sizeof( float );
	float* data = malloc( size );

	// x1
	for( i = 0; i < spec ->np; i++ )
		data[i] = ( spec -> n_move + spec -> part[i].ix + spec -> part[i].x ) * spec -> dx[0];
	zdf_part_file_add_quant( &part_file, quants[0], data, spec ->np );

	// x2
	for( i = 0; i < spec ->np; i++ )
		data[i] = (spec -> part[i].iy + spec -> part[i].y ) * spec -> dx[1];
	zdf_part_file_add_quant( &part_file, quants[1], data, spec ->np );

	// ux
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].ux;
	zdf_part_file_add_quant( &part_file, quants[2], data, spec ->np );

	// uy
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uy;
	zdf_part_file_add_quant( &part_file, quants[3], data, spec ->np );

	// uz
	for( i = 0; i < spec ->np; i++ ) data[i] = spec -> part[i].uz;
	zdf_part_file_add_quant( &part_file, quants[4], data, spec ->np );

	free( data );

	zdf_close_file( &part_file );
}	


void spec_rep_charge( const t_species *spec )
{
	t_part_data *buf, *charge, *b, *c;
	size_t size;
	int i, j;
	
	// Add 1 guard cell to the upper boundary
	size = ( spec -> nx[0] + 1 ) * ( spec -> nx[1] + 1 ) * sizeof( t_part_data );
	charge = malloc( size );
	memset( charge, 0, size );
	
	// Deposit the charge
	spec_deposit_charge( spec, charge );
	
	// Compact the data to save the file (throw away guard cells)
	size = ( spec -> nx[0] ) * ( spec -> nx[1] );
	buf = malloc( size * sizeof( float ) );
	
	b = buf;
	c = charge;
	for( j = 0; j < spec->nx[1]; j++) {
		for ( i = 0; i < spec->nx[0]; i++ ) {
			b[i] = c[i];
		}
		b += spec->nx[0];
		c += spec->nx[0] + 1;
	}
	
	free( charge );

  t_zdf_grid_axis axis[2];
  axis[0].min = 0.0;
  axis[0].max = spec->box[0];
  axis[0].label = "x_1";
  axis[0].units = "c/\\omega_p";

  axis[1].min = 0.0,
  axis[1].max = spec->box[1];
  axis[1].label = "x_2";
  axis[1].units = "c/\\omega_p";

  t_zdf_grid_info info;
  info.ndims = 2;
  info.label = "charge";
  info.units = "n_e";
  info.axis  = axis;

  info.nx[0] = spec->nx[0];
  info.nx[1] = spec->nx[1];

  t_zdf_iteration iter;
  iter.n = spec->iter;
  iter.t = spec -> iter * spec -> dt;
  iter.time_units = "1/\\omega_p";

	zdf_save_grid( buf, &info, &iter, spec->name );	


	free( buf );
}	


void spec_pha_axis( const t_species *spec, int i0, int np, int quant, float *axis )
{
	int i;
	
	switch (quant) {
		case X1:
			for (i = 0; i < np; i++) 
				axis[i] = ( spec -> part[i0+i].x + spec -> part[i0+i].ix ) * spec -> dx[0];
			break;
		case X2:
			for (i = 0; i < np; i++) 
				axis[i] = ( spec -> part[i0+i].y + spec -> part[i0+i].iy ) * spec -> dx[1];
			break;
		case U1:
			for (i = 0; i < np; i++) 
				axis[i] = spec -> part[i0+i].ux;
			break;
		case U2:
			for (i = 0; i < np; i++) 
				axis[i] = spec -> part[i0+i].uy;
			break;
		case U3:
			for (i = 0; i < np; i++) 
				axis[i] = spec -> part[i0+i].uz;
			break;
	}
}

const char * spec_pha_axis_units( int quant ) {
	switch (quant) {
		case X1:
		case X2:
			return("c/\\omega_p");
			break;
		case U1:
		case U2:
		case U3:
			return("m_e c");
	}
	return("");
}


void spec_deposit_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2], float* restrict buf )
{
	const int BUF_SIZE = 1024;
	float pha_x1[BUF_SIZE], pha_x2[BUF_SIZE];


	const int nrow = pha_nx[0];

	const int quant1 = rep_type & 0x000F;
	const int quant2 = (rep_type & 0x00F0)>>4;

	const float x1min = pha_range[0][0];
	const float x2min = pha_range[1][0];

	const float rdx1 = pha_nx[0] / ( pha_range[0][1] - pha_range[0][0] );
	const float rdx2 = pha_nx[1] / ( pha_range[1][1] - pha_range[1][0] );

	for ( int i = 0; i<spec->np; i+=BUF_SIZE ) {
		int np = ( i + BUF_SIZE > spec->np )? spec->np - i : BUF_SIZE;

		spec_pha_axis( spec, i, np, quant1, pha_x1 );
	    spec_pha_axis( spec, i, np, quant2, pha_x2 );

		for ( int k = 0; k < np; k++ ) {

			float nx1 = ( pha_x1[k] - x1min ) * rdx1;
			float nx2 = ( pha_x2[k] - x2min ) * rdx2;

			int i1 = (int)(nx1 + 0.5f);
			int i2 = (int)(nx2 + 0.5f);

			float w1 = nx1 - i1 + 0.5f;
			float w2 = nx2 - i2 + 0.5f;

			int idx = i1 + nrow*i2;

			if ( i2 >= 0 && i2 < pha_nx[1] ) {

				if (i1 >= 0 && i1 < pha_nx[0]) {
					buf[ idx ] += (1.0f-w1)*(1.0f-w2)*spec->q;
				}

				if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
					buf[ idx + 1 ] += w1*(1.0f-w2)*spec->q;
				}
			}

			idx += nrow;
			if ( i2+1 >= 0 && i2+1 < pha_nx[1] ) {

				if (i1 >= 0 && i1 < pha_nx[0]) {
					buf[ idx ] += (1.0f-w1)*w2*spec->q;
				}

				if (i1+1 >= 0 && i1+1 < pha_nx[0] ) {
					buf[ idx + 1 ] += w1*w2*spec->q;
				}
			}

		}

	}
}

void spec_rep_pha( const t_species *spec, const int rep_type,
			  const int pha_nx[], const float pha_range[][2] )
{

	char const * const pha_ax_name[] = {"x1","x2","x3","u1","u2","u3"};
	char pha_name[64];

	// Allocate phasespace buffer
	float* restrict buf = malloc( pha_nx[0] * pha_nx[1] * sizeof( float ));
	memset( buf, 0, pha_nx[0] * pha_nx[1] * sizeof( float ));

	// Deposit the phasespace
	spec_deposit_pha( spec, rep_type, pha_nx, pha_range, buf );

	// save the data in hdf5 format
	int quant1 = rep_type & 0x000F;
	int quant2 = (rep_type & 0x00F0)>>4;

    const char * pha_ax1_units = spec_pha_axis_units(quant1);
    const char * pha_ax2_units = spec_pha_axis_units(quant2);

	sprintf( pha_name, "%s%s", pha_ax_name[quant1-1], pha_ax_name[quant2-1] );

    t_zdf_grid_axis axis[2];
    axis[0] = (t_zdf_grid_axis) {
    	.min = pha_range[0][0],
    	.max = pha_range[0][1],
    	.label = (char *) pha_ax_name[ quant1 - 1 ],
    	.units = (char *) pha_ax1_units
    };

    axis[1] = (t_zdf_grid_axis) {
    	.min = pha_range[1][0],
    	.max = pha_range[1][1],
    	.label = (char *) pha_ax_name[ quant2 - 1 ],
    	.units = (char *) pha_ax2_units
    };

    t_zdf_grid_info info = {
    	.ndims = 2,
    	.label = pha_name,
    	.units = "a.u.",
    	.axis  = axis
    };

    info.nx[0] = pha_nx[0];
    info.nx[1] = pha_nx[1];

    t_zdf_iteration iter = {
    	.n = spec->iter,
    	.t = spec -> iter * spec -> dt,
    	.time_units = "1/\\omega_p"
    };

	zdf_save_grid( buf, &info, &iter, spec->name );

	// Free temp. buffer
	free( buf );

}

void spec_report( const t_species *spec, const int rep_type, 
				  const int pha_nx[], const float pha_range[][2] )
{
	
	switch (rep_type & 0xF000) {
		case CHARGE:
			spec_rep_charge( spec );
			break;

		case PHA:
			spec_rep_pha( spec, rep_type, pha_nx, pha_range );
			break;

		case PARTICLES:
			spec_rep_particles( spec );
			break;
	}
	
	
}
