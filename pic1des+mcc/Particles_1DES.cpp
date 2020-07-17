//
//  Particles_1DES.cpp
//  es1d2v
//
//  Created by Jun Hasegawa on 2013/03/05.
//  Last Updated by Jun Hasegawa on 2019/07/12.
//
//  Copyright (c) 2019 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#include "Particle_1DES.h"
#include "Particles_1DES.h"
#include "Field_1DES.h"

#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <random>

#include "phys_const.h"

Particles_1DES::Particles_1DES( long np, double m, double q, double x0, double x1, double kt, double vd, long np0, double kt0, double vd0 )
: pid_(0), nr_(0), xmin_(0), xmax_(0), m_(m), q_(q), np0_(np0), kt0_(kt0), vd0_(vd0), n0_(0), n1_(0)
{
    // check input parameters
	if ( m == 0 ) throw std::runtime_error( "bad parameters" );
    
    // set left boundary ( x = 0 )
    boundary_.push_back( 0 );
	
	// prepare random generator (C++11)
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());
	//std::mt19937 engine(1);	//
	std::uniform_real_distribution<> rand(0,1);
	std::normal_distribution<> norm(0,sqrt(kt/m));	// Maxwell distribution
	
	// initialize particles
	for ( long i = 0; i < np; i++ ) {
        
        // initial position
		double r0 = rand(engine);
		double x = r0*(x1-x0) + x0;
        
		// initial velocity
		double um = norm(engine);
		
        // get direction
        //double r2 = genrand_real2();
		double r2 = rand(engine);
		double theta = si::pi*(2*r2-1);
        double vx = um*cos(theta);
        double vy = um*sin(theta);

        vx += vd;   // add drift velocity
        
        Particle_1DES* p = new Particle_1DES( pid_++, m, q, x, vx, vy );
        ps_.push_back( p );
    }
	
	j0_.push_back( 0 );
	j1_.push_back( 0 );
    
}

Particles_1DES::~Particles_1DES()
{
    std::vector<Particle_1DES*>::iterator p;
    
    // delete all particles
	for ( p = ps_.begin(); p != ps_.end(); p++ ) {
        delete *p;
	}
}

void Particles_1DES::add_particle( double m, double q, double x, double vx, double vy )
{
    Particle_1DES* p = new Particle_1DES( pid_++, m, q, x, vx, vy );
    ps_.push_back( p );
}

void Particles_1DES::add_particle( double x, double vx, double vy )
{
    Particle_1DES* p = new Particle_1DES( pid_++, m_, q_, x, vx, vy );
    ps_.push_back( p );
}

void Particles_1DES::advance( double dt )
{
	long np = ps_.size();
    
    // advance all particles
    for ( long i = 0; i < np; i++ ) {
		double x = ps_[i]->x();
        // get region number
        int k = get_region( x );
        
		// get field values at particle position
        double ex = field_[k]->ex( x );
        double bz = field_[k]->bz();
        
		// accelerate and move particle
        ps_[i]->accelerate( dt, ex, bz );
        ps_[i]->move( dt );
    }
	
    // remonve/reinject/reflect particles that move beyond boundaries
    long n = 0;
	double q0 = 0;
	double q1 = 0;
    for ( long i = 0; i < np; i++ ) {
        double x = ps_[i]->x();
        if ( x < xmin_ ) {	// remove particles if they go beyond left boundary
            q0 += ps_[i]->q();
			n0_++;
            if ( field_[0]->bc0() == 0 ) { // remove particles
                delete ps_[i];
            }
            else if ( field_[0]->bc0() == 1 ) {    // periodic case: reinject particles from opposite boundary
                ps_[i]->set_x( xmax_ + x );
                ps_[n] = ps_[i];
                n++;
            }
            else {  // synmetric case: reflect particles
                ps_[i]->set_x( 2*xmax_ - x );
                ps_[i]->set_vx( -ps_[i]->vx() );
                ps_[i]->set_vy( -ps_[i]->vy() );
                ps_[n] = ps_[i];
                n++;
            }
            continue;
        }
        else if ( x > xmax_ ) {
			q1 += ps_[i]->q();
			n1_++;
            if ( field_[0]->bc1() == 0 ) { // remove particles
                delete ps_[i];
            }
            else if ( field_[0]->bc1() == 1 ) {
                ps_[i]->set_x( x - xmax_ );
                ps_[n] = ps_[i];
                n++;
            }
            else {
                ps_[i]->set_x( 2*xmax_ - x );
                ps_[i]->set_vx( -ps_[i]->vx() );
                ps_[i]->set_vy( -ps_[i]->vy() );
                ps_[n] = ps_[i];
                n++;
            }
			continue;
        }
        ps_[n] = ps_[i];
        n++;
    }
	
	// record conduction currents
	j0_.push_back( q0/dt );
	j1_.push_back( q1/dt );
	
    ps_.resize(n);
}

void Particles_1DES::accelerate( double dt )
{
    long np = ps_.size();
    for ( long i = 0; i < np; i++ ) {
		double x = ps_[i]->x(); // particle position
		int k = get_region( x );	// get region number

        // get field values at particle position
        double ex = field_[k]->ex( x );
        double bz = field_[k]->bz();

        // accelerate particle
        ps_[i]->accelerate( dt, ex, bz );
	}
}


void Particles_1DES::weight( int s ) const
{
    long np = ps_.size();

    for ( long i = 0; i < np; i++ ) {
		int k = get_region( ps_[i]->x() );	// get region number
		field_[k]->set_source( s, ps_[i]->q(), ps_[i]->x() );
	}
}

/*
void Particles_1DES::inject( double dt  )
{
	// prepare random generator (C++11)
	static std::random_device seed_gen;
	static std::mt19937 engine(seed_gen());
	static std::uniform_real_distribution<> rand(0,1);
	static std::normal_distribution<> norm(0,vt0_/sqrt(2));	// Maxwell distribution

	// initialize particles
	for ( long i = 0; i < np0_; i++ ) {
		// initial velocity
		double um = norm(engine);
		
		// get direction
		double r2 = rand(engine);
		double phi = si::pi*(2*r2-1);
		double vx = um*cos(phi);
		double vy = um*sin(phi);
		
		vx += vd0_;   // add drift velocity
		if ( vx < 0 ) vx = -vx;
		
		// get initial position
        double rdt = rand(engine) * dt;
        double x = vx * rdt;
		
        Particle_1DES* p = new Particle_1DES( pid_++, m_, q_, x, vx, vy );
		
		// push back injected particles
		int k = get_region( x );
		double ex = field_[k]->ex( x );
        double bz = field_[k]->bz();
        p->accelerate( -0.5*dt, ex, bz );
		
        ps_.push_back( p );
    }    
}
 */

void Particles_1DES::generate_uni( double dt )
{
	// prepare random generator (C++11)
	static std::random_device seed_gen;
	static std::mt19937 engine(seed_gen());
	static std::uniform_real_distribution<> rand(0,1);
	static std::normal_distribution<> norm(0,sqrt(kt0_/m_));	// Maxwell distribution
	
	// initialize particles
	for ( long i = 0; i < np0_; i++ ) {
		// initial velocity
		double um = norm(engine);
		
		// get direction
		double r0 = rand(engine);
		double phi = si::pi*(2*r0-1);
		double vx = um*cos(phi);
		double vy = um*sin(phi);
		
		vx += vd0_;   // add drift velocity
		if ( vx < 0 ) vx = -vx;
		
		// get initial position
		double r1 = rand(engine);
		double x = r1*(xmax_-xmin_) + xmin_;
		
		Particle_1DES* p = new Particle_1DES( pid_++, m_, q_, x, vx, vy );
		
		// advance generated particle
		int k = get_region( x );
		double ex = field_[k]->ex( x );
		double bz = field_[k]->bz();
		double r2 = rand(engine);	// distribute particles uniformly on temporal axis
		p->accelerate( (0.5-r2)*dt, ex, bz );
		p->move( (1-r2)*dt );
		
		ps_.push_back( p );
	}
}

double Particles_1DES::pdist_exp( double x, double a )
{
	double e0 = exp(a*xmin_);
	double e1 = exp(a*xmax_);
	double ans = (exp(a*x)-e0)/(e1-e0);
	return ans;
}

void Particles_1DES::generate_exp( double dt, double a )
{
	// prepare random generator (C++11)
	static std::random_device seed_gen;
	static std::mt19937 engine(seed_gen());
	static std::uniform_real_distribution<> rand(0,1);
	static std::normal_distribution<> norm(0,sqrt(kt0_/m_));	// Maxwell distribution
	
	// initialize particles
	for ( long i = 0; i < np0_; i++ ) {
		// initial velocity
		double um = norm(engine);
		
		// get direction
		double r0 = rand(engine);
		double phi = si::pi*(2*r0-1);
		double vx = um*cos(phi);
		double vy = um*sin(phi);
		
		vx += vd0_;   // add drift velocity
		if ( vx < 0 ) vx = -vx;
		
		// get initial position based on a exponential function
		double r1 = rand(engine);
		double x0 = xmin_;
		double x1 = xmax_;
		double x = 0.5*(x0 + x1);
		for (;;) {
			double f0 = pdist_exp(x0,a) - r1;
			double f1 = pdist_exp(x,a) - r1;
			if ( f0*f1 < 0 ) x1 = x;
			else x0 = x;
			if ( fabs(x0-x1)/x0 < 1e-5 ) break;
			x = 0.5*(x0 + x1);
		}
		
		Particle_1DES* p = new Particle_1DES( pid_++, m_, q_, x, vx, vy );
		
		// advance generated particle
		int k = get_region( x );
		double ex = field_[k]->ex( x );
		double bz = field_[k]->bz();
		double r2 = rand(engine);	// distribute particles uniformly on temporal axis
		p->accelerate( (0.5-r2)*dt, ex, bz );
		p->move( (1-r2)*dt );
		
		ps_.push_back( p );
	}
}

/*
void Particles_1DES::generate_seed( double dt, double x0, long n )
{
	Particle_1DES* p = new Particle_1DES( pid_++, m_, q_, x0, 0, 0 );
	
	ps_.push_back( p );
}
*/

void Particles_1DES::snapshot( std::ostream& out, bool title ) const
{
	if ( title ) out << "#id\tm\tq\tx\tvx\tvy\n";
	out << scientific << setprecision(6);
    long np = ps_.size();
    for ( long i = 0; i < np; i++ ) {
		out << ps_[i]->id() << '\t'
			<< ps_[i]->m() << '\t'
            << ps_[i]->q() << '\t'
            << ps_[i]->x() << '\t'
            << ps_[i]->vx() << '\t'
            << ps_[i]->vy() << '\n';
	}
}

void Particles_1DES::add_field( Field_1DES* field )
{
	field_.push_back( field );
	field->add_particles();
	
    xmax_ += field->lx();	// update right boundary position
	boundary_.push_back( xmax_ );
	nr_++;
}

int Particles_1DES::get_region( double x ) const
{
    int k = 0;
	for ( int i = 1; i < nr_; i++ ) {
		if( x < boundary_[i] ) break;
		k++;
	}
	return k;
}


