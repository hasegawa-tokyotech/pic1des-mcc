//
//  World_1DES.cpp
//  es1d2v
//
//  Created by Jun Hasegawa on 2013/08/01.
//  Last Updated by Jun Hasegawa on 2019/07/12.
//
//  Copyright (c) 2013-2019 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <random>

#include "World_1DES.h"

#include "Particles_1DES.h"
#include "Field_1DES.h"

//#define PIC_USE_EXPDIST	// comment out if use uniform initial distribution

World_1DES::World_1DES( const char* fn )
: xd_(0), dt_(0), tstep_(0), tout_(0), ns_(0), nr_(0), nd_(0), xmin_(0), xmax_(0)
{
	using namespace std;

	ifstream fin( fn, ios::in );
	if ( !fin ) throw runtime_error( "can't open input file" );
	
	string s;
	fin >> s;
	if ( s != "#pic_1des" ) throw runtime_error( "illegal file format" );
	
	for (;;) {
		char c;
		char buf[255];
		fin.get( c );	// read first character
		if ( fin.eof() ) break;
		if ( c == '#' ) {
			fin.getline( buf, 255 );	// skip comment
			continue;
		}
		else if ( c == '$' ) {
			string kind;
			fin >> kind;
			if ( kind == "world" ) {
				fin >> buf >> tstep_
					>> buf >> tout_
					>> buf >> dt_
					>> buf >> ng_;
			}
			else if ( kind == "particles" ) {
				long np;
				double m, q, x0, x1, kt, vd, np0, kt0, vd0;
				fin >> buf >> np
					>> buf >> m
					>> buf >> q
					>> buf >> x0
					>> buf >> x1
					>> buf >> kt
					>> buf >> vd
					>> buf >> np0
					>> buf >> kt0
					>> buf >> vd0;
				particles_.push_back( new Particles_1DES( np, m, q, x0, x1, kt, vd, np0, kt0, vd0 ) );
				ns_++;
			}
			else if ( kind == "field" ) {
				static double x0 = 0;
				int nx, ns, bc0, bc1;
				double lx, n0, phi_n, bz;
				fin >> buf >> nx
					>> buf >> lx
					>> buf >> ns
					>> buf >> n0
					>> buf >> phi_n
					>> buf >> bz
                    >> buf >> bc0
                    >> buf >> bc1;
				field_.push_back( new Field_1DES( nx, lx, x0, n0, phi_n, bz, bc0, bc1 ) );
				nr_++;
				x0 += lx;	// set new left boundary
				xmax_ += lx;	// set maximum value of x
				if ( nr_ >= 2 ) {	// set potentials at boundaries
					double bias = field_[nr_-2]->phin();
					field_[nr_-1]->set_phi0( bias );
					field_[nr_-1]->set_phin( phi_n + bias );
				}
			}
			else {
				throw runtime_error( "unknown input parameter" );
			}
		}
		else continue;
	}
		
	for ( int i = 0; i < ns_; i++ )
		for ( int j = 0; j < nr_; j++ )
			particles_[i]->add_field( field_[j] );

	
	//dt_ = 0.1;	// set time step to 0.1/wp

}

World_1DES::~World_1DES()
{
	for ( int i = 0; i < ns_; i++ )
		delete particles_[i];
	
	for ( int i = 0; i < nr_; i++ )
		delete field_[i];
}

void World_1DES::run()
{
	using namespace std;
	
	const double mfp = 0.1*xmax_;	// mean free path of electrons

	// generate seed electrons
	//particles_[0]->generate_seed( dt_, xd_, 1 );
	particles_[0]->add_particle( xmax_, 0, 0 );

	// push back initial particles
	for ( int i = 0; i < ns_; i++ )
		particles_[i]->accelerate( -0.5*dt_ );

#ifdef PIC_USE_EXPDIST
	double alpha = -2.0;	// collisional ionization coefficient for H2 gas at 1 Pa
#endif
	
	int count = 0;	// counter for file output
    int count_eq = 0;   // counter for eqilibrium condition
	int neq = 3;	// number of establishment of ji=jeq needed for termination
	//int ts = 0;
	
	vector<double> j0(ns_,0), j1(ns_,0);	// containers for accumulating charges lost at left and right boudaries
	double jeq = particles_[1]->np0()*particles_[1]->q()/dt_;	// current density expected from number of generated particles

	// prepare output file for currents
	ofstream fout0( "current.txt", ios::out );
	if ( !fout0 ) throw runtime_error( "cannot open file: current.txt" );
	fout0 << "time";
	for ( int j = 0; j < ns_; j++ )
		fout0 << '\t' << "j0(pa" << j << ')' << '\t' << "j1(pa" << j << ')';
	fout0 << '\n';

	
	//for ( int i = 0; i <= tstep_; i++ ) {	// main loop
	for ( int i = 0; ; i++ ) {	// main loop

		// clear source
		for ( int j = 0; j < nr_; j++ )
			field_[j]->clear_source();
		 
		// set source
		for ( int k = 0; k < ns_; k++ )
			particles_[k]->weight( k );
			
		// solve field equation
		for ( int j = 0; j < nr_; j++ )
			field_[j]->solve();

		// output data to files
		if ( i%tout_ == 0 ) {
			char fn[255];
			for ( int j = 0; j < ns_; j++ ) {
				sprintf( fn, "pa%01d_%03d.txt", j, count );
				ofstream fout( fn, std::ios::out );
				particles_[j]->snapshot( fout );
				fout.close();
			}
			for ( int j = 0; j < nr_; j++ ) {
				sprintf( fn, "fi%01d_%03d.txt", j, count );
				ofstream fout( fn, std::ios::out );
				field_[j]->snapshot( fout );
				fout.close();
			}
			
			// output progress to console
			cout << '[' << setw(4) << setfill('0') << count << setfill(' ') << ']';
			fout0 << i*dt_;
			for ( int j = 0; j < ns_; j++ )
				cout << "\tnp" << j << "= " << setw(8) << right << particles_[j]->np();
            for ( int j = 0; j < ns_; j++ ) {
				//cout << scientific << setprecision(3);
				//cout << "\tj0[" << j << "]=\t" << setw(8) << right << j0[j]/tout_;
				//cout << "\tj1[" << j << "]=\t" << setw(8) << right << j1[j]/tout_;
				fout0 << '\t' << j0[j]/tout_ << '\t' << j1[j]/tout_;
            }
            fout0 << '\n';

            cout << scientific << setprecision(3);
            cout << "\tj0[0]=\t" << setw(8) << right << j0[0]/tout_;
            cout << "\tj1[1]=\t" << setw(8) << right << j1[1]/tout_;
			
			// detection of equilibrium condition ji = jeq
			double dj = fabs(j1[1]/tout_/jeq - 1);
            if ( dj < 1e-2 ) {
                if ( ++count_eq == neq ) {  // if ji = jeq is established three times, terminate calculation
                    cout << " *";
                    cout << "\nequilibrium is established: calculation stopped" << endl;
                    tstep_ = i;	// update total number of iteration
                    break;	// stop iteration
                }
                cout << " *";
			}
            for ( int j = 0; j < ns_; j++ ) j0[j] = j1[j] = 0;	// clear current values
			
            cout << endl;
            cout.flush();
			count++;
		}
			
        // move and accelerate particles
		for ( int j = 0; j < ns_; j++ )
			particles_[j]->advance( dt_ );
		
		
		if ( particles_[0]->np()==0 && particles_[1]->np()==0 )
			break;	// if no particles exist, stop calculation
		
		// inject particles
		//for ( int j = 0; j < ns_; j++ )
			//particles_[j]->inject( dt_ );

		// electron impact ionization
		long np = particles_[0]->np();	// number of electrons

		for ( int i = 0; i < np; i++ ) {
			//double mfp = 1.0/(ng_*sigma);	// mean free path
			
			double x = particles_[0]->x(i);
			double vx = particles_[0]->vx(i);
			double p = 1 - exp(-fabs(vx)*dt_/mfp);
			//if ( p > 1 ) throw runtime_error( "wrong probability value in World_1DES::run()" );
			// prepare random generator (C++11)
			static std::random_device seed_gen;
			static std::mt19937 engine(seed_gen());
			static std::uniform_real_distribution<> rand(0,1);
			double r1 = rand(engine);
			if ( r1 < p ) {
				double r2 = rand(engine);
				double xc = x + r2*(vx*dt_);	// position of collision event
				
				// generate an ion-electron pair
				particles_[0]->add_particle( xc, 0, 0 );	// generate an electron
				particles_[1]->add_particle( xc, 0, 0 );	// generate an ion
				
				// change position and velocity of projectile electron
				particles_[0]->set_x( i, xc );	// make vx = 0
				particles_[0]->set_vx( i, 0 );	// make vx = 0
			}
		}
		
		// generate plasma particles
		for ( int j = 0; j < ns_; j++ )
#ifdef PIC_USE_EXPDIST
			particles_[j]->generate_exp( dt_, alpha );
#else
			//particles_[j]->generate_uni( dt_ );
#endif
		
		// generate seed electrons
		//particles_[0]->generate_seed( dt_, xd_, 1 );
		//particles_[0]->add_particle( xmax_, 0, 0 );


		// accumulate charges lost at boundaries
		for ( int j = 0; j < ns_; j++ ) {
			j0[j] += particles_[j]->j0(i);
			j1[j] += particles_[j]->j1(i);
		}
	}
	
	/*
	// output conduction current
	ofstream fout( "current.txt", ios::out );
	if ( !fout ) throw runtime_error( "cannot open file: current.txt" );
	
	fout << "time";
	for ( int j = 0; j < ns_; j++ ) {
		fout << '\t' << "j0(pa" << j << ')'
			 << '\t' << "j1(pa" << j << ')';
	}
	fout << '\n';
	for ( int i = 0; i <= tstep_; i++ ) {
		fout << i*dt_;
		for ( int j = 0; j < ns_; j++ ) {
			fout << '\t' << particles_[j]->j0(i)
				 << '\t' << particles_[j]->j1(i);
		}
		fout << '\n';
	}
	*/
	
}
