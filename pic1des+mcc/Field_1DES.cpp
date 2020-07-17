//
//  Field_1DES.cpp
//  es1d2v
//
//  Created by Jun Hasegawa on 2013/03/06.
//  Last Updated by Jun Hasegawa on 2019/07/12.
//
//  Copyright (c) 2019 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#include <fstream>
#include <stdexcept>
#include <iomanip>

#include "Field_1DES.h"
#include "Particles_1DES.h"

const int ns_max = 5;

Field_1DES::Field_1DES( int nx, double lx, double x0, double n0, double phi_n, double bz, int bc0, int bc1 )
: nx_(nx), lx_(lx), x0_(x0), x1_(x0+lx), dx_(lx/nx), n0_(n0), bz_(bz), bc0_(bc0), bc1_(bc1), ns_(0), sig0_(0), sig1_(0)
{
    // resize field values
    rho_.resize( nx+1 );
    phi_.resize( nx+1 );
    ex_.resize( nx+1 );
	
	// resize density matrices
	rhok_.resize( ns_max, nx+1 );
		
    // clear field values
    for ( int i = 0; i <= nx; i++ ) {
        rho_.at(i) = 0;
        phi_.at(i) = 0;
        ex_.at(i) = 0;
    }

    // set boundary conditions
    phi_.at(0) = 0; // potential at left boundary
    phi_.at(nx) = phi_n;    // potential at right boundary
    
    // set coefficients for Poisson solver
    a_.resize( nx );
    b_.resize( nx );
    c_.resize( nx );
    r_.resize( nx );
    for ( int i = 1; i <= nx-1; i++ ) {
		a_[i] = 1;
		b_[i] = -2;
		c_[i] = 1;
	}
    f_ = -dx_*dx_/n0;
}

void Field_1DES::solve()
{
	//if ( bc0_ == 0 && bc1_ == 0 ) {	// both boundaries are conducting walls
		// solve Poisson equation
		r_[1] 	= f_*rho_[1] - phi_[0];
		r_[nx_-1] = f_*rho_[nx_-1] - phi_[nx_];
		for ( int j = 2; j <= nx_-2 ; j++ ) {
			r_[j] = f_*rho_[j];
		}
		
		solve_tridiagonal_eq();	// solve tridiagonal equation
		
		// get electric field
		for ( int j = 1; j < nx_; j++ ) {
			ex_[j] = -(phi_[j+1] - phi_[j-1])/(2*dx_);
		}
		sig0_ = -(phi_[1] - phi_[0])/dx_*n0_ - rho_[0]*dx_*0.5;
		sig1_ = (phi_[nx_] - phi_[nx_-1])/dx_*n0_ - rho_[nx_]*dx_*0.5;
		ex_[0] = sig0_/n0_;
		ex_[nx_] = -sig1_/n0_;
	//}
	
	/*
	else if ( bc0_ == 0 && bc1_ == 2 ) {	// the right is a conducting wall, the left is symmetric
		// get potential
		phi_[nx_-1] = phi_[nx_] - f_*rho_[nx_];
		for ( int j = nx_-1; j > 1; j-- ) {
			phi_[j-1] = 2*phi_[j] - phi_[j+1] + f_*rho_[j];
		}
		// get electric field
		for ( int j = 1; j < nx_; j++ ) {
			ex_[j] = -(phi_[j+1] - phi_[j-1])/(2*dx_);
		}
		double sigma0 = -(phi_[1] - phi_[0])/dx_*n0_ - rho_[0]*dx_*0.5;
		ex_[0] = sigma0/n0_;
		ex_[nx_] = 0;	// due to symmetry
	}
	else {
		throw runtime_error( "unsupported boundary conditions in Field_1DES" );
	}*/
}

void Field_1DES::solve_tridiagonal_eq()
{
	int n = int(a_.size());
    std::vector<double> gam(n);
	
    if ( b_[1] == 0.0 ) throw std::runtime_error( "error 1 in tridag" );
	
    double bet;
    phi_[1] = r_[1]/(bet = b_[1]);
	for ( int j = 2; j < n; j++ ) {
		gam[j] = c_[j-1]/bet;
		bet = b_[j] - a_[j]*gam[j];
		if ( bet == 0.0 ) throw std::runtime_error( "error 2 in tridag" );
		phi_[j] = (r_[j] - a_[j]*phi_[j-1])/bet;
	}
	
    for ( int j = n-2; j >= 1; j-- )
		phi_[j] -= gam[j+1]*phi_[j+1];
}

void Field_1DES::set_source( int s, double q, double x )
{
    double qdx = q/dx_;
    double posx = (x-x0_)/dx_;
    int j = static_cast<int>( posx );
    
	// bilinear weighting
	double qdx0 = qdx*(j+1 - posx);
	double qdx1 = qdx*(posx - j);
	
	rho_[j] += qdx0;
	rho_[j+1] += qdx1;
	rhok_[s][j] += qdx0;
	rhok_[s][j+1] += qdx1;
}

void Field_1DES::clear_source()
{
    for ( int i = 0; i <= nx_; i++ )
        rho_[i] = 0;
	
	rhok_.fill(0);	// fill all elements with zero
}

double Field_1DES::ex( double x ) const
{
	double posx = (x-x0_)/dx_;
	int j = static_cast<int>( posx );
    
    // bilinear weighting
	return (j+1 - posx)*ex_[j] + (posx - j)*ex_[j+1];
}

void Field_1DES::snapshot( std::ostream& out, bool title ) const
{
    out << "#x\trho\tphi\tex\tbz";
	out << scientific << setprecision(6);
	for ( int j = 0; j < ns_; j++ )
		out << "\trhok[" << j << ']';
	out << '\n';
    for ( int i = 0; i <= nx_; i++ ) {
        out << i*dx_ << '\t'
            << rho_[i] << '\t'
            << phi_[i] << '\t'
            << ex_[i] << '\t'
            << bz_;
		for ( int j = 0; j < ns_; j++ )
			out << '\t' << rhok_[j][i];
		out << '\n';
    }
}

