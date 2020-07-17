//
//  Field_1DES.cpp
//  es1d2v
//
//  Created by Jun Hasegawa on 2013/03/06.
//  Last Updated by Jun Hasegawa on 2019/07/12.
//
//  Copyright (c) 2019 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#ifndef Field_1DES_h
#define Field_1DES_h

#include "Matrix.h"

#include <iostream>
#include <fstream>
#include <vector>

class Field_1DES {

public:
    Field_1DES( int nx, double lx, double x0, double n0, double phi_n, double bz = 0, int bc0 = 0, int bc1 = 0 );
    // nx: number of cells
    // lx: field size
    // x0: left boundary
    // n0: number density of super particles
    // phi_n: potential of right boundary with respect to left boundary
    // bz: magnetic field
    // bc0: boundary condition for left boundary
    // bc1: boundary condition for right boundary
    
    virtual void solve();
    virtual void clear_source();
    virtual void set_source( int s, double q, double x );
    
    virtual double ex( double x ) const;
    virtual double bz() const { return bz_; }
    
    virtual bool out_of_region( double x ) const { return x < x0_ || x > x1_; }

    virtual void set_bz( double bz ) { bz_ = bz; }
    virtual void set_phi0( double phi_0 ) { phi_.at(0) = phi_0; }
    virtual void set_phin( double phi_n ) { phi_.at(nx_) = phi_n; }
    
    virtual void snapshot( std::ostream& out, bool title = true ) const;
    
    // member functions for direct access to private members
    long nx() const { return nx_; }
    double x0() const { return x0_; }
    double lx() const { return lx_; }
    double dx() const { return dx_; }
    double n0() const { return n0_; }
    double bc0() const { return bc0_; }
    double bc1() const { return bc1_; }
    double phi0() const { return phi_.at(0); }
    double phin() const { return phi_.at(nx_); }
    
    void set_ns( int ns ) { ns_ = ns; }
    virtual void add_particles() { ns_++; }
    
private:
    void solve_tridiagonal_eq();
    
    long nx_;   // number of cells
    double lx_; // field size
    double x0_; // position of left boundary
    double x1_; // position of right boundary
    double dx_; // grid spacing
    double n0_; // number density of super particles
    double bz_; // background magnetic field
    int bc0_;   // boundary condition for left boundary
    int bc1_;   // boundary condition for right boundary
    
    int ns_;    // number of particle groups
    double	sig0_;	// surface charge density at left boundary
	double	sig1_;	// surface charge density at right boundary
	
    std::vector<double> rho_;   // charge density
    std::vector<double> phi_;   // scalar potential
    std::vector<double> ex_;    // electric field in x
    Matrix<double> rhok_;	// charge density of each species
    
    // coefficients for matrix solver
    std::vector<double> a_, b_, c_, r_;
    double f_;
};

#endif
