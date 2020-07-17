//
//  Particle_1DES.h
//  es1d2v
//
//  Created by Jun Hasegawa on 2013/03/05.
//  Last Updated by Jun Hasegawa on 2013/06/27.
//
//  Copyright (c) 2013 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#ifndef Particle_1DES_h
#define Particle_1DES_h

class Particle_1DES {
public:
    Particle_1DES( long id, double m, double q, double x, double vx, double vy )
    : id_(id), m_(m), q_(q), x_(x), vx_(vx), vy_(vy) {}
	
    virtual void move( double dt );
    // dt : time step
    
    virtual void accelerate( double dt, double ex, double bz );
    // dt : time step
    // ex : electric field
    // bz : magnetic field
    
    // member functions to access private members
    double id() const { return id_; }
    double m() const { return m_; }
    double q() const { return q_; }
    double x() const { return x_; }
    double vx() const { return vx_; }
    double vy() const { return vy_; }
	
	void set_m( double m ) { m_ = m; }
	void set_q( double q ) { q_ = q; }
	void set_x( double x ) { x_ = x; }
	void set_vx( double vx ) { vx_ = vx; }
	void set_vy( double vy ) { vy_ = vy; }

private:
    long id_;  // particle ID
    double m_;  // mass
    double q_;  // charge
    double x_;  // position
    double vx_; // velocity in x
    double vy_; // velocity in y
};

#endif
