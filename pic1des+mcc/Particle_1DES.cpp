//
//  PICParticle_1DES.cpp
//  es1d2v
//
//  Created by Jun Hasegawa on 2012/11/06.
//  Last Updated by Jun Hasegawa on 2013/03/06.
//
//  Copyright (c) 2013 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#include "Particle_1DES.h"

#include <cmath>

void Particle_1DES::move( double dt )
{
    x_ += vx_*dt;
}

void Particle_1DES::accelerate( double dt, double ex, double bz )
{
    double qm = q_/m_;
    double wcdt = qm*bz*dt;
    double sin_wcdt = std::sin( wcdt );
    double cos_wcdt = std::cos( wcdt );
    
    // half acceleration
    double dv = qm*ex*(0.5*dt);
    vx_ += dv;
    
    // v x B rotation
    double vxx =  cos_wcdt * vx_ + sin_wcdt * vy_;	// rotation
    double vyy = -sin_wcdt * vx_ + cos_wcdt * vy_;
    vx_ = vxx;
    vy_ = vyy;
    
    // half acceleration
    vx_ += dv;
}

