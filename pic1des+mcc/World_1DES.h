//
//  World_1DES.h
//  es1d2v
//
//  Created by Jun on 2013/08/01.
//	Last Updated by Jun Hasegawa on 2019/07/12
//
//  Copyright (c) 2013-2019 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#ifndef World_1DES_h
#define World_1DES_h

#include <vector>
#include <string>

class Particles_1DES;
//class Field_1DES;

#include "Field_1DES.h"

class World_1DES {
public:
	World_1DES( const char* fn );
	~World_1DES();
	
	void run();
	
	// member functions for direct access private members
	double dt() const { return dt_; }
	long tstep() const { return tstep_; }
	long tout() const { return tout_; }
	int ns() const { return ns_; }
	int nr() const { return nr_; }
	int nd() const { return nd_; }
	
	double lx() const { return field_[0]->lx(); }
    double phin() const { return field_[0]->phin(); }

private:
	double xd_;	// total size of simulation space (unused)
	double dt_;	// step time
    long tstep_;	// total number of time steps
	long tout_;	// interval for file output
	
	double ng_;	// background gas density
		
	int ns_;	// number of particle species
	int nr_;	// number of field regions
	int nd_;	// number of detectors
	
	double xmin_;
	double xmax_;
	
	// containers for Class pointers
	std::vector<Particles_1DES*> particles_;
	std::vector<Field_1DES*> field_;
};

#endif
