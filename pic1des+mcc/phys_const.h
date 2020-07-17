#ifndef phys_const_h
#define phys_const_h

namespace si{
	const double c   = 2.9979246e8;		// light speed [m/s]
	const double e   = 1.6021773e-19;	// charge of electron [C]
	const double mu0 = 1.2566371e-6;	// magnetic permeability for vacuum [H/m]
	const double ep0 = 8.8541878e-12;	// dielectric constant for vacuum [F/m]
	const double h   = 6.6260755e-34;	// Planck constant [Js]
	const double h2  = 1.05457266e-34;	// h/2pi [Js]
	const double me  = 9.1093897e-31;	// mass of electron [kg]
	const double mp  = 1.6726231e-27;	// mass of proton [kg]
	const double mn  = 1.6749286e-27;	// mass of neutron [kg] 
	const double amu = 1.6605402e-27;	// atomic mass unit [kg]
	const double na  = 6.0221367e23;	// Avogadro number [1/mol]
	const double kb  = 1.380658e-23;	// Boltzmann constant [J/K]
	const double sg  = 5.67051e-8;		// Stefan-Boltzmann coefficient [W/m^2/K^4] 
	const double a0  = 5.2917725e-11;	// Bohr radius [m] 
	const double vb  = 2.1876914e6;		// Bohr velocity [m/s] 
	const double pi  = 3.1415926535;	// circumference
	const double af  = 7.29735308e-3;	// fine structure constant
}

namespace cgs{
	const double c   = 2.9979246e10;	// light speed [cm/s]
	const double e   = 4.80325e-10;		// charge of electron [esu]
	const double mu0 = 1;				// magnetic permeability for vacuum [H/m]
	const double ep0 = 1;				// dielectric constant for vacuum [F/m]
	const double h   = 6.6260755e-27;	// Planck constant [erg s]
	const double h2  = 1.05457266e-27;	// h/2pi [erg s]
	const double me  = 9.1093897e-28;	// mass of electron [g]
	const double mp  = 1.6726231e-24;	// mass of proton [g]
	const double mn  = 1.6749286e-24;	// mass of neutron [g] 
	const double amu = 1.6605402e-24;	// atomic mass unit [g]
	const double na  = 6.0221367e23;	// Avogadro number [1/mol]
	const double kb  = 1.380658e-16;	// Boltzmann constant [erg/K]
	const double a0  = 5.2917725e-9;	// Bohr radius [cm] 
	const double vb  = 2.1876914e8;		// Bohr velocity [cm/s] 
	const double pi  = 3.1415926535;	// circumference
	const double af  = 7.29735308e-3;	// fine structure constant
}

#endif

