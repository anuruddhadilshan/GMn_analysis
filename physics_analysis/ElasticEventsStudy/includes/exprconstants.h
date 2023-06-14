#ifndef EXPRCONSTANTS_H
#define	EXPRCONSTANTS_H

#include "constants.h"

namespace Exprconstants
{
	// Target constants: fixed for the entire experiment. //
	constexpr double targ_length = 15.0; //cm

	// LH2
	constexpr double lH2_density = 0.0723; // g/cc - liquid hydrogen density (operation conditions: 19 K, 25 psi)

	// LD2
	constexpr double lD2_density = 0.167; // g/cc - liquid deuterium density (operation conditions: 22 K, 22 psi)

	// Luminosity calculations. //
	// Luminosity = (Particle rate of the beam).("Atomic density" of the target).(Length of the target cell) => units s^-1.cm^-2
	// For simplicity, we can calculate a constant factor for both LD2 and LH2, 
	// where when multiplied by the beam current in microampheres, we get the luminosity in s^-1.cm^-2.
	constexpr double lH2_lumifactor = (1.0E-6/Constants::e_chrg)*(2*(lH2_density/Constants::H2_u)*Constants::N_A)*targ_length; // cm^-2
	constexpr double lD2_lumifactor = (1.0E-6/Constants::e_chrg)*(2*(lD2_density/Constants::D2_u)*Constants::N_A)*targ_length; // cm^-2
}

#endif