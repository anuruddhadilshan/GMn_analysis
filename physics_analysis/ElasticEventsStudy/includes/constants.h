#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants
{
	// electron; from particle-data-group (https://pdg.lbl.gov/) on April 18, 2023.
	constexpr double e_chrg = 1.602176634E-19; // exact C
	constexpr double e_mass = 0.51099895000E-3; // +/- 15E-14 GeV

	// neutron; from particle-data-group (https://pdg.lbl.gov/) on April 18, 2023.
	constexpr double n_mass = 0.93956542052; // +/- 54E-11 GeV

	// proton; from particle-data-group (https://pdg.lbl.gov/) on April 18, 2023.
	constexpr double p_mass = 0.93827208816; // +/- 29E-11 GeV

	// Avogadro number
	constexpr double N_A = 6.02214076E23; // mol^-1

	// Molar mass of Hydrogen molecule = 1H_2
	constexpr double H2_u = 2.01588; // g/mol

	// Molar mass of Deuterium molecule = 2D_2
	constexpr double D2_u = 4.0282035556; // g/mol
}

#endif