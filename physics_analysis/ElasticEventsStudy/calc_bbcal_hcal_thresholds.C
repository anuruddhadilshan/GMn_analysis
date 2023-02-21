// A script to find the minimum energy deposited on HCal and PS+SH by elastic events at a given kinematic setting.
// Copied from S.Seed's "GMnElasPeak.C" script and adapted to just spit out the HCal and PS+SH thresholds to be used for good event selection.

#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "beam_variables.h"

const double M_p = 0.938272; // GeV
const double M_n = 0.939565; // GeV
const double M_e = 0.00051; // GeV
const double ff = 0.05; // Added arbitrary smearing factor to account for beam energy fluctuations and fermi motion in downstream estimates
const double hcal_width = 1.70434; // m
const double hcal_sampfrac = 0.0795; // m
const double hcal_threshconv = 6.914; // MeV/mV
const double bbcal_threshconv = 7.2; // MeV/mV

void calc_bbcal_hcal_thresholds(const int runnum)
{
  
  readin_beamvariables(runnum);

  double hcal_minang = hcaltheta - (hcal_width/2)/hcaldist; //approx with arclength
  double hcal_maxang = hcaltheta + (hcal_width/2)/hcaldist; //approx with arclength

  double sh_ypos[7] = {-0.2565, -0.171, -.0855, 0.0, 0.0855, 0.171, 0.2565}; //Relative positions of shower columns.
  double effective_BBang[7] = {0.};
  
  //double deltaBBang = 0.;
  double sh_faceDist = 3.1 + bbdist; //1.2m to GEMs, another 1.9m to BBCal from the BigBite magnet.

  double eprimeEnergy_protonscat[7] = {0.};
  double eprimeEnergy_neutronscat[7] = {0.};
  double nu_p[7] = {0.}; // Keeping around in arrays for possible readout in future versions. Just KE where KE = Ebeam - E_e'
  double nu_n[7] = {0.};

  for (int shcol=0; shcol<7; shcol++)
  {
  	effective_BBang[shcol] = (sh_ypos[shcol]/sh_faceDist) + bbtheta;

  	// If the electron is scattered off the proton in D2.
    eprimeEnergy_protonscat[shcol] = Ebeam/( 1. + (2.*Ebeam/M_p)*pow(sin(effective_BBang[shcol]/2.), 2.) ); // For elastic scattering.
    nu_p[shcol] = Ebeam - eprimeEnergy_protonscat[shcol];

    // If the electron is scattered off the neutron in D2
    eprimeEnergy_neutronscat[shcol] = Ebeam/( 1. + ( 2.*Ebeam/M_n )*pow( sin(effective_BBang[shcol]/2. ), 2.) ); // For elastic scattering.
    nu_n[shcol] = Ebeam - eprimeEnergy_neutronscat[shcol];
  }

  // Check for the lowerst KE i.e the lowest energy deposited in HCal and BBCal(SH+PS)
  double minEnergyinHCal_fromneutron{nu_p[0]};
  double minEnergyinHCal_fromproton{nu_n[0]};
  double minEnergyinBBCal_fromprotonscat{eprimeEnergy_protonscat[0]};
  double minEnergyinBBCal_fromneutronscat{eprimeEnergy_neutronscat[0]};

  for (int shcol=0; shcol<6; shcol++)
  {
    if ( nu_p[shcol] > nu_p[shcol+1] ) minEnergyinHCal_fromproton = nu_p[shcol+1];
    if ( nu_n[shcol] > nu_n[shcol+1] ) minEnergyinHCal_fromneutron = nu_n[shcol+1];
    if ( eprimeEnergy_protonscat[shcol] > eprimeEnergy_protonscat[shcol+1]) minEnergyinBBCal_fromprotonscat = eprimeEnergy_protonscat[shcol+1];
    if ( eprimeEnergy_neutronscat[shcol] > eprimeEnergy_neutronscat[shcol+1]) minEnergyinBBCal_fromneutronscat = eprimeEnergy_neutronscat[shcol+1];
  }

  double minEnergyinHCal{0.};

  if ( minEnergyinHCal_fromproton < minEnergyinHCal_fromneutron )
  {
  	minEnergyinHCal = minEnergyinHCal_fromproton;
  	std::cout << "Lowest energy deposited in HCal = " << minEnergyinHCal << " GeV\n"; 
  }

  else 
  {
  	minEnergyinHCal = minEnergyinHCal_fromneutron;
  	std::cout << "Lowest energy deposited in HCal = " << minEnergyinHCal << " GeV\n";
  }

  double minEnergySampledinHCal{minEnergyinHCal*hcal_sampfrac};

  std::cout << "\nLowest energy sampled by HCal = " << minEnergySampledinHCal << " GeV\n";

  double minEnergySampledinHCalwithsmearing{minEnergySampledinHCal*(1-ff)};

  std::cout << "\nLowest energy sampled in HCal with estimated smearing = " << minEnergySampledinHCalwithsmearing <<" GeV\n";

  std::cout <<"##############################################################################################################\n";

  double minEnergyinBBCal{0.};

  if ( minEnergyinBBCal_fromprotonscat < minEnergyinBBCal_fromneutronscat )
  {
  	minEnergyinBBCal = minEnergyinBBCal_fromprotonscat;
  	std::cout << "\nLowerst energy deposted in BBCal = " << minEnergyinBBCal << " GeV\n";
  }

  else 
  {
  	minEnergyinBBCal = minEnergyinBBCal_fromneutronscat;
  	std::cout << "\nLowerst energy deposited in BBCal = " << minEnergyinBBCal << " GeV\n";
  }

  double minEnergyinBBCal_withsmearing{minEnergyinBBCal*(1-ff)};

  std::cout << "\nLowest energy in BBCal with estimated smearing = " << minEnergyinBBCal_withsmearing << " GeV\n";

}