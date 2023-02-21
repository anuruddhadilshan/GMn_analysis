// Script to define tight cuts to filter out good/elastic events
// E/P cut, trigger coincidence cut, shower and preshower energy...

#ifndef DEFINE_TIGHTCUTS_H
#define DEFINE_TIGHTCUTS_H

#include <TMath.h>
#include <iomanip>
#include "beam_variables.h"

//All the cuts are defined as global variables so that they can be accessed directly from the analysis script.

//// Vertex Z cut ////
// Manually define a tight vertex cut to only accept the events originating from along the length of the target.
double vzcut{0.075};

void vertexZ_cut(TH1D* h1_bbtrackvertz)
{
	TCanvas* c1 = new TCanvas();
	c1->cd();
	double bbtrackvertz_histmax{h1_bbtrackvertz->GetMaximum()};
	TLine* leftvzcut = new TLine(-vzcut,0,-vzcut,bbtrackvertz_histmax);
	TLine* rightvzcut = new TLine(vzcut,0,vzcut,bbtrackvertz_histmax);
	leftvzcut->SetLineColorAlpha(kBlue,0);
	leftvzcut->SetLineWidth(2);
	leftvzcut->SetLineStyle(9);
	rightvzcut->SetLineColorAlpha(kBlue,0);
	rightvzcut->SetLineWidth(2);
	rightvzcut->SetLineStyle(9);
	h1_bbtrackvertz->Draw();
	leftvzcut->Draw("same");
	rightvzcut->Draw("same");
	std::cout << "\ndefine_tightcuts()->vertexZ_cut() = " << leftvzcut << " < vertex_z < " << rightvzcut <<'\n'; 
}

//// BBCal and HCal coin cut left and right boundaries ////
// This cut will help to remove accidental events that falls within the trigger window.
double trigdiff_leftcut{0.}; 
double trigdiff_rightcut{0.};
	
void bbCal_HCal_coincut(TH1D* h1_bbcal_hcal_tdiff)
{
	TCanvas* c2 = new TCanvas();
	c2->cd();	
	h1_bbcal_hcal_tdiff->Fit("gaus");
	TF1* fitfunc_trigdiff = h1_bbcal_hcal_tdiff->GetFunction("gaus");
	double trigdiff_mean{fitfunc_trigdiff->GetParameter(1)};
	double trigdiff_stddev{fitfunc_trigdiff->GetParameter(2)};
	trigdiff_leftcut = trigdiff_mean-trigdiff_stddev*2;
	trigdiff_rightcut = trigdiff_mean+trigdiff_stddev*2;
	double trigdiff_histmax{h1_bbcal_hcal_tdiff->GetMaximum()};
	//std::cout << fitfunc_trigdiff->GetParameter(0) <<" "<< fitfunc_trigdiff->GetParameter(1) <<" " << fitfunc_trigdiff->GetParameter(2) <<'\n';
	TLine* ll_trigdiff = new TLine(trigdiff_leftcut,0,trigdiff_leftcut,trigdiff_histmax);
	TLine* lr_trigdiff = new TLine(trigdiff_rightcut,0,trigdiff_rightcut,trigdiff_histmax);
	ll_trigdiff->SetLineColorAlpha(kBlue,0);
	ll_trigdiff->SetLineWidth(2);
	ll_trigdiff->SetLineStyle(9);
	lr_trigdiff->SetLineColorAlpha(kBlue,0);
	lr_trigdiff->SetLineWidth(2);
	lr_trigdiff->SetLineStyle(9);
	h1_bbcal_hcal_tdiff->Draw();
	ll_trigdiff->Draw("same");
	lr_trigdiff->Draw("same");
	std::cout << "\ndefine_tightcuts()->bbCal_HCal_coincut() = " << ll_trigdiff << " < HCal_time - BBCal_time < " << lr_trigdiff << '\n';
}

//// E over P cut left and right boundaries ////
double eoverp_leftcut{0.};
double eoverp_rightcut{0.};

void e_over_p_cut(TH1D* h1_EoverP)
{
	TCanvas* c3 = new TCanvas();
	c3->cd();
	h1_EoverP->Fit("gaus");
	TF1* fitfunc_eoverp = h1_EoverP->GetFunction("gaus");
	double eoverp_mean{fitfunc_eoverp->GetParameter(1)};
	double eoverp_stddev{fitfunc_eoverp->GetParameter(2)};
	eoverp_leftcut = eoverp_mean - eoverp_stddev*3;
	eoverp_rightcut = eoverp_mean + eoverp_stddev*3;
	double eoverp_histmax{h1_EoverP->GetMaximum()};
	TLine* ll_eoverp = new TLine(eoverp_leftcut,0,eoverp_leftcut,eoverp_histmax);
	TLine* lr_eoverp = new TLine(eoverp_rightcut,0,eoverp_rightcut,eoverp_histmax);
	ll_eoverp->SetLineColorAlpha(kBlue,0);
	ll_eoverp->SetLineWidth(2);
	ll_eoverp->SetLineStyle(9);
	lr_eoverp->SetLineColorAlpha(kBlue,0);
	lr_eoverp->SetLineWidth(2);
	lr_eoverp->SetLineStyle(9);
	h1_EoverP->Draw();
	ll_eoverp->Draw("same");
	lr_eoverp->Draw("same");	
	std::cout << "\ndefine_tightcuts()->e_over_p_cut() = " << ll_eoverp << " < E/P < " << lr_eoverp << '\n';
}

//// HCal and BBCal energy cut ////
// The energy deposited on BBCal and HCal by elastically scatterred electrons and hadrons will be simulated for a given kinematic setting.
// Then a cut on HCal and BBCal (SH+PS) energy will be applied.
double hcal_energycut{0.};
double bbcal_energycut{0.};

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

  //readin_beamvariables(runnum);	
  double hcal_minang = hcaltheta - (hcal_width/2)/hcaldist; //approx with arclength
  double hcal_maxang = hcaltheta + (hcal_width/2)/hcaldist; //approx with arclength

  double sh_ypos[7] = {-0.2565, -0.171, -.0855, 0.0, 0.0855, 0.171, 0.2565}; //Relative positions of shower columns.
  double effective_BBang[7] = {0.};
  
  double sh_faceDist = 3.1 + bbdist; //1.2m to GEMs, another 1.9m to BBCal from the BigBite magnet.

  double eprimeEnergy_protonscat[7] = {0.};
  double eprimeEnergy_neutronscat[7] = {0.};
  double nu_p[7] = {0.}; // Just KE where KE = Ebeam - E_e'
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

  //std::cout << "\nLowest energy sampled by HCal = " << minEnergySampledinHCal << " GeV\n";

  double minEnergySampledinHCalwithsmearing{minEnergySampledinHCal*(1-ff)};
  hcal_energycut = minEnergySampledinHCalwithsmearing;

  //std::cout << "\nSimulated lowest energy (from elastic scattering processes) sampled in HCal with estimated smearing of " << ff*100 <<"% = " << minEnergySampledinHCalwithsmearing <<" GeV\n";
  std::cout << "\ndefine_tightcuts()->calc_bbcal_hcal_thresholds(): Simulated lowest energy (from elastic scattering processes) sampled in HCal with estimated smearing of " << ff*100 <<"% = " << minEnergySampledinHCalwithsmearing <<" GeV\n";

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
  bbcal_energycut = minEnergyinBBCal_withsmearing;

  //std::cout << "\nSimulated lowest energy in BBCal (from elastic scattering processes) with estimated smearing of " << ff*100 << "% = " << minEnergyinBBCal_withsmearing << " GeV\n";
  std::cout << "\ndefine_tightcuts()->calc_bbcal_hcal_thresholds(): Simulated lowest energy in BBCal (from elastic scattering processes) with estimated smearing of " << ff*100 << "% = " << minEnergyinBBCal_withsmearing << " GeV\n";
}

void hcal_energy_cut(TH1D* h1_HCal_e)
{
	TCanvas* c4 = new TCanvas();
	c4->cd();
	double hcale_histmax{h1_HCal_e->GetMaximum()};
	//double min_hcal_e{0.10989167};
	TLine* ll_minhcale = new TLine(hcal_energycut,0,hcal_energycut,hcale_histmax);
	ll_minhcale->SetLineColorAlpha(kBlue,0);
	ll_minhcale->SetLineWidth(2);
	ll_minhcale->SetLineStyle(9);
	h1_HCal_e->Draw();
	ll_minhcale->Draw("same");
}

void ps_plus_sh_energy_cut(TH1D* h1_sh_ps_e)
{
	TCanvas* c5 = new TCanvas();
	c5->cd();
	double shpse_histmax{h1_sh_ps_e->GetMaximum()};
	//double min_shpse{1.7};
	TLine* ll_minshpse = new TLine(bbcal_energycut,0,bbcal_energycut,shpse_histmax);
	ll_minshpse->SetLineColorAlpha(kBlue,0);
	ll_minshpse->SetLineWidth(2);
	ll_minshpse->SetLineStyle(9);
	h1_sh_ps_e->Draw();
	ll_minshpse->Draw("same");
}

void define_tightcuts(const int runnum, TH1D* h1_bbtrackvertz, TH1D* h1_bbcal_hcal_tdiff, TH1D* h1_EoverP, TH1D* h1_HCal_e, TH1D* h1_sh_ps_e)
{
	
	gStyle->SetOptFit(1);

	/*vertexZ_cut(h1_bbtrackvertz);
	
	//bbCal_HCal_coincut(h1_bbcal_hcal_tdiff);	
	
	e_over_p_cut(h1_EoverP);
	
	calc_bbcal_hcal_thresholds(runnum);

	hcal_energy_cut(h1_HCal_e);
	
	ps_plus_sh_energy_cut(h1_sh_ps_e);*/

	std::cout << "\nFunction define_tightcuts() finished.\n";
	
}

#endif