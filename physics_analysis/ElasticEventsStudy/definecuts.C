// Written by Anuruddha Rathnayake 
// Dec 16, 2022.

// 1) This script analyses the replayed data from GMn/nTPE and makes some basic histograms.
// 2) Namely: W^2, vertexZ, E/P, PS energy, SH+PS energy, HCal energy, HCal time - BBCal time, HCal block ADC time in main cluster.
// 3) Does several fittings / calculations and find out the optimal (possibly a lot of room to improve) cut thresholds for different variables from the histograms.
// 4) Make a root file of the historams and also make a PDF with all the histograms with cut thresholds marked on them (the cut thresholds are still not being marked on the PDFs due to some issue with the Print() option when making histograms.s

#include <iostream>
#include <string>
#include <fstream>
#include "TChain.h"
#include "TFile.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "beam_variables.h"
#include "calc_q_W2.h"

const double target_mass{0.5*(0.938272+0.939565)}; //Average of neutron and proton mass.
constexpr int n_histos{9}; // The number of histograms that needs to be printed into the PDF.

void readin_definecutsconfigfile(const char* configfilename, TChain* C)
{

  std::cout <<'\n'<<"--- Reading configuration file " << configfilename << " --- \n";
    	
  ifstream configfile(configfilename);
	TString currentline;

	while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") )
	{
  		if( !currentline.BeginsWith("#") )
   		{
   			C->Add(currentline);
   			cout << "Loaded root file: " << currentline << '\n';
   		}
  }

}

//Function to apply global cuts and print out some iformation into the screen.
/*{
	nevents_TChain = C->GetEntries();
	std::cout<<"\nNumber of events in the TChain = " << nevents_TChain <<'\n';
	std::cout<<"--- Begin applying global cuts = " << globalcut << " ---" <<'\n';
	C->Draw(">>elist",globalcut);
	std::cout<<"--- Done applying global cuts ---" <<'\n';
	nevents_eventlist = elist->GetN();	
	std::cout<<"Number of events in the eventlist = " << nevents_eventlist <<'\n';
	std::cout<<"Percentage of events accepted for physics analysis = " << (nevents_eventlist/(double)nevents_TChain)*100 <<"%\n";
}*/

void get_TDC_times(int ndata_tdc, double tdc_elemID[1000], double tdc_time[1000], double& bbcal_time, double& hcal_time, double& coin_time, double&rf_time)
{
	for(int ihit=0; ihit<ndata_tdc; ihit++)
	{
      if(tdc_elemID[ihit]==5) bbcal_time=tdc_time[ihit];
      if(tdc_elemID[ihit]==0) hcal_time=tdc_time[ihit];
      if(tdc_elemID[ihit]==1) coin_time=tdc_time[ihit];
      if(tdc_elemID[ihit]==4) rf_time=tdc_time[ihit];
    }
}

void print_analysis_percentage(double currentevent_ana_percentage, int& previousevent_ana_percentage_int)
{
	int currentevent_ana_percentage_int{(int)currentevent_ana_percentage};
	if (currentevent_ana_percentage_int%10==0 && currentevent_ana_percentage_int>previousevent_ana_percentage_int)
	{
		std::cout << currentevent_ana_percentage_int <<"%\n";
	}

	previousevent_ana_percentage_int = currentevent_ana_percentage_int;
}

//// W^2 cut. Manually defined. ////
double W2_min{0.48};
double W2_max{1.28};

void w2_cut(TH1D* h1_W2, TCanvas* c)
{
	c->cd();
	double w2_histmax{h1_W2->GetMaximum()};
	TLine* leftw2cut = new TLine(W2_min,0,W2_min,w2_histmax);
	TLine* rightw2cut = new TLine(W2_max,0,W2_max,w2_histmax);
	leftw2cut->SetLineColorAlpha(kBlue,0);
	leftw2cut->SetLineWidth(2);
	leftw2cut->SetLineStyle(9);
	rightw2cut->SetLineColorAlpha(kBlue,0);
	rightw2cut->SetLineWidth(2);
	rightw2cut->SetLineStyle(9);
	h1_W2->Draw();
	leftw2cut->Draw();
	rightw2cut->Draw();
}

//// preshower cut. Manually define. ////
double pse_min{0.2};

void preshower_cut(TH1D* h1_bb_ps_e, TCanvas* c)
{
	c->cd();
	double pse_histmax(h1_bb_ps_e->GetMaximum());
	TLine* min_pse_line = new TLine(pse_min,0,pse_min,pse_histmax);
	min_pse_line->SetLineColorAlpha(kBlue,0);
	min_pse_line->SetLineWidth(2);
	min_pse_line->SetLineStyle(9);
	h1_bb_ps_e->Draw();
	min_pse_line->Draw();
}

//// Vertex Z cut ////
// Manually define a tight vertex cut to only accept the events originating from along the length of the target.
double vzcut{0.075};

void vertexZ_cut(TH1D* h1_bbtrackvertz, TCanvas* c)
{
	c->cd();
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
	leftvzcut->Draw();
	rightvzcut->Draw();
}

//// E over P cut left and right boundaries ////
double eoverp_leftcut{0.};
double eoverp_rightcut{0.};

void e_over_p_cut(TH1D* h1_EoverP, TCanvas* c)
{
	c->cd();
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
	ll_eoverp->Draw();
	lr_eoverp->Draw();	
}

//// BBCal and HCal coin cut left and right boundaries ////
// This cut will help to remove accidental events that falls within the trigger window.
double trigdiff_leftcut{0.}; 
double trigdiff_rightcut{0.};
	
void bbCal_HCal_coincut(TH1D* h1_bbcal_hcal_tdiff, TCanvas* c)
{
	c->cd();	
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
	ll_trigdiff->Draw();
	lr_trigdiff->Draw();
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
//const double hcal_sampfrac = 0.0795; // m -> From MC simulations.
const double hcal_sampfrac = 0.04; // From S.Seeds -> Looking at real data. (0.049 is the value he gave but putting in 0.04 with a bit of a "safety margin").
const double hcal_threshconv = 6.914; // MeV/mV
const double bbcal_threshconv = 7.2; // MeV/mV

void calc_bbcal_hcal_thresholds(const int kine_num, TH1D* h1_HCal_e, TH1D* h1_sh_ps_e, TCanvas* c1, TCanvas* c2)
{

  lookup_kinematic_info(kine_num);	
  
  double hcal_minang = hcaltheta - (hcal_width/2)/hcaldist; //approx with arclength
  double hcal_maxang = hcaltheta + (hcal_width/2)/hcaldist; //approx with arclength

  double sh_ypos[7] = {-0.2565, -0.171, -.0855, 0.0, 0.0855, 0.171, 0.2565}; //Relative positions of shower columns.
  double effective_BBang[7] = {0.};
  
  //double deltaBBang = 0.;
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
  	std::cout << "\nLowest energy deposited in HCal = " << minEnergyinHCal << " GeV\n"; 
  }
  else 
  {
  	minEnergyinHCal = minEnergyinHCal_fromneutron;
  	std::cout << "\nLowest energy deposited in HCal = " << minEnergyinHCal << " GeV\n";
  }

  double minEnergySampledinHCal{minEnergyinHCal*hcal_sampfrac};

  //std::cout << "\nLowest energy sampled by HCal = " << minEnergySampledinHCal << " GeV\n";

  double minEnergySampledinHCalwithsmearing{minEnergySampledinHCal*(1-ff)};
  hcal_energycut = minEnergySampledinHCalwithsmearing;

  std::cout << "\nSimulated lowest energy (from elastic scattering processes) sampled in HCal with estimated smearing of " << ff*100 <<"% = " << minEnergySampledinHCalwithsmearing <<" GeV\n";

  c1->cd();
  double hcale_histmax{h1_HCal_e->GetMaximum()};
  TLine* ll_minhcale = new TLine(hcal_energycut,0,hcal_energycut,hcale_histmax);
  ll_minhcale->SetLineColorAlpha(kBlue,0);
  ll_minhcale->SetLineWidth(2);
  ll_minhcale->SetLineStyle(9);
  h1_HCal_e->Draw();
  ll_minhcale->Draw();

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

  std::cout << "\nLowest energy in BBCal (from elastic scattering processes) with estimated smearing of " << ff*100 << "% = " << minEnergyinBBCal_withsmearing << " GeV\n";

  c2->cd();
  double shpse_histmax{h1_sh_ps_e->GetMaximum()};
  TLine* ll_minshpse = new TLine(bbcal_energycut,0,bbcal_energycut,shpse_histmax);
  ll_minshpse->SetLineColorAlpha(kBlue,0);
  ll_minshpse->SetLineWidth(2);
  ll_minshpse->SetLineStyle(9);
  h1_sh_ps_e->Draw();
  ll_minshpse->Draw();  

}


double hcal_clusblk_ADCtime_leftcut{0.};
double hcal_clusblk_ADCtime_rightcut{0.};

void hcal_clusblk_ADCtime_cut(TH1D* h1_hcal_clusblk_ADCtime, TCanvas* C)
{
	C->cd();
	TF1* gaus_fit = new TF1("gaus_fit", "gaus", 20, 50); // <<<<< ADJUST THE RANGE HERE TO GET THE FIT RIGHT.
	h1_hcal_clusblk_ADCtime->Fit(gaus_fit,"R");
	double mean{gaus_fit->GetParameter(1)};
	double sigma{gaus_fit->GetParameter(2)};
	double hist_max{h1_hcal_clusblk_ADCtime->GetMaximum()};
	double ll_xpos{mean - 3*sigma};
	double rl_xpos{mean + 3*sigma};
	hcal_clusblk_ADCtime_leftcut = ll_xpos;
	hcal_clusblk_ADCtime_rightcut = rl_xpos;
	TLine* ll = new TLine(ll_xpos,0,ll_xpos,hist_max);
	TLine* rl = new TLine(rl_xpos,0,rl_xpos,hist_max);
	ll->SetLineColorAlpha(kBlue,0);
  ll->SetLineWidth(2);
  ll->SetLineStyle(9);
  rl->SetLineColorAlpha(kBlue,0);
  rl->SetLineWidth(2);
  rl->SetLineStyle(9);
 
	h1_hcal_clusblk_ADCtime->Draw();
  ll->Draw();
  rl->Draw();
}

void make_pdf(TCanvas* c[n_histos],TString output_file_name)
{
	TString pdffilename = output_file_name;
	TString openfilename = pdffilename+"(";
	TString closefilename = pdffilename+")";

	double lmargin=0.15;
  	double rmargin=0.15;
    double bmargin=0.15;
    double tmargin=0.09;

    for (int icanvas = 0; icanvas < n_histos; icanvas++)
	{
		if(icanvas == 0) c[icanvas]->Print(openfilename);
		else	if (icanvas == n_histos-1) c[icanvas]->Print(closefilename);
		else c[icanvas]->Print(pdffilename);
	}
}

void make_configfile(TString outputconfigfilename)
{
	ofstream configfile;
	configfile.open(outputconfigfilename,ios::out);

	configfile << "W2_min " << W2_min <<'\n';
	configfile << "W2_max " << W2_max <<'\n';
	configfile << "vzcut " << vzcut <<'\n';
	configfile << "pse_min " << pse_min <<'\n';
	configfile << "eoverp_leftcut " << eoverp_leftcut <<'\n';
	configfile << "eoverp_rightcut " << eoverp_rightcut <<'\n';
	configfile << "hcal_energycut " << hcal_energycut <<'\n';
	configfile << "bbcal_energycut " << bbcal_energycut <<'\n';
	configfile << "trigdiff_leftcut " << trigdiff_leftcut <<'\n';
	configfile << "trigdiff_rightcut " << trigdiff_rightcut <<'\n';
	configfile << "hcal_clusblk_ADCtime_leftcut " << hcal_clusblk_ADCtime_leftcut <<'\n';
	configfile << "hcal_clusblk_ADCtime_rightcut " << hcal_clusblk_ADCtime_rightcut <<'\n';

	configfile.close();
}


//// Start Main Function ////

void definecuts(const int kine_num, const char* configfilename = "setup_definecuts.cfg", TString outputfilename  = "setup_analysis")
{

	lookup_kinematic_info(kine_num);

	TChain* C = new TChain("T"); //Initialize the TChain to chain the root files for analysis.
	
	readin_definecutsconfigfile(configfilename,C); // TChain and all the cut thresholds are defined inside the "readconfigfile.h".

	//Defining event list.
	/*TEventList* elist = new TEventList("elist","");
	long nevents_TChain{C->GetEntries()};
	long nevents_eventlist{0};
	apply_globalcuts(C,elist,nevents_TChain,nevents_eventlist); // Make the event list passing the global cuts out of the TChain.*/

	std::cout<<"\n--- Beginning making histograms for cut definition ---" <<'\n';

	const int MAXNTRACKS{10};
	const int MAXNTDC{1000};

	//variables needed are BigBite track px,py,pz and sbs hcal x,y,e
 	double ntrack{0.};
 	double vz[MAXNTRACKS];
	double epx[MAXNTRACKS];
	double epy[MAXNTRACKS];
	double epz[MAXNTRACKS];
	double ep[MAXNTRACKS];
	double sh_e[MAXNTRACKS];
	double ps_e[MAXNTRACKS];
	double hcal_e {0.};
	double bbsh_e {0.};
	double bbps_e {0.};
	double tdc_time[MAXNTDC];
	double tdc_elemID[MAXNTDC];
	int ndata_tdc{0};
	double hcal_clusblk_ADC_time[15]; // Maximum number of blocks in a cluster is 15 as per S.Seeds.
	

	C->SetBranchStatus("*",0);
	C->SetBranchStatus("bb.tr.n",1);
	C->SetBranchStatus("bb.tr.vz",1);
	C->SetBranchStatus("bb.tr.px",1);
	C->SetBranchStatus("bb.tr.py",1);
	C->SetBranchStatus("bb.tr.pz",1);
	C->SetBranchStatus("bb.tr.p",1);
	C->SetBranchStatus("bb.sh.e",1);
	C->SetBranchStatus("bb.ps.e",1);
	C->SetBranchStatus("sbs.hcal.e",1);
	C->SetBranchStatus("bb.tdctrig.tdc",1);
	C->SetBranchStatus("bb.tdctrig.tdcelemID",1);
	C->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
	C->SetBranchStatus("sbs.hcal.clus_blk.atime",1);

	C->SetBranchAddress("bb.tr.n",&ntrack);
	C->SetBranchAddress("bb.tr.vz",vz);
	C->SetBranchAddress("bb.tr.px",epx);
	C->SetBranchAddress("bb.tr.py",epy);
	C->SetBranchAddress("bb.tr.pz",epz);
	C->SetBranchAddress("bb.tr.p",ep);
	C->SetBranchAddress("bb.sh.e",&bbsh_e);
	C->SetBranchAddress("bb.ps.e",&bbps_e);
	C->SetBranchAddress("sbs.hcal.e",&hcal_e);
	C->SetBranchAddress("bb.tdctrig.tdc",tdc_time);
	C->SetBranchAddress("bb.tdctrig.tdcelemID",tdc_elemID);
	C->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID",&ndata_tdc);
	C->SetBranchAddress("sbs.hcal.clus_blk.atime",hcal_clusblk_ADC_time);

	TString outputrootfilename = outputfilename+".root";
	TFile* fout = new TFile(outputrootfilename,"RECREATE");

	TH1D* h1_W2 = new TH1D("h1_W2","W^{2} of all the events passing the global cut; W^2 (GeV^2/c^4); Entries",250,-1,4);
	TH1D* h1_bb_trackvertz = new TH1D("h1_bb_trackvertz","BB Track vertex Z position; vertex Z (m); Entries",600,-0.15,0.15);
	TH1D* h1_eprime_EoverP = new TH1D("h1_eprime_EoverP","Scattered electron E/P; E/P; Entries",1000,-1,3);
	TH1D* h1_bb_ps_e = new TH1D("h1_bb_ps_e","Energy Deopsited in BB Pre Shower; Energy (GeV); Entries",500,0,5);
	TH1D* h1_bb_sh_plus_ps_e = new TH1D("h1_bb_sh_plus_ps_e","Total SH and PS cluster energy sum; SH+PS Energy (GeV)' Entries",500,0,5);
	TH2D* h2_bb_sh_ps_energycorrelation = new TH2D("h2_bb_sh_ps_energycorrelation","Energy deposited in BB PS vs SH;SH Energy (GeV);PS Energy (GeV)",500,0,5,500,0,5);
	TH1D* h1_hcal_e = new TH1D("h1_hcal_e","HCal Energy Deopsited; HCal Energy (GeV); Entries",250,0,0.5);
	TH1D* h1_bbcal_hcal_tdiff = new TH1D("h1_bbcal_hcal_tdiff","HCal time - BBCal time; HCal_{time}-BBCal_{time} (ns); Entries",300,400,700);
	TH1D* h1_hcal_clusblk_ADCtime = new TH1D("h1_hcal_clusblk_ADCtime","ADC time of the highest energy block in the largest cluster; ADC Time (ns)",300,-100,200);
	//Defining the constant 4-vecors of the incoming beam electrons and the target 
	TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
	TLorentzVector Ptarg(0,0,0,target_mass);

	//Variables to keep track of printout to the terminal to keep track of analysis progress.
	int previousevent_ana_percentage_int{0};

	long nevents_TChain{C->GetEntries()};
	long nevent {0};

	//while (C->GetEntry(elist->GetEntry(nevent++)))
	while (C->GetEntry(nevent++))
	{
		
		//double ana_percentage{(nevent/(double)nevents_eventlist)*100}; //Percentage of events analyzed in the Event List.
		double ana_percentage{(nevent/(double)nevents_TChain)*100}; //Percentage of events analyzed in the Event List.
		
		// Calculating q and W^2
		calcq(Pbeam,epx,epy,epz,ep); // Calculates the q vector of the scattering reaction.
		double W2{calcW2(Ptarg)};    // Calculates the invariant mass squared (W^2) of the virtual photon - nucleon system.
		h1_W2->Fill(W2);          // Fill the 1D histogram of W^2 of all the events analyzed.

		h1_bb_trackvertz->Fill(vz[0]);
		
		h1_bb_ps_e->Fill(bbps_e);

		double ps_sh_e{bbsh_e+bbps_e};
		h1_bb_sh_plus_ps_e->Fill(ps_sh_e);

		h2_bb_sh_ps_energycorrelation->Fill(bbsh_e,bbps_e);

		h1_eprime_EoverP->Fill(ps_sh_e/ep[0]);
		
		h1_hcal_e->Fill(hcal_e);
		
 		// Calculate BBCal and HCal coincidence times
		double bbcal_time{0.};
		double hcal_time{0.};
		double coin_time{0.};
		double rf_time{0.};
		get_TDC_times(ndata_tdc,tdc_elemID,tdc_time,bbcal_time,hcal_time,coin_time,rf_time);
		double bbcalHcal_time_diff{hcal_time-bbcal_time};
		h1_bbcal_hcal_tdiff->Fill(bbcalHcal_time_diff);

		h1_hcal_clusblk_ADCtime->Fill(hcal_clusblk_ADC_time[0]);
						
		print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
	}

	std::cout << "\n--- Making histograms for the cut definitions is complete ---\n";

	TCanvas* c[n_histos]; //Define canvasses to print histos.
	for(int i=0; i<n_histos; i++) c[i] = new TCanvas(); 

	gStyle->SetOptFit(1);
 
 	w2_cut(h1_W2,c[0]);
 	TString w2cut_pngName = outputfilename+"_W2cut.png";
 	c[0]->SaveAs(w2cut_pngName);

 	vertexZ_cut(h1_bb_trackvertz,c[1]);
 	TString vzcut_pngName = outputfilename+"_VZcut.png";
 	c[1]->SaveAs(vzcut_pngName);

	preshower_cut(h1_bb_ps_e,c[2]);
	TString pscut_pngName = outputfilename+"_PScut.png";
 	c[2]->SaveAs(pscut_pngName);

 	c[3]->cd();
 	h2_bb_sh_ps_energycorrelation->Draw("COLZ");
 	TString pssh_corr_pngName = outputfilename+"_SH_PS_correlation.png";
 	c[3]->SaveAs(pssh_corr_pngName);

	e_over_p_cut(h1_eprime_EoverP,c[4]);
	TString EoverPcut_pngName = outputfilename+"_EoverPcut.png";
 	c[4]->SaveAs(EoverPcut_pngName);

	calc_bbcal_hcal_thresholds(kine_num,h1_hcal_e,h1_bb_sh_plus_ps_e,c[5],c[6]);
	TString hcalEcut_pngName = outputfilename+"_hcalEcut.png";
 	c[5]->SaveAs(hcalEcut_pngName);
 	TString bbcalEcut_pngName = outputfilename+"_bbcalEcut.png";
 	c[6]->SaveAs(bbcalEcut_pngName);

	//bbCal_HCal_coincut(h1_bbcal_hcal_tdiff,c[7]);
	TString coincut_pngName = outputfilename+"_coincut.png";
 	c[7]->SaveAs(coincut_pngName);

 	hcal_clusblk_ADCtime_cut(h1_hcal_clusblk_ADCtime, c[8]);
 	TString hcal_clusblkADCtimecut_pngName = outputfilename+"_hcal_clusblk_ADCtime.png";
 	c[8]->SaveAs(hcal_clusblkADCtimecut_pngName);
 	
	/*for(int i=0; i<n_histos; i++)
	{
		//c[i]->Modified();
		c[i]->Update();
	}*/

	TString outputpdffilename = outputfilename+".pdf";
	make_pdf(c,outputpdffilename); // Dec 19 - TLines will not show up in the saved PDF files. The reason is unkown and still to be resolved. Saving all the canvases in png works as expected.
	std::cout << "\n--- PDF file " <<  outputpdffilename << " was created with all the histograms ---" <<'\n';
	std::cout << "\n--- PNG files were also created with the corresponding histograms ---" <<'\n';
	
	TString outputconfigfilename = outputfilename+"_foranaconfigfile.cfg";
	make_configfile(outputconfigfilename);
	std::cout << "\n--- Cut threshold inputs for analysis were put in " << outputconfigfilename << ". Copy and paste the contents at the end of analysis script's configuation file. \n";

	//elist->Delete();
	fout->Write();
	std::cout << "\n--- Output histograms are written in into the  " << outputrootfilename << '\n';
	
}

//// End Main Function ////

