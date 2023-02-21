// Written by Anuruddha Rathanayke
// Dec 20, 2022

// 1) Reads in confiugaration file and get the root files to be analyzed and the cuts to be used.
// 2) Run the analysis with all the cuts applied.
// 3) Outputs a root file with all the histogram results.

#include <iostream>
#include <string>
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
#include "readconfigfile.h"
#include "calc_q_W2.h"
#include "calc_HCalintersect.h"


const double target_mass{0.5*(0.938272+0.939565)}; //Average of neutron and proton mass.

//Function to apply global cuts and print out some iformation into the screen.
void apply_globalcuts(TChain* C, TEventList* elist, long& nevents_TChain, long& nevents_eventlist)
{
	nevents_TChain = C->GetEntries();
	std::cout<<"\nNumber of events in the TChain = " << nevents_TChain <<'\n';
	std::cout<<"--- Begin applying global cuts = " << globalcut << " ---" <<'\n';
	C->Draw(">>elist",globalcut);
	std::cout<<"--- Done applying global cuts ---" <<'\n';
	nevents_eventlist = elist->GetN();	
	std::cout<<"Number of events in the eventlist = " << nevents_eventlist <<'\n';
	std::cout<<"Percentage of events accepted for physics analysis = " << (nevents_eventlist/(double)nevents_TChain)*100 <<"%\n";
}

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

void print_analysis_percentage(double cuurentevent_ana_percentage, int& previousevent_ana_percentage_int)
{
	int cuurentevent_ana_percentage_int{(int)cuurentevent_ana_percentage};
	if (cuurentevent_ana_percentage_int%10==0 && cuurentevent_ana_percentage_int>previousevent_ana_percentage_int)
	{
		std::cout << cuurentevent_ana_percentage_int <<"%\n";
	}

	previousevent_ana_percentage_int = cuurentevent_ana_percentage_int;
}

//If using multiple root files with multiple run numbers under the same settings, just input the run number of one of the root files. 
//runnum variable is only used to access the constant parameters for the given kinematic setting.

void elastix_4(const int runnum, const char* configfilename = "setup_analysis.cfg", const char* outputrootfilename = "elastix_4_output.root")
{
	
	readin_beamvariables(runnum); // Reads in kinematic variables that are a constant to the given configuration, from the "beam_variables.h" file.

	TChain* C = new TChain("T"); //Initialize the TChain to chain the root files for analysis.
	
	readin_anaconfigfile(configfilename,C); // Reads in things such as root files to be analyzed, the global cuts to be used, and etc. *Also initializes the TChain (pointer C)*.

	//Defining event list.
	TEventList* elist = new TEventList("elist","");
	long nevents_TChain{0};
	long nevents_eventlist{0};
	apply_globalcuts(C,elist,nevents_TChain,nevents_eventlist); // Make the event list passing the global cuts out of the TChain.

	std::cout<<"\n--- Beginning Physics Analysis ---" <<'\n';

	const int MAXNTRACKS{10};
	const int MAXNTDC{1000};

	//variables needed are BigBite track px,py,pz and sbs hcal x,y,e
 	double ntrack{0.};
 	double vz[MAXNTRACKS];
	double epx[MAXNTRACKS];
	double epy[MAXNTRACKS];
	double epz[MAXNTRACKS];
	double ep[MAXNTRACKS];
	double xhcal {0.};
	double yhcal {0.};
	double hcal_e {0.};
	double bbsh_e {0.};
	double bbps_e {0.};
	double tdc_time[MAXNTDC];
	double tdc_elemID[MAXNTDC];
	int ndata_tdc{0};

	C->SetBranchStatus("*",0);
	C->SetBranchStatus("bb.tr.n",1);
	C->SetBranchStatus("bb.tr.vz",1);
	C->SetBranchStatus("bb.tr.px",1);
	C->SetBranchStatus("bb.tr.py",1);
	C->SetBranchStatus("bb.tr.pz",1);
	C->SetBranchStatus("bb.tr.p",1);
	C->SetBranchStatus("bb.sh.e",1);
	C->SetBranchStatus("bb.ps.e",1);
	C->SetBranchStatus("sbs.hcal.x",1);
	C->SetBranchStatus("sbs.hcal.y",1);
	C->SetBranchStatus("sbs.hcal.e",1);
	C->SetBranchStatus("bb.tdctrig.tdc",1);
	C->SetBranchStatus("bb.tdctrig.tdcelemID",1);
	C->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);

	C->SetBranchAddress("bb.tr.n",&ntrack);
	C->SetBranchAddress("bb.tr.vz",vz);
	C->SetBranchAddress("bb.tr.px",epx);
	C->SetBranchAddress("bb.tr.py",epy);
	C->SetBranchAddress("bb.tr.pz",epz);
	C->SetBranchAddress("bb.tr.p",ep);
	C->SetBranchAddress("bb.sh.e",&bbsh_e);
	C->SetBranchAddress("bb.ps.e",&bbps_e);
	C->SetBranchAddress("sbs.hcal.x",&xhcal);
	C->SetBranchAddress("sbs.hcal.y",&yhcal);
	C->SetBranchAddress("sbs.hcal.e",&hcal_e);
	C->SetBranchAddress("bb.tdctrig.tdc",tdc_time);
	C->SetBranchAddress("bb.tdctrig.tdcelemID",tdc_elemID);
	C->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID",&ndata_tdc);

	TFile* fout = new TFile(outputrootfilename,"RECREATE");

	TH1D* h1_W2 = new TH1D("h1_W2","W^2 of all the events passing the global cut; W^2 (GeV^2/c^4)",250,-1,4);
	TH1D* h1_bbtrackvertz = new TH1D("h1_bbtrackvertz","BB Track vertex Z position; vertex Z (m)",600,-0.15,0.15);
	TH1D* h1_bbshower_e = new TH1D("h1_bbshower_e","Energy Deopsited in BB Shower; Energy (GeV)",500,0,5);
	TH1D* h1_bbpreshower_e = new TH1D("h1_bbpreshower_e","Energy Deopsited in BB Pre Shower; Energy (GeV)",500,0,5);
	TH1D* h1_sh_ps_e = new TH1D("h1_sh_ps_e","Total SH and PS cluster energy sum; SH+PS Energy (GeV)",500,0,5);
	TH1D* h1_eprime_EoverP = new TH1D("h1_eprime_EoverP","Scattered electron E/P; E/P",1000,-1,3);
	TH1D* h1_hcal_e = new TH1D("h1_hcal_e","HCal Energy Deopsited; HCal Energy (GeV)",250,0,0.5);
	TH1D* h1_bbcal_hcal_tdiff = new TH1D("h1_bbcal_hcal_tdiff","HCal time - BBCal time; HCal_{time}-BBCal_{time} (ns)",300,400,700);
	TH2D* h2_dxdy_all = new TH2D("h2_dxdy_all","HCal delta plots - All events;Y_{HCal}-Y_{expected} (m);X_{HCal}-X_{expected} (m)",125,-2,2,125,-4,6);
	TH2D* h2_dxdy_withcuts = new TH2D("h2_dxdy_withcuts","HCal delta plots - good events with cuts;Y_{HCal}-Y_{expected} (m);X_{HCal}-X_{expected} (m)",125,-2,2,125,-4,6);
	TH1D* h1_dx_all = new TH1D("h1_dx_all","dx - All events;X_{HCal}-X_{expected} (m)",125,-4,6);
	TH1D* h1_dx_withcuts = new TH1D("h1_dx_withcuts","dx - good events;X_{HCal}-X_{expected} (m)",125,-4,6);
	TH1D* h1_dy_all = new TH1D("h1_dy_all","dy - All events;Y_{HCal}-Y_{expected} (m)",125,-2,2);
	TH1D* h1_dy_withcuts = new TH1D("h1_dy_withcuts","dy - good events;Y_{HCal}-Y_{expected} (m)",125,-2,2);
	

	//Defining the constant 4-vecors of the incoming beam electrons and the target 
	TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
	TLorentzVector Ptarg(0,0,0,target_mass);

	make_HCal_vectors(hcaldist,sbstheta); //Creates constant vectors that definces the postion of the HCal w.r.t Hall coordinate system.

	//double W2_min{W2_mean-W2_sigma};
	//double W2_max{W2_mean+W2_sigma};

	//Variables to keep track of printout to the terminal to keep track of analysis progress.
	int previousevent_ana_percentage_int{0};

	//const double tdiffmax = 20; // Max deviation from coin via tdctrig cut

	long nevent {0};
	long elastic_yield{0};

	while (C->GetEntry(elist->GetEntry(nevent++)))
	{
		
		double ana_percentage{(nevent/(double)nevents_eventlist)*100}; //Percentage of events analyzed in the Event List.

		////**** Fill certain histograms before applying any tight cuts ****////
		// Calculating q and W^2
		calcq(Pbeam,epx,epy,epz,ep); // Calculates the q vector of the scattering reaction.
		double W2{calcW2(Ptarg)};    // Calculates the invariant mass squared (W^2) of the virtual photon - nucleon system.
		h1_W2->Fill(W2);          // Fill the 1D histogram of W^2 of all the events analyzed.
		h1_bbshower_e->Fill(bbsh_e);
		h1_bbpreshower_e->Fill(bbps_e);
		double ps_sh_e{bbsh_e+bbps_e};
		h1_sh_ps_e->Fill(ps_sh_e);
		double e_over_p{ps_sh_e/ep[0]};
		h1_eprime_EoverP->Fill(e_over_p);
		h1_hcal_e->Fill(hcal_e);
		h1_bbtrackvertz->Fill(vz[0]);
		// Calculate BBCal and HCal coincidence times
		double bbcal_time{0.};
		double hcal_time{0.};
		double coin_time{0.};
		double rf_time{0.};
		get_TDC_times(ndata_tdc,tdc_elemID,tdc_time,bbcal_time,hcal_time,coin_time,rf_time);
		double bbcalHcal_time_diff{hcal_time-bbcal_time};
		h1_bbcal_hcal_tdiff->Fill(bbcalHcal_time_diff);
 		
 		// Calculating the dx and dy for hits on HCal
		double xexpected_hcal{0.};
		double yexpected_hcal{0.};
		calc_expected_xyonHCal(vz,xexpected_hcal,yexpected_hcal); // Calculating the intersection point on HCal.
		double delta_x{xhcal-xexpected_hcal};
		double delta_y{yhcal-yexpected_hcal};
		// Fill dx dy histos with all the events.
		h2_dxdy_all->Fill(delta_y,delta_x);
		h1_dx_all->Fill(delta_x);
		h1_dy_all->Fill(delta_y);
		

		
		////****Elastic and "good event" cuts****////

		// The W^2 cut. The primary cut(?) to select elastic scattering D(e,e'n/p) events.
		if (W2<W2_min || W2>W2_max)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue;			
		}

		// Tighten the aleady implemented vertex cut (in elist) to strictly select events originating from the target
		// IMPORTANT! - vzcut is defined manulally inside "definecuts.C" to match the GMn/nTPE target lenght of 15 cm. It does not get automatically calculated like the other cuts in "definecuts()".
		if (abs(vz[0])>vzcut)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue; //The target is 15 cm long very accurately. Here the cut accepts events from a 16 cm length "target".
		}
		 
		// BB track E/P cut
		if (e_over_p<eoverp_leftcut || e_over_p>eoverp_rightcut)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue;
		}

		// Preshower energy cut. This is for pion rejection. This is already loosely implemented in the eventlis. We tighten the cut further more in here.
		if (bbps_e<pse_min)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue;
		}

		// Total preshower and shower cluster energy cut.
		if (ps_sh_e<bbcal_energycut)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue;
		}

		// HCal energy cut
		if (hcal_e<hcal_energycut)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue;
		}

		// BB and SBS trigger coincidence cut
		/*if (abs(bbcalHcal_time_diff-510)>tdiffmax)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue;
		}*/
		if (bbcalHcal_time_diff<trigdiff_leftcut || bbcalHcal_time_diff>trigdiff_rightcut)
		{
			print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
			continue;
		}		

		////****End cuts****////


		//Filling histograms with the events passing the cuts.
		h2_dxdy_withcuts->Fill(delta_y,delta_x);
		h1_dx_withcuts->Fill(delta_x);
		h1_dy_withcuts->Fill(delta_y);

		print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);
		elastic_yield++;
	}
	
	elist->Delete();
	fout->Write();
	std::cout << "\n--- Analysis completed and the ouput histograms are written to the root file: " << outputrootfilename << " ---" <<'\n';
	std::cout << "\nElastic yield from the analysis = " << elastic_yield <<'\n';
	std::cout << "Percentage elastic yield of the run = " << (elastic_yield/(double)nevents_TChain)*100 <<"%\n";

}