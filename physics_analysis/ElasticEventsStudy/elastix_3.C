//Plots W2, HCal dx dy plots with and without W2 cuts. 

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
void apply_globalcuts(TEventList* elist, long& nevents_TChain, long& nevents_eventlist)
{
	nevents_TChain = C->GetEntries();
	std::cout<<"\nNumber of events in the TChain = " << nevents_TChain <<'\n';
	std::cout<<"--- Begin applying global cuts = "<< globalcut << " ---" <<'\n';
	C->Draw(">>elist",globalcut);
	std::cout<<"--- Done applying global cuts ---"<<'\n';
	nevents_eventlist = elist->GetN();	
	std::cout<<"Number of events in the eventlist = "<< nevents_eventlist <<'\n';
	std::cout<<"Percentage of events accepted for physics analysis = "<< (nevents_eventlist/(double)nevents_TChain)*100 <<"%\n";
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

void elastix_3(const int runnum, const char* configfilename = "setup_analysis.cfg", const char* outputrootfilename = "elastix_3_output.root")
{
	
	readin_beamvariables(runnum); // Reads in kinematic variables that are a constant to the given configuration, from the "beam_variables.h" file.

	readin_configfile(runnum,configfilename); // Reads in things such as root files to be analyzed, the global cuts to be used, and etc. *Also initializes the TChain (pointer C)*.

	//Defining event list.
	TEventList* elist = new TEventList("elist","");
	long nevents_TChain{0};
	long nevents_eventlist{0};
	apply_globalcuts(elist,nevents_TChain,nevents_eventlist); // Make the event list passing the global cuts out of the TChain.

	std::cout<<"\n--- Beginning Physics Analysis ---" <<'\n';

	int MAXNTRACKS=10;

	//variables needed are BigBite track px,py,pz and sbs hcal x,y,e
 	double ntrack{0.};

	double vz[MAXNTRACKS];

	double epx[MAXNTRACKS];
	double epy[MAXNTRACKS];
	double epz[MAXNTRACKS];
	double ep[MAXNTRACKS];
	double sh_e[MAXNTRACKS];
	double ps_e[MAXNTRACKS];

	double xhcal {0.};
	double yhcal {0.};
	double hcal_e {0.};
	double bbsh_e {0.};
	double bbps_e {0.};

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

	TFile* fout = new TFile(outputrootfilename,"RECREATE");

	TH1D* h1W2_all = new TH1D("h1W2_all","W^2 of all the events passing the global cut; W^2 (GeV^2/c^4)",250,-1,4);
	TH2D* h2dxdy_all = new TH2D("h2dxdy_all","All events;#Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
	TH2D* h2dxdy_Wcut = new TH2D("h2dxdy_Wcut","|W^{2}-0.88|<0.4;#Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
	TH1D* h1dx_all = new TH1D("h1dx_all","All events; #Deltax (m)",125,-4,6);
	TH1D* h1dx_Wcut = new TH1D("h1dx_Wcut","|W^{2}-0.88|<0.4;#Deltax (m)",125,-4,6);
	TH1D* h1dy_all = new TH1D("h1dy_all","All events; #Deltay (m)",125,-2,2);
	TH1D* h1dy_Wcut = new TH1D("h1dy_Wcut","|W^2-0.88|<0.4;#Deltay (m)",125,-2,2);
	TH1D* h1bbshower_e = new TH1D("h1bbshower_e","Energy Deopsited in BB Shower; Energy (GeV)",500,0,5);
	TH1D* h1bbpreshower_e = new TH1D("h1bbpreshower_e","Energy Deopsited in BB Pre Shower; Energy (GeV)",500,0,5);
	TH1D* h1eprime_EoverP = new TH1D("h1eprime_EoverP","Scattered electron E/P",1000,-1,3);
	TH1D* h1Hcal_e = new TH1D("h1HCal_e","HCal Energy Deopsited",250,0,0.5);

	//Defining the constant 4-vecors of the incoming beam electrons and the target 
	TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
	TLorentzVector Ptarg(0,0,0,target_mass);

	make_HCal_vectors(hcaldist,sbstheta); //Creates constant vectors that definces the postion of the HCal w.r.t Hall coordinate system.

	double W2_min{W2_mean-W2_sigma};
	double W2_max{W2_mean+W2_sigma};

	//Variables to keep track of printout to the terminal to keep track of analysis progress.
	//double ana_percentage{0.};
	int previousevent_ana_percentage_int{0};

	long nevent {0};
	long elastic_yield{0};

	while (C->GetEntry(elist->GetEntry(nevent++)))
	{
		// Calculating q and W^2
		calcq(Pbeam,epx,epy,epz,ep); // Calculates the q vector of the scattering reaction.
		double W2{calcW2(Ptarg)};    // Calculates the invariant mass squared (W^2) of the virtual photon - nucleon system.
		h1W2_all->Fill(W2);          //Fill the 1D histogram of W^2 of all the events analyzed.
		h1bbshower_e->Fill(bbsh_e);
		h1bbpreshower_e->Fill(bbps_e);
		h1Hcal_e->Fill(hcal_e);
		h1eprime_EoverP->Fill((bbsh_e+bbps_e)/ep[0]);
 		
 		// Calculating the dx and dy for hits on HCal
		double xexpected_hcal{0.};
		double yexpected_hcal{0.};
		calc_expected_xyonHCal(vz,xexpected_hcal,yexpected_hcal); // Calculating the intersection point on HCal.
		double delta_x{xhcal-xexpected_hcal};
		double delta_y{yhcal-yexpected_hcal};
		h2dxdy_all->Fill(delta_y,delta_x);
		h1dx_all->Fill(delta_x);
		h1dy_all->Fill(delta_y);
		

		if (W2>W2_min && W2<W2_max)
		{
			h2dxdy_Wcut->Fill(delta_y,delta_x);
			h1dx_Wcut->Fill(delta_x);
			h1dy_Wcut->Fill(delta_y);
		}

		double ana_percentage{(nevent/(double)nevents_eventlist)*100}; //Percentage of events analyzed in the Event List.
		print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);

	}
	
	elist->Delete();
	fout->Write();
	std::cout << "\n--- Analysis completed and the ouput histograms are written to the root file: " << outputrootfilename << " ---" <<'\n';

}