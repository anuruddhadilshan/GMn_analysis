///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Written by Anuruddha Rathanayke																	 						     //
// Feb 8 , 2023																												     //
// This is an extention of elastix_5.C script. With an addition of the function of outputting a TTree with the analysis results. // 
// 																															     //
// 1) Reads in confiugaration file and get the *parsed* root files to be analyzed and the cuts to be used.						 //
// 2) Run the analysis with all the cuts applied and makes the "dx", "dy" and "dx vs dy" plots. 								 //
// 3) Outputs a root file with analysis results histograms and a TTree.														     //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <chrono>
#include "TChain.h"
#include "TFile.h"
#include "TCut.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include <Math/Vector4D.h>
#include "TVector3.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include "includes/beam_variables.h"
#include "includes/readconfigfile.h"
#include "includes/HCalConstants.h"
#include "includes/calc_HCalintersect.h"
#include "includes/fiducialcut.h"
#include "includes/constants.h"
#include "includes/exprconstants.h"
#include "includes/calc_thetapq.h"


const double target_mass{0.5*(Constants::n_mass + Constants::p_mass)}; //Average of neutron and proton mass.

void apply_globalcuts(TChain*, TEventList*, long&, long&);
void get_TDC_times(int, double [1000], double [1000], double&, double&, double&, double&);
void print_analysis_percentage(double, int&);


void gmn_ana_dxdy(const int kine_num, const double sbsfieldscale, const char* configfilename = "setup_dxdy_analysis.cfg", const char* outputrootfilename = "gmn_ana_output")
{
	auto total_time_start = std::chrono::high_resolution_clock::now();

	lookup_kinematic_info(kine_num); // Reads in kinematic variables that are a constant to the given configuration, from the "beam_variables.h" file.

	TChain* C = new TChain("T"); //Initialize the TChain to chain the root files for analysis.
	
	readin_anaconfigfile(configfilename,C); // Reads in things such as root files to be analyzed, the global cuts to be used, and etc. *Also initializes the TChain (pointer C)*.

	//Defining event list.
	TEventList* elist = new TEventList("elist","");
	long nevents_TChain{0};
	long nevents_eventlist{0};
	//apply_globalcuts(C,elist,nevents_TChain,nevents_eventlist); // Make the event list passing the global cuts out of the TChain.
	nevents_TChain = C->GetEntries();

	double average_proton_deflection{getavg_proton_deflection(kine_num,sbsfieldscale)}; // This value is needed for the fiducial cut. NEEDS to be RE EVALUATED.

	std::cout<<"\n--- Beginning Physics Analysis ---" <<'\n';

	const int MAXNTRACKS{10};
	const int MAXNTDC{1000};

	// Define the variables needed to save information from ROOT tree for each event.
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
	double hcal_clusblk_ADC_time[15]; // Maximum number of blocks in a cluster is 15 as per S.Seeds.

	// Initialize the TTree variables.
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
	C->SetBranchStatus("sbs.hcal.clus_blk.atime",1);

	// Associate the TTree variables with the local varibales defined above.
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
	C->SetBranchAddress("sbs.hcal.clus_blk.atime",hcal_clusblk_ADC_time);


	TFile* fout = new TFile(Form("%s.root",outputrootfilename),"RECREATE");

	TTree* tree = new TTree("T","Analysis outputs from the gmn_ana_dxdy.C script");

	// Define variables that are not already define that needs to be written to the analys output tree,
	double Q2{0.};
	double W2{0.};
	double theta_pq{0.};
	double xexpected_hcal{0.};
	double yexpected_hcal{0.};
	double delta_x{0.};
	double delta_y{0.};
	double bbcalHcal_time_diff{0};
	double ps_sh_e{0.};
	double e_over_p{0.};

	// Output tree branch definitions.
	tree->Branch("ekine.Q2", &Q2);
	tree->Branch("ekine.W2", &W2);
	tree->Branch("theta_pq", &theta_pq);
	tree->Branch("bb.tr.vz0", &vz[0]);
	tree->Branch("bb.sh.e", &bbsh_e);
	tree->Branch("bb.ps.e", &bbps_e);
	tree->Branch("ekine_E_over_P", &e_over_p);
	tree->Branch("sbs.hcal.e", &hcal_e);
	tree->Branch("hcal_bbcal_TDC_trigtime_diff", &bbcalHcal_time_diff);
	tree->Branch("sbs.hcal.clus_blk.atime", &hcal_clusblk_ADC_time[0]);
	tree->Branch("sbs.hcal.xdet", &xhcal);
	tree->Branch("sbs.hcal.ydet", &yhcal);
	tree->Branch("sbs.hcal.xexpect", &xexpected_hcal);
	tree->Branch("sbs.hcal.yexpect", &yexpected_hcal);
	tree->Branch("sbs.hcal.dx", &delta_x);
	tree->Branch("sbs.hcal.dy", &delta_y);

	// Histogram definitions.
	TH1D* h1_W2_all = new TH1D("h1_W2_all","W^2 of all the events passing the global cut; W^2 (GeV^2/c^4)",250,-1,4);
	TH1D* h1_Q2_all = new TH1D("h1_Q2_all", "Q^2 of all the events passing the global cut; Q^2 (=-q^2) (GeV^2/c^2)",200,-2,18);
	TH1D* h1_W2_withcuts = new TH1D("h1_W2_withcuts","W^2 of events passing all the cuts; W^2 (GeV^2/c^4)",250,-1,4);
	TH1D* h1_Q2_withcuts = new TH1D("h1_Q2_withcuts", "Q^2 of events passing all the cuts; Q^2 (=-q^2) (GeV^2/c^2)",200,-2,18);
	TH1D* h1_thetapq_all = new TH1D("h1_thetapq_all", "#theta_{pq} for all the events without eleastic cuts; #theta_{pq}", 1000, 0, 10);
	TH2D* h2_thetapq_vs_W2_all = new TH2D("h2_thetapq_vs_W2_all", "#theta_{pq} vs W^{2} for all the events without elastic cuts; W^{2}; #theta_{pq} degrees", 100, 0.4, 1.4, 1000, 0, 10);
	TH1D* h1_bbtrackvertz_withcuts = new TH1D("h1_bbtrackvertz_withcuts","BB Track vertex Z position; vertex Z (m)",600,-0.15,0.15);
	TH1D* h1_bbshower_e_withcuts = new TH1D("h1_bbshower_e_withcuts","Energy Deopsited in BB Shower; Energy (GeV)",500,0,5);
	TH1D* h1_bbpreshower_e_withcuts = new TH1D("h1_bbpreshower_e_withcuts","Energy Deopsited in BB Pre Shower; Energy (GeV)",500,0,5);
	TH1D* h1_sh_ps_e_withcuts = new TH1D("h1_sh_ps_e_withcuts","Total SH and PS cluster energy sum; SH+PS Energy (GeV)",500,0,5);
	TH1D* h1_eprime_EoverP_withcuts = new TH1D("h1_eprime_EoverP_withcuts","Scattered electron E/P; E/P",1000,-1,3);
	TH1D* h1_hcal_e_withcuts = new TH1D("h1_hcal_e_withcuts","HCal Energy Deopsited; HCal Energy (GeV)",750,0,1.5);
	//TH1D* h1_bbcal_hcal_tdiff_withcuts = new TH1D("h1_bbcal_hcal_tdiff_withcuts","HCal time - BBCal time; HCal_{time}-BBCal_{time} (ns)",300,400,700);
	TH1D* h1_hcal_clusblk_ADCtime_withcuts = new TH1D("h1_hcal_clusblk_ADCtime_withcuts","ADC time of the highest energy block in the largest cluster; ADC Time (ns)",300,-100,200);
	TH2D* h2_dxdy_all = new TH2D("h2_dxdy_all","HCal delta plot - All events;Y_{HCal}-Y_{expected} (m);X_{HCal}-X_{expected} (m)",125,-2,2,125,-4,6);
	TH2D* h2_dxdy_withcuts_beforefiducial = new TH2D("h2_dxdy_withcuts_beforefiducial","HCal delta plots - good events with all cuts except fiducial cut;Y_{HCal}-Y_{expected} (m);X_{HCal}-X_{expected} (m)",125,-2,2,125,-4,6);
	TH2D* h2_dxdy_withcuts = new TH2D("h2_dxdy_withcuts","HCal delta plot - good events with all cuts;Y_{HCal}-Y_{expected} (m);X_{HCal}-X_{expected} (m)",125,-2,2,125,-4,6);
	TH2D* h2_dxdy_failfiducial = new TH2D("h2_dxdy_failfiducial","HCal delta plot - events that fail the fiducial cut;Y_{HCal}-Y_{expected} (m);X_{HCal}-X_{expected} (m)",125,-2,2,125,-4,6);
	//TH2D* h2_dxdy_ = new TH2D("h2_dxdy_withcuts_afterfiducial","HCal delta plots - good events with all cuts;Y_{HCal}-Y_{expected} (m);X_{HCal}-X_{expected} (m)",125,-2,2,125,-4,6);
	TH1D* h1_dx_all = new TH1D("h1_dx_all","dx - All events;X_{HCal}-X_{expected} (m)",1000,-4,6);
	TH1D* h1_dx_withcuts = new TH1D("h1_dx_withcuts","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	TH1D* h1_dy_all = new TH1D("h1_dy_all","dy - All events;Y_{HCal}-Y_{expected} (m)",400,-2,2);
	TH1D* h1_dy_withcuts = new TH1D("h1_dy_withcuts","dy - with cuts;Y_{HCal}-Y_{expected} (m)",400,-2,2);
	//TH1D* h1_hcalx = new TH1D("h1_hcalx","HCal X - with cuts",125,-4,6);
	//TH1D* h1_hcaly = new TH1D("h1_hcaly","HCal Y - with cuts",125,-2,2);
	//TH2D* h2_hcal_xy = new TH2D("h2_hcal_xy","HCal XY plots - good events with cuts;Y_{HCal} (m);X_{HCal} (m)",125,-2,2,125,-4,6);
	TH2D* h2_hcal_xyexpected = new TH2D("h2_hcal_xyexpected","HCal XY_expected Surviving Fiducial Cut;Y_{Expected} (m);X_{Expected} (m)",125,-2,2,125,-4,6);
	TH2D* h2_hcal_xyexpected_failfiducial = new TH2D("h2_hcal_xyexpected_failfiducial","HCal XY_expected Fail Fiducial Cut;Y_{Expected} (m);X_{Expected} (m)",125,-2,2,125,-4,6);
	TH1D* h1_thetapq = new TH1D("h1_thetapq", "#theta_{pq} after elastic cuts; #theta_{pq}", 1000, 0, 10);
	TH2D* h2_thetapq_vs_W2 = new TH2D("h2_thetapq_vs_W2", "#theta_{pq} vs W^{2}; W^{2}; #theta_{pq} degrees", 100, 0.4, 1.4, 1000, 0, 10);

	//Defining the constant 4-vecors of the incoming beam electrons and the target 
	ROOT::Math::PxPyPzEVector Pbeam(0,0,Ebeam,Ebeam);
	ROOT::Math::PxPyPzEVector Ptarg(0,0,0,target_mass);

	HCalVectors hcal; // Defined an object named "hcal" 
	hcal.make_HCal_vectors(hcaldist, hcaltheta); //Creates constant vectors that definces the postion of the HCal w.r.t Hall coordinate system.
	//make_HCal_vectors(hcaldist, hcaltheta);

	//const double tdiffmax = 20; // Max deviation from coin via tdctrig cut.

	//Variables to keep track of printout to the terminal.
	int previousevent_ana_percentage_int{0};
	
	long nevent {0};
	long elastic_yield{0};

	//variables to keep track the number of events get cut out by each cut.
	long n_w2cut{0};
	long n_vzcut{0};
	long n_eoverpcut{0};
	long n_psecut{0};
	long n_psandshecut{0};
	long n_hcalecut{0};
	long n_HCalandBBCaltimecorrcut{0};
	long n_fiducialcut{0};
	long n_hcal_clusblk_atime_cut{0};

	while (C->GetEntry(nevent++))
	{
		
		double ana_percentage{(nevent/(double)nevents_TChain)*100}; //Percentage of events analyzed in the Event List.
		print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int);

		// Calculating q 
		ROOT::Math::PxPyPzEVector kprime; //Four vector of the scattered electron.
		kprime.SetPxPyPzE(epx[0],epy[0],epz[0],ep[0]); 
		ROOT::Math::PxPyPzEVector q;     //Four momentum trasnferred to the scattered nucleon.
		q = Pbeam - kprime; 
		Q2 = -q.M2();
		W2 = (Ptarg+q).M2();    // Calculates the invariant mass squared (W^2) of the virtual photon - nucleon system.

		h1_W2_all->Fill(W2);
		h1_Q2_all->Fill(Q2);

		// Calculating the dx and dy for hits on HCal
		hcal.calc_expected_xyonHCal(q, vz); // Calculating the intersection point on HCal.
		//calc_expected_xyonHCal(q, vz, xexpected_hcal, yexpected_hcal);
		xexpected_hcal = hcal.return_xexpected();
		yexpected_hcal = hcal.return_yexpected();
		delta_x = xhcal-xexpected_hcal;
		delta_y = yhcal-yexpected_hcal;
		// Fill dx dy histos with all the events.
		h2_dxdy_all->Fill(delta_y,delta_x);
		h1_dx_all->Fill(delta_x);
		h1_dy_all->Fill(delta_y);

		// Calculate theta_pq
		Calc_thetapq calc_thetapq;
		theta_pq = calc_thetapq.return_thetaPq(hcal, xhcal, yhcal);

		h1_thetapq_all->Fill(theta_pq);
		h2_thetapq_vs_W2_all->Fill(W2,theta_pq);
						
		// Calculate BBCal and HCal coincidence times
		double bbcal_time{0.};
		double hcal_time{0.};
		double coin_time{0.};
		double rf_time{0.};
		get_TDC_times(ndata_tdc,tdc_elemID,tdc_time,bbcal_time,hcal_time,coin_time,rf_time);
		bbcalHcal_time_diff = hcal_time - bbcal_time;

		
		////****Elastic and "good event" cuts****////

		// Count the number of events that will get removed from each cut //
		bool b_w2cut{false};
		if (W2<W2_min || W2>W2_max) 
		{ 	
			b_w2cut = true;
			n_w2cut++;
		}

		bool b_vzcut{false};
		if (abs(vz[0])>vzcut)
		{
			b_vzcut = true;
			n_vzcut++;
		}

		ps_sh_e = bbsh_e+bbps_e;
		e_over_p = ps_sh_e/ep[0];
		bool b_eoverpcut{false};
		if (e_over_p < eoverp_leftcut || e_over_p > eoverp_rightcut)
		{
			b_eoverpcut = true;
			n_eoverpcut++;	
		} 

		bool b_psecut{false};
		if (bbps_e<pse_min)
		{
			b_psecut = true;
			n_psecut++;
		}

		bool b_psshecut{false};
		if (ps_sh_e<bbcal_energycut)
		{
			b_psshecut = true;
			n_psandshecut++;
		}

		bool b_hcalectu{false};
		if (hcal_e<hcal_energycut)
		{
			b_hcalectu = true;
			n_hcalecut++;
		}

		bool b_hcalbbcaltimecorrcut{false};
		if (bbcalHcal_time_diff < trigdiff_leftcut || bbcalHcal_time_diff > trigdiff_rightcut)
		{
			b_hcalbbcaltimecorrcut = true;
			n_HCalandBBCaltimecorrcut++;
		}

		bool b_hcal_clusblk_ADC_cut{false};
		if ( hcal_clusblk_ADC_time[0] <  hcal_clusblk_ADCtime_leftcut || hcal_clusblk_ADC_time[0] > hcal_clusblk_ADCtime_rightcut )
		{
			b_hcal_clusblk_ADC_cut = true;
			n_hcal_clusblk_atime_cut++;
		}
		
		// Fiducial cut
		bool is_fiducialcut_pass{fiducial_cut(average_proton_deflection,xexpected_hcal,yexpected_hcal)};
		if (!is_fiducialcut_pass) n_fiducialcut++;
		////

		// Now applying the cuts in the order that I am confident about them. i.e: most confident ones first so that what get cuts out at the beginning are definite bad events and then I can play around with the thresholds o the last ones.
		// The W^2 cut. The primary cut(?) to select elastic scattering D(e,e'n/p) events.
		if (b_w2cut) continue;
		
		// Tighten the aleady implemented vertex cut to strictly select events originating from the target region.
		// IMPORTANT! - vzcut is defined manulally inside "definecuts.C" to match the GMn/nTPE target lenght of 15 cm. It does not get automatically calculated like the other cuts in "definecuts()".
		if (b_vzcut) continue; //The target is 15 cm long very accurately. Here the cut accepts events from a 16 cm length "target".

		// BB track E/P cut
		if (b_eoverpcut) continue;

		// Preshower energy cut. This is for pion rejection. This is already loosely implemented in the eventlis. We tighten the cut further more in here.
		if (b_psecut) continue;

		// BB and SBS trigger coincidence cut
		/*if (abs(bbcalHcal_time_diff-510)>tdiffmax) print_analysis_percentage(ana_percentage,previousevent_ana_percentage_int) continue;
		}*/
		//if (b_hcalbbcaltimecorrcut) continue; // Cut is removed as it was found that HCal f1 TDC timing is not reliable.

		// HCal ADC time cut
		if (b_hcal_clusblk_ADC_cut) continue; // Using HCal ADC time cut instead of BBCal and HCal coincidence cut.

		// Total preshower and shower cluster energy cut.
		if (b_psshecut) continue;

		// HCal energy cut
		if (b_hcalectu) continue;

		h2_dxdy_withcuts_beforefiducial->Fill(delta_y,delta_x);

		if (!is_fiducialcut_pass) 
		{
			h2_dxdy_failfiducial->Fill(delta_y,delta_x);
			h2_hcal_xyexpected_failfiducial->Fill(yexpected_hcal,xexpected_hcal);
			continue;
		}
		
		////****End cuts****////


		//Filling histograms with the events passing the cuts.
		h1_W2_withcuts->Fill(W2);     // Fill the 1D histogram of W^2 of all the events analyzed.
		h1_Q2_withcuts->Fill(Q2);
		h1_bbtrackvertz_withcuts->Fill(vz[0]);
		h1_bbshower_e_withcuts->Fill(bbsh_e);
		h1_bbpreshower_e_withcuts->Fill(bbps_e);
		h1_sh_ps_e_withcuts->Fill(ps_sh_e);
		h1_eprime_EoverP_withcuts->Fill(e_over_p);
		h1_hcal_e_withcuts->Fill(hcal_e);
		//h1_bbcal_hcal_tdiff_withcuts->Fill(bbcalHcal_time_diff);
		h1_hcal_clusblk_ADCtime_withcuts->Fill(hcal_clusblk_ADC_time[0]);
		h2_dxdy_withcuts->Fill(delta_y,delta_x);
		h1_dx_withcuts->Fill(delta_x);
		h1_dy_withcuts->Fill(delta_y);
		//h1_hcalx->Fill(xhcal);
		//h1_hcaly->Fill(yhcal);
		//h2_hcal_xy->Fill(yhcal,xhcal);
		h2_hcal_xyexpected->Fill(yexpected_hcal,xexpected_hcal);
		h1_thetapq->Fill(theta_pq);
		h2_thetapq_vs_W2->Fill(W2,theta_pq);

		// Fill the tree
		tree->Fill();

		elastic_yield++;
	}
	
	//elist->Delete();
	tree->Write();
	fout->Write();

	//draw_fiducialcut(h2_dxdy_withcuts_beforefiducial,Form("%s_h2_dxdy_beforefid.png",outputrootfilename));
	draw_fiducialcut(h2_dxdy_withcuts,Form("%s_h2_dxdy_allcuts.png",outputrootfilename));
	//draw_fiducialcut(h2_hcal_xy,Form("%s_h2_hcal_xy.png",outputrootfilename));
	//draw_fiducialcut(h2_dxdy_failfiducial,Form("%s_h2_dxdy_failfiducial.png",outputrootfilename));
	//draw_fiducialcut(h2_hcal_xyexpected,Form("%s_h2_hcal_xyexpected.png",outputrootfilename));
	//draw_fiducialcut(h2_hcal_xyexpected_failfiducial,Form("%s_h2_hcal_xyexpected_failfiducial.png",outputrootfilename));

	std::cout << "\n--- Percentage of events cut out by different cuts ---\n";
	std::cout << "W^{2} cut = " << (n_w2cut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "Vertex Z cut = " << (n_vzcut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "E over P cut = " << (n_eoverpcut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "PS energy cut = " << (n_psecut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "PS and SH cluster energy sum cut = " << (n_psandshecut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "HCal cluster energy sum cut = " << (n_hcalecut/(double)nevents_TChain)*100 <<"%\n";
	//std::cout << "HCal and BBCal time correlation cut = " << (n_HCalandBBCaltimecorrcut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "HCal ADC time cut = " << (n_hcal_clusblk_atime_cut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "Fiducial cut = "        << (n_fiducialcut/(double)nevents_TChain)*100 <<"%\n";
	std::cout << "---                                                ---\n";

	std::cout << "\n   ##################################################  \n";
	std::cout << "\n--- Analysis completed and the ouput histograms are written to the root file: " << outputrootfilename << " ---" <<'\n';
	std::cout << "\nElastic yield from the analysis = " << elastic_yield <<'\n';
	std::cout << "\nPercentage elastic yield of the run = " << (elastic_yield/(double)nevents_TChain)*100 <<"%\n";

	auto total_time_end = std::chrono::high_resolution_clock::now();
	auto total_time_duration = std::chrono::duration_cast<std::chrono::seconds>(total_time_end - total_time_start);
	std::cout << "\nTotal time taken for analysis: " << total_time_duration.count() << " seconds.\n";

}
//// End Main Function ////

//Function to apply global cuts and print out some iformation into the screen.
void apply_globalcuts(TChain* C, TEventList* elist, long& nevents_TChain, long& nevents_eventlist)
{
	nevents_TChain = C->GetEntries();
	std::cout<<"\nNumber of events in the TChain = " << nevents_TChain <<'\n';
	std::cout<<"--- Begin applying global cuts = " << globalcut << " ---" <<'\n';
	C->Draw(">>elist",globalcut); // "globalcut" is read in from the configuration file.
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
