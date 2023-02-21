// Selects events peaking at the 3rd or 4th timesamples
// Applies a "good signal shape cut" to only consider events that has a desired APV signal shape
// Fit those events using the APV function and extract t0, tau, and max ADC values and make 1D histograms. Also make (trigger_time - t0) plots
// Apply Time Deconvolution to the above same events and find at time sample the events are peaking


#include <iostream>
#include <string>

void copyADCsamples(Int_t ientry, const Double_t adcsamples[46080], Double_t raw_6ADCsamples[6])
{
	const Int_t ntimesamples {6};

	for(Int_t itimesamp = 0; itimesamp < ntimesamples; ++itimesamp)
	{
		raw_6ADCsamples[itimesamp] = adcsamples[ientry+itimesamp];
	}
}		

//Function to determine whether the signal we are going to analyse has a "good" shape.
//I.e. TS1 < TS2 < TS3 < TS4 > TS5 > TS6 or TS1 < TS2 < TS3 > TS4 > TS5 > TS6
bool goodSignalShapeCut(const Double_t raw_6ADCsamples[6], Double_t maxTS)
{
	if (maxTS == 2)
	{
		if (raw_6ADCsamples[0]<raw_6ADCsamples[1] && raw_6ADCsamples[1]<raw_6ADCsamples[2] && raw_6ADCsamples[2]>raw_6ADCsamples[3] && raw_6ADCsamples[3]>raw_6ADCsamples[4] && raw_6ADCsamples[4]>raw_6ADCsamples[5]) return true;
		else false;
	}

	if (maxTS == 3)
	{
		if (raw_6ADCsamples[0]<raw_6ADCsamples[1] && raw_6ADCsamples[1]<raw_6ADCsamples[2] && raw_6ADCsamples[2]<raw_6ADCsamples[3] && raw_6ADCsamples[3]>raw_6ADCsamples[4] && raw_6ADCsamples[4]>raw_6ADCsamples[5]) return true;
		else false;
	}

	return false;
}

void deconvolutionADCsamples(const Double_t raw_6ADCsamples[6], Double_t deconv_6ADCsamples[6])
{
	const Double_t delta_t {24.0}; // nanoseconds
	const Double_t tau {56.0};     // nanoseconds

	Double_t x {delta_t/tau};
	Double_t w_1 {exp(x-1)/x};
	Double_t w_2 {-2*exp(-1)/x};
	Double_t w_3 {exp(-x-1)/x};

	const Int_t ntimesamples = 6;

	for (Int_t ntimesample = 2; ntimesample < ntimesamples; ++ntimesample)
	{
		deconv_6ADCsamples[ntimesample] = w_1*raw_6ADCsamples[ntimesample] + w_2*raw_6ADCsamples[ntimesample-1] + w_3*raw_6ADCsamples[ntimesample-2];
	}
}

Double_t apv_function(Double_t* x, Double_t* par)
{	
	Double_t arg=0;

	if (par[2] != 0)
	{	
		arg = (x[0]-par[1])/par[2];
	}	
	Double_t fitval = par[0]*arg*TMath::Exp(-arg+1);
	return fitval;
}		

void apply_fit(TGraph* g, TF1* fitfunction, Double_t maxADCvalue)
{
	fitfunction->SetParameters(maxADCvalue,0.,56.);
	fitfunction->SetParNames("MaxADC","t0","tau");
	
	g->Fit("apv_function","RQ");
}

Double_t get_bbcal_trigtime(Double_t bb_sh_nclus, Int_t Ndata_bb_tdctrig_tdcelemID, Double_t bb_tdctrig_tdcelemID[6], Double_t bb_tdctrig_tdc[6])
{
	Double_t trig_time {0.};

	if (bb_sh_nclus == 0) return trig_time;

	for (Int_t ihit=0; ihit < Ndata_bb_tdctrig_tdcelemID; ++ihit)
	{	
		if(bb_tdctrig_tdcelemID[ihit]==5) trig_time = bb_tdctrig_tdc[ihit];
	}

	return trig_time;
}

//Define as global variables so they can be accesse by other functions.
Int_t max_timesample {0};
Double_t max_ADC {0.};	
void findmaxtimesample_and_ADC(const Double_t sixADCsamples[6])
{
	const Int_t ntimesamples = 6;
	max_timesample = 0;
	max_ADC = sixADCsamples[0];

	for (Int_t ntimesample=0; ntimesample < ntimesamples-1; ++ntimesample)
	{
		if (sixADCsamples[ntimesample] < sixADCsamples[ntimesample+1]) 
		{	
			max_timesample = ntimesample+1;
			max_ADC = sixADCsamples[ntimesample+1];
		}	
	}
	//if (max_ADC > 100) return max_timesample;
	//else return 6;
	//return max_timesample;
}

void fill_histos_noDeconvolution(TF1* fitfunction, Double_t bbcal_trigtime, TH1F* h1_maxADC, TH1F* h1_t0, TH1F* h1_tau, TH1F* h1_trigtime, TH1F* h1_t0_minus_trigtime, TH2F* h2_t0vsTrigTime)
{	
	Double_t maxADC {fitfunction->GetParameter(0)};
	Double_t t0 {fitfunction->GetParameter(1)};
	Double_t tau {fitfunction->GetParameter(2)};

	h1_maxADC->Fill(maxADC);
	h1_t0->Fill(t0);
	h1_tau->Fill(tau);
	h1_trigtime->Fill(bbcal_trigtime);
	h1_t0_minus_trigtime->Fill(bbcal_trigtime - t0);
	h2_t0vsTrigTime->Fill(bbcal_trigtime, t0);
	//h1_deconv_maxTS->Fill(findmaxtimesample(sixADCsamples));
}

void fill_histos_withDeconvolution(TH1D* h1_deconv_maxTS, TH1F* h1_deconv_maxADC, TH1D* h1_deconv_maxTS_ADCbiggerthan50, TH1D* h1_deconv_maxTS_ADCbiggerthan100)
{
	h1_deconv_maxTS->Fill(max_timesample);
	h1_deconv_maxADC->Fill(max_ADC);
	if (max_ADC > 50) h1_deconv_maxTS_ADCbiggerthan50->Fill(max_timesample);
	if (max_ADC > 100) h1_deconv_maxTS_ADCbiggerthan100->Fill(max_timesample);
}

void make_TMultiGraphs(Int_t nevents_display, TGraph* raw_gr[nevents_display], TGraph* deconv_gr[nevents_display], TMultiGraph* mul_gr[nevents_display])
{
	for(Int_t ievent_display = 0; ievent_display < nevents_display; ++ievent_display)
	{	
		raw_gr[ievent_display]->SetTitle("Original ADC signal");
		raw_gr[ievent_display]->SetMarkerStyle(21);
		raw_gr[ievent_display]->SetMarkerColor(1);

		deconv_gr[ievent_display]->SetTitle("Deconvoluted ADC signal");
		deconv_gr[ievent_display]->SetMarkerStyle(22);
		deconv_gr[ievent_display]->SetMarkerColor(2);

		TString name_mul_gr = "mg"+std::to_string(ievent_display);
		TString title_mul_gr = "Event Number: "+std::to_string(ievent_display);
		mul_gr[ievent_display] = new TMultiGraph(name_mul_gr,title_mul_gr);
		mul_gr[ievent_display]->Add(raw_gr[ievent_display]);
		mul_gr[ievent_display]->Add(deconv_gr[ievent_display]);
		mul_gr[ievent_display]->GetXaxis()->SetTitle("Time (ns)");
		mul_gr[ievent_display]->GetYaxis()->SetTitle("ADC");
	}
}

void make_pdfs_for_displayed_events(Int_t nevents_display, TMultiGraph* gr[nevents_display], TString output_file_name)
{
	TCanvas* c[nevents_display];

	TString pdffilename = output_file_name;
	TString openfilename = pdffilename+"(";
	TString closefilename = pdffilename+")";

	Double_t lmargin=0.15;
  	Double_t rmargin=0.15;
    Double_t bmargin=0.15;
    Double_t tmargin=0.09;

	for (Int_t ievent_display = 0; ievent_display < nevents_display; ++ievent_display)
	{
		TString event = "Event Number: "+std::to_string(ievent_display);
		c[ievent_display] = new TCanvas(event, event);
		c[ievent_display]->SetGrid();
		gr[ievent_display]->Draw("AP");
		//c[ievent_display]->BuildLegend();
		
		if(ievent_display == 0) c[ievent_display]->Print(openfilename,"Title:"+event);
		else	if (ievent_display == nevents_display-1) c[ievent_display]->Print(closefilename,"Title:"+event);
		else c[ievent_display]->Print(pdffilename,"Title:"+event);
	}
}

void fit_test3(const Int_t run_number_int = 11590, Int_t nevents_display = 10, const TString output_file_name = "output_with_fits.pdf")
{
	TChain* tchain_T = new TChain("T");
	const TString run_number_string = std::to_string(run_number_int);
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg0_0*.root");
	
	// Disable all the unused branches to save time and enable only the branches used for analysis
	tchain_T->SetBranchStatus("*",false);
	tchain_T->SetBranchStatus("Ndata.bb.gem.m3.strip.ADCsamples",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.ADCsamples",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.nstripsfired",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.ontrack",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.isampmax",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.ADCmax",true);
	tchain_T->SetBranchStatus("bb.tdctrig.tdc",true);
	tchain_T->SetBranchStatus("bb.tdctrig.tdcelemID",true);
	tchain_T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",true);
	tchain_T->SetBranchStatus("bb.sh.nclus",true);

	//Define local the variables to store data from Tree
	Int_t ndata_stripADCsamples{0}; // Variable to hold the number of ADC samples for the event
	Double_t adcsamples[46080]{0.}; // Array to hold ADC smaples with max number of elements needed for a U-V GEM
	Double_t nstripsfired{0.};      // Number of strips fired for the event
	Double_t ontrack[7680]{0};      // The strip is on a track or not. Max number of elements needed for a U-V GEM
	Double_t maxtimesamp[7680] {0}; // Max time sample of the hit.
	Double_t maxADCvalue[7680] {0}; // Max ADC value of the hit.
	Int_t Ndata_bb_tdctrig_tdcelemID{0};
	Double_t bb_tdctrig_tdcelemID[6]{0.};
	Double_t bb_tdctrig_tdc[6]{0.};
	Double_t bb_sh_nclus{0.};

	//Associate local variables with the Tree variables
	tchain_T->SetBranchAddress("Ndata.bb.gem.m3.strip.ADCsamples",&ndata_stripADCsamples);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ADCsamples",&adcsamples);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.nstripsfired",&nstripsfired);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ontrack",&ontrack);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.isampmax",&maxtimesamp);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ADCmax",&maxADCvalue);
	tchain_T->SetBranchAddress("bb.tdctrig.tdc",&bb_tdctrig_tdc);
	tchain_T->SetBranchAddress("bb.tdctrig.tdcelemID",&bb_tdctrig_tdcelemID);
	tchain_T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID",&Ndata_bb_tdctrig_tdcelemID);
	tchain_T->SetBranchAddress("bb.sh.nclus",&bb_sh_nclus);

	//Define histograms to store t0, tau and max_ADC and other parameters that we get out fitting the 6 time sample data 
	TH1F* h1_maxADC = new TH1F("h1_maxADC","Fit Result Distribution for maxADCs; ADC",3500,-500,3000);
	TH1F* h1_t0 = new TH1F("h1_t0","Fit Result Distribution for t0; Time(ns)",400,-100,300);
	TH1F* h1_tau = new TH1F("h1_tau","Fit Result Distribution for APV Time Constant tau; Time(ns) ",500,-200,300); 
	TH1F* h1_trigtime = new TH1F("h1_trigtime", "BBCal Trigger time distribution; Time(ns)", 100, 320 ,420);
	TH1F* h1_t0_minus_trigtime = new TH1F("h1_t0_minus_trigtime","TriggerTime-t0  Distribution; Time(ns)",700,-100,600);
	TH2F* h2_t0vsTrigTime = new TH2F("h2_t0vsTrigTime","t0 and Trigger Time correlation; BBCal Trigger Time (ns); t0 (ns)",100,320,420,100,-50,50);

	TH1D* h1_deconv_maxTS = new TH1D();
	TH1F* h1_deconv_maxADC = new TH1F("h1_deconv_maxADC", "Distribution of maxADC values after Time Deconvolution; ADC",3500,-500,3000);
	TH1D* h1_deconv_maxTS_ADCbiggerthan50 = new TH1D("h1_deconv_maxTS_ADCbiggerthan50","Distribution of Tims Samples with Max ADC > 50 after Time Deconvolution ; Time(ns)",6,-0.5,5.5);
	TH1D* h1_deconv_maxTS_ADCbiggerthan100 = new TH1D("h1_deconv_maxTS_ADCbiggerthan50","Distribution of Tims Samples with Max ADC > 100 after Time Deconvolution ; Time(ns)",6,-0.5,5.5);



	Long64_t nevents = tchain_T->GetEntries();
	constexpr Int_t ntimesamples {6};
	Double_t timesamples[6] {12.0,36.0,60.0,84.0,108.0,132.0};
	Long64_t nevents_analyzed {0}; 
	Long64_t nevents_NOT_analyzed {0};
	TGraph* graph_noDeconv[nevents_display];
	TGraph* graph_Deconv[nevents_display];
	Int_t ievent_display {0};

	for(Long64_t nevent = 0; nevent < nevents; ++nevent)
	{
		tchain_T->GetEntry(nevent);

		Int_t istrip {0};
		Double_t bbcal_trigtime {get_bbcal_trigtime(bb_sh_nclus, Ndata_bb_tdctrig_tdcelemID, bb_tdctrig_tdcelemID, bb_tdctrig_tdc)};

		for(Int_t ientry = 0; ientry < ndata_stripADCsamples; ientry += ntimesamples)
		{
			
			// All the strips that are not on tracks gets regejected from this steps.
			if(ontrack[istrip] == 0.)
			{	
				++istrip;
				continue;
			} 

			bool third_or_fourth_TimeSample_peaking {(maxtimesamp[istrip]==2) || (maxtimesamp[istrip]==3)};	
			if(!third_or_fourth_TimeSample_peaking)
			{
				++nevents_NOT_analyzed; //This is a hit on track but we do not include it.
				++istrip;
				continue;
			}

			//Copy the 6 ADC samples in to a an array of lenght 6.
			Double_t raw_6ADCsamples[ntimesamples] {0.};
			copyADCsamples(ientry, adcsamples, raw_6ADCsamples);

			//std::cout << maxtimesamp[istrip] <<'\n';
			Double_t maxTS = maxtimesamp[istrip];

			if(!goodSignalShapeCut(raw_6ADCsamples, maxTS))
			{
				++nevents_NOT_analyzed; //This is a hit on a track but we do not include it.
				++istrip;
				continue;
			}
			//// The strips left after this step are the ones we are going to use for applying fits. ////

			// Apply the APV function fit to all the selected hits.
			TGraph* istrip_graph = new TGraph(ntimesamples, timesamples, raw_6ADCsamples);
			TF1* fitfunction = new TF1("apv_function", apv_function, 0, 150.0, 3);
			apply_fit(istrip_graph, fitfunction, maxADCvalue[istrip]);
			fill_histos_noDeconvolution(fitfunction, bbcal_trigtime, h1_maxADC, h1_t0, h1_tau, h1_trigtime, h1_t0_minus_trigtime, h2_t0vsTrigTime);

			// Apply time deconvolution to all the selected hits.
			Double_t deconv_6ADCsamples[ntimesamples] {0.};
			deconvolutionADCsamples(raw_6ADCsamples, deconv_6ADCsamples);
			findmaxtimesample_and_ADC(deconv_6ADCsamples);
			fill_histos_withDeconvolution(h1_deconv_maxTS, h1_deconv_maxADC, h1_deconv_maxTS_ADCbiggerthan50, h1_deconv_maxTS_ADCbiggerthan100);
			
			//findmaxtimesample_and_fill_histos_withDeconvolution();
			
			//Copy the TGraphs for the events that are going to be displayed
			if(ievent_display < nevents_display && istrip%100==0)
			{	
				graph_noDeconv[ievent_display] = istrip_graph;
				graph_Deconv[ievent_display] = new TGraph(ntimesamples, timesamples, deconv_6ADCsamples);
				++ievent_display;
			}

			++nevents_analyzed;
			++istrip;
		}

		if(nevent%100 == 0) std::cout <<"Current Event Number = " << nevent <<'\n';
	}

	TMultiGraph* multi_graph[nevents_display];
	make_TMultiGraphs(nevents_display, graph_noDeconv, graph_Deconv, multi_graph);
	make_pdfs_for_displayed_events(nevents_display, multi_graph, output_file_name);

	h1_t0->Fit("gaus","Q");
	h1_tau->Fit("gaus","Q");

	TCanvas* a1 = new TCanvas();
	h1_maxADC->Draw();

	TCanvas* a2 = new TCanvas();
	h1_t0->Draw();

	TCanvas* a3 = new TCanvas();
	h1_tau->Draw();

	TCanvas* a4 = new TCanvas();
	h1_trigtime->Draw();

	TCanvas* a5 = new TCanvas();
	h1_t0_minus_trigtime->Draw();

	TCanvas* a6 = new TCanvas();
	h2_t0vsTrigTime->Draw("colz");

	TCanvas* a7 = new TCanvas();
	h1_deconv_maxTS->Draw();

	TCanvas* a8 = new TCanvas();
	h1_deconv_maxADC->Draw();

	TCanvas* a9 = new TCanvas();
	h1_deconv_maxTS_ADCbiggerthan50->Draw();

	TCanvas* a10 = new TCanvas();
	h1_deconv_maxTS_ADCbiggerthan100->Draw();

	std::cout <<"Total number of hits on tracks accepted for the analysis= "<< nevents_analyzed <<'\n';
	std::cout <<"Percentage hits on tracks accepted for the analysis= "<< static_cast<double>(nevents_analyzed)*100/(static_cast<double>(nevents_analyzed)+static_cast<double>(nevents_NOT_analyzed)) <<" %\n";

}