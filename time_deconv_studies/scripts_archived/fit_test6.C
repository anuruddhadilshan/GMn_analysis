// Selects events peaking at the 2nd,3rd,4th, and 5th timesamples
// Applies a "good signal shape cut" to only consider events that has a desired APV signal shape and analyze them seperately
// Fit those events using the APV function and extract t0, tau, and max ADC values and make 1D histograms. Also make (trigger_time - t0) plots.
// Apply Time Deconvolution to the above same events and find at what time sample the events are peaking.
// Insert analyzeOnTrackOrNot = 1 to analyze *only* the events that are ON tracks and "0" to analyze *only* the events that are NOT ON tracks.
// Insert analyzeOnTrackOrNot = 2 to analyze *both* the events that are ON and NOT ON tracks.

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

	// Applying Time Deconvolution to the time sample 0 *assuming* that the ADC value of the time sampel -1 and -2 to be zero.
	//deconv_6ADCsamples[0] = w_1*raw_6ADCsamples[0];

	// Applying Time Deconvolution to the time sample 1 *assuming* that the ADC value of the time sampel -1 to be zero.
	//deconv_6ADCsamples[1] = w_1*raw_6ADCsamples[1] + w_2*raw_6ADCsamples[0];
	
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
	Double_t fitval = std::max(0.0,par[0]*arg*TMath::Exp(-arg+1));
	return fitval;
}		

void apply_fit(TGraph* g, TF1* fitfunction, Double_t maxADCvalue)
{
	fitfunction->SetParameters(maxADCvalue,0.,56.);
	fitfunction->SetParNames("MaxADC","t0","tau");
	fitfunction->FixParameter(2,56.);
	
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
		if (max_ADC < sixADCsamples[ntimesample+1]) 
		{	
			max_timesample = ntimesample+1;
			max_ADC = sixADCsamples[ntimesample+1];
		}	
	}
}

Double_t findMinADC(const Double_t sixADCsamples[6])
{
	const Int_t ntimesamples = 6;
	Double_t min_ADC = sixADCsamples[0];

	for (Int_t ntimesample=0; ntimesample < ntimesamples-1; ++ntimesample) 
	{
		if (min_ADC > sixADCsamples[ntimesample+1]) 
		{	
			min_ADC = sixADCsamples[ntimesample+1];
		}	
	}

	return min_ADC;
}

Double_t adc_weighted_time {0.};	
void find_ADCweighted_time(const Double_t timesamples[6], const Double_t sixADCsamples[6])
{	
	const Int_t ntimesamples {6};
	Double_t sum_adcsamples {0.};
	Double_t sum_adcsampleval_times_timesampleval {0.};

	for (Int_t i=0; i<ntimesamples; ++i)
	{
		sum_adcsamples += sixADCsamples[i];
		sum_adcsampleval_times_timesampleval += sixADCsamples[i]*timesamples[i];
	}

	adc_weighted_time = sum_adcsampleval_times_timesampleval/sum_adcsamples;
}

//Find ADC weighted time using only the nearest neighbours.
Double_t adc_weighted_time_NN {0.};
void find_ADCweighted_time_NN(Double_t maxTS,const Double_t timesamples[6],const Double_t sixADCsamples[6])
{
	const Int_t ntimesamples {6};
	Double_t sum_adcsamples {0.};
	Double_t sum_adcsampleval_times_timesampleval {0.};
	const Int_t maxTS_int {static_cast<int>(maxTS)};

	if (maxTS_int==0) 
	{
		sum_adcsamples = sixADCsamples[0]+sixADCsamples[1]+sixADCsamples[2];
		sum_adcsampleval_times_timesampleval = sixADCsamples[0]*timesamples[0]+sixADCsamples[1]*timesamples[1]+sixADCsamples[2]*timesamples[2];
	}
	else if (maxTS_int==5) 
	{
		sum_adcsamples = sixADCsamples[3]+sixADCsamples[4]+sixADCsamples[5];
		sum_adcsampleval_times_timesampleval = sixADCsamples[3]*timesamples[3]+sixADCsamples[4]*timesamples[4]+sixADCsamples[5]*timesamples[5];
	}
	else
	{
		for (Int_t i=maxTS_int-1; i<=maxTS_int+1; ++i)
		{
			sum_adcsamples += sixADCsamples[i];
			sum_adcsampleval_times_timesampleval += sixADCsamples[i]*timesamples[i];
		}
	}

	adc_weighted_time_NN = sum_adcsampleval_times_timesampleval/sum_adcsamples;
}

void fill_ADCweighted_time_histos(Double_t maxTS,TH1F* h1_ADCweightedTime_TS2_peakingevnts,TH1F* h1_ADCweightedTime_TS3_peakingevnts,TH1F* h1_ADCweightedTime_TS4_peakingevnts,TH1F* h1_ADCweightedTime_TS5_peakingevnts,TH1F* h1_ADCweightedTimeNN_TS2_peakingevnts,TH1F* h1_ADCweightedTimeNN_TS3_peakingevnts,TH1F* h1_ADCweightedTimeNN_TS4_peakingevnts,TH1F* h1_ADCweightedTimeNN_TS5_peakingevnts)
{
	if (maxTS==1) 
	{
		h1_ADCweightedTime_TS2_peakingevnts->Fill(adc_weighted_time);
		h1_ADCweightedTimeNN_TS2_peakingevnts->Fill(adc_weighted_time_NN);
	}
	else if (maxTS==2)
	{
		h1_ADCweightedTime_TS3_peakingevnts->Fill(adc_weighted_time);
		h1_ADCweightedTimeNN_TS3_peakingevnts->Fill(adc_weighted_time_NN);
	}
	else if (maxTS==3)
	{ 
		h1_ADCweightedTime_TS4_peakingevnts->Fill(adc_weighted_time);
		h1_ADCweightedTimeNN_TS4_peakingevnts->Fill(adc_weighted_time_NN);
	}
	else if (maxTS==4)
	{
		h1_ADCweightedTime_TS5_peakingevnts->Fill(adc_weighted_time);
		h1_ADCweightedTimeNN_TS5_peakingevnts->Fill(adc_weighted_time_NN);
	}
}

Double_t correction_slope {0.24};
Double_t correction_intercept {-87.26};
void fill_histos_noDeconvolution(TF1* fitfunction, Double_t bbcal_trigtime,Double_t maxADC_from_tree, TH1F* h1_maxADC, TH1F* h1_t0, TH1F* h1_tau, TH1F* h1_t0_minus_trigtime, TH2F* h2_t0vsTrigTime, TH1F* h1_ADCweightedTime, TH1F* h1_ADCweightedTime_minus_t0, TH2F* h2_ADCweightedTimevsTrigTime, TH2F* h2_t0vsADCweightedTime,TH1F* h1_t0_minus_correctionTimesTrigTime,TH2F* h2_t0CorrectedvsTrigTime,TH2F* h2_tau_vs_maxADC)
{	
	Double_t maxADC {fitfunction->GetParameter(0)};
	Double_t t0 {fitfunction->GetParameter(1)};
	Double_t tau {fitfunction->GetParameter(2)};

	h1_maxADC->Fill(maxADC);
	h1_t0->Fill(t0);
	h1_tau->Fill(tau);
	//h1_trigtime->Fill(bbcal_trigtime);
	h1_t0_minus_trigtime->Fill(bbcal_trigtime - t0);
	h2_t0vsTrigTime->Fill(bbcal_trigtime, t0);
	h1_ADCweightedTime->Fill(adc_weighted_time);
	h1_ADCweightedTime_minus_t0->Fill(adc_weighted_time-t0);
	h2_ADCweightedTimevsTrigTime->Fill(bbcal_trigtime, adc_weighted_time);
	h2_t0vsADCweightedTime->Fill(adc_weighted_time,t0);
	Double_t t0_corrected {t0-(correction_slope*bbcal_trigtime+correction_intercept)};
 	h1_t0_minus_correctionTimesTrigTime->Fill(t0_corrected);
 	h2_t0CorrectedvsTrigTime->Fill(bbcal_trigtime,t0_corrected);	
 	//h2_tau_vs_maxADC->Fill(maxADC_from_tree,tau);
 	h2_tau_vs_maxADC->Fill(maxADC,tau);
}

void fill_histos_withDeconvolution(TH1F* h1_deconv_maxTS, TH1F* h1_deconv_maxADC, TH1F* h1_deconv_maxTS_ADCbiggerthan50, TH1F* h1_deconv_maxTS_ADCbiggerthan300, TH1F* h1_deconv_ADCweightedTime)
{
	h1_deconv_maxTS->Fill(max_timesample);
	h1_deconv_maxADC->Fill(max_ADC);
	if (max_ADC > 50) h1_deconv_maxTS_ADCbiggerthan50->Fill(max_timesample);
	if (max_ADC > 300) h1_deconv_maxTS_ADCbiggerthan300->Fill(max_timesample);
	h1_deconv_ADCweightedTime->Fill(adc_weighted_time);
}

void make_TMultiGraphs(Int_t nevents_display,TGraph* raw_gr[nevents_display],TGraph* deconv_gr[nevents_display],TMultiGraph* mul_gr[nevents_display])
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

void make_pdf_for_displayed_events(Int_t nevents_display,TMultiGraph* gr[nevents_display],TString output_file_name)
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

void divideTheCanvasIntoTwoAndDrawTheHistograms(TCanvas* c, TH1F* h1, TH1F* h2)
{
	c->Divide(2,1);
	c->cd(1);
	h1->Draw();
	c->cd(2);
	h2->Draw();
}

void divideTheCanvasIntoTwoAndDrawTheHistograms(TCanvas* c, TH2F* h1, TH2F* h2)
{
	c->Divide(2,1);
	c->cd(1);
	h1->Draw("COLZ");
	c->cd(2);
	h2->Draw("COLZ");
}

const Int_t n_canvases_for_histos {21};
void make_pdf_for_histograms(TCanvas* c[n_canvases_for_histos], TString output_file_name)
{
	TString pdffilename = "analysis_output_histograms_"+output_file_name;
	TString openfilename = pdffilename+"(";
	TString closefilename = pdffilename+")";

	Double_t lmargin=0.15;
  	Double_t rmargin=0.15;
    Double_t bmargin=0.15;
    Double_t tmargin=0.09;

	for (Int_t icanvas = 0; icanvas < n_canvases_for_histos; ++icanvas)
	{
		if(icanvas == 0) c[icanvas]->Print(openfilename);
		else	if (icanvas == n_canvases_for_histos-1) c[icanvas]->Print(closefilename);
		else c[icanvas]->Print(pdffilename);
	}
}


void fit_test6(const Int_t run_number_int = 11590, Double_t analyzeOnTrackOrNot = 1, Int_t nevents_display = 10, const TString output_file_name = "output_with_fits.pdf")
{
	TChain* tchain_T = new TChain("T");
	const TString run_number_string = std::to_string(run_number_int);
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg0_0*.root");
	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg1_1*.root");
	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg2_2*.root");
	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg3_3*.root");
	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg4_4*.root");
	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg5_5*.root");
	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg6_6*.root");
	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg7_7*.root");
	//tchain_T->Add("/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_11590_stream0_seg0_9.root");

	//tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/tests_3/e1209019_replayed_11594_stream0_seg0_0_firstevent0_nevent1000.root");
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

	//Define histograms to store t0, tau and max_ADC and other parameters that we get fitting the 6 time sample data 
	TH1F* h1_maxADC = new TH1F("h1_maxADC","Fit Results for maxADCs - any ts peaking, no shape cut; ADC",3500,-500,3000); 
	TH1F* h1_maxADC_goodshape = new TH1F("h1_maxADC_goodshape","Fit Results for maxADCs - 3,4 ts peaking, with shape cut; ADC",3500,-500,3000); 
	TH1F* h1_t0 = new TH1F("h1_t0","Fit Result for t0 - any ts peaking, no shape cut; Time(ns)",400,-100,300);
	TH1F* h1_t0_goodshape = new TH1F("h1_t0_goodshape","Fit Results for t0 - 3,4 ts peaking, with shape cut; Time(ns)",400,-100,300);
	TH1F* h1_tau = new TH1F("h1_tau","Fit Results fotr APV Time Constant tau - any ts peaking, no shape cut; Time(ns) ",500,-200,300);
	TH1F* h1_tau_goodshape = new TH1F("h1_tau_goodshape","Fit Results for APV Time Constant tau - 3,4 ts peaking, with shape cut; Time(ns) ",500,-200,300);
	TH1F* h1_trigtime = new TH1F("h1_trigtime", "BBCal Trigger time distribution; Time(ns)", 100, 320 ,420);
	TH1F* h1_maxADC_TS_before_deconv = new TH1F("h1_maxADC_TS_before_deconv","Maximum ADC Time Sample Distribution BEFORE Deconvolution",6,-0.5,5.5);
	TH1F* h1_t0_minus_trigtime = new TH1F("h1_t0_minus_trigtime","TriggerTime-t0 : any ts peaking, no shape cut ; Time(ns)",700,-100,600);
	TH1F* h1_t0_minus_trigtime_goodshape = new TH1F("h1_t0_minus_trigtime_goodshape","TriggerTime-t0 : 3,4 ts peaking, with shape cut ; Time(ns)",700,-100,600);
	TH2F* h2_t0vsTrigTime = new TH2F("h2_t0vsTrigTime","t0 and Trigger Time correlation - any ts peaking, no shape cut; BBCal Trigger Time (ns); t0 (ns)",100,320,420,100,-50,50);
	TH2F* h2_t0vsTrigTime_goodshape = new TH2F("h2_t0vsTrigTime_goodshape","t0 and Trigger Time correlation - 3,4 ts peaking, with shape cut; BBCal Trigger Time (ns); t0 (ns)",100,320,420,100,-50,50);
	TH1F* h1_t0_minus_correctionTimesTrigTime = new TH1F("h1_t0_minus_correctionTimesTrigTime","t0-CorrectedTrigTime : any ts peaking, no shape cut ; Time(ns)",300,-100,200);
	TH1F* h1_t0_minus_correctionTimesTrigTime_goodshape = new TH1F("h1_t0_minus_correctionTimesTrigTime_goodshape","t0-CorrectedTrigTime : 3,4 ts peaking, with shape cut ; Time(ns)",300,-100,200);
	TH1F* h1_ADCweightedTime = new TH1F("h1_ADCweightedTime", "ADC Weighted Time Before Time Deconvolution - any ts peaking, no shape cut; Time(ns)",151,0,150);
	TH1F* h1_ADCweightedTime_goodshape = new TH1F("h1_ADCweightedTime_goodshape", "ADC Weighted Time Before Time Deconvolution - 3,4 ts peaking, with shape cut; Time(ns)",151,0,150);
	TH1F* h1_ADCweightedTime_minus_t0 = new TH1F("h1_ADCweightedTime_minus_t0", "ADCWeightedTime-t0 - any ts peaking, no shape cut; Time(ns)",151,0,150);
	TH1F* h1_ADCweightedTime_minus_t0_goodshape = new TH1F("h1_ADCweightedTime_minus_t0_goodshape", "ADCWeightedTime-t0 - 3,4 ts peaking, with shape cut ; Time(ns)",151,0,150);
	TH2F* h2_ADCweightedTimevsTrigTime = new TH2F("h2_ADCweightedTimevsTrigTime","ADC weighted time and Trigger Time correlation - any ts peaking, no shape cut; BBCal Trigger Time (ns); ADC weighted time (ns)",100,320,420,151,0,150);
	TH2F* h2_ADCweightedTimevsTrigTime_goodshape = new TH2F("h2_ADCweightedTimevsTrigTime_goodshape","ADC Weighted Time and Trigger Time correlation - 3,4 ts peaking, with shape cut; BBCal Trigger Time (ns); ADC weighted time (ns)",100,320,420,151,0,150);
	TH2F* h2_t0vsADCweightedTime = new TH2F("h2_t0vsADCweightedTime","t0 and ADC Weighted Time correlation - any ts peaking, no shape cut;ADC weighted time (ns);t0 (ns)",151,0,150,200,-100,100);
	TH2F* h2_t0vsADCweightedTime_goodshape = new TH2F("h2_t0vsADCweightedTime_goodshape","t0 and ADC Weighted Time correlation - 3,4 ts peaking, with shape cut;ADC weighted time (ns);t0 (ns)",151,0,150,200,-100,100);
	TH2F* h2_t0CorrectedvsTrigTime = new TH2F("h2_t0CorrectedvsTrigTime","t0 Corrected and Trigger Time correlation - any ts peaking, no shape cut; BBCal Trigger Time (ns); t0 (ns)",100,320,420,100,-50,50);
	TH2F* h2_t0CorrectedvsTrigTime_goodshape = new TH2F("h2_t0CorrectedvsTrigTime_goodshape","t0 Corrected and Trigger Time correlation - 3,4 ts peaking, with shape cut; BBCal Trigger Time (ns); t0 (ns)",100,320,420,100,-50,50);
	TH2F* h2_tau_vs_maxADC = new TH2F("h2_tau_vs_maxADC","Tau vs maxADC - any ts peaking, no shape cut;Max ADC;tau(ns)",3500,-500,3000,400,-100,300);
	TH2F* h2_tau_vs_maxAD_goodshape = new TH2F("h2_tau_vs_maxADC_goodshape","Tau vs maxADC - any ts peaking, 3,4 ts peaking with shape cut;Max ADC;tau(ns)",3500,-500,3000,400,-100,300);

	TH1F* h1_ADCweightedTime_TS2_peakingevnts = new TH1F("h1_ADCweightedTime_TS2_PeakingEvents","ADC Weighted Time for Events Peaking in TS-2",151,0,150);
	TH1F* h1_ADCweightedTime_TS3_peakingevnts = new TH1F("h1_ADCweightedTime_TS3_PeakingEvents","ADC Weighted Time for Events Peaking in TS-3",151,0,150);
	TH1F* h1_ADCweightedTime_TS4_peakingevnts = new TH1F("h1_ADCweightedTime_TS4_PeakingEvents","ADC Weighted Time for Events Peaking in TS-4",151,0,150);
	TH1F* h1_ADCweightedTime_TS5_peakingevnts = new TH1F("h1_ADCweightedTime_TS5_PeakingEvents","ADC Weighted Time for Events Peaking in TS-5",151,0,150);
	TH1F* h1_ADCweightedTimeNN_TS2_peakingevnts = new TH1F("h1_ADCweightedTimeNN_TS2_PeakingEvents","ADC Weighted Time Using Nearest Neighbours for Events Peaking in TS-2",151,0,150);
	TH1F* h1_ADCweightedTimeNN_TS3_peakingevnts = new TH1F("h1_ADCweightedTimeNN_TS3_PeakingEvents","ADC Weighted Time Using Nearest Neighbours for Events Peaking in TS-3",151,0,150);
	TH1F* h1_ADCweightedTimeNN_TS4_peakingevnts = new TH1F("h1_ADCweightedTimeNN_TS4_PeakingEvents","ADC Weighted Time Using Nearest Neighbours for Events Peaking in TS-4",151,0,150);
	TH1F* h1_ADCweightedTimeNN_TS5_peakingevnts = new TH1F("h1_ADCweightedTimeNN_TS5_PeakingEvents","ADC Weighted Time Using Nearest Neighbours for Events Peaking in TS-5",151,0,150);

	TH1F* h1_deconv_maxTS = new TH1F("h1_deconv_maxTS","Max Time Sample after Time Decon - any ts peaking, no shape cut; Time(ns)",6,-0.5,5.5);
	TH1F* h1_deconv_maxTS_goodshape = new TH1F("h1_deconv_maxTS_goodshape","Max Time Sample after Time Decon - 3,4 ts peaking, with shape cut; Time(ns)",6,-0.5,5.5);
	TH1F* h1_deconv_maxADC = new TH1F("h1_deconv_maxADC","Max Time Sample ADC after Time Decon - any ts peaking, no shape cut; ADC",3500,-500,3000);
	TH1F* h1_deconv_maxADC_goodshape = new TH1F("h1_deconv_maxADC_goodshape","Max Time Sample ADC after Time Decon - 3,4 ts peaking, with shape cut; ADC",3500,-500,3000);
	TH1F* h1_deconv_maxTS_ADCbiggerthan50 = new TH1F("h1_deconv_maxTS_ADCbiggerthan50","Tims Samples with Max ADC > 50 after Time Decon - any ts peaking, no shape cut ; Time(ns)",6,-0.5,5.5);
	TH1F* h1_deconv_maxTS_ADCbiggerthan50_goodshape = new TH1F("h1_deconv_maxTS_ADCbiggerthan50_goodshape","Tims Samples with Max ADC > 50 after Time Decon - 3,4 ts peaking, with shape cut ; Time(ns)",6,-0.5,5.5);
	TH1F* h1_deconv_maxTS_ADCbiggerthan300 = new TH1F("h1_deconv_maxTS_ADCbiggerthan300","Time Samples with Max ADC > 300 after Time Decon -  any ts peaking, no shape cut; Time(ns)",6,-0.5,5.5);
	TH1F* h1_deconv_maxTS_ADCbiggerthan300_goodshape = new TH1F("h1_deconv_maxTS_ADCbiggerthan300_goodshape","Time Samples with Max ADC > 300 after Time Decon -  3,4 ts peaking, with shape cut ; Time(ns)",6,-0.5,5.5);
	TH1F* h1_deconv_ADCweightedTime = new TH1F("h1_deconv_ADCweightedTime","ADC Weighted Time after Time Decon - any ts peaking, no shape cut; Time(ns)",151,0,150);
	TH1F* h1_deconv_ADCweightedTime_goodshape = new TH1F("h1_deconv_ADCweightedTime_goodshape","ADC Weighted Time after Time Decon - 3,4 ts peaking, with shape cut; Time(ns)",151,0,150);

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
		h1_trigtime->Fill(bbcal_trigtime);

		for(Int_t ientry = 0; ientry < ndata_stripADCsamples; ientry += ntimesamples)
		{
			// Decide whether to analyze only the strips that are not on tracks, on track or both on track and not on track events.
			if (!(analyzeOnTrackOrNot==2))
			{			
				if (!(ontrack[istrip] == analyzeOnTrackOrNot))
				{	
					++nevents_NOT_analyzed; 
					++istrip;
					continue;
				}
			}

			Double_t maxTS {maxtimesamp[istrip]};
			/*bool first_or_sixth_TimeSample_peaking {(maxTS==0) || (maxTS==5)};	
			if(first_or_sixth_TimeSample_peaking)
			{
				++nevents_NOT_analyzed; 
				++istrip;
				continue;
			}*/

			//Copy the 6 ADC samples in to a an array of lenght 6.
			Double_t raw_6ADCsamples[ntimesamples] {0.};
			copyADCsamples(ientry, adcsamples, raw_6ADCsamples);

			/*if(findMinADC(raw_6ADCsamples)<50.0)
			{
				++nevents_NOT_analyzed;
				++istrip;
				continue;
			}*/
			//// The strips left after this step are the ones we are going to use for applying fits. ////

			bool third_or_fourth_TimeSample_peaking {(maxTS==2) || (maxTS==3)};
			bool good_signal_shape {goodSignalShapeCut(raw_6ADCsamples, maxTS)};
			bool sigShape_and_3rdOr4thTS_peaking {third_or_fourth_TimeSample_peaking && good_signal_shape};

			h1_maxADC_TS_before_deconv->Fill(maxTS);

			// Apply the APV function fit to all the selected hits.
			TGraph* istrip_graph = new TGraph(ntimesamples, timesamples, raw_6ADCsamples);
			TF1* fitfunction = new TF1("apv_function", apv_function, 0, 150.0, 3);
			apply_fit(istrip_graph, fitfunction, maxADCvalue[istrip]);
			find_ADCweighted_time(timesamples,raw_6ADCsamples);
			find_ADCweighted_time_NN(maxTS,timesamples,raw_6ADCsamples);
			fill_ADCweighted_time_histos(maxTS,h1_ADCweightedTime_TS2_peakingevnts,h1_ADCweightedTime_TS3_peakingevnts,h1_ADCweightedTime_TS4_peakingevnts,h1_ADCweightedTime_TS5_peakingevnts,h1_ADCweightedTimeNN_TS2_peakingevnts,h1_ADCweightedTimeNN_TS3_peakingevnts,h1_ADCweightedTimeNN_TS4_peakingevnts,h1_ADCweightedTimeNN_TS5_peakingevnts);
			fill_histos_noDeconvolution(fitfunction,bbcal_trigtime,maxADCvalue[istrip],h1_maxADC,h1_t0,h1_tau,h1_t0_minus_trigtime,h2_t0vsTrigTime,h1_ADCweightedTime,h1_ADCweightedTime_minus_t0,h2_ADCweightedTimevsTrigTime,h2_t0vsADCweightedTime,h1_t0_minus_correctionTimesTrigTime,h2_t0CorrectedvsTrigTime,h2_tau_vs_maxADC);
			if (sigShape_and_3rdOr4thTS_peaking) fill_histos_noDeconvolution(fitfunction,bbcal_trigtime,maxADCvalue[istrip],h1_maxADC_goodshape,h1_t0_goodshape,h1_tau_goodshape,h1_t0_minus_trigtime_goodshape,h2_t0vsTrigTime_goodshape,h1_ADCweightedTime_goodshape,h1_ADCweightedTime_minus_t0_goodshape,h2_ADCweightedTimevsTrigTime_goodshape,h2_t0vsADCweightedTime_goodshape,h1_t0_minus_correctionTimesTrigTime_goodshape,h2_t0CorrectedvsTrigTime_goodshape,h2_tau_vs_maxAD_goodshape);

			// Apply time deconvolution to all the selected hits.
			Double_t deconv_6ADCsamples[ntimesamples] {0.};
			deconvolutionADCsamples(raw_6ADCsamples, deconv_6ADCsamples);
			findmaxtimesample_and_ADC(deconv_6ADCsamples);
			find_ADCweighted_time(timesamples,deconv_6ADCsamples);
			fill_histos_withDeconvolution(h1_deconv_maxTS,h1_deconv_maxADC,h1_deconv_maxTS_ADCbiggerthan50,h1_deconv_maxTS_ADCbiggerthan300,h1_deconv_ADCweightedTime);
			if (sigShape_and_3rdOr4thTS_peaking) fill_histos_withDeconvolution(h1_deconv_maxTS_goodshape,h1_deconv_maxADC_goodshape,h1_deconv_maxTS_ADCbiggerthan50_goodshape,h1_deconv_maxTS_ADCbiggerthan300_goodshape,h1_deconv_ADCweightedTime_goodshape);
								
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
	make_pdf_for_displayed_events(nevents_display, multi_graph, output_file_name);

	gStyle->SetOptFit(1);
	//h1_t0->Fit("gaus");
	//h1_t0_goodshape->Fit("gaus");
	h1_tau->Fit("gaus");
	h1_tau_goodshape->Fit("gaus");
	h1_maxADC->Fit("landau");
	h1_maxADC_goodshape->Fit("landau");

	TF1 *straight_line = new TF1("straight_line","[0]*x+[1]",350,385);
	straight_line->SetParameters(0.3,-200);
	straight_line->SetLineColor(kRed);
	h2_t0vsTrigTime->Fit(straight_line);
	h2_t0vsTrigTime_goodshape->Fit(straight_line);
	
	TCanvas* c_hist[n_canvases_for_histos];
	for (Int_t i=0; i < n_canvases_for_histos; ++i) c_hist[i] = new TCanvas();
	
	//Divide The Canvas Into Two And Draw The Histograms	
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[0],h1_maxADC,h1_maxADC_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[1],h1_t0,h1_t0_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[2],h1_tau,h1_tau_goodshape);
	c_hist[3]->cd(); h1_trigtime->Draw();
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[4],h1_t0_minus_trigtime,h1_t0_minus_trigtime_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[5],h2_t0vsTrigTime,h2_t0vsTrigTime_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[6],h2_t0CorrectedvsTrigTime,h2_t0CorrectedvsTrigTime_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[7],h1_ADCweightedTime,h1_ADCweightedTime_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[8],h1_ADCweightedTime_minus_t0,h1_ADCweightedTime_minus_t0_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[9],h2_ADCweightedTimevsTrigTime,h2_ADCweightedTimevsTrigTime_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[10],h2_t0vsADCweightedTime,h2_t0vsADCweightedTime_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[11],h1_t0_minus_correctionTimesTrigTime,h1_t0_minus_correctionTimesTrigTime_goodshape);
	c_hist[12]->Divide(2,2); 
	c_hist[12]->cd(1); h1_ADCweightedTime_TS2_peakingevnts->Draw();
	c_hist[12]->cd(2); h1_ADCweightedTime_TS3_peakingevnts->Draw();
	c_hist[12]->cd(3); h1_ADCweightedTime_TS4_peakingevnts->Draw();
	c_hist[12]->cd(4); h1_ADCweightedTime_TS5_peakingevnts->Draw();
	c_hist[13]->Divide(2,2); 
	c_hist[13]->cd(1); h1_ADCweightedTimeNN_TS2_peakingevnts->Draw();
	c_hist[13]->cd(2); h1_ADCweightedTimeNN_TS3_peakingevnts->Draw();
	c_hist[13]->cd(3); h1_ADCweightedTimeNN_TS4_peakingevnts->Draw();
	c_hist[13]->cd(4); h1_ADCweightedTimeNN_TS5_peakingevnts->Draw();
	c_hist[14]->cd(); h1_maxADC_TS_before_deconv->Draw();
	//divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[15],h1_deconv_maxTS,h1_deconv_maxTS_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[15],h1_maxADC_TS_before_deconv,h1_deconv_maxTS);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[16],h1_deconv_maxTS_ADCbiggerthan50,h1_deconv_maxTS_ADCbiggerthan50_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[17],h1_deconv_maxTS_ADCbiggerthan300,h1_deconv_maxTS_ADCbiggerthan300_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[18],h1_deconv_maxADC,h1_deconv_maxADC_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[19],h1_deconv_ADCweightedTime,h1_deconv_ADCweightedTime_goodshape);
	divideTheCanvasIntoTwoAndDrawTheHistograms(c_hist[20],h2_tau_vs_maxADC,h2_tau_vs_maxAD_goodshape);
	
	make_pdf_for_histograms(c_hist, output_file_name);

	std::cout <<"Total number of hits accepted for the analysis= "<< nevents_analyzed <<'\n';
	std::cout <<"Percentage hits accepted for the analysis= "<< static_cast<double>(nevents_analyzed)*100/(static_cast<double>(nevents_analyzed)+static_cast<double>(nevents_NOT_analyzed)) <<" %\n";

}