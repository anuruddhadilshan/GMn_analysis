// Fit events using the APV function and extract t0, tau, and max ADC values and make 1D histograms.
// Insert analyzeOnTrackOrNot = 2 to analyze *both* the events that are ON and NOT ON tracks.
// Insert goosSignalShape = 1 to analyze only the events that are peaking in the 3rd and 4th time samples and ts1<ts2<ts3>ts4>ts5>ts6 / ts1<ts2<ts3<ts4>ts5>ts6. Inser 0 to analyze any event.

#include <iostream>
#include <string>

constexpr Int_t n_canvases_for_histos {3}; // No of histos to be drawn on canvases

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
		else return false;
	}
	if (maxTS == 3)
	{
		if (raw_6ADCsamples[0]<raw_6ADCsamples[1] && raw_6ADCsamples[1]<raw_6ADCsamples[2] && raw_6ADCsamples[2]<raw_6ADCsamples[3] && raw_6ADCsamples[3]>raw_6ADCsamples[4] && raw_6ADCsamples[4]>raw_6ADCsamples[5]) return true;
		else return false;
	}
	return false;
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

void fill_histos(const Double_t maxADC,const Double_t t0,const Double_t tau,TH1F* h1_maxADC,TH1F* h1_t0,TH1F* h1_tau)
{	 
	h1_maxADC->Fill(maxADC);
	h1_t0->Fill(t0);
	h1_tau->Fill(tau);
}

void make_pdf_for_displayed_events(const Int_t nevents_display,const Int_t nevents_recorded,TGraph* gr[nevents_display],TString output_file_name)
{
	TCanvas* c[nevents_recorded];

	TString pdffilename = output_file_name;
	TString openfilename = pdffilename+"(";
	TString closefilename = pdffilename+")";

	Double_t lmargin=0.15;
  	Double_t rmargin=0.15;
    Double_t bmargin=0.15;
    Double_t tmargin=0.09;

	for (Int_t ievent_display = 0; ievent_display < nevents_recorded; ++ievent_display)
	{
		TString event = "Event Number: "+std::to_string(ievent_display);
		c[ievent_display] = new TCanvas();
		c[ievent_display]->SetGrid();
		gr[ievent_display]->SetMarkerStyle(21);
		gr[ievent_display]->SetMarkerColor(1);
		gr[ievent_display]->SetTitle(event);
		gr[ievent_display]->GetXaxis()->SetTitle("Time (ns)");
		gr[ievent_display]->GetYaxis()->SetTitle("ADC");
		gr[ievent_display]->Draw("AP");
		//c[ievent_display]->BuildLegend();
		
		if(ievent_display == 0) c[ievent_display]->Print(openfilename,"Title:"+event);
		else if (ievent_display == nevents_recorded-1) c[ievent_display]->Print(closefilename,"Title:"+event);
		else c[ievent_display]->Print(pdffilename,"Title:"+event);
	}
}

void draw_hists_on_Canvas(TCanvas* c[n_canvases_for_histos],TH1F* h1_maxADC,TH1F* h1_t0,TH1F* h1_tau)
{
	for (Int_t i=0; i < n_canvases_for_histos; ++i) c[i] = new TCanvas();

	c[0]->cd(); h1_maxADC->Draw();
	c[1]->cd(); h1_t0->Draw();
	c[2]->cd(); h1_tau->Draw();
}

void make_pdf_for_histograms(TCanvas* c[n_canvases_for_histos],TString output_file_name)
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


void fit_test7(const Int_t run_number_int = 11590,const Double_t analyzeOnTrackOrNot = 1,const Int_t goodSignalShape = 0,const Int_t nevents_display = 10,const TString output_file_name = "fit_results_output")
{
	TChain* tchain_T = new TChain("T");
	const TString run_number_string = std::to_string(run_number_int);
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg0_0*.root");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg1_1*.root");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg2_2*.root");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg3_3*.root");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg4_4*.root");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg5_5*.root");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg6_6*.root");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_"+run_number_string+"_stream0_seg7_7*.root");

	
	// Disable all the unused branches to save time and enable only the branches used for analysis
	tchain_T->SetBranchStatus("*",false);
	tchain_T->SetBranchStatus("Ndata.bb.gem.m3.strip.ADCsamples",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.ADCsamples",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.nstripsfired",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.ontrack",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.isampmax",true);
	tchain_T->SetBranchStatus("bb.gem.m3.strip.ADCmax",true);
	

	//Define local the variables to store data from Tree
	Int_t ndata_stripADCsamples{0}; // Variable to hold the number of ADC samples for the event
	Double_t adcsamples[46080]{0.}; // Array to hold ADC smaples with max number of elements needed for a U-V GEM
	Double_t nstripsfired{0.};      // Number of strips fired for the event
	Double_t ontrack[7680]{0};      // The strip is on a track or not. Max number of elements needed for a U-V GEM
	Double_t maxtimesamp[7680] {0}; // Max time sample of the hit.
	Double_t maxADCvalue[7680] {0}; // Max ADC value of the hit.
	

	//Associate local variables with the Tree variables
	tchain_T->SetBranchAddress("Ndata.bb.gem.m3.strip.ADCsamples",&ndata_stripADCsamples);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ADCsamples",&adcsamples);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.nstripsfired",&nstripsfired);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ontrack",&ontrack);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.isampmax",&maxtimesamp);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ADCmax",&maxADCvalue);
	
	//Define histograms to store t0, tau and max_ADC and other parameters that we get fitting the 6 time sample data 
	TH1F* h1_maxADC = new TH1F("h1_maxADC","Fit Results for maxADCs - any ts peaking, no shape cut; ADC",3500,-500,3000); 
	TH1F* h1_t0 = new TH1F("h1_t0","Fit Result for t0 - any ts peaking, no shape cut; Time(ns)",400,-100,300);
	TH1F* h1_tau = new TH1F("h1_tau","Fit Results fotr APV Time Constant tau - any ts peaking, no shape cut; Time(ns) ",500,-200,300);
		
	Long64_t nevents = tchain_T->GetEntries();
	constexpr Int_t ntimesamples {6};
	Double_t timesamples[6] {12.0,36.0,60.0,84.0,108.0,132.0};
	Long64_t nevents_analyzed {0}; 
	Long64_t nevents_NOT_analyzed {0};
	TGraph* graph_min100_pos10[nevents_display];
	TGraph* graph_pos10_pos35[nevents_display];
	TGraph* graph_pos35_pos60[nevents_display];

	Int_t ievent_display_min100_pos10 {0};
	Int_t ievent_display_pos10_pos35 {0};
	Int_t ievent_display_pos35_pos60 {0};

	for(Long64_t nevent = 0; nevent < nevents; ++nevent)
	{
		tchain_T->GetEntry(nevent);

		Int_t istrip {0};
		
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
			
			//Copy the 6 ADC samples in to a an array of lenght 6.
			Double_t raw_6ADCsamples[ntimesamples] {0.};
			copyADCsamples(ientry, adcsamples, raw_6ADCsamples);

			// Decide whether the events have the "good signal shape" as described above.
			if (goodSignalShape==1)
			{
				if (!goodSignalShapeCut(raw_6ADCsamples,maxTS))
				{
					++nevents_NOT_analyzed; 
					++istrip;
					continue;
				}
			}

			//// The strips left after this step are the ones we are going to use for applying fits. ////

			// Apply the APV function fit to all the selected hits.
			TGraph* istrip_graph = new TGraph(ntimesamples,timesamples,raw_6ADCsamples);
			TF1* fitfunction = new TF1("apv_function",apv_function,0,150.0,3);

			apply_fit(istrip_graph,fitfunction,maxADCvalue[istrip]);
			
			Double_t maxADC {fitfunction->GetParameter(0)};
			Double_t t0 {fitfunction->GetParameter(1)};
			Double_t tau {fitfunction->GetParameter(2)};
			
			fill_histos(maxADC,t0,tau,h1_maxADC,h1_t0,h1_tau);
											
			//Copy the TGraphs for the events that are going to be displayed
			if((ievent_display_min100_pos10<nevents_display || ievent_display_pos10_pos35<nevents_display || ievent_display_pos35_pos60<nevents_display) && istrip%100==0)
			{	
				if(t0>=-100 && t0<=10 && ievent_display_min100_pos10<nevents_display) 
				{
					/*TString name_gr = "g"+std::to_string(ievent_display_min100_pos10);
					TString title_gr = "Event Number: "+std::to_string(ievent_display_min100_pos10);
					graph_min100_pos10[ievent_display_min100_pos10] = new TGraph(name_gr,title_gr);*/
					graph_min100_pos10[ievent_display_min100_pos10] = istrip_graph;
					++ievent_display_min100_pos10;
				}
				else if(t0>10 && t0<=35 && ievent_display_pos10_pos35<nevents_display)
				{
					/*TString name_gr = "g"+std::to_string(ievent_display_pos10_pos35);
					TString title_gr = "Event Number: "+std::to_string(ievent_display_pos10_pos35);
					graph_pos10_pos35[ievent_display_pos10_pos35] = new TGraph(name_gr,title_gr);*/
					graph_pos10_pos35[ievent_display_pos10_pos35] = istrip_graph;
					++ievent_display_pos10_pos35;
				}
				else if(t0>35 && t0<=60 && ievent_display_pos35_pos60<nevents_display)
				{
					/*TString name_gr = "g"+std::to_string(ievent_display_pos35_pos60);
					TString title_gr = "Event Number: "+std::to_string(ievent_display_pos35_pos60);
					graph_pos35_pos60[ievent_display_pos35_pos60] = new TGraph(name_gr,title_gr);*/
					graph_pos35_pos60[ievent_display_pos35_pos60] = istrip_graph;
					++ievent_display_pos35_pos60;
				}
			}

			++nevents_analyzed;
			++istrip;
		}

		if(nevent%100 == 0) std::cout <<"Current Event Number = " << nevent <<'\n';
	}

	//gStyle->SetOptFit(1);
	//h1_t0->Fit("gaus");
	//h1_tau->Fit("gaus");
	//h1_maxADC->Fit("landau");
		
	TCanvas* c_hist[n_canvases_for_histos];
	draw_hists_on_Canvas(c_hist,h1_maxADC,h1_t0,h1_tau);
	TString hist_pdf_filename = output_file_name+"_histo_output.pdf";
	make_pdf_for_histograms(c_hist,hist_pdf_filename);

	TString min100_pos10_evnts_filename = output_file_name+"_First_TimeBin.pdf";
	make_pdf_for_displayed_events(nevents_display,ievent_display_min100_pos10,graph_min100_pos10,min100_pos10_evnts_filename);

	TString pos10_pos35_evnts_filename = output_file_name+"_Second_TimeBin.pdf";
	make_pdf_for_displayed_events(nevents_display,ievent_display_pos10_pos35,graph_pos10_pos35,pos10_pos35_evnts_filename);

	TString pos35_pos60_evnts_filename = output_file_name+"_Third_TimeBin.pdf";
	make_pdf_for_displayed_events(nevents_display,ievent_display_pos35_pos60,graph_pos35_pos60,pos35_pos60_evnts_filename);
	
	std::cout <<"Total number of hits accepted for the analysis= "<< nevents_analyzed <<'\n';
	std::cout <<"Percentage hits accepted for the analysis= "<< static_cast<double>(nevents_analyzed)*100/(static_cast<double>(nevents_analyzed)+static_cast<double>(nevents_NOT_analyzed)) <<" %\n";
	std::cout <<"Number of events displayed with t0 between -100 and 10 ns = " << ievent_display_min100_pos10 <<'\n'; 
	std::cout <<"Number of events displayed with t0 between 10 and 35 ns = " << ievent_display_pos10_pos35 <<'\n'; 
	std::cout <<"Number of events displayed with t0 between 35 and 60 ns = " << ievent_display_pos35_pos60 <<'\n'; 
	
} 