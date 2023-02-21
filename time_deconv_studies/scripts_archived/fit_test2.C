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
bool goodSignalShapeCut(Double_t raw_6ADCsamples[6], Double_t maxTS)
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

void apply_fit(TGraph* g, TF1* fitfunction)
{
	fitfunction->SetParameters(500.,0.,56.);
	fitfunction->SetParNames("MaxADC","t0","tau");
	
	g->Fit("apv_function","RQ");
}

void fill_1D_histos(TF1* fitfunction, TH1F* h1_maxADC, TH1F* h1_t0, TH1F* h1_tau)
{	
	h1_maxADC->Fill(fitfunction->GetParameter(0));
	h1_t0->Fill(fitfunction->GetParameter(1));
	h1_tau->Fill(fitfunction->GetParameter(2));
}
void edit_TGraphs(Int_t nevents_display, TGraph* gr[nevents_display])
{
	for(Int_t ievent_display = 0; ievent_display < nevents_display; ++ievent_display)
	{	
		TString event = std::to_string(ievent_display);
		gr[ievent_display]->SetTitle("ADC signal No Deconvolution Evnt: "+event);
		gr[ievent_display]->SetMarkerStyle(21);
		gr[ievent_display]->SetMarkerColor(1);
		gr[ievent_display]->GetXaxis()->SetTitle("Time (ns)");
		gr[ievent_display]->GetYaxis()->SetTitle("ADC");
	}
}

void make_canvases_and_pdfs(Int_t nevents_display, TGraph* gr[nevents_display], TString output_file_name)
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

void fit_test2(const Int_t run_number_int = 11590, Int_t nevents_display = 10, const TString output_file_name = "output_with_fits.pdf")
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


	//Define local the variables to store data from Tree
	Int_t ndata_stripADCsamples{0}; // Variable to hold the number of ADC samples for the event
	Double_t adcsamples[46080]{0.}; // Array to hold ADC smaples with max number of elements needed for a U-V GEM
	Double_t nstripsfired{0.};      // Number of strips fired for the event
	Double_t ontrack[7680]{0};      // The strip is on a track or not. Max number of elements needed for a U-V GEM
	Double_t maxtimesamp[7680] {0}; // Max time sample of the hit.

	//Associate local variables with the Tree variables
	tchain_T->SetBranchAddress("Ndata.bb.gem.m3.strip.ADCsamples",&ndata_stripADCsamples);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ADCsamples",&adcsamples);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.nstripsfired",&nstripsfired);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.ontrack",&ontrack);
	tchain_T->SetBranchAddress("bb.gem.m3.strip.isampmax",&maxtimesamp);
	

	//Define 1D histograms to store t0, tau and max_ADC fit parametrs
	TH1F* h1_maxADC = new TH1F("h1_maxADC","Fit Result Distribution for maxADCs",3500,-500,3000);
	TH1F* h1_t0 = new TH1F("h1_t0","Fit Result Distribution for t0",400,-100,300);
	TH1F* h1_tau = new TH1F("h1_tau","Fit Result Distribution for APV Time Constant tau",500,-200,300); 

	Long64_t nevents = tchain_T->GetEntries();
	constexpr Int_t ntimesamples {6};
	Double_t timesamples[6] {12.0,36.0,60.0,84.0,108.0,132.0};
	Long64_t nevents_analyzed {0}; 
	Long64_t nevents_NOT_analyzed {0};
	TGraph* gr[nevents_display];
	Int_t ievent_display {0};


	for(Long64_t nevent = 0; nevent < nevents; ++nevent)
	{
		tchain_T->GetEntry(nevent);

		Int_t istrip {0};

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
			//The strips left after this step are the ones we are going to use for applying fits.

			TGraph* istrip_graph = new TGraph(ntimesamples, timesamples, raw_6ADCsamples);
			TF1* fitfunction = new TF1("apv_function", apv_function, 0, 150.0, 3);
			apply_fit(istrip_graph, fitfunction);
			
			//Fill the 1D histograms
			fill_1D_histos(fitfunction, h1_maxADC, h1_t0, h1_tau);

			//std::cout << fitfunction->GetChisquare() <<'\n';
			
			if(ievent_display < nevents_display && istrip%100==0)
			{	
				gr[ievent_display] = istrip_graph;
				apply_fit(istrip_graph, fitfunction);
				++ievent_display;
			}

			++nevents_analyzed;
			++istrip;
		}

		if(nevent%100 == 0) std::cout <<"Current Event Number = " << nevent <<'\n';
	}

	edit_TGraphs(nevents_display, gr);
	make_canvases_and_pdfs(nevents_display, gr, output_file_name);

	TCanvas* a1 = new TCanvas();
	h1_maxADC->Draw();

	TCanvas* a2 = new TCanvas();
	h1_t0->Draw();

	TCanvas* a3 = new TCanvas();
	h1_tau->Draw();

	std::cout <<"Total number of hits on tracks accepted for the analysis= "<< nevents_analyzed <<'\n';
	std::cout <<"Percentage hits on tracks accepted for the analysis= "<< static_cast<double>(nevents_analyzed)*100/(static_cast<double>(nevents_analyzed)+static_cast<double>(nevents_NOT_analyzed)) <<" %\n";

}