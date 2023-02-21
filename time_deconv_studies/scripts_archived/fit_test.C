#include <iostream>


void copyADCsamples(Int_t ientry, const Double_t adcsamples[46080], Double_t raw_6ADCsamples[6])
{
	const Int_t ntimesamples {6};

	for(Int_t itimesamp = 0; itimesamp < ntimesamples; ++itimesamp)
	{
		raw_6ADCsamples[itimesamp] = adcsamples[ientry+itimesamp];
	}
}	

Double_t apv_function(Double_t* x, Double_t* par)
{	
	Double_t arg=0;

	if(par[2] != 0)
	{	
		arg = (x[0]-par[1])/par[2];
	}	
	Double_t fitval = par[0]*arg*TMath::Exp(-arg+1);
	return fitval;
}	

void apply_fit(TGraph* g, TF1* func)
{
	//func = new TF1("apv_func", apv_function, 0, 150.0, 3);
	func->SetParameters(500,0,56);
	func->SetParNames("MaxADC","t0","tau");
	g->Fit("apv_func","R");


}	

void fit_test(Double_t maxADC = 500.0, Double_t t0 = 0, Double_t tau = 56.0)
{
	TChain* tchain_T = new TChain("T");
	tchain_T->Add("/volatile/halla/sbs/adr/Rootfiles/replay_with_stripinfo/e1209019_fullreplay_11590_stream0_seg0_0.root");

	tchain_T->SetBranchStatus("*",false);
	tchain_T->SetBranchStatus("Ndata.bb.gem.m0.strip.ADCsamples",true);
	tchain_T->SetBranchStatus("bb.gem.m0.strip.ADCsamples",true);

	Int_t ndata_stripADCsamples{0};
	Double_t adcsamples[46080]{0.};

	tchain_T->SetBranchAddress("Ndata.bb.gem.m0.strip.ADCsamples",&ndata_stripADCsamples);
	tchain_T->SetBranchAddress("bb.gem.m0.strip.ADCsamples",&adcsamples);


	tchain_T->GetEntry(193);
	constexpr Int_t ntimesamples{6};
	const Int_t istrip = 2500;
	const Int_t ientry = ntimesamples*istrip;
	Double_t raw_6ADCsamples[ntimesamples] {0.};

	copyADCsamples(ientry, adcsamples, raw_6ADCsamples);

	Double_t timesample[6]{12.0,36.0,60.0,84.0,108.0,132.0};

	TGraph* g = new TGraph(ntimesamples, timesample, raw_6ADCsamples);

	TF1* func = new TF1("apv_func", apv_function, 0, 150.0, 3);
	//TF1* func = new TF1("apv_func","([0]*(x-[1])/[2])*TMath::Exp(-(x-[1])/[2]+1)",12.0,132.0);

	//TF1* func;
	//func->SetParameters(maxADC,t0,tau);
	//func->SetParNames("MaxADC","t0","tau");
	apply_fit(g,func);

	//g->Fit("apv_func","R");
	std::cout << func->GetChisquare() <<'\n';
	//g->Fit("gaus");
	g->Draw("A*");

}

