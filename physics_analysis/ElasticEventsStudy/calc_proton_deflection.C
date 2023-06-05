// This scripts takes the results root files from the scripts "gmn_ana_dxdy.C" and "gmn_MC_ana_dxdy.C" and calculates the avg. proton deflection.
// Anuruddha Rathnayake - 04/24/2023

void calc_proton_deflection(const char* rootfile)
{
	// Load the root file and get the TTree object "T"
	TFile* anarootfile = new TFile(rootfile, "READ");
	//TTree* T = (TTree*)anarootfile->Get("T"); // Load the TTree, "T" which is made during the analysis by "gmn_ana_dxdy.C" script.

	// Define Canvas
	TCanvas* C1 = new TCanvas("C1", "Calculate prton deflection", 800, 800);

	// Load the 1D dx histogram.
	TH1D* h1_dx = (TH1D*)(anarootfile->Get("h1_dx_withcuts"));
	//TH1D* h1_dx = new TH1D("h1_dx","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	//T->Draw("sbs.hcal.dx>>h1_dx");
	h1_dx->Draw("HIST");
	//h1_dx->Fit("gaus");
	C1->Update();

	// Figure out a sub range to fit the proton peak.
	double p_peak_lowlimit{0.};
	double p_peak_highlimit{0.};
	TF1* p_peak_gaussfit;
	string fits_good = "n";
	bool b_fits_good = false;

	do 
	{
		std::cout << "Lower limit for the proton peak gaus fit: ";
		std::cin >> p_peak_lowlimit;
		std::cout << "Higher limit for the proton peak gaus fit: ";
		std::cin >> p_peak_highlimit;
		
		p_peak_gaussfit = new TF1("p_peak_gaussfit", "gaus", p_peak_lowlimit, p_peak_highlimit);
		
		h1_dx->Fit(p_peak_gaussfit, "R");
		p_peak_gaussfit->Draw("SAME");
		C1->Update();

		std::cout << "\nIs the sub-range fit good? y/n: ";
		std::cin >> fits_good;
		
		if (fits_good == "n" || fits_good == "no") 
		{
			std::cout << "Let us try again. \n ";
			std::cout << '\n';
			b_fits_good = false;
			h1_dx->Draw("HIST");
			C1->Update();
		}
		else if (fits_good == "y" || fits_good == "yes") 
		{
			std::cout << '\n';
			b_fits_good = true;
		}
	}
	while(!b_fits_good);

	std::cout << "The average proton deflection: " << TMath::Abs(p_peak_gaussfit->GetParameter(1)) << " meters\n";
}