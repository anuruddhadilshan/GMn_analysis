// Anuruddha Rathnayake - 4/20/2023
// ROOT macro that compares dx and dy plots of simulated and real data.



void compare_data_and_sim(const char* realdat_rootfile, const char* simudat_rootfile)
{
	
	// Load the root files and get the TTree object "T"
	TFile* real = new TFile(realdat_rootfile, "READ");
	TTree* T_real = (TTree*)real->Get("T"); // Load the TTree, "T" which is made during the analysis by "gmn_ana_dxdy.C" script.
	TFile* simu = new TFile(simudat_rootfile, "READ");
	TTree* T_simu = (TTree*)simu->Get("T"); // Load the TTree, "T" which is made during the analysis by "gmn_ana_dxdy.C" script.

	// Define the histograms.
	TH1D* h1_realdat_dx = new TH1D("h1_realdat_dx", "dx; X_{HCal}-X_{expected} (m)", 1000, -4, 6);
	TH1D* h1_realdat_dy = new TH1D("h1_realdat_dy", "dy; Y_{HCal}-Y_{expected} (m)", 400, -2, 2);
	TH1D* h1_simudat_dx = new TH1D("h1_simudat_dx", "dx; X_{HCal}-X_{expected} (m)", 1000, -4, 6);
	TH1D* h1_simudat_dy = new TH1D("h1_simudat_dy", "dy; Y_{HCal}-Y_{expected} (m)", 400, -2, 2);

	//T_real->Draw("sbs.hcal.dx>>h1_data_dx");
	//T_real->Draw("sbs.hcal.dy>>h1_data_dy");

	// Fill the real data histograms.
	T_real->SetBranchStatus("*", 0);
	T_real->SetBranchStatus("sbs.hcal.dx", 1);
	T_real->SetBranchStatus("sbs.hcal.dy", 1);
	double dx_real{0.}; 
	double dy_real{0.};
	T_real->SetBranchAddress("sbs.hcal.dx", &dx_real);
	T_real->SetBranchAddress("sbs.hcal.dy", &dy_real );
	long nentries_realdat = T_real->GetEntries();
	
	for(long i = 0; i < nentries_realdat; i++)
	{
		T_real->GetEntry(i);
		h1_realdat_dx->Fill(dx_real);
		h1_realdat_dy->Fill(dy_real);
	}
	//

	// Fill simulated data histograms.
	T_simu->SetBranchStatus("*", 0);
	T_simu->SetBranchStatus("sbs.hcal.dx", 1);
	T_simu->SetBranchStatus("sbs.hcal.dy", 1);
	T_simu->SetBranchStatus("mc.weight", 1);
	double dx_simu{0.};
	double dy_simu{0.};
	double weight{0.};
	T_simu->SetBranchAddress("sbs.hcal.dx", &dx_simu);
	T_simu->SetBranchAddress("sbs.hcal.dy", &dy_simu);
	T_simu->SetBranchAddress("mc.weight", &weight);
	long nentries_simudat = T_simu->GetEntries();
	
	for(long i = 0; i < nentries_simudat; i++)
	{
		T_simu->GetEntry(i);
		h1_simudat_dx->Fill(dx_simu, weight);
		h1_simudat_dy->Fill(dy_simu, weight);
	}
	//

	// Re-scalling all the histograms so they can be compared with each other.
	h1_realdat_dx->Scale(1./h1_realdat_dx->Integral(), "width");
	h1_realdat_dy->Scale(1./h1_realdat_dy->Integral(), "width");
	h1_simudat_dx->Scale(1./h1_simudat_dx->Integral(), "width");
	h1_simudat_dy->Scale(1./h1_simudat_dy->Integral(), "width");

	TCanvas* C1 = new TCanvas("C1", "Compare dx distributions", 600, 600);
	C1->SetTitle("Compare dx distributions");
	h1_realdat_dx->SetLineColor(4);
	h1_simudat_dx->SetLineColor(8);
	TLegend* dxplot_legend = new TLegend(0.66,0.7,0.9,0.9);
	dxplot_legend->AddEntry(h1_realdat_dx, "Experiment data", "l");
	dxplot_legend->AddEntry(h1_simudat_dx, "Simulated data", "l");
	h1_realdat_dx->Draw("HIST");
	h1_simudat_dx->Draw("HIST+SAME");
	dxplot_legend->Draw();
	C1->Update();

	TCanvas* C2 = new TCanvas("C2", "Compare dy distributions", 600, 600);
	C2->SetTitle("Compare dy distributions");
	h1_realdat_dy->SetLineColor(4);
	h1_simudat_dy->SetLineColor(8);
	TLegend* dyplot_legend = new TLegend(0.66,0.7,0.9,0.9);
	dyplot_legend->AddEntry(h1_realdat_dy, "Experiment data", "l");
	dyplot_legend->AddEntry(h1_simudat_dy, "Simulated data", "l");
	h1_realdat_dy->Draw("HIST");
	h1_simudat_dy->Draw("HIST+SAME");
	dyplot_legend->Draw();
	C2->Update();

	string save_plots = "no";

	std::cout << "Save the plots (yes/no): ";
	std::cin >> save_plots;

	if(save_plots == "yes" || save_plots == "y")
	{
		std::cout << "File name: ";
		TString plot_file_name;
		std::cin >> plot_file_name;

		C1->SaveAs(Form("%s_dx.png", plot_file_name.Data()));
		C2->SaveAs(Form("%s_dy.png", plot_file_name.Data()));
	}
	
}

