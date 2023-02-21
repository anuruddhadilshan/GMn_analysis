// Macro to fit the dx plots with two gausians for n and p peaks and a polynormial for inelastic background.
// Only to be used for LD2 data.



void fit_dx_plots(const char* rootfile, TString  outputname = "sbsx_sbsfield_ypercent_LD2", const int n_polynormial = 4)
{
	// Load the root file and get the histogram "h1_dx_withcuts"
	TFile* anarootfile = new TFile(rootfile, "READ");
	TTree* T = (TTree*)anarootfile->Get("T"); // Load the TTree, "T" which is made during the analysis by "gmn_ana_dxdy.C" script.

	// Draw the root file on a Canvas. 
	TCanvas* C1 = new TCanvas("C1", outputname, 800, 800);
	TH1D* h1_dx_1 = new TH1D("h1_dx_1","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	T->Draw("sbs.hcal.dx>>h1_dx_1");
	C1->Update(); // Updates the Canvas and displays it.

	// Select the sub-ranges for the neutron and proton gaus fits and the background fit (n_polynormal'th order polynomial).
	double n_gausfit_lowlimit{0.};
	double n_gausfit_highlimit{0.};
	double p_gausfit_lowlimit{0.};
	double p_gausfit_highlimit{0.};
	double backg_polfit_lowlimit{0.};
	double backg_polfit_highlimit{0.};

	std::cout << "Lower limit for neutron peak gaus fit: ";
	std::cin >> n_gausfit_lowlimit;
	std::cout << "Higher limit for the neutron peak gaus fit: ";
	std::cin >> n_gausfit_highlimit;
	std::cout << "Lower limit for the proton peak gaus fit: ";
	std::cin >> p_gausfit_lowlimit;
	std::cout << "Higher limit for the proton peak gaus fit: ";
	std::cin >> p_gausfit_highlimit;
	std::cout << "Lower limit for the backgroun polynomial fit: ";
	std::cin >> backg_polfit_lowlimit;
	std::cout << "Higher limit for the background polynomial fit: ";
	std::cin >> backg_polfit_highlimit;

	// Apply the fits
	TF1* n_fit = new TF1("n_fit", "gaus", n_gausfit_lowlimit, n_gausfit_highlimit);
	TF1* p_fit = new TF1("p_fit", "gaus", p_gausfit_lowlimit, p_gausfit_highlimit);
	TF1* backg_fit = new TF1("backg_fit", Form("pol%i", n_polynormial), backg_polfit_lowlimit, backg_polfit_highlimit);

	h1_dx_1->Fit(n_fit, "R");
	h1_dx_1->Fit(p_fit, "R+");
	h1_dx_1->Fit(backg_fit, "R+");
	C1->Update();

	// Get the output parameters to an array. Define a "total_fit" function "gaus(0)+gaus(3)+poln(6)".
	double par[6+n_polynormial+1];
	n_fit->GetParameters(&par[0]);
	p_fit->GetParameters(&par[3]);
	backg_fit->GetParameters(&par[6]);

	// Get lower and upper limits of which the total_fit parameter should be applied.
	double totalfit_lowlimit{0.};
	double totalfit_highlimit{0.};

	std::cout << "Lower limit for the total fit: ";
	std::cin >> totalfit_lowlimit;
	std::cout << "Higher limit for the total fit: ";
	std::cin >> totalfit_highlimit;

	TF1* total_fit = new TF1("total_fit", Form("gaus(0)+gaus(3)+pol%i(6)", n_polynormial), totalfit_lowlimit, totalfit_highlimit );
	// Set the starting values of the parameters of the "total_fit" function to the ones we get after individual sub-range fitting.
	total_fit->SetParameters(par);
	total_fit->SetLineColor(7);

	//// Apply the total fit. ////
	TString totalfit_plot = outputname + "_total_fit";
	TCanvas* C2 = new TCanvas("C2", totalfit_plot, 800, 800);
	TH1D* h1_dx_2 = new TH1D("h1_dx_2","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	T->Draw("sbs.hcal.dx>>h1_dx_2");
	h1_dx_2->Fit(total_fit, "R+");

	TLegend* total_fit_legend = new TLegend(0.45,0.7,0.9,0.9);
	//total_fit_legend->AddEntry(h1_dx_2, "Histogram of X_{expected} - X_{detected}", "l");
	total_fit_legend->AddEntry(total_fit, Form("Total fit function = gaus(0)+gaus(3)+pol%i(6)", n_polynormial), "l");
	total_fit_legend->Draw();

	C2->Update();
	////

	// Get the resulting fit parameters from the "total_fit".                                                                
	double par_from_totalfit[6+n_polynormial+1];
	total_fit->GetParameters(par_from_totalfit);

	// Define 3 TF1 functions for the neutron_function, proton_function, and backgruond_function and use the output parameter values from the "total_fit" to define them.
	TF1* n_function = new TF1("n_function", "gaus", totalfit_lowlimit, totalfit_highlimit);
	n_function->SetParameters(&par_from_totalfit[0]);
	n_function->SetLineColor(3);

	TF1* p_function = new TF1("p_function", "gaus", totalfit_lowlimit, totalfit_highlimit);
	p_function->SetParameters(&par_from_totalfit[3]);
	p_function->SetLineColor(6);

	TF1* backg_function = new TF1("backg_function", Form("pol%i", n_polynormial), totalfit_lowlimit, totalfit_highlimit);
	backg_function->SetParameters(&par_from_totalfit[6]);
	backg_function->SetLineColor(9);

	// Draw them on a canvas. The integrals of the nuetron_function and proton_function will be the respective neutron and proton yields.
	TString results_plot = outputname + "_results_plot";
	TCanvas* C3 = new TCanvas("C3", results_plot, 800, 800);
	TH1D* h1_dx_3 = new TH1D("h1_dx_3","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	T->Draw("sbs.hcal.dx>>h1_dx_3");
	n_function->Draw("same");
	p_function->Draw("same");
	backg_function->Draw("same");
	
	TLegend* results_plot_legend = new TLegend(0.45,0.7,0.9,0.9);
	results_plot_legend->AddEntry(n_function, "neutron function - gaussian", "l");
	results_plot_legend->AddEntry(p_function, "proton function - gaussian", "l");
	results_plot_legend->AddEntry(backg_function, Form("background function - %i^{th} order polynomial", n_polynormial), "l");
	results_plot_legend->Draw();

	std::cout << "\nArea under the proton function plus neutron function = " << n_function->Integral(totalfit_lowlimit, totalfit_highlimit) + p_function->Integral(totalfit_lowlimit, totalfit_highlimit) << '\n';

	C1->SaveAs(Form("%s.png", outputname.Data()));
	C2->SaveAs(Form("%s.png", totalfit_plot.Data()));
	C3->SaveAs(Form("%s.png", results_plot.Data()));

}