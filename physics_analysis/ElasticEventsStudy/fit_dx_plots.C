// Macro to fit the dx plots with two gausians for n and p peaks and a polynormial for inelastic background.
// Only to be used for LD2 data.

void fit_sub_ranges(TTree*, TCanvas*, const int, double [100]);
void do_total_fit(double&, double&, TTree*, TCanvas*, const int, double [100], double[100]);
void draw_result_functions(TTree*, TCanvas*,const int, double[100], double, double);
void compare_backgfunc_withdycut_dxdat(TTree*, TString, TCanvas*, TCanvas*, int, double, double, double[100]);

void fit_dx_plots(const char* rootfile, TString  outputname = "sbsx_sbsfield_ypercent_LD2", const int n_polynomial = 4)
{
	// Load the root file and get the TTree object "T"
	TFile* anarootfile = new TFile(rootfile, "READ");
	TTree* T = (TTree*)anarootfile->Get("T"); // Load the TTree, "T" which is made during the analysis by "gmn_ana_dxdy.C" script.

	// Define Canvas, and parameter array that goes into the sub-rang fitting of the dx plot.
	TCanvas* C1 = new TCanvas("C1", outputname, 800, 800);
	double par[6+n_polynomial+1];

	fit_sub_ranges(T, C1, n_polynomial, par);
	////  

	
	// Define Canvas, dx histo, fit function, and etc that goes into the total fit.
	double par_from_totalfit[6+n_polynomial+1];
	TString totalfit_plot = outputname + "_totalfit_plot";
	TCanvas* C2 = new TCanvas("C2", totalfit_plot, 800, 800);
	double totalfit_lowlimit{0.};
	double totalfit_highlimit{0.};

	do_total_fit(totalfit_lowlimit, totalfit_highlimit, T, C2, n_polynomial, par, par_from_totalfit);
	////


	// Draw the resulting functions from the total fit.
	TString results_plot = outputname + "_results_plot";
	TCanvas* C3 = new TCanvas("C3", results_plot, 800, 800);

	draw_result_functions(T, C3, n_polynomial, par_from_totalfit, totalfit_lowlimit, totalfit_highlimit);	
	////

	// Compare the background function with dx data with a dy cut to select obvious background events.
	
	TCanvas* C4; // = new TCanvas("C4", dy_plot, 800, 800);
	TCanvas* C5; // = new TCanvas("C5", compare_backg_plot, 800, 800);
	compare_backgfunc_withdycut_dxdat(T, outputname, C4, C5, n_polynomial, totalfit_lowlimit, totalfit_highlimit, par_from_totalfit);
	
	C1->SaveAs(Form("%s.png", outputname.Data()));
	C2->SaveAs(Form("%s.png", totalfit_plot.Data()));
	C3->SaveAs(Form("%s.png", results_plot.Data()));
}


void fit_sub_ranges(TTree* T, TCanvas* C, const int n_polynomial, double par[6+n_polynomial+1])
{
	C->cd();
	TH1D* h1_dx_1 = new TH1D("h1_dx_1","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	T->Draw("sbs.hcal.dx>>h1_dx_1");
	C->Update(); // Updates the Canvas and displays it.

	TF1* n_fit;
	TF1* p_fit;
	TF1* backg_fit;

	double n_gausfit_lowlimit{0.};
	double n_gausfit_highlimit{0.};
	double p_gausfit_lowlimit{0.};
	double p_gausfit_highlimit{0.};
	double backg_polfit_lowlimit{0.};
	double backg_polfit_highlimit{0.};

	// Select the sub-ranges for the neutron and proton gaus fits and the background fit (n_polynormal'th order polynomial).
	string fits_good = "n";
	bool b_fits_good{false};

	do 
	{
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

		n_fit = new TF1("n_fit", "gaus", n_gausfit_lowlimit, n_gausfit_highlimit);
		p_fit = new TF1("p_fit", "gaus", p_gausfit_lowlimit, p_gausfit_highlimit);
		backg_fit = new TF1("backg_fit", Form("pol%i", n_polynomial), backg_polfit_lowlimit, backg_polfit_highlimit);

		h1_dx_1->Fit(n_fit, "R");
		h1_dx_1->Fit(p_fit, "R+");
		h1_dx_1->Fit(backg_fit, "R+");
		C->Update();

		std::cout << "\nIs the sub-range fit good? y/n: ";
		std::cin >> fits_good;
		
		if (fits_good == "n" || fits_good == "no") 
		{
			std::cout << "Let us try again. \n ";
			std::cout << '\n';
			b_fits_good = false;
		}
		else if (fits_good == "y" || fits_good == "yes") 
		{
			std::cout << '\n';
			b_fits_good = true;
		}
	}
	while(!b_fits_good);

	// Get the output parameters to an array. Define a "total_fit" function "gaus(0)+gaus(3)+poln(6)".
	n_fit->GetParameters(&par[0]);
	p_fit->GetParameters(&par[3]);
	backg_fit->GetParameters(&par[6]);
}

void do_total_fit(double& totalfit_lowlimit, double& totalfit_highlimit, TTree* T, TCanvas* C, const int n_polynomial, double par[6+n_polynomial+1], double par_from_totalfit[6+n_polynomial+1])
{
	// Get lower and upper limits of which the total_fit parameter should be applied.
	std::cout << "Lower limit for the total fit: ";
	std::cin >> totalfit_lowlimit;
	std::cout << "Higher limit for the total fit: ";
	std::cin >> totalfit_highlimit;

	TF1* total_fit = new TF1("total_fit", Form("gaus(0)+gaus(3)+pol%i(6)", n_polynomial), totalfit_lowlimit, totalfit_highlimit );
	// Set the starting values of the parameters of the "total_fit" function to the ones we get after individual sub-range fitting.
	total_fit->SetParameters(par);
	total_fit->SetLineColor(7);

	C->cd();
	TH1D* h1_dx_2 = new TH1D("h1_dx_2","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	T->Draw("sbs.hcal.dx>>h1_dx_2");
	h1_dx_2->Fit(total_fit, "R+"); //// Apply the total fit. ////

	TLegend* total_fit_legend = new TLegend(0.45,0.7,0.9,0.9);
	total_fit_legend->AddEntry(total_fit, Form("Total fit function = gaus(0)+gaus(3)+pol%i(6)", n_polynomial), "l");
	total_fit_legend->Draw();

	C->Update();

	total_fit->GetParameters(par_from_totalfit); // Get the resulting fit parameters from the "total_fit".
}

void draw_result_functions(TTree* T, TCanvas* C, const int n_polynomial, double par_from_totalfit[6+n_polynomial+1], double totalfit_lowlimit, double totalfit_highlimit)
{
	// Define 3 TF1 functions for the neutron_function, proton_function, and backgruond_function and use the output parameter values from the "total_fit" to define them.
	TF1* n_function = new TF1("n_function", "gaus", totalfit_lowlimit, totalfit_highlimit);
	n_function->SetParameters(&par_from_totalfit[0]);
	n_function->SetLineColor(3);

	TF1* p_function = new TF1("p_function", "gaus", totalfit_lowlimit, totalfit_highlimit);
	p_function->SetParameters(&par_from_totalfit[3]);
	p_function->SetLineColor(6);

	TF1* backg_function = new TF1("backg_function", Form("pol%i", n_polynomial), totalfit_lowlimit, totalfit_highlimit);
	backg_function->SetParameters(&par_from_totalfit[6]);
	backg_function->SetLineColor(2);

	TH1D* h1_dx_3 = new TH1D("h1_dx_3","dx - with cuts;X_{HCal}-X_{expected} (m)",1000,-4,6);
	C->cd();
	T->Draw("sbs.hcal.dx>>h1_dx_3");

	n_function->Draw("same");
	p_function->Draw("same");
	backg_function->Draw("same");
	
	TLegend* results_plot_legend = new TLegend(0.45,0.7,0.9,0.9);
	results_plot_legend->AddEntry(n_function, "neutron function - gaussian", "l");
	results_plot_legend->AddEntry(p_function, "proton function - gaussian", "l");
	results_plot_legend->AddEntry(backg_function, Form("background function - %i^{th} order polynomial", n_polynomial), "l");
	results_plot_legend->Draw();

	C->Update();

	std::cout << "\nArea under the proton function plus neutron function = " << n_function->Integral(totalfit_lowlimit, totalfit_highlimit) + p_function->Integral(totalfit_lowlimit, totalfit_highlimit) << '\n';
}

void compare_backgfunc_withdycut_dxdat(TTree* T, TString outputname, TCanvas* C_1, TCanvas* C_2, int n_polynomial, double totalfit_lowlimit, double totalfit_highlimit, double par_from_totalfit[6+n_polynomial+1])
{
	TH1D* h1_dy = new TH1D("h1_dy", "dy distribution;Y_{HCal}-Y_{expected}", 400, -2, 2);
	TString dy_plot = outputname + "_dy_plot";
	C_1 = new TCanvas("C4", dy_plot, 800, 800);
	//C_1->cd();
	T->Draw("sbs.hcal.dy>>h1_dy");
	C_1->Update();

	// Manually define a dy cut region by looking at the above h1_dy plot.
	double dy_cut_lowlimit{0.};
	double dy_cut_highlimit{0.};
	TLine* ll;
	TLine* rl;

	TString compare_backg_plot = outputname + "_compare_backg_plot";
	C_2 = new TCanvas("C5", compare_backg_plot, 800, 800);
	TH1D* h1_dx_outside_QE_dy_cut = new TH1D("h1_dx_outside_QE_dy_cut", "dx distribution of data outside QE dy cut;X_{HCal}-X_{expected} (m)",1000,-4,6);

	string dy_cut_good = "n";
	bool b_dy_cut_good{false};

	do
	{
		std::cout << "\nLower limit for the dy cut: ";
		std::cin >> dy_cut_lowlimit;
		std::cout << "Higher limit for the dy cut: "; 
		std::cin >> dy_cut_highlimit;

		double dy_hist_max{h1_dy->GetMaximum()};
		ll = new TLine(dy_cut_lowlimit, 0, dy_cut_lowlimit, dy_hist_max);
		rl = new TLine(dy_cut_highlimit, 0, dy_cut_highlimit, dy_hist_max);
		ll->SetLineColorAlpha(kBlue,0);
		ll->SetLineWidth(2);
		ll->SetLineStyle(9);
		rl->SetLineColorAlpha(kBlue,0);
		rl->SetLineWidth(2);
		rl->SetLineStyle(9);
		C_1->cd();
		ll->Draw();
		rl->Draw();
		C_1->Update();
		
		/*TString compare_backg_plot = outputname + "_compare_backg_plot";
		C_2 = new TCanvas("C5", compare_backg_plot, 800, 800);

		TH1D* h1_dx_outside_QE_dy_cut = new TH1D("h1_dx_outside_QE_dy_cut", "dx distribution of data outside QE dy cut;X_{HCal}-X_{expected} (m)",1000,-4,6);*/
		C_2->cd();
		T->Draw("sbs.hcal.dx>>h1_dx_outside_QE_dy_cut", Form("sbs.hcal.dy<%f||sbs.hcal.dy>%f", dy_cut_lowlimit, dy_cut_highlimit));
		C_2->Update();

		std::cout << "\nThe dx distribution looks like mostly inelastic background? y/n: ";
		std::cin >> dy_cut_good;

		if (dy_cut_good == "n" || dy_cut_good == "no") 
		{
			std::cout << "Let us try again. \n ";
			std::cout << '\n';
			b_dy_cut_good = false;
		}
		else if (dy_cut_good == "y" || dy_cut_good == "yes") 
		{
			std::cout << '\n';
			b_dy_cut_good = true;
		}
	
	}
	while(!b_dy_cut_good);

	//
	TH1* h1_dx_outside_QE_dy_cut_clone = (TH1*)h1_dx_outside_QE_dy_cut->Clone("h1_dx_outside_QE_dy_cut_clone");
	h1_dx_outside_QE_dy_cut_clone->Scale(1./h1_dx_outside_QE_dy_cut_clone->Integral(), "width");
	//h1_dx_outside_QE_dy_cut_clone->Draw("HIST");

	TF1* backg_function = new TF1("backg_function", Form("pol%i", n_polynomial), totalfit_lowlimit, totalfit_highlimit);
	backg_function->SetParameters(&par_from_totalfit[6]);
	TH1* backg_function_histo = (TH1*)backg_function->GetHistogram();
	backg_function_histo->SetLineColor(2);
	backg_function_histo->Scale(1./backg_function_histo->Integral(), "width");
	//backg_function_histo->Draw("HIST+SAME");

	double dx_dis_max{h1_dx_outside_QE_dy_cut_clone->GetMaximum()};
	double backg_func_max{backg_function_histo->GetMaximum()};

	if(backg_func_max > dx_dis_max) 
	{
		backg_function_histo->Draw("HIST");
		h1_dx_outside_QE_dy_cut_clone->Draw("HIST+SAME");
	}
	else 
	{
		h1_dx_outside_QE_dy_cut_clone->Draw("HIST");
		backg_function_histo->Draw("HIST+SAME");
	}

	TLegend* comparebackgfunc_plot_legend = new TLegend(0.66,0.7,0.9,0.9);
	comparebackgfunc_plot_legend->AddEntry(h1_dx_outside_QE_dy_cut_clone, "dx with dy cut", "l");
	comparebackgfunc_plot_legend->AddEntry(backg_function_histo, Form("%i th order background pol", n_polynomial), "l");
	comparebackgfunc_plot_legend->Draw();

	C_1->SaveAs(Form("%s.png", dy_plot.Data()));
	C_2->SaveAs(Form("%s.png", compare_backg_plot.Data()));	
}	
