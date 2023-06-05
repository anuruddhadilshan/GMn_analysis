#ifndef FIDUCIALCUT_H
#define FIDUCIALCUT_H

#include "HCalConstants.h"

// Find Average Proton Deflection from the SBS magnet for a given kinematic setting and SBS magnet field scale.
// Found by seeing how much the LH2 data gets shifted by the SBS magnet for the given kinematic setting and SBS magnet field scale.
std::vector<std::vector<double>> avg_proton_deflection_onHCal = { //{kinematicsetting(sbs #), sbs fieldscale, avg proton deflection from the analysis of the LH2 data}
	{4, 0, 0},
	{4, 30.0, -0.66863820},
	{4, 50.0, -1.1193552},
	{7, 85.0, -0.693155},
	{8, 0, 0},
	{8, 50.0, -0.661306},
	{8, 70.0, -0.920279},
	{8, 100.0, -1.29957},
	{11, 0, 0},
	{11, 100.0, -0.7144345},
	{14, 0, 0},
	{14, 70.0, -0.8045086}
};

double getavg_proton_deflection(int kine_num, double sbsfieldscale)
{
	for (int i=0; i<avg_proton_deflection_onHCal.size(); ++i)
	{
		if( avg_proton_deflection_onHCal[i][0] == kine_num && avg_proton_deflection_onHCal[i][1] == sbsfieldscale )
		{
			double average_proton_deflection = avg_proton_deflection_onHCal[i][2];
			std::cout << "\nAverage proton deflection by the SBS magnet, calculated using LH2 data from SBS" << kine_num << " kinematic setting and the SBS field scale " << sbsfieldscale << "% = " << average_proton_deflection << " m\n";
			std::cout << "### IMPORTANT: The above deflection will be used to estimate the proton deflection for the fiducial cut ###\n";
			return average_proton_deflection;
		}
	}

	std::cout << "\n### ERROR: Either the kinematic setting number or the percentage SBS magnet field scale entered are invalid ####\n";
	return 0;
}
////

//// Fiducial Cut ////

// Define the boundaries of HCal dimensions with respect to the HCal coordinate system.
const double hcal_active_xlow = hcal_topXpos;
const double hcal_active_xhigh = hcal_botXpos;
const double hcal_active_ylow = hcal_rightYpos;
const double hcal_active_yhigh = hcal_leftYpos;

// Define the boundaries of HCal active area with some "safety margins" included.
// Safety margin = Exclude the two outer columns and two outer blocks.
const double hcal_active_xlow_safe =  hcal_active_xlow + hcalblk_h;
const double hcal_active_xhigh_safe = hcal_active_xhigh - hcalblk_h;
const double hcal_active_ylow_safe = hcal_active_ylow + hcalblk_w;
const double hcal_active_yhigh_safe = hcal_active_yhigh - hcalblk_w;

bool fiducial_cut(double avg_proton_deflection, double xexpected_hcal, double yexpected_hcal)
{
	// neutron coordinates: The kinematic calculation to find "xexpected_hcal" and "yexpected_hcal" does not take into account the SBS magnet deflection.
	// So essentially what we get is the neutron position from that calculation.
	double neutron_xpos_hcal{xexpected_hcal};
	double neutron_ypos_hcal{yexpected_hcal};
	// First we see whether neutron falls within the HCal. If not we return a "false" value as we do not need to proceed if the neutron already is outside the region we consider.
	if( neutron_xpos_hcal>=hcal_active_xhigh_safe || neutron_xpos_hcal<=hcal_active_xlow_safe || neutron_ypos_hcal>=hcal_active_yhigh_safe || neutron_ypos_hcal<= hcal_active_ylow_safe ) return false;

	// proton coordinates: Just add the expected average deflection expected to the neutron's x pos value.
	double proton_xpos_hcal{xexpected_hcal+avg_proton_deflection};
	double proton_ypos_hcal{yexpected_hcal};  
	// We are here only if we have established already above that the neutron is within the HCal region we consider.
	// Now see whether the proton will fall within the considered region in HCal.
	if( proton_xpos_hcal>=hcal_active_xhigh_safe || proton_xpos_hcal<=hcal_active_xlow_safe || proton_ypos_hcal>=hcal_active_yhigh_safe || proton_ypos_hcal<= hcal_active_ylow_safe ) return false;

	// We are here only if both the neutron and proton are within the HCal region considered. So we return a "true" value only at this stage.
	return true;
}

////

void draw_fiducialcut(TH2D* h2_dxdy,const char* h2_dxdy_filename)
{
	TCanvas* C = new TCanvas();

	TLine* hcal_active_horizontal_low = new TLine(hcal_active_ylow,hcal_active_xlow,hcal_active_yhigh,hcal_active_xlow);
	TLine* hcal_active_horizontal_high = new TLine(hcal_active_ylow,hcal_active_xhigh,hcal_active_yhigh,hcal_active_xhigh);
	TLine* hcal_active_vertical_left = new TLine(hcal_active_yhigh,hcal_active_xlow,hcal_active_yhigh,hcal_active_xhigh);
	TLine* hcal_active_vertical_right = new TLine(hcal_active_ylow,hcal_active_xlow,hcal_active_ylow,hcal_active_xhigh);
	TLine* hcal_active_horizontalsafemarg_low = new TLine(hcal_active_ylow_safe,hcal_active_xlow_safe,hcal_active_yhigh_safe,hcal_active_xlow_safe);
	TLine* hcal_active_horizontalsafemarg_high = new TLine(hcal_active_ylow_safe,hcal_active_xhigh_safe,hcal_active_yhigh_safe,hcal_active_xhigh_safe);
	TLine* hcal_active_verticalsafemarg_left = new TLine(hcal_active_yhigh_safe,hcal_active_xlow_safe,hcal_active_yhigh_safe,hcal_active_xhigh_safe);
	TLine* hcal_active_verticalsafemarg_right = new TLine(hcal_active_ylow_safe,hcal_active_xlow_safe,hcal_active_ylow_safe,hcal_active_xhigh_safe);

	hcal_active_horizontal_low->SetLineColorAlpha(kGreen,0);
	hcal_active_horizontal_low->SetLineWidth(2);
	hcal_active_horizontal_high->SetLineColorAlpha(kGreen,0);
	hcal_active_horizontal_high->SetLineWidth(2);
	hcal_active_vertical_left->SetLineColorAlpha(kGreen,0);
	hcal_active_vertical_left->SetLineWidth(2);
	hcal_active_vertical_right->SetLineColorAlpha(kGreen,0);
	hcal_active_vertical_right->SetLineWidth(2);

	hcal_active_horizontalsafemarg_low->SetLineColorAlpha(kRed,0);
	hcal_active_horizontalsafemarg_low->SetLineWidth(2);
	hcal_active_horizontalsafemarg_low->SetLineStyle(9);
	hcal_active_horizontalsafemarg_high->SetLineColorAlpha(kRed,0);
	hcal_active_horizontalsafemarg_high->SetLineWidth(2);
	hcal_active_horizontalsafemarg_high->SetLineStyle(9);
	hcal_active_verticalsafemarg_left->SetLineColorAlpha(kRed,0);
	hcal_active_verticalsafemarg_left->SetLineWidth(2);
	hcal_active_verticalsafemarg_left->SetLineStyle(9);
	hcal_active_verticalsafemarg_right->SetLineColorAlpha(kRed,0);
	hcal_active_verticalsafemarg_right->SetLineWidth(2);
	hcal_active_verticalsafemarg_right->SetLineStyle(9);

	h2_dxdy->Draw("COLZ");
	hcal_active_horizontal_low->Draw();
	hcal_active_horizontal_high->Draw();
	hcal_active_vertical_left->Draw();
	hcal_active_vertical_right->Draw();
	hcal_active_horizontalsafemarg_low->Draw();
	hcal_active_horizontalsafemarg_high->Draw();
	hcal_active_verticalsafemarg_right->Draw();
	hcal_active_verticalsafemarg_left->Draw();

	C->SaveAs(h2_dxdy_filename);
}

#endif