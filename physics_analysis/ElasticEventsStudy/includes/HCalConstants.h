#ifndef HCALCONSTANTS_H
#define HCALCONSTANTS_H

const double hcalblk_h = 0.15494; // Height of all HCAL blocks in m from MC database
const double hcalblk_w = 0.15875; // Width of all HCAL blocks in m from MC database
const int Nhcal_rows = 24;
const int Nhcal_columns = 12;
const double hcal_width = hcalblk_w*Nhcal_columns; //1.905; // m
const double hcal_height = hcalblk_h*Nhcal_rows; // 3.71856 m
const double hcal_topXpos = -2.355005; // Distance to the top of the HCal w.r.t HCal origin.
const double hcal_botXpos = 1.454995; // Distance to the bottom of the HCal w.r.t HCal origin.
const double hcal_leftYpos = 0.92964; // Distance to the left side of HCal (when looking from up-stream to down-stream direction) w.r.t HCal origin.
const double hcal_rightYpos = -0.92964; // Distance to the right side of HCal (when looking from up-stream to down-stream direction) w.r.t HCal origin.
const double hcal_height_abovebeamline = -0.2897; // Vertical distance (X) of the HCal origin above beamline.
//const double hcal_height_abovebeamline = -0.75; // Vertical distance (X) of the HCal origin above beamline.


#endif