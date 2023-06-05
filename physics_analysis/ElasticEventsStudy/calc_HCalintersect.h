#ifndef CALC_HCALINTERSECT_H
#define CALC_HCALINTERSECT_H

#include "HCalConstants.h"

//Defining HCal vectors
TVector3 hcal_origin; // Vector pointing to the origin of HCal starting from the origin of the Hall coordinate system.
TVector3 hcalXaxis;				
TVector3 hcalYaxis;
TVector3 hcalZaxis;

void make_HCal_vectors(double hcaldist, double hcaltheta)
{
	//In the Tree, HCal cluster position is given in the HCal coordinate system. Where x - vertically down, y - beam left, z - beam downstream.
	//In this analysis we will recontsruct the position vector of the hadron striking the HCal with the help of the recontructed kinematics of the scattered electron.
	//This vecotr will be define in the Hall coordinate system as the 4 momentum of the scattered electron in the Tree is given in the Hall coordinate system.
	//Thereofre, we need to take the dot product between the recontructed postion vector of the Hadron with the unit vectors of the HCal x axis and y axis, 
	//to get the predicted x and y cluster posion respectively, on the face of the HCal in HCal coordinates.
	//Defining unit vectors for HCal coordinate system (Using the Hall coordinate system).
	hcalXaxis.SetXYZ(0,-1,0);
	hcalZaxis.SetXYZ(-sin(hcaltheta),0,cos(hcaltheta));
	hcalYaxis = hcalZaxis.Cross(hcalXaxis).Unit();

	//3D vector defining the position of the HCal origin w.r.t the Hall origin. Hall coordinate system = x-beam left, y-verically up, z-beam downstream.
	hcal_origin = hcaldist*hcalZaxis + hcal_height_abovebeamline*hcalXaxis;
}


// Make the HCal vectors for MC simulation data analysis with the HCal vertical offset excluded.
void make_HCal_vectors_forsim(double hcaldist, double hcaltheta)
{
	//In the Tree, HCal cluster position is given in the HCal coordinate system. Where x - vertically down, y - beam left, z - beam downstream.
	//In this analysis we will recontsruct the position vector of the hadron striking the HCal with the help of the recontructed kinematics of the scattered electron.
	//This vecotr will be define in the Hall coordinate system as the 4 momentum of the scattered electron in the Tree is given in the Hall coordinate system.
	//Thereofre, we need to take the dot product between the recontructed postion vector of the Hadron with the unit vectors of the HCal x axis and y axis, 
	//to get the predicted x and y cluster posion respectively, on the face of the HCal in HCal coordinates.
	//Defining unit vectors for HCal coordinate system (Using the Hall coordinate system).
	hcalXaxis.SetXYZ(0,-1,0);
	hcalZaxis.SetXYZ(-sin(hcaltheta),0,cos(hcaltheta));
	hcalYaxis = hcalZaxis.Cross(hcalXaxis).Unit();

	//3D vector defining the position of the HCal origin w.r.t the Hall origin. Hall coordinate system = x-beam left, y-verically up, z-beam downstream.
	hcal_origin = hcaldist*hcalZaxis;
}



// Vector calculation to find the vector that starts from the orgin of the Hall cordinate system and points to spot on the face of the HCal that the nucleon strikes.
/*void calc_expected_xyonHCal(double vz[10], double &xexpected_hcal, double &yexpected_hcal)
{
	TVector3 qdirection = q.Vect().Unit(); //Unit vector defining the direction of the spacial component of the q.
	TVector3 vertex(0,0,vz[0]); // Vertex vector (the vector going from the origin of the Hall coordinate system(~middle of the target) to the (e,e'p/n) interaction point)
	TVector3 vertextoHCalorigin = hcal_position - vertex;
	double hcalZdistancefromvertex = vertextoHCalorigin.Dot(hcalZaxis);
	double hadronvectormagnitude = hcalZdistancefromvertex/qdirection.Dot(hcalZaxis);
	TVector3 hcalintersect = hadronvectormagnitude*qdirection + vertex;

	xexpected_hcal = hcalintersect.Dot(hcalXaxis);
	yexpected_hcal = hcalintersect.Dot(hcalYaxis);
}*/

// Vector calculation to find the vector that starts from the orgin of the Hall cordinate system and points to spot on the face of the HCal that the nucleon strikes.
void calc_expected_xyonHCal(TLorentzVector q, double vz[10], double &xexpected_hcal, double &yexpected_hcal)
{
	TVector3 vertex(0,0,vz[0]); // Vertex vector (the vector going from the origin of the Hall coordinate system(~middle of the target) to the (e,e'p/n) interaction point)
	TVector3 vertextoHCalorigin = hcal_origin - vertex;
	double hcalZdistancefromvertex = vertextoHCalorigin.Dot(hcalZaxis);
	TVector3 qdirection = q.Vect().Unit(); //Unit vector defining the direction of the spacial component of the q.
	double hadronvectormagnitude = hcalZdistancefromvertex/qdirection.Dot(hcalZaxis);
	TVector3 hcalintersect = hadronvectormagnitude*qdirection + vertex;

	xexpected_hcal = (hcalintersect-hcal_origin).Dot(hcalXaxis);
	yexpected_hcal = (hcalintersect-hcal_origin).Dot(hcalYaxis);
}

#endif