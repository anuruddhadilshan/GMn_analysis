#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "HCalConstants.h"

#ifndef CALC_HCALINTERSECT_H
#define CALC_HCALINTERSECT_H

//Defining HCal vectors
ROOT::Math::XYZVector hcal_origin; // Vector pointing to the origin of HCal starting from the origin of the Hall coordinate system.
ROOT::Math::XYZVector hcalXaxis;				
ROOT::Math::XYZVector hcalYaxis;
ROOT::Math::XYZVector hcalZaxis;

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
	ROOT::Math::XYZVector qdirection = q.Vect().Unit(); //Unit vector defining the direction of the spacial component of the q.
	ROOT::Math::XYZVector vertex(0,0,vz[0]); // Vertex vector (the vector going from the origin of the Hall coordinate system(~middle of the target) to the (e,e'p/n) interaction point)
	ROOT::Math::XYZVector vertextoHCalorigin = hcal_position - vertex;
	double hcalZdistancefromvertex = vertextoHCalorigin.Dot(hcalZaxis);
	double hadronvectormagnitude = hcalZdistancefromvertex/qdirection.Dot(hcalZaxis);
	ROOT::Math::XYZVector hcalintersect = hadronvectormagnitude*qdirection + vertex;

	xexpected_hcal = hcalintersect.Dot(hcalXaxis);
	yexpected_hcal = hcalintersect.Dot(hcalYaxis);
}*/

// Vector calculation to find the vector that starts from the orgin of the Hall cordinate system and points to spot on the face of the HCal that the nucleon strikes.
void calc_expected_xyonHCal(ROOT::Math::PxPyPzEVector& q, double vz[10], double& xexpected_hcal, double& yexpected_hcal)
{
	ROOT::Math::XYZVector vertex(0,0,vz[0]); // Vertex vector (the vector going from the origin of the Hall coordinate system(~middle of the target) to the (e,e'p/n) interaction point)
	ROOT::Math::XYZVector vertextoHCalorigin = hcal_origin - vertex;
	double hcalZdistancefromvertex = vertextoHCalorigin.Dot(hcalZaxis);
	ROOT::Math::XYZVector qdirection = q.Vect().Unit(); //Unit vector defining the direction of the spacial component of the q.
	double hadronvectormagnitude = hcalZdistancefromvertex/qdirection.Dot(hcalZaxis);
	ROOT::Math::XYZVector hcalintersect = hadronvectormagnitude*qdirection + vertex;

	xexpected_hcal = (hcalintersect-hcal_origin).Dot(hcalXaxis);
	yexpected_hcal = (hcalintersect-hcal_origin).Dot(hcalYaxis);
}


// A single class that does all the HCal vector calculations.

class HCalVectors
{
	ROOT::Math::XYZVector m_hcal_origin; // Vector pointing to the origin of HCal starting from the origin of the Hall coordinate system.
	ROOT::Math::XYZVector m_hcalXaxis;				
	ROOT::Math::XYZVector m_hcalYaxis;
	ROOT::Math::XYZVector m_hcalZaxis;
	
public:
	void make_HCal_vectors(double hcaldist, double hcaltheta)
	{
		m_hcalXaxis.SetXYZ(0,-1,0);
		m_hcalZaxis.SetXYZ(-sin(hcaltheta),0,cos(hcaltheta));
		m_hcalYaxis = m_hcalZaxis.Cross(m_hcalXaxis).Unit();

		//3D vector defining the position of the HCal origin w.r.t the Hall origin. Hall coordinate system = x-beam left, y-verically up, z-beam downstream.
		m_hcal_origin = hcaldist*m_hcalZaxis + hcal_height_abovebeamline*m_hcalXaxis;
	}

	void make_HCal_vectors_forsim(double hcaldist, double hcaltheta)
	{
		m_hcalXaxis.SetXYZ(0,-1,0);
		m_hcalZaxis.SetXYZ(-sin(hcaltheta),0,cos(hcaltheta));
		m_hcalYaxis = m_hcalZaxis.Cross(m_hcalXaxis).Unit();

		//3D vector defining the position of the HCal origin w.r.t the Hall origin. Hall coordinate system = x-beam left, y-verically up, z-beam downstream.
		m_hcal_origin = hcaldist*m_hcalZaxis;
	}

private: 
	ROOT::Math::XYZVector m_vertex;
	ROOT::Math::XYZVector m_vertextoHCalorigin;
	double m_hcalZdistancefromvertex {0.};
	ROOT::Math::XYZVector m_qdirection;
	double m_hadronvectormagnitude {0.};
	ROOT::Math::XYZVector m_hcalintersect;
	double m_xexpected_hcal {0.};
	double m_yexpected_hcal {0.};

public:
	void calc_expected_xyonHCal(ROOT::Math::PxPyPzEVector& q, double vz[10])
	{
		m_vertex.SetXYZ(0,0,vz[0]); // Vertex vector (the vector going from the origin of the Hall coordinate system(~middle of the target) to the (e,e'p/n) interaction point)
		m_vertextoHCalorigin = m_hcal_origin - m_vertex;
		m_hcalZdistancefromvertex = m_vertextoHCalorigin.Dot(m_hcalZaxis);
		m_qdirection = q.Vect().Unit(); //Unit vector defining the direction of the spacial component of the q.
		m_hadronvectormagnitude = m_hcalZdistancefromvertex/m_qdirection.Dot(m_hcalZaxis);
		m_hcalintersect = m_hadronvectormagnitude*m_qdirection + m_vertex;

		m_xexpected_hcal = (m_hcalintersect-m_hcal_origin).Dot(m_hcalXaxis);
		m_yexpected_hcal = (m_hcalintersect-m_hcal_origin).Dot(m_hcalYaxis);
	}

	double return_xexpected()
	{
		return m_xexpected_hcal;
	}

	double return_yexpected()
	{
		return m_yexpected_hcal;
	}

	ROOT::Math::XYZVector return_q3directionvect()
	{
		return m_qdirection;
	}

	ROOT::Math::XYZVector return_HCalZaxis()
	{
		return m_hcalZaxis;
	}

	ROOT::Math::XYZVector return_VertextoHCalOrigin()
	{
		return m_vertextoHCalorigin;
	}

	double return_hcalZdistancefromvertex()
	{
		return m_hcalZdistancefromvertex;
	}

};

#endif