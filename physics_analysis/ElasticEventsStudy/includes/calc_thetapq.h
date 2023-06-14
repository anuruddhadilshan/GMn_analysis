#include "HCalConstants.h"
#include "calc_HCalintersect.h"
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include "TLorentzVector.h"
#include "TVector3.h"


#ifndef CALC_THETAPQ_H
#define CALC_THETAPQ_H

class Calc_thetapq
{
	ROOT::Math::XYZVector m_P3vect;
	double m_hcalZdistancefromvertex {0.}; 
	ROOT::Math::XYZVector m_qdisvect;
	double m_cos_theta_P {0.};
	double m_theta_Pq_degrees {0.};

public:
	void make_P3vect(HCalVectors& hcal_vect, double xhcal, double yhcal)
	{
		ROOT::Math::XYZVector hcalorigintocluspos(yhcal, -xhcal, 0); // Transformed to HCal coordinate sys -> Hall coordinate sys.
		m_P3vect = hcal_vect.return_VertextoHCalOrigin() + hcalorigintocluspos;
	}

	void make_qdisvect(HCalVectors& hcal_vect)
	{
		m_hcalZdistancefromvertex = hcal_vect.return_hcalZdistancefromvertex();
		double qdisvect_mag = m_hcalZdistancefromvertex / (hcal_vect.return_q3directionvect().Dot(hcal_vect.return_HCalZaxis()));
		m_qdisvect = qdisvect_mag*hcal_vect.return_q3directionvect();
	}

	double return_thetaPq(HCalVectors& hcal_vect, double xhcal, double yhcal)
	{
		make_P3vect(hcal_vect, xhcal, yhcal);
		make_qdisvect(hcal_vect);
		double cos_thetapq = (m_P3vect.Dot(m_qdisvect))/(m_P3vect.r()*m_qdisvect.r());
		m_theta_Pq_degrees = acos(cos_thetapq)*TMath::RadToDeg();
		return m_theta_Pq_degrees;
	}
};

#endif