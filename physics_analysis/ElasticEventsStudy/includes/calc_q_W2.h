#ifndef CALC_Q_W2_H
#define CALC_Q_W2_H

//Defining four vectors for the scattered particles
TLorentzVector kprime; //Four vector of the scattered electron.
TLorentzVector q;      //Four momentum trasnferred to the scattered nucleon.

void calcq(TLorentzVector Pbeam, double epx[10], double epy[10], double epz[10], double ep[10])
{
	kprime.SetPxPyPzE(epx[0],epy[0],epz[0],ep[0]); 
	q = Pbeam - kprime; 
}

double calcW2(TLorentzVector Ptarg)
{
	return (Ptarg+q).M2(); //This is the invariant mass squared (W^2) of the scattering reaction.
}

#endif