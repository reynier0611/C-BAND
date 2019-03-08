#ifndef __PPARTICLE_H__
#define __PPARTICLE_H__

#include <string>
#include "TVector3.h"
#include "TLorentzVector.h"

class TF1;

class PParticle
{
	public:
		PParticle();
		PParticle(TVector3 mom0,double z_m0);
		~PParticle();
		
		bool pointsToBand();

	private:
		double mp;
		double mn; 

		TVector3 mom;
        	double z_m;

};

#endif
