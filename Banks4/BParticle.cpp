#include "BParticle.h"
#include "TVector3.h"
// ==============================================================
void   BParticle::init(const char *bankName, hipo::reader &r){
	initBranches(bankName,r);
	pid_order     = getEntryOrder("pid"    );
	px_order      = getEntryOrder("px"     );
	py_order      = getEntryOrder("py"     );
	pz_order      = getEntryOrder("pz"     );
	vx_order      = getEntryOrder("vx"     );
	vz_order      = getEntryOrder("vy"     );
	vz_order      = getEntryOrder("vz"     );
	charge_order  = getEntryOrder("charge" );
	beta_order    = getEntryOrder("beta"   );
	chi2pid_order = getEntryOrder("chi2pid");
	status_order  = getEntryOrder("status" );
}
// ==============================================================
BParticle::~BParticle(){}
// ==============================================================
TVector3 BParticle::getV3v(int index)
{
	TVector3 * V3v = new TVector3();
        float vx = BParticle::getVx(index);
        float vy = BParticle::getVy(index);
        float vz = BParticle::getVz(index);
        V3v -> SetXYZ(vx,vy,vz);
        return * V3v;
}
// ==============================================================
TVector3 BParticle::getV3P(int index)
{
	TVector3 * V3P = new TVector3();
	float px = BParticle::getPx(index);
	float py = BParticle::getPy(index);
	float pz = BParticle::getPz(index);
	V3P -> SetXYZ(px,py,pz);
	return * V3P;
}
// ==============================================================
