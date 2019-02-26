#include "particle.h"
#include "TVector3.h"
// ==============================================================
void   particle::init(const char *bankName, hipo::reader &r){
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
particle::~particle(){}
// ==============================================================
TVector3 particle::getV3v(int index)
{
	TVector3 * V3v = new TVector3();
        float vx = particle::getVx(index);
        float vy = particle::getVy(index);
        float vz = particle::getVz(index);
        V3v -> SetXYZ(vx,vy,vz);
        return * V3v;
}
// ==============================================================
TVector3 particle::getV3P(int index)
{
	TVector3 * V3P = new TVector3();
	float px = particle::getPx(index);
	float py = particle::getPy(index);
	float pz = particle::getPz(index);
	V3P -> SetXYZ(px,py,pz);
	return * V3P;
}
// ==============================================================
