#include "BEvent.h"
#include "TVector3.h"
// ==============================================================
void   BEvent::init(const char *bankName, hipo::reader &r){
	initBranches(bankName,r);
	Category_order  = getEntryOrder("category"  );
	Topo_order      = getEntryOrder("topology"   );
	BCG_order       = getEntryOrder("beamCharge");
	LiveTime_order  = getEntryOrder("liveTime"  );
	StartTime_order = getEntryOrder("startTime" );
	RFTime_order    = getEntryOrder("RFTime" );
	Helic_order     = getEntryOrder("helicity"  );
	HelicRaw_order  = getEntryOrder("helicityRaw"  );
	ProcTime_order  = getEntryOrder("procTime"  );
}
// ==============================================================
BEvent::~BEvent(){}
// ==============================================================
