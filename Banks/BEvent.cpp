#include "BEvent.h"
#include "TVector3.h"
// ==============================================================
void   BEvent::init(const char *bankName, hipo::reader &r){
	initBranches(bankName,r);
	NRUN_order    = getEntryOrder("NRUN"   );
	NEVENT_order  = getEntryOrder("NEVENT" );
	EVNTime_order = getEntryOrder("EVNTime");
	TYPE_order    = getEntryOrder("TYPE"   );
	EvCAT_order   = getEntryOrder("EvCAT"  );
	NPGP_order    = getEntryOrder("NPGP"   );
	TRG_order     = getEntryOrder("TRG"    );
	BCG_order     = getEntryOrder("BCG"    );
	LT_order      = getEntryOrder("LT"     );
	STTime_order  = getEntryOrder("STTime" );
	RFTime_order  = getEntryOrder("RFTime" );
	Helic_order   = getEntryOrder("Helic"  );
	PTIME_order   = getEntryOrder("PTIME"  );
}
// ==============================================================
BEvent::~BEvent(){}
// ==============================================================
