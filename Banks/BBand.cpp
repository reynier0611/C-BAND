#include "BBand.h"
#include "TVector3.h"
// ==============================================================
void   BBand::init(const char *bankName, hipo::reader &r){
	initBranches(bankName,r);
	id_order           = getEntryOrder("id"          );
        sector_order       = getEntryOrder("sector"      );
        layer_order        = getEntryOrder("layer"       );
        component_order    = getEntryOrder("component"   );
        meantimeTdc_order  = getEntryOrder("meantimeTdc" );
        meantimeFadc_order = getEntryOrder("meantimeFadc");
        difftimeTdc_order  = getEntryOrder("difftimeTdc" );
        difftimeFadc_order = getEntryOrder("difftimeFadc");
        adcLcorr_order     = getEntryOrder("adcLcorr"    );
        adcRcorr_order     = getEntryOrder(" adcRcorr"   );
        tFadcLcorr_order   = getEntryOrder("tFadcLcorr"  );
        tFadcRcorr_order   = getEntryOrder("tFadcRcorr"  );
        tTdcLcorr_order    = getEntryOrder("tTdcLcorr"   );
        tTdcRcorr_order    = getEntryOrder("tTdcRcorr"   );
        x_order            = getEntryOrder("x"           );
        y_order            = getEntryOrder("y"           );
        z_order            = getEntryOrder("z"           );
        ux_order           = getEntryOrder("ux"          );
        uy_order           = getEntryOrder("uy"          );
        uz_order           = getEntryOrder("uz"          );
}
// ==============================================================
BBand::~BBand(){}
// ==============================================================
int BBand::getBarKey(int index)
{
	int s = getSector   (index);
	int l = getLayer    (index);
	int c = getComponent(index);
	return s*100+l*10+c;
}
