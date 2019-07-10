#include "BBand.h"
#include "TVector3.h"
// ==============================================================
BBand::BBand(hipo::schema __schema) : hipo::bank(__schema){

	id_order           = __schema.getEntryOrder("id"          );
	sector_order       = __schema.getEntryOrder("sector"      );
	layer_order        = __schema.getEntryOrder("layer"       );
	component_order    = __schema.getEntryOrder("component"   );
	meantimeTdc_order  = __schema.getEntryOrder("meantimeTdc" );
	meantimeFadc_order = __schema.getEntryOrder("meantimeFadc");
	difftimeTdc_order  = __schema.getEntryOrder("difftimeTdc" );
	difftimeFadc_order = __schema.getEntryOrder("difftimeFadc");
	adcLcorr_order     = __schema.getEntryOrder("adcLcorr"    );
	adcRcorr_order     = __schema.getEntryOrder("adcRcorr"    );
	tFadcLcorr_order   = __schema.getEntryOrder("tFadcLcorr"  );
	tFadcRcorr_order   = __schema.getEntryOrder("tFadcRcorr"  );
	tTdcLcorr_order    = __schema.getEntryOrder("tTdcLcorr"   );
	tTdcRcorr_order    = __schema.getEntryOrder("tTdcRcorr"   );
	x_order            = __schema.getEntryOrder("x"           );
	y_order            = __schema.getEntryOrder("y"           );
	z_order            = __schema.getEntryOrder("z"           );
	ux_order           = __schema.getEntryOrder("ux"          );
	uy_order           = __schema.getEntryOrder("uy"          );
	uz_order           = __schema.getEntryOrder("uz"          );

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
