#include "BScintillator.h"
#include "TVector3.h"
// ==============================================================
void   BScintillator::init(const char *bankName, hipo::reader &r){
	initBranches(bankName,r);

	index_order     = getEntryOrder("index"    );
	pindex_order    = getEntryOrder("pindex"   );
	detector_order  = getEntryOrder("detector" );
	sector_order    = getEntryOrder("sector"   );
	layer_order     = getEntryOrder("layer"    );
	component_order = getEntryOrder("component");
	energy_order    = getEntryOrder("energy"   );
	time_order      = getEntryOrder("time"     );
	path_order      = getEntryOrder("path"     );
	chi2_order      = getEntryOrder("chi2"     );
	x_order         = getEntryOrder("x"        );
	y_order         = getEntryOrder("y"        );
	z_order         = getEntryOrder("z"        );
	hx_order        = getEntryOrder("hx"       );
	hy_order        = getEntryOrder("hy"       );
	hz_order        = getEntryOrder("hz"       );
	status_order    = getEntryOrder("status"   );
}
// ==============================================================
BScintillator::~BScintillator(){}
// ==============================================================
