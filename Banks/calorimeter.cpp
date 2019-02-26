#include "calorimeter.h"
#include "TVector3.h"
// ==============================================================
void   calorimeter::init(const char *bankName, hipo::reader &r){
	initBranches(bankName,r);
	pindex_order   = getEntryOrder("pindex"   );
	detector_order = getEntryOrder("detector" );
	sector_order   = getEntryOrder("sector"   );
	layer_order    = getEntryOrder("layer"    );
	energy_order   = getEntryOrder("energy"   );
	time_order     = getEntryOrder("time"     );
	path_order     = getEntryOrder("path"     );
	x_order        = getEntryOrder("x"        );
	y_order        = getEntryOrder("y"        );
	z_order        = getEntryOrder("z"        );
	lu_order       = getEntryOrder("lu"       );
	lv_order       = getEntryOrder("lv"       );
	lw_order       = getEntryOrder("lw"       );
}
// ==============================================================
calorimeter::~calorimeter(){}
// ==============================================================
float calorimeter::getPcalE (int index)
{
	int nCal = getSize();
	for(int i = 0 ; i < nCal ; i++){
		if(calorimeter::getIndex(i)==index&&calorimeter::getLayer(i)==1){
			return calorimeter::getEnergy(i);
		}
	}
	return 0;
}
// ==============================================================
float calorimeter::getECinE (int index)
{
	int nCal = getSize();
        for(int i = 0 ; i < nCal ; i++){
                if(calorimeter::getIndex(i)==index&&calorimeter::getLayer(i)==4){
                        return calorimeter::getEnergy(i);
                }
        }
        return 0;
}
// ==============================================================
float calorimeter::getECoutE(int index)
{
	int nCal = getSize();
        for(int i = 0 ; i < nCal ; i++){
                if(calorimeter::getIndex(i)==index&&calorimeter::getLayer(i)==7){
                        return calorimeter::getEnergy(i);
                }
        }
        return 0;
}
// ==============================================================
float calorimeter::getTotE  (int index)
{
	return calorimeter::getPcalE(index) + calorimeter::getECinE(index) + calorimeter::getECoutE(index);
}
// ==============================================================
