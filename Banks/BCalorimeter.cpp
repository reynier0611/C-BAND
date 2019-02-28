#include "BCalorimeter.h"
#include "TVector3.h"
// ==============================================================
void   BCalorimeter::init(const char *bankName, hipo::reader &r){
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
BCalorimeter::~BCalorimeter(){}
// ==============================================================
float BCalorimeter::getPcalE (int index)
{
	int nCal = getSize();
	for(int i = 0 ; i < nCal ; i++){
		if(BCalorimeter::getIndex(i)==index&&BCalorimeter::getLayer(i)==1){
			return BCalorimeter::getEnergy(i);
		}
	}
	return 0;
}
// ==============================================================
float BCalorimeter::getECinE (int index)
{
	int nCal = getSize();
        for(int i = 0 ; i < nCal ; i++){
                if(BCalorimeter::getIndex(i)==index&&BCalorimeter::getLayer(i)==4){
                        return BCalorimeter::getEnergy(i);
                }
        }
        return 0;
}
// ==============================================================
float BCalorimeter::getECoutE(int index)
{
	int nCal = getSize();
        for(int i = 0 ; i < nCal ; i++){
                if(BCalorimeter::getIndex(i)==index&&BCalorimeter::getLayer(i)==7){
                        return BCalorimeter::getEnergy(i);
                }
        }
        return 0;
}
// ==============================================================
float BCalorimeter::getTotE  (int index)
{
	return BCalorimeter::getPcalE(index) + BCalorimeter::getECinE(index) + BCalorimeter::getECoutE(index);
}
// ==============================================================
