#ifndef BCALORIMETER_H
#define BCALORIMETER_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BCalorimeter : public hipo::bank {

	private:
		int pindex_order   ;
		int detector_order ;
		int sector_order   ;
		int layer_order    ;
		int energy_order   ;
		int time_order     ;
		int path_order     ;
		int x_order        ;
		int y_order        ;
		int z_order        ;
		int lu_order       ; 
		int lv_order       ;
		int lw_order       ;

	public:

		BCalorimeter(){};

		BCalorimeter(const char *bankName, hipo::reader &r) : hipo::bank(bankName,r){
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

		~BCalorimeter();

		void  init(const char *bankName, hipo::reader &r);
		int   getIndex   (int index) { return getInt  (pindex_order    ,index);}
		int   getDetector(int index) { return getInt  (detector_order  ,index);}
		int   getSector  (int index) { return getInt  (sector_order    ,index);}
		int   getLayer   (int index) { return getInt  (layer_order     ,index);}
		float getEnergy  (int index) { return getFloat(energy_order    ,index);}
		float getLU      (int index) { return getFloat(lu_order        ,index);}
		float getLV      (int index) { return getFloat(lv_order        ,index);}
		float getLW      (int index) { return getFloat(lw_order        ,index);}

		float getPcalE (int index);
		float getECinE (int index);
		float getECoutE(int index);
		float getTotE  (int index);
};

#endif
