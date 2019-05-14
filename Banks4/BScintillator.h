#ifndef BSCINTILLATOR_H
#define BSCINTILLATOR_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BScintillator : public hipo::bank {

	private:

		int index_order     ;
		int pindex_order    ;
		int detector_order  ;
		int sector_order    ;
		int layer_order     ;
		int component_order ;
		int energy_order    ;
		int time_order      ;
		int path_order      ;
		int chi2_order      ;
		int x_order         ;
		int y_order         ;
		int z_order         ;
		int hx_order        ;
		int hy_order        ;
		int hz_order        ;
		int status_order    ;

	public:

		BScintillator(){};

		BScintillator(const char *bankName, hipo::reader &r) : hipo::bank(bankName,r){

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

		~BScintillator();

		void  init(const char *bankName, hipo::reader &r);

		int   getIndex	  (int index) { return getInt   ( index_order     ,index);}
		int   getPindex   (int index) { return getInt   ( pindex_order    ,index);}
		int   getDetector (int index) { return getInt   ( detector_order  ,index);}
		int   getSector   (int index) { return getInt   ( sector_order    ,index);}
		int   getLayer    (int index) { return getInt   ( layer_order     ,index);}
		int   getComponent(int index) { return getInt   ( component_order ,index);}
		float getEnergy   (int index) { return getFloat ( energy_order    ,index);}
		float getTime     (int index) { return getFloat ( time_order      ,index);}
		float getPath     (int index) { return getFloat ( path_order      ,index);}
		float getChi2     (int index) { return getFloat ( chi2_order      ,index);}
		float getX        (int index) { return getFloat ( x_order         ,index);}
		float getY        (int index) { return getFloat ( y_order         ,index);}
		float getZ        (int index) { return getFloat ( z_order         ,index);}
		float getHx       (int index) { return getFloat ( hx_order        ,index);}
		float getHy       (int index) { return getFloat ( hy_order        ,index);}
		float getHz       (int index) { return getFloat ( hz_order        ,index);}
		int   getStatus   (int index) { return getInt   ( status_order    ,index);}

};

#endif
