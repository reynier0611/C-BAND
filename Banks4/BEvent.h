#ifndef BEVENT_H
#define BEVENT_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BEvent : public hipo::bank {

	private:

		int Category_order  ;
		int Topo_order      ;
		int BCG_order       ;
		int LiveTime_order  ;
		int StartTime_order ;
		int RFTime_order    ;
		int Helic_order     ;
		int HelicRaw_order  ;
		int ProcTime_order  ;


	public:

		BEvent(){};

		BEvent(const char *bankName, hipo::reader &r) : hipo::bank(bankName,r){

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

		~BEvent();

		void  init(const char *bankName, hipo::reader &r);


		int   getCategory(int index) { return getInt   ( Category_order   ,index);}
		int   getTopo    (int index) { return getInt   ( Topo_order       ,index);}
		float getBCG     (int index) { return getFloat ( BCG_order        ,index);}
		float getLT      (int index) { return getFloat ( LiveTime_order   ,index);}
		float getSTTime  (int index) { return getFloat ( StartTime_order  ,index);}
		float getRFTime  (int index) { return getFloat ( RFTime_order     ,index);}
		int   getHelic   (int index) { return getInt   ( Helic_order      ,index);}
		int   getHelicRaw(int index) { return getInt   ( HelicRaw_order   ,index);}
		float getPTime   (int index) { return getFloat ( ProcTime_order   ,index);}
};

#endif
