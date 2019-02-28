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

		int NRUN_order    ;
		int NEVENT_order  ;
		int EVNTime_order ;
		int TYPE_order    ;
		int EvCAT_order   ;
		int NPGP_order    ;
		int TRG_order     ;
		int BCG_order     ;
		int LT_order      ;
		int STTime_order  ;
		int RFTime_order  ;
		int Helic_order   ;
		int PTIME_order   ;

	public:

		BEvent(){};

		BEvent(const char *bankName, hipo::reader &r) : hipo::bank(bankName,r){

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

		~BEvent();

		void  init(const char *bankName, hipo::reader &r);

		int   getNrun    (int index) { return getInt   ( NRUN_order    ,index);}
		int   getNEvent  (int index) { return getInt   ( NEVENT_order  ,index);}
		float getEvnTime (int index) { return getFloat ( EVNTime_order ,index);}
		int   getType    (int index) { return getInt   ( TYPE_order    ,index);}
		int   getEvCat   (int index) { return getInt   ( EvCAT_order   ,index);}
		int   getNPGP    (int index) { return getInt   ( NPGP_order    ,index);}
		int   getTrg     (int index) { return getInt   ( TRG_order     ,index);}
		float getBCG     (int index) { return getFloat ( BCG_order     ,index);}
		float getLT      (int index) { return getFloat ( LT_order      ,index);}
		float getSTTime  (int index) { return getFloat ( STTime_order  ,index);}
		float getRFTime  (int index) { return getFloat ( RFTime_order  ,index);}
		int   getHelic   (int index) { return getInt   ( Helic_order   ,index);}
		float getPTime   (int index) { return getFloat ( PTIME_order   ,index);}
};

#endif
