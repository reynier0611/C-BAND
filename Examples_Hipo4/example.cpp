#include <cstdlib>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TCanvas.h"

#include "reader.h"
#include "bank.h"
#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BEvent.h"
#include "BBand.h"
#include "BConfig.h"
#include "BScaler.h"

using namespace std;


// ========================================================================================================================================
int main(int argc, char** argv) {

	std::cout << " reading file  (HIPO) "  << __cplusplus << std::endl;

	char inputFile[256];
	char outputFile[256];

	if(argc>1) {
		sprintf(inputFile,"%s",argv[1]);
	} else {
		std::cout << " *** please provide a file name..." << std::endl;
		exit(0);
	}

	// ----------------------------------------------------------------------------------
	// Useful variables
	double mp      = 0.93827; //GeV (proton mass      )
	double mPiC    = 0.13957; //GeV (charged pion mass)
	double mD      = 1.8756;  //GeV (deuteron mass    )
	double rad2deg = 180./3.14159;

	// ----------------------------------------------------------------------------------
	// Event selection cuts
	double cut_ep      =     2; //GeV
	double cut_chi2pid =     5;
	double cut_min_vz  =   -15; //cm
	double cut_max_vz  =    10; //cm
	double cut_W       =     0; //GeV
	double cut_uvw     =    15; //cm
	double cut_Epcal   = 0.060; //GeV (60 MeV)
	double cut_tof_e   =    10; //ns

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	double Ebeam = 10.6; //GeV
	TVector3 V3_Ebeam(0,0,Ebeam);

	double mtar    = mD;

	//Reader for the inputFile
	hipo::reader reader;
	reader.open(inputFile);

	//Read Dictionary of Hipo File
	hipo::dictionary  factory;
	reader.readDictionary(factory);
	//factory.show();

	//Create Banks. Constructor needs now the Schema information from the dictionary. This is one of the
	//main differences than the Hipo3 file reading !!!

	BParticle     particles   (factory.getSchema("REC::Particle"    ));
	BEvent     		event   		(factory.getSchema("REC::Event"    		));
  BCalorimeter  calorimeter (factory.getSchema("REC::Calorimeter" ));
	BScintillator scintillator(factory.getSchema("REC::Scintillator"));
	BBand         band_hits   (factory.getSchema("BAND::hits"       ));
	BConfig       runconfig   (factory.getSchema("RUN::config"       ));

	//In general one can also use hipo::bank directly
	//hipo::bank REC_particle(factory.getSchema("REC::Particle"    ));

  //One also needs a hipo::event object which is called from the reader for each event to get
	//the information for each bank
	hipo::event readevent;

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;
		//Reader has to load information about event in hipo::event class
		reader.read(readevent);

		//Load explicitly all information for each bank for the event
		readevent.getStructure(event);
		readevent.getStructure(particles);
	  readevent.getStructure(calorimeter);
		readevent.getStructure(scintillator);
		readevent.getStructure(band_hits);
		readevent.getStructure(runconfig);

		//Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
		//The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.

		//particles.show();
		//calorimeter.show();
		//event.show();
		//scintillator.show();

		// Particle bank
		int pid0       = particles.getPid    (0);	// electron candidate id assigned by clas
		TVector3 V3_ev = particles.getV3v    (0);	// electron candidate vertex vector
		TVector3 V3_ep = particles.getV3P    (0);	// electron candidate momentum vector
		float chr0     = particles.getCharge (0);	// electron candidate charge
		float eBeta    = particles.getBeta   (0);	// electron candidate beta = v/c
		float chi2pid  = particles.getChi2pid(0);	// electron candidate goodness of pid fit
		int eStatus    = particles.getStatus (0);	// electron candidate status

		// Calorimeter bank
		float Epcal = calorimeter.getPcalE(0);
		float Ee    = calorimeter.getTotE (0);
		float lU    = calorimeter.getLU   (0);	// electron candidate distance on U-side [cm?]
		float lV    = calorimeter.getLV   (0);	// electron candidate distance on V-side [cm?]
		float lW    = calorimeter.getLW   (0);	// electron candidate distance on W-side [cm?]

		if(Ee==0) continue;

		// Event bank
		double t_vtx   = event.getSTTime(0);  // event start time (STTime is mapped to startTime entry of bank)

		// Scintillator bank
		double t_e     = scintillator.getTime(0);

		// calculated variables
		double ep     = V3_ep.Mag();		// electron candidate momentum magnitude [GeV]
		double tof_e  = t_e - t_vtx;		// electron candidate time-of-flight [ns]

		// Transfer variables
		TVector3 V3_q = V3_Ebeam - V3_ep;
		double Q2     = 4*ep*Ebeam*pow(TMath::Sin(V3_ep.Theta()/2.),2);       // Q-squared [GeV^2]
		double nu     = Ebeam - ep;                                              // Transfer energy [GeV]
		double W2     = mp*mp-Q2+2*nu*mp;
		double xB     = Q2/2./mp/nu;

		// -------------------------------------------------------------------------


		// -------------------------------------------------------------------------
		// Electron PID from Dan Carman
		// - (DONE)	pid=11 from EB
		// - (DONE)	p > 2 GeV
		// - (DONE)	p < Ebeam
		// - (DONE)	TOF > 10 ns //(need bank 330) event.json
		// - (DONE)	vz: -15 to 10 cm
		// - (DONE)	W > 0 GeV
		// - 		Sampling fraction +/-3sigma about E/p vs. p mean
		// - (DONE)	~15 cm fiducial cuts on U, V, W to contain full shower (need bank 332 lu, lv, lw)
		// - (DONE)	abs(chisq PID) < 5 (goodness of PID from EB)
		// - (DONE)	PCAL > 60 MeV (to remove min-i) (bank 332 layers 1(PCAL), 4(EC inner), 7(EC outter))
		// -------------------------------------------------------------------------
		// Only keep events for which the first particle is an electron
		if(             (pid0!=11              )||
				(chr0!=-1              )||
				(chi2pid>=cut_chi2pid  )||
				(ep<=cut_ep            )||
				(ep>=Ebeam             )||
				(V3_ev.Z()>cut_max_vz  )||
				(V3_ev.Z()<cut_min_vz  )||
				(lU<cut_uvw            )||
				(lV<cut_uvw            )||
				(lW<cut_uvw            )||
				(Epcal<cut_Epcal       )||
				(TMath::Sqrt(W2)<=cut_W)||
				(tof_e<cut_tof_e       )
		  ) continue;

//This is also new for Hipo4: It was band_hits.getSize() before
			int nHits = band_hits.getRows();
			if( nHits != 1 ) continue;

			float  meantimeTdc       = band_hits.getMeantimeTdc (0);
			float  meantimeFadc      = band_hits.getMeantimeFadc(0);
		// -------------------------------------------------------------------------
	  cout << "Event " << runconfig.getEvent(0) << " , BAND mean time FADC " << meantimeFadc  << " , mean time TDC " <<   meantimeTdc  << endl;


	} // End of while loop over events

	// -------------------------------------------------------------------------------------------

	return 0;
}
