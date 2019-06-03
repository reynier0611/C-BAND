#include <cstdlib>
#include <iostream>

#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

#include "reader.h"
#include "node.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

// Loading files from include directory
#include "constants.h"

using namespace std;
double parA_L[600];		// loaded
double parB_L[600];		// loaded
double parA_R[600];		// loaded
double parB_R[600];		// loaded
double TDC_TDIFF[600];		// loaded
double FADC_TDIFF[600];		// loaded
double TDC_P2P[600];		// loaded
double TDC_L2L[600];		// loaded
double FADC_P2P[600];		// loaded
double FADC_L2L[600];		// loaded
double TDC_VEFF[600];		// loaded
double FADC_VEFF[600];		// loaded
double FADC_ATTEN_LENGTH[600];	// loaded
double globPos[600][3];		
double BARLENGTHS[]  = {163.7,201.9,51.2,51.2,201.9};

double getTriggerPhase( long timeStamp );
void LoadTimeWalk();
void LoadLROffsets();
void LoadPaddleOffsets();
void LoadLayerOffsets();
void LoadVelocityMap();
void LoadAttenuation();
void CreateGeo();

int main(int argc, char** argv) {


	// Load calibration constants from include DIR:
	cout << "Loading constants...\n";
	LoadTimeWalk();
	LoadLROffsets();
	LoadPaddleOffsets();
	LoadLayerOffsets();
	LoadVelocityMap();
	LoadAttenuation();
	CreateGeo();
	cout << "...Done!\n";


	TApplication *myapp = new TApplication("myapp",0,0);

	// check number of arguments
	if( argc != 3 ){
		cerr << "Incorrect number of arguments. Instead use:\n\t./psMomentum [inputFile] [outputFile]\n";
		return -1;
	}

	// Event selection cuts for electron
	double cut_ep      =     2; //GeV
	double cut_chi2pid =     5;
	double cut_min_vz  =   -15; //cm
	double cut_max_vz  =    10; //cm
	double cut_W       =     0; //GeV
	double cut_uvw     =    15; //cm
	double cut_Epcal   = 0.060; //GeV (60 MeV)
	double cut_tof_e   =    10; //ns

	// Create output tree
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TTree * outTree = new TTree("skim","CLAS and BAND Physics");
	// BAND variables
	int nHits;
	int sector, layer, component;
	double adcLcorr, adcRcorr;
	double meantimeFadc, meantimeTdc;
	double difftimeFadc, difftimeTdc;
	double x,y,z;
	double dL, ToF;
	double beta, pN_mag, theta_n, phi_n;
	double En, pN_cosTheta, phi_en;
	double CosTheta_qn, Xp, Wp, As, theta_qn;
	// Raw BAND variables 
	int nADC, nTDC;
	double adcLraw, adcRraw;
	double tTdcLraw, tTdcRraw;
	double tFadcLraw, tFadcRraw;
	// By hand variables using tables
	double byHand_adcL, byHand_adcR;
	double byHand_meantimeFadc, byHand_meantimeTdc;
	double byHand_difftimeFadc, byHand_difftimeTdc;
	double byHand_x, byHand_y, byHand_z;
	double byHand_dL;
	// CLAS variables
	double theta_e, phi_e;
	double theta_q, phi_q;
	double Q2, nu, xB, W2, q;
	double start_time;

	// Branches for BAND
	outTree->Branch("nHits",		&nHits		,	"nHits/I");
	outTree->Branch("sector",		&sector		,	"sector/I");
	outTree->Branch("layer",		&layer		,	"layer/I");
	outTree->Branch("component",		&component	,	"component/I");
	outTree->Branch("adcLcorr",		&adcLcorr	,	"adcLcorr/D");
	outTree->Branch("adcRcorr",		&adcRcorr	,	"adcRcorr/D");
	outTree->Branch("meantimeFadc",		&meantimeFadc	,	"meantimeFadc/D");
	outTree->Branch("meantimeTdc",		&meantimeTdc	,	"meantimeTdc/D");
	outTree->Branch("difftimeFadc",		&difftimeFadc	,	"difftimeFadc/D");
	outTree->Branch("difftimeTdc",		&difftimeTdc	,	"difftimeTdc/D");
	outTree->Branch("x",			&x		,	"x/D");
	outTree->Branch("y",			&y		,	"y/D");
	outTree->Branch("z",			&z		,	"z/D");
	outTree->Branch("dL",			&dL		,	"dL/D");
	outTree->Branch("ToF",			&ToF		,	"ToF/D");
	outTree->Branch("beta",			&beta		,	"beta/D");
	outTree->Branch("pN_mag",		&pN_mag		,	"pN_mag/D");
	outTree->Branch("theta_n",		&theta_n	,	"theta_n/D");
	outTree->Branch("phi_n",		&phi_n		,	"phi_n/D");
	outTree->Branch("En",			&En		,	"En/D");
	outTree->Branch("pN_cosTheta",		&pN_cosTheta	,	"pN_cosTheta/D");
	outTree->Branch("phi_en",		&phi_en		,	"phi_en/D");
	outTree->Branch("CosTheta_qn",		&CosTheta_qn	,	"CosTheta_qn/D");
	outTree->Branch("Xp",			&Xp		,	"Xp/D");
	outTree->Branch("Wp",			&Wp		,	"Wp/D");
	outTree->Branch("As",			&As		,	"As/D");
	outTree->Branch("theta_qn",		&theta_qn	,	"theta_qn/D");
	// Raw branches for BAND
	outTree->Branch("nADC",			&nADC		,	"nADC/I");
	outTree->Branch("nTDC",			&nTDC		,	"nTDC/I");
	outTree->Branch("adcLraw",		&adcLraw	,	"adcLraw/D");
	outTree->Branch("adcRraw",		&adcRraw	,	"adcRraw/D");
	outTree->Branch("tTdcLraw",		&tTdcLraw	,	"tTdcLraw/D");
	outTree->Branch("tTdcRraw",		&tTdcRraw	,	"tTdcRraw/D");
	outTree->Branch("tFadcLraw",		&tFadcLraw	, 	"tFadcLraw/D");
	outTree->Branch("tFadcRraw",		&tFadcRraw	,	"tFadcRraw/D");
	// By-hand branches for BAND
	outTree->Branch("byHand_adcL",		&byHand_adcL	,	"byHand_adcL/D");
	outTree->Branch("byHand_adcR",		&byHand_adcR	, 	"byHand_adcR/D");
	outTree->Branch("byHand_meantimeFadc",	&byHand_meantimeFadc	,	"byHand_meantimeFadc/D");
	outTree->Branch("byHand_meantimeTdc",	&byHand_meantimeTdc	,	"byHand_meantimeTdc/D");
	outTree->Branch("byHand_difftimeFadc",	&byHand_difftimeFadc	,	"byHand_difftimeFadc/D");
	outTree->Branch("byHand_difftimeTdc",	&byHand_difftimeTdc	,	"byHand_difftimeTdc/D");
	outTree->Branch("byHand_x",		&byHand_x	,	"byHand_x/D");
	outTree->Branch("byHand_y",		&byHand_y	,	"byHand_y/D");
	outTree->Branch("byHand_z",		&byHand_z	,	"byHand_z/D");
	outTree->Branch("byHand_dL",		&byHand_dL	,	"byHand_dL/D");
	// Branches for CLAS
	outTree->Branch("theta_e",		&theta_e	,	"theta_e/D");
	outTree->Branch("phi_e",		&phi_e		,	"phi_e/D");
	outTree->Branch("theta_q",		&theta_q	,	"theta_q/D");
	outTree->Branch("phi_q",		&phi_q		,	"phi_q/D");
	outTree->Branch("Q2",			&Q2		,	"Q2/D");
	outTree->Branch("nu",			&nu		,	"nu/D");
	outTree->Branch("xB",			&xB		,	"xB/D");
	outTree->Branch("W2",			&W2		,	"W2/D");
	outTree->Branch("q",			&q		,	"q/D");
	outTree->Branch("STTime",		&start_time	,	"STTime/D");


	// Declaring histograms
	TH1F * hToF_sig 	= new TH1F("hToF_sig","hToF",40,20,60);
	TH1F * hToF_full 	= new TH1F("hToF_full","hToF",450,-100,350);


	// Load input file
	TString inputFile = argv[1];
	hipo::reader reader;
	reader.open(inputFile);
	// Banks for EMC-SRC physics
	BEvent        event       ("REC::Event"       ,reader);
	BParticle     particles   ("REC::Particle"    ,reader);
	BCalorimeter  calo        ("REC::Calorimeter" ,reader);
	BScintillator scintillator("REC::Scintillator",reader);
	BBand         band_hits   ("BAND::hits"       ,reader);
	hipo::bank BAND_ADC  ("BAND::adc"  ,reader);
	hipo::bank BAND_TDC  ("BAND::tdc"  ,reader);
	hipo::bank RUN_config("RUN::config",reader);

	// Setup initial vector for beam
	TVector3 e0(0,0,Ebeam);
	TLorentzVector e0_4(e0,Ebeam);

	// Loop over events in hipo fil
	int event_counter = 0;
	while(reader.next()==true){
		nHits,sector,layer,component,adcLcorr,adcRcorr			= 0.;
		meantimeFadc,meantimeTdc,difftimeFadc,difftimeTdc 		= 0.;
		x,y,z,dL,ToF,beta,pN_mag,theta_n,phi_n,En,pN_cosTheta 		= 0.;
		phi_en,CosTheta_qn,Xp,Wp,As,theta_qn 				= 0.;
		theta_e,phi_e,theta_q,phi_q,Q2,nu,xB,W2,q 			= 0.;

		nADC, nTDC 							= 0.;
		adcLraw, adcRraw, tTdcLraw, tTdcRraw, tFadcLraw, tFadcRraw 	= 0.;
		byHand_adcL, byHand_adcR, byHand_meantimeFadc			= 0.;
		byHand_meantimeTdc, byHand_difftimeFadc, byHand_difftimeTdc	= 0.;
		byHand_x, byHand_y, byHand_z, byHand_dL				= 0.;

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;

		// Debugging print option
		//band_hits.show();
		//event.show();
		//scintillator.show();	

		// Particle bank
		int pid0		= particles.getPid    (0);       // electron candidate id assigned by clas
		TVector3 eP_vertex	= particles.getV3v    (0);       // electron candidate vertex vector
		TVector3 eP 		= particles.getV3P    (0);       // electron candidate momentum vector
		float chr0     		= particles.getCharge (0);       // electron candidate charge
		float eBeta    		= particles.getBeta   (0);       // electron candidate beta = v/c
		float chi2pid  		= particles.getChi2pid(0);       // electron candidate goodness of pid fit
		int eStatus    		= particles.getStatus (0);       // electron candidate status
		// Calorimeter bank
		float Epcal = calo.getPcalE(0); 
		float Ee    = calo.getTotE (0);
		float lU    = calo.getLU   (0);	// electron candidate distance on U-side [cm?]
		float lV    = calo.getLV   (0);	// electron candidate distance on V-side [cm?]
		float lW    = calo.getLW   (0);	// electron candidate distance on W-side [cm?]

		// Event vertex time calibrated from FTOF
		double t_vtx   = event.getSTTime(0);
		double t_e     = scintillator.getTime(0);
		double tof_e  = t_e - t_vtx;		// electron candidate time-of-flight [ns]

		// Only keep events for which the first particle is an electron
		if(             (pid0!=11              )||
				(chr0!=-1              )||
				(chi2pid>=cut_chi2pid  )||
				(eP.Mag()<=cut_ep      )||
				(eP.Mag()>=Ebeam       )||
				(eP_vertex.Z()>cut_max_vz  )||
				(eP_vertex.Z()<cut_min_vz  )||
				(lU<cut_uvw            )||
				(lV<cut_uvw            )||
				(lW<cut_uvw            )||
				(Epcal<cut_Epcal       )||
				//(TMath::Sqrt(W2)<=cut_W)||
				(tof_e<cut_tof_e       )
		  ) continue;


		// Calculate kinematic variables
		theta_e 	= eP.Theta();
		Q2 		= 4. * e0_4.Energy() * eP.Mag() * pow( TMath::Sin(theta_e/2.) , 2 );	// momentum transfer
		nu 		= e0_4.Energy() - eP.Mag();						// energy transfer
		xB		= Q2 / (2.*mP*nu);
		q		= sqrt( Q2 + pow(nu,2) ); 						// 3 vector q-magnitude
		phi_q		= eP.Phi() + M_PI;	  						// q vector phi angle
		W2     		= mP*mP - Q2 + 2*nu*mP;						// invariant jet mass based on e kinematics
		if (phi_q > M_PI) phi_q -= 2.*M_PI;
		theta_q 	= acos(  ( e0_4.Energy() - eP.Mag() * TMath::Cos(theta_e) ) / q );	// q vector theta angle
		phi_e		= eP.Phi();


		// Looking in BAND
		nHits = band_hits.getSize();
		if( nHits == 1){
			for(int hit = 0; hit < nHits; hit++) {

				int    barKey           = band_hits.getBarKey  		(hit); 
				sector          	= band_hits.getSector      	(hit);
				layer			= band_hits.getLayer		(hit);
				component		= band_hits.getComponent	(hit);

				adcLcorr         	= band_hits.getAdcLcorr    	(hit);
				adcRcorr         	= band_hits.getAdcRcorr    	(hit);
				meantimeFadc    	= band_hits.getMeantimeFadc	(hit);
				meantimeTdc		= band_hits.getMeantimeTdc	(hit);	

				difftimeFadc		= band_hits.getDifftimeFadc	(hit);
				difftimeTdc		= band_hits.getDifftimeTdc	(hit);

				x			= band_hits.getX		(hit);
				y			= band_hits.getY		(hit);
				z			= band_hits.getZ		(hit);

				//double cutMeVee = 5.;
				//if( adcLcorr < cutMeVee*2000. || adcRcorr < cutMeVee*2000. ) continue;

				//if( sector == 3 || sector == 4 ) continue; 	// Don't want short bars for now because they have a systematic
				// shift from long bars due to length of bar
				start_time = t_vtx;					
				ToF = meantimeFadc-t_vtx;	// [ns]
				dL  = sqrt( x*x + y*y + z*z );	// [cm]

				// Real ToF = (measured ToF - gamma peak position) + (where gamma peak should really be)
				ToF = ToF - 4.99481e+01 + dL/cAir;
				hToF_sig->Fill( ToF );
				hToF_full->Fill( ToF );

				// Calculate beta
				beta = dL/(cAir*ToF);	

				// Build neutron 4 vector
				pN_mag = 1./sqrt( 1./(beta*beta) - 1.) * mN; 	// [GeV]
				En = sqrt( pN_mag*pN_mag + mN*mN );		// [GeV]

				pN_cosTheta = z/dL;
				theta_n = acos( pN_cosTheta );
				phi_n = atan2( y , x );

				TVector3 pN( pN_mag*sin(theta_n)*sin(phi_n) , pN_mag*sin(theta_n)*cos(phi_n) , pN_mag*cos(theta_n) );
				TLorentzVector pN_4( pN , En );

				// Now look at electron-neutron quantities to build the physics
				phi_en = phi_n - phi_e;
				if (phi_en < -M_PI) phi_en += 2.*M_PI;
				if (phi_en >  M_PI) phi_en -= 2.*M_PI;
				CosTheta_qn = cos(theta_n)*cos(theta_q) - sin(theta_n)*sin(theta_q)*cos(phi_en);
				theta_qn = acos(CosTheta_qn);
				Xp = Q2/(2.*( nu*(mD-En) + pN_mag*q*CosTheta_qn));
				Wp = sqrt((mD*mD) - Q2 + (mP*mP) + 2.*mD*(nu-En) -2.* nu * En + 2.*q*pN_mag*CosTheta_qn);
				As = (En - pN_mag*CosTheta_qn)/mN;



				//if (Q2 < 2)       				continue;
				//if (acos(CosTheta_qn)*180./M_PI < 110.)		continue;
				//if (Xp > 0.8)      				continue;
				//if (Wp < 1.8)     	 			continue;

				//if( beta > 1 ) 					continue;
				//if( ToF > 100 && ToF < 250 )			continue;
				//if( ToF < 0 )	 				continue; 	// these are below our photon peak, so ignore
				//if( ToF < 3*1.19734e+00) 			continue; 	// this is our photon peak, so ignore
				//if( pN_mag < 0.2)				continue;
				//if (pN_mag < 0.25) 				continue;
				//if (pN_mag > 0.6)  	    			continue;


			}// end loop over band hits in an event

			// Now check the corresponding ADC banks -- we should only have 2 ADCs, 2 TDCs:
			nADC = BAND_ADC.getSize();
			nTDC = BAND_TDC.getSize();
			if( nADC == 2 && nTDC == 2){
				//RUN_config.show();
				//long timestamp = RUN_config.getLong(4,0);
				//double phaseCorr = getTriggerPhase(timestamp);
				int adc_barKey, tdc_barKey;

				// Get the raw ADC information, uncorrected
				//BAND_ADC.show();
				for(int aIdx = 0 ; aIdx < nADC ; aIdx++){
					int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
					int   ADC_sector    = BAND_ADC.getInt  (0,aIdx);
					int   ADC_layer     = BAND_ADC.getInt  (1,aIdx);
					int   ADC_component = BAND_ADC.getInt  (2,aIdx);
					if( ADC_order == 0 ){
						adcLraw = (float)(BAND_ADC.getInt(4,aIdx));
						tFadcLraw = BAND_ADC.getFloat(5,aIdx);
						adc_barKey = ADC_sector*100 + ADC_layer*10 + ADC_component;
					}
					if( ADC_order == 1 ){
						adcRraw = (float)(BAND_ADC.getInt(4,aIdx));
						tFadcRraw = BAND_ADC.getFloat(5,aIdx);
						adc_barKey = ADC_sector*100 + ADC_layer*10 + ADC_component;
					}
				}

				// Get the raw TDC information, uncorrected
				//BAND_TDC.show();
				for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
					int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
					int   TDC_sector    = BAND_TDC.getInt  (0,tIdx);
					int   TDC_layer     = BAND_TDC.getInt  (1,tIdx);
					int   TDC_component = BAND_TDC.getInt  (2,tIdx);
					TDC_order -= 2;
					if( TDC_order == 0 ){
						tTdcLraw = (float)(BAND_TDC.getInt(4,tIdx))*0.02345;
						tdc_barKey = TDC_sector*100 + TDC_layer*10 + TDC_component;
					}
					if( TDC_order == 1 ){
						tTdcRraw = (float)(BAND_TDC.getInt(4,tIdx))*0.02345;
						tdc_barKey = TDC_sector*100 + TDC_layer*10 + TDC_component;
					}

				}

				// Correct everything by hand using tables in include DIR
				if( 	adcLraw != 0. && adcRraw != 0. && 			// ADC non zero
						tFadcLraw != 0. && tFadcRraw != 0. &&			// ADC time non zero
						tTdcLraw != 0. && tTdcRraw != 0. && 			// TDC time non zero
						adc_barKey == sector*100+layer*10+component && 		// matching bar ID for ADC
						tdc_barKey == sector*100+layer*10+component ){		/// matching bar ID for TDC

					int barID = sector*100+layer*10+component;					

					// TDC phase correction:
					//tTdcLraw -= phaseCorr;
					//tTdcRraw -= phaseCorr;
					// TDC time walk
					tTdcLraw = tTdcLraw - (parA_L[barID]/sqrt(adcLraw) + parB_L[barID]);
					tTdcRraw = tTdcRraw - (parA_R[barID]/sqrt(adcRraw) + parB_R[barID]);
					//cout << sector << " " << layer << " " << component << " " 
					//	<< parA_L[barID] << " " << parB_L[barID] << " "
					//	<< parA_R[barID] << " " << parB_R[barID] << " "
					//	<< TDC_TDIFF[barID] << " " << FADC_TDIFF[barID] << " "
					//	<< TDC_P2P[barID] << " " << TDC_L2L[barID] << " " 
					//	<< FADC_P2P[barID] << " " << FADC_L2L[barID] << " "
					//	<< TDC_VEFF[barID] << " " << FADC_VEFF[barID] << " "
					//	<< FADC_ATTEN_LENGTH[barID] << " "
					//	<< globPos[barID][0] << " " << globPos[barID][1] << " " << globPos[barID][2] << "\n";
													
					// TDiff:
					byHand_difftimeTdc = (tTdcLraw-tTdcRraw) - TDC_TDIFF[barID];
					byHand_difftimeFadc = (tFadcLraw-tFadcRraw) - FADC_TDIFF[barID];
					// Meantime:
					byHand_meantimeTdc = (tTdcLraw+tTdcRraw)/2. - TDC_TDIFF[barID]/2. - TDC_P2P[barID] - TDC_L2L[barID];
					byHand_meantimeFadc = (tFadcLraw+tFadcRraw)/2. - FADC_TDIFF[barID]/2. - FADC_P2P[barID] - FADC_L2L[barID];
					// ADC attenuation corr:
					double xpos_tdc = (-1./2) * byHand_difftimeTdc * TDC_VEFF[barID];
					double xpos_fadc = (-1./2) * byHand_difftimeFadc * FADC_VEFF[barID];
					double sectorLen = BARLENGTHS[sector-1];
					double mu_cm = FADC_ATTEN_LENGTH[barID];
					byHand_adcL = adcLraw * exp( (sectorLen/2.-xpos_fadc) / mu_cm );
					byHand_adcR = adcRraw * exp( (sectorLen/2.+xpos_fadc) / mu_cm );
					// Path length:
					byHand_x = (xpos_tdc+xpos_fadc)/2.;
					byHand_x += globPos[barID][0];
					byHand_y = globPos[barID][1];
					byHand_z = globPos[barID][2];
					byHand_dL = sqrt( byHand_x*byHand_x + byHand_y*byHand_y + byHand_z*byHand_z );

				}
			} // end if for raw ADC==2 and TDC =2

		} // end if for nHits == 1 for BAND

		outTree->Fill();


	}// end file



	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}

double getTriggerPhase( long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}

void LoadTimeWalk(){
	ifstream f;
	int sector, layer, component, barId;
	double parA, parB, temp;

	f.open("../include/time_walk_corr_left.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> temp;
		f >> temp;
		parA_L[barId] = parA;
		parB_L[barId] = parB;
	}
	f.close();

	f.open("../include/time_walk_corr_right.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parA;
		f >> parB;
		f >> temp;
		f >> temp;
		parA_R[barId] = parA;
		parB_R[barId] = parB;
	}
	f.close();
	return;
}
void LoadLROffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double tdc_off, fadc_off, temp;

	f.open("../include/lr_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> tdc_off;
		f >> fadc_off;
		f >> temp;
		f >> temp;
		TDC_TDIFF[barId] = tdc_off;
		FADC_TDIFF[barId] = fadc_off;
	}
	f.close();
	return;
}
void LoadPaddleOffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double offset_fadc, offset_tdc, temp;

	f.open("../include/paddle_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_fadc;
		f >> temp;
		FADC_P2P[barId] = offset_fadc;
	}
	f.close();
	f.open("../include/paddle_offsets_tdc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_tdc;
		f >> temp;
		TDC_P2P[barId] = offset_tdc;
	}
	f.close();

	return;
}
void LoadLayerOffsets(){
	ifstream f;
	int sector, layer, component, barId;
	double offset_fadc, offset_tdc, temp;

	f.open("../include/layer_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_fadc;
		f >> temp;
		FADC_L2L[barId] = offset_fadc;
	}
	f.close();
	f.open("../include/layer_offsets_tdc.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> offset_tdc;
		f >> temp;
		TDC_L2L[barId] = offset_tdc;
	}
	f.close();

	return;
}
void LoadVelocityMap(){
	ifstream f;
	int sector, layer, component, barId;
	double veff_tdc, veff_fadc, temp;

	f.open("../include/effective_velocity.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> veff_tdc;
		f >> veff_fadc;
		f >> temp;
		f >> temp;
		TDC_VEFF[barId] = veff_tdc;
		FADC_VEFF[barId] = veff_fadc;
	}
	f.close();

	return;
}
void LoadAttenuation(){
	ifstream f;
	int sector, layer, component, barId;
	double atten_len, temp;

	f.open("../include/attenuation_lengths.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> atten_len;
		f >> temp;
		FADC_ATTEN_LENGTH[barId] = atten_len;
	}
	f.close();

	return;
}

void CreateGeo(){
	// GEOMETRY PARAMETERS
	int sectNum = 5;	// Number of sectors (blocks)
	int layNum = 6; 	// Number of layer
	int compNum = 7;    // Maximum Number of components in a sector
	// Number of components per layer and sector [layer][sector] index
	int compNumSecLay[6][5] = { {3,7,6,6,2}, {3,7,6,6,2}, {3,7,6,6,2}, {3,7,6,6,2}, {3,7,5,5,0}, {3,7,6,6,2} };

	double distVetoLead = 17.463;			// distance from veto to lead wall (cm)

	double thickness = 7.2;				// thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};  // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6
	double zOffset = 100;                               // distance from center first layer to target.
	double surveyBox[4][3] = {  	{-24.05,-21.10,-302.69},
					{-24.05, 22.81,-302.69},
					{ 24.10,-21.06,-302.57},
					{ 24.37, 22.81,-302.64}  	};
	double lenLG = 8.9; // [cm] -- length of LG
	double lenET = 16.; // [cm] -- length of PMT tube for ET PMTs
	double lenHam = 13.3; // [cm] -- length of PMT tube for Hamamatsu PMTs

	double avgX = ( (surveyBox[0][0] + surveyBox[2][0]) + (surveyBox[1][0] + surveyBox[3][0]) )/2.;
	double avgY = ( (surveyBox[0][1] + thickness*3.) + (surveyBox[1][1] - thickness*3.) 
			+ (surveyBox[2][1] + thickness*3.) + (surveyBox[3][1] - thickness*3.) )/4.;
	double avgZ = ( surveyBox[0][2] + surveyBox[1][2] + surveyBox[2][2] + surveyBox[3][2] )/4.;

	double globPt[] = {avgX,avgY,avgZ}; // single global position

	double barLengthSector[] = {164, 202, 51, 51, 202} ;           // Bar length in each layer (cm)


	for( int layer = 1 ; layer < layNum + 1 ; layer++){
		double localZ = 0.;
		localZ += (layerGap[1])/2.; // taking this thickness because wrapping material
		// isn't 'squeezed' by the weight of the detector
		for( int i = 1 ; i < layer ; i++ ){
			localZ += layerGap[i-1];
		}

		for( int sector = 1 ; sector < sectNum + 1 ; sector++){
			int nBars = compNumSecLay[layer-1][sector-1];
			for( int bar = 1 ; bar < nBars + 1 ; bar++){
				int key = sector*100+layer*10+bar;
				double localY = 666666.;				
				double localX = 666666.;

				double secYOff = 666666.;

				if( sector == 1){
					secYOff = 10.;
					localX = 0.;
				}
				else if( sector == 2){
					secYOff = 3.;
					localX = 0.;
				}
				else if( sector == 3 || sector == 4){
					secYOff = -3.;
					if( sector == 3){
						localX = (surveyBox[2][0]+surveyBox[3][0])/2. + lenLG + lenET + barLengthSector[sector-1]/2.;
					}
					else if( sector == 4){
						localX = (surveyBox[0][0]+surveyBox[1][0])/2. - lenLG - lenET - barLengthSector[sector-1]/2.;
					}
				}
				else if( sector == 5){
					secYOff = -5.;
					localX = 0.;
				}
				localY = secYOff*thickness + (nBars - (bar-1) )*thickness - thickness/2.;

				globPos[key][0] =  localX + globPt[0];
				globPos[key][1] =  localY + globPt[1];
				globPos[key][2] =  localZ + globPt[2];
				
			}


		}


	}
	return;
}
