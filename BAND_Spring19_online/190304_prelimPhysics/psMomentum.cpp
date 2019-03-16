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


int main(int argc, char** argv) {

	TApplication *myapp = new TApplication("myapp",0,0);

	// check number of arguments
	if( argc != 2 ){
		cerr << "Incorrect number of arguments. Instead use:\n\t./psMomentum [inputFile]\n";
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
	TFile * outFile = new TFile("clas_band_physics.root","RECREATE");
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
		// CLAS variables
	double theta_e, phi_e;
	double theta_q, phi_q;
	double Q2, nu, xB, W2, q;
	
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

				double cutMeVee = 5.;
				if( adcLcorr < cutMeVee*2000. || adcRcorr < cutMeVee*2000. ) continue;

				//if( sector == 3 || sector == 4 ) continue; 	// Don't want short bars for now because they have a systematic
										// shift from long bars due to length of bar
										
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
		}
		outTree->Fill();


	}// end file

	
	// Fit the background part of the ToF spectrum
	TF1 * bkg = new TF1("bkg","pol0",120,200);
	hToF_full -> Fit(bkg,"QESRN");
	double bkg_per_bin = bkg->GetParameter(0);
	cout << "Background level of ToF spectrum: " << bkg_per_bin << "\n";
	double bkg_int = bkg_per_bin * (200.-120.);
	
	// Calculate the total signa with uncertainty:
	double SpB = hToF_sig->Integral();
	double Rb = (200.-120.);
	double Rs = 40.;
	double B = bkg_int * (Rs/Rb);
	double S = SpB - B;
	double Serr = sqrt( 1./SpB+ 1./B * pow( Rs/Rb , 2) );
	cout << "Signal + Background Integral: " << SpB << "\n";
	cout << "Signal and Error: " << S << " " << Serr << "\n";

	// Take our ToF histogram and subtract of BKG level
		// integral of background level
	TH1F * hToF_new = new TH1F("hToF_new","hToF",40,20,60);
	for( int bin = 1 ; bin < hToF_sig->GetXaxis()->GetNbins() ; bin++ ){
		double val = hToF_sig->GetBinContent(bin);
		double newVal = val - bkg_per_bin;
		if( newVal < 0 ) newVal = 0;
		hToF_new->SetBinContent( bin , newVal );
		// Error is 1/sqrt(N) propagated
		double err = sqrt( 1./val + 1./bkg_per_bin );
		hToF_new->SetBinError( bin , err );
	}
	TCanvas * c1 = new TCanvas("c1","ToF Background Subtraction",900,900);
	c1->Divide(1,3);
	c1->cd(1);
	hToF_full->GetYaxis()->SetRange(0,2000);
	hToF_full->SetStats(0);
	hToF_full->SetTitle("");
	hToF_full->Draw();
	//bkg->Draw("same");
	c1->cd(2);
	hToF_sig->SetStats(0);
	hToF_new->GetYaxis()->SetRange(0,1000);
	hToF_sig->SetTitle("");
	hToF_sig->Draw();
	c1->cd(3);
	hToF_new->GetYaxis()->SetRange(0,600);
	hToF_new->SetStats(0);
	hToF_new->SetTitle("");
	hToF_new->Draw("E HIST");
	c1->Update();


	outFile->cd();
	c1->Write();
	outTree->Write();
	outFile->Close();


	//myapp -> Run();
	return 0;
}


