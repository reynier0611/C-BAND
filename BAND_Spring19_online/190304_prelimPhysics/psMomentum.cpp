#include <cstdlib>
#include <iostream>

#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
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

	// Declaring histograms
	TH1F * hToF_bkg = new TH1F("hToF_bkg","hToF_bkg",400,0,400);
	TH1F * hToF 	= new TH1F("hToF","hToF",40,0,80);
	TH1F * hToF_pre = new TH1F("hToF_pre","hToF_pre",60,0,120);
	TH1F * hDL 	= new TH1F("hDL","hDL",50,250,350);
	TH1F * hQ2      = new TH1F("hQ2","hQ2",100,0,10);
        TH1F * hNu      = new TH1F("hNu","hnu",100,0,10);
        TH1F * hXB      = new TH1F("hXB","xB",100,0,2);
        TH1F * hTheta_q = new TH1F("hTheta_q","htheta_q",100,-180,180);
        TH1F * hqMag    = new TH1F("hqMag","hq",100,0,10);	
	TH1F * hBeta	= new TH1F("hBeta","hBeta",50,0,1);
	TH1F * hPsMom	= new TH1F("hPsMom","hPsMom",10,0,1000);
	TH1F * hEn	= new TH1F("hEn","hEn",50,0,2000);
	TH1F * hTheta_n = new TH1F("hTheta_n","hTheta_n",100,-180,180);
	TH1F * hPhi_n 	= new TH1F("hPhi_n","hPhi_n",100,-180,180);
	TH1F * hPhi_en 	= new TH1F("hPhi_en","hPhi_en",100,-180,180);
	TH1F * hTheta_qn= new TH1F("hTheta_qn","hTheta_qn",100,-180,180);
	TH1F * hXp	= new TH1F("hXp","hXp",100,0,1);
	TH1F * hWp	= new TH1F("hWp","hWp",50,0,5);
	TH1F * hAs	= new TH1F("hAs","hAs",10,1,2);
	TH2D * hXpAs	= new TH2D("hXpAs","hXpAs",14,0.1,0.8,25,1.3,1.55);
	TH2D * hXpQ2	= new TH2D("hXpQ2","hXpQ2",14,0.1,0.8,30,2,8);
	

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
		double theta_e 		= eP.Theta();
		double Q2 		= 4. * e0_4.Energy() * eP.Mag() * pow( TMath::Sin(theta_e/2.) , 2 );	// momentum transfer
		double nu 		= e0_4.Energy() - eP.Mag();						// energy transfer
		double xB		= Q2 / (2.*mP*nu);
		double q		= sqrt( Q2 + pow(nu,2) ); 						// 3 vector q-magnitude
		double phi_q		= eP.Phi() + M_PI;	  						// q vector phi angle
		if (phi_q > M_PI) phi_q -= 2.*M_PI;
		double theta_q 		= acos(  ( e0_4.Energy() - eP.Mag() * TMath::Cos(theta_e) ) / q );	// q vector theta angle
		double phi_e		= eP.Phi();


		// Looking in BAND
		int nHits = band_hits.getSize();
		if( nHits != 1 ) continue;
		for(int hit = 0; hit < nHits; hit++) {

			int    sector          	= band_hits.getSector      	(hit);
			float adcLcorr         	= band_hits.getAdcLcorr    	(hit);
			float adcRcorr         	= band_hits.getAdcRcorr    	(hit);
			float  meantimeFadc    	= band_hits.getMeantimeFadc	(hit);
			float x			= band_hits.getX		(hit);
			float y			= band_hits.getY		(hit);
			float z			= band_hits.getZ		(hit);

			if( adcLcorr < 4000 || adcRcorr < 4000 ) continue;

			if( sector == 3 || sector == 4 ) continue; 	// Don't want short bars for now because they have a systematic
									// shift from long bars due to length of bar
									
			double ToF = meantimeFadc-t_vtx;	// [ns]
			double dL  = sqrt( x*x + y*y + z*z );	// [cm]


			// Subtract off gamma peak position in ToF:
			hToF_bkg -> Fill( ToF + dL/cAir -4.99481e+01);
			ToF = ToF - (4.99481e+01);

			// Fill spectrum that will be used for background subtraction
			//	At t = 100ns from gamma peak, highest p of neutrons could be
			//	100MeV, which we aren't interested in. So let's take from 100ns
			//	for +100ns more to get our background level estimation
			
			double beta = dL/(cAir*(ToF+dL/cAir));	
	
			// Build neutron 4 vector
			double pN_mag = 1./sqrt( 1./(beta*beta) - 1.) * mN; 	// [GeV]
			double En = sqrt( pN_mag*pN_mag + mN*mN );		// [GeV]

			double pN_cosTheta = z/dL;
			double theta_n = acos( pN_cosTheta );
			double phi_n = atan2( y , x );

			TVector3 pN( pN_mag*sin(theta_n)*sin(phi_n) , pN_mag*sin(theta_n)*cos(phi_n) , pN_mag*cos(theta_n) );
			TLorentzVector pN_4( pN , En );
			
			// Now look at electron-neutron quantities to build the physics
			double phi_en = phi_n - phi_e;
			if (phi_en < -M_PI) phi_en += 2.*M_PI;
			if (phi_en >  M_PI) phi_en -= 2.*M_PI;
			double CosTheta_qn = cos(theta_n)*cos(theta_q) - sin(theta_n)*sin(theta_q)*cos(phi_en);
			double Xp = Q2/(2.*( nu*(mD-En) + pN_mag*q*CosTheta_qn));
			double Wp = sqrt((mD*mD) - Q2 + (mP*mP) + 2.*mD*(nu-En) -2.* nu * En + 2.*q*pN_mag*CosTheta_qn);
			double As = (En - pN_mag*CosTheta_qn)/mN;



			if (Q2 < 2)       				continue;
			if (acos(CosTheta_qn)*180./M_PI < 110.)		continue;
			if (Xp > 0.8)      				continue;
			if (Wp < 1.8)     	 			continue;

			hToF_pre -> Fill(ToF+dL/cAir);
			if( beta > 1 ) 					continue;
			if( ToF > 100 && ToF < 250 )			continue;
			if( ToF < 0 )	 				continue; 	// these are below our photon peak, so ignore
			if( ToF < 3*1.19734e+00) 			continue; 	// this is our photon peak, so ignore
			if( pN_mag < 0.2)				continue;
			/*
			if (pN_mag < 0.25) 				continue;
			if (pN_mag > 0.6)  	    			continue;
			*/

			hXpQ2		-> Fill(Xp,Q2);
			hXpAs		-> Fill(Xp,As);
			if (Xp < 0.2)  	   	 			continue;
			// Fill histograms
				// only CLAS
			hQ2		-> Fill(Q2);
			hNu		-> Fill(nu);
			hXB		-> Fill(xB);
			hTheta_q	-> Fill(theta_q*180./M_PI);
			hqMag		-> Fill(q);
				// only BAND
			hToF		-> Fill(ToF+dL/cAir);
			hBeta		-> Fill(beta);
			hDL		-> Fill(dL);
			hPsMom		-> Fill(pN_mag * 1E3); // convert to MeV
			hEn		-> Fill(En * 1E3);
			hTheta_n	-> Fill(theta_n * 180./M_PI);
			hPhi_n		-> Fill(phi_n * 180./M_PI);
				// combo
			hPhi_en		-> Fill(phi_en * 180./M_PI);
			hTheta_qn	-> Fill( acos(CosTheta_qn) * 180./M_PI );
			hXp		-> Fill(Xp);
			hWp		-> Fill(Wp);
			hAs		-> Fill(As);
		
		

		}// end loop over hits in event


	}// end file

	// From hToF_bkg let's calculate our background levels. But not all of these pass the cuts used above
	double bkg_ev = hToF_bkg->Integral();
	double n_bins = (250.-100.);
	double bkgPerBin = bkg_ev / n_bins;

	

	TCanvas *cbkg = new TCanvas("cbkg","Background Estimation",900,900);
	hToF_bkg->Draw();
	cbkg->Update();
	

	TCanvas *c0 = new TCanvas("c0","BAND Physics",900,900);
	c0->Divide(2,2);
	c0->cd(1);
	hToF_pre ->SetLineColor(2);
	hToF_pre -> Draw();
	hToF->Draw("same");
	c0->cd(2);
	hDL->Draw();
	c0->cd(3);	
	hBeta->Draw();
	c0->cd(4);
	hPsMom->Draw();
	c0->Update();

	TCanvas *c1 = new TCanvas("c1","CLAS Physics",900,900);
	c1->Divide(2,2);
	c1->cd(1);
	hQ2->Draw();
	c1->cd(2);
	hXB->Draw();
	c1->cd(3);
	hTheta_q->Draw();
	c1->cd(4);
	hNu->Draw();
	c1->Update();

	TCanvas * c2 = new TCanvas("c2","CLAS and BAND",900,900);
	c2->Divide(2,3);
	c2->cd(1);
	hTheta_qn->Draw();
	c2->cd(2);
	hPhi_en->Draw();
	c2->cd(3);
	hXp->Draw();
	c2->cd(4);
	hWp->Draw();
	c2->cd(5);
	hAs->Draw();
	c2->cd(6);
	hXpAs->Draw("col");
	c2->Update();

	TCanvas * c3 = new TCanvas("c3");
	hPsMom->SetStats(0);
	hPsMom->SetTitle("");
	hPsMom->Draw();
	c3->Update();

	TCanvas * c4 = new TCanvas("c4");
	hXpQ2->SetStats(0);
	hXpQ2->SetTitle("");
	hXpQ2->Draw("col");
	c4->Update();

	
	TCanvas * c5 = new TCanvas("c5");
	hXpAs->SetStats(0);
	hXpAs->SetTitle("");
	hXpAs->Draw("col");
	c5->Update();

	myapp -> Run();
	return 0;
}


