#include <cstdlib>
#include <iostream>

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
#include "BBand.h"

using namespace std;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc);

// ========================================================================================================================================
int main(int argc, char** argv) {

	std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

	TString inputFile;

	if(argc==2) {
		inputFile = argv[1];
	}
	else {
		cout << "=========================\nRun this code as:\n./code path/to/input/file\n=========================" << endl;
		exit(0);
	}

	// ----------------------------------------------------------------------------------
	// Useful variables
	double mp      = 0.93827; //GeV
	double mD      = 1.8756;  //GeV (deuteron mass    )
	double rad2deg = 180./3.14159;
	double c = 29.9792;

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
	// Declaring histograms
	// 1D histograms
	TH1F * h1_tdiff_fadc    = new TH1F("h1_tdiff_fadc"   ,"",100 ,-40,   40);	PrettyTH1F(h1_tdiff_fadc    ,"L-R for FADC"                 ,"Counts",1);
	TH1F * h1_tdiff_tdc     = new TH1F("h1_tdiff_tdc"    ,"",100 ,-40,   40);	PrettyTH1F(h1_tdiff_tdc     ,"L-R for TDC"                  ,"Counts",1);
	TH1F * h1_adcL          = new TH1F("h1_adcL"         ,"",1000,  0,30000);	PrettyTH1F(h1_adcL          ,"L ADC"                        ,"Counts",1);
	TH1F * h1_adcR          = new TH1F("h1_adcR"         ,"",1000,  0,30000);	PrettyTH1F(h1_adcR          ,"R ADC"                        ,"Counts",1);

	TH1F * h1_meantime_fadc = new TH1F("h1_meantime_fadc","",400 ,-50,  350);	PrettyTH1F(h1_meantime_fadc ,"(L+R)/2 (FADC)"               ,"Counts",1);
	TH1F * h1_meantime_tdc  = new TH1F("h1_meantime_tdc" ,"",400 ,-50,  350);	PrettyTH1F(h1_meantime_tdc  ,"(L+R)/2 corrected (FADC)"     ,"Counts",1);
	TH1F * h1_ToF_fadc      = new TH1F("h1_ToF_fadc"     ,"",400 ,-50,  350);	PrettyTH1F(h1_ToF_fadc      ,"(L+R)/2 - RF (FADC)"          ,"Counts",1);
	TH1F * h1_ToF_tdc       = new TH1F("h1_ToF_tdc"      ,"",400 ,-50,  350);	PrettyTH1F(h1_ToF_tdc       ,"(L+R)/2 corrected - RF (FADC)","Counts",1);

	TH1F * h1_bar_nPho      = new TH1F("h1_bar_nPho"     ,"",652 ,110,  652);	PrettyTH1F(h1_bar_nPho      ,"Bar ID","Number Hits Between 5-15 ns"  ,1);	

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	double Ebeam = 10.6; //GeV
	TVector3 V3_Ebeam(0,0,Ebeam);
	TLorentzVector V4_Ebeam(V3_Ebeam,Ebeam);

	double mtar    = mD;
	TLorentzVector V4_mtar(0,0,0,mtar);

	hipo::reader reader;
	reader.open(inputFile);

	hipo::bank    bank_scintillator("REC::Scintillator",reader);
	hipo::bank    bank_event       ("REC::Event"       ,reader);
	BParticle     particles        ("REC::Particle"    ,reader);
	BCalorimeter  calo             ("REC::Calorimeter" ,reader);
	BBand         band_hits        ("BAND::hits"       ,reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;

		//band_hits.show();
		//bank_event.show();
		//bank_scintillator.show();	

		// Particle bank
		int pid0       = particles.getPid    (0);       // electron candidate id assigned by clas
		TVector3 V3_ev = particles.getV3v    (0);       // electron candidate vertex vector
		TVector3 V3_ep = particles.getV3P    (0);       // electron candidate momentum vector
		float chr0     = particles.getCharge (0);       // electron candidate charge
		float eBeta    = particles.getBeta   (0);       // electron candidate beta = v/c
		float chi2pid  = particles.getChi2pid(0);       // electron candidate goodness of pid fit
		int eStatus    = particles.getStatus (0);       // electron candidate status

		// Calorimeter bank     
		float Epcal = calo.getPcalE(0);
		float Ee    = calo.getTotE (0);
		float lU    = calo.getLU   (0); // electron candidate distance on U-side [cm?]
		float lV    = calo.getLV   (0); // electron candidate distance on V-side [cm?]
		float lW    = calo.getLW   (0); // electron candidate distance on W-side [cm?]

		if(Ee==0) continue;

		// Event bank
		double t_vtx   = bank_event.getFloat( 9,0);

		// Scintillator bank
		double t_e     = bank_scintillator.getFloat(7,0);

		// calculated variables
		double ep     = V3_ep.Mag();            // electron candidate momentum magnitude [GeV]
		TLorentzVector V4_ep(V3_ep,ep);
		double tof_e  = t_e - t_vtx;            // electron candidate time-of-flight [ns]

		// Transfer variables
		TVector3 V3_q = V3_Ebeam - V3_ep;
		double Q2     = 4*ep*Ebeam*pow(TMath::Sin(V3_ep.Theta()/2.),2);       // Q-squared [GeV^2]
		double omega  = Ebeam - ep;                                              // Transfer energy [GeV]
		double W2     = mp*mp-Q2+2*omega*mp;
		double xB     = Q2/2./mp/omega;

		TLorentzVector V4_q(V3_q,omega);

		// -------------------------------------------------------------------------
		// Electron PID from Dan Carman
		// - (DONE)     pid=11 from EB
		// - (DONE)     p > 2 GeV
		// - (DONE)     p < Ebeam
		// - (DONE)     TOF > 10 ns //(need bank 330) event.json
		// - (DONE)     vz: -15 to 10 cm
		// - (DONE)     W > 0 GeV
		// -            Sampling fraction +/-3sigma about E/p vs. p mean
		// - (DONE)     ~15 cm fiducial cuts on U, V, W to contain full shower (need bank 332 lu, lv, lw)
		// - (DONE)     abs(chisq PID) < 5 (goodness of PID from EB)
		// - (DONE)     PCAL > 60 MeV (to remove min-i) (bank 332 layers 1(PCAL), 4(EC inner), 7(EC outter))
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

		int nHits = band_hits.getSize();
		if( nHits > 1 ) continue;
		for(int hit = 0; hit < nHits; hit++) {

			int    sector            = band_hits.getSector      (hit);
			int    layer             = band_hits.getLayer       (hit);
			int    component         = band_hits.getComponent   (hit);
			int    barKey            = band_hits.getBarKey      (hit); 
			float  meantimeTdc       = band_hits.getMeantimeTdc (hit);
			float  meantimeFadc      = band_hits.getMeantimeFadc(hit);
			float difftimeTdc        = band_hits.getDifftimeTdc (hit);
			float difftimeFadc       = band_hits.getDifftimeFadc(hit);
			float adcLcorr           = band_hits.getAdcLcorr    (hit);
			float adcRcorr           = band_hits.getAdcRcorr    (hit);

			double meantimeFadc_corr = (double)meantimeFadc ;

			if( adcLcorr < 4000 || adcRcorr < 4000 ) continue;
			//if( sector < 3 || sector > 4 ) continue;
			if( sector == 3 || sector == 4 ) continue;

			float tFadcLcorr = band_hits.getTfadcLcorr(hit);
			float tFadcRcorr = band_hits.getTfadcRcorr(hit);
			float tTdcLcorr  = band_hits.getTtdcLcorr (hit);
			float tTdcRcorr  = band_hits.getTtdcRcorr (hit);
			float x          = band_hits.getX         (hit);
			float y          = band_hits.getY         (hit);
			float z          = band_hits.getZ         (hit);
			float ux         = band_hits.getUx        (hit);
			float uy         = band_hits.getUy        (hit);
			float uz         = band_hits.getUz        (hit);

			//if( sector == 3 || sector == 4 ) continue;

			double dL = TMath::Sqrt( pow(x,2) + pow(y,2) + pow(z,2) );

			// Fill histograms
			h1_tdiff_fadc    -> Fill(difftimeFadc           );
			h1_tdiff_tdc     -> Fill(difftimeTdc            );
			h1_adcL          -> Fill(adcLcorr               );
			h1_adcR          -> Fill(adcRcorr               );
			h1_meantime_fadc -> Fill(meantimeFadc           );
			h1_meantime_tdc  -> Fill(meantimeFadc_corr      );
			h1_ToF_fadc      -> Fill(meantimeFadc-t_vtx     );
			h1_ToF_tdc       -> Fill(meantimeFadc_corr-t_vtx);

			if( TMath::Abs( meantimeFadc-t_vtx-40 - 10 ) < 5 ){
				h1_bar_nPho -> Fill( barKey );
			}

		}// end loop over hits in event


	}// end file

	TCanvas * c0 = new TCanvas();
	c0 -> Divide(2, 4);
	c0 -> cd(1);	h1_adcL          -> Draw();
	c0 -> cd(2);	h1_adcR          -> Draw();
	c0 -> cd(3);	h1_tdiff_fadc    -> Draw();
	c0 -> cd(4);	h1_tdiff_tdc     -> Draw();
	c0 -> cd(5);	h1_meantime_fadc -> Draw();
	c0 -> cd(6);	h1_meantime_tdc  -> Draw();
	c0 -> cd(7);	h1_ToF_fadc      -> Draw();
	c0 -> cd(8);	h1_ToF_tdc       -> Draw();
	c0 -> Modified();
	c0 -> Update();

	c0 -> Print("results_mean_time.pdf");

	return 0;
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color) {
	h1 -> GetXaxis() -> SetTitle(titx);
	h1 -> GetYaxis() -> SetTitle(tity);
	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(3);
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2,TString titx,TString tity) {
	h2 -> GetXaxis() -> SetTitle(titx);
	h2 -> GetYaxis() -> SetTitle(tity);
}
// ========================================================================================================================================
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc){
	// Calculates missing energy
	// Takes as input: missing momentum, transfer energy, struck nucleon energy, and struck nucleon mass.
	double Tb   = omega + M_tar - Enuc - TMath::Sqrt(pow(omega + M_tar - Enuc,2)- Pmiss*Pmiss );	// Kinetic energy of A-1 system
	double Tnuc = Enuc - Mnuc;                                                             			// Kinetic energy of struck neutron
	return omega - Tnuc - Tb;
}
