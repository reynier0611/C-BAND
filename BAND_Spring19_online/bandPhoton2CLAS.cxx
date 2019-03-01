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

using namespace std;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc);

// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

	// ----------------------------------------------------------------------------------
        // Useful variables
        double mp      = 0.93827; //GeV (proton mass      )
        double mPiC    = 0.13957; //GeV (charged pion mass)
        double mD      = 1.8756;  //GeV (deuteron mass    )
        double rad2deg = 180./3.14159;
	double c = 29.9792;

        // ----------------------------------------------------------------------------------
        // Getting input arguments
        TString inputFile;
        double Ebeam, mtar;

        if(argc==3){
                if(atoi(argv[1])==1){
                        cout << "Will assume this hipo file corresponds to: Ebeam =  6.4 GeV, target = H (i.e. RGA)" << endl;
                        Ebeam = 6.4; //GeV
                        mtar  = mp;
                }
                else if(atoi(argv[1])==2){
                        cout << "Will assume this hipo file corresponds to: Ebeam = 10.6 GeV, target = D (i.e. RGB)" << endl;
                        Ebeam = 10.6; //GeV
                        mtar  = mD;
                }
                inputFile = argv[2];
        }
        else {
                cout << "=========================\nRun this code as:\n./code A path/to/input/file\n" << endl;
                cout << "where: A = 1 -> Ebeam =  6.4 GeV, target = H (i.e. RGA)" << endl;
                cout << "         = 2 -> Ebeam = 10.6 GeV, target = D (i.e. RGB)" << endl;
                cout << "=========================" << endl;
                exit(0);
        }

        TVector3 V3_Ebeam(0,0,Ebeam);
        TLorentzVector V4_Ebeam(V3_Ebeam,Ebeam);
        TLorentzVector V4_mtar(0,0,0,mtar);

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
	TH1F * h1_ToF_fadc_0 = new TH1F("h1_ToF_fadc_0" ,"",400 ,-50,  350);	PrettyTH1F(h1_ToF_fadc_0 ,"(L+R)/2 - RF (FADC)" ,"Counts",1);
	TH1F * h1_ToF_fadc_1 = new TH1F("h1_ToF_fadc_1" ,"",400 ,-50,  350);	PrettyTH1F(h1_ToF_fadc_1 ,"(L+R)/2 - RF (FADC)" ,"Counts",2);

	TH1F * h1_e_px0  = new TH1F("h1_e_px0"  ,"h1_e_px0"  ,100,  -2,  2);       PrettyTH1F(h1_e_px0   ,"electron p_{x} [GeV]" ,"Counts",62);
	TH1F * h1_e_py0  = new TH1F("h1_e_py0"  ,"h1_e_py0"  ,100,  -2,  2);       PrettyTH1F(h1_e_py0   ,"electron p_{y} [GeV]" ,"Counts",62);
	TH1F * h1_e_pz0  = new TH1F("h1_e_pz0"  ,"h1_e_pz0"  ,100,   0, 10);       PrettyTH1F(h1_e_pz0   ,"electron p_{z} [GeV]" ,"Counts",62);
	TH1F * h1_e_p0   = new TH1F("h1_e_p0"   ,"h1_e_p0"   ,100,   0, 10);       PrettyTH1F(h1_e_p0    ,"electron |p| [GeV]"   ,"Counts",62);
	TH1F * h1_e_th0  = new TH1F("h1_e_th0"  ,"h1_e_th0"  ,100,   0, 30);       PrettyTH1F(h1_e_th0   ,"#theta_e [deg]"       ,"Counts",62);
	TH1F * h1_e_phi0 = new TH1F("h1_e_phi0" ,"h1_e_phi0" ,100,-190,190);       PrettyTH1F(h1_e_phi0  ,"#phi_e [deg]"         ,"Counts",62);

	TH1F * h1_e_px1  = new TH1F("h1_e_px1"  ,"h1_e_px1"  ,100,  -2,  2);       PrettyTH1F(h1_e_px1   ,"electron p_{x} [GeV]" ,"Counts", 2);
	TH1F * h1_e_py1  = new TH1F("h1_e_py1"  ,"h1_e_py1"  ,100,  -2,  2);       PrettyTH1F(h1_e_py1   ,"electron p_{y} [GeV]" ,"Counts", 2);
	TH1F * h1_e_pz1  = new TH1F("h1_e_pz1"  ,"h1_e_pz1"  ,100,   0, 10);       PrettyTH1F(h1_e_pz1   ,"electron p_{z} [GeV]" ,"Counts", 2);
	TH1F * h1_e_p1   = new TH1F("h1_e_p1"   ,"h1_e_p1"   ,100,   0, 10);       PrettyTH1F(h1_e_p1    ,"electron |p| [GeV]"   ,"Counts", 2);
	TH1F * h1_e_th1  = new TH1F("h1_e_th1"  ,"h1_e_th1"  ,100,   0, 30);       PrettyTH1F(h1_e_th1   ,"#theta_e [deg]"       ,"Counts", 2);
	TH1F * h1_e_phi1 = new TH1F("h1_e_phi1" ,"h1_e_phi1" ,100,-190,190);       PrettyTH1F(h1_e_phi1  ,"#phi_e [deg]"         ,"Counts", 2);

	h1_e_px0  -> SetMinimum(1e-02);
        h1_e_py0  -> SetMinimum(1e-02);
        h1_e_pz0  -> SetMinimum(1e-02);
        h1_e_p0   -> SetMinimum(1e-02);
        h1_e_th0  -> SetMinimum(1e-02);
        h1_e_phi0 -> SetMinimum(1e-02);

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	BEvent        event       ("REC::Event"       ,reader);
	BParticle     particles   ("REC::Particle"    ,reader);
	BCalorimeter  calo        ("REC::Calorimeter" ,reader);
	BScintillator scintillator("REC::Scintillator",reader);
	BBand         band_hits   ("BAND::hits"       ,reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;

		//band_hits.show();
		//event.show();
		//scintillator.show();	

		// Particle bank
		int pid0       = particles.getPid    (0);       // electron candidate id assigned by clas
		TVector3 V3_ev = particles.getV3v    (0);       // electron candidate vertex vector
		TVector3 V3_ep = particles.getV3P    (0);       // electron candidate momentum vector
		float chr0     = particles.getCharge (0);       // electron candidate charge
		float eBeta    = particles.getBeta   (0);       // electron candidate beta = v/c
		float chi2pid  = particles.getChi2pid(0);       // electron candidate goodness of pid fit
		int eStatus    = particles.getStatus (0);       // electron candidate status

		// Event bank
		double t_vtx   = event.getSTTime(0);

		// -------------------------------------------------------------------------
		// Only keep events for which the first particle is an electron
		if(             (pid0!=11              )||
				(chr0!=-1              )
		  ) continue;

		h1_e_px0  -> Fill( V3_ep.X    ()         );
		h1_e_py0  -> Fill( V3_ep.Y    ()         );
		h1_e_pz0  -> Fill( V3_ep.Z    ()         );
		h1_e_p0   -> Fill( V3_ep.Mag  ()         );
		h1_e_th0  -> Fill( V3_ep.Theta()*rad2deg );
		h1_e_phi0 -> Fill( V3_ep.Phi  ()*rad2deg );

		// Looking in BAND
		int nHits = band_hits.getSize();
		if( nHits > 1 ) continue;
		for(int hit = 0; hit < nHits; hit++) {

			int    sector            = band_hits.getSector      (hit);
			float adcLcorr           = band_hits.getAdcLcorr    (hit);
			float adcRcorr           = band_hits.getAdcRcorr    (hit);
			float  meantimeFadc      = band_hits.getMeantimeFadc(hit);

			if( adcLcorr < 4000 || adcRcorr < 4000 ) continue;
			if( sector     == 3 || sector     == 4 ) continue;

			double corr_meanT = meantimeFadc-t_vtx;

			// Fill histograms
			h1_ToF_fadc_0      -> Fill( corr_meanT );

			// Only select events in the photon peak
			if( (corr_meanT>43) && (corr_meanT<55) ){
				h1_ToF_fadc_1 -> Fill( corr_meanT );
				h1_e_px1      -> Fill( V3_ep.X    ()         );
                		h1_e_py1      -> Fill( V3_ep.Y    ()         );
                		h1_e_pz1      -> Fill( V3_ep.Z    ()         );
                		h1_e_p1       -> Fill( V3_ep.Mag  ()         );
                		h1_e_th1      -> Fill( V3_ep.Theta()*rad2deg );
                		h1_e_phi1     -> Fill( V3_ep.Phi  ()*rad2deg );
			}

		}// end loop over hits in event


	}// end file

	// ------------------------------------------------
	// Printing results in canvases
	TCanvas * c0 = new TCanvas("c0");
	h1_ToF_fadc_0 -> Draw(      );
	h1_ToF_fadc_1 -> Draw("same");
	c0 -> Modified();
	c0 -> Update();

	TCanvas * c1 = new TCanvas("c1");
	c1 -> Divide(2,2);
	c1 -> cd(1);	gPad -> SetLogy();	h1_e_px0  -> Draw();	h1_e_px1  -> Draw("same");
        c1 -> cd(2);	gPad -> SetLogy();	h1_e_py0  -> Draw();	h1_e_py1  -> Draw("same");
        c1 -> cd(3);	gPad -> SetLogy();	h1_e_pz0  -> Draw();	h1_e_pz1  -> Draw("same");
        c1 -> cd(4);	gPad -> SetLogy();	h1_e_p0   -> Draw();	h1_e_p1   -> Draw("same");
	c1 -> Modified();
        c1 -> Update();

	TCanvas * c2 = new TCanvas("c2");
	c2 -> Divide(2,2);
        c2 -> cd(1);	gPad -> SetLogy();	h1_e_th0  -> Draw();	h1_e_th1  -> Draw("same");
        c2 -> cd(4);	gPad -> SetLogy();	h1_e_phi0 -> Draw();	h1_e_phi1 -> Draw("same");
	c2 -> Modified();
        c2 -> Update();

	// ------------------------------------------------
	// Saving results to a pdf file
	c0 -> Print("results_mean_time.pdf");

	myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color) {
	h1 -> GetXaxis() -> SetTitle(titx);
	h1 -> GetYaxis() -> SetTitle(tity);
	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);
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
