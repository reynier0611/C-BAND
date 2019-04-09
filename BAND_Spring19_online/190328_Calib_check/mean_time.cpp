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
#include "TF1.h"

#include "reader.h"
#include "node.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

double p2p   [600] = {0};
double l2l   [600] = {0};
double twAL  [600] = {0};
double twBL  [600] = {0};
double twAR  [600] = {0};
double twBR  [600] = {0};
double LRtdc [600] = {0};
double LRfadc[600] = {0};

// Forward-declaring functions
void LoadPaddleCorrectionPar();
void LoadT_WalkCorrectionPar();
void Load_LminRCorrectionPar();
double TDC_corr( int barId );
double TimeWalk_corr( int barId , TString LR , double ADC );
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

	TString inputFile;

	if(argc==2) {
		inputFile = argv[1];
	}
	else {
		cout << "=========================\nRun this code as:\n./code path/to/input/file\n=========================" << endl;
		exit(0);
	}

	TString outRootName = "out";
	// ----------------------------------------------------------------------------------
        // Load time-walk correction parameters only if data is not already corrected
        cout << "Is your data already timewalk corrected?\n1) no, 2) yes, 3) I don't know?" << endl;
        int optionTW;
        cin >> optionTW;
        if(optionTW==3){
                cout << "Bitch please! Get your shit together and then come back. I cannot answer this for you." << endl;
                exit(0);
        }
        if(optionTW==1){
                LoadT_WalkCorrectionPar();
		outRootName + "_TW";
	}
	// ----------------------------------------------------------------------------------
        // Load L-R correction parameters only if data is not already corrected
	cout << "Is your data already L-R time corrected?\n1) no, 2) yes, 3) I don't know?" << endl;
	int optionLR;
	cin >> optionLR;
	if(optionLR==3){
                cout << "Bitch please! Get your shit together and then come back. I cannot answer this for you." << endl;
                exit(0);
        }
	if(optionLR==1){
		Load_LminRCorrectionPar();
		outRootName + "_LR";
	}
	// ----------------------------------------------------------------------------------
        // Load TDC paddle-to-paddle correction parameters only if data is not already corrected
        cout << "Is your data already TDC paddle offset corrected?\n1) no, 2) yes, 3) I don't know?" << endl;
        int option;
        cin >> option;
        if(option==3){
                cout << "Bitch please! Get your shit together and then come back. I cannot answer this for you." << endl;
                exit(0);
        }
        if(option==1){
                LoadPaddleCorrectionPar();
		outRootName + "_p2p";
	}
	// ---------------------------------------------------------------------------------- 

	outRootName + ".root";

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
	TH1F * h1_ToF_tdc_uncorr  = new TH1F("h1_ToF_tdc_uncorr" ,"uncorrected",1800 ,200,  600);	PrettyTH1F(h1_ToF_tdc_uncorr  ,"(L+R)/2 corrected - RF - pathlength (TDC)","Counts",4);
	TH1F * h1_ToF_tdc_corr    = new TH1F("h1_ToF_tdc_corr"   ,"corrected"  ,1800 ,200,  600);	PrettyTH1F(h1_ToF_tdc_corr    ,"(L+R)/2 corrected - RF - pathlength (TDC)","Counts",8);

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	double Ebeam = 10.6; //GeV
	TVector3 V3_Ebeam(0,0,Ebeam);
	TLorentzVector V4_Ebeam(V3_Ebeam,Ebeam);

	double mtar    = mD;
	TLorentzVector V4_mtar(0,0,0,mtar);

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

		// Calorimeter bank     
		float Epcal = calo.getPcalE(0);
		float Ee    = calo.getTotE (0);
		float lU    = calo.getLU   (0); // electron candidate distance on U-side [cm?]
		float lV    = calo.getLV   (0); // electron candidate distance on V-side [cm?]
		float lW    = calo.getLW   (0); // electron candidate distance on W-side [cm?]

		//if(Ee==0) continue;

		// Event bank
		double t_vtx   = event.getSTTime(0);

		// Scintillator bank
		double t_e     = scintillator.getTime(0);

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
				(chr0!=-1              )//||
				//(chi2pid>=cut_chi2pid  )||
				//(ep<=cut_ep            )||
				//(ep>=Ebeam             )||
				//(V3_ev.Z()>cut_max_vz  )||
				//(V3_ev.Z()<cut_min_vz  )||
				//(lU<cut_uvw            )||
				//(lV<cut_uvw            )||
				//(lW<cut_uvw            )||
				//(Epcal<cut_Epcal       )||
				//(TMath::Sqrt(W2)<=cut_W)||
				//(tof_e<cut_tof_e       )
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

			//if( adcLcorr < 4000 || adcRcorr < 4000 ) continue;
			if( TMath::Sqrt(adcLcorr*adcRcorr) < 6000 ) continue;
			//if( sector < 3 || sector > 4 ) continue; // Use this line to not include long bars
			if( sector == 3 || sector == 4 ) continue; // Use this line to not include short bars

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

			double meantimeTdc_byHand = (tTdcLcorr+tTdcRcorr)/2.;
                        //                                                 L-R corr         P2P & L2L corr                              time-walk corrections
			double meantimeTdcCorr = meantimeTdc_byHand - (LRtdc[barKey]/2.) + TDC_corr(barKey) + (TimeWalk_corr(barKey,"L",adcLcorr) + TimeWalk_corr(barKey,"R",adcRcorr))/2.;

			//if( sector == 3 || sector == 4 ) continue;

			double dL = TMath::Sqrt( pow(x,2) + pow(y,2) + pow(z,2) );

			// Fill histograms
			h1_ToF_tdc_uncorr      -> Fill(meantimeTdc     -t_vtx -(442.741-287.925) - dL/30.);
			h1_ToF_tdc_corr        -> Fill(meantimeTdcCorr -t_vtx                    - dL/30.);

		}// end loop over hits in event


	}// end file

	// Fitting the photo-peak
	TF1 * f_photopeak_1 = new TF1("f_photopeak_1","gaus(0) + pol2(3)",270,290);
	f_photopeak_1 -> SetParameter(0,100);
	f_photopeak_1 -> SetParameter(1,278);
	f_photopeak_1 -> SetParameter(2,  1);
	f_photopeak_1 -> SetParLimits(2,0,5);
	TF1 * f_photopeak_2 = new TF1("f_photopeak_2","gaus(0) + pol2(3)",270,290);
	f_photopeak_2 -> SetParameter(0,100);
        f_photopeak_2 -> SetParameter(1,278);
        f_photopeak_2 -> SetParameter(2,  1);
	f_photopeak_2 -> SetParLimits(2,0,5);

	h1_ToF_tdc_uncorr -> Fit("f_photopeak_1","R");
	h1_ToF_tdc_corr      -> Fit("f_photopeak_2","R");

	// Making plots
	TCanvas * c2 = new TCanvas("c2","c2",1200,800);
	c2 -> Divide(1,2);
	c2 -> cd(1);	gPad -> SetBottomMargin(0.2);	h1_ToF_tdc_uncorr -> Draw();
	c2 -> cd(2);	gPad -> SetBottomMargin(0.2);	h1_ToF_tdc_corr      -> Draw();
	c2 -> Modified();
	c2 -> Update();

	// Saving plots to a pdf file
	c2 -> Print("results_mean_time.pdf");

	// Saving plots to a root file
	TFile * output = new TFile(outRootName,"recreate");
	h1_ToF_tdc_uncorr -> Write("h1_ToF_tdc_uncorr");
	h1_ToF_tdc_corr      -> Write("h1_ToF_tdc_corr"     );
	output -> Close();

	cout << "ROOT FILE HAS BEEN CREATED" << endl;

	myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color) {
	h1 -> GetXaxis() -> SetTitle(titx);
	h1 -> GetYaxis() -> SetTitle(tity);
	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);

	h1 -> GetXaxis() -> CenterTitle();
        h1 -> GetYaxis() -> CenterTitle();

	h1 -> GetXaxis() -> SetTitleSize(0.07);
	h1 -> GetXaxis() -> SetLabelSize(0.07);
	h1 -> GetXaxis() -> SetNdivisions(107);

	h1 -> GetYaxis() -> SetTitleSize(0.07);
        h1 -> GetYaxis() -> SetLabelSize(0.07);
        h1 -> GetYaxis() -> SetNdivisions(107);
	h1 -> GetYaxis() -> SetTitleOffset(0.63);
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
// ========================================================================================================================================
void LoadPaddleCorrectionPar(){
        ifstream f;
        int layer, sector, component, barId;
        double parameter, temp;

        f.open("input_params/TDC_paddle_offsets.txt");
        while(!f.eof()){
		f >> sector;
                f >> layer;
                f >> component;
                barId = 100*sector + 10*layer + component;
                f >> parameter;
                f >> temp;
                p2p[barId] = parameter;
        }
        f.close();

        f.open("input_params/TDC_layer_offsets.txt");
        while(!f.eof()){
		f >> sector;
                f >> layer;
                f >> component;
                barId = 100*sector + 10*layer + component;
                f >> parameter;
                f >> temp;
                l2l[barId] = parameter;
        }
        f.close();
}
// ========================================================================================================================================
void LoadT_WalkCorrectionPar(){
	ifstream f;
        int layer, sector, component, barId;
        double parA, parB, temp;

	f.open("input_params/timeWalkPar_L.txt");
        while(!f.eof()){
                f >> sector;
                f >> layer;
                f >> component;
                barId = 100*sector + 10*layer + component;
                f >> parA;
		f >> parB;
                f >> temp;
		f >> temp;
                twAL[barId] = parA;
        	twBL[barId] = parB;
	}
        f.close();

        f.open("input_params/timeWalkPar_R.txt");
        while(!f.eof()){
                f >> sector;
                f >> layer;
                f >> component;
                barId = 100*sector + 10*layer + component;
		f >> parA;
                f >> parB;
                f >> temp;
                f >> temp;
		twAR[barId] = parA;
                twBR[barId] = parB;
        }
        f.close();
}
// ========================================================================================================================================
void Load_LminRCorrectionPar(){
	ifstream f;
        int layer, sector, component, barId;
        double par_TDC, par_FADC, temp;

	f.open("input_params/lr_offsets_261.txt");
        while(!f.eof()){
                f >> sector;
                f >> layer;
                f >> component;
                barId = 100*sector + 10*layer + component;
                f >> par_TDC;
                f >> par_FADC;
                f >> temp;
                f >> temp;
                LRtdc [barId] = TMath::Abs(par_TDC );
                LRfadc[barId] = TMath::Abs(par_FADC);
        }
        f.close();
}
// ========================================================================================================================================
double TDC_corr( int barId ){
        return -( p2p[barId] + l2l[barId] );
}
// ========================================================================================================================================
double TimeWalk_corr( int barId , TString LR , double ADC ){
	if(LR=="L"||LR=="l"||LR=="Left"||LR=="left"){
		return -( twAL[barId]/TMath::Sqrt(ADC) + twBL[barId] );
	}
	else if(LR=="R"||LR=="r"||LR=="Right"||LR=="right"){
                return -( twAR[barId]/TMath::Sqrt(ADC) + twBR[barId] );
        }
	return 0;
}
