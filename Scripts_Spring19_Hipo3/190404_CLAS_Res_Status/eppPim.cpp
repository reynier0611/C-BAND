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
	TH1F * h1_e_vz  = new TH1F("h1_e_vz"  ,"h1_e_vz"  ,100, -50, 50);	PrettyTH1F(h1_e_vz  ,"v_{z} [cm]"           ,"Counts",62);
	TH1F * h1_e_tof = new TH1F("h1_e_tof" ,"h1_e_tof" ,100,  20, 27);	PrettyTH1F(h1_e_tof ,"electron TOF [ns]"    ,"Counts",62);

	TH1F * h1_e_px  = new TH1F("h1_e_px"  ,"h1_e_px"  ,100,  -2,  2);	PrettyTH1F(h1_e_px  ,"electron p_{x} [GeV]" ,"Counts",62);
	TH1F * h1_e_py  = new TH1F("h1_e_py"  ,"h1_e_py"  ,100,  -2,  2);	PrettyTH1F(h1_e_py  ,"electron p_{y} [GeV]" ,"Counts",62);
	TH1F * h1_e_pz  = new TH1F("h1_e_pz"  ,"h1_e_pz"  ,100,   0, 10);	PrettyTH1F(h1_e_pz  ,"electron p_{z} [GeV]" ,"Counts",62);
	TH1F * h1_e_p   = new TH1F("h1_e_p"   ,"h1_e_p"   ,100,   0, 10);	PrettyTH1F(h1_e_p   ,"electron |p| [GeV]"   ,"Counts",62);

	TH1F * h1_e_th  = new TH1F("h1_e_th"  ,"h1_e_th"  ,100,   0, 30);	PrettyTH1F(h1_e_th  ,"#theta_e [deg]"       ,"Counts",62);
	TH1F * h1_e_phi = new TH1F("h1_e_phi" ,"h1_e_phi" ,100,-190,190);	PrettyTH1F(h1_e_phi ,"#phi_e [deg]"         ,"Counts",62);
	TH1F * h1_e_lu  = new TH1F("h1_e_lu"  ,"h1_e_lu"  ,100,   0,500);	PrettyTH1F(h1_e_lu  ,"distance on U-side"   ,"Counts",62);
	TH1F * h1_e_lv  = new TH1F("h1_e_lv"  ,"h1_e_lv"  ,100,   0,500);	PrettyTH1F(h1_e_lv  ,"distance on V-side"   ,"Counts",62);
	TH1F * h1_e_lw  = new TH1F("h1_e_lw"  ,"h1_e_lw"  ,100,   0,500);	PrettyTH1F(h1_e_lw  ,"distance on W-side"   ,"Counts",62);

	TH1F * h1_p_vz  = new TH1F("h1_p_vz"  ,"h1_p_vz"  ,100, -50, 50);	PrettyTH1F(h1_p_vz  ,"v_{z} [cm]"           ,"Counts",2);
	TH1F * h1_p_num = new TH1F("h1_p_num" ,"h1_p_num" , 20,   0, 10);	PrettyTH1F(h1_p_num ,"proton number"        ,"Counts",62);

	TH1F * h1_p_px  = new TH1F("h1_p_px"  ,"h1_p_px"  ,100,  -2,  2);	PrettyTH1F(h1_p_px  ,"proton p_{x} [GeV]"   ,"Counts",62);
	TH1F * h1_p_py  = new TH1F("h1_p_py"  ,"h1_p_py"  ,100,  -2,  2);	PrettyTH1F(h1_p_py  ,"proton p_{y} [GeV]"   ,"Counts",62);
	TH1F * h1_p_pz  = new TH1F("h1_p_pz"  ,"h1_p_pz"  ,100,   0,  9);	PrettyTH1F(h1_p_pz  ,"proton p_{z} [GeV]"   ,"Counts",62);
	TH1F * h1_p_p   = new TH1F("h1_p_p"   ,"h1_p_p"   ,100,   0,  9);	PrettyTH1F(h1_p_p   ,"proton |p| [GeV]"     ,"Counts",62);

	TH1F * h1_p_th  = new TH1F("h1_p_th"  ,"h1_p_th"  ,100,   0, 80);	PrettyTH1F(h1_p_th  ,"#theta_p [deg]"       ,"Counts",62);
	TH1F * h1_p_phi = new TH1F("h1_p_phi" ,"h1_p_phi" ,100,-190,190);	PrettyTH1F(h1_p_phi ,"#phi_p [deg]"         ,"Counts",62);

	TH1F * h1_pmiss = new TH1F("h1_pmiss" ,"P_{miss}" ,100,   0,  5);	PrettyTH1F(h1_pmiss ,"Pm [GeV]"             ,"Counts",62);
	TH1F * h1_pmx   = new TH1F("h1_pmx"   ,"h1_pmx"   ,100,  -2,  2);	PrettyTH1F(h1_pmx   ,"Pmx [GeV]"            ,"Counts",62);
	TH1F * h1_pmy   = new TH1F("h1_pmy"   ,"h1_pmy"   ,100,  -2,  2);	PrettyTH1F(h1_pmy   ,"Pmy [GeV]"            ,"Counts",62);
	TH1F * h1_pmz   = new TH1F("h1_pmz"   ,"h1_pmz"   ,100,  -3,  3);	PrettyTH1F(h1_pmz   ,"Pmz [GeV]"            ,"Counts",62);
	TH1F * h1_pm_th = new TH1F("h1_pm_th" ,"h1_pm_th" ,100,   0,180);	PrettyTH1F(h1_pm_th ,"#theta_{Pm} [deg]"    ,"Counts",62);
        TH1F * h1_pm_ph = new TH1F("h1_pm_ph" ,"h1_pm_ph" ,100,-190,190);	PrettyTH1F(h1_pm_ph ,"#phi_{Pm} [deg]"      ,"Counts",62);
	TH1F * h1_Mm    = new TH1F("h1_Mm"    ,"h1_Mm"    ,100,  -3,  3);	PrettyTH1F(h1_Mm    ,"m_{miss} [GeV]"       ,"Counts",62);
	TH1F * h1_MmSqr = new TH1F("h1_MmSqr" ,"h1_MmSqr" ,100,  -3,  3);	PrettyTH1F(h1_MmSqr ,"m_{miss}^{2} [GeV^{2}]","Counts",62);
	TH1F * h1_Em    = new TH1F("h1_Em"    ,"h1_Em"    ,100,  -3,  3);	PrettyTH1F(h1_Em    ,"Em [GeV]"             ,"Counts",62);

	TH1F * h1_W     = new TH1F("h1_W"     ,"h1_W"     ,100,   0,  4);	PrettyTH1F(h1_W     ,"W [GeV]"              ,"Counts",62);
	TH1F * h1_xB    = new TH1F("h1_xB"    ,"h1_xB"    ,100,   0,  4);	PrettyTH1F(h1_xB    ,"x_B"                  ,"Counts",62);

	TH1F * h1_dlt_vz_ep  = new TH1F("h1_dlt_vz_ep"  ,"electron - proton",100, -20, 20);	PrettyTH1F(h1_dlt_vz_ep  ,"#Delta v_{z} [cm]"    ,"Counts",2);
	TH1F * h1_dlt_vz_epip= new TH1F("h1_dlt_vz_epip","electron - #pi+"  ,100, -20, 20);	PrettyTH1F(h1_dlt_vz_epip,"#Delta v_{z} [cm]"    ,"Counts",3);
	TH1F * h1_dlt_vz_epim= new TH1F("h1_dlt_vz_epim","electron - #pi-"  ,100, -20, 20);	PrettyTH1F(h1_dlt_vz_epim,"#Delta v_{z} [cm]"    ,"Counts",62);

	// Goodness of PID cut
	TH1F * h1_chi2pid_e   = new TH1F("h1_chi2pid_e"  ,";elecron #chi^{2} PID;Counts" ,100,-30,30);
	TH1F * h1_chi2pid_p1  = new TH1F("h1_chi2pid_p1" ,";proton 1 #chi^{2} PID;Counts",100,-30,30);
	TH1F * h1_chi2pid_p2  = new TH1F("h1_chi2pid_p2" ,";proton 2 #chi^{2} PID;Counts",100,-30,30);
	TH1F * h1_chi2pid_pim = new TH1F("h1_chi2pid_pim",";#pi^{-} #chi^{2} PID;Counts" ,100,-30,30);

	// 2D histograms
	TH2F * h2_e_Ep_p_0   = new TH2F("h2_e_Ep_p_0"   ,"h2_e_Ep_p_0"  ,100,1   , 10,100,  0,0.4);
	TH2F * h2_e_Ep_p_1   = new TH2F("h2_e_Ep_p_1"   ,"h2_e_Ep_p_1"  ,100,1   , 10,100,  0,0.4);	
	TH2F * h2_e_th_phi   = new TH2F("h2_e_th_phi"   ,"h2_e_th_phi"  ,100,-190,190,100,  0, 30);
	TH2F * h2_p_th_phi   = new TH2F("h2_p_th_phi"   ,"h2_p_th_phi"  ,100,-190,190,100,  0, 80);
	TH2F * h2_beta_p_pos = new TH2F("h2_beta_p_pos" ,"h2_beta_p_pos",100,0   ,  6,100,0.1,1.1);
	TH2F * h2_beta_p_neg = new TH2F("h2_beta_p_neg" ,"h2_beta_p_neg",100,0   ,  6,100,0.1,1.1);
	TH2F * h2_beta_p_p   = new TH2F("h2_beta_p_p"   ,"h2_beta_p_p"  ,100,0   ,  6,100,0.1,1.1);
	TH2F * h2_beta_p_pip = new TH2F("h2_beta_p_pip" ,"h2_beta_p_pip",100,0   ,  6,100,0.1,1.1);
	TH2F * h2_beta_p_pim = new TH2F("h2_beta_p_pim" ,"h2_beta_p_pim",100,0   ,  6,100,0.1,1.1);
	TH2F * h2_e_vz_phi   = new TH2F("h2_e_vz_phi"   ,"h2_e_vz_phi"  ,100,-190,190,100,-50, 50);
	TH2F * h2_p_vz_phi   = new TH2F("h2_p_vz_phi"   ,"h2_p_vz_phi"  ,100,-190,190,100,-50, 50);
	TH2F * h2_e_tof_p    = new TH2F("h2_e_tof_p"    ,"h2_e_tof_p"   ,100,1   ,  7,100, 20, 27);
	TH2F * h2_p_dtT_p_0  = new TH2F("h2_p_dtT_p_0"  ,"h2_p_dtT_p_0" ,100,0   ,  6,100,- 4,  4);
	TH2F * h2_p_dtT_p_1  = new TH2F("h2_p_dtT_p_1"  ,"h2_p_dtT_p_1" ,100,0   ,  6,100,- 4,  4);
	TH2F * h2_p_tof_det  = new TH2F("h2_p_tof_det"  ,"h2_p_tof_det" , 72,0   , 36,100,  0, 40);
	TH2F * h2_p_dtT_det  = new TH2F("h2_p_dtT_det"  ,"h2_p_dtT_det" , 72,0   , 36,100,-20, 20);
	TH2F * h2_Em_Pm      = new TH2F("h2_Em_Pm"      ,"h2_Em_Pm"     ,100,0   ,  3,100,-2 ,  2);
	TH2F * h2_pe_pp      = new TH2F("h2_pe_pp"      ,"h2_pe_pp"     ,100,0   ,  6,100, 1 , 10);
	TH2F * h2_pe_the     = new TH2F("h2_pe_the"     ,"h2_pe_the"    ,100,0   , 30,100, 1 , 10);

	PrettyTH2F(h2_e_Ep_p_0  ,"p_{e} [GeV]"   ,"E_{e}/p_{e}"         );
	PrettyTH2F(h2_e_Ep_p_1  ,"p_{e} [GeV]"   ,"E_{e}/p_{e}"         );
	PrettyTH2F(h2_e_th_phi  ,"#phi_e [deg]"  ,"#theta_e [deg]"      );
	PrettyTH2F(h2_p_th_phi  ,"#phi_p [deg]"  ,"#theta_p [deg]"      );
	PrettyTH2F(h2_beta_p_pos,"p [GeV]"       ,"#beta"               );
	PrettyTH2F(h2_beta_p_neg,"p [GeV]"       ,"#beta"               );
	PrettyTH2F(h2_beta_p_p  ,"p [GeV]"       ,"#beta"               );
	PrettyTH2F(h2_beta_p_pip,"p [GeV]"       ,"#beta"               );
	PrettyTH2F(h2_beta_p_pim,"p [GeV]"       ,"#beta"               );
	PrettyTH2F(h2_e_vz_phi  ,"#phi_e [deg]"  ,"e v_{z} [cm]"        );
	PrettyTH2F(h2_p_vz_phi  ,"#phi_p [deg]"  ,"p v_{z} [cm]"        );
	PrettyTH2F(h2_e_tof_p   ,"p_{e} [GeV]"   ,"electron TOF [ns]"   );
	PrettyTH2F(h2_p_dtT_p_0 ,"p_{e} [GeV]"   ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_p_dtT_p_1 ,"p_{e} [GeV]"   ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_p_tof_det ,"detector id"   ,"candidate p tof [ns]");
	PrettyTH2F(h2_p_dtT_det ,"detector id"   ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_Em_Pm     ,"Pm [GeV]"      ,"Em [GeV]"            );
	PrettyTH2F(h2_pe_pp     ,"p p [GeV]"     ,"e p [GeV]"           );
	PrettyTH2F(h2_pe_the    ,"#theta_e [deg]","p_{e} [GeV]"         );
	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	BEvent        event       ("REC::Event"       ,reader);
	BParticle     particles   ("REC::Particle"    ,reader);
	BCalorimeter  calo        ("REC::Calorimeter" ,reader);
	BScintillator scintillator("REC::Scintillator",reader);

	int event_counter = 0;
	int ctr_triggers = 0;
	int ctr_electrons = 0;
	int ctr_selectedReaction = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;
		ctr_triggers++;

		//particles.show();
		//calo.show();
		//event.show();
		//bank_scintillator.show();	

		// Particle bank
		int pid0       = particles.getPid    (0);	// electron candidate id assigned by clas
		TVector3 V3_ev = particles.getV3v    (0);	// electron candidate vertex vector
		TVector3 V3_ep = particles.getV3P    (0);	// electron candidate momentum vector
		float chr0     = particles.getCharge (0);	// electron candidate charge
		float eBeta    = particles.getBeta   (0);	// electron candidate beta = v/c
		float chi2pid  = particles.getChi2pid(0);	// electron candidate goodness of pid fit
		int eStatus    = particles.getStatus (0);	// electron candidate status

		// Calorimeter bank	
		float Epcal = calo.getPcalE(0); 
		float Ee    = calo.getTotE (0);
		float lU    = calo.getLU   (0);	// electron candidate distance on U-side [cm?]
		float lV    = calo.getLV   (0);	// electron candidate distance on V-side [cm?]
		float lW    = calo.getLW   (0);	// electron candidate distance on W-side [cm?]

		if(Ee==0) continue;

		// Event bank
		double t_vtx   = event.getSTTime(0);

		// Scintillator bank
		double t_e     = scintillator.getTime(0);

		// calculated variables
		double ep     = V3_ep.Mag();		// electron candidate momentum magnitude [GeV]
		TLorentzVector V4_ep(V3_ep,ep);
		double tof_e  = t_e - t_vtx;		// electron candidate time-of-flight [ns]

		// Transfer variables
		TVector3 V3_q = V3_Ebeam - V3_ep;
		double Q2     = 4*ep*Ebeam*pow(TMath::Sin(V3_ep.Theta()/2.),2);       // Q-squared [GeV^2]
		double omega  = Ebeam - ep;                                              // Transfer energy [GeV]
		double W2     = mp*mp-Q2+2*omega*mp;
		double xB     = Q2/2./mp/omega; 

		TLorentzVector V4_q(V3_q,omega);

		// -------------------------------------------------------------------------
		// Fill some histograms before cutting on good electrons
		h2_e_Ep_p_0 -> Fill(ep,Ee/ep);
		h2_e_vz_phi -> Fill(rad2deg*V3_ep.Phi(),V3_ev.Z());

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

		ctr_electrons++;

		// -------------------------------------------------------------------------
		// Filling electron histograms
		h1_e_lu  -> Fill(lU                   );
		h1_e_lv  -> Fill(lV                   );
		h1_e_lw  -> Fill(lW                   );
		h1_e_px  -> Fill(V3_ep.X()            );
		h1_e_py  -> Fill(V3_ep.Y()            );
		h1_e_pz  -> Fill(V3_ep.Z()            );
		h1_e_p   -> Fill(V3_ep.Mag()          );
		h1_e_vz  -> Fill(V3_ev.Z()            );
		h1_W     -> Fill(TMath::Sqrt(W2)      );
		h1_xB    -> Fill(xB                   );
		h1_e_th  -> Fill(rad2deg*V3_ep.Theta());
		h1_e_phi -> Fill(rad2deg*V3_ep.Phi()  );
		h1_e_tof -> Fill(tof_e                );

		h2_e_Ep_p_1 -> Fill(ep           ,Ee/ep       );
		h2_e_th_phi -> Fill(rad2deg*V3_ep.Phi(),rad2deg*V3_ep.Theta());
		h2_e_tof_p  -> Fill(ep           ,tof_e       );
		h2_pe_the   -> Fill(rad2deg*V3_ep.Theta(),ep  );

		h1_chi2pid_e-> Fill(chi2pid);
		// -------------------------------------------------------------------------
		// Loop over the remaining particles and require a proton, a pi+, and a pi-

		// Proton variables
		int nProtons = 0;
		int tmp_fast_p1_idx  = 0;	// index of the proton 1
		int tmp_fast_p2_idx  = 0;        // index of the proton 2

		// pi- variables
		int nPim = 0;
		int tmp_fast_pim_idx  = 0;	// index of the fastest proton

		int nParticles = particles.getSize();	

		for(int par = 1; par < nParticles; par++) {
			// Particle bank

			int pidi       = particles.getPid    (par);       // id assigned by clas
			TVector3 V3_hv = particles.getV3v    (par);       // vertex vector
			TVector3 V3_hp = particles.getV3P    (par);       // momentum vector
			float chri     = particles.getCharge (par);       // charge
			float beta_p   = particles.getBeta   (par);       // beta = v/c
			double pp     = V3_hp.Mag();

			if     (chri== 1) h2_beta_p_pos -> Fill(V3_hp.Mag(), beta_p  );
			else if(chri==-1) h2_beta_p_neg -> Fill(V3_hp.Mag(), beta_p  );

			// -------------------------------------------------------------------------
			// Proton PID from Dan Carman
			// - (DONE)	pid=2212 from EB
			// -		bank 331 status word (not entirely sure what he means by this???)
			// -		p > 0.5 GeV (need to check if we want this for our event selection)
			// -		charge != 0 (I implemented specifically 
			// - (DONE)	p < Ebeam
			// -		Tof > 10 ns
			// -		abs(chi^2pid) < 5
			// -		Delta_t < 5 ns
			// -------------------------------------------------------------------------
			// Count number of protons, pi+, and pi-

			// Protons
			if(             (pidi==2212         )&&
					(chri==1            )&&
					(V3_hp.Mag()<Ebeam  )
			  ) {
				nProtons++;
				if(nProtons==1) tmp_fast_p1_idx=par;
				if(nProtons==2) tmp_fast_p2_idx=par;
			}
			// pi-
			if(             (pidi==-211         )&&
					(chri==-1           )&&
					(V3_hp.Mag()<Ebeam  )
			  ) {
				nPim++;	
				tmp_fast_pim_idx=par;
			}
			// ----------------------------------------------------------------------

		}

		if(nProtons==2&&nPim==1&&nParticles==4) {

			ctr_selectedReaction++;

			TVector3 V3_p1v  = particles.getV3v    (tmp_fast_p1_idx);		// proton vertex vector [cm]
			TVector3 V3_p1p  = particles.getV3P    (tmp_fast_p1_idx);		// proton momentum vector [GeV]
			float beta_p1    = particles.getBeta   (tmp_fast_p1_idx);		// proton beta = v/c
			float chi2pid_p1 = particles.getChi2pid(tmp_fast_p1_idx);

			TVector3 V3_p2v  = particles.getV3v    (tmp_fast_p2_idx);		// proton 2 vertex vector [cm]
			TVector3 V3_p2p  = particles.getV3P    (tmp_fast_p2_idx);		// proton 2 momentum vector [GeV]
			float beta_p2    = particles.getBeta   (tmp_fast_p2_idx);		// proton 2 beta = v/c
			float chi2pid_p2 = particles.getChi2pid(tmp_fast_p2_idx);

			TVector3 V3_pimv = particles.getV3v    (tmp_fast_pim_idx);		// pi- vertex vector [cm]
			TVector3 V3_pimp = particles.getV3P    (tmp_fast_pim_idx);		// pi- momentum vector [GeV]
			float beta_pim   = particles.getBeta   (tmp_fast_pim_idx);		// pi- beta = v/c
			float chi2pid_pim= particles.getChi2pid(tmp_fast_pim_idx);

			TLorentzVector V4_p1p, V4_p2p, V4_pimp;
			V4_p1p .SetXYZM( V3_p1p .X() , V3_p1p .Y() , V3_p1p .Z() , mp   );
			V4_p2p .SetXYZM( V3_p2p .X() , V3_p2p .Y() , V3_p2p .Z() , mp   );
			V4_pimp.SetXYZM( V3_pimp.X() , V3_pimp.Y() , V3_pimp.Z() , mPiC );

			double Ep1  = V4_p1p .E();
			double Ep2  = V4_p2p .E();
			double EPim = V4_pimp.E();

			// Missing momentum components
			TVector3       V3_Pm = V3_p1p + V3_p2p + V3_pimp - V3_q;
			TLorentzVector V4_Pm = V4_Ebeam + V4_mtar - V4_ep - V4_p1p - V4_p2p - V4_pimp;	

			double Mmiss  = V4_Pm.Mag();
			double Mmiss2 = V4_Pm.M2();	

			// Filling histograms
			h1_p_vz        -> Fill(V3_p1v.Z()              );
			h1_dlt_vz_ep   -> Fill(V3_p1v.Z()  - V3_ev.Z() );
			h1_dlt_vz_epip -> Fill(V3_p2v.Z()  - V3_ev.Z() );
			h1_dlt_vz_epim -> Fill(V3_pimv.Z() - V3_ev.Z() );
			h1_p_px        -> Fill(V3_p1p.X()              );
			h1_p_py        -> Fill(V3_p1p.Y()              );
			h1_p_pz        -> Fill(V3_p1p.Z()              );
			h1_p_p         -> Fill(V3_p1p.Mag()            );
			h1_p_th        -> Fill(rad2deg*V3_p1p.Theta()  );
			h1_p_phi       -> Fill(rad2deg*V3_p1p.Phi()    );
			h1_pmx         -> Fill(V3_Pm.X()               );
			h1_pmy         -> Fill(V3_Pm.Y()               );
			h1_pmz         -> Fill(V3_Pm.Z()               );
			h1_pm_th       -> Fill(rad2deg*V3_Pm.Theta()   );
                        h1_pm_ph       -> Fill(rad2deg*V3_Pm.Phi()     );
			h1_pmiss       -> Fill(V3_Pm.Mag()             );
			h1_Mm          -> Fill(Mmiss                   );
			h1_MmSqr       -> Fill(Mmiss2                  );
			//h1_Em       -> Fill(Emiss        );

			h2_p_th_phi   -> Fill(rad2deg*V3_p1p.Phi(), rad2deg*V3_p1p.Theta());
			h2_p_vz_phi   -> Fill(rad2deg*V3_p1p.Phi(), V3_p1v.Z()         );
			h2_beta_p_p   -> Fill(V3_p1p .Mag(), beta_p1      );
			h2_beta_p_pip -> Fill(V3_p2p .Mag(), beta_p2    );
			h2_beta_p_pim -> Fill(V3_pimp.Mag(), beta_pim    );
			//h2_p_dtT_p_1 -> Fill(pp          , delta_tP    );
			//h2_Em_Pm     -> Fill(Pm          , Emiss       );
			h2_pe_pp      -> Fill(V3_p1p .Mag(), ep          );

			h1_chi2pid_p1  -> Fill(chi2pid_p1 );
			h1_chi2pid_p2  -> Fill(chi2pid_p2 );
			h1_chi2pid_pim -> Fill(chi2pid_pim);

		}


		h1_p_num -> Fill((double)(nProtons));

	} // End of while loop over events

	// -------------------------------------------------------------------------------------------
	// Drawing histograms
	TCanvas * c0 = new TCanvas();
	c0 -> Divide(2,1);
	c0 -> cd(1);
	h1_e_vz -> Draw(      );
	h1_p_vz -> Draw("same");
	c0 -> cd(2);
	h1_dlt_vz_ep   -> Draw(      );
	h1_dlt_vz_epip -> Draw("same");
	h1_dlt_vz_epim -> Draw("same");
	c0 -> Modified();
	c0 -> Update();

	TCanvas * c1 = new TCanvas();
	c1 -> Divide(2, 2);
	c1 -> cd(1);	h1_pmx  -> Draw();
	c1 -> cd(2);	h1_pmy  -> Draw();
	c1 -> cd(3);	h1_pmz  -> Draw();
	c1 -> cd(4);	h1_pmiss-> Draw();
	c1 -> Modified();
	c1 -> Update(); 

	TCanvas * c2 = new TCanvas();
	c2 -> Divide(2, 2);
	c2 -> cd(1);	h1_e_th     -> Draw(      );
	c2 -> cd(2);	h2_e_th_phi -> Draw("COLZ");
	c2 -> cd(4);	h1_e_phi    -> Draw(      );
	c2 -> Modified();
	c2 -> Update();

	TCanvas * c3 = new TCanvas();
	c3 -> Divide(2, 1);
	c3 -> cd(1);	h1_W  -> Draw();
	c3 -> cd(2);	h1_xB -> Draw();
	c3 -> Modified();
	c3 -> Update();

	TCanvas * c5 = new TCanvas();
	c5 -> Divide(2,1);
	c5 -> cd(1);	gPad -> SetLogz();	h2_e_Ep_p_0 -> Draw("COLZ");
	c5 -> cd(2);	gPad -> SetLogz();	h2_e_Ep_p_1 -> Draw("COLZ");
	c5 -> Modified();
	c5 -> Update();

	TCanvas * c6 = new TCanvas();
	c6 -> Divide(2, 2);
	c6 -> cd(1);	h1_e_px -> Draw();
	c6 -> cd(2);	h1_e_py -> Draw();
	c6 -> cd(3);	h1_e_pz -> Draw();
	c6 -> cd(4);	h1_e_p  -> Draw();
	c6 -> Modified();
	c6 -> Update();

	TCanvas * c7 = new TCanvas();
	c7 -> Divide(2, 2);
	c7 -> cd(1);	h1_p_px -> Draw();
	c7 -> cd(2);	h1_p_py -> Draw();
	c7 -> cd(3);	h1_p_pz -> Draw();
	c7 -> cd(4);	h1_p_p  -> Draw();
	c7 -> Modified();
	c7 -> Update();

	TCanvas * c8 = new TCanvas();
	c8->Divide(3, 2);
	c8->cd(1);	gPad -> SetLogz();	h2_beta_p_pos -> Draw("COLZ");
	c8->cd(2);	gPad -> SetLogz();	h2_beta_p_p   -> Draw("COLZ");
	c8->cd(3);	gPad -> SetLogz();	h2_beta_p_pip -> Draw("COLZ");
	c8->cd(4);	gPad -> SetLogz();	h2_beta_p_neg -> Draw("COLZ");
	c8->cd(5);	gPad -> SetLogz();	h2_beta_p_pim -> Draw("COLZ");
	c8 -> Modified();
	c8 -> Update(); 

	TCanvas * c9 = new TCanvas();
	c9 -> Divide(2, 2);
	c9 -> cd(1);	h1_p_th    -> Draw(      );
	c9 -> cd(2);	h2_p_th_phi-> Draw("COLZ");
	c9 -> cd(4);	h1_p_phi   -> Draw(      );
	c9 -> Modified();
	c9 -> Update();

	TCanvas * c10 = new TCanvas();
	c10 -> Divide(2, 2);
	c10 -> cd(1);	h1_e_lu -> Draw();
	c10 -> cd(2);	h1_e_lv -> Draw();
	c10 -> cd(3);	h1_e_lw -> Draw();
	c10 -> Modified();
	c10 -> Update();

	TCanvas * c11 = new TCanvas();
	h1_p_num -> Draw();
	c11 -> Modified();
	c11 -> Update();
	
	TCanvas * c12 = new TCanvas();
        c12 -> Divide(2,1);
        c12 -> cd(1);
        h1_Mm -> Draw();
        c12 -> cd(2);
        h1_MmSqr -> Draw();
        c12 -> Modified();
        c12 -> Update();

	TCanvas * c13 = new TCanvas();
	c13 -> Divide(2,1);
	c13 -> cd(1);	h2_e_vz_phi-> Draw("COLZ");
	c13 -> cd(2);	h2_p_vz_phi-> Draw("COLZ");
	c13 -> Modified();
	c13 -> Update();

	TCanvas * c14 = new TCanvas();
	c14 -> Divide(2,1);
	c14 -> cd(1);	h1_e_tof  -> Draw(      );
	c14 -> cd(2);	h2_e_tof_p-> Draw("COLZ");
	c14 -> Modified();
	c14 -> Update();

	TCanvas * c15 = new TCanvas();
	c15 -> Divide(2, 1);
	c15 -> cd(1);	gPad -> SetLogz();	h2_p_dtT_p_0 -> Draw("COLZ");
	c15 -> cd(2);	gPad -> SetLogz();	h2_p_dtT_p_1 -> Draw("COLZ");
	c15 -> Modified();
	c15 -> Update();

	TCanvas * c16 = new TCanvas();
	h2_p_tof_det -> Draw("COLZ");
	c16 -> Modified();
	c16 -> Update();

	TCanvas * c17 = new TCanvas();
	h2_p_dtT_det -> Draw("COLZ");
	c17 -> Modified();
	c17 -> Update();

	TCanvas * c18 = new TCanvas();
	h1_Em -> Draw();
	c18 -> Modified();
	c18 -> Update();

	TCanvas * c19 = new TCanvas();
	h2_Em_Pm -> Draw("COLZ");
	c19 -> Modified();
	c19 -> Update();

	TCanvas * c20 = new TCanvas();
	h2_pe_pp -> Draw("COLZ");
	c20 -> Modified();
	c20 -> Update();

	TCanvas * c21 = new TCanvas();
        h1_pm_th -> Draw();
        c21 -> Modified();
        c21 -> Update();

        TCanvas * c22 = new TCanvas();
        h1_pm_ph -> Draw();
        c22 -> Modified();
        c22 -> Update();

	TCanvas * c23 = new TCanvas();
        h2_pe_the -> Draw("COLZ");
        c23 -> Modified();
        c23 -> Update();

	TCanvas * c24 = new TCanvas();
	c24 -> Divide(2,2);
	c24 -> cd(1);	h1_chi2pid_e   -> Draw();
	c24 -> cd(2);	h1_chi2pid_p1  -> Draw();
	c24 -> cd(3);	h1_chi2pid_p2  -> Draw();
	c24 -> cd(4);	h1_chi2pid_pim -> Draw();
	c24 -> Modified();
        c24 -> Update();

	// -------------------------------------------------------------------------------------------
	// Printing some info to screen
	cout << "Triggers: " << ctr_triggers << endl;
        cout << "Electrons: " << ctr_electrons << endl;
        cout << "Selected reaction: " << ctr_selectedReaction << endl;

	// -------------------------------------------------------------------------------------------
	// Saving plots to system

	c0  -> Print("results_eppPi.pdf(");
	c1  -> Print("results_eppPi.pdf" );
	c2  -> Print("results_eppPi.pdf" );
	c3  -> Print("results_eppPi.pdf" );
	c5  -> Print("results_eppPi.pdf" );
	c6  -> Print("results_eppPi.pdf" );
	c7  -> Print("results_eppPi.pdf" );
	c8  -> Print("results_eppPi.pdf" );
	c9  -> Print("results_eppPi.pdf" );
	c10 -> Print("results_eppPi.pdf" );
	c11 -> Print("results_eppPi.pdf" );
	c12 -> Print("results_eppPi.pdf" );
	c13 -> Print("results_eppPi.pdf" );
	c14 -> Print("results_eppPi.pdf" );
	c15 -> Print("results_eppPi.pdf" );
	c16 -> Print("results_eppPi.pdf" );
	c17 -> Print("results_eppPi.pdf" );
	c18 -> Print("results_eppPi.pdf" );
	c19 -> Print("results_eppPi.pdf" );
	c20 -> Print("results_eppPi.pdf" );
	c21 -> Print("results_eppPi.pdf" );
	c22 -> Print("results_eppPi.pdf" );
	c23 -> Print("results_eppPi.pdf" );
	c24 -> Print("results_eppPi.pdf)");

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
