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
#include "TFile.h"
#include "TTree.h"
#include "TVectorT.h"

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
void PrettyTH1(TH1F * h1,int color);
void PrettyTH2(TH2F * h2);
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
	TString inputFile, Reac, ReacPDF, par_names[4];
	double Ebeam, mtar;
	int nFSPart = 0;

	if(argc==4){
		// First argument
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
		// Second argument
		if(atoi(argv[2])==1){
			Reac = "(e,e'p )";	ReacPDF = "ep";
			nFSPart = 2;
			par_names[0] = "proton";
		}
		else if(atoi(argv[2])==2){
			Reac = "(e,e'π˖)";	ReacPDF = "ePip";
			nFSPart = 2;
			par_names[0] = "#pi^{+}";
		}
		else if(atoi(argv[2])==3){
			Reac = "(e,e'pπ˖π˗)";	ReacPDF = "epPipPim";
			nFSPart = 4;
			par_names[0] = "proton";
			par_names[1] = "#pi^{+}";
			par_names[2] = "#pi^{-}";
		}
		else if(atoi(argv[2])==4){	
			Reac = "e,e'pp π˗";	ReacPDF = "eppPim";
			nFSPart = 4;
			par_names[0] = "proton 1";
			par_names[1] = "proton 2";
			par_names[2] = "#pi^{-}";
		}
		cout << "Will be looking for the reaction: " << Reac << endl;
		// Third argument
		inputFile = argv[3];
	}
	else {
		cout << "=========================\nRun this code as:\n./code A B path/to/input/file\n" << endl;
		cout << "where: A = 1 -> Ebeam =  6.4 GeV, target = H (i.e. RGA)" << endl;
		cout << "         = 2 -> Ebeam = 10.6 GeV, target = D (i.e. RGB)" << endl << endl;
		cout << "       B = 1 -> (e,e'p )" << endl;
		cout << "         = 2 -> (e,e'π˖)" << endl;
		cout << "         = 3 -> (e,e'pπ˖π˗)" << endl;
		cout << "         = 4 -> (e,e'pp π˗) -> (study CLAS resolution)" << endl;
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
	// Electron histograms
	// 1D histograms
	TH1F * h1_e_vz      = new TH1F("h1_e_vz"     ,"v_{z} [cm];Counts"              ,100, -50, 50);	PrettyTH1(h1_e_vz     ,4);
	TH1F * h1_e_tof     = new TH1F("h1_e_tof"    ,"electron TOF [ns];Counts"       ,100,  20, 27);	PrettyTH1(h1_e_tof    ,4);
	TH1F * h1_e_px      = new TH1F("h1_e_px"     ,"electron p_{x} [GeV];Counts"    ,100,  -2,  2);	PrettyTH1(h1_e_px     ,4);
	TH1F * h1_e_py      = new TH1F("h1_e_py"     ,"electron p_{y} [GeV];Counts"    ,100,  -2,  2);	PrettyTH1(h1_e_py     ,4);
	TH1F * h1_e_pz      = new TH1F("h1_e_pz"     ,"electron p_{z} [GeV];Counts"    ,100,   0, 10);	PrettyTH1(h1_e_pz     ,4);
	TH1F * h1_e_p       = new TH1F("h1_e_p"      ,"electron |p| [GeV];Counts"      ,100,   0, 10);	PrettyTH1(h1_e_p      ,4);
	TH1F * h1_e_th      = new TH1F("h1_e_th"     ,"#theta_{e} [deg];Counts"        ,100,   0, 30);	PrettyTH1(h1_e_th     ,4);
	TH1F * h1_e_phi     = new TH1F("h1_e_phi"    ,"#phi_{e} [deg];Counts"          ,100,-190,190);	PrettyTH1(h1_e_phi    ,4);
	TH1F * h1_e_lu      = new TH1F("h1_e_lu"     ,"distance on U-side;Counts"      ,100,   0,500);	PrettyTH1(h1_e_lu     ,4);
	TH1F * h1_e_lv      = new TH1F("h1_e_lv"     ,"distance on V-side;Counts"      ,100,   0,500);	PrettyTH1(h1_e_lv     ,4);
	TH1F * h1_e_lw      = new TH1F("h1_e_lw"     ,"distance on W-side;Counts"      ,100,   0,500);	PrettyTH1(h1_e_lw     ,4);
	
	// 2D histograms
	TH2F * h2_e_Ep_p_0  = new TH2F("h2_e_Ep_p_0" ,"before cuts;p_{e} [GeV];E_{e}/p_{e}"         ,100,1   , 10,100,  0,0.4 );	PrettyTH2(h2_e_Ep_p_0);
        TH2F * h2_e_Ep_p_1  = new TH2F("h2_e_Ep_p_1" ,"after cuts;p_{e} [GeV];E_{e}/p_{e}"          ,100,1   , 10,100,  0,0.4 );	PrettyTH2(h2_e_Ep_p_1);
	TH2F * h2_e_th_phi  = new TH2F("h2_e_th_phi" ,";electron #phi [deg];electron #theta [deg]"  ,100,-190,190,100,  0, 30 );	PrettyTH2(h2_e_th_phi);
	TH2F * h2_e_vz_phi  = new TH2F("h2_e_vz_phi" ,";electron #phi [deg];electron v_{z} [cm]"    ,100,-190,190,100,-50, 50 );	PrettyTH2(h2_e_vz_phi);
	TH2F * h2_e_tof_p   = new TH2F("h2_e_tof_p"  ,";p_{e} [GeV];electron TOF [ns]"              ,100,1   ,  7,100, 20, 27 );	PrettyTH2(h2_e_tof_p );
	TH2F * h2_pe_the    = new TH2F("h2_pe_the"   ,";#theta_e [deg];p_{e} [GeV]"                 ,100,0   , 30,100, 1 , 10 );	PrettyTH2(h2_pe_the  );

	// Missing histograms
	// 1D histograms
	TH1F * h1_pmiss     = new TH1F("h1_pmiss"    ,";p_{miss} [GeV];Counts"         ,100,   0,  5);	PrettyTH1(h1_pmiss    ,4);
        TH1F * h1_pmx       = new TH1F("h1_pmx"      ,";p_{miss,x} [GeV];Counts"       ,100,  -2,  2);	PrettyTH1(h1_pmx      ,4);
        TH1F * h1_pmy       = new TH1F("h1_pmy"      ,";p_{miss,y} [GeV];Counts"       ,100,  -2,  2);	PrettyTH1(h1_pmy      ,4);
        TH1F * h1_pmz       = new TH1F("h1_pmz"      ,";p_{miss,z} [GeV];Counts"       ,100,  -3,  3);	PrettyTH1(h1_pmz      ,4);
	TH1F * h1_Mm        = new TH1F("h1_Mm"       ,";m_{miss} [GeV];Counts"         ,100,  -3,  3);	PrettyTH1(h1_Mm       ,4);
	TH1F * h1_MmSqr     = new TH1F("h1_MmSqr"    ,";m_{miss}^{2} [GeV^{2}];Counts" ,100,  -3,  3);	PrettyTH1(h1_MmSqr    ,4);
	TH1F * h1_pm_th     = new TH1F("h1_pm_th"    ,";#theta_{Pm} [deg];Counts"      ,100,   0,180);	PrettyTH1(h1_pm_th    ,4);
        TH1F * h1_pm_ph     = new TH1F("h1_pm_ph"    ,";#phi_{Pm} [deg];Counts"        ,100,-190,190);	PrettyTH1(h1_pm_ph    ,4);
        TH1F * h1_Em        = new TH1F("h1_Em"       ,";E_{miss} [GeV];Counts"         ,100,  -3,  3);	PrettyTH1(h1_Em       ,4);
	TH1F * h1_W         = new TH1F("h1_W"        ,";W [GeV];Counts"                ,100,   0,  4);	PrettyTH1(h1_W        ,4);
	TH1F * h1_xB        = new TH1F("h1_xB"       ,";x_{B};Counts"                  ,100,   0,  4);	PrettyTH1(h1_xB       ,4);
	TH1F * h1_chi2pid_e = new TH1F("h1_chi2pid_e",";elecron #chi^{2} PID;Counts"   ,100, -30, 30);	PrettyTH1(h1_chi2pid_e,4);

	// 2D histograms
	TH2F * h2_Em_Pm     = new TH2F("h2_Em_Pm"    ,";p_{miss} [GeV];E_{miss} [GeV]" ,100,0 , 3,100,-2 , 2);	PrettyTH2(h2_Em_Pm);	

	// Histograms for particles other than the electron
	// 1D histograms	
	TH1F ** h1_chi2pid_i  = new TH1F*[nFSPart-1];
	TH1F ** h1_i_vz       = new TH1F*[nFSPart-1];
        TH1F ** h1_i_dlt_vz   = new TH1F*[nFSPart-1];
	TH1F ** h1_i_px       = new TH1F*[nFSPart-1];
        TH1F ** h1_i_py       = new TH1F*[nFSPart-1];
        TH1F ** h1_i_pz       = new TH1F*[nFSPart-1];
        TH1F ** h1_i_p        = new TH1F*[nFSPart-1];
        TH1F ** h1_i_th       = new TH1F*[nFSPart-1];
        TH1F ** h1_i_phi      = new TH1F*[nFSPart-1];

	// 2D histograms
	TH2F ** h2_i_th_phi   = new TH2F*[nFSPart-1];
        TH2F ** h2_beta_i_p   = new TH2F*[nFSPart-1];
        TH2F ** h2_i_vz_phi   = new TH2F*[nFSPart-1];
	TH2F ** h2_pe_pi      = new TH2F*[nFSPart-1];

	for(int i = 0 ; i < nFSPart-1 ; i++){
		h1_chi2pid_i[i] = new TH1F(Form("h1_chi2pid_%i",i),";"+par_names[i]+" #chi^{2} PID;Counts"    ,100, -30, 30); 	PrettyTH1(h1_chi2pid_i[i],i);
		h1_i_vz     [i] = new TH1F(Form("h1_%i_vz"     ,i),";"+par_names[i]+" v_{z} [cm];Counts"      ,100, -50, 50);	PrettyTH1(h1_i_vz     [i],i);
	        h1_i_dlt_vz [i] = new TH1F(Form("h1_%i_dlt_vz" ,i),";"+par_names[i]+"#Delta v_{z} [cm];Counts",100, -20, 20);	PrettyTH1(h1_i_dlt_vz [i],i);
		h1_i_px     [i] = new TH1F(Form("h1_%i_px"     ,i),";"+par_names[i]+" p_{x} [GeV];Counts"     ,100,  -2,  2);	PrettyTH1(h1_i_px     [i],i);
	        h1_i_py     [i] = new TH1F(Form("h1_%i_py"     ,i),";"+par_names[i]+" p_{y} [GeV];Counts"     ,100,  -2,  2);	PrettyTH1(h1_i_py     [i],i);
	        h1_i_pz     [i] = new TH1F(Form("h1_%i_pz"     ,i),";"+par_names[i]+" p_{z} [GeV];Counts"     ,100,   0,  9);	PrettyTH1(h1_i_pz     [i],i);
	        h1_i_p      [i] = new TH1F(Form("h1_%i_p"      ,i),";"+par_names[i]+" |p| [GeV];Counts"       ,100,   0,  9);	PrettyTH1(h1_i_p      [i],i);
        	h1_i_th     [i] = new TH1F(Form("h1_%i_th"     ,i),";"+par_names[i]+" #theta [deg];Counts"    ,100,   0, 80);	PrettyTH1(h1_i_th     [i],i);
        	h1_i_phi    [i] = new TH1F(Form("h1_%i_phi"    ,i),";"+par_names[i]+" #phi [deg];Counts"      ,100,-190,190);	PrettyTH1(h1_i_phi    [i],i);

		// 2D histograms
		h2_i_th_phi [i] = new TH2F(Form("h2_%i_th_phi" ,i),";"+par_names[i]+" #phi [deg];"+par_names[i]+" #theta [deg]" ,100,-190,190,100,  0, 80);	
		h2_beta_i_p [i] = new TH2F(Form("h2_beta_%i_p" ,i),";"+par_names[i]+" p [GeV];"   +par_names[i]+" #beta"        ,100,0   ,  6,100,0.1,1.1);
        	h2_i_vz_phi [i] = new TH2F(Form("h2_%i_vz_phi" ,i),";"+par_names[i]+" #phi [deg];"+par_names[i]+" v_{z} [cm]"   ,100,-190,190,100,-50, 50);	
		h2_pe_pi    [i] = new TH2F(Form("h2_pe_p%i"    ,i),";"+par_names[i]+" p [GeV];"   +par_names[i]+"E [GeV]"       ,100,0   ,  6,100, 1 , 10);

		PrettyTH2(h2_i_th_phi[i]);
		PrettyTH2(h2_beta_i_p[i]);
		PrettyTH2(h2_i_vz_phi[i]);
		PrettyTH2(h2_pe_pi   [i]);
	}

	// Other histograms
	// 2D histograms
	TH2F * h2_beta_p_pos = new TH2F("h2_beta_p_pos" ,"p (positive) [GeV];#beta (positive)",100,0 , 6,100,0.1,1.1);	PrettyTH2(h2_beta_p_pos);
	TH2F * h2_beta_p_neg = new TH2F("h2_beta_p_neg" ,"p (negative) [GeV];#beta (negative)",100,0 , 6,100,0.1,1.1);	PrettyTH2(h2_beta_p_neg);

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
        // Setting up output root file
        TFile * out = new TFile("outRoot_"+ReacPDF+".root","RECREATE");
        TTree * tree = new TTree("T","(e,"+ReacPDF+")");

        double br_pex, br_pey, br_pez, br_e_chi2pid;
	double br_vex, br_vey, br_vez;
	double br_qx , br_qy , br_qz ;
        double br_Q2 , br_Nu , br_W2 , br_xB;
	tree -> Branch("br_e_chi2pid", &br_e_chi2pid , "br_e_chi2pid/D");
	tree -> Branch("br_pex"      , &br_pex       , "br_pex/D"      );
	tree -> Branch("br_pey"      , &br_pey       , "br_pey/D"      );
	tree -> Branch("br_pez"      , &br_pez       , "br_pez/D"      );
	tree -> Branch("br_vex"      , &br_vex       , "br_vex/D"      );
        tree -> Branch("br_vey"      , &br_vey       , "br_vey/D"      );
        tree -> Branch("br_vez"      , &br_vez       , "br_vez/D"      );
	tree -> Branch("br_qx"       , &br_qx        , "br_qx/D"       );
	tree -> Branch("br_qy"       , &br_qy        , "br_qy/D"       );
	tree -> Branch("br_qz"       , &br_qz        , "br_qz/D"       );
	tree -> Branch("br_Q2"       , &br_Q2        , "br_Q2/D"       );
	tree -> Branch("br_Nu"       , &br_Nu        , "br_Nu/D"       );
	tree -> Branch("br_W2"       , &br_W2        , "br_W2/D"       );
	tree -> Branch("br_xB"       , &br_xB        , "br_xB/D"       );

	double br_vix[nFSPart], br_viy[nFSPart], br_viz[nFSPart];
	double br_pix[nFSPart], br_piy[nFSPart], br_piz[nFSPart], br_i_chi2pid[nFSPart];
	for(int par = 0 ; par < nFSPart-1 ; par++){
		tree -> Branch(Form("br_%i_chi2pid",par), &br_i_chi2pid[par] , Form("br_%i_chi2pid/D",par));
		tree -> Branch(Form("br_p%ix"      ,par), &br_pix      [par] , Form("br_p%ix/D"      ,par));
		tree -> Branch(Form("br_p%iy"      ,par), &br_piy      [par] , Form("br_p%iy/D"      ,par));
		tree -> Branch(Form("br_p%iz"      ,par), &br_piz      [par] , Form("br_p%iz/D"      ,par));
		tree -> Branch(Form("br_v%ix"      ,par), &br_vix      [par] , Form("br_v%ix/D"      ,par));
                tree -> Branch(Form("br_v%iy"      ,par), &br_viy      [par] , Form("br_v%iy/D"      ,par));
                tree -> Branch(Form("br_v%iz"      ,par), &br_viz      [par] , Form("br_v%iz/D"      ,par));
	}

	double br_thpm, br_phipm;
	double br_Pmx , br_Pmy  , br_Pmz, br_Mm, br_Mm2;
	tree -> Branch("br_Pmx"  , &br_Pmx   , "br_Pmx/D"  );
	tree -> Branch("br_Pmy"  , &br_Pmy   , "br_Pmy/D"  );
	tree -> Branch("br_Pmz"  , &br_Pmz   , "br_Pmz/D"  );
	tree -> Branch("br_Mm"   , &br_Mm    , "br_Mm/D"   );
	tree -> Branch("br_Mm2"  , &br_Mm2   , "br_Mm2/D"  );
	tree -> Branch("br_thpm" , &br_thpm  , "br_thpm/D" );
	tree -> Branch("br_phipm", &br_phipm , "br_phipm/D");

	TVectorT<double> part_mass(nFSPart-1);
	TVectorT<double> part_char(nFSPart-1);
	TVectorT<double> part_pid (nFSPart-1);

	if(ReacPDF=="ePip"){
		part_mass[0] = mPiC;
		part_char[0] = 1   ;
		part_pid [0] = 211 ;
	}
        if( (ReacPDF=="ep")||(ReacPDF=="epPipPim")||(ReacPDF=="eppPim") ){
       		part_mass[0] = mp  ; 
        	part_char[0] = 1   ;
		part_pid [0] = 2212;
	}
        if(ReacPDF=="epPipPim"){
    		part_mass[1] = mPiC;
		part_char[1] = 1   ;
		part_pid [1] = 211 ;
		part_mass[2] = mPiC;       
        	part_char[2] = -1  ;
		part_pid [2] = -211;
	}
        if(ReacPDF=="eppPim"){
		part_mass[1] = mp  ;
                part_char[1] = 1   ;
		part_pid [1] = 2212;
		part_mass[2] = mPiC;
		part_char[2] = -1  ;
		part_pid [2] = -211;
	}

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
		// Only keep events for which the first particle is an electron
		if(             (pid0!=11                   )||
				(chr0!=-1                   )||
				(abs(chi2pid)>=cut_chi2pid  )||
				(ep<=cut_ep                 )||
				(ep>=Ebeam                  )||
				(V3_ev.Z()>cut_max_vz       )||
				(V3_ev.Z()<cut_min_vz       )||
				(lU<cut_uvw                 )||
				(lV<cut_uvw                 )||
				(lW<cut_uvw                 )||
				(Epcal<cut_Epcal            )||
				(TMath::Sqrt(W2)<=cut_W     )||
				(tof_e<cut_tof_e            )
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
		// Loop over the remaining particles in the bank
		// Proton variables
		int nProtons = 0;
		int tmp_fast_p1_idx  = 0;	// index of the proton 1
		int tmp_fast_p2_idx  = 0;	// index of the proton 2
		// π˗ variables
		int nPim = 0;
		int tmp_fast_pim_idx  = 0;
		// π˖ variables
		int nPip = 0;
		int tmp_fast_pip_idx = 0;

		const int nParticles = particles.getSize();	

		for(int par = 1; par < nParticles; par++) {
			int pidi       = particles.getPid    (par);	// PID assigned by CLAS
			float chi2pidi = particles.getChi2pid(par);     // χ² PID
			float chri     = particles.getCharge (par);	// charge
			TVector3 V3_pi = particles.getV3P    (par);	// momentum vector
			float betai    = particles.getBeta   (par);     // β = v/c

			if     (chri== 1) h2_beta_p_pos -> Fill(V3_pi.Mag(), betai );
			else if(chri==-1) h2_beta_p_neg -> Fill(V3_pi.Mag(), betai );

			// -------------------------------------------------------------------------
			// PID for particles other than electrons
			// Count number of protons, π˖, and π˗
			if((abs(chi2pidi)<cut_chi2pid)&&(V3_pi.Mag()<Ebeam)){
				// -------
				// Protons
				if(pidi==2212){
					nProtons++;
					if(nProtons==1) tmp_fast_p1_idx=par;
					if(nProtons==2) tmp_fast_p2_idx=par;
				}
				// π˗
				if(pidi==-211) {
					nPim++;	
					tmp_fast_pim_idx=par;
				}
				// π˖
				if(pidi==211) {
					nPip++; 
					tmp_fast_pip_idx=par;
				}
				// ----------------------------------------------------------------------
			}
		}

		if(      	((ReacPDF=="ep"      )&&(nProtons==1)                      &&(nParticles==2))||
				((ReacPDF=="ePip"    )&&(nPip    ==1)                      &&(nParticles==2))||
				((ReacPDF=="epPipPim")&&(nProtons==1)&&(nPim==1)&&(nPip==1)&&(nParticles==4))||
				((ReacPDF=="eppPim"  )&&(nProtons==2)&&(nPim==1)           &&(nParticles==4))
		  ){
			ctr_selectedReaction++;

			int par_idx[4] = {0};
			TVector3 V3_par_p[3];
			TVector3 V3_par_v[3];
			TLorentzVector V4_par_p[3];
			float beta_par   [3];
			float chi2pid_par[3];
			double m_par     [3];
			double E_par     [3];

			if(ReacPDF=="ePip"){
				par_idx[0] = tmp_fast_pip_idx; // first particle (other than the e') is a π˖
				m_par  [0] = mPiC;
			}

			if( (ReacPDF=="ep")||(ReacPDF=="epPipPim")||(ReacPDF=="eppPim") ){
				par_idx[0] = tmp_fast_p1_idx; // first particle (other than the e') is a proton
				m_par  [0] = mp;
			}

			if(ReacPDF=="epPipPim"){
				par_idx[1] = tmp_fast_pip_idx; // second particle (other than the e') is a π˖
				par_idx[2] = tmp_fast_pim_idx; // third particle (other than the e') is a π˗
				m_par  [1] = mPiC;
				m_par  [2] = mPiC;
			}

			if(ReacPDF=="eppPim"){
				par_idx[1] = tmp_fast_p2_idx ; // second particle (other than the e') is a proton
				par_idx[2] = tmp_fast_pim_idx; // third particle (other than the e') is a π˗
				m_par  [1] = mp;
				m_par  [2] = mPiC;
			}

			// Filling particle variables
			for(int par = 1; par < nParticles; par++){
				V3_par_p   [par-1] = particles.getV3P    (par_idx[par-1]);
				V3_par_v   [par-1] = particles.getV3v    (par_idx[par-1]);
				beta_par   [par-1] = particles.getBeta   (par_idx[par-1]);
				chi2pid_par[par-1] = particles.getChi2pid(par_idx[par-1]);
				V4_par_p   [par-1].SetXYZM( V3_par_p[par-1].X() , V3_par_p[par-1].Y() , V3_par_p[par-1].Z() , m_par[par-1] );
				E_par      [par-1] = V4_par_p[par-1].E();
			}

			// Missing momentum components
			TVector3       V3_Pm = - V3_q;
			TLorentzVector V4_Pm = V4_Ebeam + V4_mtar - V4_ep;
			for(int par = 1; par < nParticles; par++){
				V3_Pm += V3_par_p[par-1];
				V4_Pm -= V4_par_p[par-1];
			}

			double Mmiss  = V4_Pm.Mag();
			double Mmiss2 = V4_Pm.M2();
			double Emiss  = V4_Pm.E();

			// Filling histograms
			// 1D histograms
			h1_pmx         -> Fill(V3_Pm.X()               );
			h1_pmy         -> Fill(V3_Pm.Y()               );
			h1_pmz         -> Fill(V3_Pm.Z()               );
			h1_pm_th       -> Fill(rad2deg*V3_Pm.Theta()   );
			h1_pm_ph       -> Fill(rad2deg*V3_Pm.Phi()     );
			h1_pmiss       -> Fill(V3_Pm.Mag()             );
			h1_Mm          -> Fill(Mmiss                   );
			h1_MmSqr       -> Fill(Mmiss2                  );
			h1_Em          -> Fill(Emiss                   );

			// 2D histograms
			h2_Em_Pm       -> Fill(V3_Pm.Mag(), Emiss      );

			for(int par = 0; par < nParticles-1; par++){
				// 1D histograms
				h1_chi2pid_i[par] -> Fill(chi2pid_par[par]);
				h1_i_vz     [par] -> Fill(V3_par_v[par].Z    ());
				h1_i_dlt_vz [par] -> Fill(V3_par_v[par].Z()-V3_ev.Z());
                		h1_i_px     [par] -> Fill(V3_par_p[par].X    ());
                		h1_i_py     [par] -> Fill(V3_par_p[par].Y    ());
                		h1_i_pz     [par] -> Fill(V3_par_p[par].Z    ());
                		h1_i_p      [par] -> Fill(V3_par_p[par].Mag  ());
                		h1_i_th     [par] -> Fill(V3_par_p[par].Theta());
                		h1_i_phi    [par] -> Fill(V3_par_p[par].Phi  ());

				// 2D histograms
				h2_i_th_phi [par] -> Fill(rad2deg*V3_par_p[par].Phi(),rad2deg*V3_par_p[par].Theta());
               			h2_beta_i_p [par] -> Fill(V3_par_p[par].Mag()        ,beta_par[par]                );
                		h2_i_vz_phi [par] -> Fill(rad2deg*V3_par_p[par].Phi(),V3_par_v[par].Mag()          );
				h2_pe_pi    [par] -> Fill(V3_par_p[par].Mag()        ,ep                           );
			}

			// ----------------------------------------------------------------------
			// Filling tree
			br_e_chi2pid = chi2pid  ;
			br_pex       = V3_ep.X();
			br_pey       = V3_ep.Y();
			br_pez       = V3_ep.Z();
			br_vex       = V3_ev.X();
			br_vey       = V3_ev.Y();
			br_vez       = V3_ev.Z();
			br_qx        = V3_q .X();
			br_qy        = V3_q .Y();
			br_qz        = V3_q .Z();
        		br_Q2        = Q2       ;
        		br_Nu        = omega    ;
        		br_W2        = W2       ;
        		br_xB        = xB       ;

			for(int par = 0; par < nParticles-1; par++){
				br_i_chi2pid[par] = chi2pid_par [par];
				br_pix      [par] = V3_par_p[par].X();
				br_piy      [par] = V3_par_p[par].Y();
				br_piz      [par] = V3_par_p[par].Z();
				br_vix      [par] = V3_par_v[par].X();
				br_viy      [par] = V3_par_v[par].Y();
				br_viz      [par] = V3_par_v[par].Z();
			}

			br_Pmx   = V3_Pm.X()    ;
			br_Pmy   = V3_Pm.Y()    ;
			br_Pmz   = V3_Pm.Z()    ;
			br_Mm    = Mmiss        ;
			br_Mm2   = Mmiss2       ;
			br_thpm  = V3_Pm.Theta();
			br_phipm = V3_Pm.Phi  ();

			tree -> Fill();
		}    

	} // End of while loop over events

	// -------------------------------------------------------------------------------------------
	// Drawing histograms
	TCanvas * c0 = new TCanvas();
	c0 -> Divide(2,1);
	c0 -> cd(1);
	h1_e_vz -> Draw(      );
	for(int par = 0; par < nFSPart-1; par++) h1_i_vz[par] -> Draw("same");
	c0 -> cd(2);
	h1_i_dlt_vz[0] -> Draw(      );
	for(int par = 0; par < nFSPart-1; par++) h1_i_dlt_vz[par] -> Draw("same");
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

	TCanvas ** c7 = new TCanvas*[nFSPart-1];
	for(int par = 0; par < nFSPart-1; par++){
		c7[par] = new TCanvas(Form("c7_%i",par));
		c7[par] -> Divide(2,2);
		c7[par] -> cd(1);	h1_i_px[par] -> Draw();
		c7[par] -> cd(2);	h1_i_py[par] -> Draw();
		c7[par] -> cd(3);	h1_i_pz[par] -> Draw();
		c7[par] -> cd(4);	h1_i_p [par] -> Draw();
		c7[par] -> Modified();
        	c7[par] -> Update();
	}

	TCanvas * c8 = new TCanvas();
	c8->Divide(3,2);
	c8->cd(1);	gPad -> SetLogz();	h2_beta_p_pos -> Draw("COLZ");
	c8->cd(2);	gPad -> SetLogz();	h2_beta_p_neg -> Draw("COLZ");
	for(int par = 0; par < nFSPart-1; par++){
		c8->cd(par+3);
		gPad -> SetLogz();
		h2_beta_i_p[par] -> Draw("COLZ");
	}
	c8 -> Modified();
	c8 -> Update(); 

	TCanvas ** c9 = new TCanvas*[nFSPart-1];
        for(int par = 0; par < nFSPart-1; par++){	
	        c9[par] = new TCanvas(Form("c9_%i",par));
                c9[par] -> Divide(2,2);
		c9[par] -> cd(1);	h1_i_th    [par] -> Draw();
		c9[par] -> cd(2);	h2_i_th_phi[par] -> Draw("COLZ");
		c9[par] -> cd(4);	h1_i_phi   [par] -> Draw();
		c9[par] -> Modified();
                c9[par] -> Update();
        }

	TCanvas * c10 = new TCanvas();
	c10 -> Divide(2, 2);
	c10 -> cd(1);	h1_e_lu -> Draw();
	c10 -> cd(2);	h1_e_lv -> Draw();
	c10 -> cd(3);	h1_e_lw -> Draw();
	c10 -> Modified();
	c10 -> Update();

	TCanvas * c12 = new TCanvas();
	c12 -> Divide(2,1);
	c12 -> cd(1);
	h1_Mm -> Draw();
	c12 -> cd(2);
	h1_MmSqr -> Draw();
	c12 -> Modified();
	c12 -> Update();

	TCanvas * c13 = new TCanvas();
	c13 -> Divide(2,2);
	c13 -> cd(1);	h2_e_vz_phi-> Draw("COLZ");
	for(int par = 0; par < nFSPart-1; par++){
		c13 -> cd(par+2);
		h2_i_vz_phi[par] -> Draw("COLZ");
	}
	c13 -> Modified();
	c13 -> Update();

	TCanvas * c14 = new TCanvas();
	c14 -> Divide(2,1);
	c14 -> cd(1);	h1_e_tof  -> Draw(      );
	c14 -> cd(2);	h2_e_tof_p-> Draw("COLZ");
	c14 -> Modified();
	c14 -> Update();

	TCanvas * c18 = new TCanvas();
	h1_Em -> Draw();
	c18 -> Modified();
	c18 -> Update();

	TCanvas * c19 = new TCanvas();
	h2_Em_Pm -> Draw("COLZ");
	c19 -> Modified();
	c19 -> Update();

	TCanvas * c20 = new TCanvas();	
	for(int par = 0; par < nFSPart-1; par++){
                c20 -> cd(par+1);
                h2_pe_pi[par] -> Draw("COLZ");
        }
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
	for(int par = 0; par < nFSPart-1; par++){
        	c24 -> cd(par+2);
        	h1_chi2pid_i[par] -> Draw();
        } 
	c24 -> Modified();
	c24 -> Update();

	// -------------------------------------------------------------------------------------------
	// Printing some info to screen
	cout << "Triggers: " << ctr_triggers << endl;
	cout << "Electrons: " << ctr_electrons << endl;
	cout << "Selected reaction: " << ctr_selectedReaction << endl;

	// -------------------------------------------------------------------------------------------
	// Saving plots to system
	TString out_pdf_name = "results_" + ReacPDF + ".pdf";
	c0  -> Print(out_pdf_name + "(");
	c1  -> Print(out_pdf_name );
	c2  -> Print(out_pdf_name );
	c3  -> Print(out_pdf_name );
	c5  -> Print(out_pdf_name );
	c6  -> Print(out_pdf_name );
	for(int par = 0; par < nFSPart-1; par++) c7[par] -> Print(out_pdf_name );
	c8  -> Print(out_pdf_name );
	for(int par = 0; par < nFSPart-1; par++) c7[par] -> Print(out_pdf_name );
	c10 -> Print(out_pdf_name );	
	c12 -> Print(out_pdf_name );
	c13 -> Print(out_pdf_name );
	c14 -> Print(out_pdf_name );	
	c18 -> Print(out_pdf_name );
	c19 -> Print(out_pdf_name );
	c20 -> Print(out_pdf_name );
	c21 -> Print(out_pdf_name );
	c22 -> Print(out_pdf_name );
	c23 -> Print(out_pdf_name );
	c24 -> Print(out_pdf_name + ")");

	out -> cd();
	part_pid .Write("particles_pid"   );
	part_mass.Write("particles_mass"  );
	part_char.Write("particles_charge");
	tree -> Write();
	out -> Close();

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
void PrettyTH1(TH1F * h1,int color){
	double color_val = 1;
	if     (color==0) color_val =  2;
	else if(color==1) color_val = 62;
	else if(color==2) color_val =  8;
	else if(color==3) color_val = 93;
	h1 -> SetLineColor(color_val);
        h1 -> SetLineWidth(2);
	h1 -> GetXaxis() -> CenterTitle();
	h1 -> GetYaxis() -> CenterTitle();
	h1 -> GetXaxis() -> SetNdivisions(109);
        h1 -> GetYaxis() -> SetNdivisions(109);
	h1 -> GetXaxis() -> SetTitleSize(0.04);
	h1 -> GetYaxis() -> SetTitleSize(0.04);
	h1 -> GetXaxis() -> SetLabelSize(0.04);
        h1 -> GetYaxis() -> SetLabelSize(0.04);
}
// ========================================================================================================================================
void PrettyTH2(TH2F * h2){
	h2 -> GetXaxis() -> CenterTitle();
        h2 -> GetYaxis() -> CenterTitle();
        h2 -> GetXaxis() -> SetNdivisions(109);
        h2 -> GetYaxis() -> SetNdivisions(109);
        h2 -> GetXaxis() -> SetTitleSize(0.04);
        h2 -> GetYaxis() -> SetTitleSize(0.04);
        h2 -> GetXaxis() -> SetLabelSize(0.04);
        h2 -> GetYaxis() -> SetLabelSize(0.04);
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
