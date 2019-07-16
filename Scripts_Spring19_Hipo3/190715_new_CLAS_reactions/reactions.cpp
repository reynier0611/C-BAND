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
	TString inputFile, Reac, ReacPDF;
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
		}
		else if(atoi(argv[2])==2){
			Reac = "(e,e'π˖)";	ReacPDF = "ePip";
			nFSPart = 2;
		}
		else if(atoi(argv[2])==3){
			Reac = "(e,e'pπ˖π˗)";	ReacPDF = "epPipPim";
			nFSPart = 4;
		}
		else if(atoi(argv[2])==4){	
			Reac = "e,e'pp π˗";	ReacPDF = "eppPim";
			nFSPart = 4;
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
		cout << "         = 4 -> (e,e'pp π˗)" << endl;
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
        // Setting up output root file
        TFile * out = new TFile("outRoot_"+ReacPDF+".root","RECREATE");
        TTree * tree = new TTree("T","(e,"+ReacPDF+")");

        double br_pex, br_pey, br_pez;
	double br_qx , br_qy , br_qz ;
        double br_Q2 , br_Nu , br_W2 , br_xB;
	tree -> Branch("br_pex", &br_pex , "br_pex/D");
	tree -> Branch("br_pey", &br_pey , "br_pey/D");
	tree -> Branch("br_pez", &br_pez , "br_pez/D");
	tree -> Branch("br_qx" , &br_qx  , "br_qx/D" );
	tree -> Branch("br_qy" , &br_qy  , "br_qy/D" );
	tree -> Branch("br_qz" , &br_qz  , "br_qz/D" );
	tree -> Branch("br_Q2" , &br_Q2  , "br_Q2/D" );
	tree -> Branch("br_Nu" , &br_Nu  , "br_Nu/D" );
	tree -> Branch("br_W2" , &br_W2  , "br_W2/D" );
	tree -> Branch("br_xB" , &br_xB  , "br_xB/D" );

	double br_pix[nFSPart], br_piy[nFSPart], br_piz[nFSPart];

	for(int par = 0 ; par < nFSPart-1 ; par++){
		tree -> Branch(Form("br_p%ix",par), &br_pix[par] , Form("br_p%ix/D",par));
		tree -> Branch(Form("br_p%iy",par), &br_piy[par] , Form("br_p%iy/D",par));
		tree -> Branch(Form("br_p%iz",par), &br_piz[par] , Form("br_p%iz/D",par));
	}

	double br_Pmx, br_Pmy, br_Pmz, br_Mm, br_Mm2;
	tree -> Branch("br_Pmx", &br_Pmx , "br_Pmx/D");
	tree -> Branch("br_Pmy", &br_Pmy , "br_Pmy/D");
	tree -> Branch("br_Pmz", &br_Pmz , "br_Pmz/D");
	tree -> Branch("br_Mm" , &br_Mm  , "br_Mm/D" );
	tree -> Branch("br_Mm2", &br_Mm2 , "br_Mm2/D");

	TVectorT<double> part_mass(nFSPart-1);
	TVectorT<double> part_char(nFSPart-1);

	if(ReacPDF=="ePip"){
		part_mass[0] = mPiC;
		part_char[0] = 1   ;
	}
        if( (ReacPDF=="ep")||(ReacPDF=="epPipPim")||(ReacPDF=="eppPim") ){
       		part_mass[0] = mp; 
        	part_char[0] = 1 ;
	}
        if(ReacPDF=="epPipPim"){
    		part_mass[1] = mPiC;
		part_char[1] = 1   ;
		part_mass[2] = mPiC;       
        	part_char[2] = -1  ;
	}
        if(ReacPDF=="eppPim"){
		part_mass[1] = mp  ;
                part_char[1] = 1   ;
		part_mass[2] = mPiC;
		part_char[2] = -1  ;
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
			if((chi2pidi<cut_chi2pid)&&(V3_pi.Mag()<Ebeam)){
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

			/*
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
			*/
			h1_Mm          -> Fill(Mmiss                   );
			h1_MmSqr       -> Fill(Mmiss2                  );
			/*
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
			 */
		
			// Filling tree
			br_pex = V3_ep.X();
			br_pey = V3_ep.Y();
			br_pez = V3_ep.Z();
			br_qx  = V3_q .X();
			br_qy  = V3_q .Y();
			br_qz  = V3_q .Z();
        		br_Q2  = Q2       ;
        		br_Nu  = omega    ;
        		br_W2  = W2       ;
        		br_xB  = xB       ;

			for(int par = 0; par < nParticles-1; par++){
				br_pix[par] = V3_par_p[par].X();
				br_piy[par] = V3_par_p[par].Y();
				br_piz[par] = V3_par_p[par].Z();
			}

			br_Pmx = V3_Pm.X();
			br_Pmy = V3_Pm.Y();
			br_Pmz = V3_Pm.Z();
			br_Mm  = Mmiss    ;
			br_Mm2 = Mmiss2   ;

			tree -> Fill();
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
	TString out_pdf_name = "results_" + ReacPDF + ".pdf";
	c0  -> Print(out_pdf_name + "(");
	c1  -> Print(out_pdf_name );
	c2  -> Print(out_pdf_name );
	c3  -> Print(out_pdf_name );
	c5  -> Print(out_pdf_name );
	c6  -> Print(out_pdf_name );
	c7  -> Print(out_pdf_name );
	c8  -> Print(out_pdf_name );
	c9  -> Print(out_pdf_name );
	c10 -> Print(out_pdf_name );
	c11 -> Print(out_pdf_name );
	c12 -> Print(out_pdf_name );
	c13 -> Print(out_pdf_name );
	c14 -> Print(out_pdf_name );
	c15 -> Print(out_pdf_name );
	c16 -> Print(out_pdf_name );
	c17 -> Print(out_pdf_name );
	c18 -> Print(out_pdf_name );
	c19 -> Print(out_pdf_name );
	c20 -> Print(out_pdf_name );
	c21 -> Print(out_pdf_name );
	c22 -> Print(out_pdf_name );
	c23 -> Print(out_pdf_name );
	c24 -> Print(out_pdf_name + ")");

	out -> cd();
	part_mass.Write("particles_mass");
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
