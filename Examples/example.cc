#include <cstdlib>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "reader.h"
#include "node.h"
#include "bank.h"
//#include "particle.h"
//#include "detector.h"

using namespace std;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);
double fn_Emiss(double Pmiss, double omega, double M_tar, double Enuc, double Mnuc);

// ========================================================================================================================================
int main(int argc, char** argv) {

	std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

	char inputFile[256];
	char outputFile[256];

	if(argc>1) {
		sprintf(inputFile,"%s",argv[1]);
	} else {
		std::cout << " *** please provide a file name..." << std::endl;
		exit(0);
	}

	/*
	   hipo::reader  reader;
	   reader.open(inputFile);

	   while(reader.next()==true){
	   cout << "One more" << endl;

	   int size = particles.getSize();
	   for(int i = 0; i < size; i++){
	   int pid = particles.getPid(i);
	   if(pid==11){
	   particles.getVector4(i,electron,0.0005);
	   clas12::vector4 w2 = beam + target - electron;
	   printf("w2 = %12.5f\n",w2.m());
	   }
	   }

	   }
	 */

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
	// Declaring histograms
	// 1D histograms
	TH1F * h1_e_vz  = new TH1F("h1_e_vz"  ,"h1_e_vz"  ,100, -50, 50);	PrettyTH1F(h1_e_vz  ,"v_{z} [cm]"           ,"Counts",4);
	TH1F * h1_e_tof = new TH1F("h1_e_tof" ,"h1_e_tof" ,100,  20, 27);	PrettyTH1F(h1_e_tof ,"electron TOF [ns]"    ,"Counts",4);

	TH1F * h1_e_px  = new TH1F("h1_e_px"  ,"h1_e_px"  ,100,  -2,  2);	PrettyTH1F(h1_e_px  ,"electron p_{x} [GeV]" ,"Counts",4);
	TH1F * h1_e_py  = new TH1F("h1_e_py"  ,"h1_e_py"  ,100,  -2,  2);	PrettyTH1F(h1_e_py  ,"electron p_{y} [GeV]" ,"Counts",4);
	TH1F * h1_e_pz  = new TH1F("h1_e_pz"  ,"h1_e_pz"  ,100,   0, 10);	PrettyTH1F(h1_e_pz  ,"electron p_{z} [GeV]" ,"Counts",4);
	TH1F * h1_e_p   = new TH1F("h1_e_p"   ,"h1_e_p"   ,100,   0, 10);	PrettyTH1F(h1_e_p   ,"electron |p| [GeV]"   ,"Counts",4);

	TH1F * h1_e_th  = new TH1F("h1_e_th"  ,"h1_e_th"  ,100,   0, 30);	PrettyTH1F(h1_e_th  ,"#theta_e [deg]"       ,"Counts",4);
	TH1F * h1_e_phi = new TH1F("h1_e_phi" ,"h1_e_phi" ,100,-190,190);	PrettyTH1F(h1_e_phi ,"#phi_e [deg]"         ,"Counts",4);
	TH1F * h1_e_lu  = new TH1F("h1_e_lu"  ,"h1_e_lu"  ,100,   0,500);	PrettyTH1F(h1_e_lu  ,"distance on U-side"   ,"Counts",4);
	TH1F * h1_e_lv  = new TH1F("h1_e_lv"  ,"h1_e_lv"  ,100,   0,500);	PrettyTH1F(h1_e_lv  ,"distance on V-side"   ,"Counts",4);
	TH1F * h1_e_lw  = new TH1F("h1_e_lw"  ,"h1_e_lw"  ,100,   0,500);	PrettyTH1F(h1_e_lw  ,"distance on W-side"   ,"Counts",4);

	TH1F * h1_p_vz  = new TH1F("h1_p_vz"  ,"h1_p_vz"  ,100, -50, 50);	PrettyTH1F(h1_p_vz  ,"v_{z} [cm]"           ,"Counts",2);
	TH1F * h1_p_num = new TH1F("h1_p_num" ,"h1_p_num" , 20,   0, 10);	PrettyTH1F(h1_p_num ,"proton number"        ,"Counts",4);

	TH1F * h1_p_px  = new TH1F("h1_p_px"  ,"h1_p_px"  ,100,  -2,  2);	PrettyTH1F(h1_p_px  ,"proton p_{x} [GeV]"   ,"Counts",4);
	TH1F * h1_p_py  = new TH1F("h1_p_py"  ,"h1_p_py"  ,100,  -2,  2);	PrettyTH1F(h1_p_py  ,"proton p_{y} [GeV]"   ,"Counts",4);
	TH1F * h1_p_pz  = new TH1F("h1_p_pz"  ,"h1_p_pz"  ,100,   0,  9);	PrettyTH1F(h1_p_pz  ,"proton p_{z} [GeV]"   ,"Counts",4);
	TH1F * h1_p_p   = new TH1F("h1_p_p"   ,"h1_p_p"   ,100,   0,  9);	PrettyTH1F(h1_p_p   ,"proton |p| [GeV]"     ,"Counts",4);

	TH1F * h1_p_th  = new TH1F("h1_p_th"  ,"h1_p_th"  ,100,   0, 80);	PrettyTH1F(h1_p_th  ,"#theta_p [deg]"       ,"Counts",4);
	TH1F * h1_p_phi = new TH1F("h1_p_phi" ,"h1_p_phi" ,100,-190,190);	PrettyTH1F(h1_p_phi ,"#phi_p [deg]"         ,"Counts",4);

	TH1F * h1_pmiss = new TH1F("h1_pmiss" ,"P_{miss}" ,100,   0,  5);	PrettyTH1F(h1_pmiss ,"Pm [GeV]"             ,"Counts",4);
	TH1F * h1_pmx   = new TH1F("h1_pmx"   ,"h1_pmx"   ,100,  -2,  2);	PrettyTH1F(h1_pmx   ,"Pmx [GeV]"            ,"Counts",4);
	TH1F * h1_pmy   = new TH1F("h1_pmy"   ,"h1_pmy"   ,100,  -2,  2);	PrettyTH1F(h1_pmy   ,"Pmy [GeV]"            ,"Counts",4);
	TH1F * h1_pmz   = new TH1F("h1_pmz"   ,"h1_pmz"   ,100,  -3,  3);	PrettyTH1F(h1_pmz   ,"Pmz [GeV]"            ,"Counts",4);
	TH1F * h1_Mmiss = new TH1F("h1_Mmiss" ,"h1_Mmiss" ,100,   0,  4);	PrettyTH1F(h1_Mmiss ,"m_{miss} [GeV]"       ,"Counts",4);
	TH1F * h1_Em    = new TH1F("h1_Em"    ,"h1_Em"    ,100,  -3,  3);	PrettyTH1F(h1_Em    ,"Em [GeV]"             ,"Counts",4);

	TH1F * h1_W     = new TH1F("h1_W"     ,"h1_W"     ,100,   0,  4);	PrettyTH1F(h1_W     ,"W [GeV]"              ,"Counts",4);
	TH1F * h1_xB    = new TH1F("h1_xB"    ,"h1_xB"    ,100,   0,  4);	PrettyTH1F(h1_xB    ,"x_B"                  ,"Counts",4);

	TH1F * h1_dlt_vz_ep  = new TH1F("h1_dlt_vz_ep"  ,"electron - proton",100, -20, 20);	PrettyTH1F(h1_dlt_vz_ep  ,"#Delta v_{z} [cm]"    ,"Counts",2);
	TH1F * h1_dlt_vz_epip= new TH1F("h1_dlt_vz_epip","electron - #pi+"  ,100, -20, 20);	PrettyTH1F(h1_dlt_vz_epip,"#Delta v_{z} [cm]"    ,"Counts",3);
	TH1F * h1_dlt_vz_epim= new TH1F("h1_dlt_vz_epim","electron - #pi-"  ,100, -20, 20);	PrettyTH1F(h1_dlt_vz_epim,"#Delta v_{z} [cm]"    ,"Counts",4);

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

	PrettyTH2F(h2_e_Ep_p_0  ,"p_{e} [GeV]" ,"E_{e}/p_{e}"         );
	PrettyTH2F(h2_e_Ep_p_1  ,"p_{e} [GeV]" ,"E_{e}/p_{e}"         );
	PrettyTH2F(h2_e_th_phi  ,"#phi_e [deg]","#theta_e [deg]"      );
	PrettyTH2F(h2_p_th_phi  ,"#phi_p [deg]","#theta_p [deg]"      );
	PrettyTH2F(h2_beta_p_pos,"p [GeV]"     ,"#beta"               );
	PrettyTH2F(h2_beta_p_neg,"p [GeV]"     ,"#beta"               );
	PrettyTH2F(h2_beta_p_p  ,"p [GeV]"     ,"#beta"               );
	PrettyTH2F(h2_beta_p_pip,"p [GeV]"     ,"#beta"               );
	PrettyTH2F(h2_beta_p_pim,"p [GeV]"     ,"#beta"               );
	PrettyTH2F(h2_e_vz_phi  ,"#phi_e [deg]","e v_{z} [cm]"        );
	PrettyTH2F(h2_p_vz_phi  ,"#phi_p [deg]","p v_{z} [cm]"        );
	PrettyTH2F(h2_e_tof_p   ,"p_{e} [GeV]" ,"electron TOF [ns]"   );
	PrettyTH2F(h2_p_dtT_p_0 ,"p_{e} [GeV]" ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_p_dtT_p_1 ,"p_{e} [GeV]" ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_p_tof_det ,"detector id" ,"candidate p tof [ns]");
	PrettyTH2F(h2_p_dtT_det ,"detector id" ,"p #Delta t [ns]"     );
	PrettyTH2F(h2_Em_Pm     ,"Pm [GeV]"    ,"Em [GeV]"            );
	PrettyTH2F(h2_pe_pp     ,"p p [GeV]"   ,"e p [GeV]"           );
	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	//hipo::reader reader = new HipoReader();

	// RGA example
	//String dataFile = "/Users/Reynier/WORK/CLAS12/data/cooked/rga/out_clas_003842.evio.295.hipo"; // target lH2
	//String dataFile = "/Users/Reynier/WORK/CLAS12/data/cooked/rga/rga_3842_200_299.hipo"; // target lH2
	//double Ebeam = 6.4; //GeV
	//double mtar    = mp;

	// RGB example
	//String dataFile = "/Users/Reynier/WORK/CLAS12/data/cooked/rgb/hipo_6164_00000_00029.hipo"; // target deuterium
	//String dataFile = "/Users/Reynier/WORK/CLAS12/data/cooked/rgb/hipo_6164_00000_00100.hipo"; // target deuterium
	double Ebeam = 10.6; //GeV
	double mtar    = mD;

	hipo::reader  reader;
	reader.open(inputFile);

	hipo::bank bank_particle    ("REC::Particle"    ,reader);
	hipo::bank bank_calorimeter ("REC::Calorimeter" ,reader);
	hipo::bank bank_scintillator("REC::Scintillator",reader);
	hipo::bank bank_event       ("REC::Event"       ,reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){
		//bank_particle.show();
		//bank_event.show();
		bank_scintillator.show();

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;

		// Particle bank
		int    pid0    = bank_particle.getInt  ( 0,0);	// electron candidate id assigned by clas
		double epx     = bank_particle.getFloat( 1,0);	// electron candidate momentum x-component [GeV]
		double epy     = bank_particle.getFloat( 2,0);	// electron candidate momentum y-component [GeV]
		double epz     = bank_particle.getFloat( 3,0);	// electron candidate momentum z-component [GeV]
		double evx     = bank_particle.getFloat( 4,0);	// electron candidate vertex x coordinate [cm]
		double evy     = bank_particle.getFloat( 5,0);	// electron candidate vertex y coordinate [cm]
		double evz     = bank_particle.getFloat( 6,0);	// electron candidate vertex z coordinate [cm]
		int    chr0    = bank_particle.getInt  ( 7,0);	// electron candidate charge
		double eBeta   = bank_particle.getFloat( 8,0);	// electron candidate beta = v/c
		double chi2pid = bank_particle.getFloat( 9,0);	// electron candidate goodness of PID from EB
		double eStatus = bank_particle.getInt  (10,0);	// electron candidate status (whatever this means)

		// Event bank
		double t_vtx   = bank_event.getFloat( 9,0);

		// Scintillator bank
                double t_e     = bank_scintillator.getFloat(7,0);

		short cal_id = 0;
		double Epcal = 0;
		double Eecin = 0;
		double Eecou = 0;

		/*
		// Calorimeter bank
		int nCal = bank_calorimeter.getNode("pindex").getDataSize();

		double lU      = bank_calorimeter.getNode("lu"    ).getFloat(0);	// electron candidate distance on U-side [cm?]
		double lV      = bank_calorimeter.getNode("lv"    ).getFloat(0);	// electron candidate distance on V-side [cm?]
		double lW      = bank_calorimeter.getNode("lw"    ).getFloat(0);	// electron candidate distance on W-side [cm?]

		for(int cal=0 ; cal < nCal ; cal++) {
		cal_id = bank_calorimeter.getNode("pindex").getShort(cal);
		int layer = bank_calorimeter.getNode("layer").getByte(cal);
		if     (cal_id==0&&layer==1) Epcal = bank_calorimeter.getNode("energy").getFloat(cal);	// PCal
		else if(cal_id==0&&layer==4) Eecin = bank_calorimeter.getNode("energy").getFloat(cal);	// ECinner
		else if(cal_id==0&&layer==7) Eecou = bank_calorimeter.getNode("energy").getFloat(cal);	// ECouter
		}

		double Ee     = Epcal + Eecin + Eecou;	// electron candidate energy from calorimeter [GeV]

		// calculated variables
		double ep     = TMath::Sqrt(epx*epx + epy*epy + epz*epz);		// electron candidate momentum magnitude [GeV]
		Vector3 v3_ep = new Vector3(epx,epy,epz);					// electron candidate momentum vector [GeV]
		double th_e   = v3_ep.theta();			// electron candidate theta [rad]
		double phi_e  = v3_ep.phi();								// electron candidate phi [rad]

		double tof_e  = t_e - t_vtx;								// electron candidate time-of-flight [ns]

		// Transfer variables
		double qx = - epx;											// q vector x component [GeV]
		double qy = - epy;											// q vector y component [GeV]
		double qz = Ebeam - epz;									// q vector z component [GeV]
		double Q2     = 4*ep*Ebeam*Math.pow(Math.sin(th_e/2.),2);	// Q-squared [GeV^2]
		double nu     = Ebeam - ep;									// Transfer energy [GeV]
		//double W2     = mtar*mtar-Q2+2*nu*mtar;						// Invariant mass ^2 [GeV]
		double W2     = mp*mp-Q2+2*nu*mp;
		double xB     = Q2/2./mp/nu;								// Bjorken-x

		// -------------------------------------------------------------------------
		// Fill some histograms before cutting on good electrons
		h2_e_Ep_p_0.fill(ep,Ee/ep);
		h2_e_vz_phi.fill(rad2deg*phi_e, evz);

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
		if(     (pid0!=11            )||
				(chr0!=-1            )||
				(chi2pid>=cut_chi2pid)||
				(ep<=cut_ep          )||
				(ep>=Ebeam           )||
				(evz>cut_max_vz      )||
				(evz<cut_min_vz      )||
				(lU<cut_uvw          )||
				(lV<cut_uvw          )||
				(lW<cut_uvw          )||
				(Epcal<cut_Epcal     )||
				(Math.sqrt(W2)<=cut_W)||
				(tof_e<cut_tof_e     )
		  ) continue;

		// -------------------------------------------------------------------------
		// Filling electron histograms
		h1_e_lu .fill(lU           );
		h1_e_lv .fill(lV           );
		h1_e_lw .fill(lW           );
		h1_e_px .fill(epx          );
		h1_e_py .fill(epy          );
		h1_e_pz .fill(epz          );
		h1_e_p  .fill(ep           );
		h1_e_vz .fill(evz          );
		h1_W    .fill(Math.sqrt(W2));
		h1_xB   .fill(xB           );
		h1_e_th .fill(rad2deg*th_e );
		h1_e_phi.fill(rad2deg*phi_e);
		h1_e_tof.fill(tof_e        );

		h2_e_Ep_p_1.fill(ep           ,Ee/ep       );
		h2_e_th_phi.fill(rad2deg*phi_e,rad2deg*th_e);
		h2_e_tof_p .fill(ep           ,tof_e       );

		// Saving these events in output file
		writer.writeEvent(event);

		// -------------------------------------------------------------------------
		// Loop over the remaining particles and require a proton, a pi+, and a pi-

		// Proton variables
		int nProtons = 0;
		int tmp_fast_p_idx  = 0;	// index of the fastest proton
		double tmp_fast_p_p = 0;	// momentum of the fastest proton

		// pi+ variables
		int nPip = 0;
		int tmp_fast_pip_idx  = 0;	// index of the fastest proton
		double tmp_fast_pip_p = 0;	// momentum of the fastest proton

		// pi- variables
		int nPim = 0;
		int tmp_fast_pim_idx  = 0;	// index of the fastest proton
		double tmp_fast_pim_p = 0;	// momentum of the fastest proton

		// int size = particles.getSize(); << new way in C++
		int nParticles = bank_particle.getNode("pid").getDataSize();
		for(int par = 1; par < nParticles; par++) {
			// Particle bank
			int pidi      = bank_particle.getNode("pid"   ).getInt  (par);		// proton candidate id assigned by clas
			int chri      = bank_particle.getNode("charge").getInt  (par);		// proton candidate charge
			double beta_p = bank_particle.getNode("beta"  ).getFloat(par);		// proton candidate beta (v/c)
			double ppx    = bank_particle.getNode("px"    ).getFloat(par);		// proton candidate momentum x-component [GeV]
			double ppy    = bank_particle.getNode("py"    ).getFloat(par);		// proton candidate momentum y-component [GeV]
			double ppz    = bank_particle.getNode("pz"    ).getFloat(par);		// proton candidate momentum z-component [GeV]
			double pvx    = bank_particle.getNode("vx"    ).getFloat(par);		// proton candidate vertex x coordinate [cm]
			double pvy    = bank_particle.getNode("vy"    ).getFloat(par);		// proton candidate vertex y coordinate [cm]
			double pvz    = bank_particle.getNode("vz"    ).getFloat(par);		// proton candidate vertex z coordinate [cm]

			double pp     = Math.sqrt(ppx*ppx + ppy*ppy + ppz*ppz);		

			double mSq    = pp*pp*(1-beta_p*beta_p)/(beta_p*beta_p);
			boolean alreadyGotHit = false;
			double t_p      = -1000;
			double tof_p    = -1000;
			double delta_tP = -1000;
			int detectorID = -1;

			// Scintillator bank
			int nScin = bank_scintillator.getNode("pindex").getDataSize();
			for(int scin = 1 ; scin < nScin ; scin++ ) {
				int scint_id = bank_scintillator.getNode("pindex").getInt(scin);
				if(scint_id==par&&!alreadyGotHit) {
					alreadyGotHit=true;
					detectorID = bank_scintillator.getNode("detector").getInt  (scin);
					t_p        = bank_scintillator.getNode("time"    ).getFloat(scin);
					tof_p  = t_p - t_vtx;
					delta_tP = tof_p*(1-Math.sqrt((pp*pp+mp*mp)/(pp*pp+mSq)));
				}
			}


			if     (chri== 1) h2_beta_p_pos.fill(pp, beta_p  );
			else if(chri==-1) h2_beta_p_neg.fill(pp, beta_p  );

			if(delta_tP!=-1000&&chri>0) {
				h2_p_dtT_p_0.fill(pp, delta_tP);
				h2_p_tof_det.fill((double)(detectorID),tof_p);
				h2_p_dtT_det.fill((double)(detectorID),delta_tP);
			}
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
			if(     (pidi==2212)&&
					(chri==1   )&&
					(pp<Ebeam  )
			  ) {
				nProtons++;
				if(pp>tmp_fast_p_p) {
					tmp_fast_p_p = pp;
					tmp_fast_p_idx=par;
				}
			}
			// pi+
			if(     (pidi==211)&&
					(chri==1   )&&
					(pp<Ebeam  )
			  ) {
				nPip++;
				if(pp>tmp_fast_pip_p) {
					tmp_fast_pip_p = pp;
					tmp_fast_pip_idx=par;
				}
			}
			// pi-
			if(     (pidi==-211)&&
					(chri==-1  )&&
					(pp<Ebeam  )
			  ) {
				nPim++;
				if(pp>tmp_fast_pim_p) {
					tmp_fast_pim_p = pp;
					tmp_fast_pim_idx=par;
				}
			}
			// ----------------------------------------------------------------------
		}
		if(nProtons==1&&nPip==1&&nPim==1) {

			double beta_p  = bank_particle.getNode("beta").getFloat(tmp_fast_p_idx  );	// proton candidate beta (v/c)
			double ppx     = bank_particle.getNode("px"  ).getFloat(tmp_fast_p_idx  );	// proton candidate momentum x-component [GeV]
			double ppy     = bank_particle.getNode("py"  ).getFloat(tmp_fast_p_idx  );	// proton candidate momentum y-component [GeV]
			double ppz     = bank_particle.getNode("pz"  ).getFloat(tmp_fast_p_idx  );	// proton candidate momentum z-component [GeV]
			double pvz     = bank_particle.getNode("vz"  ).getFloat(tmp_fast_p_idx  );	// proton candidate vertex z coordinate [cm]

			double beta_pip= bank_particle.getNode("beta").getFloat(tmp_fast_pip_idx);	// proton candidate beta (v/c)
			double pip_px  = bank_particle.getNode("px"  ).getFloat(tmp_fast_pip_idx);	// proton candidate momentum x-component [GeV]
			double pip_py  = bank_particle.getNode("py"  ).getFloat(tmp_fast_pip_idx);	// proton candidate momentum y-component [GeV]
			double pip_pz  = bank_particle.getNode("pz"  ).getFloat(tmp_fast_pip_idx);	// proton candidate momentum z-component [GeV]
			double pip_vz  = bank_particle.getNode("vz"  ).getFloat(tmp_fast_pip_idx);	// proton candidate vertex z coordinate [cm]

			double beta_pim= bank_particle.getNode("beta").getFloat(tmp_fast_pim_idx);	// proton candidate beta (v/c)
			double pim_px  = bank_particle.getNode("px"  ).getFloat(tmp_fast_pim_idx);	// proton candidate momentum x-component [GeV]
			double pim_py  = bank_particle.getNode("py"  ).getFloat(tmp_fast_pim_idx);	// proton candidate momentum y-component [GeV]
			double pim_pz  = bank_particle.getNode("pz"  ).getFloat(tmp_fast_pim_idx);	// proton candidate momentum z-component [GeV]
			double pim_vz  = bank_particle.getNode("vz"  ).getFloat(tmp_fast_pim_idx);	// proton candidate vertex z coordinate [cm]

			double pp     = Math.sqrt(ppx*ppx + ppy*ppy + ppz*ppz);
			double Ep     = Math.sqrt(pp*pp+mp*mp);
			Vector3 v3_pp = new Vector3(ppx,ppy,ppz);							// proton candidate momentum vector [GeV]		
			double th_p   = v3_pp.theta();										// proton candidate theta [rad]
			double phi_p  = v3_pp.phi();										// proton candidate phi [rad]

			double pPip   = Math.sqrt(pip_px*pip_px + pip_py*pip_py + pip_pz*pip_pz);
			double EPip   = Math.sqrt(pPip*pPip+mPiC*mPiC);

			double pPim   = Math.sqrt(pim_px*pim_px + pim_py*pim_py + pim_pz*pim_pz);
			double EPim   = Math.sqrt(pPim*pPim+mPiC*mPiC);		

			// Missing momentum components
			double pmx = ppx + pip_px + pim_px - qx;
			double pmy = ppy + pip_py + pim_py - qy;
			double pmz = ppz + pip_pz + pim_pz - qz;
			double Pm = Math.sqrt(pmx*pmx + pmy*pmy + pmz*pmz);

			// Missing mass
			double E_mmiss = Ebeam + mtar - ep - Ep - EPip - EPim;
			double Mmiss = Math.sqrt(E_mmiss*E_mmiss - Pm*Pm);

			//double Emiss = fn_Emiss( Pm, nu, mtar, Ep, mp);

			if(     (Math.abs(pvz    - evz + 4.219763e-01) < 3*5.139898e+00)&&
					(Math.abs(pip_vz - evz - 2.137405e+00) < 3*4.259259e+00)&&
					(Math.abs(pim_vz - evz - 2.230033e+00) < 3*3.190657e+00)
			  ) {

				h1_p_vz       .fill(pvz          );
				h1_dlt_vz_ep  .fill(pvz    - evz );
				h1_dlt_vz_epip.fill(pip_vz - evz );
				h1_dlt_vz_epim.fill(pim_vz - evz );
				h1_p_px       .fill(ppx          );
				h1_p_py       .fill(ppy          );
				h1_p_pz       .fill(ppz          );
				h1_p_p        .fill(pp           );
				h1_p_th       .fill(rad2deg*th_p );
				h1_p_phi      .fill(rad2deg*phi_p);
				h1_pmx        .fill(pmx          );
				h1_pmy        .fill(pmy          );
				h1_pmz        .fill(pmz          );
				h1_pmiss      .fill(Pm           );
				h1_Mmiss      .fill(Mmiss        );
				//h1_Em      .fill(Emiss        );

				h2_p_th_phi .fill(rad2deg*phi_p, rad2deg*th_p);
				h2_p_vz_phi .fill(rad2deg*phi_p, pvz         );

				h2_beta_p_p  .fill(pp           , beta_p      );
				h2_beta_p_pip.fill(pPip         , beta_pip    );
				h2_beta_p_pim.fill(pPim         , beta_pim    );

				//h2_p_dtT_p_1.fill(pp           , delta_tP    );

				//h2_Em_Pm    .fill(Pm           , Emiss       );
				h2_pe_pp    .fill(pp           , ep          );
			}
		}


		h1_p_num.fill((double)(nProtons));


	*/
	} // End of while loop over events

	//writer.close();

	/*
	// -------------------------------------------------------------------------------------------
	// Drawing histograms
	TCanvas c0 = new TCanvas("c0", 1024, 768);
	c0.divide(2,1);
	c0.cd(0);
	c0.getPad().setLegend(true);
	c0.draw(h1_e_vz);
	c0.draw(h1_p_vz,"same");
	c0.cd(1);
	c0.getPad().setLegend(true);
	c0.draw(h1_dlt_vz_ep);
	c0.draw(h1_dlt_vz_epip,"same");
	c0.draw(h1_dlt_vz_epim,"same");

	TCanvas c1 = new TCanvas("c1", 1024, 768);
	c1.divide(2, 2);
	c1.cd(0);	c1.draw(h1_pmx);
	c1.cd(1);	c1.draw(h1_pmy);
	c1.cd(2);	c1.draw(h1_pmz);
	c1.cd(3);	c1.draw(h1_pmiss);

	TCanvas c2 = new TCanvas("c2", 1024, 768);
	c2.divide(2, 2);
	c2.cd(0);	c2.draw(h1_e_th);
	c2.cd(1);	c2.draw(h2_e_th_phi);
	c2.cd(3);	c2.draw(h1_e_phi);

	TCanvas c3 = new TCanvas("c3", 1024, 768);
	c3.divide(2, 1);
	c3.cd(0);	c3.draw(h1_W);
	c3.cd(1);	c3.draw(h1_xB);

	EmbeddedCanvas c5 = new EmbeddedCanvas();
	c5.divide(2, 1);
	c5.cd(0);	c5.draw(h2_e_Ep_p_0);	c5.getPad(0).getAxisZ().setLog(true);
	c5.cd(1);	c5.draw(h2_e_Ep_p_1);	c5.getPad(1).getAxisZ().setLog(true);
	CanvasMaker(c5,"c5",800,600);

	TCanvas c6 = new TCanvas("c6", 1024, 768);
	c6.divide(2, 2);
	c6.cd(0);	c6.draw(h1_e_px);
	c6.cd(1);	c6.draw(h1_e_py);
	c6.cd(2);	c6.draw(h1_e_pz);
	c6.cd(3);	c6.draw(h1_e_p );

	TCanvas c7 = new TCanvas("c7", 1024, 768);
	c7.divide(2, 2);
	c7.cd(0);	c7.draw(h1_p_px);
	c7.cd(1);	c7.draw(h1_p_py);
	c7.cd(2);	c7.draw(h1_p_pz);
	c7.cd(3);	c7.draw(h1_p_p );

	EmbeddedCanvas c8 = new EmbeddedCanvas();
	c8.divide(3, 2);
	c8.cd(0);	c8.draw(h2_beta_p_pos);	c8.getPad(0).getAxisZ().setLog(true);
	c8.cd(1);	c8.draw(h2_beta_p_p	 );	c8.getPad(1).getAxisZ().setLog(true);
	c8.cd(2);	c8.draw(h2_beta_p_pip);	c8.getPad(2).getAxisZ().setLog(true);
	c8.cd(3);	c8.draw(h2_beta_p_neg);	c8.getPad(3).getAxisZ().setLog(true);
	c8.cd(4);	c8.draw(h2_beta_p_pim);	c8.getPad(4).getAxisZ().setLog(true);
	CanvasMaker(c8,"c8",800,600);

	TCanvas c9 = new TCanvas("c9", 1024, 768);
	c9.divide(2, 2);
	c9.cd(0);	c9.draw(h1_p_th);
	c9.cd(1);	c9.draw(h2_p_th_phi);
	c9.cd(3);	c9.draw(h1_p_phi);

	TCanvas c10 = new TCanvas("c10", 1024, 768);
	c10.divide(2, 2);
	c10.cd(0);	c10.draw(h1_e_lu);
	c10.cd(1);	c10.draw(h1_e_lv);
	c10.cd(2);	c10.draw(h1_e_lw);

	TCanvas c11 = new TCanvas("c11", 1024, 768);
	c11.draw(h1_p_num);

	TCanvas c12 = new TCanvas("c12", 1024, 768);
	c12.draw(h1_Mmiss);

	TCanvas c13 = new TCanvas("c13", 1024, 768);
	c13.divide(2, 1);
	c13.cd(0);	c13.draw(h2_e_vz_phi);
	c13.cd(1);	c13.draw(h2_p_vz_phi);

	TCanvas c14 = new TCanvas("c14", 1024, 768);
	c14.divide(2, 1);
	c14.cd(0);	c14.draw(h1_e_tof  );
	c14.cd(1);	c14.draw(h2_e_tof_p);

	EmbeddedCanvas c15 = new EmbeddedCanvas();
	c15.divide(2, 1);
	c15.getPad(0).getAxisZ().setLog(true);	c15.cd(0);	c15.draw(h2_p_dtT_p_0);
	c15.getPad(1).getAxisZ().setLog(true);	c15.cd(1);	c15.draw(h2_p_dtT_p_1);
	CanvasMaker(c15,"c15",800,600);

	TCanvas c16 = new TCanvas("c16", 1024, 768);
	c16.draw(h2_p_tof_det);

	TCanvas c17 = new TCanvas("c17", 1024, 768);
	c17.draw(h2_p_dtT_det);

	TCanvas c18 = new TCanvas("c18", 1024, 768);
	c18.draw(h1_Em);

	TCanvas c20 = new TCanvas("c20", 1024, 768);
	c20.draw(h2_Em_Pm);

	TCanvas c21 = new TCanvas("c21", 1024, 768);
	c21.draw(h2_pe_pp);

	// -------------------------------------------------------------------------------------------
	// Saving plots to system
	c0 .save("/Users/Reynier/WORK/CLAS12/data/c0.png" );
	c1 .save("/Users/Reynier/WORK/CLAS12/data/c1.png" );
	c2 .save("/Users/Reynier/WORK/CLAS12/data/c2.png" );
	c3 .save("/Users/Reynier/WORK/CLAS12/data/c3.png" );
	c5 .save("/Users/Reynier/WORK/CLAS12/data/c5.png" );
	c6 .save("/Users/Reynier/WORK/CLAS12/data/c6.png" );
	c7 .save("/Users/Reynier/WORK/CLAS12/data/c7.png" );
	c8 .save("/Users/Reynier/WORK/CLAS12/data/c8.png" );
	c9 .save("/Users/Reynier/WORK/CLAS12/data/c9.png" );
	c10.save("/Users/Reynier/WORK/CLAS12/data/c10.png");
	c11.save("/Users/Reynier/WORK/CLAS12/data/c11.png");
	c12.save("/Users/Reynier/WORK/CLAS12/data/c12.png");
	c13.save("/Users/Reynier/WORK/CLAS12/data/c13.png");
	c14.save("/Users/Reynier/WORK/CLAS12/data/c14.png");
	c15.save("/Users/Reynier/WORK/CLAS12/data/c15.png");
	c16.save("/Users/Reynier/WORK/CLAS12/data/c16.png");
	c17.save("/Users/Reynier/WORK/CLAS12/data/c17.png");
	c18.save("/Users/Reynier/WORK/CLAS12/data/c18.png");

	c20.save("/Users/Reynier/WORK/CLAS12/data/c20.png");
	c21.save("/Users/Reynier/WORK/CLAS12/data/c21.png");
	*/

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
