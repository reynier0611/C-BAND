#include <cstdlib>
#include <iostream>
#include <cmath>

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

const int PCal = 1;
const int ECin = 4;
const int ECout = 7;

// Forward-declaring functions
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);

int main(int argc, char** argv) {

	// Useful variables
	double mp      = 0.93827; //GeV (proton mass      )
	double mPiC    = 0.13957; //GeV (charged pion mass)
	double mD      = 1.8756;  //GeV (deuteron mass    )
	double rad2deg = 180./3.14159;

	// Getting input arguments
	double Ebeam = 10.6;
	double mtar = mD;
	if( argc!= 3 ){
		cerr << "Wrong number of arguments. Please instead use:\n"
			<< "\t./pcal_fiducial [inputFile] [outputFile]\n";
		return -1;
	}
	TString inputFile = argv[1];
	TFile * outFile = new TFile(argv[2],"RECREATE");

	TVector3 V3_Ebeam(0,0,Ebeam);
	TLorentzVector V4_Ebeam(V3_Ebeam,Ebeam);
	TLorentzVector V4_mtar(0,0,0,mtar);

	// Declare histograms:
	TH1D ** h_EoP_noCut = new TH1D*[6];
	TH1D ** h_EoP_vw7 = new TH1D*[6];
	TH1D ** h_EoP_vw15 = new TH1D*[6];
	TH1D ** h_PIDchi_noCut = new TH1D*[6];
	TH1D ** h_PIDchi_vw7 = new TH1D*[6];
	TH1D ** h_PIDchi_vw15 = new TH1D*[6];
	TH2D ** h_EoP_u = new TH2D*[6];
	TH2D ** h_EoP_v = new TH2D*[6];
	TH2D ** h_EoP_w = new TH2D*[6];
	TH2D * h_thph_noCut = new TH2D("h_thph_noCut","Electron Theta-Phi, No Cut",80,-30,50,400,0,40);
	TH2D * h_thph_vw7 = new TH2D("h_thph_vw7","Electron Theta-Phi, Loose Cut",80,-30,50,400,0,40);
	TH2D * h_thph_vw15 = new TH2D("h_thph_vw15","Electron Theta-Phi, Tight Cut",80,-30,50,400,0,40);
	for( int i = 0 ; i < 6 ; i ++){
		h_EoP_u[i] = new TH2D(Form("h_EoP_u%i",i),Form("Electron E/p vs u Coord., Sector %i; u [cm]; E/p",(i+1)),450,0,450,400,0,0.4);
		h_EoP_v[i] = new TH2D(Form("h_EoP_v%i",i),Form("Electron E/p vs v Coord., Sector %i; v [cm]; E/p",(i+1)),450,0,450,400,0,0.4);
		h_EoP_w[i] = new TH2D(Form("h_EoP_w%i",i),Form("Electron E/p vs w Coord., Sector %i; w [cm]; E/p",(i+1)),450,0,450,400,0,0.4);

		h_PIDchi_noCut[i]	= new TH1D(Form("h_PIDchi_noCut%i",(i+1)),Form("Electron PID-Chi2, Sector %i; PID-Chi2; Counts",(i+1)),200,-6,6);
		h_PIDchi_noCut[i]->SetLineWidth(2);
		h_PIDchi_vw7[i]		= new TH1D(Form("h_PIDchi_vw7%i",(i+1)),Form("Electron PID-Chi2, Sector %i; PID-Chi2; Counts",(i+1)),200,-6,6);
		h_PIDchi_vw7[i]->SetLineWidth(2);
		h_PIDchi_vw7[i]->SetLineColor(2);
		h_PIDchi_vw15[i]	= new TH1D(Form("h_PIDchi_vw15%i",(i+1)),Form("Electron PID-Chi2, Sector %i; PID-Chi2; Counts",(i+1)),200,-6,6);
		h_PIDchi_vw15[i]->SetLineWidth(2);
		h_PIDchi_vw15[i]->SetLineColor(8);

		h_EoP_noCut[i] 	= new TH1D(Form("h_EoP_noCut_%i",(i+1)),Form("Electron E/p, Sector %i; E/p; Counts",(i+1)),400,0,0.4);
		h_EoP_noCut[i]->SetLineWidth(2);
		h_EoP_vw7[i] 	= new TH1D(Form("h_EoP_vw7_%i",(i+1))  ,Form("Electron E/p, Sector %i; E/p; Counts",(i+1)),400,0,0.4);
		h_EoP_vw7[i]->SetLineWidth(2);
		h_EoP_vw7[i]->SetLineColor(2);
		h_EoP_vw15[i] 	= new TH1D(Form("h_EoP_vw15_%i",(i+1)) ,Form("Electron E/p, Sector %i; E/p; Counts",(i+1)),400,0,0.4);
		h_EoP_vw15[i]->SetLineWidth(2);
		h_EoP_vw15[i]->SetLineColor(8);
	}

	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	BEvent        event       ("REC::Event"       ,reader);
	BParticle     particles   ("REC::Particle"    ,reader);
	BCalorimeter  calo        ("REC::Calorimeter" ,reader);
	BScintillator scintillator("REC::Scintillator",reader);

	int event_counter = 0;
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if( event_counter % 10000 == 0 ) cout << "Working on event: " << event_counter << endl;
		if( event_counter > 100000) break;
		event_counter++;

		// Get electron-particle quantities
		int pid0       = particles.getPid    (0);	// electron candidate id assigned by clas
		TVector3 V3_ev = particles.getV3v    (0);	// electron candidate vertex vector
		TVector3 V3_ep = particles.getV3P    (0);	// electron candidate momentum vector
		float chr0     = particles.getCharge (0);	// electron candidate charge
		float eBeta    = particles.getBeta   (0);	// electron candidate beta = v/c
		float chi2pid  = particles.getChi2pid(0);	// electron candidate goodness of pid fit
		int eStatus    = particles.getStatus (0);	// electron candidate status

		if( pid0!=11 || chr0 != -1 ) continue;

		// Get electron particle momentum
		double ep	= V3_ep.Mag();		// electron candidate momentum magnitude [GeV]

		// Fill histograms
		//calo.show();
		for( int hit = 0 ; hit < calo.getSize() ; hit++ ){
			if( calo.getIndex(hit) == 0 ){ // all the hits here correspond to the ParticleRow 0 == electron
				int secIdx = calo.getSector(hit) - 1;
				double uC = calo.getLU(hit);
				double vC = calo.getLV(hit);
				double wC = calo.getLW(hit);
				double En = calo.getTotE(hit);
				if( En == 0 ) continue;
				

				switch ( calo.getLayer(hit) ){
					//case ECin: continue;
					//case ECout: continue;
					case PCal:
						h_EoP_noCut[secIdx] -> Fill( En / ep );
						h_PIDchi_noCut[secIdx] -> Fill( chi2pid );
						h_thph_noCut -> Fill( V3_ep.Phi()*180./M_PI , V3_ep.Theta()*180./M_PI );
						if( wC < (418.-7.5) && vC < (418.-7.5) ){
							h_EoP_vw7[secIdx]  -> Fill( En / ep );
							h_PIDchi_vw7[secIdx] -> Fill( chi2pid );
							h_thph_vw7 -> Fill( V3_ep.Phi()*180./M_PI, V3_ep.Theta()*180./M_PI );
						}
						if( wC < (418.-15) && vC < (418.-15) ){
							h_EoP_vw15[secIdx] -> Fill( En / ep );
							h_PIDchi_vw15[secIdx] -> Fill( chi2pid );
							h_thph_vw15 -> Fill( V3_ep.Phi()*180./M_PI, V3_ep.Theta()*180./M_PI );
						}

						h_EoP_u[secIdx] -> Fill( uC , En / ep );
						h_EoP_v[secIdx] -> Fill( vC , En / ep );
						h_EoP_w[secIdx] -> Fill( wC , En / ep );
					//default: cerr << "Something happened...\n";
				}
			}
		}

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


	} // End of while loop over events

	outFile->cd();
	TCanvas * pcalU = new TCanvas("pcalU");
	TCanvas * pcalV = new TCanvas("pcalV");
	TCanvas * pcalW = new TCanvas("pcalW");
	TCanvas * epCuts = new TCanvas("epCuts");
	TCanvas * chi2pid = new TCanvas("chi2pid");
	pcalU->Divide(2,3);
	pcalV->Divide(2,3);
	pcalW->Divide(2,3);
	epCuts->Divide(2,3);
	chi2pid->Divide(2,3);
	for( int i = 0 ; i < 6 ; i ++){
		pcalU->cd(i+1);
		h_EoP_u[i] -> Draw("col");
		h_EoP_u[i] -> Write();

		pcalV->cd(i+1);
		h_EoP_v[i] -> Draw("col");
		h_EoP_v[i] -> Write();

		pcalW->cd(i+1);
		h_EoP_w[i] -> Draw("col");
		h_EoP_w[i] -> Write();

		epCuts->cd(i+1);
		h_EoP_noCut[i]->Write();
		h_EoP_vw7[i]->Write();
		h_EoP_vw15[i]->Write();
		int fracLoose = 100* (h_EoP_vw7[i]->Integral() / h_EoP_noCut[i]->Integral());
		int fracTight = 100* (h_EoP_vw15[i]->Integral() / h_EoP_noCut[i]->Integral());
		h_EoP_noCut[i]->SetStats(0); h_EoP_vw7[i]->SetStats(0); h_EoP_vw15[i]->SetStats(0);
		h_EoP_noCut[i]->SetTitle(Form("E/p Projection, Sector %i, No Cut 100%, Loose Cut %i%, Tight Cut %i%",(i+1),fracLoose,fracTight));
		h_EoP_noCut[i]->Draw();
		h_EoP_vw7[i]->Draw("same");
		h_EoP_vw15[i]->Draw("same");

		chi2pid->cd(i+1);
		h_PIDchi_noCut[i]->Write();
		h_PIDchi_vw7[i]->Write();
		h_PIDchi_vw15[i]->Write();
		h_PIDchi_noCut[i]->SetStats(0); h_PIDchi_vw7[i]->SetStats(0); h_PIDchi_vw15[i]->SetStats(0);
		h_PIDchi_noCut[i]->SetTitle(Form("Electron PID-Chi2, Sector %i, No Cut 100%, Loose Cut %i%, Tight Cut %i%",(i+1),fracLoose,fracTight));
		h_PIDchi_noCut[i]->Draw();
		h_PIDchi_vw7[i]->Draw("same");
		h_PIDchi_vw15[i]->Draw("same");

	}
	h_thph_noCut->Write();
	h_thph_vw7->Write();
	h_thph_vw15->Write();
	TCanvas * thetaphi = new TCanvas("thetaphi");
	thetaphi->Divide(3,1);
	thetaphi->cd(1);
	h_thph_noCut->Draw("col");
	thetaphi->cd(2);
	h_thph_vw7->Draw("col");
	thetaphi->cd(3);
	h_thph_vw15->Draw("col");
	thetaphi->Update();
	thetaphi->Modified();
	thetaphi->Write();
	
	

	pcalU->Update(); 
	pcalV->Update(); 
	pcalW->Update(); 
	epCuts->Update();
	chi2pid->Update();
	pcalU->Modified();
	pcalV->Modified();
	pcalW->Modified();
	epCuts->Modified();
	chi2pid->Modified();
	pcalU->Write();
	pcalV->Write();
	pcalW->Write();
	epCuts->Write();
	chi2pid->Write();
	
	outFile->Close();

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
