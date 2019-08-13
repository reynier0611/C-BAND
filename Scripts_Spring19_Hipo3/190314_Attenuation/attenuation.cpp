#include <cstdlib>
#include <iostream>

#include "TRint.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TProfile.h"

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
void PrettyTH1F(TH1F * h1, int color);
void PrettyTH2F(TH2F * h2);
double getTriggerPhase( long timeStamp );
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);
void LoadEffectiveVelocity();
void LoadEfrainsAttenuationLengths();

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

double mean_c = 14; //cm/ns speed of light in the bar
double eff_c    [600] = {0};
double mu_efrain[600] = {0};
// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetTitleSize(0.3,"t");
	gStyle->SetOptStat(0);

	TString inputFile;

	if(argc==2) {
		inputFile = argv[1];
	}
	else {
		cout << "=========================\nRun this code as:\n./code path/to/input/file\n=========================" << endl;
		exit(0);
	}

	LoadEffectiveVelocity        ();
	LoadEfrainsAttenuationLengths();
	
	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	double par[nHistos][2] = {{0}};

	TH2F ** h2_T_vDt = new TH2F * [nHistos];
	TH1F ** h1_ADC_L = new TH1F * [nHistos];
	TH1F ** h1_ADC_R = new TH1F * [nHistos];

	for(int i = 0 ; i < nHistos ; i++){
		h2_T_vDt[i] = new TH2F(Form("h2_T_vDt_%i",i),";position [m];ln(R)",200,-1.5,1.5,200,-1,1);
		PrettyTH2F(h2_T_vDt[i]);
	
		h1_ADC_L[i] = new TH1F(Form("h1_ADC_L_%i",i),";ADC",100,0,50000);	PrettyTH1F(h1_ADC_L[i] , 2);
        	h1_ADC_R[i] = new TH1F(Form("h1_ADC_R_%i",i),";ADC",100,0,50000);	PrettyTH1F(h1_ADC_R[i] ,62);
	}

	TH1F * h1_E_minus_R_mu = new TH1F("h1_E_minus_R_mu",";#mu (Efrain - Rey) [m]",40,-2,2); // Attenuation lengths from Efrain - from Rey

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	hipo::bank BAND_ADC  ("BAND::adc"  ,reader);
	hipo::bank BAND_TDC  ("BAND::tdc"  ,reader);
	hipo::bank RUN_config("RUN::config",reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000==0) cout << "event: " << event_counter << endl;
		event_counter++;	

		//BAND_ADC.show();
		//BAND_TDC.show();
		//RUN_config.show();

		int nADC = BAND_ADC.getSize();
		int nTDC = BAND_TDC.getSize();

		// Skip events with no entries
		if(nADC==0||nTDC==0) continue;

		// Skip events with less than 200 entries (the laser lights all bars at the same time)
		//if(nADC<200||nTDC<200) continue;

		// Skip all laser events
		if(nADC>100||nTDC>100) continue;

		//long timestamp = RUN_config.getLong(4,0);
		//double phaseCorr = getTriggerPhase(timestamp);

		for(int aIdx1 = 0 ; aIdx1 < nADC ; aIdx1++){

			bool matchedBarADC = false;

			int   ADC1_sector    = BAND_ADC.getInt  (0,aIdx1);
			int   ADC1_layer     = BAND_ADC.getInt  (1,aIdx1);
			int   ADC1_component = BAND_ADC.getInt  (2,aIdx1);
			int   ADC1_order     = BAND_ADC.getInt  (3,aIdx1);
			float ADC1_adc       = (float)(BAND_ADC.getInt(4,aIdx1));
			float ADC1_time      = BAND_ADC.getFloat(5,aIdx1);

			for(int aIdx2 = 0 ; aIdx2 < nADC ; aIdx2++){

				int   ADC2_sector    = BAND_ADC.getInt  (0,aIdx2);
				int   ADC2_layer     = BAND_ADC.getInt  (1,aIdx2);
				int   ADC2_component = BAND_ADC.getInt  (2,aIdx2);
				int   ADC2_order     = BAND_ADC.getInt  (3,aIdx2);
				float ADC2_adc       = (float)(BAND_ADC.getInt(4,aIdx2));
				float ADC2_time      = BAND_ADC.getFloat(5,aIdx2);

				// Matching each ADC to the corresponding ADC from the other side of the bar
				if(             (!matchedBarADC                       )&&
						(ADC1_sector   ==ADC2_sector          )&&
						(ADC1_layer    ==ADC2_layer           )&&
						(ADC1_component==ADC2_component       )&&
						(ADC1_order+1==ADC2_order             )
				  ){
					matchedBarADC = true;

					bool matchedBarTDC = false;

					for(int tIdx1 = 0 ; tIdx1 < nTDC ; tIdx1++){
						int   TDC1_sector    = BAND_TDC.getInt  (0,tIdx1);
						int   TDC1_layer     = BAND_TDC.getInt  (1,tIdx1);
						int   TDC1_component = BAND_TDC.getInt  (2,tIdx1);
						int   TDC1_order     = BAND_TDC.getInt  (3,tIdx1);
						float TDC1_tdc       = (float)(BAND_TDC.getInt(4,tIdx1));
						float TDC1_time      = TDC1_tdc*0.02345;

						//TDC1_time -= phaseCorr;

						// Matching these ADCs to a TDC
						if(             (ADC1_sector   ==TDC1_sector   )&&
								(ADC1_layer    ==TDC1_layer    )&&
								(ADC1_component==TDC1_component)
						  ){

							for(int tIdx2 = 0 ; tIdx2 < nADC ; tIdx2++){

								int   TDC2_sector    = BAND_TDC.getInt  (0,tIdx2);
								int   TDC2_layer     = BAND_TDC.getInt  (1,tIdx2);
								int   TDC2_component = BAND_TDC.getInt  (2,tIdx2);
								int   TDC2_order     = BAND_TDC.getInt  (3,tIdx2);
								float TDC2_tdc       = (float)(BAND_TDC.getInt(4,tIdx2));
								float TDC2_time      = TDC2_tdc*0.02345; 

								//TDC2_time -= phaseCorr;

								if(             (!matchedBarTDC                       )&&
										(TDC1_sector   ==TDC2_sector          )&&
										(TDC1_layer    ==TDC2_layer           )&&
										(TDC1_component==TDC2_component       )&&
										(TDC1_order+1==TDC2_order             )
								  ){
									matchedBarTDC = true;

									int barKey = 100*ADC1_sector + 10*ADC1_layer + ADC1_component;

									double lnR = TMath::Log( ADC1_adc/ADC2_adc );
									double deltaT = TDC1_time-TDC2_time;
									//double x = deltaT*mean_c/100./2.;
									double x = deltaT*eff_c[barKey]/100./2.;

									h2_T_vDt[barKey] -> Fill(x,lnR);
									h2_T_vDt[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",TDC1_sector,TDC1_layer,TDC1_component));
								
									h1_ADC_L[barKey] -> Fill(ADC1_adc);	h1_ADC_L[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",TDC1_sector,TDC1_layer,TDC1_component));
									h1_ADC_R[barKey] -> Fill(ADC2_adc);	h1_ADC_R[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",TDC1_sector,TDC1_layer,TDC1_component));
								}			
							}
						}
					}		
				}
			}		
		}

	}// end file

	// -------------------------------------------------------------------------------------------------
	// Fitting functions
	TProfile ** p_T_vDt = new TProfile*[nHistos];
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_T_vDt[i]->Integral());
		if(notEmpty){
			int nBins = h2_T_vDt[i] -> GetYaxis() -> GetNbins();
                        p_T_vDt[i] = h2_T_vDt[i]->ProfileX(Form("px_%i",i), 1, nBins );

			TF1 * f_atten = new TF1("f_atten","pol1",-1,1);

			p_T_vDt[i] -> Fit("f_atten","QR");

			par[i][0] = -2/(f_atten->GetParameter(1));		//slope
			par[i][1] = 2*(f_atten->GetParError(1))/pow(par[i][0],2);	//error
		
			p_T_vDt[i] -> SetMarkerStyle(20);
                        p_T_vDt[i] -> SetMarkerColor( 1);
                        p_T_vDt[i] -> SetMarkerSize (.5);
		}
	}

	ofstream output;
	output.open("output/attenuation_lengths.txt");

	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0","c0",900,900);

	TCanvas *** cSLC = new TCanvas**[5];
	TCanvas *** cADC = new TCanvas**[5];

	TLegend * leg_adc = new TLegend(0.6,0.6,0.85,0.75);
	leg_adc -> SetLineColor(0);
	leg_adc -> AddEntry(h1_ADC_L[0],"Left PMT" );
	leg_adc -> AddEntry(h1_ADC_R[0],"Right PMT");

	for(int is = 0 ; is < 5 ; is++){
		cSLC[is] = new TCanvas*[6];
		cADC[is] = new TCanvas*[6];
		for(int il = 0 ; il < 5 ; il++){
			cSLC[is][il] = new TCanvas(Form("S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cSLC[is][il] -> Divide(2,4);

			cADC[is][il] = new TCanvas(Form("cADC_S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cADC[is][il] -> Divide(2,4);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h2_T_vDt[identifier]->Integral());
				if(notEmpty){
					int min, max;
					cSLC[is][il] -> cd(cIdx+1);
					gPad -> SetBottomMargin(0.26);
					gPad -> SetLeftMargin(0.14);
					min = getMinNonEmptyBin(h2_T_vDt[identifier]);
					max = getMaxNonEmptyBin(h2_T_vDt[identifier]);
					h2_T_vDt[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_T_vDt[identifier] -> Draw("COLZ");

					p_T_vDt[identifier] -> Draw("same");

					TLatex * tex_mu = new TLatex(-1.4,-0.1,Form("#mu = %.3f #pm %.4f m",par[identifier][0],par[identifier][1]));
					tex_mu -> SetTextSize(0.08);
					tex_mu -> Draw("same");

					output << is+1 << "\t" << il+1 << "\t" << cIdx+1 << "\t" << par[identifier][0] << "\t" << par[identifier][1] << endl;

					h1_E_minus_R_mu -> Fill(mu_efrain[identifier]/100.-par[identifier][0]);
				
					// ADC spectra
					cADC[is][il] -> cd(cIdx+1);
					gPad -> SetBottomMargin(0.26);

					h1_ADC_L[identifier] -> Draw(); 
					h1_ADC_R[identifier] -> Draw("same");
					leg_adc -> Draw("same");
				}
			}
			cSLC[is][il] -> Modified();	cSLC[is][il] -> Update();
			cADC[is][il] -> Modified();	cADC[is][il] -> Update();
		}
	}

	TCanvas * c2 = new TCanvas("c2");
	h1_E_minus_R_mu -> Draw();
	c2 -> Modified();
	c2 -> Update();

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_T_vDt[i]->Integral());
		if(!notEmpty){
			delete h2_T_vDt[i];
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c0 -> Print("results_attenuation.pdf(");
	c0 -> Print("results_adc_spectra.pdf(");

	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			cSLC[is][il] -> Print("results_attenuation.pdf");
			cADC[is][il] -> Print("results_adc_spectra.pdf");
		}
	}
	c0 -> Print("results_attenuation.pdf)");
	c0 -> Print("results_adc_spectra.pdf)");

	output.close();

	myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2) {
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetYaxis() -> CenterTitle();

	h2 -> SetTitleSize(0.5);

	h2 -> GetXaxis() -> SetLabelSize(0.09);
	h2 -> GetYaxis() -> SetLabelSize(0.09);
	h2 -> GetXaxis() -> SetTitleSize(0.09);
	h2 -> GetYaxis() -> SetTitleSize(0.09);

	h2 -> GetYaxis() -> SetTitleOffset(0.60);

	h2 -> GetYaxis() -> SetNdivisions(109);
	h2 -> GetXaxis() -> SetNdivisions(109);
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1, int color) {
        h1 -> GetXaxis() -> CenterTitle();
        h1 -> GetYaxis() -> CenterTitle();

	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);

        h1 -> SetTitleSize(0.5);

        h1 -> GetXaxis() -> SetLabelSize(0.09);
        h1 -> GetYaxis() -> SetLabelSize(0.09);
        h1 -> GetXaxis() -> SetTitleSize(0.09);
        h1 -> GetYaxis() -> SetTitleSize(0.09);

        h1 -> GetYaxis() -> SetTitleOffset(0.60);

        h1 -> GetYaxis() -> SetNdivisions(109);
        h1 -> GetXaxis() -> SetNdivisions(109);
}
// ========================================================================================================================================
double getTriggerPhase( long timeStamp ) {
	double tPh = 0.;

	long Period = 4.0;
	long Cycles = 6.0;
	long Phase  = 3.0;

	if( timeStamp != -1 ) 
		tPh = (double)(Period *( (timeStamp + Phase) % Cycles ));

	return tPh;
}
// ========================================================================================================================================
int getMinNonEmptyBin(TH2F * h2){
	int nBinX = h2 -> GetXaxis() -> GetNbins();
	int nBinY = h2 -> GetYaxis() -> GetNbins();

	int binContent = 0;

	for(int i = 1 ; i <= nBinY ; i++){
		for(int j = 1 ; j <= nBinX ; j++){
			binContent = h2 -> GetBinContent(j,i);
			if(binContent!=0) return i;
		}
	}
	return 1;
}
// ========================================================================================================================================
int getMaxNonEmptyBin(TH2F * h2){
	int nBinX = h2 -> GetXaxis() -> GetNbins();
	int nBinY = h2 -> GetYaxis() -> GetNbins();

	int binContent = 0;

	for(int i = nBinY ; i > 0 ; i--){
		for(int j = 1 ; j <= nBinX ; j++){
			binContent = h2 -> GetBinContent(j,i);
			if(binContent!=0) return i;
		}
	}
	return nBinY;
}
// ========================================================================================================================================
void LoadEffectiveVelocity(){
	ifstream f;
	int layer, sector, component, barId;
	double veff, temp;
	f.open("input/effective_velocity.txt");
	while(!f.eof()){
		f >> layer;
                f >> sector;
                f >> component;
                barId = 100*layer + 10*sector + component;
		f >> veff;
		f >> temp;
		f >> temp;
		f >> temp;
		eff_c[barId] = veff;
	}
	f.close();
}
// ========================================================================================================================================
void LoadEfrainsAttenuationLengths(){
//mu_efrain[600] = {0};
	ifstream f;
	int layer, sector, component, barId;
	double attL, temp;
	f.open("input/attenuation_lengths.txt");
	while(!f.eof()){
		f >> layer;
		f >> sector;
		f >> component;
		barId = 100*layer + 10*sector + component;
		f >> attL;
		f >> temp;
		mu_efrain[barId] = attL;
	}
	f.close();
}
// ========================================================================================================================================
