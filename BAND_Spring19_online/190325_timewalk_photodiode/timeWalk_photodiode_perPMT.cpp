#include <cstdlib>
#include <iostream>

#include "TLine.h"
#include "TProfile.h"
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
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"

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
void PrettyTGraphErrors(TGraphErrors * gP, int color);
double getTriggerPhase( long timeStamp );
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetTitleSize(0.2,"t");
	gStyle->SetOptStat(0);

	TString inputFile;

	if(argc==2) {
		inputFile = argv[1];
	}
	else {
		cout << "=========================\nRun this code as:\n./code path/to/input/file\n=========================" << endl;
		exit(0);
	}

	cout << "****************************************" << endl;
	cout << "WARNING: This code may not work unless"   << endl;
	cout << "         you run it on run 261."          << endl;
	cout << "****************************************" << endl;

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	double ref_time = 0;

	double parA_L[nHistos][2] = {{0}};
	double parB_L[nHistos][2] = {{0}};
	double parA_R[nHistos][2] = {{0}};
	double parB_R[nHistos][2] = {{0}};

	TH2F ** h2_tdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_tdc_adc_R = new TH2F * [nHistos];

	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_adc_L[i] = new TH2F(Form("h2_tdc_adc_L_%i",i),";ADC_{R};t_{TDC,L} - ref [ns]",400,100,18000,300,142,180);
		h2_tdc_adc_R[i] = new TH2F(Form("h2_tdc_adc_R_%i",i),";ADC_{R};t_{TDC,R} - ref [ns]",400,100,18000,300,142,180);
		PrettyTH2F(h2_tdc_adc_L[i]);
		PrettyTH2F(h2_tdc_adc_R[i]);
	}

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	hipo::bank BAND_ADC  ("BAND::adc"  ,reader);
	hipo::bank BAND_TDC  ("BAND::tdc"  ,reader);
	//hipo::bank RUN_config("RUN::config",reader);

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
		if(nADC<200||nTDC<200) continue;

		// Skip all laser events
		//if(nADC>100||nTDC>100) continue;

		//long timestamp = RUN_config.getLong(4,0);
		//double phaseCorr = getTriggerPhase(timestamp);

		ref_time = -100000;

		for(int tIdx1 = 0 ; tIdx1 < nTDC ; tIdx1++){
			int   TDC1_sector    = BAND_TDC.getInt  (0,tIdx1);
			int   TDC1_layer     = BAND_TDC.getInt  (1,tIdx1);
			int   TDC1_component = BAND_TDC.getInt  (2,tIdx1);
			int   TDC1_order     = BAND_TDC.getInt  (3,tIdx1);
			float TDC1_tdc       = (float)(BAND_TDC.getInt(4,tIdx1));
			float TDC1_time      = TDC1_tdc*0.02345;

			//TDC1_time -= phaseCorr;

			// Reference is PMT in bar V16-A (L), corresponding to Sector: 3, Layer: 6, Component: 6, Order
			if( TDC1_sector==3 && TDC1_layer==6 && TDC1_component==6 ){
				ref_time = TDC1_time;
			}

		}

		if(ref_time == -100000) continue;

		for(int aIdx = 0 ; aIdx < nADC ; aIdx++){

			int   ADC_sector    = BAND_ADC.getInt  (0,aIdx);
			int   ADC_layer     = BAND_ADC.getInt  (1,aIdx);
			int   ADC_component = BAND_ADC.getInt  (2,aIdx);
			int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
			float ADC_adc       = (float)(BAND_ADC.getInt(4,aIdx));
			float ADC_time      = BAND_ADC.getFloat(5,aIdx);

			bool already_matched = false;

			for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
				int   TDC_sector    = BAND_TDC.getInt  (0,tIdx);
				int   TDC_layer     = BAND_TDC.getInt  (1,tIdx);
				int   TDC_component = BAND_TDC.getInt  (2,tIdx);
				int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
				float TDC_tdc       = (float)(BAND_TDC.getInt(4,tIdx));
				float TDC_time      = TDC_tdc*0.02345;

				if(             (!already_matched            )&& // Avoid multi-hits
						(ADC_sector   ==TDC_sector   )&&
						(ADC_layer    ==TDC_layer    )&&
						(ADC_component==TDC_component)&&
						(ADC_order +2 ==TDC_order    )
				  ){
					already_matched = true;

					int barKey = 100*ADC_sector + 10*ADC_layer + ADC_component;

					//TDC_time -= phaseCorr;
					double delta_time = TDC_time - ref_time;

					// Left PMTs
					if(ADC_order==0){
						h2_tdc_adc_L[barKey] -> SetTitle(Form("Left PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
						h2_tdc_adc_L[barKey] -> Fill(ADC_adc,delta_time);
					}

					// Right PMTs
					if(ADC_order==1){
						h2_tdc_adc_R[barKey] -> SetTitle(Form("Right PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
						h2_tdc_adc_R[barKey] -> Fill(ADC_adc,delta_time);
					}

				}

			}
		}	

	}// end file

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// Fitting functions
	TProfile ** p_tdc_adc_L = new TProfile*[nHistos];
	TProfile ** p_tdc_adc_R = new TProfile*[nHistos];

	// Fitting paddle-to-paddle parameters
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_tdc_adc_L[i]->Integral());
		if(notEmpty){
			int nBins = h2_tdc_adc_L[i] -> GetYaxis() -> GetNbins();
			p_tdc_adc_L[i] = h2_tdc_adc_L[i]->ProfileX(Form("px_L_%i",i), 1, nBins );
			TF1 * f_twCorr_L = new TF1("f_twCorr_L","[0]/TMath::Sqrt(x)+[1]",100,18000);
			f_twCorr_L -> SetParameters(200,200);
			p_tdc_adc_L[i] -> Fit("f_twCorr_L","QR");

			parA_L[i][0] = f_twCorr_L->GetParameter(0);
			parA_L[i][1] = f_twCorr_L->GetParError (0);
			parB_L[i][0] = f_twCorr_L->GetParameter(1);
			parB_L[i][1] = f_twCorr_L->GetParError (1);

			p_tdc_adc_L[i] -> SetMarkerStyle(20);
			p_tdc_adc_L[i] -> SetMarkerColor( 1);
			p_tdc_adc_L[i] -> SetMarkerSize (.5);
		}

		notEmpty = (h2_tdc_adc_R[i]->Integral());
		if(notEmpty){
			int nBins = h2_tdc_adc_R[i] -> GetYaxis() -> GetNbins();
			p_tdc_adc_R[i] = h2_tdc_adc_R[i]->ProfileX(Form("px_R_%i",i), 1, nBins );
			TF1 * f_twCorr_R = new TF1("f_twCorr_R","[0]/TMath::Sqrt(x)+[1]",100,18000);
			f_twCorr_R -> SetParameters(200,200);
			p_tdc_adc_R[i] -> Fit("f_twCorr_R","QR");

			parA_R[i][0] = f_twCorr_R->GetParameter(0);
			parA_R[i][1] = f_twCorr_R->GetParError (0);
			parB_R[i][0] = f_twCorr_R->GetParameter(1);
			parB_R[i][1] = f_twCorr_R->GetParError (1);

			p_tdc_adc_R[i] -> SetMarkerStyle(20);
			p_tdc_adc_R[i] -> SetMarkerColor( 1);
			p_tdc_adc_R[i] -> SetMarkerSize (.5);
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0","c0",700,900);
	TLatex * tex_f = new TLatex(0.1,0.5,"y = A/#sqrt{x} + B");
	tex_f -> Draw();
	c0 -> Modified();
	c0 -> Update(); 

	TCanvas *** cSLC = new TCanvas**[5];

	for(int is = 0 ; is < 5 ; is++){
		cSLC[is] = new TCanvas*[6];
		for(int il = 0 ; il < 6 ; il++){
			cSLC[is][il] = new TCanvas(Form("S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),700,900);
			cSLC[is][il] -> Divide(2,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h2_tdc_adc_L[identifier]->Integral());
				if(notEmpty){
					int min, max;
					double min_val, max_val;

					cSLC[is][il] -> cd(2*cIdx+1);
					gPad -> SetLeftMargin(0.16);
					gPad -> SetBottomMargin(0.20);
					gPad -> SetTopMargin(0);
					min = getMinNonEmptyBin(h2_tdc_adc_L[identifier]);
					max = getMaxNonEmptyBin(h2_tdc_adc_L[identifier]);
					min_val = h2_tdc_adc_L[identifier] -> GetYaxis() -> GetBinCenter(min);
					max_val = h2_tdc_adc_L[identifier] -> GetYaxis() -> GetBinCenter(max);
					h2_tdc_adc_L[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_tdc_adc_L[identifier] -> Draw("COLZ");
					p_tdc_adc_L[identifier] -> Draw("same");

					TLatex * tex_AL = new TLatex(10500,max_val-0.30*(max_val-min_val),Form("A = %.3f #pm %.3f",parA_L[identifier][0],parA_L[identifier][1]));
					tex_AL -> SetTextSize(0.08);
					tex_AL -> Draw("same");

					TLatex * tex_BL = new TLatex(10500,max_val-0.42*(max_val-min_val),Form("B = %.3f #pm %.3f",parB_L[identifier][0],parB_L[identifier][1]));
					tex_BL -> SetTextSize(0.08);
					tex_BL -> Draw("same");
				}
				notEmpty = (h2_tdc_adc_R[identifier]->Integral());
				if(notEmpty){
					int min, max;
                                        double min_val, max_val;

					cSLC[is][il] -> cd(2*cIdx+2);
					gPad -> SetLeftMargin(0.16);
					gPad -> SetBottomMargin(0.20);
					gPad -> SetTopMargin(0);
					min = getMinNonEmptyBin(h2_tdc_adc_R[identifier]);
					max = getMaxNonEmptyBin(h2_tdc_adc_R[identifier]);
					min_val = h2_tdc_adc_R[identifier] -> GetYaxis() -> GetBinCenter(min);
					max_val = h2_tdc_adc_R[identifier] -> GetYaxis() -> GetBinCenter(max);
					h2_tdc_adc_R[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_tdc_adc_R[identifier] -> Draw("COLZ");
					p_tdc_adc_R[identifier] -> Draw("same");

					TLatex * tex_AR = new TLatex(10500,max_val-0.30*(max_val-min_val),Form("A = %.3f #pm %.3f",parA_R[identifier][0],parA_R[identifier][1]));
					tex_AR -> SetTextSize(0.08);
					tex_AR -> Draw("same");

					TLatex * tex_BR = new TLatex(10500,max_val-0.42*(max_val-min_val),Form("B = %.3f #pm %.3f",parB_R[identifier][0],parB_R[identifier][1]));
					tex_BR -> SetTextSize(0.08);
					tex_BR -> Draw("same");
				}
			}
			cSLC[is][il] -> Modified();     cSLC[is][il] -> Update();
		}
	}
	// -------------------------------------------------------------------------------------------------
	// Saving fit values to ccdb tables
	ofstream tabL, tabR;
	tabL.open("timeWalkPar_L.txt");
	tabR.open("timeWalkPar_R.txt");

	for(int is = 1 ; is <= 5 ; is++){
		for(int il = 1 ; il <= 6 ; il++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;
				tabL << is << "\t" << il << "\t" << ic << "\t" << "\t";
				//if(il==6) tabL << "0\t0\t0\t0" << endl;
				//else
				tabL << parA_L[idx][0] << "\t" << parB_L[idx][0] << "\t" << parA_L[idx][1]  << "\t" << parB_L[idx][1] << endl;
				// ---
				tabR << is << "\t" << il << "\t" << ic << "\t" << "\t";
				//if(il==6) tabR << "0\t0\t0\t0" << endl;
				//else
				tabR << parA_R[idx][0] << "\t" << parB_R[idx][0] << "\t" << parA_R[idx][1]  << "\t" << parB_R[idx][1] << endl;
			}
		}
	}
	tabL.close();
	tabR.close();

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_tdc_adc_L[i]->Integral()+h2_tdc_adc_R[i]->Integral());
		if(!notEmpty){
			delete h2_tdc_adc_L[i];
			delete h2_tdc_adc_R[i];
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c0 -> Print("results_timewalk_corr_perPMT.pdf(");

	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 6 ; il++){
			cSLC[is][il] -> Print("results_timewalk_corr_perPMT.pdf");
		}
	}
	c0 -> Print("results_timewalk_corr_perPMT.pdf)");

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

	h2 -> GetYaxis() -> SetTitleOffset(0.80);

	h2 -> GetYaxis() -> SetNdivisions(109);
	h2 -> GetXaxis() -> SetNdivisions(107);
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
	h1 -> GetXaxis() -> SetNdivisions(107);
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
void PrettyTGraphErrors(TGraphErrors * gP, int color){
	gP -> SetTitle("");
	gP -> SetMarkerColor(color);
	gP -> SetMarkerStyle(20);
	gP -> GetXaxis() -> CenterTitle();
	gP -> GetXaxis() -> SetTitle("Bar ID");
	gP -> GetYaxis() -> CenterTitle();
	gP -> GetYaxis() -> SetTitle("resolution/#sqrt{2} [ps]");

	gP -> GetYaxis() -> SetTitleOffset(1.30);

	gP -> SetMinimum(150);
	gP -> SetMaximum(400);
}
