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

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	double ref_time = 0;

	double parA_L[nHistos][2] = {{0}};
	double parB_L[nHistos][2] = {{0}};

	TH2F ** h2_dMeantime_adc = new TH2F * [nHistos];

	for(int i = 0 ; i < nHistos ; i++){
		h2_dMeantime_adc[i] = new TH2F(Form("h2_dMeantime_adc_%i",i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC} - ref",500,0,20000,800,100,200);
		PrettyTH2F(h2_dMeantime_adc[i]);
	}

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
		if(nADC<200||nTDC<200) continue;

		// Skip all laser events
		//if(nADC>100||nTDC>100) continue;

		long timestamp = RUN_config.getLong(4,0);
		double phaseCorr = getTriggerPhase(timestamp);

		ref_time = -100000;

		for(int tIdx1 = 0 ; tIdx1 < nTDC ; tIdx1++){
			int   TDC1_sector    = BAND_TDC.getInt  (0,tIdx1);
			int   TDC1_layer     = BAND_TDC.getInt  (1,tIdx1);
			int   TDC1_component = BAND_TDC.getInt  (2,tIdx1);
			int   TDC1_order     = BAND_TDC.getInt  (3,tIdx1);
			float TDC1_tdc       = (float)(BAND_TDC.getInt(4,tIdx1));
			float TDC1_time      = TDC1_tdc*0.02345;

			TDC1_time -= phaseCorr;

			// Reference is PMT in bar V16-A (L), corresponding to Sector: 3, Layer: 6, Component: 6, Order
			if( TDC1_sector==3 && TDC1_layer==6 && TDC1_component==6 ){
				ref_time = TDC1_time;
			}

		}

		if(ref_time == -100000) continue;

		// Looping over entries of the event for a second time
		for(int aIdx1 = 0 ; aIdx1 < nADC ; aIdx1++){

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
				if(             (ADC1_sector   ==ADC2_sector          )&&
						(ADC1_layer    ==ADC2_layer           )&&
						(ADC1_component==ADC2_component       )&&
						(ADC1_order+1==ADC2_order             )
				  ){

					bool already_matched = false;
					for(int tIdx1 = 0 ; tIdx1 < nTDC ; tIdx1++){

						int   TDC1_sector    = BAND_TDC.getInt  (0,tIdx1);
						int   TDC1_layer     = BAND_TDC.getInt  (1,tIdx1);
						int   TDC1_component = BAND_TDC.getInt  (2,tIdx1);
						int   TDC1_order     = BAND_TDC.getInt  (3,tIdx1);
						float TDC1_tdc       = (float)(BAND_TDC.getInt(4,tIdx1));
						float TDC1_time      = TDC1_tdc*0.02345;

						TDC1_time -= phaseCorr;

						// Matching these ADCs to a TDC
						if(             (ADC1_sector   ==TDC1_sector   )&&
								(ADC1_layer    ==TDC1_layer    )&&
								(ADC1_component==TDC1_component)&&
								(ADC1_order+2  ==TDC1_order    )
						  ){

							for(int tIdx2 = 0 ; tIdx2 < nTDC ; tIdx2++){

								int   TDC2_sector    = BAND_TDC.getInt  (0,tIdx2);
								int   TDC2_layer     = BAND_TDC.getInt  (1,tIdx2);
								int   TDC2_component = BAND_TDC.getInt  (2,tIdx2);
								int   TDC2_order     = BAND_TDC.getInt  (3,tIdx2);
								float TDC2_tdc       = (float)(BAND_TDC.getInt(4,tIdx2));
								float TDC2_time      = TDC2_tdc*0.02345;

								TDC2_time -= phaseCorr;

								if(             (!already_matched                     )&& // Avoid multi-hits
										(TDC1_sector   ==TDC2_sector          )&&
										(TDC1_layer    ==TDC2_layer           )&&
										(TDC1_component==TDC2_component       )&&
										(TDC1_order+1==TDC2_order             )
								  ){	

									already_matched = true;

									int barKey = 100*ADC1_sector + 10*ADC1_layer + ADC1_component;			

									double meantime = ( TDC1_time + TDC2_time )/2.;
									double delta_meantime = meantime - ref_time;
									double adc_geometric_mean = TMath::Sqrt(ADC1_adc*ADC2_adc);
									
									h2_dMeantime_adc[barKey] -> Fill(adc_geometric_mean,delta_meantime);
									h2_dMeantime_adc[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",TDC1_sector,TDC1_layer,TDC1_component));
								}
							}
						}
					}
				}
			}
		}

	}// end file


	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// Fitting functions
	TProfile ** p_dMeantime_adc = new TProfile*[nHistos];
	// Fitting paddle-to-paddle parameters
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_dMeantime_adc[i]->Integral());
		if(notEmpty){
			int nBins = h2_dMeantime_adc[i] -> GetYaxis() -> GetNbins();
			p_dMeantime_adc[i] = h2_dMeantime_adc[i]->ProfileX(Form("px_%i",i), 1, nBins );
			TF1 * f_twCorr = new TF1("f_twCorr","[0]/TMath::Sqrt(x)+[1]",2000,20000);
			f_twCorr -> SetParameters(200,200);
			f_twCorr -> SetParLimits(0,100,300);
			p_dMeantime_adc[i] -> Fit("f_twCorr","QR");

			parA_L[i][0] = f_twCorr->GetParameter(0);
			parA_L[i][1] = f_twCorr->GetParError (0);
			parB_L[i][0] = f_twCorr->GetParameter(1);
                        parB_L[i][1] = f_twCorr->GetParError (1);

			p_dMeantime_adc[i] -> SetMarkerStyle(20);
			p_dMeantime_adc[i] -> SetMarkerColor( 1);
			p_dMeantime_adc[i] -> SetMarkerSize (.5);
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0","c0",900,900);

	TCanvas *** c_tw = new TCanvas**[5];

	for(int is = 0 ; is < 5 ; is++){
		c_tw[is] = new TCanvas*[6];


		for(int il = 0 ; il < 5 ; il++){
			c_tw[is][il] = new TCanvas(Form("S%iL%i_paddle",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			c_tw[is][il] -> Divide(2,4);	

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h2_dMeantime_adc[identifier]->Integral());
				if(notEmpty){
					int min, max;
					double min_val, max_val;
					c_tw[is][il] -> cd(cIdx+1);
					gPad -> SetBottomMargin(0.26);
					gPad -> SetLeftMargin(0.14);
					
					min = getMinNonEmptyBin(h2_dMeantime_adc[identifier]);
					max = getMaxNonEmptyBin(h2_dMeantime_adc[identifier]);
					min_val = h2_dMeantime_adc[identifier] -> GetYaxis() -> GetBinCenter(min);
					max_val = h2_dMeantime_adc[identifier] -> GetYaxis() -> GetBinCenter(max);

					h2_dMeantime_adc[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_dMeantime_adc[identifier] -> Draw("COLZ");
					p_dMeantime_adc[identifier] -> Draw("same");		

					TLatex * tex_f = new TLatex(500,min_val+0.30*(max_val-min_val),"y = A/#sqrt{x} + B");
					tex_f -> SetTextSize(0.08);
                                        tex_f -> Draw("same");

					TLatex * tex_A = new TLatex(500,min_val+0.18*(max_val-min_val),Form("A = %.4f #pm %.4f",parA_L[identifier][0],parA_L[identifier][1]));
					tex_A -> SetTextSize(0.08);
					tex_A -> Draw("same");

					TLatex * tex_B = new TLatex(500,min_val+0.06*(max_val-min_val),Form("B = %.4f #pm %.4f",parB_L[identifier][0],parB_L[identifier][1]));
                                        tex_B -> SetTextSize(0.08);
                                        tex_B -> Draw("same");
				} 
			}
			c_tw[is][il] -> Modified();	c_tw[is][il] -> Update();		
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
				if(il==6) tabL << "0\t0\t0\t0" << endl;
				else tabL << parA_L[idx][0] << "\t" << parB_L[idx][0] << "\t" << parA_L[idx][1] << "\t" << parB_L[idx][1] << endl;
			}
		}
	}
	tabL.close();
	tabR.close();

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_dMeantime_adc[i]->Integral());
		if(!notEmpty){
			delete h2_dMeantime_adc[i];
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c0 -> Print("results_timewalk_corr.pdf(");

	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			c_tw[is][il] -> Print("results_timewalk_corr.pdf");
		}
	}
	c0 -> Print("results_timewalk_corr.pdf)");

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
