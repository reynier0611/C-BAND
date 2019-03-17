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

#include "reader.h"
#include "node.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

double p2p[600] = {0};
double l2l[600] = {0};

// Forward-declaring functions
void LoadPaddleCorrectionPar();
double TDC_corr( int barId );
void PrettyTH1F(TH1F * h1, int color);
void PrettyTH2F(TH2F * h2);
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
	// Load TDC paddle-to-paddle correction parameters only if data is not already corrected
	cout << "Is your data already TDC paddle offset corrected?\n1) no, 2) yes, 3) I don't know?" << endl;
	int option;
	cin >> option;
	if(option==3){
		cout << "Bitch please! Get your shit together and then come back. I cannot answer this for you." << endl;
		exit(0);
	}
	if(option==1)
		LoadPaddleCorrectionPar();

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	TH2F ** h2_dMeantime_adc_paddle = new TH2F * [nHistos];

	TH2F * h2_dMeantime_adc_before = new TH2F("h2_dMeantime_adc_before","Before;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,15000,500,660,720);
	TH2F * h2_dMeantime_adc_after  = new TH2F("h2_dMeantime_adc_after" ,"After ;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,15000,500,660,720);

	for(int i = 0 ; i < nHistos ; i++){
		h2_dMeantime_adc_paddle[i] = new TH2F(Form("h2_dMeantime_adc_paddle_%i",i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,15000,500,660,720);		
		PrettyTH2F(h2_dMeantime_adc_paddle[i]);	
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
									double meantime_corr = meantime + TDC_corr(barKey);	
									double adc_geometric_mean = TMath::Sqrt(ADC1_adc*ADC2_adc);

									h2_dMeantime_adc_paddle[barKey] -> Fill(adc_geometric_mean,meantime_corr);
									h2_dMeantime_adc_paddle[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",TDC1_sector,TDC1_layer,TDC1_component));


									h2_dMeantime_adc_before -> Fill(adc_geometric_mean,meantime     );
									h2_dMeantime_adc_after  -> Fill(adc_geometric_mean,meantime_corr);

								}
							}
						}
					}
				}
			}
		}


	}// end file
	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0","c0",900,900);

	TCanvas *** cSLC_paddle = new TCanvas**[5];

	for(int is = 0 ; is < 5 ; is++){
		cSLC_paddle[is] = new TCanvas*[6];

		for(int il = 0 ; il < 5 ; il++){
			cSLC_paddle[is][il] = new TCanvas(Form("S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cSLC_paddle[is][il] -> Divide(2,4);	

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h2_dMeantime_adc_paddle[identifier]->Integral());
				if(notEmpty){
					int min, max;
					cSLC_paddle[is][il] -> cd(cIdx+1);
					gPad -> SetBottomMargin(0.26);
					gPad -> SetLeftMargin(0.14);
					//min = getMinNonEmptyBin(h2_dMeantime_adc_paddle[identifier]);
					//max = getMaxNonEmptyBin(h2_dMeantime_adc_paddle[identifier]);
					//h2_dMeantime_adc_paddle[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_dMeantime_adc_paddle[identifier] -> Draw("COLZ");
				} 
			}
			cSLC_paddle[is][il] -> Modified();	cSLC_paddle[is][il] -> Update();		
		}
	}

	TCanvas * c1 = new TCanvas("c1","c1",900,900);
	c1 -> Divide(1,2);
	c1 -> cd(1);	h2_dMeantime_adc_before -> Draw("COLZ");
	c1 -> cd(2);	h2_dMeantime_adc_after  -> Draw("COLZ");
	c1 -> Modified();
	c1 -> Update();

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_dMeantime_adc_paddle[i]->Integral());
		if(!notEmpty){
			delete h2_dMeantime_adc_paddle[i];
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c1 -> Print("check_parameters_extracted.pdf(");

	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			cSLC_paddle[is][il] -> Print("check_parameters_extracted.pdf");
		}
	}
	c0 -> Print("check_parameters_extracted.pdf)");

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

	h2 -> GetYaxis() -> SetTitleOffset(0.70);

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
void LoadPaddleCorrectionPar(){
	ifstream f;
	int layer, sector, component, barId;
	double parameter, temp;

	f.open("TDC_paddle_offsets.txt");
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

	f.open("TDC_layer_offsets.txt");
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
double TDC_corr( int barId ){
	return -( p2p[barId] + l2l[barId] );
}
