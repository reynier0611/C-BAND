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

double p2p[600] = {0};
double l2l[600] = {0};

// Forward-declaring functions
void LoadPaddleCorrectionPar();
double TDC_corr( int barId );
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

	cout << "*************************************************************" << endl;
	cout << "                  ___                  ___" << endl;
	cout << "\\    /\\    / /\\   |  \\ |\\  | | |\\  |  /" << endl;
	cout << " \\  /  \\  / /--\\  |-\\/ | \\ | | | \\ | |  --\\" << endl;
	cout << "  \\/    \\/ /    \\ |  \\ |  \\| | |  \\|  \\___/" << endl << endl;
	cout << "Run this code on laser data on a file that does not have" << endl;
	cout << "paddle-to-paddle corrections yet, but already has time-walk" << endl;
	cout << "and L-R corrections" << endl << endl;
	cout << "*************************************************************" << endl;

	LoadPaddleCorrectionPar();

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	double ref_meantime[10] = {0};

	TH2F ** h2_dMeantime_adc_paddle = new TH2F * [nHistos];

	TH2F * h2_dMeantime_adc_before = new TH2F("h2_dMeantime_adc_before","Before;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,25000,500,-15,15);
	TH2F * h2_dMeantime_adc_after  = new TH2F("h2_dMeantime_adc_after" ,"After ;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,25000,500,-15,15);

	for(int i = 0 ; i < nHistos ; i++){
		h2_dMeantime_adc_paddle[i] = new TH2F(Form("h2_dMeantime_adc_paddle_%i",i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,15000,500,-15,15);
		PrettyTH2F(h2_dMeantime_adc_paddle[i]);
	}

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	//hipo::bank BAND_ADC  ("BAND::adc"  ,reader);
	//hipo::bank BAND_TDC  ("BAND::tdc"  ,reader);
	//hipo::bank RUN_config("RUN::config",reader);
	BBand         band_hits   ("BAND::hits"       ,reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000==0) cout << "event: " << event_counter << endl;
		event_counter++;	

		//BAND_ADC.show();
		//BAND_TDC.show();
		//RUN_config.show();

		int nHits = band_hits.getSize();

		// Skip events with no entries
		if(nHits==0) continue;

		// Skip events with less than 200 entries (the laser lights all bars at the same time)
		if(nHits<50) continue;

		// Skip all laser events
		//if(nHits>100) continue;

		for(int hit = 0; hit < nHits; hit++) {
			int    sector            = band_hits.getSector    (hit);
			int    layer             = band_hits.getLayer     (hit);
			int    component         = band_hits.getComponent (hit);
			float  tTdcLcorr         = band_hits.getTtdcLcorr (hit);
			float  tTdcRcorr         = band_hits.getTtdcRcorr (hit);
			float  meantimeTdc       = band_hits.getMeantimeTdc (hit);

			if( (sector==2)&&(component==1) ){
				ref_meantime[layer] = meantimeTdc;
			}
		}

		// Looping over entries of the event for a second time

		for(int hit = 0; hit < nHits; hit++) {

			int    sector            = band_hits.getSector      (hit);
			int    layer             = band_hits.getLayer       (hit);
			int    component         = band_hits.getComponent   (hit);
			int    barKey            = band_hits.getBarKey      (hit);
			float  meantimeTdc       = band_hits.getMeantimeTdc (hit);
			float adcLcorr           = band_hits.getAdcLcorr    (hit);                        
			float adcRcorr           = band_hits.getAdcRcorr    (hit);
			float tTdcLcorr          = band_hits.getTtdcLcorr   (hit);
			float tTdcRcorr          = band_hits.getTtdcRcorr   (hit);

			if(adcLcorr>6000&&adcRcorr>6000){

				double meantime_corr = meantimeTdc + TDC_corr(barKey);
				double adc_geometric_mean = TMath::Sqrt(adcLcorr*adcRcorr);	

				h2_dMeantime_adc_paddle[barKey] -> Fill(adc_geometric_mean,meantime_corr-ref_meantime[layer]);
				h2_dMeantime_adc_paddle[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",sector,layer,component));

				h2_dMeantime_adc_before -> Fill(adc_geometric_mean,meantimeTdc  -ref_meantime[layer]);
				h2_dMeantime_adc_after  -> Fill(adc_geometric_mean,meantime_corr-ref_meantime[layer]);
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
			cSLC_paddle[is][il] -> Modified();      cSLC_paddle[is][il] -> Update();
		}
	}

	TCanvas * c1 = new TCanvas("c1","c1",900,900);
	c1 -> Divide(1,2);
	c1 -> cd(1);    h2_dMeantime_adc_before -> Draw("COLZ");
	c1 -> cd(2);    h2_dMeantime_adc_after  -> Draw("COLZ");
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
void PrettyTGraphErrors(TGraphErrors * gP, int color){
	gP -> SetTitle("");
	gP -> SetMarkerColor(color);
	gP -> SetMarkerStyle(20);
	gP -> GetXaxis() -> CenterTitle();
	gP -> GetXaxis() -> SetTitle("Bar ID");
	gP -> GetYaxis() -> CenterTitle();
	gP -> GetYaxis() -> SetTitle("resolution/#sqrt{2} [ps]");

	gP -> GetYaxis() -> SetTitleOffset(1.30);

	gP -> SetMinimum( 80);
	gP -> SetMaximum(200);
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

