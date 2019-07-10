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
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

double p2p_tdc[600] = {0};
double l2l_tdc[600] = {0};
double p2p_fadc[600] = {0};
double l2l_fadc[600] = {0};

// Forward-declaring functions
void LoadPaddleCorrectionPar();
void PrettyTH1F(TH1F * h1, int color);
void PrettyTH2F(TH2F * h2);
void PrettyTGraphErrors(TGraphErrors * gP, int color);
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};

// ========================================================================================================================================
int main(int argc, char** argv) {

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

	double ref_meanTime_tdc[5] 	= {0.};
	double ref_meanTime_fadc[5] 	= {0.};

	TH2F * h2_dMeantime_tdc_before = new TH2F("h2_dMeantime_tdc_before","TDC Before;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,25000,600,-3,3);
	TH2F * h2_dMeantime_tdc_after  = new TH2F("h2_dMeantime_tdc_after" ,"TDC After ;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}",500,0,25000,600,-3,3);
	TH2F * h2_dMeantime_fadc_before = new TH2F("h2_dMeantime_fadc_before","FADC Before;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{FADC}",500,0,25000,600,-3,3);
	TH2F * h2_dMeantime_fadc_after  = new TH2F("h2_dMeantime_fadc_after" ,"FADC After ;#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{FADC}",500,0,25000,600,-3,3);


	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);
	
	//Read Dictionary of Hipo File  // new hipo4
        hipo::dictionary  factory;      // new hipo4
        reader.readDictionary(factory); // new hipo4
        //factory.show();               // new hipo4

	BBand         band_hits   (factory.getSchema("BAND::hits" ));

	//One also needs a hipo::event object which is called from the reader for each event to get
        //the information for each bank
        hipo::event readevent;  // new hipo4

	int event_counter = 0;
	// Loop over events
	while(reader.next()==true){
		//Reader has to load information about event in hipo::event class
                reader.read(readevent); // new hipo4

                //Load explicitly all information for each bank for the event
                readevent.getStructure(band_hits  );   // new hipo4

                //Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
                //The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.

		if(event_counter%1000==0) cout << "event: " << event_counter << endl;
		event_counter++;	

		// Grab number of hits and skip if not our optimal size
		int nHits = band_hits.getRows();
		if(nHits==0) continue;
		if(nHits<120) continue;

		// Clear out the reference times for each new event
		fill( begin(ref_meanTime_tdc), 	end(ref_meanTime_tdc), 	0. );
		fill( begin(ref_meanTime_fadc), end(ref_meanTime_fadc), 0. );


		// First loop through all the hits and get the reference time for each layer
		for(int hit = 0; hit < nHits; hit++) {
			int    sector		= band_hits.getSector     (hit);
                        int    layer		= band_hits.getLayer      (hit);
                        int    component	= band_hits.getComponent  (hit);
		
			// Reference bar in each layer is sector 2, component 1
			if( sector == 2 && component == 1 ){
				float  adcLcorr		= band_hits.getAdcLcorr   (hit);
				float  adcRcorr		= band_hits.getAdcRcorr   (hit);
				if(adcLcorr< 6000||adcRcorr< 6000) continue; 
				if(adcLcorr>20000||adcRcorr>20000) continue; 

				// If all passed, grab the reference times for this event. 
				// This reference time has an offset already for L-R and TW implicitly
				double meantimeTdc	= band_hits.getMeantimeTdc(hit);
				double meantimeFadc	= band_hits.getMeantimeFadc(hit);

				ref_meanTime_tdc[layer-1]	= meantimeTdc;
				ref_meanTime_fadc[layer-1]	= meantimeFadc;
			}

		}

		// Looping over all bars now that we have our reference time
		for(int hit = 0; hit < nHits; hit++) {
			// Grab all the info for this bar hit
			int    sector            = band_hits.getSector       (hit);
                        int    layer             = band_hits.getLayer        (hit);
                        int    component         = band_hits.getComponent    (hit);
                        int    barKey            = band_hits.getBarKey       (hit);
                        float adcLcorr           = band_hits.getAdcLcorr     (hit);                        
                        float adcRcorr           = band_hits.getAdcRcorr     (hit);
                        double  meantimeTdc      = band_hits.getMeantimeTdc  (hit);
			double  meantimeFadc	 = band_hits.getMeantimeFadc (hit);

			// If we don't have any reference time, we can't use this event
			if( ref_meanTime_tdc[layer-1] == 0. ) continue;
			if( ref_meanTime_fadc[layer-1] == 0.) continue;

			// Skip for bad ADCs
			if(adcLcorr< 6000||adcRcorr< 6000) continue; 
			if(adcLcorr>20000||adcRcorr>20000) continue; 

			// Get this time minus our reference time to do P2P offsets
			double delta_meantime_tdc 	= meantimeTdc  - ref_meanTime_tdc[5-1];
			double delta_meantime_fadc	= meantimeFadc - ref_meanTime_fadc[5-1];
			double adc_geometric_mean 	= sqrt(adcLcorr*adcRcorr);

			if( (sector == 2 && component == 1) ){
			}
			else{
				// Now store all this in histograms:
				h2_dMeantime_tdc_after -> Fill(adc_geometric_mean, delta_meantime_tdc - ( p2p_tdc[barKey] + l2l_tdc[barKey] )  );
				h2_dMeantime_tdc_before -> Fill( adc_geometric_mean, delta_meantime_tdc );
				h2_dMeantime_fadc_after -> Fill(adc_geometric_mean, delta_meantime_fadc - ( p2p_fadc[barKey] + l2l_fadc[barKey] )  );
				h2_dMeantime_fadc_before -> Fill( adc_geometric_mean, delta_meantime_fadc );
			}
		}
	}// end file


	TCanvas * c1 = new TCanvas("c1","c1",900,900);
	c1 -> Divide(1,2);
	c1 -> cd(1);    h2_dMeantime_tdc_before -> Draw("COLZ");
	c1 -> cd(2);    h2_dMeantime_tdc_after  -> Draw("COLZ");
	c1 -> Modified();
	c1 -> Update();
	TCanvas * c2 = new TCanvas("c2","c2",900,900);
	c2 -> Divide(1,2);
	c2 -> cd(1);    h2_dMeantime_fadc_before -> Draw("COLZ");
	c2 -> cd(2);    h2_dMeantime_fadc_after  -> Draw("COLZ");
	c2 -> Modified();
	c2 -> Update();

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c1 -> Print("check_parameters_extracted.pdf(");
	c2 -> Print("check_parameters_extracted.pdf)");


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
		p2p_tdc[barId] = parameter;
	}
	f.close();
	f.open("FADC_paddle_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parameter;
		f >> temp;
		p2p_fadc[barId] = parameter;
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
		l2l_tdc[barId] = parameter;
	}
	f.close();
	f.open("FADC_layer_offsets.txt");
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> parameter;
		f >> temp;
		l2l_fadc[barId] = parameter;
	}
	f.close();
}
// ========================================================================================================================================

