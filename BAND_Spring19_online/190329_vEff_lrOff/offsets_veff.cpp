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

double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};
// ========================================================================================================================================
int main(int argc, char** argv) {

        TApplication *myapp = new TApplication("myapp",0,0);

	std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

	TString inputFile;

	if(argc==2) {
		inputFile = argv[1];
	}
	else {
		cout << "=========================\nRun this code as:\n./code path/to/input/file\n=========================" << endl;
		exit(0);
	}

	// ----------------------------------------------------------------------------------
	// Useful variables
	double mp      = 0.93827; //GeV
	double mD      = 1.8756;  //GeV (deuteron mass    )
	double rad2deg = 180./3.14159;
	double c = 29.9792;

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	TH1F ** h1_tdc_diff = new TH1F * [nHistos];
	TH1F ** h1_ftdc_diff = new TH1F * [nHistos];
	for( int i = 0 ; i < nHistos ; i++ ){
		h1_tdc_diff[i] = new TH1F(Form("h1_tdc_diff_%i",i),"",1000,-50,50);
		h1_ftdc_diff[i] = new TH1F(Form("h1_ftdc_diff_%i",i),"",1000,-50,50);
	}

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	BBand         band_hits   ("BAND::hits"       ,reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;

		int nHits = band_hits.getSize();
		if( nHits > 10 ) continue;
		//cout << "**" << nHits <<"\n";
		for(int hit = 0; hit < nHits; hit++) {

			int    sector            = band_hits.getSector      (hit);
			int    layer             = band_hits.getLayer       (hit);
			int    component         = band_hits.getComponent   (hit);
			int    barKey            = band_hits.getBarKey      (hit); 
			float  meantimeTdc       = band_hits.getMeantimeTdc (hit);
			float  meantimeFadc      = band_hits.getMeantimeFadc(hit);
			float difftimeTdc        = band_hits.getDifftimeTdc (hit);
			float difftimeFadc       = band_hits.getDifftimeFadc(hit);
			float adcLcorr           = band_hits.getAdcLcorr    (hit);
			float adcRcorr           = band_hits.getAdcRcorr    (hit);

			float tFadcLcorr = band_hits.getTfadcLcorr(hit);
			float tFadcRcorr = band_hits.getTfadcRcorr(hit);
			float tTdcLcorr  = band_hits.getTtdcLcorr (hit);
			float tTdcRcorr  = band_hits.getTtdcRcorr (hit);
			float x          = band_hits.getX         (hit);
			float y          = band_hits.getY         (hit);
			float z          = band_hits.getZ         (hit);
			float ux         = band_hits.getUx        (hit);
			float uy         = band_hits.getUy        (hit);
			float uz         = band_hits.getUz        (hit);
			//cout << barKey << " " << tTdcLcorr - tFadcLcorr << " " << adcLcorr << "\n";
			h1_tdc_diff[barKey]->Fill( difftimeTdc );
			h1_ftdc_diff[barKey]->Fill( difftimeFadc );

		}// end loop over hits in event
	}// end file

	// Get the width for each histogram:
	ofstream lr_offsets,effective_velocity;
	lr_offsets.open("lr_offsets.txt");
	effective_velocity.open("effective_velocity.txt");
	for( int i = 0 ; i < nHistos ; i++ ){
		if( h1_tdc_diff[i]->Integral() == 0 ) continue;

		int sector = (i/100);
		int layer = (i - sector*100)/10;
		int component = (i - sector*100 - layer*10);
		double barLength = bandlen[sector-1];

		// For FADC
		int ftdc_low_Bin = h1_ftdc_diff[i]->FindFirstBinAbove(3);
		int ftdc_high_Bin = h1_ftdc_diff[i]->FindLastBinAbove(3);
		double ftdc_left = h1_ftdc_diff[i]->GetXaxis()->GetBinCenter(ftdc_low_Bin);
		double ftdc_right = h1_ftdc_diff[i]->GetXaxis()->GetBinCenter(ftdc_high_Bin);
		double ftdc_width = ftdc_right - ftdc_left;
		double ftdc_vEff = (2*barLength) / ftdc_width;
		double ftdc_lr_off = (ftdc_left+ftdc_right)/2.;

		// For TDC
		int tdc_low_Bin = h1_tdc_diff[i]->FindFirstBinAbove(3);
		int tdc_high_Bin = h1_tdc_diff[i]->FindLastBinAbove(3);
		double tdc_left = h1_tdc_diff[i]->GetXaxis()->GetBinCenter(tdc_low_Bin);
		double tdc_right = h1_tdc_diff[i]->GetXaxis()->GetBinCenter(tdc_high_Bin);
		double tdc_width = tdc_right - tdc_left;
		double tdc_vEff = (2*barLength) / tdc_width;
		double tdc_lr_off = (tdc_left+tdc_right)/2.;

		if( isinf(tdc_vEff) ) tdc_vEff = 0.;
		if( isinf(ftdc_vEff) ) ftdc_vEff = 0.;
		
		effective_velocity << sector << "\t" << layer << "\t" << component << "\t" << tdc_vEff << "\t" << ftdc_vEff << "\t" << 0.000 << "\t" << 0.000 << "\n";
		lr_offsets << sector << "\t" << layer << "\t" << component << "\t" << tdc_lr_off << "\t" << ftdc_lr_off << "\t" << 0.000 << "\t" << 0.000 << "\n";
	}
	lr_offsets.close();
	effective_velocity.close();

	//myapp -> Run();
	return 0;
}
// ========================================================================================================================================
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color) {
	h1 -> GetXaxis() -> SetTitle(titx);
	h1 -> GetYaxis() -> SetTitle(tity);
	h1 -> SetLineColor(color);
	h1 -> SetLineWidth(2);

	h1 -> GetXaxis() -> CenterTitle();
        h1 -> GetYaxis() -> CenterTitle();

	h1 -> GetXaxis() -> SetTitleSize(0.07);
	h1 -> GetXaxis() -> SetLabelSize(0.07);
	h1 -> GetXaxis() -> SetNdivisions(107);

	h1 -> GetYaxis() -> SetTitleSize(0.07);
        h1 -> GetYaxis() -> SetLabelSize(0.07);
        h1 -> GetYaxis() -> SetNdivisions(107);
	h1 -> GetYaxis() -> SetTitleOffset(0.63);
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2,TString titx,TString tity) {
	h2 -> GetXaxis() -> SetTitle(titx);
	h2 -> GetYaxis() -> SetTitle(tity);
}
// ========================================================================================================================================
