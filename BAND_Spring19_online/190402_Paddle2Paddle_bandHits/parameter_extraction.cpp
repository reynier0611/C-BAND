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

	cout << "*************************************************************" << endl;
	cout << "                  ___                  ___" << endl;
	cout << "\\    /\\    / /\\   |  \\ |\\  | | |\\  |  /" << endl;
	cout << " \\  /  \\  / /--\\  |-\\/ | \\ | | | \\ | |  --\\" << endl;
	cout << "  \\/    \\/ /    \\ |  \\ |  \\| | |  \\|  \\___/" << endl << endl;
	cout << "Run this code on laser data on a file that does not have" << endl;
        cout << "paddle-to-paddle corrections yet, but already has time-walk" << endl;
        cout << "and L-R corrections" << endl << endl;
	cout << "*************************************************************" << endl;

	// ----------------------------------------------------------------------------------
	// Declaring histograms
	const int nHistos = 600;

	double ref_meantime[10] = {0};
	double layer_ref_meantime = 0;

	double par_pad[nHistos][2] = {{0}};
	double par_lay[nHistos][2] = {{0}};
	double par_res[nHistos][2] = {{0}};

	TH2F ** h2_dMeantime_adc_paddle = new TH2F * [nHistos];
	TH2F ** h2_dMeantime_adc_layer  = new TH2F * [nHistos];
	TH1F ** h1_resolutions          = new TH1F * [nHistos];

	for(int i = 0 ; i < nHistos ; i++){
		h2_dMeantime_adc_paddle[i] = new TH2F(Form("h2_dMeantime_adc_paddle_%i",i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}(bar - ref)",400,0,25000,401,-15,15);
		h2_dMeantime_adc_layer [i] = new TH2F(Form("h2_dMeantime_adc_layer_%i" ,i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}(bar - ref)",400,0,25000,401,-15,15);
		PrettyTH2F(h2_dMeantime_adc_paddle[i]);
		PrettyTH2F(h2_dMeantime_adc_layer [i]);

		h1_resolutions         [i] = new TH1F(Form("h1_resolutions_%i"         ,i),";((t_{L}+t_{R})/2)_{TDC}(bar - ref)",401,-15,15);
		PrettyTH1F(h1_resolutions[i],62);
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
			int    sector            = band_hits.getSector     (hit);
                        int    layer             = band_hits.getLayer      (hit);
                        int    component         = band_hits.getComponent  (hit);
                        float  tTdcLcorr         = band_hits.getTtdcLcorr  (hit);
                        float  tTdcRcorr         = band_hits.getTtdcRcorr  (hit);
			float  meantimeTdc       = band_hits.getMeantimeTdc(hit);

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

			if(adcLcorr< 6000||adcRcorr< 6000) continue; // To avoid area sensitive to time-walk corrections
			if(adcLcorr>20000||adcRcorr>20000) continue; // To avoid area sensitive to adc overflow

			double delta_meantime = meantimeTdc - ref_meantime[layer];
			double adc_geometric_mean = TMath::Sqrt(adcLcorr*adcRcorr);

			h2_dMeantime_adc_paddle[barKey] -> Fill(adc_geometric_mean,delta_meantime);
			h2_dMeantime_adc_paddle[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",sector,layer,component));

			// 1-D projection to extract resolutions
			h1_resolutions[barKey] -> Fill(delta_meantime);
			h1_resolutions[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",sector,layer,component));	
			
			if(sector==2&&component==1){
				double delta_meantime_layer = meantimeTdc - ref_meantime[5];
				h2_dMeantime_adc_layer[barKey] -> Fill(adc_geometric_mean,delta_meantime_layer);
				h2_dMeantime_adc_layer[barKey] -> SetTitle(Form("Sector: %i, Layer: %i, Component: %i",sector,layer,component));
			}
		}
	}// end file

	// -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// Fitting functions
	TProfile ** p_dMeantime_adc_paddle = new TProfile*[nHistos];
	// Fitting paddle-to-paddle parameters
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_dMeantime_adc_paddle[i]->Integral());
		if(notEmpty){
			int nBins = h2_dMeantime_adc_paddle[i] -> GetXaxis() -> GetNbins();
			p_dMeantime_adc_paddle[i] = h2_dMeantime_adc_paddle[i]->ProfileX(Form("px_%i",i), 1, nBins );
			TF1 * f_line = new TF1("f_line","pol0");
			p_dMeantime_adc_paddle[i] -> Fit("f_line","Q");

			par_pad[i][0] = f_line->GetParameter(0);
			par_pad[i][1] = f_line->GetParError (0);
		}
	}

	TProfile ** p_dMeantime_adc_layer = new TProfile*[5];
	// Fitting layer-to-layer parameters
	for(int i = 0 ; i < 5 ; i++){
		TF1 * f_line = new TF1("f_line","pol0");

		int idx = 201+10*(i+1);
		int nBins = h2_dMeantime_adc_layer[i] -> GetXaxis() -> GetNbins();
		p_dMeantime_adc_layer[i] = h2_dMeantime_adc_layer[idx] ->ProfileX(Form("lx_%i",i), 1, nBins );
		p_dMeantime_adc_layer[i] -> Fit("f_line","Q");

		par_lay[i][0] = f_line->GetParameter(0);
		par_lay[i][1] = f_line->GetParError (0);

	}

	// Fitting 1-d projections to get resolutions
	for(int i = 0 ; i < nHistos ; i++){
                int notEmpty = (h1_resolutions[i]->Integral());
		if(notEmpty){
			double mean = h1_resolutions[i]->GetMean();
			TF1 * f_gaus = new TF1("f_gaus","gaus",mean-5,mean+5);
			f_gaus -> SetParameters(300,0,200);
			f_gaus -> SetParLimits(1,-15, 15);
			f_gaus -> SetParLimits(2,  0,600);

			h1_resolutions[i] -> Fit("f_gaus","QR");

			if((i!=211)&&(i!=221)&&(i!=231)&&(i!=241)&&(i!=251)){ // Don't fill parameters for reference bars
				par_res[i][0] = f_gaus->GetParameter(2);
				par_res[i][1] = f_gaus->GetParError (2);
			}
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0","c0",900,900);

	TCanvas *** cSLC_paddle = new TCanvas**[5];
	TCanvas *** cSLC_resol  = new TCanvas**[5];

	for(int is = 0 ; is < 5 ; is++){
		cSLC_paddle[is] = new TCanvas*[6];
		cSLC_resol [is] = new TCanvas*[6];

		for(int il = 0 ; il < 5 ; il++){
			cSLC_paddle[is][il] = new TCanvas(Form("S%iL%i_paddle",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cSLC_paddle[is][il] -> Divide(2,4);	

			cSLC_resol [is][il] = new TCanvas(Form("S%iL%i_resol" ,is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cSLC_resol [is][il] -> Divide(2,4); 

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
					p_dMeantime_adc_paddle[identifier] -> Draw("same");		

					TLatex * tex_pad = new TLatex(1000,par_pad[identifier][0]-5,Form("y = %f #pm %f ns",par_pad[identifier][0],par_pad[identifier][1]));
					tex_pad -> SetTextSize(0.08);
					tex_pad -> Draw("same");

					// Resolutions
					cSLC_resol [is][il] -> cd(cIdx+1);
					gPad -> SetBottomMargin(0.26);
                                        gPad -> SetLeftMargin(0.14);
					h1_resolutions[identifier] -> Draw();

					TLatex * tex_lay = new TLatex(-10,(h1_resolutions[identifier]->GetMaximum()/2.),Form("#sigma = %.1f #pm %.1f ps",1000*par_res[identifier][0],1000*par_res[identifier][1]));
					tex_lay -> SetTextSize(0.08);
                                        tex_lay -> Draw("same");
				} 
			}
			cSLC_paddle[is][il] -> Modified();	cSLC_paddle[is][il] -> Update();		
			cSLC_resol [is][il] -> Modified();	cSLC_resol [is][il] -> Update();
		}
	}

	TCanvas * cSLC_layer = new TCanvas("cSLC_layer","cSLC_layer",900,900);
	cSLC_layer -> Divide(2,3);
	for(int i = 0 ; i < 5 ; i++){
		cSLC_layer -> cd(i+1);
		gPad -> SetBottomMargin(0.26);
		gPad -> SetLeftMargin(0.14);

		int idx = 201 + 10*(i+1);
		h2_dMeantime_adc_layer[idx] -> Draw("COLZ");
		p_dMeantime_adc_layer[i] -> Draw("same");

		TLatex * tex_lay = new TLatex(1000,par_lay[i][0]-5,Form("y = %f #pm %f ns",par_lay[i][0],par_lay[i][1]));
                tex_lay -> SetTextSize(0.08);
                tex_lay -> Draw("same");
	}
	cSLC_layer -> Modified();
	cSLC_layer -> Update();

	// -------------------------------------------------------------------------------------------------
	// Making a plot with resolutions
	int ctr = 0;
	double x_val[600] = {0};
	double y_val[600] = {0};
	for(int il = 1 ; il < 6 ; il++){
		for(int is = 1 ; is <= 5 ; is++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int I_BarKey = 100*is + 10*il + ic;

				x_val[ctr] = (double)(I_BarKey);
				y_val[ctr] = 1000*par_res[I_BarKey][0]/TMath::Sqrt(2);
				ctr++;	
			}
		}
	}

	TGraphErrors * gResolutions = new TGraphErrors(ctr,x_val,y_val,0,0);
	PrettyTGraphErrors(gResolutions,62);

	TCanvas * cres = new TCanvas("cres","cres",900,900);
	gResolutions -> Draw("AP");
	cres -> Modified();
	cres -> Update();

	// -------------------------------------------------------------------------------------------------
	// Saving fit values to ccdb tables
	ofstream tab_pad, tab_lay;
	tab_pad.open("TDC_paddle_offsets.txt");
	tab_lay.open("TDC_layer_offsets.txt" );

	for(int il = 1 ; il <= 6 ; il++){
		for(int is = 1 ; is <= 5 ; is++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;

				tab_pad << is << "\t" << il << "\t" << ic << "\t";
				if(il==6||(is==2&&ic==1)) tab_pad << "0.0000\t0.0000" << endl;
				else tab_pad << par_pad[idx][0] << "\t" << par_pad[idx][1] << endl;

				tab_lay << is << "\t" << il << "\t" << ic << "\t";
				if(il==6||il==5) tab_lay << "0.0000\t0.0000" << endl;
				else tab_lay << par_lay[il-1][0] << "\t" << par_lay[il-1][1] << endl;
			}
		}
	}
	tab_pad.close();
	tab_lay.close();

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_dMeantime_adc_paddle[i]->Integral());
		if(!notEmpty){
			delete h2_dMeantime_adc_paddle[i];
		}

		notEmpty = (h2_dMeantime_adc_layer[i]->Integral());
		if(!notEmpty){
			delete h2_dMeantime_adc_layer[i];
		}
	}

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	c0 -> Print("results_paddle_corr.pdf(");
	cSLC_layer -> Print("results_paddle_corr.pdf");
	c0 -> Print("results_paddle_corr.pdf");
	cres -> Print("results_resolutions.pdf(");

	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			cSLC_paddle[is][il] -> Print("results_paddle_corr.pdf");
			cSLC_resol [is][il] -> Print("results_resolutions.pdf");
		}
	}
	c0 -> Print("results_paddle_corr.pdf)");
	c0 -> Print("results_resolutions.pdf)");

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
