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

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);
// ========================================================================================================================================
int main(int argc, char** argv) {

	gStyle->SetOptStat(0);

	if( argc != 2 ){
		cerr << "\nIncorrect number of arguments. Please use:\n" << "\t./efrain_check [inputFile]\n\n";
		return -1;
	}

	// Set up histograms for paddle-to-paddle and layer-to-layer
	const int nHistos = 600;

	double ref_meanTime_tdc[5] 	= {0.};
	double ref_meanTime_fadc[5] 	= {0.};
	double par_pad_tdc[nHistos][2] = {{0}};
	double par_pad_fadc[nHistos][2] = {{0}};
	double par_lay_tdc[5][2] = {{0.}};
	double par_lay_fadc[5][2] = {{0.}};

	TH2F ** h2_dMeantime_tdc_adc_paddle = new TH2F * [nHistos];
	TH2F ** h2_dMeantime_tdc_adc_layer  = new TH2F * [nHistos];
	TH2F ** h2_dMeantime_fadc_adc_paddle = new TH2F * [nHistos];
	TH2F ** h2_dMeantime_fadc_adc_layer  = new TH2F * [nHistos];

	for(int i = 0 ; i < nHistos ; i++){
		h2_dMeantime_tdc_adc_paddle[i] = new TH2F(Form("h2_dMeantime_tdc_adc_paddle_%i",i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}(bar - ref)",400,0,25000,300,-20,20);
		h2_dMeantime_tdc_adc_layer [i] = new TH2F(Form("h2_dMeantime_tdc_adc_layer_%i" ,i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}(bar - ref)",400,0,25000,300,-20,20);
		h2_dMeantime_fadc_adc_paddle[i] = new TH2F(Form("h2_dMeantime_fadc_adc_paddle_%i",i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}(bar - ref)",400,0,25000,300,-20,20);
		h2_dMeantime_fadc_adc_layer [i] = new TH2F(Form("h2_dMeantime_fadc_adc_layer_%i" ,i),";#sqrt{ADC_{L}ADC_{R}};((t_{L}+t_{R})/2)_{TDC}(bar - ref)",400,0,25000,300,-20,20);
	}

	// Opening input HIPO file
	hipo::reader reader;
	TString inputFile = argv[1];
	reader.open(inputFile);
	
	//Read Dictionary of Hipo File  // new hipo4
        hipo::dictionary  factory;      // new hipo4
        reader.readDictionary(factory); // new hipo4
        //factory.show();               // new hipo4

	BBand         band_hits   (factory.getSchema("BAND::hits"));

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
			double delta_meantime_tdc 	= meantimeTdc  - ref_meanTime_tdc[layer-1];
			double delta_meantime_fadc	= meantimeFadc - ref_meanTime_fadc[layer-1];
			double adc_geometric_mean 	= sqrt(adcLcorr*adcRcorr);

			// Now store all this in histograms:
			h2_dMeantime_tdc_adc_paddle[barKey] -> Fill(adc_geometric_mean,delta_meantime_tdc);
			h2_dMeantime_tdc_adc_paddle[barKey] -> SetTitle(Form("TDC: Sector: %i, Layer: %i, Component: %i",sector,layer,component));
			h2_dMeantime_fadc_adc_paddle[barKey] -> Fill(adc_geometric_mean,delta_meantime_fadc);
			h2_dMeantime_fadc_adc_paddle[barKey] -> SetTitle(Form("FADC: Sector: %i, Layer: %i, Component: %i",sector,layer,component));
			
			// If we have one of our reference bars, we can fill layer-by-layer offsets
			if( sector == 2 && component == 1){
				if( delta_meantime_tdc != 0. ) cout << "ERROR!\n";
				if( delta_meantime_fadc != 0.) cout << "ERROR!\n";
				double delta_meantime_layer_tdc  = meantimeTdc  - ref_meanTime_tdc[5-1];
				double delta_meantime_layer_fadc = meantimeFadc - ref_meanTime_fadc[5-1];

				if( delta_meantime_layer_tdc  != (ref_meanTime_tdc[layer-1]-ref_meanTime_tdc[5-1])   ) cout << "ERROR!\n";
				if( delta_meantime_layer_fadc != (ref_meanTime_fadc[layer-1]-ref_meanTime_fadc[5-1]) ) cout << "ERROR!\n";

				h2_dMeantime_tdc_adc_layer[barKey] -> Fill(adc_geometric_mean,delta_meantime_layer_tdc);
				h2_dMeantime_tdc_adc_layer[barKey] -> SetTitle(Form("TDC: Layer %i offset to Layer 5",layer));
				h2_dMeantime_fadc_adc_layer[barKey] -> Fill(adc_geometric_mean,delta_meantime_layer_fadc);
				h2_dMeantime_fadc_adc_layer[barKey] -> SetTitle(Form("FADC: Layer %i offset to Layer 5",layer));
			}
		}
	}// end file


	// Fitting functions
	TProfile ** p_dMeantime_tdc_adc_paddle = new TProfile*[nHistos];
	TProfile ** p_dMeantime_fadc_adc_paddle = new TProfile*[nHistos];
	// Fitting paddle-to-paddle parameters
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = (h2_dMeantime_tdc_adc_paddle[i]->Integral());
		if(notEmpty){
			int nBins;
			// For TDC
			nBins = h2_dMeantime_tdc_adc_paddle[i] -> GetXaxis() -> GetNbins();
			p_dMeantime_tdc_adc_paddle[i] = h2_dMeantime_tdc_adc_paddle[i]->ProfileX(Form("tdc_px_%i",i), 1, nBins );
			TF1 * f_line_tdc = new TF1("f_line_tdc","pol0",10000,18000);
			p_dMeantime_tdc_adc_paddle[i] -> Fit("f_line_tdc","QESR");

			par_pad_tdc[i][0] = f_line_tdc->GetParameter(0);
			par_pad_tdc[i][1] = f_line_tdc->GetParError (0);

			int nBins_fadc;
			// For FADC
			nBins_fadc = h2_dMeantime_fadc_adc_paddle[i] -> GetXaxis() -> GetNbins();
			p_dMeantime_fadc_adc_paddle[i] = h2_dMeantime_fadc_adc_paddle[i]->ProfileX(Form("fadc_px_%i",i), 1, nBins_fadc );
			TF1 * f_line_fadc = new TF1("f_line_fadc","pol0",10000,18000);
			p_dMeantime_fadc_adc_paddle[i] -> Fit("f_line_fadc","QESR");

			par_pad_fadc[i][0] = f_line_fadc->GetParameter(0);
			par_pad_fadc[i][1] = f_line_fadc->GetParError (0);

			cout << "Paddle: " << i << " " << par_pad_tdc[i][0] << " " << par_pad_fadc[i][0] << "\n";

		}
	}

	
	TProfile ** p_dMeantime_tdc_adc_layer = new TProfile*[5];
	TProfile ** p_dMeantime_fadc_adc_layer = new TProfile*[5];
	// Fitting layer-to-layer parameters
	for(int i = 0 ; i < 5 ; i++){
		int idx = 201+10*(i+1);
		int nBins;

		// For TDC:
		nBins = h2_dMeantime_tdc_adc_layer[idx] -> GetXaxis() -> GetNbins();
		p_dMeantime_tdc_adc_layer[i] = h2_dMeantime_tdc_adc_layer[idx] ->ProfileX(Form("tdc_lx_%i",i), 1, nBins );
		TF1 * f_line_tdc = new TF1("f_line_tdc","pol0",10000,18000);
		p_dMeantime_tdc_adc_layer[i] -> Fit("f_line_tdc","QESR");
		
		par_lay_tdc[i][0] = f_line_tdc->GetParameter(0);
		par_lay_tdc[i][1] = f_line_tdc->GetParError (0);

		int nBins_fadc;
		// For FADC:
		nBins_fadc = h2_dMeantime_fadc_adc_layer[idx] -> GetXaxis() -> GetNbins();
		p_dMeantime_fadc_adc_layer[i] = h2_dMeantime_fadc_adc_layer[idx] ->ProfileX(Form("fadc_lx_%i",i), 1, nBins );
		TF1 * f_line_fadc = new TF1("f_line_fadc","pol0",10000,18000);
		p_dMeantime_fadc_adc_layer[i] -> Fit("f_line_fadc","QESR");
		
		par_lay_fadc[i][0] = f_line_fadc->GetParameter(0);
		par_lay_fadc[i][1] = f_line_fadc->GetParError (0);
			
		cout << "Layer: " << i << " " << par_lay_tdc[i][0] << " " << par_lay_fadc[i][0] << "\n";
	}


	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0","c0",900,900);

	TCanvas *** cSLC_paddle = new TCanvas**[5];

	for(int is = 0 ; is < 5 ; is++){
		cSLC_paddle[is] = new TCanvas*[6];

		for(int il = 0 ; il < 5 ; il++){
			cSLC_paddle[is][il] = new TCanvas(Form("S%iL%i_paddle",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			cSLC_paddle[is][il] -> Divide(2,8);	

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){

				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h2_dMeantime_tdc_adc_paddle[identifier]->Integral());
				if(notEmpty){
					int min,max;

					// Histogram on left is TDC: meantime - ref
					cSLC_paddle[is][il] -> cd(2*cIdx+1);
					gPad -> SetBottomMargin(0.26);
					gPad -> SetLeftMargin(0.14);
						// set range
					min = getMinNonEmptyBin(h2_dMeantime_tdc_adc_paddle[identifier]);
					max = getMaxNonEmptyBin(h2_dMeantime_tdc_adc_paddle[identifier]);
					h2_dMeantime_tdc_adc_paddle[identifier] -> GetYaxis() -> SetRange(min,max);
						// draw
					h2_dMeantime_tdc_adc_paddle[identifier] -> Draw("COLZ");
					p_dMeantime_tdc_adc_paddle[identifier] -> Draw("same");	
					
					// Histogram on right is FADC: meantime - ref
					cSLC_paddle[is][il] -> cd(2*cIdx+2);
					gPad -> SetBottomMargin(0.26);
					gPad -> SetLeftMargin(0.14);
						// set range
					min = getMinNonEmptyBin(h2_dMeantime_fadc_adc_paddle[identifier]);
					max = getMaxNonEmptyBin(h2_dMeantime_fadc_adc_paddle[identifier]);
					h2_dMeantime_fadc_adc_paddle[identifier] -> GetYaxis() -> SetRange(min,max);
						// draw
					h2_dMeantime_fadc_adc_paddle[identifier] -> Draw("COLZ");
					p_dMeantime_fadc_adc_paddle[identifier] -> Draw("same");		
				} 
			}
			// Last histogram is the layer histogram
			if( is == 1 ){
				int min,max;
				int refID = 200+10*(il+1)+1;
				// Draw TDC on left
				cSLC_paddle[is][il] -> cd(15);
				gPad -> SetBottomMargin(0.26);
				gPad -> SetLeftMargin(0.14);
					// set range
				min = getMinNonEmptyBin(h2_dMeantime_tdc_adc_layer[refID]);
				max = getMaxNonEmptyBin(h2_dMeantime_tdc_adc_layer[refID]);
				h2_dMeantime_tdc_adc_layer[refID] -> GetYaxis() -> SetRange(min,max);
					// draw
				h2_dMeantime_tdc_adc_layer[refID] -> Draw("COLZ");
				p_dMeantime_tdc_adc_layer[il] -> Draw("same");		
		
				cSLC_paddle[is][il] -> cd(16);
				gPad -> SetBottomMargin(0.26);
				gPad -> SetLeftMargin(0.14);
					// set range
				min = getMinNonEmptyBin(h2_dMeantime_fadc_adc_layer[refID]);
				max = getMaxNonEmptyBin(h2_dMeantime_fadc_adc_layer[refID]);
				h2_dMeantime_fadc_adc_layer[refID] -> GetYaxis() -> SetRange(min,max);
					// draw
				h2_dMeantime_fadc_adc_layer[refID] -> Draw("COLZ");
				p_dMeantime_fadc_adc_layer[il] -> Draw("same");		
			}

			cSLC_paddle[is][il] -> Modified();	cSLC_paddle[is][il] -> Update();		
		}
	}

	// Saving plots to a pdf file
	c0 -> Print("results_paddle_corr.pdf(");
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			cSLC_paddle[is][il] -> Print("results_paddle_corr.pdf");
		}
	}
	c0 -> Print("results_paddle_corr.pdf)");

	// Create some txt files for the offsets
	ofstream tab_pad_tdc, tab_lay_tdc;
	ofstream tab_pad_fadc, tab_lay_fadc;
	tab_pad_tdc.open("TDC_paddle_offsets.txt");
	tab_lay_tdc.open("TDC_layer_offsets.txt" );
	tab_pad_fadc.open("FADC_paddle_offsets.txt");
	tab_lay_fadc.open("FADC_layer_offsets.txt" );

	for(int il = 1 ; il <= 6 ; il++){
		for(int is = 1 ; is <= 5 ; is++){
			for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
				int idx = 100*is + 10*il + ic;

				// TDC paddle to paddle:
				tab_pad_tdc << is << "\t" << il << "\t" << ic << "\t";
				if( il==6 || (is==2 && ic==1)) tab_pad_tdc << "0.0000\t0.0000\n";
				else tab_pad_tdc << par_pad_tdc[idx][0] << "\t" << par_pad_tdc[idx][1] << "\n";
				// TDC layer by layer:
				tab_lay_tdc << is << "\t" << il << "\t" << ic << "\t";
				if( il==6 || il==5 ) tab_lay_tdc << "0.0000\t0.0000\n";
				else tab_lay_tdc << par_lay_tdc[il-1][0] << "\t" << par_lay_tdc[il-1][1] << "\n";

				// FADC paddle to paddle:
				tab_pad_fadc << is << "\t" << il << "\t" << ic << "\t";
				if( il == 6 || (is==2 && ic==1) ) tab_pad_fadc << "0.0000\t0.0000\n";
				else tab_pad_fadc << par_pad_fadc[idx][0] << "\t" << par_pad_fadc[idx][1] << "\n";
				// FADC layer by layer:
				tab_lay_fadc << is << "\t" << il << "\t" << ic << "\t";
				if( il == 6 || il == 5) tab_lay_fadc << "0.0000\t0.0000\n";
				else tab_lay_fadc << par_lay_fadc[il-1][0] << "\t" << par_lay_fadc[il-1][1] << "\n";
			}
		}
	}
	tab_pad_tdc.close();
	tab_lay_tdc.close();
	tab_pad_fadc.close();
	tab_lay_fadc.close();


	return 0;
}


int getMinNonEmptyBin(TH2F * h2){
        int nBinX = h2 -> GetXaxis() -> GetNbins();
        int nBinY = h2 -> GetYaxis() -> GetNbins();

        int binContent = 0;
	int nextBinContent = 0;

        for(int i = 1 ; i <= nBinY ; i++){
                for(int j = 1 ; j <= nBinX ; j++){
                        binContent     = h2 -> GetBinContent(j  ,i);
                        nextBinContent = h2 -> GetBinContent(j+1,i);
			if(binContent!=0&&nextBinContent!=0) return i;
			else{
				binContent     = 0;
				nextBinContent = 0;
			}
                }
        }
        return 1;
}
int getMaxNonEmptyBin(TH2F * h2){
        int nBinX = h2 -> GetXaxis() -> GetNbins();
        int nBinY = h2 -> GetYaxis() -> GetNbins();

        int binContent = 0;
	int nextBinContent = 0;

        for(int i = nBinY ; i > 0 ; i--){
                for(int j = 1 ; j <= nBinX ; j++){
                        binContent = h2 -> GetBinContent(j,i);
                        nextBinContent = h2 -> GetBinContent(j-1,i);
			if(binContent!=0) return i;
			 else{
                                binContent     = 0;
                                nextBinContent = 0;
                        }
                }
        }
        return nBinY;
}
