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

/*
int S1_comp[] = {3,3,3,3,3,3};
int S2_comp[] = {7,7,7,7,7,7};
int S3_comp[] = {6,6,6,6,5,6};
int S4_comp[] = {6,6,6,6,5,6};
int S5_comp[] = {2,2,2,2,0,2};
*/

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

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

	double parL[nHistos][4] = {{0}};
	double parR[nHistos][4] = {{0}};

	TH2F ** h2_tdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_tdc_adc_R = new TH2F * [nHistos];

	TH2F * h2_tdc_adc_tot_L = new TH2F("h2_tdc_adc_tot_L","Left PMTs;ADC;t_{TDC}-t_{FADC}" ,200,0,30000,200,350,450);
	TH2F * h2_tdc_adc_tot_R = new TH2F("h2_tdc_adc_tot_R","Right PMTs;ADC;t_{TDC}-t_{FADC}",200,0,30000,200,350,450);

	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_adc_L[i] = new TH2F(Form("h2_tdc_adc_L_%i",i),";ADC;t_{TDC}-t_{FADC}",200,0,30000,200,350,450);
		h2_tdc_adc_R[i] = new TH2F(Form("h2_tdc_adc_R_%i",i),";ADC;t_{TDC}-t_{FADC}",200,0,30000,200,350,450);
	}

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	hipo::bank BAND_ADC("BAND::adc",reader);
	hipo::bank BAND_TDC("BAND::tdc",reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;	

		//BAND_ADC.show();
		//BAND_TDC.show();

		int nADC = BAND_ADC.getSize();
		int nTDC = BAND_TDC.getSize();

		// Skip events with no entries
		if(nADC==0||nTDC==0) continue;

		for(int aIdx = 0 ; aIdx < nADC ; aIdx++){

			int   ADC_sector    = BAND_ADC.getInt  (0,aIdx);
			int   ADC_layer     = BAND_ADC.getInt  (1,aIdx);
			int   ADC_component = BAND_ADC.getInt  (2,aIdx);
			int   ADC_order     = BAND_ADC.getInt  (3,aIdx);
			float ADC_adc       = (float)(BAND_ADC.getInt(4,aIdx));
			float ADC_time      = BAND_ADC.getFloat(5,aIdx);

			for(int tIdx = 0 ; tIdx < nTDC ; tIdx++){
				int   TDC_sector    = BAND_TDC.getInt  (0,tIdx);
				int   TDC_layer     = BAND_TDC.getInt  (1,tIdx);
				int   TDC_component = BAND_TDC.getInt  (2,tIdx);
				int   TDC_order     = BAND_TDC.getInt  (3,tIdx);
				float TDC_tdc       = (float)(BAND_TDC.getInt(4,0));
				float TDC_time      = TDC_tdc*0.02345;

				if(		(ADC_sector   ==TDC_sector   )&&
						(ADC_layer    ==TDC_layer    )&&
						(ADC_component==TDC_component)&&
						(ADC_order +2 ==TDC_order    )
				  ){
					int barKey = 100*ADC_sector + 10*ADC_layer + ADC_component;

					// Left PMTs
					if(ADC_order==0){
						h2_tdc_adc_L[barKey] -> SetTitle(Form("Left PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
						h2_tdc_adc_L[barKey] -> Fill(ADC_adc,TDC_time-ADC_time);
						h2_tdc_adc_tot_L     -> Fill(ADC_adc,TDC_time-ADC_time);
					}

					// Right PMTs
					if(ADC_order==1){
						h2_tdc_adc_R[barKey] -> SetTitle(Form("Right PMTs, Sector:%i, Layer:%i, Component:%i",ADC_sector,ADC_layer,ADC_component));
						h2_tdc_adc_R[barKey] -> Fill(ADC_adc,TDC_time-ADC_time);
						h2_tdc_adc_tot_R     -> Fill(ADC_adc,TDC_time-ADC_time);
					}

				}

			}
		}

	}// end file

	// -------------------------------------------------------------------------------------------------
	// Fitting functions
	for(int i = 0 ; i < nHistos ; i++){
                int notEmpty = h2_tdc_adc_L[i] -> Integral();
                if(notEmpty){
			TF1 * f_timewalk_L = new TF1("f_timewalk_L","[0]/TMath::Sqrt(x)+[1]");
			TF1 * f_timewalk_R = new TF1("f_timewalk_R","[0]/TMath::Sqrt(x)+[1]");

                        h2_tdc_adc_L[i] -> Fit("f_timewalk_L");
                        h2_tdc_adc_R[i] -> Fit("f_timewalk_R");
                
			parL[i][0] = f_timewalk_L -> GetParameter(0);	//A
			parL[i][1] = f_timewalk_L -> GetParameter(1);	//B
			parL[i][2] = f_timewalk_L -> GetParError (0);	//eA
			parL[i][3] = f_timewalk_L -> GetParError (1);	//eB

			parR[i][0] = f_timewalk_R -> GetParameter(0);   //A
                        parR[i][1] = f_timewalk_R -> GetParameter(1);   //B
                        parR[i][2] = f_timewalk_R -> GetParError (0);   //eA
                        parR[i][3] = f_timewalk_R -> GetParError (1);   //eB
		}
        }

	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	TCanvas * c0 = new TCanvas("c0");
	c0 -> Divide(2,1);
	c0 -> cd(1);	h2_tdc_adc_tot_L -> Draw("COLZ");
	c0 -> cd(2);	h2_tdc_adc_tot_R -> Draw("COLZ");
	c0 -> Modified();
	c0 -> Update();

	TCanvas ** cS1 = new TCanvas*[5];
	TCanvas ** cS2 = new TCanvas*[5];
	TCanvas ** cS3 = new TCanvas*[5];
	TCanvas ** cS4 = new TCanvas*[5];
	TCanvas ** cS5 = new TCanvas*[5];

	for(int i = 0 ; i < 5 ; i++){
		cS1[i] = new TCanvas(Form("S1L%i",i+1),Form("Sector 1, Layer %i",i+1),600,800);	cS1[i] -> Divide(2,7);
		cS2[i] = new TCanvas(Form("S2L%i",i+1),Form("Sector 2, Layer %i",i+1),600,800);	cS2[i] -> Divide(2,7);
		cS3[i] = new TCanvas(Form("S3L%i",i+1),Form("Sector 3, Layer %i",i+1),600,800);	cS3[i] -> Divide(2,7);
		cS4[i] = new TCanvas(Form("S4L%i",i+1),Form("Sector 4, Layer %i",i+1),600,800);	cS4[i] -> Divide(2,7);
		cS5[i] = new TCanvas(Form("S5L%i",i+1),Form("Sector 5, Layer %i",i+1),600,800);	cS5[i] -> Divide(2,7);

		// Sector 1
		for(int cIdx = 0 ; cIdx < 7 ; cIdx++){
			int identifier = 100+10*(i+1)+(cIdx+1);
			int notEmpty = h2_tdc_adc_L[identifier] -> Integral();
			if(notEmpty){
				cS1[i] -> cd(2*cIdx+1);	h2_tdc_adc_L[identifier] -> Draw("COLZ");
				cS1[i] -> cd(2*cIdx+2);	h2_tdc_adc_R[identifier] -> Draw("COLZ");
			}
		}
		cS1[i] -> Modified();	cS1[i] -> Update();

		// Sector 2
                for(int cIdx = 0 ; cIdx < 7 ; cIdx++){
                        int identifier = 200+10*(i+1)+(cIdx+1);
                        int notEmpty = h2_tdc_adc_L[identifier] -> Integral();
                        if(notEmpty){
                                cS2[i] -> cd(2*cIdx+1);	h2_tdc_adc_L[identifier] -> Draw("COLZ");
                                cS2[i] -> cd(2*cIdx+2);	h2_tdc_adc_R[identifier] -> Draw("COLZ");
                        }
                }
                cS2[i] -> Modified();   cS2[i] -> Update();

		// Sector 3
                for(int cIdx = 0 ; cIdx < 7 ; cIdx++){
                        int identifier = 300+10*(i+1)+(cIdx+1);
                        int notEmpty = h2_tdc_adc_L[identifier] -> Integral();
                        if(notEmpty){
                                cS3[i] -> cd(2*cIdx+1);	h2_tdc_adc_L[identifier] -> Draw("COLZ");
                                cS3[i] -> cd(2*cIdx+2);	h2_tdc_adc_R[identifier] -> Draw("COLZ");
                        }
                }
                cS3[i] -> Modified();   cS3[i] -> Update();

		// Sector 4
                for(int cIdx = 0 ; cIdx < 7 ; cIdx++){
                        int identifier = 400+10*(i+1)+(cIdx+1);
                        int notEmpty = h2_tdc_adc_L[identifier] -> Integral();
                        if(notEmpty){
                                cS4[i] -> cd(2*cIdx+1);	h2_tdc_adc_L[identifier] -> Draw("COLZ");
                                cS4[i] -> cd(2*cIdx+2);	h2_tdc_adc_R[identifier] -> Draw("COLZ");
                        }
                }
                cS4[i] -> Modified();   cS4[i] -> Update();

		// Sector 5
                for(int cIdx = 0 ; cIdx < 7 ; cIdx++){
                        int identifier = 500+10*(i+1)+(cIdx+1);
                        int notEmpty = h2_tdc_adc_L[identifier] -> Integral();
                        if(notEmpty){
                                cS5[i] -> cd(2*cIdx+1);	h2_tdc_adc_L[identifier] -> Draw("COLZ");
                                cS5[i] -> cd(2*cIdx+2);	h2_tdc_adc_R[identifier] -> Draw("COLZ");
                        }
                }
                cS5[i] -> Modified();   cS5[i] -> Update();
		
	}

	// -------------------------------------------------------------------------------------------------
	// Saving fit values to ccdb tables
	for(int is = 1 ; is <= 5 ; is++){	
		for(int il = 1 ; il <= 6 ; il++){
                	for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
                        	int idx = 100*is + 10*il + ic;
                	        cout << is << "\t" << il << "\t" << ic << "\t" << "\t";
				cout << parL[idx][0] << "\t" << parL[idx][1] << "\t" << parL[idx][2]  << "\t" << parL[idx][3] << endl;
        	        }
	        }
	}

	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	for(int i = 0 ; i < nHistos ; i++){
		int notEmpty = h2_tdc_adc_L[i] -> Integral();
		if(!notEmpty){
                	delete h2_tdc_adc_L[i];
                	delete h2_tdc_adc_R[i];
        	}
	}

	myapp -> Run();
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

