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

#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "Rtypes.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

// Forward-declaring functions
void PrettyTH2F(TH2F * h2);
double getTriggerPhase( long timeStamp );
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);
int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE);
void walkCorr(	std::vector<double> *widths		,
		std::vector<double> *widthsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				);
void fitGraph( 	std::vector<double> widths	,
		std::vector<double> times	,
		std::vector<double> widthsErr	,
		std::vector<double> timesErr	,
		double &par0			,
		double &par1			,
		TCanvas *c			,
		int cd				,
		TString name			);
double wlk( double *x , double *p);

const int nADCbins = 400;
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// ========================================================================================================================================
int main(int argc, char** argv) {
	/*
#ifdef WITHRINT
TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
TApplication *myapp = new TApplication("myapp",0,0);
#endif
*/
	//gStyle->SetTitleSize(0.3,"t");
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
        gStyle->SetStatX(0.8);
        gStyle->SetStatY(0.9);
        gStyle->SetStatW(0.2);
        gStyle->SetStatH(0.2);

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
	TH2F ** h2_ftdc_adc_L = new TH2F * [nHistos];
	TH2F ** h2_ftdc_adc_R = new TH2F * [nHistos];
	TH2F ** h2_tdc_adc_L_Cut = new TH2F * [nHistos];
	TH2F ** h2_tdc_adc_R_Cut = new TH2F * [nHistos];
	TH2F ** h2_ftdc_adc_L_Cut = new TH2F * [nHistos];
	TH2F ** h2_ftdc_adc_R_Cut = new TH2F * [nHistos];
	TH2F ** h2_tdcMean = new TH2F * [nHistos];
	TH2F ** h2_ftdcMean = new TH2F * [nHistos];
	TH2F ** h2_tdcMean_Cut = new TH2F * [nHistos];
	TH2F ** h2_ftdcMean_Cut = new TH2F * [nHistos];
	TH2F ** h2_tdiff_L	= new TH2F * [nHistos];
	TH2F ** h2_tdiff_R	= new TH2F * [nHistos];
	int hasEvents[nHistos] = {0};

	//TH2F * h2_tdc_adc_tot_L = new TH2F("h2_tdc_adc_tot_L","Left PMTs;ADC;t_{TDC}-t_{FADC} [ns]" ,2000,1000,32500,600,70,130);
	//TH2F * h2_tdc_adc_tot_R = new TH2F("h2_tdc_adc_tot_R","Right PMTs;ADC;t_{TDC}-t_{FADC} [ns]",2000,1000,32500,600,70,130);
	//TH2F * h2_tdc_adc_tot_L = new TH2F("h2_tdc_adc_tot_L","Left PMTs;ADC;t_{TDC}-t_{FADC} [ns]" ,2000,1000,32500,100,-5,5);
	//TH2F * h2_tdc_adc_tot_R = new TH2F("h2_tdc_adc_tot_R","Right PMTs;ADC;t_{TDC}-t_{FADC} [ns]",2000,1000,32500,100,-5,5);
	//TH2F * h2_tdc_adc_tot_L = new TH2F("h2_tdc_adc_tot_L","Left PMTs;ADC;t_{TDC}-t_{FADC} [ns]" ,200,1000,32500,600,1185,1245);
	//TH2F * h2_tdc_adc_tot_R = new TH2F("h2_tdc_adc_tot_R","Right PMTs;ADC;t_{TDC}-t_{FADC} [ns]",200,1000,32500,600,1185,1245);
	//h2_tdc_adc_tot_L -> GetYaxis() -> SetTitleOffset(1.40);
	//h2_tdc_adc_tot_L -> GetXaxis() -> SetNdivisions(509);
	//h2_tdc_adc_tot_R -> GetYaxis() -> SetTitleOffset(1.40);
	//h2_tdc_adc_tot_R -> GetXaxis() -> SetNdivisions(509);

	for(int i = 0 ; i < nHistos ; i++){
		h2_tdc_adc_L[i] = new TH2F(Form("h2_tdc_adc_L_%i",i),";ADC_{L};	  t_{L,TDC} - Ref [ns]",		nADCbins,1000,33000,	1200,-30,30);
		h2_tdc_adc_R[i] = new TH2F(Form("h2_tdc_adc_R_%i",i),";ADC_{R};	  t_{R,TDC} - Ref [ns]",		nADCbins,1000,33000,	1200,-30,30);
		h2_ftdc_adc_L[i] = new TH2F(Form("h2_ftdc_adc_L_%i",i),";ADC_{L};t_{L,FTDC} - Ref [ns]",		nADCbins,1000,33000,	1200,-30,30);
		h2_ftdc_adc_R[i] = new TH2F(Form("h2_ftdc_adc_R_%i",i),";ADC_{R};t_{R,FTDC} - Ref [ns]",		nADCbins,1000,33000,	1200,-30,30);
		h2_tdc_adc_L_Cut[i] = new TH2F(Form("h2_tdc_adc_L_Cut_%i",i),";ADC_{L};	  t_{L,TDC} - Ref [ns]",	nADCbins,1000,33000,	1200,-30,30);
		h2_tdc_adc_R_Cut[i] = new TH2F(Form("h2_tdc_adc_R_Cut_%i",i),";ADC_{R};	  t_{R,TDC} - Ref [ns]",	nADCbins,1000,33000,	1200,-30,30);
		h2_ftdc_adc_L_Cut[i] = new TH2F(Form("h2_ftdc_adc_L_Cut_%i",i),";ADC_{L};t_{L,FTDC} - Ref [ns]",	nADCbins,1000,33000,	1200,-30,30);
		h2_ftdc_adc_R_Cut[i] = new TH2F(Form("h2_ftdc_adc_R_Cut_%i",i),";ADC_{R};t_{R,FTDC} - Ref [ns]",	nADCbins,1000,33000,	1200,-30,30);
		h2_tdcMean[i]	= new TH2F(Form("h2_tdcMean_%i",i),";GM ADC;   TDC Mean T - Ref [ns]",        		nADCbins,1000,33000, 1200,-30,30);
		h2_ftdcMean[i]   = new TH2F(Form("h2_ftdcMean_%i",i),";GM ADC;   FTDC Mean T - Ref [ns]",        	nADCbins,1000,33000, 1200,-30,30);
		h2_tdcMean_Cut[i]	= new TH2F(Form("h2_tdcMean_Cut_%i",i),";GM ADC;   TDC Mean T - Ref [ns]",      nADCbins,1000,33000, 1200,-30,30);
		h2_ftdcMean_Cut[i]   = new TH2F(Form("h2_ftdcMean_Cut_%i",i),";GM ADC;   FTDC Mean T - Ref [ns]",       nADCbins,1000,33000, 1200,-30,30);
		// laser or cosmic standalone
		h2_tdiff_L[i]	= new TH2F(Form("h2_tdiff_L_%i",i),";ADC L; TDC L - FADC L [ns]",			nADCbins,1000,33000, 200,1220,1240);
		h2_tdiff_R[i]	= new TH2F(Form("h2_tdiff_R_%i",i),";ADC R; TDC R - FADC R [ns]",			nADCbins,1000,33000, 200,1220,1240);
		// laser with CLAS
		//h2_tdiff_L[i]	= new TH2F(Form("h2_tdiff_L_%i",i),";ADC L; TDC L - FADC L [ns]",			nADCbins,1000,33000, 200,375,395);
		//h2_tdiff_R[i]	= new TH2F(Form("h2_tdiff_R_%i",i),";ADC R; TDC R - FADC R [ns]",			nADCbins,1000,33000, 200,375,395);

		//PrettyTH2F(h2_tdc_adc_L[i]);
		//PrettyTH2F(h2_tdc_adc_R[i]);
		//PrettyTH2F(h2_ftdc_adc_L[i]);
		//PrettyTH2F(h2_ftdc_adc_R[i]);
	}

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	//Read Dictionary of Hipo File  // new hipo4
        hipo::dictionary  factory;      // new hipo4
        reader.readDictionary(factory); // new hipo4
        //factory.show();               // new hipo4

	//hipo::bank BAND_ADC  (factory.getSchema("BAND::adc"  ));
	//hipo::bank BAND_TDC  (factory.getSchema("BAND::tdc"  ));
	//hipo::bank RUN_config(factory.getSchema("RUN::config"));
	BBand	band_hits(factory.getSchema("BAND::hits"));

	//One also needs a hipo::event object which is called from the reader for each event to get
        //the information for each bank
        hipo::event readevent;  // new hipo4

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		//Reader has to load information about event in hipo::event class
                reader.read(readevent); // new hipo4

                //Load explicitly all information for each bank for the event
                readevent.getStructure(band_hits);   // new hipo4

                //Now everything is loaded and can be used as before with HIPO3 files. There is only one difference:
                //The number of hits in each bank is determined by the function "getRows()" and not by "getSize" as before.

		if(event_counter%10000==0) cout << "event: " << event_counter << endl;
		event_counter++;	

		//BAND_ADC.show();
		//BAND_TDC.show();
		//RUN_config.show();

		int nHits = band_hits.getRows();
		// Skip events that don't have laser data (cosmics)
		if( nHits != 135 ) continue;
		//if( nHits < 2 || nHits > 10 ) continue;
		double diff_timesL[600] = {0.};
		double diff_timesR[600] = {0.};
		double mean_Ttdc[600] 	= {0.};
		double mean_Tfadc[600] 	= {0.};
		double TtdcL[600] 	= {0.};
		double TtdcR[600] 	= {0.};
		double TfadcL[600] 	= {0.};
		double TfadcR[600] 	= {0.};

		double adcL[600] = {0.};
		double adcR[600] = {0.};
		for( int hit = 0 ; hit < nHits ; hit++){
			int sector		= band_hits.getSector(hit);
			int layer		= band_hits.getLayer(hit);
			int component		= band_hits.getComponent(hit);
			float meanTimeTdc	= band_hits.getMeantimeTdc(hit);
			float meanTimeFadc	= band_hits.getMeantimeFadc(hit);
			float difftimeTdc	= band_hits.getDifftimeTdc(hit);
			float difftimeFadc	= band_hits.getDifftimeFadc(hit);
			float adcLcorr		= band_hits.getAdcLcorr(hit); 
			float adcRcorr		= band_hits.getAdcRcorr(hit);
			float TfadcLcorr	= band_hits.getTfadcLcorr(hit);
			float TfadcRcorr	= band_hits.getTfadcRcorr(hit);
			float TtdcLcorr		= band_hits.getTtdcLcorr(hit);
			float TtdcRcorr		= band_hits.getTtdcRcorr(hit);
			float x			= band_hits.getX(hit);
			float y			= band_hits.getY(hit);
			float z			= band_hits.getZ(hit);
			int ID 			= band_hits.getBarKey(hit);

			if( ID == 251 || ID == 252 ){
				if( ID == 251 ){
					//if( fabs(difftimeFadc)>.75 || fabs(difftimeFadc)>.75 ) continue;
					//if( sqrt(adcLcorr*adcRcorr) < 12000 ) continue;
				}
	
				//if( TtdcLcorr - TfadcLcorr < 1215 ) continue;
				//if( TtdcRcorr - TfadcRcorr < 1215 ) continue;

				// Looking at Ttdc - Reference and Tfadc - Reference for both left and right sides
				TtdcL[ID] 	= TtdcLcorr;
				TtdcR[ID]	= TtdcRcorr;
				TfadcL[ID]	= TfadcLcorr;
				TfadcR[ID]	= TfadcRcorr;
				adcL[ID]	= adcLcorr;
				adcR[ID]	= adcRcorr;
			}
		}// end read-in event
		for( int ID = 0 ; ID < 600 ; ID++){
			// Left PMTs
			int sector = round(ID/100);
			int layer = round((ID-sector*100)/10);
			int component = ID-sector*100-layer*10;

			// Reference based on one that reaches the largest ADC
			int refS = 2;
			int refL = 5;
			int refC = 1;
			int refID = refS*100 + refL*10 + refC;

			double refT_tdc  = 	( TtdcL[refID]  + TtdcR[refID]  )/2.;
			double refT_fadc = 	( TfadcL[refID] + TfadcR[refID] )/2.;

			if( refT_tdc == 0. || refT_fadc == 0. || TtdcL[ID] == 0. || TtdcR[ID] == 0. || TfadcL[ID] == 0. || TfadcR[ID] == 0. ) continue;
			hasEvents[ID]++;

			//cout << TtdcL[ID] - refT_tdc << " " << TtdcR[ID] - refT_tdc << " " << TfadcL[ID] - refT_fadc << " " << TfadcR[ID] - refT_fadc << "\n";

			// TDC_L - Ref bar time (no ref ADC cut)
			h2_tdc_adc_L[ID]	-> Fill( adcL[ID] , TtdcL[ID] - refT_tdc );
			h2_tdc_adc_L[ID] 	-> SetTitle(Form("Left TDC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
			// TDC_R - Ref bar time (no ref ADC cut)
			h2_tdc_adc_R[ID]	-> Fill( adcR[ID] , TtdcR[ID] - refT_tdc );
			h2_tdc_adc_R[ID] 	-> SetTitle(Form("Right TDC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
			// FADC_L - Ref bar time (no ref ADC cut)
			h2_ftdc_adc_L[ID]	-> Fill( adcL[ID] , TfadcL[ID] - refT_fadc );
			h2_ftdc_adc_L[ID] 	-> SetTitle(Form("Left FADC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
			// FADC_R - Ref bar time (no ref ADC cut)
			h2_ftdc_adc_R[ID]	-> Fill( adcR[ID] , TfadcR[ID] - refT_fadc );
			h2_ftdc_adc_R[ID] 	-> SetTitle(Form("Right FADC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
			// TDC - Ref bar time (no ref ADC cut)
			h2_tdcMean[ID]		-> Fill( sqrt(adcL[ID]*adcR[ID]) , (TtdcL[ID]+TtdcR[ID])/2. - refT_tdc );
			h2_tdcMean[ID] 		-> SetTitle(Form("TDC Mean, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
			// FADC - Ref bar time (no ref ADc cut)
			h2_ftdcMean[ID]		-> Fill( sqrt(adcL[ID]*adcR[ID]) , (TfadcL[ID]+TfadcR[ID])/2. - refT_fadc );
			h2_ftdcMean[ID] 	-> SetTitle(Form("FADC Mean, Sector:%i, Layer:%i, Component:%i",sector,layer,component));

			//cout << TtdcL[ID] - TfadcL[ID] << " " << TtdcR[ID] - TfadcR[ID] << "\n";

			// Correct myself for the overall 24ns trigger jitter
			//if( TtdcL[ID] - TfadcL[ID] < 1215 ) continue;
			//if( TtdcR[ID] - TfadcR[ID] < 1215 ) continue;

			// TDiff for L pmt
			h2_tdiff_L[ID]		-> Fill( adcL[ID]	, TtdcL[ID] - TfadcL[ID] );
			h2_tdiff_L[ID] 		-> SetTitle(Form("TDiff Left, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
			// TDiff for R pmt
			h2_tdiff_R[ID]		-> Fill( adcR[ID]	, TtdcR[ID] - TfadcR[ID] );
			h2_tdiff_R[ID] 		-> SetTitle(Form("TDiff Right, Sector:%i, Layer:%i, Component:%i",sector,layer,component));

			double refADC = sqrt( adcL[refID] * adcR[refID] );
			// Cut on high ADC
			if( fabs( refADC - 20000 ) < 100 ){
				// TDC_L - Ref bar time (ref ADC cut)
				h2_tdc_adc_L_Cut[ID] 	-> Fill( adcL[ID] , TtdcL[ID] - refT_tdc 	); 
				h2_tdc_adc_L_Cut[ID]	-> SetTitle(Form("Cut, Left TDC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				// TDC_R - Ref bar time (ref ADC cut)
				h2_tdc_adc_R_Cut[ID] 	-> Fill( adcR[ID] , TtdcR[ID] - refT_tdc 	); 
				h2_tdc_adc_R_Cut[ID]  	-> SetTitle(Form("Cut, Right TDC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				// FADC_L - Ref bar time (ref ADC cut)
				h2_ftdc_adc_L_Cut[ID]	-> Fill( adcL[ID] , TfadcL[ID] - refT_fadc 	);
				h2_ftdc_adc_L_Cut[ID]	-> SetTitle(Form("Cut, Left FADC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				// FADC_R - Ref bar time (ref ADC cut)
				h2_ftdc_adc_R_Cut[ID]	-> Fill( adcR[ID] , TfadcR[ID] - refT_fadc 	);
				h2_ftdc_adc_R_Cut[ID] 	-> SetTitle(Form("Cut, Right FADC, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				// TDC - Ref bar time (ref ADC cut)
				h2_tdcMean_Cut[ID]	-> Fill( sqrt(adcL[ID]*adcR[ID]) , (TtdcL[ID]+TtdcR[ID])/2. - refT_tdc );
				h2_tdcMean_Cut[ID] 	-> SetTitle(Form("TDC Mean, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
				// FADC - Ref bar time (ref ADc cut)
				h2_ftdcMean_Cut[ID]	-> Fill( sqrt(adcL[ID]*adcR[ID]) , (TfadcL[ID]+TfadcR[ID])/2. - refT_fadc );
				h2_ftdcMean_Cut[ID] 	-> SetTitle(Form("FADC Mean, Sector:%i, Layer:%i, Component:%i",sector,layer,component));
			}
		}// end loop over all bars for single event

	}// end file

	// Now that we have our histograms full, we want to, for each histogram, take projections onto Y for a given slice in X
	// to make a TGraph for each PMT
	TApplication theApp("App",&argc,argv);

	TFile * outFile = new TFile("testing.root","RECREATE");
	outFile->cd();

	TCanvas *c1 = new TCanvas("c1","Testing",900,900);
	c1->Divide(2,5);

	std::vector<double>	adcs,adcsErr,times,timesErr,res,resErr;
	double par0, par1;
	int ADCcut = 25000;

	// Draw TDC Mean Plot
	c1->cd(1);
	TF1 * model = new TF1("timeWalk",wlk,0,33000,2);
	h2_tdcMean[252]->Fit("timeWalk","QES");
	h2_tdcMean[252]->Draw("COLZ");
	// Fit the TDC mean plot
	walkCorr(	&adcs		,
			&adcsErr	,
			&times		,
			&timesErr	,
			&res		,
			&resErr		,
			ADCcut		,
			h2_tdcMean[252]	);
	// Draw the fit of the projections
	fitGraph(adcs, times, adcsErr, timesErr, par0, par1 ,c1, 2, Form("Test Bar %i",252) );
	// Draw the resulting fit on the original histogram
	TF1 * model2 = new TF1("timeWalk2",wlk,0,33000,2);	
	model2->SetParameter(0,par0);
	model2->SetParameter(1,par1);
	model2->SetLineColor(4);
	c1->cd(1);
	model2->Draw("same");
	c1->Update();


	adcs.clear(); adcsErr.clear(); times.clear(); timesErr.clear(); res.clear(); resErr.clear();
	// Draw the TDIFF_L Plot
	c1->cd(3);
	h2_tdiff_L[252]->Draw("COLZ");
	// Fit the TDIFF_L Plot
	walkCorr(       &adcs           ,
                        &adcsErr        ,
                        &times          ,
                        &timesErr       ,
                        &res            ,
                        &resErr         ,
                        ADCcut          ,
                        h2_tdiff_L[252] );
	// Draw the fit of the projections
	fitGraph(adcs, times, adcsErr, timesErr, par0, par1 ,c1, 4, Form("Test Bar %i",252) );
	TF1 * model3 = new TF1("timeWalk3",wlk,0,33000,2);
	model3->SetParameter(0,par0);
	model3->SetParameter(1,par1);
	model3->SetLineColor(4);
	c1->cd(3);
	model3->Draw("same");
	c1->Update();
	
	
	adcs.clear(); adcsErr.clear(); times.clear(); timesErr.clear(); res.clear(); resErr.clear();
	// Draw the TDIFF_R Plot
	c1->cd(5);
	h2_tdiff_R[252]->Draw("COLZ");
	// Fit the TDIFF_R Plot
	walkCorr(       &adcs           ,
                        &adcsErr        ,
                        &times          ,
                        &timesErr       ,
                        &res            ,
                        &resErr         ,
                        ADCcut          ,
                        h2_tdiff_R[252] );
	// Draw the fit of the projecitons
	fitGraph(adcs, times, adcsErr, timesErr, par0, par1 ,c1, 6, Form("Test Bar %i",252) );
	TF1 * model4 = new TF1("timeWalk4",wlk,0,33000,2);
	model4->SetParameter(0,par0);
	model4->SetParameter(1,par1);
	model4->SetLineColor(4);
	c1->cd(5);
	model4->Draw("same");
	c1->Update();


	adcs.clear(); adcsErr.clear(); times.clear(); timesErr.clear(); res.clear(); resErr.clear();
	// Draw the TDC_L - Ref Plot
	c1->cd(7);
	h2_tdc_adc_L[252]->Draw("COLZ");
	// Fit the TDC_L - Ref Plot
	walkCorr(       &adcs           ,
                        &adcsErr        ,
                        &times          ,
                        &timesErr       ,
                        &res            ,
                        &resErr         ,
                        ADCcut          ,
                        h2_tdc_adc_L[252] );
	// Draw the fit of the projecitons
	fitGraph(adcs, times, adcsErr, timesErr, par0, par1 ,c1, 8, Form("Test Bar %i",252) );
	TF1 * model5 = new TF1("timeWalk5",wlk,0,33000,2);
        model5->SetParameter(0,par0);
        model5->SetParameter(1,par1);
        model5->SetLineColor(4);
        c1->cd(7);
        model5->Draw("same");
        c1->Update();


	adcs.clear(); adcsErr.clear(); times.clear(); timesErr.clear(); res.clear(); resErr.clear();
	// Draw the TDC_R - Ref Plot
	c1->cd(9);
	h2_tdc_adc_R[252]->Draw("COLZ");	
	walkCorr(       &adcs           ,
                        &adcsErr        ,
                        &times          ,
                        &timesErr       ,
                        &res            ,
                        &resErr         ,
                        ADCcut          ,
                        h2_tdc_adc_R[252] );
	// Draw the fit of the projecitons
	fitGraph(adcs, times, adcsErr, timesErr, par0, par1 ,c1, 10, Form("Test Bar %i",252) );
        TF1 * model6 = new TF1("timeWalk5",wlk,0,33000,2);
        model6->SetParameter(0,par0);
        model6->SetParameter(1,par1);
        model6->SetLineColor(4);
        c1->cd(9);
        model6->Draw("same");
        c1->Update();

	
	outFile->Close();
	
	theApp.Run();

	// -------------------------------------------------------------------------------------------------
	// Fitting functions
	/*
	   for(int i = 0 ; i < nHistos ; i++){
	   cout << "fitting " << i << "\n";
	   int notEmpty = (h2_tdc_adc_L[i]->Integral()+h2_tdc_adc_R[i]->Integral());
	   if(notEmpty){
	   TF1 * f_timewalk_L = new TF1("f_timewalk_L","[0]/TMath::Sqrt(x)+[1]");
	   TF1 * f_timewalk_R = new TF1("f_timewalk_R","[0]/TMath::Sqrt(x)+[1]");

	   h2_ff_Ltdc_adc_L[i] -> Fit("f_timewalk_L","Q");
	   h2_tdc_adc_R[i] -> Fit("f_timewalk_R","Q");

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
	   */
	// -------------------------------------------------------------------------------------------------
	// Printing results onto canvases
	/*
	TCanvas * c0 = new TCanvas("c0");
	c0 -> Divide(2,1);
	//c0 -> cd(1);	h2_tdc_adc_tot_L -> Draw("COLZ");
	//c0 -> cd(2);	h2_tdc_adc_tot_R -> Draw("COLZ");
	c0 -> Modified();
	c0 -> Update();

	TCanvas **** cSLC = new TCanvas***[5];

	for(int is = 0 ; is < 5 ; is++){
		cSLC[is] = new TCanvas**[6];
		for(int il = 0 ; il < 5 ; il++){
			cSLC[is][il] = new TCanvas*[ slc[il][is] ];

			for( int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++ ){
				cSLC[is][il][cIdx] = new TCanvas(Form("S%iL%iC%i",is,il,cIdx),Form("Sector %i, Layer %i, Component %i",is+1,il+1,cIdx+1),900,900);
				cSLC[is][il][cIdx]->Divide(2,4);

				int ID = 100*(is+1)+10*(il+1)+(cIdx+1);
				int min,max;

				// Draw TDC_L - Ref with no cut
				cSLC[is][il][cIdx]->cd(1);
				min 	= 	getMinNonEmptyBin(h2_tdc_adc_L[ID]);
				max 	=	getMaxNonEmptyBin(h2_tdc_adc_L[ID]);
				h2_tdc_adc_L[ID]-> GetYaxis() -> SetRange(min,max);
				h2_tdc_adc_L[ID]	->	Draw("COLZ");
				// Draw TDC_L - Ref with ref adc cut
				//cSLC[is][il][cIdx]->cd(2);
				//min     =       getMinNonEmptyBin(h2_tdc_adc_L_Cut[ID]);
				//max     =       getMaxNonEmptyBin(h2_tdc_adc_L_Cut[ID]);
				//h2_tdc_adc_L_Cut[ID]-> GetYaxis() -> SetRange(min,max);
				//h2_tdc_adc_L_Cut[ID]	->	Draw("COLZ");

				// Draw TDC_R - Ref with no cut
				cSLC[is][il][cIdx]->cd(2);
				min     =       getMinNonEmptyBin(h2_tdc_adc_R[ID]);
				max     =       getMaxNonEmptyBin(h2_tdc_adc_R[ID]);
				h2_tdc_adc_R[ID]-> GetYaxis() -> SetRange(min,max);
				h2_tdc_adc_R[ID]	->	Draw("COLZ");
				// Draw TDC_R - Ref with ref adc cut
				//cSLC[is][il][cIdx]->cd(4);
				//min     =       getMinNonEmptyBin(h2_tdc_adc_R_Cut[ID]);
				//max     =       getMaxNonEmptyBin(h2_tdc_adc_R_Cut[ID]);
				//h2_tdc_adc_R_Cut[ID]-> GetYaxis() -> SetRange(min,max);
				//h2_tdc_adc_R_Cut[ID]	->	Draw("COLZ");

				// Draw FADC_L - Ref with no cut
				cSLC[is][il][cIdx]->cd(3);
				min     =       getMinNonEmptyBin(h2_ftdc_adc_L[ID]);
				max     =       getMaxNonEmptyBin(h2_ftdc_adc_L[ID]);
				h2_ftdc_adc_L[ID]-> GetYaxis() -> SetRange(min,max);
				h2_ftdc_adc_L[ID]        ->      Draw("COLZ");
				// Draw FADC_L - Ref with ref adc cut
				//cSLC[is][il][cIdx]->cd(6);
				//min     =       getMinNonEmptyBin(h2_ftdc_adc_L_Cut[ID]);
				//max     =       getMaxNonEmptyBin(h2_ftdc_adc_L_Cut[ID]);
				//h2_ftdc_adc_L_Cut[ID]-> GetYaxis() -> SetRange(min,max);
				//h2_ftdc_adc_L_Cut[ID]        ->      Draw("COLZ");

				// Draw FADC_R - Ref with no cut
				cSLC[is][il][cIdx]->cd(4);
				min     =       getMinNonEmptyBin(h2_ftdc_adc_R[ID]);
				max     =       getMaxNonEmptyBin(h2_ftdc_adc_R[ID]);
				h2_ftdc_adc_R[ID]-> GetYaxis() -> SetRange(min,max);
				h2_ftdc_adc_R[ID]        ->      Draw("COLZ");
				// Draw FADC_R - Ref with ref adc cut
				//cSLC[is][il][cIdx]->cd(8);
				//min     =       getMinNonEmptyBin(h2_ftdc_adc_R_Cut[ID]);
				//max     =       getMaxNonEmptyBin(h2_ftdc_adc_R_Cut[ID]);
				//h2_ftdc_adc_R_Cut[ID]-> GetYaxis() -> SetRange(min,max);
				//h2_ftdc_adc_R_Cut[ID]        ->      Draw("COLZ");


				// Draw TDC mean time - Ref with no cut
				cSLC[is][il][cIdx]->cd(5);
				min     =       getMinNonEmptyBin(h2_tdcMean[ID]);
				max     =       getMaxNonEmptyBin(h2_tdcMean[ID]);
				h2_tdcMean[ID]->GetYaxis() -> SetRange(min,max);
				h2_tdcMean[ID]		->	Draw("COLZ");
				// Draw TDC mean time - Ref with ref adc cut
				//cSLC[is][il][cIdx]->cd(10);
				//min     =       getMinNonEmptyBin(h2_tdcMean_Cut[ID]);
				//max     =       getMaxNonEmptyBin(h2_tdcMean_Cut[ID]);
				//h2_tdcMean_Cut[ID] -> GetYaxis() -> SetRange(min,max);
				//h2_tdcMean_Cut[ID]	->	Draw("COLZ");

				// Draw FADC mean time - Ref with no cut
				cSLC[is][il][cIdx]->cd(6);
				min     =       getMinNonEmptyBin(h2_ftdcMean[ID]);
				max     =       getMaxNonEmptyBin(h2_ftdcMean[ID]);
				h2_ftdcMean[ID] -> GetYaxis() -> SetRange(min,max);
				h2_ftdcMean[ID]		->	Draw("COLZ");
				// Draw FADC mean time - Ref with ref adc cut
				//cSLC[is][il][cIdx]->cd(12);
				//min     =       getMinNonEmptyBin(h2_ftdcMean_Cut[ID]);
				//max     =       getMaxNonEmptyBin(h2_ftdcMean_Cut[ID]);
				//h2_ftdcMean_Cut[ID] -> GetYaxis() -> SetRange(min,max);
				//h2_ftdcMean_Cut[ID]	->	Draw("COLZ");

				// Draw TDIFF for left PMT
				cSLC[is][il][cIdx]->cd(7);
				min     =       getMinNonEmptyBin(h2_tdiff_L[ID]);
				max     =       getMaxNonEmptyBin(h2_tdiff_L[ID]);
				h2_tdiff_L[ID] -> GetYaxis() -> SetRange(min,max);
				h2_tdiff_L[ID]		->	Draw("COLZ");
				// Draw TDIFF for right PMT
				cSLC[is][il][cIdx]->cd(8);
				min     =       getMinNonEmptyBin(h2_tdiff_R[ID]);
				max     =       getMaxNonEmptyBin(h2_tdiff_R[ID]);
				h2_tdiff_R[ID] -> GetYaxis() -> SetRange(min,max);
				h2_tdiff_R[ID]		->	Draw("COLZ");


				cSLC[is][il][cIdx]->Modified(); cSLC[is][il][cIdx]->Update();
			}// end loop over component
			*/
			/*	
			//cSLC[is][il] = new TCanvas(Form("S%iL%i",is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);
			//cSLC[is][il] -> Divide(2,7);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
			cout << "\t plotting " << is << " " << il << " " << cIdx <<"\n";
			int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
			int notEmpty = (h2_tdc_adc_L[identifier]->Integral()+h2_tdc_adc_R[identifier]->Integral());
			if(notEmpty){
			int min, max;

			cSLC[is][il] -> cd(2*cIdx+1);
			gPad -> SetBottomMargin(0.26);
			min = getMinNonEmptyBin(h2_tdc_adc_L[identifier]);
			max = getMaxNonEmptyBin(h2_tdc_adc_L[identifier]);
			h2_tdc_adc_L[identifier] -> GetYaxis() -> SetRange(min,max);
			h2_tdc_adc_L[identifier] -> Draw("COLZ");

			cSLC[is][il] -> cd(2*cIdx+2);
			gPad -> SetBottomMargin(0.26);
			min = getMinNonEmptyBin(h2_tdc_adc_R[identifier]);
			max = getMaxNonEmptyBin(h2_tdc_adc_R[identifier]);
			h2_tdc_adc_R[identifier] -> GetYaxis() -> SetRange(min,max);
			h2_tdc_adc_R[identifier] -> Draw("COLZ");
			}
			}
			cSLC[is][il] -> Modified();	cSLC[is][il] -> Update();
			*/
		//}// end loop over layer
	//}// end loop over sector

	// -------------------------------------------------------------------------------------------------
	// Saving fit values to ccdb tables
	/*
	   ofstream tabL, tabR;
	   tabL.open("timeWalkPar_L.txt");
	   tabR.open("timeWalkPar_R.txt");

	   for(int is = 1 ; is <= 5 ; is++){	
	   for(int il = 1 ; il <= 6 ; il++){
	   for(int ic = 1 ; ic <= slc[il-1][is-1] ; ic++){
	   int idx = 100*is + 10*il + ic;

	//if(parL[idx][0]!=0){
	tabL << is << "\t" << il << "\t" << ic << "\t" << "\t";
	if(il==6) tabL << "0\t0\t0\t0" << endl;
	else tabL << parL[idx][0] << "\t" << parL[idx][1] << "\t" << parL[idx][2]  << "\t" << parL[idx][3] << endl;
	//}
	// ---
	//if(parR[idx][0]!=0){
	tabR << is << "\t" << il << "\t" << ic << "\t" << "\t";
	if(il==6) tabR << "0\t0\t0\t0" << endl;
	else tabR << parR[idx][0] << "\t" << parR[idx][1] << "\t" << parR[idx][2]  << "\t" << parR[idx][3] << endl;
	//}
	}
	}
	}
	tabL.close();
	tabR.close();
	*/
	// -------------------------------------------------------------------------------------------------
	// Deleting empty histograms
	/*
	   for(int i = 0 ; i < nHistos ; i++){
	   int notEmpty = (h2_tdc_adc_L[i]->Integral()+h2_tdc_adc_R[i]->Integral());
	   if(!notEmpty){
	   delete h2_tdc_adc_L[i];
	   delete h2_tdc_adc_R[i];
	   }
	   }

	// -------------------------------------------------------------------------------------------------
	// Saving plots to a pdf file
	*/
	/*
	c0 -> Print("results_timeWalk.pdf(");
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			for( int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				cSLC[is][il][cIdx] -> Print("results_timeWalk.pdf");
			}
		}
	}
	c0 -> Print("results_timeWalk.pdf)");
	//myapp -> Run();
	*/
	return 0;
}
// ========================================================================================================================================
void PrettyTH2F(TH2F * h2) {
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetYaxis() -> CenterTitle();

	h2 -> SetTitleSize(0.5);

	//h2 -> GetXaxis() -> SetLabelSize(0.13);
	//h2 -> GetYaxis() -> SetLabelSize(0.13);
	//h2 -> GetXaxis() -> SetTitleSize(0.13);
	//h2 -> GetYaxis() -> SetTitleSize(0.13);

	//h2 -> GetYaxis() -> SetTitleOffset(0.38);
	h2 -> GetYaxis() -> SetNdivisions(509);

	h2 -> GetXaxis() -> SetNdivisions(509);
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


void walkCorr(	std::vector<double> *widths		,
		std::vector<double> *widthsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2F * hist				){
	int currBin = 0;
	while( currBin < nADCbins){
		cout << "\tCurrently on bin: " << currBin << "\n";
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , true , 0, xPt, yPt, yEr, ySig, ySigEr );
		cout << "\tAfter projection: " << xPt << " " << yPt << "\n";
		currBin += step ;
		if( xPt > widthCut) continue;
		widths		->push_back(xPt);
		widthsErr	->push_back(0);
		times		->push_back(yPt);
		timesErr	->push_back(yEr);
		res		->push_back(ySig);
		resErr		->push_back(ySigEr);
	}
}

int doProj( TH2F * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE){
	int cnt = 0;
	int step = 0;
	TCanvas * trash = new TCanvas("trash");
	int thres = hist->GetEntries() / 5;
	while ( cnt < thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
		if( (bin+step) >= nADCbins) break;
		step+=1;
	}
	double scale = (33000.-1000.)/400.;
	cout << "\t\tTaking projection between " << bin*scale << " " << (bin+step)*scale <<"\n";
	char temp[100];
	sprintf(temp,"slice_%d_%d",flag,bin);
	TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
	
	// Getting the mean x value in this range:
	hist->GetXaxis()->SetRange(bin,bin+step);
	x = hist->GetMean(1);
	hist->GetXaxis()->SetRange();
	TFitResultPtr f = pj->Fit("gaus","QES");
	y = f->Parameter(1);
	yE = f->ParError(1);
	sig = f->Parameter(2);
	sigE = f->ParError(2);
	
	if( write ){
		pj->Write();
	}

	delete trash;
	delete pj;

	return step;
}

void fitGraph( 	std::vector<double> widths	,
		std::vector<double> times	,
		std::vector<double> widthsErr	,
		std::vector<double> timesErr	,
		double &par0			,
		double &par1			,
		TCanvas *c			,
		int cd				,
		TString name			){
	
	int dim = widths.size();
	TGraphErrors *g = new TGraphErrors(dim, &widths[0], &times[0], &widthsErr[0], &timesErr[0]);
	TF1 * model = new TF1("timeWalk",wlk,0,33000,2);
	model->SetLineColor(4);
	TFitResultPtr ptr = g->Fit(model,"QES");
	par0 = ptr->Parameter(0);
	par1 = ptr->Parameter(1);

	c->cd(cd);
	g->SetMarkerStyle(20);
	//g->GetXaxis()->SetLimits(17,24);
	g->Draw("AP");
	g->SetTitle("Timewalk of "+name);
	g->GetXaxis()->SetTitle("ADC of PMT [a.u.]");
	g->GetYaxis()->SetTitleOffset(1.6);
	g->GetYaxis()->SetTitle(name+" Time - Ref Time [ns]");
	g->GetXaxis()->SetNdivisions(509);
	c->Update();

}
double wlk( double *x , double *p){
	double var = *x;
	return p[0] + p[1] / sqrt(var);
}
