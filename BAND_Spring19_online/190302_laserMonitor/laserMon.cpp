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
#include "TGraphErrors.h"

#include "reader.h"
#include "node.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

using namespace std;

// Define struct for each pmt to store the ADC and TDC per event
struct Bar{
	int id;
	int sector, layer, component;
	std::vector<double> adcL,adcR;
	std::vector<double> meanT_tdc,meanT_ftdc;
	std::vector<double> diffT_tdc,diffT_ftdc;
	std::vector<double> events;
};
// Define operator to search for pmt ID in our vector later on
struct find_bar{
	int uniqueID;
	find_bar(int id): uniqueID(id){ }
	bool operator()( Bar const& bar) const{
		return bar.id == uniqueID;
	}
};

void PrettyTGraph( TGraphErrors * h , int ID, int s, int l, int c);

int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// ========================================================================================================================================
int main(int argc, char** argv) {

#ifdef WITHRINT
	//TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	//TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle->SetTitleSize(0.1,"t");
	gStyle->SetOptStat(0);

	if( argc != 2){
		cerr << "Wrong number of arguments given. Instead use:\n\t./laserMon [inputFile]\n";
		return -1;
	}

	TString inputFile = argv[1];

	// Declaring histograms
	const int nHistos = 600;
	TGraphErrors ** evNo_adc = new TGraphErrors * [nHistos];
	TGraphErrors ** evNo_ftdc = new TGraphErrors * [nHistos];
	TGraphErrors ** evNo_tdc = new TGraphErrors * [nHistos];

	// Define our PMT storage vector
	std::vector< Bar > Bars;
	int refBar = 251;
	Bar RefBar;
	RefBar.id = refBar;
	RefBar.sector = 2;
	RefBar.layer = 5;
	RefBar.component = 1;

	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	// Prep reader for BAND banks
	//hipo::bank BAND_HITS("BAND::hits",reader);
	BBand	band_hits("BAND::hits"	,	reader);
	//hipo::bank BAND_ADC("BAND::adc",reader);
	//hipo::bank BAND_TDC("BAND::tdc",reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000==0) cout << "event: " << event_counter << endl;
		event_counter++;	

		int nHits = band_hits.getSize();
		// Skip events that don't have laser data (cosmics)
		if( nHits != 135 ) continue;
		//band_hits.show();
		for( int hit = 0 ; hit < nHits ; hit++){
			/*
			int sector		= BAND_HITS.getInt(1,nHits);
			int layer		= BAND_HITS.getInt(2,nHits);
			int component		= BAND_HITS.getInt(3,nHits);
			float meanTimeTdc	= BAND_HITS.getFloat(4,nHits);	// this has been corrected for paddle offsets
			float meanTimeFadc	= BAND_HITS.getFloat(5,nHits);	// this has been corrected for paddle offsets
			float difftimeTdc	= BAND_HITS.getFloat(6,nHits);	// this has been corrected for L-R offset
			float difftimeFadc	= BAND_HITS.getFloat(7,nHits);	// this has been corrected for L-R offset
			float adcLcorr		= BAND_HITS.getFloat(8,nHits);	// this has been corrected for attenuation
			float adcRcorr		= BAND_HITS.getFloat(9,nHits);	// this has been corrected for atteunation
			float x			= BAND_HITS.getFloat(14,nHits);
			float y			= BAND_HITS.getFloat(15,nHits);
			float z			= BAND_HITS.getFloat(16,nHits);
			*/
			int sector		= band_hits.getSector(hit);
			int layer		= band_hits.getLayer(hit);
			int component		= band_hits.getComponent(hit);
			float meanTimeTdc	= band_hits.getMeantimeTdc(hit);
			float meanTimeFadc	= band_hits.getMeantimeFadc(hit);
			float difftimeTdc	= band_hits.getDifftimeTdc(hit);
			float difftimeFadc	= band_hits.getDifftimeFadc(hit);
			float adcLcorr		= band_hits.getAdcLcorr(hit); 
			float adcRcorr		= band_hits.getAdcRcorr(hit);
			float x			= band_hits.getX(hit);
			float y			= band_hits.getY(hit);
			float z			= band_hits.getZ(hit);
			int ID 			= band_hits.getBarKey(hit);


			if( ID == refBar && event_counter % 50 == 0 ){
				RefBar.meanT_tdc.push_back(		(double)	meanTimeTdc	);
				RefBar.meanT_ftdc.push_back(		(double)	meanTimeFadc	);
				RefBar.diffT_tdc.push_back(		(double)	difftimeTdc	);
				RefBar.diffT_ftdc.push_back(		(double)	difftimeFadc	);
				RefBar.events.push_back( 		(double)	event_counter 	);
				continue;
			}
		
			std::vector< Bar >::iterator it;
			it = std::find_if( Bars.begin(), Bars.end(), find_bar( ID ) );
			// If I've already started saving this bar
			if( it != Bars.end() ){
				if( event_counter % 50 == 0 ){
					int idx = std::distance( Bars.begin() , it );
					Bars.at(idx).adcL.push_back( 		(double)	adcLcorr	);
					Bars.at(idx).adcR.push_back( 		(double)	adcRcorr	);
					Bars.at(idx).meanT_tdc.push_back(	(double)	meanTimeTdc	);
					Bars.at(idx).meanT_ftdc.push_back(	(double)	meanTimeFadc	);
					Bars.at(idx).diffT_tdc.push_back(	(double)	difftimeTdc	);
					Bars.at(idx).diffT_ftdc.push_back(	(double)	difftimeFadc	);
					Bars.at(idx).events.push_back( 		(double)	event_counter 	);
				}
			}
			else{	// Otherwise create and add new pmt
				//cout << "\t\tFirst hit for this PMT!\n"
				if( event_counter % 50 == 0 ){
					Bar thisBar;

					thisBar.sector = sector;
					thisBar.layer = layer;
					thisBar.component = component;
					thisBar.id = ID;

					thisBar.adcL.push_back( 		(double)	adcLcorr	);
					thisBar.adcR.push_back( 		(double)	adcRcorr	);
					thisBar.meanT_tdc.push_back(		(double)	meanTimeTdc	);
					thisBar.meanT_ftdc.push_back(		(double)	meanTimeFadc	);
					thisBar.diffT_tdc.push_back(		(double)	difftimeTdc	);
					thisBar.diffT_ftdc.push_back(		(double)	difftimeFadc	);
					thisBar.events.push_back( 		(double)	event_counter 	);

					Bars.push_back( thisBar );
				}
			}
		}
	}// end file

	// Create canvases for each layer and each sector
	TCanvas *** cSLC = new TCanvas**[5];
	for( int sector = 0 ; sector < 5 ; sector++ ){
		cSLC[sector] = new TCanvas*[6];
		for( int layer = 0 ; layer < 6 ; layer++ ){
			cSLC[sector][layer] = new TCanvas(Form("S%iL%i",sector,layer),
									Form("Sector %i, Layer %i",sector+1,layer+1),900,900);
			cSLC[sector][layer] -> Divide(3,7);
		}
		//break;
	}


	// Fill corresponding canvases
	for( int i = 0 ; i < Bars.size() ; i++){

		Bar thisBar = Bars.at(i);

		// Subtract off a reference time for each bar
		std::vector<double> meanT_ftdc;
		std::vector<double> meanT_tdc;
		std::vector<double> gMean;
		for( int j = 0 ; j < RefBar.meanT_ftdc.size() ; j++){
			if( fabs(thisBar.meanT_ftdc[j] - RefBar.meanT_ftdc[j]) > 100 ) continue;
			meanT_ftdc.push_back( thisBar.meanT_ftdc[j] - RefBar.meanT_ftdc[j] );
			meanT_tdc.push_back(  thisBar.meanT_tdc[j] - RefBar.meanT_tdc[j] );
			gMean.push_back( sqrt( thisBar.adcL[j] * thisBar.adcR[j] ) );
		}

		int ID = thisBar.id;

		evNo_adc[ID]	= new TGraphErrors( thisBar.events.size() , &thisBar.events[0] , &gMean[0] , 0 , 0 );
		evNo_ftdc[ID] 	= new TGraphErrors( thisBar.events.size() , &thisBar.events[0] , &meanT_ftdc[0] , 0 , 0 );
		evNo_tdc[ID]	= new TGraphErrors( thisBar.events.size() , &thisBar.events[0] , &meanT_tdc[0] , 0 , 0 );

		PrettyTGraph( evNo_adc[ID] 	,	0, thisBar.sector, thisBar.layer, thisBar.component );
		PrettyTGraph( evNo_ftdc[ID] 	,	1, thisBar.sector, thisBar.layer, thisBar.component );
		PrettyTGraph( evNo_tdc[ID] 	,	2, thisBar.sector, thisBar.layer, thisBar.component );
		
		
		cSLC[thisBar.sector-1][thisBar.layer-1] -> cd( 3* (thisBar.component-1) + 1);
		gPad -> SetBottomMargin(0.26);
		evNo_adc[ID] -> Draw();
		
		cSLC[thisBar.sector-1][thisBar.layer-1] -> cd( 3* (thisBar.component-1) + 2);
		gPad -> SetBottomMargin(0.26);
		evNo_ftdc[ID] -> Draw();
	
		cSLC[thisBar.sector-1][thisBar.layer-1] -> cd( 3* (thisBar.component-1) + 3);
		gPad -> SetBottomMargin(0.26);
		evNo_tdc[ID] -> Draw();

		cSLC[thisBar.sector-1][thisBar.layer-1] -> Modified();
		cSLC[thisBar.sector-1][thisBar.layer-1] -> Update();
		
	}

	TCanvas * c0 = new TCanvas("c0");
	c0 -> Print("monitoring_plots.pdf(");
	for(int is = 0 ; is < 5 ; is++){
                for(int il = 0 ; il < 6 ; il++){
                        cSLC[is][il] -> Print("monitoring_plots.pdf");
                }
        }
	c0 -> Print("monitoring_plots.pdf)");
	//myapp -> Run();
	return 0;
}
void PrettyTGraph( TGraphErrors * h , int ID, int s, int l, int c){
	if( ID == 0){
		h->SetTitle(Form("Sector %i, Layer %i, Component %i",s,l,c));
		h->GetXaxis()->SetTitle("Event (time)");
		h->GetYaxis()->SetTitle("ADC from Laser");
		//h->GetHistogram()->SetMaximum(15000);
		//h->GetHistogram()->SetMinimum(0);
	}
	else if( ID == 1 ){
		h->SetTitle(Form("Sector %i, Layer %i, Component %i",s,l,c));
		h->GetXaxis()->SetTitle("Event (time)");
		h->GetYaxis()->SetTitle("ADC time [ns] from Laser");
		h->GetHistogram()->SetMaximum(2);
		h->GetHistogram()->SetMinimum(-2);
	}
	else{
		h->SetTitle(Form("Sector %i, Layer %i, Component %i",s,l,c));
		h->GetXaxis()->SetTitle("Event (time)");
		h->GetYaxis()->SetTitle("TDC time [ns] from Laser");
		h->GetHistogram()->SetMaximum(2);
		h->GetHistogram()->SetMinimum(-2);
	}
}
