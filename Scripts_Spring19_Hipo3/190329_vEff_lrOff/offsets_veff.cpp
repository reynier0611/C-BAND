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
#include "TLine.h"
#include "TProfile.h"
#include "TF1.h"
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
void PrettyTH1F(TH1F * h1,TString titx,TString tity,int color);
void PrettyTH2F(TH2F * h2,TString titx,TString tity);
int getMinNonEmptyBin(TH2F * h2);
int getMaxNonEmptyBin(TH2F * h2);

double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};
int slc[6][5] = {{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,6,6,2},{3,7,5,5,0},{3,7,6,6,2}};
// ========================================================================================================================================
int main(int argc, char** argv) {


	//TApplication *myapp = new TApplication("myapp",0,0);

	std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

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
        cout << "Run this code on cosmic or production data" << endl << endl;
        cout << "*************************************************************" << endl;

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
	TH2F ** h2_empty = new TH2F * [nHistos];
	TH2F ** h2_meantime = new TH2F * [nHistos];
	for( int i = 0 ; i < nHistos ; i++ ){
		h1_tdc_diff [i] = new TH1F(Form("h1_tdc_diff_%i" ,i),"",400,-20,20);
		h1_ftdc_diff[i] = new TH1F(Form("h1_ftdc_diff_%i",i),"",400,-20,20);
		h2_empty    [i] = new TH2F(Form("h2_empty_%i"    ,i),";TDC L-R;FADC L-R",400,-20,20,400,-20,20);
		h2_meantime [i] = new TH2F(Form("h2_meantime_%i" ,i),";#sqrt{ADC_{L}*ADC_{R}};#frac{L+R}{2} TDC - FADC",500,0,25000,1000,1045,1075);
	
		PrettyTH2F(h2_empty    [i],"","");
		PrettyTH2F(h2_meantime [i],"","");
	}

	double par_pad[nHistos][2] = {0};
	double val_slope = 0;
	int ctr_slope = 0;
	TH1D * h1_slopes       = new TH1D("h1_slopes"      ,";slope",100,0.8,1.5);	h1_slopes       -> SetLineColor( 1);	h1_slopes       -> SetLineWidth(2);	h1_slopes -> SetLineStyle(2);
	TH1D * h1_slopes_short = new TH1D("h1_slopes_short",";slope",100,0.8,1.5);	h1_slopes_short -> SetLineColor( 2);	h1_slopes_short -> SetLineWidth(2);
	TH1D * h1_slopes_long  = new TH1D("h1_slopes_long" ,";slope",100,0.8,1.5);	h1_slopes_long  -> SetLineColor(62);	h1_slopes_long  -> SetLineWidth(2);
	TH1D * h1_slopes_med   = new TH1D("h1_slopes_med"  ,";slope",100,0.8,1.5);	h1_slopes_med   -> SetLineColor( 8);	h1_slopes_med   -> SetLineWidth(2);

	TH1D * h1_vEff_tdc_s   = new TH1D("h1_vEff_tdc_s"  ,"short bars only;effective speed of light [cm/ns]" ,60,10,17);	h1_vEff_tdc_s  -> SetLineColor(1);
	TH1D * h1_vEff_fadc_s  = new TH1D("h1_vEff_fadc_s" ,"short bars only;effective speed of light [cm/ns]" ,60,10,17);	h1_vEff_fadc_s -> SetLineColor(2);
	TH1D * h1_vEff_tdc_m   = new TH1D("h1_vEff_tdc_m"  ,"medium bars only;effective speed of light [cm/ns]",60,10,17);	h1_vEff_tdc_m  -> SetLineColor(1);
        TH1D * h1_vEff_fadc_m  = new TH1D("h1_vEff_fadc_m" ,"medium bars only;effective speed of light [cm/ns]",60,10,17);	h1_vEff_fadc_m -> SetLineColor(2);
	TH1D * h1_vEff_tdc_l   = new TH1D("h1_vEff_tdc_l"  ,"long bars only;effective speed of light [cm/ns]"  ,60,10,17);	h1_vEff_tdc_l  -> SetLineColor(1);
        TH1D * h1_vEff_fadc_l  = new TH1D("h1_vEff_fadc_l" ,"long bars only;effective speed of light [cm/ns]"  ,60,10,17);	h1_vEff_fadc_l -> SetLineColor(2);
	TH1D * h1_vEff_tdc     = new TH1D("h1_vEff_tdc"    ,"all bars;effective speed of light [cm/ns]"        ,60,10,17);	h1_vEff_tdc    -> SetLineColor(1);
        TH1D * h1_vEff_fadc    = new TH1D("h1_vEff_fadc"   ,"all bars;effective speed of light [cm/ns]"        ,60,10,17);	h1_vEff_fadc   -> SetLineColor(2);

	// ----------------------------------------------------------------------------------
	// Opening input HIPO file
	hipo::reader reader;
	reader.open(inputFile);

	BBand band_hits   ("BAND::hits"       ,reader);

	int event_counter = 0;
	// ----------------------------------------------------------------------------------
	// Loop over events and print them on the screen
	while(reader.next()==true){

		if(event_counter%1000000==0) cout << "event: " << event_counter << endl;
		event_counter++;

		int nHits = band_hits.getSize();
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

		
			//band_hits.show();

			if( adcLcorr <  6000 || adcRcorr <  6000 ) continue;
			if( adcLcorr > 20000 || adcRcorr > 20000 ) continue;

			h1_tdc_diff [barKey]->Fill( difftimeTdc );
			h1_ftdc_diff[barKey]->Fill( difftimeFadc );
			h2_empty    [barKey]->Fill( difftimeTdc,difftimeFadc);
			h2_meantime [barKey]->Fill( TMath::Sqrt(adcLcorr*adcRcorr) ,(meantimeTdc-meantimeFadc));

		}// end loop over hits in event
	}// end file

	// Get the width for each histogram:
	ofstream lr_offsets,effective_velocity;
	lr_offsets.open("output/TEST_lr_offsets.txt");
	effective_velocity.open("output/TEST_effective_velocity.txt");
	double ftdc_LEFTS[600] = {0.};
	double ftdc_RIGHTS[600] = {0.};
	double ftdc_HEIGHTS[600] = {0.};
	double tdc_LEFTS[600] = {0.};
	double tdc_RIGHTS[600] = {0.};
	double tdc_HEIGHTS[600] = {0.};
	for( int sector = 1 ; sector < 6 ; sector++){
		for( int layer = 1 ; layer < 7 ; layer++ ){
			for( int component = 1 ; component <= slc[layer-1][sector-1] ; component++ ){
				int i = sector*100 + layer*10 + component;

				double tdc_vEff,tdc_lr_off,ftdc_vEff,ftdc_lr_off;
				if( h1_tdc_diff[i]->Integral() == 0 ){
					tdc_vEff = 0.;
					tdc_lr_off = 0.;
					ftdc_vEff = 0.;
					ftdc_lr_off = 0.;
				}
				else{
					double barLength = bandlen[sector-1];
			
					// Taking just 50% is problematic for some bars as the L-R distribution is asymmetric. What I
					// will do now is first find the rough middle of the distribution given by the mean and then 
					// taking an average in the middle, and finding first bin about 50% for that average. This will
					// also help with FADC
					
					// For FADC
					int meanBin_ftdc = h1_ftdc_diff[i]->FindBin( h1_ftdc_diff[i]->GetMean() );
					double avg_ftdc = h1_ftdc_diff[i]->Integral( meanBin_ftdc-5, meanBin_ftdc+5, "width");
	
					// Now get the bin that correponds to 20% of max for this
					int ftdc_low_bin 	= h1_ftdc_diff[i]->FindFirstBinAbove( avg_ftdc * 0.2 );
					int ftdc_high_bin	= h1_ftdc_diff[i]->FindLastBinAbove ( avg_ftdc * 0.2 );
					double ftdc_left	= h1_ftdc_diff[i]->GetXaxis()->GetBinCenter( ftdc_low_bin  );
					double ftdc_right	= h1_ftdc_diff[i]->GetXaxis()->GetBinCenter( ftdc_high_bin );
					double ftdc_width	= ftdc_right - ftdc_left;
					ftdc_LEFTS[i] 	= ftdc_left;
					ftdc_RIGHTS[i]	= ftdc_right;
					ftdc_HEIGHTS[i] = avg_ftdc * 0.2;
					ftdc_vEff 		= (2*barLength) / ftdc_width;
					ftdc_lr_off 		= (ftdc_left+ftdc_right)/2.;

					// For TDC
					int meanBin_tdc = h1_tdc_diff[i]->FindBin( h1_tdc_diff[i]->GetMean() );
					double avg_tdc = h1_tdc_diff[i]->Integral( meanBin_tdc-5 , meanBin_tdc+5, "width");
					
					int tdc_low_bin		= h1_tdc_diff[i]->FindFirstBinAbove( avg_tdc * 0.2 );
					int tdc_high_bin	= h1_tdc_diff[i]->FindLastBinAbove ( avg_tdc * 0.2 );
					double tdc_left		= h1_tdc_diff[i]->GetXaxis()->GetBinCenter( tdc_low_bin  );
					double tdc_right	= h1_tdc_diff[i]->GetXaxis()->GetBinCenter( tdc_high_bin );
					double tdc_width = tdc_right - tdc_left;
					tdc_LEFTS[i]	= tdc_left;
					tdc_RIGHTS[i]	= tdc_right;
					tdc_HEIGHTS[i]  = avg_tdc * 0.2;
					tdc_vEff = (2*barLength) / tdc_width;
					tdc_lr_off = (tdc_left+tdc_right)/2.;

					if( isinf(tdc_vEff) ) tdc_vEff = 0.;
					if( isinf(ftdc_vEff) ) ftdc_vEff = 0.;
				}

				effective_velocity << sector << "\t" << layer << "\t" << component << "\t" << tdc_vEff << "\t" << ftdc_vEff << "\t" << 0.000 << "\t" << 0.000 << "\n";
				lr_offsets << sector << "\t" << layer << "\t" << component << "\t" << tdc_lr_off << "\t" << ftdc_lr_off << "\t" << 0.000 << "\t" << 0.000 << "\n";
			
				h1_vEff_tdc  -> Fill(tdc_vEff );
				h1_vEff_fadc -> Fill(ftdc_vEff);

				if(i>=300&&i<500){ // Short bars only
					h1_vEff_tdc_s  -> Fill(tdc_vEff );
                                	h1_vEff_fadc_s -> Fill(ftdc_vEff);
				}
				else if(i<200){ // Medium bars only
					h1_vEff_tdc_m  -> Fill(tdc_vEff );
                                        h1_vEff_fadc_m -> Fill(ftdc_vEff);
				}
				else{ // Long bars only
					h1_vEff_tdc_l  -> Fill(tdc_vEff );
                                        h1_vEff_fadc_l -> Fill(ftdc_vEff);
				}
			}
		}
	}
	lr_offsets.close();
	effective_velocity.close();

	// --------------------------------------------------------------------------------------------------------
	// Fitting lines to extract the slope and try to determine the needed TDC -> time conversion factor
        TProfile ** p_empty = new TProfile*[nHistos];
        // Fitting paddle-to-paddle parameters
        for(int i = 0 ; i < nHistos ; i++){
                int notEmpty = (h2_empty[i]->Integral());
                if(notEmpty){
                        int nBins = h2_empty[i] -> GetXaxis() -> GetNbins();
                        p_empty[i] = h2_empty[i]->ProfileX(Form("px_%i",i), 1, nBins );
                        TF1 * f_line = new TF1("f_line","pol1");
                        p_empty[i] -> Fit("f_line","Q");
                        par_pad[i][0] = f_line->GetParameter(1);
                        par_pad[i][1] = f_line->GetParError (1);
                
			val_slope += par_pad[i][0];
        		ctr_slope++;

			h1_slopes -> Fill(par_pad[i][0]);

			if( i >= 300 && i < 500 )
				h1_slopes_short -> Fill(par_pad[i][0]);
			else if( i < 200 )
				h1_slopes_med -> Fill(par_pad[i][0]);
			else
				h1_slopes_long -> Fill(par_pad[i][0]);
		}
        }

	double avg_m = val_slope/((double)(ctr_slope));
	double std_TDC2ns  = 0.02345;
	double std_FADC2ns = 0.0625;
	cout << "Average slope = " << avg_m << endl;
	cout << "TDC  to ns conversion factor (assuming TDC  t conversion factor is " << std_TDC2ns << "): " << avg_m*std_TDC2ns  << endl;
	cout << "FADC to ns conversion factor (assuming FADC t conversion factor is " << std_FADC2ns  << "): " << std_FADC2ns/avg_m << endl;

	TCanvas * c1 = new TCanvas("c1","c1");
	h1_slopes_long  -> Draw();
	h1_slopes_short -> Draw("same");
	h1_slopes_med   -> Draw("same");
	//h1_slopes       -> Draw("same");
        TLegend * leg = new TLegend (0.6,0.5,0.8,0.7);
	leg -> SetLineColor(0);
	leg -> AddEntry(h1_slopes_short,"short" );
	leg -> AddEntry(h1_slopes_med  ,"medium");
	leg -> AddEntry(h1_slopes_long ,"long"  );
	leg -> Draw("same");
	c1 -> Modified();
        c1 -> Update();

	// --------------------------------------------------------------------------------------------------------
	// Histograms with all the extracted speeds of light
	TCanvas * c2 = new TCanvas("c2","c2");
	c2 -> Divide(2,2);
	c2 -> cd(1);	h1_vEff_tdc_s -> Draw();	h1_vEff_fadc_s -> Draw("same");
	c2 -> cd(2);	h1_vEff_tdc_m -> Draw();	h1_vEff_fadc_m -> Draw("same");
	c2 -> cd(3);	h1_vEff_tdc_l -> Draw();	h1_vEff_fadc_l -> Draw("same");
	c2 -> cd(4);	h1_vEff_tdc   -> Draw();	h1_vEff_fadc   -> Draw("same");
	
	c2 -> cd(2);
	TLegend * leg1 = new TLegend(0.15,0.6,0.30,0.87);
	leg1 -> SetLineColor(0);
	leg1 -> AddEntry( h1_vEff_tdc  , "TDC"  );
	leg1 -> AddEntry( h1_vEff_fadc , "FADC" );
	leg1 -> Draw("same");
	c2 -> Modified();
	c2 -> Update();

	// --------------------------------------------------------------------------------------------------------
	// Create plots for all the L-R distributions of TDC and FADC
	TCanvas *** c_tdc_diff  = new TCanvas**[5];
	TCanvas *** c_empty     = new TCanvas**[5];
	TCanvas *** c_meantime  = new TCanvas**[5];

	int min, max;

	for(int is = 0 ; is < 5 ; is++){
		c_tdc_diff [is] = new TCanvas*[6];
		c_empty    [is] = new TCanvas*[6];
		c_meantime [is] = new TCanvas*[6];

		for(int il = 0 ; il < 5 ; il++){
			c_tdc_diff [is][il] = new TCanvas(Form("c_tdc_diff_S%iL%i" ,is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);	c_tdc_diff [is][il] -> Divide(2,7);
			c_empty    [is][il] = new TCanvas(Form("c_empty_S%iL%i"    ,is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);	c_empty    [is][il] -> Divide(2,4);
			c_meantime [is][il] = new TCanvas(Form("c_meantime_S%iL%i" ,is,il),Form("Sector %i, Layer %i",is+1,il+1),900,900);	c_meantime [is][il] -> Divide(2,4);

			for(int cIdx = 0 ; cIdx < slc[il][is] ; cIdx++){
				int identifier = 100*(is+1)+10*(il+1)+(cIdx+1);
				int notEmpty = (h1_tdc_diff[identifier]->Integral()+h1_ftdc_diff[identifier]->Integral());
				if(notEmpty){
					// --------------------------------------------------------------
					c_tdc_diff[is][il] -> cd(2*cIdx+1);
					gPad -> SetBottomMargin(0.26);

					// Draw TDC
					PrettyTH1F( h1_tdc_diff[identifier] , "TDC L-R [ns]","", 2);
					h1_tdc_diff[identifier]->Draw();

					// Draw two lines:
					TLine * lLow 	= new TLine( tdc_LEFTS[identifier],  0, tdc_LEFTS[identifier],  tdc_HEIGHTS[identifier] );
					TLine * lHigh 	= new TLine( tdc_RIGHTS[identifier], 0, tdc_RIGHTS[identifier], tdc_HEIGHTS[identifier] );
					TLine * across	= new TLine( tdc_LEFTS[identifier],  tdc_HEIGHTS[identifier], tdc_RIGHTS[identifier], tdc_HEIGHTS[identifier]);
					lLow->SetLineColor(4);
					lHigh->SetLineColor(4);
					across->SetLineColor(4);
					lLow  -> Draw("same");
					lHigh -> Draw("same");
					across->Draw("Same");

					// --------------------------------------------------------------
					c_tdc_diff[is][il] -> cd(2*cIdx+2);
					gPad -> SetBottomMargin(0.26);

					// Draw FADC
					PrettyTH1F( h1_ftdc_diff[identifier] , "FADC L-R [ns]","", 1);
					h1_ftdc_diff[identifier]->Draw();
					// make two lines
					TLine * lLow2 	= new TLine( ftdc_LEFTS[identifier],  0, ftdc_LEFTS[identifier],  ftdc_HEIGHTS[identifier] );
					TLine * lHigh2 	= new TLine( ftdc_RIGHTS[identifier], 0, ftdc_RIGHTS[identifier], ftdc_HEIGHTS[identifier] );
					TLine * across2 = new TLine( ftdc_LEFTS[identifier],  ftdc_HEIGHTS[identifier], ftdc_RIGHTS[identifier], ftdc_HEIGHTS[identifier] );
					lLow2->SetLineColor(8);
					lHigh2->SetLineColor(8);
					across2->SetLineColor(8);
					lLow2  -> Draw("same");
					lHigh2 -> Draw("same");
					across2->Draw("same");

					// --------------------------------------------------------------
					c_empty[is][il] -> cd(cIdx+1);
					gPad -> SetTopMargin   (0.00);
                                        gPad -> SetBottomMargin(0.26);
                                        gPad -> SetLeftMargin  (0.26);	
					double fadc_len = ftdc_RIGHTS[identifier] - ftdc_LEFTS[identifier];
					double tdc_len = tdc_RIGHTS[identifier] - tdc_LEFTS[identifier];
					h2_empty[identifier]->SetTitle(Form("FADC Length: %f / TDC Length: %f",fadc_len,tdc_len));
					h2_empty[identifier]->SetStats(0);
					h2_empty[identifier]->Draw("col");
					p_empty [identifier]->Draw("same");
					TLine * slope1 = new TLine(-20,-20,20,20);
					slope1->SetLineColor(1);
					slope1->Draw("same");

					// --------------------------------------------------------------
					c_meantime[is][il] -> cd(cIdx+1);
					gPad -> SetTopMargin   (0.00);
					gPad -> SetBottomMargin(0.26);
					gPad -> SetLeftMargin  (0.26);
					h2_meantime[identifier]->SetStats(0);
					min = getMinNonEmptyBin(h2_meantime[identifier]);
					max = getMaxNonEmptyBin(h2_meantime[identifier]);
					h2_meantime[identifier] -> GetYaxis() -> SetRange(min,max);
					h2_meantime[identifier] -> Draw("col");

				}
			}		
			c_tdc_diff [is][il] -> Modified();	c_tdc_diff [is][il] -> Update();	
			c_empty    [is][il] -> Modified();	c_empty    [is][il] -> Update();
			c_meantime [is][il] -> Modified();	c_meantime [is][il] -> Update();
		}
	}

	// Saving to pdf
	TCanvas * c0 = new TCanvas("c0","",900,900);
	c0 -> Print("output/results_offsets_veff.pdf(");
	c1 -> Print("output/results_offsets_veff.pdf");
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			c_tdc_diff[is][il] -> Print("output/results_offsets_veff.pdf");
		}
	}
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			c_empty[is][il] -> Print("output/results_offsets_veff.pdf");
		}
	}
	for(int is = 0 ; is < 5 ; is++){
		for(int il = 0 ; il < 5 ; il++){
			c_meantime[is][il] -> Print("output/results_offsets_veff.pdf");
		}
	}
	c0 -> Print("output/results_offsets_veff.pdf)");

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
	//h2 -> GetXaxis() -> SetTitle(titx);
	//h2 -> GetYaxis() -> SetTitle(tity);
	h2 -> GetXaxis() -> SetTitleSize(0.10);
        h2 -> GetYaxis() -> SetTitleSize(0.10);
	h2 -> GetXaxis() -> SetLabelSize(0.10);
        h2 -> GetYaxis() -> SetLabelSize(0.10);
	h2 -> GetXaxis() -> CenterTitle();
        h2 -> GetYaxis() -> CenterTitle();
	h2 -> GetXaxis() -> SetNdivisions(105);
        h2 -> GetYaxis() -> SetNdivisions(107);
}
// ========================================================================================================================================
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
// ========================================================================================================================================
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
// ========================================================================================================================================
