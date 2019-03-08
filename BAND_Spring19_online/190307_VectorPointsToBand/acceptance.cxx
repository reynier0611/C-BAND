#include <cstdlib>
#include <iostream>

#include "TStyle.h"
#include "TRint.h"
#include "TApplication.h"
#include "TRandom.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLine.h"
#include "TGraph.h"
#include "TBox.h"
#include "TLegend.h"

using namespace std;

bool pointsToBand(double theta,double phi,double z_m);
void prettyTH2D( TH2D * h2 , int color );
// ========================================================================================================================
int main(int argc, char ** argv) {

#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	gStyle -> SetOptStat(0);

	const int nMC = 1000000;

	double thickness  = 7.2;                                // thickness of each bar (cm)

	// Distance from ideal target to upstream end of BAND
	// (from BAND survey report, 02/18/2019)
	double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

	// Distance from ideal target to downstream end of layer 5
	double zDown = zUpst + 5*thickness;

	double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

	double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};
	
	// Sector boundaries
        double topSec1  = globalY + 13*thickness;
        double topSec2  = globalY + 10*thickness;
        double topSec34 = globalY +  3*thickness;
        double topSec5  = globalY -  3*thickness;
        double downSec5 = globalY -  5*thickness;

	// -------------------------------------------------------------------
	// Defining histograms
	TH2D * h2_th_phi_1 = new TH2D("h2_th_phi_1",";#phi [deg];#theta [deg]",100,-190,190,100, 150,180);
	TH2D * h2_th_phi_2 = new TH2D("h2_th_phi_2",";#phi [deg];#theta [deg]",100,-190,190,100, 150,180);
	TH2D * h2_x_y_1    = new TH2D("h2_x_y_1"   ,";x [cm];y [cm]"          ,100,-110,110,100,- 50,150);
	TH2D * h2_x_y_2    = new TH2D("h2_x_y_2"   ,";x [cm];y [cm]"          ,100,-110,110,100,- 50,150);

	prettyTH2D( h2_th_phi_1 ,62 );
	prettyTH2D( h2_th_phi_2 ,92 );
	prettyTH2D( h2_x_y_1    ,62 );
	prettyTH2D( h2_x_y_2    ,92 );

	// -------------------------------------------------------------------
	TLine ** trk = new TLine * [nMC];

	int nLines = 0;

	double theta_min = 150;
	double theta_max = 180;
	double y_min_gen = zDown/TMath::Cos(theta_min)*TMath::Sin(theta_min);
	double y_max_gen = zDown/TMath::Cos(theta_max)*TMath::Sin(theta_max);

	// --------------------------------------------------------------------
	// Doing monte carlo generation in theta and phi
	for( int mc = 0 ; mc < nMC ; mc++ ){

		if(mc%100000 == 0) cout << mc << "\t out of:" << nMC << endl;
	
		double genTh = gRandom -> Uniform(theta_min*TMath::DegToRad(),theta_max*TMath::DegToRad());
		double genPh = gRandom -> Uniform(     -180*TMath::DegToRad(),      180*TMath::DegToRad());

		double rho   = zDown/TMath::Cos(genTh);
		double xDown = rho*TMath::Sin(genTh)*TMath::Cos(genPh);
		double yDown = rho*TMath::Sin(genTh)*TMath::Sin(genPh);	

		h2_th_phi_1 -> Fill(genPh*TMath::RadToDeg(),genTh*TMath::RadToDeg());
		h2_x_y_1    -> Fill(xDown                  ,yDown                  );

		if(pointsToBand(genTh,genPh,0)){
			h2_th_phi_2 -> Fill(genPh*TMath::RadToDeg(),genTh*TMath::RadToDeg());
			h2_x_y_2    -> Fill(xDown                  ,yDown                  );

		
			trk[nLines] = new TLine(0,0,zDown,yDown);		
			if(
				(yDown < topSec34 && yDown >= topSec5  && TMath::Abs(xDown) < bandlen[1]/2. && TMath::Abs(xDown) > bandlen[1]/2.-bandlen[2])
			)
				trk[nLines] -> SetLineColor(92);
			else
				trk[nLines] -> SetLineColor(94);

			nLines++;
		}
	}

	// --------------------------------------------------------------------
	// Plotting results
	
	TCanvas * c1 = new TCanvas("c1","c1",1300,900);
        h2_x_y_1 -> Draw();
        h2_x_y_2 -> Draw("same");

        TLine * lh = new TLine(-5, 0,5,0);      lh -> SetLineWidth(2);  lh -> Draw("same");
        TLine * lv = new TLine( 0,-5,0,5);      lv -> SetLineWidth(2);  lv -> Draw("same");

	TLegend * leg1 = new TLegend(0.70,0.75,0.87,0.87);
	leg1 -> AddEntry( h2_x_y_1 , "generated" );
	leg1 -> AddEntry( h2_x_y_2 , "on BAND"   );
	leg1 -> AddEntry( lh       , "beam pipe" );
	leg1 -> Draw("same");

        c1 -> Modified();
        c1 -> Update();

	// -------------------

	TCanvas * c2 = new TCanvas("c2","c2",1300,900);	
	h2_th_phi_1 -> Draw();	
	h2_th_phi_2 -> Draw("same");

	TLegend * leg2 = new TLegend(0.24,0.17,0.37,0.30);
        leg2 -> AddEntry( h2_x_y_1 , "generated" );
        leg2 -> AddEntry( h2_x_y_2 , "on BAND"   );	
	leg2 -> Draw("same");

	c2 -> Modified();
	c2 -> Update();

	// -------------------
	// Visualization of tracks
	double z[] = {-300,0  };
	double y[] = {-180,180};
	TGraph * gDummy = new TGraph(2,z,y);
	gDummy -> SetMinimum(-50);
	gDummy -> SetMaximum(120);
	gDummy -> GetXaxis() -> SetNdivisions(409);
        gDummy -> GetYaxis() -> SetNdivisions(309);
	gDummy -> GetXaxis() -> SetTitle("z [cm]");
	gDummy -> GetYaxis() -> SetTitle("y [cm]");
	gDummy -> GetXaxis() -> CenterTitle();
        gDummy -> GetYaxis() -> CenterTitle();
	gDummy -> SetTitle("");	

	gStyle->SetFillColorAlpha(62,0.4);
	TBox * TopBox = new TBox(zUpst,globalY+3*thickness,zDown,globalY+13*thickness);	TopBox -> SetLineColor(62);	TopBox -> SetLineWidth(3);
	TBox * BotBox = new TBox(zUpst,globalY-5*thickness,zDown,globalY- 3*thickness);	BotBox -> SetLineColor(62);	BotBox -> SetLineWidth(3);
	TBox * TotBox = new TBox(zUpst,globalY-5*thickness,zDown,globalY+13*thickness);	TotBox -> SetLineColor(62);	TotBox -> SetLineWidth(3);

	TCanvas * c3 = new TCanvas("c3","c3",1300,900);
	gDummy -> Draw("AP");
	for(int i = 0 ; i < nLines ; i++) trk[i] -> Draw("same");
	
	TopBox -> Draw("samel");
	TotBox -> Draw("samel");
	BotBox -> Draw("samel");

	TGraph * gTrkL = new TGraph();	gTrkL -> SetLineColor(94);	gTrkL -> SetLineWidth(4);	gTrkL -> SetFillColor(0);
	TGraph * gTrkS = new TGraph();	gTrkS -> SetLineColor(92);	gTrkS -> SetLineWidth(4);	gTrkS -> SetFillColor(0);

	TLegend * leg3 = new TLegend(0.6,0.6,0.87,0.87);
	leg3 -> SetLineColor(0);
	leg3 -> AddEntry( TopBox , "BAND position"        );
	leg3 -> AddEntry( gTrkL  , "tracks on long bars"  );
	leg3 -> AddEntry( gTrkS  , "tracks on short bars" );
	leg3 -> Draw("same");

	c3 -> Modified();
        c3 -> Update();

	// --------------------------------------------------------------------
	// Saving plots
	c1 -> Print("c1.png");
	c2 -> Print("c2.png");
	c3 -> Print("c3.png");

	myapp -> Run();
	return 0;

}
// ========================================================================================================================
bool pointsToBand(double theta,double phi,double z_m){
	double z = z_m*100; // from m to cm

	// Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
	double thickness  = 7.2;                                // thickness of each bar (cm)
	double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6

	// Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java
	double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

	// Distance from ideal target to upstream end of BAND
	// (from BAND survey report, 02/18/2019)
	double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

	// Distance from ideal target to downstream end of layer 5
	double zDown = (zUpst + 5*thickness) - z_m;

	double rho   = zDown/TMath::Cos(theta);
	double xDown = rho*TMath::Sin(theta)*TMath::Cos(phi);
	double yDown = rho*TMath::Sin(theta)*TMath::Sin(phi);

	double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
	double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

	// Sector boundaries
	double topSec1  = globalY + 13*thickness;
	double topSec2  = globalY + 10*thickness;
	double topSec34 = globalY +  3*thickness;
	double topSec5  = globalY -  3*thickness;
	double downSec5 = globalY -  5*thickness;

	if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

	if(             (yDown < topSec1  && yDown >= topSec2  && TMath::Abs(xDown) < bandlen[0]/2. )||
			(yDown < topSec2  && yDown >= topSec34 && TMath::Abs(xDown) < bandlen[1]/2. )||
			(yDown < topSec34 && yDown >= topSec5  && TMath::Abs(xDown) < bandlen[1]/2. && TMath::Abs(xDown) > bandlen[1]/2.-bandlen[2])||
			(yDown < topSec5  && yDown >= downSec5 && TMath::Abs(xDown) < bandlen[4]/2. )
	  )
		return 1;

	return 0;
}
// ========================================================================================================================
void prettyTH2D( TH2D * h2 , int color ){
	h2 -> SetMarkerStyle(7);
	h2 -> SetMarkerColor(color);
	h2 -> GetXaxis() -> CenterTitle();
	h2 -> GetYaxis() -> CenterTitle();
	h2 -> GetXaxis() -> SetNdivisions(409);
	h2 -> GetYaxis() -> SetNdivisions(409);
	h2 -> SetLineColor(color);
	h2 -> SetLineWidth(3);
}
