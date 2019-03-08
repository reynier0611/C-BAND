#include "PParticle.h"
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TFile.h"
// ===================================================================================================================================
PParticle::PParticle()
{
}
// ===================================================================================================================================
PParticle::PParticle(TVector3 mom0,double z_m0)
{
	mp = 0.938272; // GeV -> Proton Mass
        mn = 0.939566; // GeV -> Neutron Mass

	mom = mom0;
	z_m = z_m0;

}
// ===================================================================================================================================
PParticle::~PParticle()
{
}
// ===================================================================================================================================
bool PParticle::pointsToBand(){
        
        double theta = mom.Theta();
        double phi   = mom.Phi();
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
