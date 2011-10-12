//fraction of events with SUSY_genbHT>350 in the mGL, mLSP plane
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TChain.h>
#include <TCut.h>
#include <TLegend.h>
#include <TFile.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TGraphPainter.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
using namespace std;



void plotFractionEvents_SUSY_genbHT_ge350()
{

  TChain ch("reducedTree");
  ch.Add("/cu2/ra2b/reducedTrees/specialNoHtCut/reducedTree.SSVHPT_PF2PATjets_JES0_JER0_PFMET_METunc0_PUunc0_BTagEff0_HLTEff0.SMS-T1bbbb_Mgluino-100to1200_mLSP-50to1150_7TeV-Pythia6Z.root");

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);

  TFile * fout = new TFile("htCurveT1bbbb.root","RECREATE");

  TCanvas c1("c1","c1");

  int nm = 60; double xmlow = 0; double xmhigh = 1500;
  double mbinwidth = 25; //scan points are generated 25 GeV apart
  int nm_allowed = 60; //no scan points for masses >= 1225

  int nHT = 100; double xHTlow = 0; double xHThigh = 1000;
  double HT_binwidth = double ( xHThigh - xHTlow ) / (double) nHT;

  TH2D h2d("genbHTfrac","", nm, xmlow, xmhigh, nm, xmlow, xmhigh);


  double x20[ 100 ] = { }; 
  double y20[ 100 ] = { }; 
  double f20[ 100 ] = { };


  for ( int ix=1; ix<=nm_allowed; ix++ ) {

    double xxlow = xmlow + ( ix - 1) * mbinwidth;

    if ( xxlow < 100 ) continue; //no scan points below this mass
    if ( xxlow > 1200 ) continue; //no scan points above this mass

    double tempDelta = 10.; double tempM = xmlow; double tempF = 99.;

    for ( int iy=1; iy<=nm_allowed; iy++ ) {
      
      //one value of (ix,iy) corresponds to one point in the m0, m12 plane
      //e.g. one scan point
      
      double xylow = xmlow + ( iy - 1) * mbinwidth;
      double xhigh = xmlow + ix * mbinwidth;


      if ( xylow < 50 ) continue; //no scan points below this mass
      if ( xylow > 1200 ) continue; //no scan points above this mass

      if ( xxlow <= xylow ) continue; //upper plane forbidden

      if ( xxlow - xylow == 25 ) continue; //no scan point here


      stringstream masscut_ss; masscut_ss<<"m0=="<<xxlow<<" && m12=="<<xylow;
      string masscut_str = masscut_ss.str();


      TH1D h_HT("h_HT","",nHT, xHTlow,xHThigh);
      ch.Draw("SUSY_genbHT>>h_HT",masscut_str.c_str());

      double numEvents_scanPoint = h_HT.GetEntries();
      if ( numEvents_scanPoint == 0 ) {
        cout<<"there are no events at the scan point with m0="<<xxlow<<" and m12="<<xylow<<" !"<<endl;
        continue;
      }

      else{

	int ibin_HT350 = 1 + (float) 350 / (float) HT_binwidth;

//      //check that it's 350!
//      int HT_min = h_HT.GetBinLowEdge( ibin_HT350 );
//      if ( HT_min != 350 ) {
//	cout<<"something went wrong"<<endl;
//	continue;
//      }

      
	double numEvents_scanPoint_SUSYgenbHTge350 = h_HT.Integral( ibin_HT350 , nHT + 1 ); 

	double fractionEvents_scanPoint_SUSYgenbHTge350 = numEvents_scanPoint_SUSYgenbHTge350 / numEvents_scanPoint;
	cout<<fractionEvents_scanPoint_SUSYgenbHTge350<<endl;

	if ( fractionEvents_scanPoint_SUSYgenbHTge350 > 0 && fractionEvents_scanPoint_SUSYgenbHTge350 <= 1 ) {
	  h2d.SetBinContent( ix, iy, fractionEvents_scanPoint_SUSYgenbHTge350 );
	}
	else if ( fractionEvents_scanPoint_SUSYgenbHTge350 == 0 ) {
	  h2d.SetBinContent( ix, iy, 0.0001 );
	}
	else cout<<"fractionEvents_scanPoint_SUSYgenbHTge350 has an unphysiscal value: "<<fractionEvents_scanPoint_SUSYgenbHTge350<<endl;


	double thisDelta = abs ( fractionEvents_scanPoint_SUSYgenbHTge350 - 0.2 ) ;
	if ( thisDelta < tempDelta ) {
	  tempDelta = thisDelta;
	  tempM = xylow;
	  tempF = fractionEvents_scanPoint_SUSYgenbHTge350;
	}

      }

    }//end loop over y (m12) bins


    //determine line coordinates
    if ( tempF > 0.15 && tempF < 0.25 ) {
      x20[ ix - 1 ] = xxlow;
      y20[ ix - 1 ] = tempM;
      f20[ ix - 1 ] = tempF;
//      cout<<"for m0="<<xxlow<<" the bin with content closest to 0.2 is m12="<<tempM<<endl;
//      cout<<"this bin's content is: "<<tempF<<endl;
    }

  }//end loop over xy (m0) bins


  for (int im=0; im<nm; im++) cout<<x20[im]<<" "<<y20[im]<<endl;

  //get rid of 0's
  double x20_f[ 100 ];
  double y20_f[ 100 ];
  double f20_f[ 100 ];
  int i0 = 0;

  for (int i=0; i<100; i++) {
    if ( x20[ i ] == 0 || y20[ i ] == 0 ) continue;
    else{
      cout<<i0<<endl;
      x20_f[ i0 ] = x20[ i ];
      y20_f[ i0 ] = y20[ i ];
      f20_f[ i0 ] = f20[ i ];
      i0++;
    }
  }
  //extend the line to upper right of last bin
  x20_f[ i0 ] =  x20_f[ i0 - 1 ] + 25;
  y20_f[ i0 ] =  y20_f[ i0 - 1 ] + 25;


  TGraph line( i0, x20_f, y20_f );
  line.SetMarkerStyle( 1 );
  line.SetName("htCurveT1bbb");

  h2d.GetXaxis()->SetTitle("m_{0}");
  h2d.GetYaxis()->SetTitle("m_{12}");
  h2d.SetTitle("Fraction of events with SUSY_genbHT>350");

  h2d.Draw("COLZ");
  line.Draw("sameCP0");
  c1.SaveAs("fraction_events_SUSY_genbHT_ge350_20pct.pdf");
  
  line.Write();
  h2d.Write();
  fout->Close();
}



