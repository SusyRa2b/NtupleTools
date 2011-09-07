#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"

#include <vector>
#include <set>
#include <iostream>

TCanvas * c1;

Long_t counter_=0;

TFile * outputfile=0;

std::set<pair<double,double> > smoothCurve(const std::set<pair<double,double> > & inputPoints, TSpline3 ** sp,
					   pair<double,double> * leftpoint=0, pair<double,double> * rightpoint=0) {
  Double_t * xx = new Double_t[ inputPoints.size() ];
  Double_t * yy = new Double_t[ inputPoints.size() ];
  std::set<pair<double,double> > xyvals_smoothed;
  int j=0;
  for (std::set<pair<double,double> >::iterator i=inputPoints.begin(); i!=inputPoints.end(); ++i) {
    xx[j] = i->first;
    yy[j] = i->second;
    if (j%2 !=0 ) xyvals_smoothed.insert(make_pair( 0.5*(xx[j]+xx[j-1]), 0.5*(yy[j]+yy[j-1]) ) );
    j++;
  }

  int mysize = inputPoints.size();
  if (leftpoint!=0) ++mysize;
  if (rightpoint!=0) ++mysize;
  Double_t * xxx = new Double_t[ mysize ];
  Double_t * yyy = new Double_t[ mysize ];
  j=0;
  if (leftpoint != 0) {
    xxx[0] = leftpoint->first;
    yyy[0] = leftpoint->second;
    j++;
  }
  for (std::set<pair<double,double> >::iterator i=inputPoints.begin(); i!=inputPoints.end(); ++i) {
    xxx[j] = i->first;
    yyy[j] = i->second;
    j++;
  }
  if (rightpoint != 0) {
    xxx[j] = rightpoint->first;
    yyy[j] = rightpoint->second;
  }

  *sp = new TSpline3(TString("spline")+counter_,xxx,yyy,mysize);
  ++counter_;

  delete [] xx;
  delete [] yy;
  return xyvals_smoothed;
}


void makeSplines(TString file) {

  TString id=file(file.Index("ge"),9);
  if (c1!=0) delete c1;
  c1=new TCanvas(id);

  TFile* fmsugra=new TFile(file);

  TH2F* hcls = (TH2F*) fmsugra->Get("hcls");

  TH2* hist = (TH2*) hcls->Clone("interpolated");
  int xMax = hist->GetNbinsX();
  int yMax = hist->GetNbinsY();

  int xMin = 0;
  
  {// for interactive root
    for(int xBin = 1; xBin <= xMax; xBin++)
      {
	for(int yBin = 1; yBin <= yMax; yBin++)
	  {
	    if(hist->GetBinContent(xBin,yBin)>0)
	      {
		xMin = xBin;
		yBin = yMax+1;
		xBin = xMax+1;
	      }
	  }
      }
  }
  
  hist->Draw("COLZ");

  { //interactive root
    for (int jj=0; jj<25; jj++) {
      for(int xBin = xMin; xBin <= xMax; xBin++)  {
	for(int yBin = 1; yBin <= yMax; yBin++) {
	  if (hist->GetBinContent(xBin,yBin) <=0 ) {
	    double av=0;
	    int n=0;
	    if (hist->GetBinContent(xBin + 1,yBin + 1)>0 && xBin!=xMax && yBin!=yMax) { av+= hist->GetBinContent(xBin + 1,yBin + 1); ++n;}
	    if (hist->GetBinContent(xBin + 1,yBin)    >0 && xBin!=xMax) { av+= hist->GetBinContent(xBin + 1,yBin);     ++n;}
	    if (hist->GetBinContent(xBin ,yBin + 1)   >0 && yBin!=yMax) { av+= hist->GetBinContent(xBin ,yBin + 1);    ++n;}
	    
	    if (n>0)  hist->SetBinContent(xBin,yBin,av/double(n));
	  }
	}
      }
      hist->Draw("COLZ");
    }
  }

  const double cutval=0.05;

  std::set<pair<double,double> > xyvals;
  //  std::vector<double> yvals;
  for (int xBin=xMin; xBin<=xMax; xBin++) {
    bool foundAtThisX=false;
    for (int yBin=yMax-1; yBin>=1; yBin--) {
      
      if (foundAtThisX) continue;

      if (hist->GetBinContent(xBin, yBin) < 0 || hist->GetBinContent(xBin, yBin+1)<0) continue;

      if (hist->GetBinContent(xBin, yBin) <cutval && hist->GetBinContent(xBin, yBin+1) >=cutval) {
	//found a transition
	//need to keep track of the x and y values of the transition
	xyvals.insert(make_pair( hist->GetXaxis()->GetBinCenter(xBin),hist->GetYaxis()->GetBinLowEdge(yBin+1) ));
	//	xvals.push_back( hist->GetXaxis()->GetBinCenter(xBin) );
	//	yvals.push_back( hist->GetYaxis()->GetBinLowEdge(yBin+1) );

	foundAtThisX=true;
      }
      else if (false && hist->GetBinContent(xBin, yBin) >=cutval && hist->GetBinContent(xBin, yBin+1) <cutval) {
	cout<<"reverse transition found..."<<endl;
	//	xvals.push_back( hist->GetXaxis()->GetBinCenter(xBin) );
	//	yvals.push_back( hist->GetYaxis()->GetBinLowEdge(yBin+1) );
      }
      
    }
  }

  pair<double,double> leftpoint = *xyvals.begin();
  pair<double,double> rightpoint = *xyvals.rbegin();

  TSpline3* spl=0;
  TSpline3* spl2=0;
  TSpline3* spl3=0;
  TSpline3* spl4=0;
  std::set<pair<double,double> > xyvals_smoothed = smoothCurve(xyvals,&spl);
  std::set<pair<double,double> > xyvals_smoothed2 = smoothCurve(xyvals_smoothed,&spl2,&leftpoint,&rightpoint);
  std::set<pair<double,double> > xyvals_smoothed3 = smoothCurve(xyvals_smoothed2,&spl3,&leftpoint,&rightpoint);
  std::set<pair<double,double> > xyvals_smoothed4 = smoothCurve(xyvals_smoothed3,&spl4,&leftpoint,&rightpoint);

  hist->SetMaximum(0.1);
  hist->SetMinimum(0);

  hist->Draw("COLZ");
  //  gr->Draw("same *");

  //  TSpline3 * spl = new TSpline3("spline",xx,yy,xyvals.size());
  //  TSpline3 * spl_smoothed = new TSpline3("spline_smoothed",xx_smoothed,yy_smoothed,xyvals_smoothed.size());

  spl->Draw("SAMEc");
  spl2->Draw("SAMEc");
  spl3->Draw("SAMEc");
  spl4->Draw("SAMEc");

  outputfile->cd();
  spl3->SetName(TString("curve3_")+id);
  spl4->SetName(TString("curve4_")+id);
  spl3->Write();
  spl4->Write();

}

void check() {

  TString which="ge2bloose";

  TFile * splinefile = new TFile("RA2b_tb40_exclusion.root");
  TFile * sourcefile = new TFile("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-"+which+".root");

  TH2F* hexcluded = (TH2F*) sourcefile->Get("hexcluded");
  hexcluded->Draw("COL");
  TSpline3* sp=  (TSpline3*)  splinefile->Get("curve4_"+which);

  sp->Draw("SAMEc");

}

void msugraExclusion() {

  outputfile = new TFile("RA2b_tb40_exclusion.root","RECREATE");
  makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge1btight.root");
  makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge1bloose.root");
  makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge2btight.root");
  makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge2bloose.root");

  outputfile->Close();

}
