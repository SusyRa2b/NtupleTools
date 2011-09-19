#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TMath.h"

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


void makeSplines(TString file, TString histname,bool binarymode=false) {

  TString id="";
  if (histname=="") {
    id=file(file.Index("ge"),9);
  }
  else if (histname.Contains("value_") ) {
    id = histname(histname.Index("value_")+6,7);
    id.Prepend("ge");
  }
  else if (histname.BeginsWith("h2d_") ) {
    id =  histname(4,histname.Length());
  }
  else {
    cout<<"histogram name structure not recognized"<<endl; 
    return;
  }
  if (c1!=0) delete c1;
  c1=new TCanvas(id);

  TFile* fmsugra=new TFile(file);

  TH2F* hcls = (TH2F*) fmsugra->Get(histname);

  TH2* hist = (TH2*) hcls->Clone("interpolated");
  int xMax = hist->GetNbinsX();
  int yMax = hist->GetNbinsY();

  int xMin = 0;
  
  {// for interactive root
    for(int xBin = 1; xBin <= xMax; xBin++)
      {
	for(int yBin = 1; yBin <= yMax; yBin++)
	  {
	    if(hist->GetBinContent(xBin,yBin)>-2)
	      {
		xMin = xBin;
		yBin = yMax+1;
		xBin = xMax+1;
		//set bins with -1 to a dummy excluded value of 0.01
		//maybe this is the wrong approach
		//		if ( hist->GetBinContent(xBin,yBin) == -1 )  hist->SetBinContent(xBin,yBin,0.01);
	      }
	  }
      }
  }
  
  hist->Draw("COLZ");

  //for binary mode, we want to interpolate over all values less than 0
  //for non-binary mode, we want to *keep* the -1 values
  const float validcut= binarymode ? -0.5 : -2;

  { //interactive root
    for (int jj=0; jj<25; jj++) {
      for(int xBin = xMin; xBin <= xMax; xBin++)  {
	for(int yBin = 1; yBin <= yMax; yBin++) {
	  if (hist->GetBinContent(xBin,yBin) <= validcut ) { //interpolate only on points with value less than the cut
	    double av=0;
	    int n=0;
	    bool neg1flag=false;
	    bool fillWithOne=false; //for binary mode
	    if (binarymode) { //back to luke's original algorithm
	      if (hist->GetBinContent(xBin + 1,yBin + 1)>0 && xBin!=xMax && yBin!=yMax) fillWithOne=true;
	      if (hist->GetBinContent(xBin + 1,yBin)    >0 && xBin!=xMax) fillWithOne=true;
	      if (hist->GetBinContent(xBin ,yBin + 1)   >0 && yBin!=yMax) fillWithOne=true;
	    }
	    else {
	      if (hist->GetBinContent(xBin + 1,yBin + 1)>0 && xBin!=xMax && yBin!=yMax) { av+= hist->GetBinContent(xBin + 1,yBin + 1); ++n;}
	      if (hist->GetBinContent(xBin + 1,yBin)    >0 && xBin!=xMax) { av+= hist->GetBinContent(xBin + 1,yBin);     ++n;}
	      if (hist->GetBinContent(xBin ,yBin + 1)   >0 && yBin!=yMax) { av+= hist->GetBinContent(xBin ,yBin + 1);    ++n;}
	      
	      
	      if ( TMath::Nint(hist->GetBinContent(xBin + 1,yBin + 1)) == -1 
		   || TMath::Nint(hist->GetBinContent(xBin + 1,yBin))==-1 || TMath::Nint(hist->GetBinContent(xBin ,yBin + 1))==-1)
		neg1flag=true;
	    }

	    if (n>0 && !neg1flag)  hist->SetBinContent(xBin,yBin,av/double(n));
	    else if (neg1flag==true) hist->SetBinContent(xBin,yBin,-1);
	    else if (fillWithOne==true) hist->SetBinContent(xBin,yBin,1);
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

      if (binarymode) {
	if (TMath::Nint(hist->GetBinContent(xBin, yBin))==1 ) {
	  if ( TMath::Nint(hist->GetBinContent(xBin, yBin+1))==0
	       || TMath::Nint(hist->GetBinContent(xBin, yBin+1))==-3
	       || TMath::Nint(hist->GetBinContent(xBin, yBin+1))==-2   ) { //found a transition
	  xyvals.insert(make_pair( hist->GetXaxis()->GetBinCenter(xBin),hist->GetYaxis()->GetBinLowEdge(yBin+1) ));
	  foundAtThisX=true;
	  }
	}
      }
      else {
	if (hist->GetBinContent(xBin, yBin) < validcut || hist->GetBinContent(xBin, yBin+1)< validcut) continue;
	
	if (hist->GetBinContent(xBin, yBin) <cutval && hist->GetBinContent(xBin, yBin+1) >=cutval) { //found a transition
	  //HACK -- if transition is at m12 greater than 500, ignore it!
	  if (hist->GetYaxis()->GetBinLowEdge(yBin+1) >500) continue;
	  //need to keep track of the x and y values of the transition
	  xyvals.insert(make_pair( hist->GetXaxis()->GetBinCenter(xBin),hist->GetYaxis()->GetBinLowEdge(yBin+1) ));
	  foundAtThisX=true;
      }
	//ignoring reverse transitions right now
	else if (false && hist->GetBinContent(xBin, yBin) >=cutval && hist->GetBinContent(xBin, yBin+1) <cutval) {
	  cout<<"reverse transition found..."<<endl;
	}
      }      

    }
  }

  pair<double,double> leftpoint = *xyvals.begin();
  pair<double,double> rightpoint = *xyvals.rbegin();

  if (!binarymode && rightpoint.first != 2000) { cout<<file<<" "<<rightpoint.first<<endl;   rightpoint = make_pair(2000,0);  }

  if (binarymode) {
    double m12atleft = leftpoint.second;
    double m12atright = rightpoint.second;
    leftpoint = make_pair(400,m12atleft);
    rightpoint = make_pair(2050,m12atright);
  }

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
  //   makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge1btight.root","");
  //   makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge1bloose.root","");
  //   makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge2btight.root","");
  //   makeSplines("/afs/cern.ch/user/o/owen/public/RA2b/clsplots-tb40-ge2bloose.root","");
  

  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans.root","cls_value_1bloose");
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans.root","cls_value_1btight");
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans.root","cls_value_2bloose");
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans.root","cls_value_2btight");
  
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_1bloose_exp",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_1btight_exp",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_2bloose_exp",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_2btight_exp",true);

  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_1bloose_exp_minus",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_1btight_exp_minus",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_2bloose_exp_minus",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_2btight_exp_minus",true);

  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_1bloose_exp_plus",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_1btight_exp_plus",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_2bloose_exp_plus",true);
  makeSplines("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root","h2d_2btight_exp_plus",true);

  outputfile->Close();

}

void compareSimple() {

  TCanvas c1("c1");
  c1.Divide(2,1);
  c1.cd(1);

  TFile f("/afs/cern.ch/user/g/gaz/public/mSugra_scans.root");

  TFile fexp("/afs/cern.ch/user/g/gaz/public/mSugra_scans_expected.root");

  TH2S* hexp = (TH2S*) fexp.Get("h2d_1btight_exp");

  TH2D* hobs = (TH2D*) f.Get("cls_value_1btight");
  hobs->Draw("COLZ");

  int xMax = hobs->GetNbinsX();
  int yMax = hobs->GetNbinsY();

  {// for interactive root
    for(int xBin = 1; xBin <= xMax; xBin++)      {
      for(int yBin = 1; yBin <= yMax; yBin++)  {
	if (TMath::Nint(hobs->GetBinContent(xBin,yBin))==-1)    hobs->SetBinContent(xBin,yBin,0);
	else if ( hobs->GetBinContent(xBin,yBin)>=0 && hobs->GetBinContent(xBin,yBin)<0.05 ) hobs->SetBinContent(xBin,yBin,0);
	else if (  hobs->GetBinContent(xBin,yBin)>=0.05 ) hobs->SetBinContent(xBin,yBin,1);

      }
    }
  }
  hobs->SetMinimum(-1);
  hobs->SetMaximum(1);
  hobs->Draw("COLZ");
   c1.cd(2);

  {// for interactive root
    for(int xBin = 1; xBin <=  hexp->GetNbinsX(); xBin++)
      {
	for(int yBin = 1; yBin <= hexp->GetNbinsY(); yBin++)
	  {

		if ( hexp->GetBinContent(xBin,yBin) == 0 )  hexp->SetBinContent(xBin,yBin,1);
		else if ( hexp->GetBinContent(xBin,yBin) == 1 )  hexp->SetBinContent(xBin,yBin,0);

	  }
      }
  }
  hexp->SetMinimum(-1);
  hexp->SetMaximum(1);
  hexp->Draw("COLZ")
    
    TH2S* hcomb=(TH2S*)hexp->Clone("hcomb");
  hcomb->Reset();

    {// for interactive root
      for(int xBin = 1; xBin <= xMax; xBin++)      {
	for(int yBin = 1; yBin <= yMax; yBin++)  {
	  if (	  TMath::Nint(hobs->GetBinContent(xBin,yBin)) == 0 && hexp->GetBinContent(xBin,yBin)==0) hcomb->SetBinContent(xBin,yBin,1);
	  else  if (	  TMath::Nint(hobs->GetBinContent(xBin,yBin)) == 1 && hexp->GetBinContent(xBin,yBin)==1) hcomb->SetBinContent(xBin,yBin,4);
	  else  if (	 TMath::Nint( hobs->GetBinContent(xBin,yBin)) == 0 && hexp->GetBinContent(xBin,yBin)==1) hcomb->SetBinContent(xBin,yBin,2);
	  else  if (	 TMath::Nint( hobs->GetBinContent(xBin,yBin)) == 1 && hexp->GetBinContent(xBin,yBin)==0) hcomb->SetBinContent(xBin,yBin,3);
	  else 	   hcomb->SetBinContent(xBin,yBin,0);
	}
	
      }
    }
    TCanvas c2("c2");
    hcomb->SetMinimum(0);
    hcomb->SetMaximum(4);
    hcomb->Draw("COLZ");

    TFile fcurves("RA2b_tb40_exclusion.root");
    TSpline3 * curve4_ge1btight = (TSpline3*) fcurves.Get("curve4_ge1btight");

    TSpline3 * curve4_ge1btight_exp = (TSpline3*) fcurves.Get("curve4_1btight_exp");

    curve4_ge1btight->Draw("samec");
    curve4_ge1btight_exp->Draw("SAMEC");
}

//something completely new to try to combine the 4 selections into one curve.

void combineCurves() {

  //input
  TFile fin("RA2b_tb40_exclusion.root");
  TSpline3 * observed[4];
  TSpline3 * expected[4];
  TSpline3 * expected_plus[4];
  TSpline3 * expected_minus[4];

  TString stubs[4];
  stubs[0]="1bloose";
  stubs[1]="1btight";
  stubs[2]="2bloose";
  stubs[3]="2btight";

  TString basename="curve4_";
  { //interactive root
  for (int i=0; i<4; i++) {
    TString observed_name = basename+"ge"+stubs[i];
    TString expected_name = basename+stubs[i]+"_exp";
    TString expected_plus_name = expected_name+"_plus";
    TString expected_minus_name = expected_name+"_minus";

    observed[i] = (TSpline3*) fin.Get(observed_name);
    expected[i] = (TSpline3*) fin.Get(expected_name);
    expected_plus[i] = (TSpline3*) fin.Get(expected_plus_name);
    expected_minus[i] = (TSpline3*) fin.Get(expected_minus_name);
  }
  }

  //now everything is loaded

  //check
//   observed[3]->Draw("C");
//   observed[2]->Draw("Csame");
//   observed[1]->Draw("Csame");
//   observed[0]->Draw("Csame");


  //output
  TFile fout("RA2b_tb40_exclusion_combination.root","RECREATE");
  const int nsteps=100; //number of steps
  TGraph * combination_observed = new TGraph(nsteps);
  combination_observed->SetName("combination_observed");

  TGraph * expectedAtBestObserved = new TGraph(nsteps);
  expectedAtBestObserved->SetName("expectedAtBestObserved");

  TGraph * expectedAtBestObservedPlus = new TGraph(nsteps);
  expectedAtBestObservedPlus->SetName("expectedAtBestObservedPlus");


  TGraph * expectedAtBestObservedMinus = new TGraph(nsteps);
  expectedAtBestObservedMinus->SetName("expectedAtBestObservedMinus");


  //need to scan the observed, looking for the best

  //first find the limits in x (m0)
  double xmin=1e9;
  double xmax=-1e9;
  {for (int i=0; i<4; i++) {
    if ( observed[i]->GetXmin() < xmin) xmin=observed[i]->GetXmin();
    if ( observed[i]->GetXmax() > xmax) xmax=observed[i]->GetXmax();
    }}

  const double stepsize = (xmax-xmin)/ double(nsteps);

  {//interactive root
  for (int i=0; i<nsteps; i++) {

    const    double thisx = xmin + i*stepsize;
    double bestUL=0;
    int whichIsBest=-1;
    for (int j=0; j<4; j++) { //see which selection is best at thisx
      if ( observed[j]->Eval(thisx) > bestUL) {
	bestUL = observed[j]->Eval(thisx);
	whichIsBest=j;
      }
    }

    combination_observed->SetPoint(i,thisx, bestUL);

    expectedAtBestObserved->SetPoint(i,thisx, expected[whichIsBest]->Eval(thisx));
    expectedAtBestObservedPlus->SetPoint(i,thisx, expected_plus[whichIsBest]->Eval(thisx));
    expectedAtBestObservedMinus->SetPoint(i,thisx, expected_minus[whichIsBest]->Eval(thisx));

  }
  }

  combination_observed->GetHistogram()->SetMinimum(100);
  combination_observed->GetHistogram()->SetMaximum(500);

  combination_observed->Draw("AC");
  expectedAtBestObserved->SetLineColor(kRed);
  expectedAtBestObserved->SetMarkerColor(kRed);
  expectedAtBestObserved->Draw("*");
  expectedAtBestObservedPlus->SetLineColor(kBlue);
  expectedAtBestObservedPlus->SetMarkerColor(kBlue);
  expectedAtBestObservedPlus->Draw("*");
  expectedAtBestObservedMinus->SetLineColor(kBlue);
  expectedAtBestObservedMinus->SetMarkerColor(kBlue);
  expectedAtBestObservedMinus->Draw("*");
  
  fout.Close();

}
