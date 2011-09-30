TList* getContours(const TH2* hist,double contourLevel,const TString filename="");

void smoothHistAcross(TH2* hist,int xMin=0);

TH2* smoothHoles(const TH2* originalHist);

//TGraph* RA2b_limit(const TString filename, const TString histname);
//TSpline3* RA2b_limit(const TString filename, const TString histname);

TSpline3* getCLs1080ObsNLOtb40(){

   Int_t nl = 10;
   Double_t xl[10]={380,500,700,900,1100,1300,1500,1700,1900,2000};
   Double_t yl[10]={520,490,390,290, 230, 210, 205, 200, 195, 190};
   Double_t exl[10];
   Double_t eyl[10];

   TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
   //  gr1->SetMarkerColor(kRed);
   //  gr1->SetMarkerStyle(21);

   TSpline3 *s = new TSpline3("grs",gr1);
   s->SetLineColor(kRed);
   s->SetLineStyle(1);
   s->SetLineWidth(3);

   return s;
}

TGraph* get_RA2b_1bloose() {
    Int_t nl = 11;
    Double_t xl[11]={385, 430, 555, 615, 715, 815, 875, 980, 1390, 1660, 2005};
    Double_t yl[11]={385, 370, 327, 300, 270, 250, 235, 215, 206,  195,  190 };
    
    TGraph* gr_obs_1bloose = new TGraph(nl,xl,yl);
    return gr_obs_1bloose;

}
TGraph* get_RA2b_1btight() {

    Int_t nl = 14;
    Double_t xl[14]={380,500,555, 600, 650, 750, 875, 975, 1035, 1185, 1390, 1560, 1630, 2005};
    Double_t yl[14]={410,405,398, 390, 370, 345, 290, 255, 245,  217,  209,  205, 197,  190};
//     Double_t exl[14];
//     Double_t eyl[14];
    
    TGraph* gr_obs_1btight = new TGraph(nl,xl,yl);

    return gr_obs_1btight;
}
TGraph* get_RA2b_2bloose() {
    Int_t nl = 12;
    Double_t xl[12]={380, 430, 500,630, 700, 730, 755, 850, 975, 1110, 1480, 2005};
    Double_t yl[12]={370, 365, 350,343,  330,318, 295, 270, 250, 232,  215,  210};
    
    TGraph* gr_obs_2bloose = new TGraph(nl,xl,yl);
  
    return gr_obs_2bloose;

}
TGraph* get_RA2b_2btight() {
    Int_t nl = 13;
    Double_t xl[13]={380,  435, 580, 665, 725, 770, 790, 810,850, 910, 1113, 1245, 2005};
    Double_t yl[13]={412,  410, 400, 383, 350, 330, 310, 257,229, 215, 200,  193, 190};
    Double_t exl[13];
    Double_t eyl[13];
    
    TGraph* gr_obs_2btight = new TGraph(nl,xl,yl);
    return gr_obs_2btight;

}


TGraph* get_RA2b_1bloose_exp() {
  //not drawn yet
  return get_RA2b_1bloose();
}
TGraph* get_RA2b_1btight_exp() {
   Int_t nl = 14;
    Double_t xl[14]={380, 600, 700, 770, 800, 875, 935, 1040, 1120, 1220, 1620, 1670, 1900, 2000};
    Double_t yl[14]={410, 405, 385, 355, 345, 301, 285, 255,   244, 225,  203,  198,  193, 190};
    
    TGraph* gr_exp = new TGraph(nl,xl,yl);

    return gr_exp;
}
TGraph* get_RA2b_2bloose_exp() {
  //not drawn yet
  return get_RA2b_1bloose();

}
TGraph* get_RA2b_2btight_exp() {
  //not drawn yet
  return get_RA2b_1bloose();

}


TGraph* get_RA2b_1bloose_exp_p() {
  //not drawn yet
  return get_RA2b_1bloose();

}
TGraph* get_RA2b_1btight_exp_p() {
    Int_t nl = 13;
    Double_t xl[13]={380, 575, 705, 755, 815, 880, 975, 1050,1150,1360,1490,1900,2000 };
    Double_t yl[13]={390, 380, 340, 304, 285, 250, 230, 216, 208, 190, 178, 170, 170};
    
    TGraphErrors* gr_exp_p = new TGraphErrors(nl,xl,yl);
    return gr_exp_p;

}
TGraph* get_RA2b_2bloose_exp_p() {
  //not drawn yet
  return get_RA2b_1bloose();

}
TGraph* get_RA2b_2btight_exp_p() {
  //not drawn yet
  return get_RA2b_1bloose();

}

TGraph* get_RA2b_1bloose_exp_m() {
  //not drawn yet
  return get_RA2b_1bloose();

}
TGraph* get_RA2b_1btight_exp_m() {
    Int_t nl = 13;
    Double_t xl[13]={380,605,675,800,880,935,1030,1135,1315,1415,1600,1915,2000 };
    Double_t yl[13]={450,438,417,385,355,326,300, 278 ,264, 250, 231, 221, 210};
    
    TGraph* gr_exp_m = new TGraph(nl,xl,yl);
    return gr_exp_m;
}
TGraph* get_RA2b_2bloose_exp_m() {
  //not drawn yet
  return get_RA2b_1bloose();

}
TGraph* get_RA2b_2btight_exp_m() {
  //not drawn yet
  return get_RA2b_1bloose();

}


TGraph* SSdilep_NLO( ){


   TGraph *graph = new TGraph(178);
   graph->SetPoint(0,5,364.256);
   graph->SetPoint(1,5.0037,364.3236);
   graph->SetPoint(2,5.667515,376.4613);
   graph->SetPoint(3,6.303097,388.0826);
   graph->SetPoint(4,6.954213,399.9881);
   graph->SetPoint(5,7.514941,410.9272);
   graph->SetPoint(6,7.646183,410.3409);
   graph->SetPoint(7,8.223984,416.1515);
   graph->SetPoint(8,8.295899,420.6026);
   graph->SetPoint(9,8.468629,424.8871);
   graph->SetPoint(10,10.88691,430.2015);
   graph->SetPoint(11,15,436.2527);
   graph->SetPoint(12,25,443.197);
   graph->SetPoint(13,35,447.5647);
   graph->SetPoint(14,45,449.7784);
   graph->SetPoint(15,55,450.4723);
   graph->SetPoint(16,57.26605,450.6156);
   graph->SetPoint(17,65,449.8667);
   graph->SetPoint(18,72.54146,448.1337);
   graph->SetPoint(19,75,447.1522);
   graph->SetPoint(20,85,443.5753);
   graph->SetPoint(21,91.66453,439.6889);
   graph->SetPoint(22,95,437.1199);
   graph->SetPoint(23,100.9103,431.6382);
   graph->SetPoint(24,105,427.0075);
   graph->SetPoint(25,113.3399,414.8456);
   graph->SetPoint(26,113.8718,412.9956);
   graph->SetPoint(27,114.4743,410.2333);
   graph->SetPoint(28,115,407.009);
   graph->SetPoint(29,115,407.009);
   graph->SetPoint(30,115,407.009);
   graph->SetPoint(31,115.5612,402.1022);
   graph->SetPoint(32,115.9149,400.0791);
   graph->SetPoint(33,116.645,397.8855);
   graph->SetPoint(34,116.7666,396.362);
   graph->SetPoint(35,125,386.6004);
   graph->SetPoint(36,135,377.4256);
   graph->SetPoint(37,143.7285,371.2169);
   graph->SetPoint(38,145,370.4765);
   graph->SetPoint(39,146.6675,369.1417);
   graph->SetPoint(40,155,365.0869);
   graph->SetPoint(41,165,360.3886);
   graph->SetPoint(42,165.6767,359.4622);
   graph->SetPoint(43,169.1272,357.7276);
   graph->SetPoint(44,175,355.2705);
   graph->SetPoint(45,185,350.2034);
   graph->SetPoint(46,195,345.0257);
   graph->SetPoint(47,197.3941,344.2251);
   graph->SetPoint(48,205,340.8926);
   graph->SetPoint(49,206.6498,339.8847);
   graph->SetPoint(50,215,336.803);
   graph->SetPoint(51,225,334.2859);
   graph->SetPoint(52,235,332.1803);
   graph->SetPoint(53,245,330.4912);
   graph->SetPoint(54,255,328.906);
   graph->SetPoint(55,265,327.7156);
   graph->SetPoint(56,275,326.7611);
   graph->SetPoint(57,281.4439,326.2159);
   graph->SetPoint(58,285,325.9328);
   graph->SetPoint(59,295,325.3143);
   graph->SetPoint(60,305,324.8972);
   graph->SetPoint(61,313.7252,324.9344);
   graph->SetPoint(62,315,325.2);
   //graph->SetPoint(63,315.2418,324.9281);
   graph->SetPoint(63,320,324.9281);
   graph->SetPoint(64,325,325.2656);
   graph->SetPoint(65,335,326.1206);
   graph->SetPoint(66,345,327.2148);
   graph->SetPoint(67,345.6352,327.3103);
   graph->SetPoint(68,355,328.3962);
   graph->SetPoint(69,365,329.5845);
   graph->SetPoint(70,375,330.7801);
   graph->SetPoint(71,385,331.9823);
   graph->SetPoint(72,395,333.2658);
   graph->SetPoint(73,405,334.7468);
   graph->SetPoint(74,415,336.3387);
   graph->SetPoint(75,418.2803,336.9342);
   graph->SetPoint(76,425,338.0564);
   graph->SetPoint(77,435,339.7781);
   graph->SetPoint(78,445,341.6835);
   graph->SetPoint(79,455,343.6595);
   graph->SetPoint(80,456.2874,343.8833);
   graph->SetPoint(81,465,345.5794);
   graph->SetPoint(82,475,347.3702);
   graph->SetPoint(83,485,349.138);
   graph->SetPoint(84,495,350.8699);
   graph->SetPoint(85,505,352.7064);
   graph->SetPoint(86,513.3109,354.2605);
   graph->SetPoint(87,515,354.6418);
   graph->SetPoint(88,516.8738,354.9615);
   graph->SetPoint(89,525,356.4621);
   graph->SetPoint(90,527.6248,356.9174);
   graph->SetPoint(91,535,358.3616);
   graph->SetPoint(92,545,360.0272);
   graph->SetPoint(93,555,360.9706);
   graph->SetPoint(94,563.2791,360.5442);
   graph->SetPoint(95,565,360.0536);
   graph->SetPoint(96,565.6864,360.458);
   graph->SetPoint(97,575,358.6743);
   graph->SetPoint(98,585,355.5729);
   graph->SetPoint(99,595,351.7765);
   graph->SetPoint(100,601.5435,349.1743);
   graph->SetPoint(101,602.3684,348.8774);
   graph->SetPoint(102,605,348.0353);
   graph->SetPoint(103,605,348.0353);
   graph->SetPoint(104,605,348.0353);
   graph->SetPoint(105,610.4187,346.9274);
   graph->SetPoint(106,615,344.9838);
   graph->SetPoint(107,615,344.9838);
   graph->SetPoint(108,615,344.9838);
   graph->SetPoint(109,619.5339,342.6968);
   graph->SetPoint(110,625,337.7304);
   graph->SetPoint(111,625,337.7304);
   graph->SetPoint(112,625,337.7304);
   graph->SetPoint(113,625.6453,337.1092);
   graph->SetPoint(114,626.6194,335.3217);
   graph->SetPoint(115,627.4621,333.9927);
   graph->SetPoint(116,635,327.6697);
   graph->SetPoint(117,645,320.4765);
   graph->SetPoint(118,655,313.9563);
   graph->SetPoint(119,665,307.91);
   graph->SetPoint(120,675,302.8112);
   graph->SetPoint(121,685,298.2886);
   graph->SetPoint(122,686.3657,297.4847);
   graph->SetPoint(123,691.7435,295.3889);
   graph->SetPoint(124,692.0528,295.3735);
   graph->SetPoint(125,695,294.1766);
   graph->SetPoint(126,702.0545,291.7193);
   graph->SetPoint(127,705,290.8017);
   graph->SetPoint(128,705,290.8017);
   graph->SetPoint(129,705,290.8017);
   graph->SetPoint(130,706.6572,290.097);
   graph->SetPoint(131,715,287.8264);
   graph->SetPoint(132,718.1839,286.8499);
   graph->SetPoint(133,720.5234,286.3812);
   graph->SetPoint(134,725,285.2411);
   graph->SetPoint(135,733.3917,283.3605);
   graph->SetPoint(136,735,283.1243);
   graph->SetPoint(137,745,282.3018);
   graph->SetPoint(138,746.1148,282.1957);
   graph->SetPoint(139,752.2693,281.9635);
   graph->SetPoint(140,755,281.7206);
   graph->SetPoint(141,765,281.4361);
   graph->SetPoint(142,775,280.6284);
   graph->SetPoint(143,785,279.3488);
   graph->SetPoint(144,791.7576,277.9269);
   graph->SetPoint(145,795,277.3645);
   graph->SetPoint(146,799.8055,276.0494);
   graph->SetPoint(147,805,274.7065);
   graph->SetPoint(148,814.4029,272.8767);
   graph->SetPoint(149,815,272.8738);
   //graph->SetPoint(150,815.2259,272.3869);
   graph->SetPoint(150,820.2259,272.3869);
   graph->SetPoint(151,825,271.7734);
   graph->SetPoint(152,835,271.6559);
   graph->SetPoint(153,837.2917,271.497);
   graph->SetPoint(154,845,271.0756);
   graph->SetPoint(155,855,270.1055);
   graph->SetPoint(156,859.854,269.3756);
   graph->SetPoint(157,865,268.2217);
   graph->SetPoint(158,875,265.7956);
   graph->SetPoint(159,882.4138,263.9001);
   graph->SetPoint(160,885,263.4806);
   graph->SetPoint(161,894.5582,260.842);
   graph->SetPoint(162,895,260.4761);
   graph->SetPoint(163,905,258.374);
   graph->SetPoint(164,915,256.9917);
   graph->SetPoint(165,919.6315,256.5705);
   graph->SetPoint(166,925,256.4248);
   graph->SetPoint(167,935,256.0899);
   graph->SetPoint(168,942.9671,255.9357);
   graph->SetPoint(169,945,255.9691);
   graph->SetPoint(170,955,256.898);
   graph->SetPoint(171,962.3684,257.9652);
   graph->SetPoint(172,965,258.3537);
   graph->SetPoint(173,975,260.1838);
   graph->SetPoint(174,985,262.0139);
   graph->SetPoint(175,994.4399,263.7415);
   graph->SetPoint(176,995,263.844);
   graph->SetPoint(177,1001,264.244);

  return graph;

}

TGraphErrors* OSdilep_NLO( ){

  
  const unsigned int nNLO=50;

  Double_t xNLO[nNLO],yNLO[nNLO],xerr[nNLO],yerr[nNLO];
  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xNLO[ierr] = 0.;
    yNLO[ierr] = 0.;
  }

  Int_t i = -1;

  xNLO[++i]=0.;    yNLO[i]=460.; //0
  xNLO[++i]=80;    yNLO[i]=410.; //1
  xNLO[++i]=100;   yNLO[i]=365.; //2
  xNLO[++i]=120;   yNLO[i]=360.; //3
  xNLO[++i]=150;   yNLO[i]=355.; //4
  xNLO[++i]=180;   yNLO[i]=325.; //5
  xNLO[++i]=200;   yNLO[i]=290.; //6
  xNLO[++i]=210;   yNLO[i]=285.; //7
  xNLO[++i]=270;   yNLO[i]=300.; //8
  xNLO[++i]=350;   yNLO[i]=310.; //9
  xNLO[++i]=450;   yNLO[i]=310.; //10
  xNLO[++i]=520;   yNLO[i]=320.; //11
  xNLO[++i]=550;   yNLO[i]=280.; //12
  xNLO[++i]=600;   yNLO[i]=265.; //13
  xNLO[++i]=700;   yNLO[i]=250.; //14
  xNLO[++i]=900;   yNLO[i]=230.; //15
  xNLO[++i]=1000;  yNLO[i]=220.; //16
  xNLO[++i]=1400;  yNLO[i]=180.; //17
  xNLO[++i]=1800;  yNLO[i]=156.; //18
  xNLO[++i]=1900;  yNLO[i]=150.; //19

  TGraphErrors* grtb10  = new TGraphErrors(i,xNLO, yNLO,xerr,yerr);
  grtb10->SetMarkerColor(kWhite);
  //grtb10->SetMarkerStyle(21);
  return grtb10;

}



TGraphErrors* RA1_NLO( ){

  const unsigned int nNLO=50;

  Double_t xNLO[nNLO],yNLO[nNLO],xerr[nNLO],yerr[nNLO];
  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xNLO[ierr] = 0.;
    yNLO[ierr] = 0.;
  }

  Double_t xl[12]={  0,100,300,500,700,900,1100,1300,1500,1700,1900,2000};
  Double_t yl[12]={540,535,520,490,400,280, 240, 220, 210, 200, 195, 190};

  TGraphErrors* grtb10  = new TGraphErrors(12,xl, yl,xerr,yerr);
  grtb10->SetMarkerColor(kWhite);
  //grtb10->SetMarkerStyle(21);
  return grtb10;


}


TGraphErrors* getRA1Observed_NLO_tanBeta10(){

  Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
    
  xl[0] = 0;
  yl[0] = 288;
  xl[1] = 100;
  yl[1] = 280;
  xl[2] = 150;
  yl[2] = 275;
  xl[3] = 200;
  yl[3] = 268;
  xl[4] = 250;
  yl[4] = 260;
  xl[5] = 300;
  yl[5] = 237;
  xl[6] = 350;
  yl[6] = 203;
  xl[7] = 400;
  yl[7] = 172;
  xl[8] = 450;
  yl[8] = 156;
  xl[9] = 500;
  yl[9] = 147;
  xl[10] = 670;
  yl[10] = 120;

  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
    
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  // s->SetLineStyle(2);
  s->SetLineWidth(3);
  
  return gr1;

}

TGraphErrors* getRA5Observed_NLO_tanBeta10(){

  int nl = 158;
   TGraphErrors *graph = new TGraphErrors(nl);
   graph->SetName("CLsNLO_ObservedTb10");
   graph->SetTitle("Graph");
   graph->SetLineWidth(3);
   for (int i=0; i<nl; i++) { graph->SetPointError(i,0,0); }
   graph->SetPoint(0,5,147.0268);
   graph->SetPoint(1,5.66862,155);
   graph->SetPoint(2,5.215185,165);
   graph->SetPoint(3,7.321363,175);
   graph->SetPoint(4,8.711564,185);
   graph->SetPoint(5,11.52623,195);
   graph->SetPoint(6,13.68512,205);
   graph->SetPoint(7,15,209.0835);
   graph->SetPoint(8,23.38369,215);
   graph->SetPoint(9,25,218.2312);
   graph->SetPoint(10,29.60877,215);
   graph->SetPoint(11,35,213.8296);
   graph->SetPoint(12,45,208.6325);
   graph->SetPoint(13,55,209.5333);
   graph->SetPoint(14,65,213.3807);
   graph->SetPoint(15,74.51391,215);
   graph->SetPoint(16,75,215.1255);
   graph->SetPoint(17,85,215.1879);
   graph->SetPoint(18,85.57727,215);
   graph->SetPoint(19,90.5631,213);
   graph->SetPoint(20,92.8214,211);
   graph->SetPoint(21,95,211);
   graph->SetPoint(22,97.6207,210);
   graph->SetPoint(23,105,200.2214);
   graph->SetPoint(24,115,203.8768);
   graph->SetPoint(25,122.0319,205);
   graph->SetPoint(26,125,205.8696);
   graph->SetPoint(27,135,207.8686);
   graph->SetPoint(28,145,206.6887);
   graph->SetPoint(29,146.963,205);
   graph->SetPoint(30,146.5303,195);
   graph->SetPoint(31,155,185.3473);
   graph->SetPoint(32,155.2025,185);
   graph->SetPoint(33,157.8652,175);
   graph->SetPoint(34,162.3215,165);
   graph->SetPoint(35,165,162.8418);
   graph->SetPoint(36,175,158.8367);
   graph->SetPoint(37,184.78,155);
   graph->SetPoint(38,185,154.9607);
   graph->SetPoint(39,185.0875,155);
   graph->SetPoint(40,193.7749,165);
   graph->SetPoint(41,195,166.1718);
   graph->SetPoint(42,205,169.5398);
   graph->SetPoint(43,215,174.8403);
   graph->SetPoint(44,225,174.0844);
   graph->SetPoint(45,235,173.9268);
   graph->SetPoint(46,245,173.3);
   graph->SetPoint(47,255,172.9746);
   graph->SetPoint(48,265,173.4697);
   graph->SetPoint(49,269.2907,175);
   graph->SetPoint(50,275,177.7745);
   graph->SetPoint(51,285,178.3641);
   graph->SetPoint(52,292.9614,185);
   graph->SetPoint(53,295,186.4867);
   graph->SetPoint(54,305,188.7906);
   graph->SetPoint(55,315,187.0713);
   graph->SetPoint(56,321.1138,185);
   graph->SetPoint(57,319.446,187);
   graph->SetPoint(58,315,187);
   graph->SetPoint(59,309.241,188);
   graph->SetPoint(60,310.12,188);
   graph->SetPoint(61,315,187.071);
   graph->SetPoint(62,325,184);
   graph->SetPoint(63,331.508,178);
   graph->SetPoint(64,335,179);
   graph->SetPoint(65,336.632,179);
   graph->SetPoint(66,338.563,179);
   graph->SetPoint(67,335,179);
   graph->SetPoint(68,331.094,178);
   graph->SetPoint(69,335,179);
   graph->SetPoint(70,343.2348,175);
   graph->SetPoint(71,345,175);
   graph->SetPoint(72,346.3166,175);
   graph->SetPoint(73,355,167.6639);
   graph->SetPoint(74,358.209,165);
   graph->SetPoint(75,355,158.0992);
   graph->SetPoint(76,353.3595,155);
   graph->SetPoint(77,355,149.5683);
   graph->SetPoint(78,365,149.4669);
   graph->SetPoint(79,375,154.471);
   graph->SetPoint(80,385,150.9277);
   graph->SetPoint(81,389.7334,145);
   graph->SetPoint(82,388.716,135);
   graph->SetPoint(83,394.9455,125);
   graph->SetPoint(84,395,124.9927);
   graph->SetPoint(85,405,119.434);
   graph->SetPoint(86,415,120.7875);
   graph->SetPoint(87,419.526,115);
   graph->SetPoint(88,425,111.8557);
   graph->SetPoint(89,427.6512,115);
   graph->SetPoint(90,435,124.8554);
   graph->SetPoint(91,435.6198,125);
   graph->SetPoint(92,445,128.0318);
   graph->SetPoint(93,453.2803,135);
   graph->SetPoint(94,455,136.4788);
   graph->SetPoint(95,457.491,135);
   graph->SetPoint(96,460.279,125);
   graph->SetPoint(97,455,120.1707);
   graph->SetPoint(98,450.4731,115);
   graph->SetPoint(99,448.9584,105);
   graph->SetPoint(100,445,102.8302);
   graph->SetPoint(101,435,99.6762);
   graph->SetPoint(102,425,102.0901);
   graph->SetPoint(103,415,100.1085);
   graph->SetPoint(104,405,98.65763);
   graph->SetPoint(105,395,97.07405);
   graph->SetPoint(106,385,100.4358);
   graph->SetPoint(107,375,101.6786);
   graph->SetPoint(108,365,97.2692);
   graph->SetPoint(109,355,98.03454);
   graph->SetPoint(110,345,102.5535);
   graph->SetPoint(111,335,104.6001);
   graph->SetPoint(112,332.1912,105);
   graph->SetPoint(113,325,106.657);
   graph->SetPoint(114,315,108.8675);
   graph->SetPoint(115,305,106.1756);
   graph->SetPoint(116,301.0824,105);
   graph->SetPoint(117,295,104.0053);
   graph->SetPoint(118,285,103.9421);
   graph->SetPoint(119,275,99.44182);
   graph->SetPoint(120,265,96.80255);
   graph->SetPoint(121,255,98.85424);
   graph->SetPoint(122,245,96.8166);
   graph->SetPoint(123,235,98.02971);
   graph->SetPoint(124,225,98.35997);
   graph->SetPoint(125,215,98.58589);
   graph->SetPoint(126,205,98.85166);
   graph->SetPoint(127,195,97.31524);
   graph->SetPoint(128,188.6431,95);
   graph->SetPoint(129,185,92.20591);
   graph->SetPoint(130,182.2059,95);
   graph->SetPoint(131,175,98.48925);
   graph->SetPoint(132,165,96.35012);
   graph->SetPoint(133,155,97.00039);
   graph->SetPoint(134,148.8757,95);
   graph->SetPoint(135,145,92.36106);
   graph->SetPoint(136,135,91.74712);
   graph->SetPoint(137,125,90.40189);
   graph->SetPoint(138,115.8813,95);
   graph->SetPoint(139,115,95.12141);
   graph->SetPoint(140,114.2028,95);
   graph->SetPoint(141,105,90.12921);
   graph->SetPoint(142,98.53317,95);
   graph->SetPoint(143,95,95.78995);
   graph->SetPoint(144,91.21521,95);
   graph->SetPoint(145,85,90.39959);
   graph->SetPoint(146,75,87.31825);
   graph->SetPoint(147,65,88.53508);
   graph->SetPoint(148,55,91.0927);
   graph->SetPoint(149,49.05011,95);
   graph->SetPoint(150,45,95.8049);
   graph->SetPoint(151,41.9025,95);
   graph->SetPoint(152,35,90.06895);
   graph->SetPoint(153,30.06895,95);
   graph->SetPoint(154,25,98.00605);
   graph->SetPoint(155,15,100.5421);
   graph->SetPoint(156,10.54208,105);
   graph->SetPoint(157,5,109.6184);
   
   //TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  graph->SetMarkerColor(kBlue);
  graph->SetMarkerStyle(21);
    
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",graph);
  s->SetLineColor(kRed);
  // s->SetLineStyle(2);
  s->SetLineWidth(3);
  
  return graph;

}


TGraphErrors* getRA6Observed_NLO_tanBeta10(){

  Int_t nNLO=47;
  Double_t xNLO[47],yNLO[47],xerr[47],yerr[47];
  for( unsigned int ierr = 0 ; ierr < nNLO ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xNLO[ierr] = 0.;
    yNLO[ierr] = 0.;
  }

  Int_t i = -1;

  xNLO[++i]=0;
  yNLO[i]=225.;

  xNLO[++i]=20;
  yNLO[i]=225.;

  xNLO[++i]=30.;
  yNLO[i]=235.;

  xNLO[++i]=40.;
  yNLO[i]=235.;

  xNLO[++i]=50.;
  yNLO[i]=235.;

  xNLO[++i]=60.;
  yNLO[i]=225.;

  xNLO[++i]=70.;
  yNLO[i]=235.;

  xNLO[++i]=80.;
  yNLO[i]=245.;

  xNLO[++i]=90.;
  yNLO[i]=245.;

  xNLO[++i]=100.;
  yNLO[i]=245.;

  xNLO[++i]=115.;
  yNLO[i]=240.;

  xNLO[++i]=110.;
  yNLO[i]=235.;

  xNLO[++i]=105.;
  yNLO[i]=230.;

  xNLO[++i]=105.;
  yNLO[i]=213.;

  xNLO[++i]=105.;
  yNLO[i]=205.;

  xNLO[++i]=95.;
  yNLO[i]=200.;

  xNLO[++i]=90.;
  yNLO[i]=185.;

  xNLO[++i]=70.;
  yNLO[i]=155.;

  xNLO[++i]=65.;
  yNLO[i]=145.;

  xNLO[++i]=75.;
  yNLO[i]=155.;

  xNLO[++i]=90.;
  yNLO[i]=170.;

  xNLO[++i]=100.;
  yNLO[i]=185.;

  xNLO[++i]=110.;
  yNLO[i]=195.;

  xNLO[++i]=120.;
  yNLO[i]=205.;

  xNLO[++i]=130.;
  yNLO[i]=215.;

  xNLO[++i]=140.;
  yNLO[i]=215.;

  xNLO[++i]=150.;
  yNLO[i]=215.;

  xNLO[++i]=160.;
  yNLO[i]=185.;

  xNLO[++i]=170.;
  yNLO[i]=175.;

  xNLO[++i]=170.;
  yNLO[i]=155.;

  xNLO[++i]=160.;
  yNLO[i]=140.;

  xNLO[++i]=170.;
  yNLO[i]=135.;

  xNLO[++i]=180.;
  yNLO[i]=125.;

  TGraphErrors* gr1 = new TGraphErrors(i++,xNLO, yNLO,xerr,yerr);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
    
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  // s->SetLineStyle(2);
  s->SetLineWidth(3);
  
  return gr1;

}


