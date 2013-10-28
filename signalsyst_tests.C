{

  TCut masscut="m0==325&&m12==1&&higgsMbb1MassDiff_correct==2";

  TChain t0("reducedTree");
  TChain tjerb("reducedTree");
  TChain tjesu("reducedTree");
  TChain tjesd("reducedTree");
  TChain tjeru("reducedTree");
  TChain tjerd("reducedTree");

  TString path = "/cu6/joshmt/reducedTrees/v71_5c/";

  t0.Add(path+"reducedTree.JES0_JER0_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71.root");
  tjerb.Add(path+"reducedTree.JES0_JERbias_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71.root");
  tjeru.Add(path+"reducedTree.JES0_JERup_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71.root");
  tjerd.Add(path+"reducedTree.JES0_JERdown_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71.root");

  tjesu.Add(path+"reducedTree.JESup_JERbias_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71.root");
  tjesd.Add(path+"reducedTree.JESdown_JERbias_PFMETTypeI_METunc0_PUunc0_hpt20.SMS-TChiHH_2b2b_2J_mChargino-130to325_mLSP-1to195_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM-v1_AODSIM_UCSB1872_v71.root");

  TCanvas *c1 = new TCanvas("c1","transparent pad",600,600);
  t0.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>h0(40,70,170)",masscut);
  tjerb.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>hjerb(40,70,170)",masscut,"sames");
  hjerb.SetLineColor(kGreen+2);
  c1->Update();
  TPaveStats *statbox_0 = (TPaveStats*)h0->GetListOfFunctions()->FindObject("stats");
  TPaveStats *statbox_jerb = (TPaveStats*)hjerb->GetListOfFunctions()->FindObject("stats");
  statbox_jerb->SetY1NDC(0.4); statbox_jerb->SetY2NDC(0.7); statbox_jerb->SetTextColor(kGreen+2);


  hjerb.Draw();
  tjeru.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>hjeru(40,70,170)",masscut,"sames");
  tjerd.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>hjerd(40,70,170)",masscut,"sames");
  TPaveStats *statbox_jerd = (TPaveStats*)hjerd->GetListOfFunctions()->FindObject("stats");
  TPaveStats *statbox_jeru = (TPaveStats*)hjeru->GetListOfFunctions()->FindObject("stats");
  statbox_jeru->SetY1NDC(0.7); statbox_jeru->SetY2NDC(1.0); statbox_jeru->SetTextColor(kRed);
  statbox_jerd->SetY1NDC(0.1); statbox_jerd->SetY2NDC(0.4); statbox_jerd->SetTextColor(kBlue);

  hjeru->SetLineColor(kRed);
  hjerd->SetLineColor(kBlue);
  c1->Update();
  c1->SaveAs("higgs_masspeak_jervaried.pdf");

  //pause here.

  hjerb.Draw(); //use JERbias for nominal
  tjesu.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>hjesu(40,70,170)",masscut,"sames");
  tjesd.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>hjesd(40,70,170)",masscut,"sames");
  hjesu->SetLineColor(kRed);
  hjesd->SetLineColor(kBlue);

  statbox_jerb = (TPaveStats*)hjerb->GetListOfFunctions()->FindObject("stats");
  statbox_jerb->SetY1NDC(0.7);  statbox_jerb->SetY2NDC(0.99);
  TPaveStats *statbox_jesu = (TPaveStats*)hjesu->GetListOfFunctions()->FindObject("stats");
  TPaveStats *statbox_jesd = (TPaveStats*)hjesd->GetListOfFunctions()->FindObject("stats");

  statbox_jesu->SetY1NDC(0.4); statbox_jesu->SetY2NDC(0.7); statbox_jesu->SetTextColor(kRed);
  statbox_jesd->SetY1NDC(0.1); statbox_jesd->SetY2NDC(0.4); statbox_jesd->SetTextColor(kBlue);

  c1->SaveAs("higgs_masspeak_jesvaried.pdf");

  // == METsig == 
  gStyle->SetOptStat(0);
  masscut="m0==325&&m12==1";
  tjerb.Draw("METsig>>hmetsig0(30,0,300)",masscut);
  hmetsig0->SetMinimum(0);
  hmetsig0->SetXTitle("MET significance");
  tjesu.Draw("METsig>>hmetsig_jesu(30,0,300)",masscut,"same");
  tjesd.Draw("METsig>>hmetsig_jesd(30,0,300)",masscut,"same");
  hmetsig_jesu->SetLineColor(kRed);
  hmetsig_jesd->SetLineColor(kBlue);

  c1->SaveAs("METsig_jesvaried.pdf");

  return;

  hmetsig0.Draw();
  tjerb.Draw("METsig>>hmetsig_jerb(30,0,300)",masscut,"same");
  hmetsig_jerb.SetLineColor(kGreen+2);

  hmetsig_jerb.Draw();
  tjeru.Draw("METsig>>hmetsig_jeru(30,0,300)",masscut,"same");
  tjerd.Draw("METsig>>hmetsig_jerd(30,0,300)",masscut,"same");
  hmetsig_jeru.SetLineColor(kRed);
  hmetsig_jerd.SetLineColor(kBlue);


  // == njets20 ==
  masscut="m0==325&&m12==1";
  t0.Draw("njets20>>hnjets0(9,0,9)",masscut);
  tjesu.Draw("njets20>>hnjets_jesu(9,0,9)",masscut,"same");
  tjesd.Draw("njets20>>hnjets_jesd(9,0,9)",masscut,"same");
  hnjets_jesu.SetLineColor(kRed);
  hnjets_jesd.SetLineColor(kBlue);

  hnjets0.Draw();
  tjerb.Draw("njets20>>hnjets_jerb(9,0,9)",masscut,"same");
  hnjets_jerb.SetLineColor(kGreen+2);

  hnjets_jerb.Draw();
  tjeru.Draw("njets20>>hnjets_jeru(9,0,9)",masscut,"same");
  tjerd.Draw("njets20>>hnjets_jerd(9,0,9)",masscut,"same");
  hnjets_jeru.SetLineColor(kRed);
  hnjets_jerd.SetLineColor(kBlue);

  // Pileup
  gStyle->SetOptStat(1111111);
  masscut="PUweight*(m0==325&&m12==1&&higgsMbb1MassDiff_correct==2)";
  TCut  masscutw="PUweightSystVar*(m0==325&&m12==1&&higgsMbb1MassDiff_correct==2)";
  t0.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>h0pu(40,70,170)",masscut);
  t0.Draw("0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>>h0pusyst(40,70,170)",masscutw,"sames");
  h0pusyst->SetLineColor(kMagenta);


  gStyle->SetOptStat(0);  
  masscut = "PUweight*(m0==325&&m12==1)";
  masscutw="PUweightSystVar*(m0==325&&m12==1)";
  t0.Draw("METsig>>hmetsig_pu(30,0,300)",masscut);
  t0.Draw("METsig>>hmetsig_pusyst(30,0,300)",masscutw,"sames");
  hmetsig_pusyst->SetLineColor(kMagenta);

  // == check LMT versus the "tweaking" implementation
  tjerb.Draw("2*(nbtag2_nomSF==1)+3*(nbtag3_nomSF==1)+4*(nbtag4_nomSF==1)>>hNomTweak(4,0.5,4.5)","m0==325&&m12==1&&njets20>=4&&njets20<=5&&METsig>=30&&jetpt2>=50 && nbtag0_nomSF==0");
  tjerb.Draw("2*(nbtag2_rawMC==1)+3*(nbtag3_rawMC==1)+4*(nbtag4_rawMC==1)>>hRaw(4,0.5,4.5)","m0==325&&m12==1&&njets20>=4&&njets20<=5&&METsig>=30&&jetpt2>=50 &&nbtag0_rawMC==0");
  tjerb.Draw("2*(nbtag2_rawMC==1)+3*(nbtag3_rawMC==1)+4*(nbtag4_rawMC==1)>>hLMT(4,0.5,4.5)","BTagWeightLMT*(m0==325&&m12==1&&njets20>=4&&njets20<=5&&METsig>=30&&jetpt2>=50 &&nbtag0_rawMC==0)");

  hRaw->Draw();
  hRaw->SetLineColor(kBlack);
  hNomTweak->SetLineColor(kGreen);
  hLMT->SetLineColor(kBlue);
  hNomTweak->Draw("SAME");
  hLMT->Draw("SAME");

}
