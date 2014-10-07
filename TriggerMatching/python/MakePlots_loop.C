#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include <iostream>
#include "TROOT.h"

int MakePlots_loop(){
 
  //gROOT->ProcessLine(".L ~/bin/style-CMSTDR.C");
  //gROOT->ProcessLine("setTDRStyle()");
  TH1::SetDefaultSumw2(true);


  TFile *f1 = new TFile("analyzeRECO_ALLEventsHWW_tightWP.root","READ");
  TTree *eventsTree = (TTree*)f1->Get("TriggerMatching/eventsTree");

  //eventsTree->MakeCode("TreeBranches.C");
  //return 0;

  //////////////////////////////////////////////////////////////
  //
  //  Input tree info
  //  copy from TreeBranches.C

  // Declaration of leaves types
  Int_t           T_Event_nTruePU;
  Int_t           T_Event_nPU;
  Int_t           T_Event_nPVoffline;
  Int_t           T_Event_HLT_Mu40;
  Int_t           T_Event_HLT_IsoMu24_IterTrk02;
  Int_t           T_Event_HLT_IsoTkMu24_IterTrk02;
  Int_t           T_Event_HLT_Mu17_Mu8;
  Int_t           T_Event_HLT_Mu17_TkMu8;
  Int_t           T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  Int_t           T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Int_t           T_Event_HLT_Mu17_NoFilters;
  Int_t           T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced;
  Int_t           T_Event_HLT_Mu30_TkMu11;
  Int_t           T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP;
  Int_t           T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP;
  Int_t           T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso;
  Int_t           T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  Int_t           T_Event_HLT_Ele17_Ele8_Gsf;
  Int_t           T_Event_HLT_Ele27_WP80_Gsf;
  vector<int>*     T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q;
  vector<int>*     T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02;
  vector<int>*     T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02;
  vector<int>*     T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8;
  vector<int>*     T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  vector<int>*     T_hltDiMuonGlb17Glb8DzFiltered0p2;
  vector<int>*     T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  vector<int>*     T_hltDiMuonGlbFiltered17TrkFiltered8;
  vector<int>*     T_hltDiMuonGlb17Trk8DzFiltered0p2;
  vector<int>*     T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4;
  vector<int>*     T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4;
  vector<int>*     T_hltL3NoFiltersfL1sMu12L3Filtered17;
  vector<int>*     T_hltL3fL1sMu25erL1f0L2f25L3Filtered30;
  vector<int>*     T_hltDiMuonGlbFiltered30TrkFiltered11;
  vector<int>*     T_hltDiMuonGlb30Trk11DzFiltered0p2;
  vector<int>*     T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8;
  vector<int>*     T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8;
  vector<int>*     T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter;
  vector<int>*     T_hltL1Mu12EG7L3MuFiltered23;
  vector<int>*     T_hltL1Mu12EG7L3IsoMuFiltered23;
  vector<int>*     T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter;
  vector<int>*     T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter;
  vector<int>*     T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter;
  vector<int>*     T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter;
  vector<int>*     T_hltEle17Ele8GsfTrackIsoLeg1Filter;
  vector<int>*     T_hltEle17Ele8GsfTrackIsoLeg2Filter;
  vector<int>*     T_hltEle27WP80GsfTrackIsoFilter;
  vector<float>*   T_mu_pt;
  vector<float>*   T_mu_eta;
  vector<float>*   T_mu_phi;
  vector<float>*   T_mu_iso;
  vector<float>*   T_ele_pt;
  vector<float>*   T_ele_eta;
  vector<float>*   T_ele_phi;

  // Set branch addresses.
  eventsTree->SetBranchAddress("T_Event_nTruePU",&T_Event_nTruePU);
  eventsTree->SetBranchAddress("T_Event_nPU",&T_Event_nPU);
  eventsTree->SetBranchAddress("T_Event_nPVoffline",&T_Event_nPVoffline);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu40",&T_Event_HLT_Mu40);
  eventsTree->SetBranchAddress("T_Event_HLT_IsoMu24_IterTrk02",&T_Event_HLT_IsoMu24_IterTrk02);
  eventsTree->SetBranchAddress("T_Event_HLT_IsoTkMu24_IterTrk02",&T_Event_HLT_IsoTkMu24_IterTrk02);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu17_Mu8",&T_Event_HLT_Mu17_Mu8);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu17_TkMu8",&T_Event_HLT_Mu17_TkMu8);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",&T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu17_NoFilters",&T_Event_HLT_Mu17_NoFilters);
  eventsTree->SetBranchAddress("T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced",&T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu30_TkMu11",&T_Event_HLT_Mu30_TkMu11);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP",&T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP);
  eventsTree->SetBranchAddress("T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP",&T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP);
  eventsTree->SetBranchAddress("T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso",&T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso);
  eventsTree->SetBranchAddress("T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",&T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
  eventsTree->SetBranchAddress("T_Event_HLT_Ele17_Ele8_Gsf",&T_Event_HLT_Ele17_Ele8_Gsf);
  eventsTree->SetBranchAddress("T_Event_HLT_Ele27_WP80_Gsf",&T_Event_HLT_Ele27_WP80_Gsf);
  eventsTree->SetBranchAddress("T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q",&T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q);
  eventsTree->SetBranchAddress("T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02",&T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02);
  eventsTree->SetBranchAddress("T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02",&T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02);
  eventsTree->SetBranchAddress("T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8",&T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8);
  eventsTree->SetBranchAddress("T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17",&T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17);
  eventsTree->SetBranchAddress("T_hltDiMuonGlb17Glb8DzFiltered0p2",&T_hltDiMuonGlb17Glb8DzFiltered0p2);
  eventsTree->SetBranchAddress("T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17",&T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17);
  eventsTree->SetBranchAddress("T_hltDiMuonGlbFiltered17TrkFiltered8",&T_hltDiMuonGlbFiltered17TrkFiltered8);
  eventsTree->SetBranchAddress("T_hltDiMuonGlb17Trk8DzFiltered0p2",&T_hltDiMuonGlb17Trk8DzFiltered0p2);
  eventsTree->SetBranchAddress("T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4",&T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4);
  eventsTree->SetBranchAddress("T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4",&T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4);
  eventsTree->SetBranchAddress("T_hltL3NoFiltersfL1sMu12L3Filtered17",&T_hltL3NoFiltersfL1sMu12L3Filtered17);
  eventsTree->SetBranchAddress("T_hltL3fL1sMu25erL1f0L2f25L3Filtered30",&T_hltL3fL1sMu25erL1f0L2f25L3Filtered30);
  eventsTree->SetBranchAddress("T_hltDiMuonGlbFiltered30TrkFiltered11",&T_hltDiMuonGlbFiltered30TrkFiltered11);
  eventsTree->SetBranchAddress("T_hltDiMuonGlb30Trk11DzFiltered0p2",&T_hltDiMuonGlb30Trk11DzFiltered0p2);
  eventsTree->SetBranchAddress("T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8",&T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8);
  eventsTree->SetBranchAddress("T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8",&T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8);
  eventsTree->SetBranchAddress("T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter",&T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter);
  eventsTree->SetBranchAddress("T_hltL1Mu12EG7L3MuFiltered23",&T_hltL1Mu12EG7L3MuFiltered23);
  eventsTree->SetBranchAddress("T_hltL1Mu12EG7L3IsoMuFiltered23",&T_hltL1Mu12EG7L3IsoMuFiltered23);
  eventsTree->SetBranchAddress("T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter",&T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter);
  eventsTree->SetBranchAddress("T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter",&T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter);
  eventsTree->SetBranchAddress("T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter",&T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter);
  eventsTree->SetBranchAddress("T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",&T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter);
  eventsTree->SetBranchAddress("T_hltEle17Ele8GsfTrackIsoLeg1Filter",&T_hltEle17Ele8GsfTrackIsoLeg1Filter);
  eventsTree->SetBranchAddress("T_hltEle17Ele8GsfTrackIsoLeg2Filter",&T_hltEle17Ele8GsfTrackIsoLeg2Filter);
  eventsTree->SetBranchAddress("T_hltEle27WP80GsfTrackIsoFilter",&T_hltEle27WP80GsfTrackIsoFilter);
  eventsTree->SetBranchAddress("T_mu_pt",&T_mu_pt);
  eventsTree->SetBranchAddress("T_mu_eta",&T_mu_eta);
  eventsTree->SetBranchAddress("T_mu_phi",&T_mu_phi);
  eventsTree->SetBranchAddress("T_mu_iso",&T_mu_iso);
  eventsTree->SetBranchAddress("T_ele_pt",&T_ele_pt);
  eventsTree->SetBranchAddress("T_ele_eta",&T_ele_eta);
  eventsTree->SetBranchAddress("T_ele_phi",&T_ele_phi);

  

  /////////////////////////////////
  //
  //   Define histograms
  //
  /////////////////////////////////
  
  int nbins = 30;//10, 20
  float xmin = 0;//10
  float xmax = 60;//70

  TH1F *h1         = new TH1F("h1","",nbins,xmin,xmax); //reco
  TH1F *h1_matched = new TH1F("h1_matched","",nbins,xmin,xmax); //reco matched to HLT

  TH1F *h2         = new TH1F("h2","",nbins,xmin,xmax); //reco
  TH1F *h2_matched = new TH1F("h2_matched","",nbins,xmin,xmax); //reco matched to HLT
  
  bool drawVector = false;
  bool drawScalar = true;
  
  int& drawWhatS = T_Event_nTruePU;
  vector<float>* drawWhatV = T_ele_eta;


  /////////////////////////////////
  //
  //   Loop on input tree
  //
  /////////////////////////////////
  
  for(int i=0; i<eventsTree->GetEntries(); i++){    
    eventsTree->GetEntry(i);

    ///////////////////////////////////////
    //    reco mu and ele selection
    ///////////////////////////////////////

    if(T_mu_pt->size() != 2 ) continue;
    // [0] is muon, [1] is electron
    
    if(fabs(T_mu_eta->at(0)) >= 2.4) continue;
    if(fabs(T_ele_eta->at(1)) >= 2.5) continue;

    if(T_mu_iso->at(0) > 0.2) continue;
   
    if ( !( (T_mu_pt->at(0) > 25 && T_ele_pt->at(1) > 15) ) ) continue;

        

    ////////////////////////////////
    // fill reco histograms
    ////////////////////////////////
    if(drawVector){
      h1->Fill(drawWhatV->at(1)); // 0 for mu, 1 for ele
    }
    if(drawScalar){
      h1->Fill(drawWhatS);
    }
    // histo2
    if(T_mu_pt->at(0) > 25){
      if(drawVector){
	h2->Fill(drawWhatV->at(1));
      }
      if(drawScalar){
	h2->Fill(drawWhatS);
      }
    }

    //////////////////////////////////////
    //    matching to HLT
    //////////////////////////////////////

   
    //bool HLTmatched = false;    
    bool Mu8LegMatched = false;
    bool Ele23LegMatched = false;
    bool Mu23LegMatched = false;
    bool Ele12LegMatched = false;
    bool SingleMu24Matched = false;
    
    
    // HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1
    if( (T_mu_pt->at(0) > 25 && T_ele_pt->at(1) > 25) ) {
      if(T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8->at(0) == 1) Mu8LegMatched = true;
      if(T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter->at(1) == 1 ) Ele23LegMatched = true;
    }
    
    // HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1
    if( T_mu_pt->at(0) > 25 && T_ele_pt->at(1) > 15){
      if(T_hltL1Mu12EG7L3IsoMuFiltered23->at(0) == 1) Mu23LegMatched = true;
      if(T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter->at(1) == 1 ) Ele12LegMatched = true;
    }
    

    // HLT_IsoTkMu24_IterTrk02_v1
    if(T_mu_pt->at(0) > 25){
      if(T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02->at(0) == 1) SingleMu24Matched = true;
    }


    //////////////////////////////////////////////////////
    // fill reco matched histograms
    //////////////////////////////////////////////////////
    if( Mu23LegMatched && Ele12LegMatched ){
      if(drawVector){
	h1_matched->Fill(drawWhatV->at(1));
      }
      if(drawScalar){
	h1_matched->Fill(drawWhatS);
      }
    }
    // histo2
    if( (Mu23LegMatched && Ele12LegMatched) || SingleMu24Matched ){
      if(drawVector){
	h2_matched->Fill(drawWhatV->at(1));
      }
      if(drawScalar){
	h2_matched->Fill(drawWhatS);
      }
    }


  } // end of loop on tree entries
  


  

  /////////////////////////
  //
  //     TGraphs
  //
  /////////////////////////
  
  TGraphAsymmErrors *gr1 = new TGraphAsymmErrors();
  gr1->Divide(h1_matched,h1);
  gr1->SetMarkerColor(4);
  gr1->SetLineColor(4);  
  gr1->SetMarkerStyle(20);


  TGraphAsymmErrors *gr2 = new TGraphAsymmErrors();
  gr2->Divide(h2_matched,h2);
  gr2->SetMarkerColor(2);
  gr2->SetLineColor(2);
  gr2->SetMarkerStyle(22);
  
  
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  c2->SetGrid(1,1);
  c2->cd();


  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr1,"p");
  mg->Add(gr2,"p");
  mg->Draw("a");

  //mg->GetXaxis()->SetTitle("electron p_{T}");
  //mg->GetXaxis()->SetTitle("electron #eta");
  mg->GetXaxis()->SetTitle("nTruePU");
  mg->GetYaxis()->SetTitle("efficiency");
  //mg->GetXaxis()->SetTitleSize(0.045);
  //mg->GetXaxis()->SetLabelSize(0.040);
  //mg->GetYaxis()->SetTitleSize(0.045);
  //mg->GetYaxis()->SetLabelSize(0.040);
  mg->GetYaxis()->SetRangeUser(0.5,1.05);
  mg->GetYaxis()->SetDecimals();
	
  gPad->Modified();
  
  //  TLegend * leg = new TLegend(0.6,0.75,0.93,0.91,"#splitline{SingleMu *fired*,}{DoubleMu *not fired*}"); 
  //TLegend * leg = new TLegend(0.6,0.75,0.93,0.91,"HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP"); 
  TLegend * leg = new TLegend(0.6,0.75,0.93,0.91); 
  //leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->AddEntry(gr1,"Mu8Ele23","p");
  leg->AddEntry(gr2,"Mu8Ele23 OR Mu24","p");
  leg->SetFillColor(0);
  leg->SetShadowColor(0);
  leg->SetLineColor(0);

  leg->Draw();



  

  return 0;
}
