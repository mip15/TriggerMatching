// -*- C++ -*-
//
// Package:    TriggerMatching/TriggerMatching
// Class:      TriggerMatching
// 
/**\class TriggerMatching TriggerMatching.cc TriggerMatching/TriggerMatching/plugins/TriggerMatching.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Mirena Ivova Rikova
//         Created:  Fri, 11 Jul 2014 18:55:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

/////

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TMath.h"
#include "TH2D.h"
#include "TH1D.h"

#include "TTree.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <vector>

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"


#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Common/interface/RefVector.h" // is it needed?

//
// class declaration
//

class TriggerMatching : public edm::EDAnalyzer {
public:
  explicit TriggerMatching(const edm::ParameterSet&);
  ~TriggerMatching();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  virtual void beginEvent();
  virtual void endEvent();
  
  virtual bool hasWasMother(const reco::GenParticle );


  // ----------member data ---------------------------
  HLTConfigProvider hltConfig;
  int triggerBit;
  std::vector<int> vec_triggerBit;
  std::vector<TString> vec_HLT_name;
  std::vector<std::string> vec_HLT_triggerObjects;

  // other
  //edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken_;
  //edm::EDGetTokenT<edm::ValueMap<bool> > electronIdToken_;

  // input tags
  edm::InputTag primaryVertexInputTag_;


  //////////////////////////////////
  //  Tree
  
  TTree* mytree_;

  int T_Event_nTruePU; //flat
  int T_Event_nPU;
  int T_Event_nPVoffline;

  //HLT path names  
  int T_Event_HLT_Mu40; // 0
  int T_Event_HLT_IsoMu24_IterTrk02; // 1
  int T_Event_HLT_IsoTkMu24_IterTrk02; // 2
  int T_Event_HLT_Mu17_Mu8; // 3
  int T_Event_HLT_Mu17_TkMu8; // 4
  int T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL; // 5
  int T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL; // 6
  int T_Event_HLT_Mu17_NoFilters; // 7
  int T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced; // 8

  int T_Event_HLT_Mu30_TkMu11; // 9
  int T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP; // 10
  int T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP; // 11
  int T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso; // 12 
  int T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL; // 13 
  int T_Event_HLT_Ele17_Ele8_Gsf; // 14 
  int T_Event_HLT_Ele27_WP80_Gsf; // 15 


  //modules for trigger matching
  std::vector<int> *T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q;
  std::vector<int> *T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02;
  std::vector<int> *T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02;
  std::vector<int> *T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8;
  std::vector<int> *T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  std::vector<int> *T_hltDiMuonGlb17Glb8DzFiltered0p2;
  std::vector<int> *T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  std::vector<int> *T_hltDiMuonGlbFiltered17TrkFiltered8;
  std::vector<int> *T_hltDiMuonGlb17Trk8DzFiltered0p2;
  //std::vector<int> *T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  //std::vector<int> *T_hltDiMuonGlb17Glb8DzFiltered0p2;
  //std::vector<int> *T_hltL3MuonRelTrkIsolationVVL; //nothing matches to this, this is a producer
  std::vector<int> *T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4;
  //std::vector<int> *T_hltDiMuonGlbFiltered17TrkFiltered8;
  //std::vector<int> *T_hltDiMuonGlb17Trk8DzFiltered0p2;
  //std::vector<int> *T_hltGlbTrkMuonRelTrkIsolationVVL; //nothing matches to this, this is a producer
  std::vector<int> *T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4;
  std::vector<int> *T_hltL3NoFiltersfL1sMu12L3Filtered17;

  //HLT_Mu30_TkMu11_v1,
  std::vector<int> *T_hltL3fL1sMu25erL1f0L2f25L3Filtered30;
  std::vector<int> *T_hltDiMuonGlbFiltered30TrkFiltered11; 
  std::vector<int> *T_hltDiMuonGlb30Trk11DzFiltered0p2;

  // HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1,
  std::vector<int> *T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8;
  std::vector<int> *T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8;  
  std::vector<int> *T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter;

  // HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1,
  std::vector<int> *T_hltL1Mu12EG7L3MuFiltered23;
  std::vector<int> *T_hltL1Mu12EG7L3IsoMuFiltered23;
  std::vector<int> *T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter;

  //HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1,
  std::vector<int> *T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter;
  std::vector<int> *T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter;

  // HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1,
  std::vector<int> *T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter;

  // HLT_Ele17_Ele8_Gsf_v1
  std::vector<int> *T_hltEle17Ele8GsfTrackIsoLeg1Filter;
  std::vector<int> *T_hltEle17Ele8GsfTrackIsoLeg2Filter;

  // HLT_Ele27_WP80_Gsf_v1
  std::vector<int> *T_hltEle27WP80GsfTrackIsoFilter;


  //reco mu properties
  std::vector<float> *T_mu_pt;
  std::vector<float> *T_mu_eta;
  std::vector<float> *T_mu_phi;
  std::vector<float> *T_mu_iso;

  //reco electron properties
  std::vector<float> *T_ele_pt;
  std::vector<float> *T_ele_eta;
  std::vector<float> *T_ele_phi;

  //will need a vector of vectors to loop:

  //std::map< int, std::vector<int>*  > test;

  //std::vector< std::vector<int> > matrix;

  //
  //////////////////////////////////



  // parameters from config file
  double muPtCut_;
  double muEtaCut_;
  double elePtCut_;
  double eleEtaCut_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerMatching::TriggerMatching(const edm::ParameterSet& iConfig)
  // :  full5x5SigmaIEtaIEtaMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap"))), 
  //electronIdToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronIDs"))),  
  //matrix(0)
{
  //now do what ever initialization is needed
 
  //take parameter values from the config
  muPtCut_= iConfig.getParameter<double>("muPtCut");
  muEtaCut_ = iConfig.getParameter<double>("muEtaCut");
  //isoCut_ = iConfig.getParameter<double>("isoCut");
  primaryVertexInputTag_  = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");

  elePtCut_ = iConfig.getParameter<double>("elePtCut");
  eleEtaCut_ = iConfig.getParameter<double>("eleEtaCut");

  // HLT paths of interest
  
  vec_HLT_name.push_back("HLT_Mu40_v"); // 0
  vec_HLT_name.push_back("HLT_IsoMu24_IterTrk02_v"); // 1
  vec_HLT_name.push_back("HLT_IsoTkMu24_IterTrk02_v"); // 2
  vec_HLT_name.push_back("HLT_Mu17_Mu8_v"); // 3
  vec_HLT_name.push_back("HLT_Mu17_TkMu8_v"); // 4
  vec_HLT_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"); // 5 
  vec_HLT_name.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"); // 6
  vec_HLT_name.push_back("HLT_Mu17_NoFilters_v"); // 7
  vec_HLT_name.push_back("HLT_DoubleMu4NoFilters_Jpsi_Displaced_v"); // 8
  vec_HLT_name.push_back("HLT_Mu30_TkMu11_v"); // 9
  vec_HLT_name.push_back("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v"); // 10
  vec_HLT_name.push_back("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v"); // 11
  vec_HLT_name.push_back("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v"); // 12
  vec_HLT_name.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"); // 13
  vec_HLT_name.push_back("HLT_Ele17_Ele8_Gsf_v"); // 14
  vec_HLT_name.push_back("HLT_Ele27_WP80_Gsf_v"); // 15



  // names of the filters in the HLT paths to be matched

  // HLT_Mu40_v // 0
  vec_HLT_triggerObjects.push_back("hltL3fL1sMu16L1f0L2f16QL3Filtered40Q"); // 0

  // HLT_IsoMu24_IterTrk02_v // 1
  vec_HLT_triggerObjects.push_back("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02");  // 1

  // HLT_IsoTkMu24_IterTrk02_v // 2
  vec_HLT_triggerObjects.push_back("hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02");  // 2

  // ("HLT_Mu17_Mu8_v"); // 3
  vec_HLT_triggerObjects.push_back("hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8"); // 3 
  vec_HLT_triggerObjects.push_back("hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17"); // 4 
  vec_HLT_triggerObjects.push_back("hltDiMuonGlb17Glb8DzFiltered0p2"); // 5

  // ("HLT_Mu17_TkMu8_v"); // 4
  vec_HLT_triggerObjects.push_back("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17"); // 6 
  vec_HLT_triggerObjects.push_back("hltDiMuonGlbFiltered17TrkFiltered8"); // 7
  vec_HLT_triggerObjects.push_back("hltDiMuonGlb17Trk8DzFiltered0p2"); // 8
  
  // ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"); // 5 
  vec_HLT_triggerObjects.push_back("hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4"); // 9

  // ("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"); // 6
  vec_HLT_triggerObjects.push_back("hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4"); // 10 
  
  // ("HLT_Mu17_NoFilters_v"); // 7
  vec_HLT_triggerObjects.push_back("hltL3NoFiltersfL1sMu12L3Filtered17"); // 11

  // ("HLT_Mu30_TkMu11_v") // 9
  vec_HLT_triggerObjects.push_back("hltL3fL1sMu25erL1f0L2f25L3Filtered30"); // 12
  vec_HLT_triggerObjects.push_back("hltDiMuonGlbFiltered30TrkFiltered11"); // 13
  vec_HLT_triggerObjects.push_back("hltDiMuonGlb30Trk11DzFiltered0p2"); // 14

  // HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP; // 10
  vec_HLT_triggerObjects.push_back("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8"); // 15
  vec_HLT_triggerObjects.push_back("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8"); // 16
  vec_HLT_triggerObjects.push_back("hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter"); // 17
  
  // HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1, // 11
  vec_HLT_triggerObjects.push_back("hltL1Mu12EG7L3MuFiltered23"); // 18
  vec_HLT_triggerObjects.push_back("hltL1Mu12EG7L3IsoMuFiltered23"); // 19
  vec_HLT_triggerObjects.push_back("hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter"); // 20
  
  // HLT_Ele23_Ele12_CaloId_TrackId_Iso; // 12 
  vec_HLT_triggerObjects.push_back("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter"); // 21
  vec_HLT_triggerObjects.push_back("hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter"); // 22
  
  // HLT_DoubleEle33_CaloIdL_GsfTrkIdVL; // 13 
  vec_HLT_triggerObjects.push_back("hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter"); // 23
  
  // HLT_Ele17_Ele8_Gsf; // 14 
  vec_HLT_triggerObjects.push_back("hltEle17Ele8GsfTrackIsoLeg1Filter"); // 24
  vec_HLT_triggerObjects.push_back("hltEle17Ele8GsfTrackIsoLeg2Filter"); // 25
 
  // HLT_Ele27_WP80_Gsf; // 15  
  vec_HLT_triggerObjects.push_back("hltEle27WP80GsfTrackIsoFilter"); // 26

}


TriggerMatching::~TriggerMatching()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerMatching::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;


  beginEvent();


  /////////////////////////                                       
  // 
  //   Pile up                                            
  // 
  /////////////////////////                
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;

  T_Event_nTruePU=-1.;
  T_Event_nPU=-1.;

  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    int BX = PVI->getBunchCrossing();
    if(BX == 0) {
      T_Event_nTruePU = PVI->getTrueNumInteractions();
      T_Event_nPU = PVI->getPU_NumInteractions();
    }
  }


  /////////////////////////////////////////////////////////
  //
  //
  
  // HLT paths names
  bool changedConfig = false;
  if (!hltConfig.init(iEvent.getRun(), iSetup, "TEST", changedConfig)) {
    //cout << "Initialization of HLTConfigProvider failed!!" << endl;
    return;
  }
  
  if (changedConfig){
    //std::cout << "the curent menu is " << hltConfig.tableName() << std::endl;
    triggerBit = -1;
    for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
      if (TString(hltConfig.triggerNames()[j]).Contains("HLT_Mu40_v")) triggerBit = j;
    }
    //if (triggerBit == -1) cout << "HLT path not found" << endl;
       
  }
  //cout<<"triggerBit = "<<triggerBit<<endl;


  // the same but in vector mode !!!
  changedConfig = false;
  if (!hltConfig.init(iEvent.getRun(), iSetup, "TEST", changedConfig)) {
    //cout << "Initialization of HLTConfigProvider failed!!" << endl;
    return;
  }
  
  if (!changedConfig){ //should be == true, must remove the ! in the end !!!
    //std::cout << "the curent menu is " << hltConfig.tableName() << std::endl;
    vec_triggerBit.clear();
    //clear module labels here !!!
    for(size_t m = 0; m < vec_HLT_name.size(); m++){
      bool foundPathInMenu = false;
      for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
	if (TString(hltConfig.triggerNames()[j]).Contains(vec_HLT_name[m])) {
	  vec_triggerBit.push_back(j);
	  foundPathInMenu = true;
	}
      }
      if (foundPathInMenu == false) {
	//cout << "HLT path not found" << endl;
	vec_triggerBit.push_back(-1);
      }
    }
  }

  //test if it's correct !!! ********
  /* 
     cout<<"vectriggerbit size = "<<vec_triggerBit.size()<<endl;
     for(size_t k = 0; k<vec_triggerBit.size(); k++){
     cout<<"trigger bit "<<k<<" = "<<vec_triggerBit[k]<<endl;
     }
  */
  // end of test !!! ********

  //
  //
  ////////////////////////////////////////////////////////////





  //open the trigger summary
  edm::InputTag triggerSummaryLabel_ = edm::InputTag("hltTriggerSummaryAOD", "", "TEST");
  edm::Handle<trigger::TriggerEvent> triggerSummary;
  iEvent.getByLabel(triggerSummaryLabel_, triggerSummary);

  //trigger results
  edm::InputTag triggerResultsLabel = edm::InputTag("TriggerResults", "", "TEST");
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(triggerResultsLabel, triggerResults);
  
  //trigger object we want to match
  //edm::InputTag filterTag = edm::InputTag("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8", "", "HLT");
  edm::InputTag filterTag = edm::InputTag("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02", "", "TEST"); 

  //now needs a vector of these !!!
  std::vector<edm::InputTag> vec_filterTag;

  for(size_t g = 0; g < vec_HLT_triggerObjects.size(); g++){
    vec_filterTag.push_back(edm::InputTag(vec_HLT_triggerObjects[g], "", "TEST") );
    
    
  }


  //find the index corresponding to the event
  size_t filterIndex = (*triggerSummary).filterIndex(filterTag);
  //now needs a vector of these !!!
  std::vector<size_t> vec_filterIndex; 
  //for(size_t j=0; j<vec_filterTag.size(); j++) { vec_filterIndex[j] = (*triggerSummary).filterIndex(vec_filterTag[j]); }
  for(size_t j=0; j<vec_filterTag.size(); j++) { vec_filterIndex.push_back((*triggerSummary).filterIndex(vec_filterTag[j])); }


  

  /////////////////////////////////////////////
  //
  //

  trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
  trigger::TriggerObjectCollection selectedObjects;
  if (filterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
    for (size_t j = 0; j < keys.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
      selectedObjects.push_back(foundObject);
    }
  }


  //now the above in vector mode !!!
  vector<int> theHLTcorr;
  trigger::TriggerObjectCollection selectedObjects_tmp;
  for(size_t t=0; t<vec_filterTag.size(); t++){
    if(vec_filterIndex[t] < (*triggerSummary).sizeFilters()){ //check if the trigger object is present
      const trigger::Keys &keys = (*triggerSummary).filterKeys(vec_filterIndex[t]);
      for (size_t j = 0; j < keys.size(); j++) {
	trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	selectedObjects_tmp.push_back(foundObject);
	theHLTcorr.push_back(t);
      }
    }
  }

  //
  //
  //////////////////////////////////////////////




  // Fill the trigger bits
	
  T_Event_HLT_Mu40 = (vec_triggerBit[0]==-1) ? -1 : triggerResults->accept(vec_triggerBit[0]);
  T_Event_HLT_IsoMu24_IterTrk02 = (vec_triggerBit[1]==-1) ? -1 : triggerResults->accept(vec_triggerBit[1]);
  T_Event_HLT_IsoTkMu24_IterTrk02 = (vec_triggerBit[2]==-1) ? -1 : triggerResults->accept(vec_triggerBit[2]);
  T_Event_HLT_Mu17_Mu8  = (vec_triggerBit[3]==-1) ? -1 : triggerResults->accept(vec_triggerBit[3]);
  T_Event_HLT_Mu17_TkMu8  = (vec_triggerBit[4]==-1) ? -1 : triggerResults->accept(vec_triggerBit[4]);
  T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL      = (vec_triggerBit[5]==-1) ? -1 : triggerResults->accept(vec_triggerBit[5]);
  T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL    = (vec_triggerBit[6]==-1) ? -1 : triggerResults->accept(vec_triggerBit[6]);
  T_Event_HLT_Mu17_NoFilters                    = (vec_triggerBit[7]==-1) ? -1 : triggerResults->accept(vec_triggerBit[7]);
  T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced = (vec_triggerBit[8]==-1) ? -1 : triggerResults->accept(vec_triggerBit[8]);
	
  T_Event_HLT_Mu30_TkMu11 = (vec_triggerBit[9]==-1) ? -1 : triggerResults->accept(vec_triggerBit[9]);
  T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP = (vec_triggerBit[10]==-1) ? -1 : triggerResults->accept(vec_triggerBit[10]);
  T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP = (vec_triggerBit[11]==-1) ? -1 : triggerResults->accept(vec_triggerBit[11]);
  T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso = (vec_triggerBit[12]==-1) ? -1 : triggerResults->accept(vec_triggerBit[12]);
  T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = (vec_triggerBit[13]==-1) ? -1 : triggerResults->accept(vec_triggerBit[13]);
  T_Event_HLT_Ele17_Ele8_Gsf = (vec_triggerBit[14]==-1) ? -1 : triggerResults->accept(vec_triggerBit[14]);
  T_Event_HLT_Ele27_WP80_Gsf= (vec_triggerBit[15]==-1) ? -1 : triggerResults->accept(vec_triggerBit[15]);


 
  edm::Handle< std::vector<reco::Vertex> > vertices_h;
  iEvent.getByLabel(primaryVertexInputTag_, vertices_h);  
  T_Event_nPVoffline = vertices_h->size();

  edm::Handle< reco::MuonCollection > muons;
  iEvent.getByLabel( "muons", muons );
  
  //edm::Handle<reco::GsfElectronCollection> electrons;
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByLabel("gedGsfElectrons",electrons);


  /////////////////////////////////////////////////////////
  //
  //    Loop on reco muons
  //
  /////////////////////////////////////////////////////////
   
  for( size_t iMuon = 0; iMuon < muons->size(); ++iMuon ) {
    const reco::Muon & recoMuon = (*muons)[iMuon];

    //////////////////////////////////////////////////////
    //     acceptance cuts 
    /////////////////////////////////////////////////////
    if ((recoMuon.pt()) < muPtCut_) continue;
    if (fabs(recoMuon.eta()) > muEtaCut_) continue;

    //////////////////////////////////////////////////////////////////////////////////////
    //    check if the muon is tight (w.r.t. any PV)
    /////////////////////////////////////////////////////////////////////////////////////

    std::vector<reco::Vertex>::const_iterator itv;
    reco::VertexRef vtx(vertices_h, 0);    
    
    if( !(muon::isTightMuon(recoMuon, *vtx)) ) continue;
    
    /////////////////////////////////////////////////////////
    //    check the muon isolation
    /////////////////////////////////////////////////////////

    float isolation = (  recoMuon.pfIsolationR04().sumChargedHadronPt+ max(0.,recoMuon.pfIsolationR04().sumNeutralHadronEt+recoMuon.pfIsolationR04().sumPhotonEt - 0.5*recoMuon.pfIsolationR04().sumPUPt)  ) / recoMuon.pt();

    //if(isolation>100) continue;
    //if(isolation > isoCut_) continue;



    //////////////////////////////////////////////////////////////
    //     check if muon is MC matched
    ///////////////////////////////////////////////////////////////

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel( "genParticles", genParticles );

 
    int theNbOfGen = genParticles->size();
    bool isMCmatched = false;
    for (int i=0 ; i < theNbOfGen; i++){
      const reco::GenParticle & genMuon = (*genParticles)[i];
      if ((fabs(genMuon.pdgId())==13)&&(genMuon.status()==1)&&(hasWasMother(genMuon))){
	float deltaR = sqrt(pow(genMuon.eta() - recoMuon.eta(),2)+ pow(acos(cos( genMuon.eta() - recoMuon.eta()  )),2)) ;	
	if(deltaR > 0.1) {continue;}
	else{isMCmatched = true;}	
      }
    }
    
    if(isMCmatched==false) continue;

    //////////////////////////////////////////////////////////////
    //     preliminary selection passed,
    //     now do the trigger matching
    //////////////////////////////////////////////////////////////

    vector<int> array_aux(50, -1);
    for(size_t mm = 0; mm <selectedObjects_tmp.size(); mm++){
      float deltaR = sqrt(pow(selectedObjects_tmp[mm].eta()-recoMuon.eta(),2)+pow(acos(cos(selectedObjects_tmp[mm].phi()-recoMuon.phi())),2));
      if(deltaR < 0.1){ //matching to trigger module successful
	array_aux[theHLTcorr[mm]] = 1; 

      }


    } //end of loop on the selected trigger objects
    
    

    int kk = 0;
    //modules for trigger matching
    T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q->push_back(array_aux[kk]); kk++;
    T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02->push_back(array_aux[kk]); kk++;
    T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02->push_back(array_aux[kk]); kk++;
    T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8->push_back(array_aux[kk]); kk++;
    T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Glb8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlbFiltered17TrkFiltered8->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Trk8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    // *T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17->push_back(array_aux[kk]); kk++;
    // *T_hltDiMuonGlb17Glb8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    //T_hltL3MuonRelTrkIsolationVVL->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4->push_back(array_aux[kk]); kk++;
    // *T_hltDiMuonGlbFiltered17TrkFiltered8->push_back(array_aux[kk]); kk++;
    // *T_hltDiMuonGlb17Trk8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    //T_hltGlbTrkMuonRelTrkIsolationVVL->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4->push_back(array_aux[kk]); kk++;
    T_hltL3NoFiltersfL1sMu12L3Filtered17->push_back(array_aux[kk]); kk++;
    T_hltL3fL1sMu25erL1f0L2f25L3Filtered30->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlbFiltered30TrkFiltered11->push_back(array_aux[kk]); kk++; 
    T_hltDiMuonGlb30Trk11DzFiltered0p2->push_back(array_aux[kk]); kk++;
    T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8->push_back(array_aux[kk]); kk++;
    T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8->push_back(array_aux[kk]); kk++;  
    T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter->push_back(array_aux[kk]); kk++;
    T_hltL1Mu12EG7L3MuFiltered23->push_back(array_aux[kk]); kk++;
    T_hltL1Mu12EG7L3IsoMuFiltered23->push_back(array_aux[kk]); kk++;
    T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter->push_back(array_aux[kk]); kk++;
    T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter->push_back(array_aux[kk]); kk++;
    T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter->push_back(array_aux[kk]); kk++;
    T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter->push_back(array_aux[kk]); kk++;
    T_hltEle17Ele8GsfTrackIsoLeg1Filter->push_back(array_aux[kk]); kk++;
    T_hltEle17Ele8GsfTrackIsoLeg2Filter->push_back(array_aux[kk]); kk++;
    T_hltEle27WP80GsfTrackIsoFilter->push_back(array_aux[kk]); kk++;



    //reco muons properties
    T_mu_pt->push_back(recoMuon.pt());
    T_mu_eta->push_back(recoMuon.eta());
    T_mu_phi->push_back(recoMuon.phi());
    T_mu_iso->push_back(isolation);


    //reco electrons  properties
    int dummy = -100;
    T_ele_pt->push_back(dummy);
    T_ele_eta->push_back(dummy);
    T_ele_phi->push_back(dummy);

  } //end of loop on reco muons




  /////////////////////////////////////////////////////////
  //
  //    Loop on reco electrons
  //
  /////////////////////////////////////////////////////////
 
  for( size_t iEle = 0; iEle < electrons->size(); ++iEle ) {
    const reco::GsfElectron & recoElectron = (*electrons)[iEle];

    //////////////////////////////////////////////////////
    //     acceptance cuts 
    /////////////////////////////////////////////////////
    if(recoElectron.pt() < elePtCut_) continue;
    if(fabs(recoElectron.eta()) > eleEtaCut_) continue;

    ///////////////////////////////////////////////////
    //    ID requirements
    ///////////////////////////////////////////////////

    reco::VertexRef vtx(vertices_h, 0);    

    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel("offlineBeamSpot", bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();

    edm::Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel("allConversions", hConversions);


    // ID and matching
    float dEtaIn_ = recoElectron.deltaEtaSuperClusterTrackAtVtx();
    float dPhiIn_ = recoElectron.deltaPhiSuperClusterTrackAtVtx();
    float hOverE_ = recoElectron.hcalOverEcal();
    float sigmaIetaIeta_ = recoElectron.sigmaIetaIeta();
    float ooEmooP_ = fabs(1.0/recoElectron.ecalEnergy() - recoElectron.eSuperClusterOverP()/recoElectron.ecalEnergy() );

    // PF isolation
    reco::GsfElectron::PflowIsolationVariables pfIso = recoElectron.pfIsolationVariables();
    float absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    float relIsoWithDBeta_ = absiso/recoElectron.pt();

    // Impact parameter
    float d0_ = (-1) * recoElectron.gsfTrack()->dxy(vtx->position() );
    float dz_ = recoElectron.gsfTrack()->dz( vtx->position() );

    // Conversion rejection
    float expectedMissingInnerHits_ = recoElectron.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    bool passConversionVeto_ = !ConversionTools::hasMatchedConversion(recoElectron,hConversions, beamspot.position());
    
    
    //edm::Handle<edm::ValueMap<float> > full5x5sieie;
    //iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_,full5x5sieie);
    //float full5x5_sigmaIetaIeta_;
    //full5x5_sigmaIetaIeta_ = (*full5x5sieie)[recoElectron.refVector()];
    


    // electron ID cuts

    // barrel electron, "tight" selection
    if(fabs(recoElectron.superCluster()->eta() ) <= 1.479 ){
      if( fabs(dEtaIn_) >= 0.0091 ) continue;
      if( fabs(dPhiIn_) >= 0.031 ) continue;
      //if( full5x5_sigmaIetaIeta_ >= 0.0106 ) continue; // must use this value
      if( sigmaIetaIeta_ >= 0.01 ) continue; // obsolete, should not use this variable, to be fixed!
      if( hOverE_ >= 0.0532 ) continue;
      if( fabs(d0_) >= 0.0126 ) continue;
      if( fabs(dz_) >= 0.0116 ) continue;
      if( fabs(ooEmooP_) >= 0.0609 ) continue;
      if( relIsoWithDBeta_ >= 0.1649 ) continue;
      if( !passConversionVeto_ ) continue;
      if( expectedMissingInnerHits_ >= 1) continue;  
    }

    // endcap electron, "tight" selection
    if( fabs(recoElectron.superCluster()->eta() ) > 1.479 &&  fabs(recoElectron.superCluster()->eta() ) < 2.5 ){
      if( fabs(dEtaIn_) >= 0.0106 ) continue;
      if( fabs(dPhiIn_) >= 0.0359 ) continue;
      //if( full5x5_sigmaIetaIeta_ >= 0.0305 ) continue; // must use this value
      if( sigmaIetaIeta_ >= 0.03 ) continue; // obsolete, should not use this variable, to be fixed!
      if( hOverE_ >= 0.0835 ) continue;
      if( fabs(d0_) >= 0.0163 ) continue;
      if( fabs(dz_) >= 0.5999 ) continue;
      if( fabs(ooEmooP_) >= 0.1126 ) continue;
      if( relIsoWithDBeta_ >= 0.2075 ) continue;
      if( !passConversionVeto_ ) continue;
      if( expectedMissingInnerHits_ >= 1 ) continue;  
    }



    //cout<<"dEtaIn:\t"<<dEtaIn_<<"\tdPhiIn:\t"<<dPhiIn_<<"\tsigmaIetaIeta:\t"<< sigmaIetaIeta_<<"\thOverE:\t"<< hOverE_<<"\td0:\t"<< d0_<<"\tdz:\t"<<dz_<<"\tooEmooP:\t"<< ooEmooP_ <<"\trelIsoWithDBeta:\t"<<relIsoWithDBeta_<<"\tpassConversionVeto:\t"<< passConversionVeto_<<"\texpectedMissingInnerHits:\t"<< expectedMissingInnerHits_<<endl;


    //////////////////////////////////////////////////////////////
    //     check if electron  is MC matched
    ///////////////////////////////////////////////////////////////

    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel( "genParticles", genParticles );
    int theNbOfGen = genParticles->size();
    bool isMCmatched = false;
    for (int i=0 ; i < theNbOfGen; i++){
      const reco::GenParticle & genElectron = (*genParticles)[i];
      if ((fabs(genElectron.pdgId())==11)&&(genElectron.status()==1)&&(hasWasMother(genElectron))){
	float deltaR = sqrt(pow(genElectron.eta() - recoElectron.eta(),2)+ pow(acos(cos( genElectron.eta() - recoElectron.eta()  )),2)) ;	
	if(deltaR > 0.1) {continue;}
	else{isMCmatched = true;}	
      }
    }
    

    if(isMCmatched==false) continue;


    // the same as for muons:

    //////////////////////////////////////////////////////////////
    //     preliminary selection passed,
    //     now do the trigger matching
    //////////////////////////////////////////////////////////////

    vector<int> array_aux(50, -1);
    for(size_t mm = 0; mm <selectedObjects_tmp.size(); mm++){
      float deltaR = sqrt(pow(selectedObjects_tmp[mm].eta()-recoElectron.eta(),2)+pow(acos(cos(selectedObjects_tmp[mm].phi()-recoElectron.phi())),2));
      if(deltaR < 0.1){ //matching to trigger module successful
	array_aux[theHLTcorr[mm]] = 1; 
      }


    } //end of loop on the selected trigger objects
    

    int kk = 0;
    //modules for trigger matching
    T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q->push_back(array_aux[kk]); kk++;
    T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02->push_back(array_aux[kk]); kk++;
    T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02->push_back(array_aux[kk]); kk++;
    T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8->push_back(array_aux[kk]); kk++;
    T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Glb8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlbFiltered17TrkFiltered8->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Trk8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    // *T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17->push_back(array_aux[kk]); kk++;
    // *T_hltDiMuonGlb17Glb8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    //T_hltL3MuonRelTrkIsolationVVL->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4->push_back(array_aux[kk]); kk++;
    // *T_hltDiMuonGlbFiltered17TrkFiltered8->push_back(array_aux[kk]); kk++;
    // *T_hltDiMuonGlb17Trk8DzFiltered0p2->push_back(array_aux[kk]); kk++;
    //T_hltGlbTrkMuonRelTrkIsolationVVL->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4->push_back(array_aux[kk]); kk++;
    T_hltL3NoFiltersfL1sMu12L3Filtered17->push_back(array_aux[kk]); kk++;
    T_hltL3fL1sMu25erL1f0L2f25L3Filtered30->push_back(array_aux[kk]); kk++;
    T_hltDiMuonGlbFiltered30TrkFiltered11->push_back(array_aux[kk]); kk++; 
    T_hltDiMuonGlb30Trk11DzFiltered0p2->push_back(array_aux[kk]); kk++;
    T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8->push_back(array_aux[kk]); kk++;
    T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8->push_back(array_aux[kk]); kk++;  
    T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter->push_back(array_aux[kk]); kk++;
    T_hltL1Mu12EG7L3MuFiltered23->push_back(array_aux[kk]); kk++;
    T_hltL1Mu12EG7L3IsoMuFiltered23->push_back(array_aux[kk]); kk++;
    T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter->push_back(array_aux[kk]); kk++;
    T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter->push_back(array_aux[kk]); kk++;
    T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter->push_back(array_aux[kk]); kk++;
    T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter->push_back(array_aux[kk]); kk++;
    T_hltEle17Ele8GsfTrackIsoLeg1Filter->push_back(array_aux[kk]); kk++;
    T_hltEle17Ele8GsfTrackIsoLeg2Filter->push_back(array_aux[kk]); kk++;
    T_hltEle27WP80GsfTrackIsoFilter->push_back(array_aux[kk]); kk++;



    //reco muons properties
    int dummy = -100;
    T_mu_pt->push_back(dummy);
    T_mu_eta->push_back(dummy);
    T_mu_phi->push_back(dummy);
    T_mu_iso->push_back(dummy);


    //reco electrons  properties
    T_ele_pt->push_back(recoElectron.pt());
    T_ele_eta->push_back(recoElectron.eta());
    T_ele_phi->push_back(recoElectron.phi());


  }// enf of loop on reco electrons

 






  mytree_->Fill();
  
  endEvent();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerMatching::beginJob()
{

  edm::Service< TFileService > fileService;

  /////////////////////////////////////
  // Tree
 
  mytree_ = fileService->make<TTree>("eventsTree","");
  
  mytree_->Branch("T_Event_nTruePU", &T_Event_nTruePU, "T_Event_nTruePU/I");
  mytree_->Branch("T_Event_nPU", &T_Event_nPU, "T_Event_nPU/I");
  mytree_->Branch("T_Event_nPVoffline", &T_Event_nPVoffline, "T_Event_nPVoffline/I");

  //HLT path names
  mytree_->Branch("T_Event_HLT_Mu40", &T_Event_HLT_Mu40, "T_Event_HLT_Mu40/I");
  mytree_->Branch("T_Event_HLT_IsoMu24_IterTrk02", &T_Event_HLT_IsoMu24_IterTrk02, "T_Event_HLT_IsoMu24_IterTrk02/I");
  mytree_->Branch("T_Event_HLT_IsoTkMu24_IterTrk02", &T_Event_HLT_IsoTkMu24_IterTrk02, "T_Event_HLT_IsoTkMu24_IterTrk02/I");
  mytree_->Branch("T_Event_HLT_Mu17_Mu8", &T_Event_HLT_Mu17_Mu8, "T_Event_HLT_Mu17_Mu8/I");
  mytree_->Branch("T_Event_HLT_Mu17_TkMu8", &T_Event_HLT_Mu17_TkMu8, "T_Event_HLT_Mu17_TkMu8/I");
  mytree_->Branch("T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, "T_Event_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL/I");
  mytree_->Branch("T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, "T_Event_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL/I");
  mytree_->Branch("T_Event_HLT_Mu17_NoFilters", &T_Event_HLT_Mu17_NoFilters, "T_Event_HLT_Mu17_NoFilters/I");
  mytree_->Branch("T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced", &T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced, "T_Event_HLT_DoubleMu4NoFilters_Jpsi_Displaced/I");

  mytree_->Branch("T_Event_HLT_Mu30_TkMu11", &T_Event_HLT_Mu30_TkMu11, "T_Event_HLT_Mu30_TkMu11/I");
  mytree_->Branch("T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP", &T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP, "T_Event_HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP/I");
  mytree_->Branch("T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP", &T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP, "T_Event_HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP/I");
  mytree_->Branch("T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso", &T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso, "T_Event_HLT_Ele23_Ele12_CaloId_TrackId_Iso/I");
  mytree_->Branch("T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, "T_Event_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL/I");
  mytree_->Branch("T_Event_HLT_Ele17_Ele8_Gsf", &T_Event_HLT_Ele17_Ele8_Gsf, "T_Event_HLT_Ele17_Ele8_Gsf/I");
  mytree_->Branch("T_Event_HLT_Ele27_WP80_Gsf", &T_Event_HLT_Ele27_WP80_Gsf, "T_Event_HLT_Ele27_WP80_Gsf/I");


 
  //modules for trigger mathing
  //  mytree_->Branch("T_Gen_Elec_Px", "std::vector<float>", &T_Gen_Elec_Px);
  
 
  
  mytree_->Branch("T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q",				"std::vector<int>",&T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q);
  mytree_->Branch("T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02",	"std::vector<int>",&T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02);
  mytree_->Branch("T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02",	"std::vector<int>",&T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02);
  mytree_->Branch("T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8",		"std::vector<int>",&T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8);
  mytree_->Branch("T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17",		"std::vector<int>",&T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17);
  mytree_->Branch("T_hltDiMuonGlb17Glb8DzFiltered0p2",					"std::vector<int>",&T_hltDiMuonGlb17Glb8DzFiltered0p2);
  mytree_->Branch("T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17",			"std::vector<int>",&T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17);
  mytree_->Branch("T_hltDiMuonGlbFiltered17TrkFiltered8",				"std::vector<int>",&T_hltDiMuonGlbFiltered17TrkFiltered8);
  mytree_->Branch("T_hltDiMuonGlb17Trk8DzFiltered0p2",					"std::vector<int>",&T_hltDiMuonGlb17Trk8DzFiltered0p2);
  //mytree_->Branch("T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17",		"std::vector<int>",&T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17);
  //mytree_->Branch("T_hltDiMuonGlb17Glb8DzFiltered0p2",				"std::vector<int>",&T_hltDiMuonGlb17Glb8DzFiltered0p2);
  //mytree_->Branch("T_hltL3MuonRelTrkIsolationVVL",					"std::vector<int>",&T_hltL3MuonRelTrkIsolationVVL);
  mytree_->Branch("T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4",		"std::vector<int>",&T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4);
  //mytree_->Branch("T_hltDiMuonGlbFiltered17TrkFiltered8",				"std::vector<int>",&T_hltDiMuonGlbFiltered17TrkFiltered8);
  //mytree_->Branch("T_hltDiMuonGlb17Trk8DzFiltered0p2",				"std::vector<int>",&T_hltDiMuonGlb17Trk8DzFiltered0p2);
  //mytree_->Branch("T_hltGlbTrkMuonRelTrkIsolationVVL",					"std::vector<int>",&T_hltGlbTrkMuonRelTrkIsolationVVL);
  mytree_->Branch("T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4",		"std::vector<int>",&T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4);
  mytree_->Branch("T_hltL3NoFiltersfL1sMu12L3Filtered17",				"std::vector<int>",&T_hltL3NoFiltersfL1sMu12L3Filtered17);


  mytree_->Branch("T_hltL3fL1sMu25erL1f0L2f25L3Filtered30",				"std::vector<int>",&T_hltL3fL1sMu25erL1f0L2f25L3Filtered30);
  mytree_->Branch("T_hltDiMuonGlbFiltered30TrkFiltered11",				"std::vector<int>",&T_hltDiMuonGlbFiltered30TrkFiltered11);
  mytree_->Branch("T_hltDiMuonGlb30Trk11DzFiltered0p2",				"std::vector<int>",&T_hltDiMuonGlb30Trk11DzFiltered0p2);
  mytree_->Branch("T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8",				"std::vector<int>",&T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8);
  mytree_->Branch("T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8",				"std::vector<int>",&T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8);
  mytree_->Branch("T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter",				"std::vector<int>",&T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter);
  mytree_->Branch("T_hltL1Mu12EG7L3MuFiltered23",				"std::vector<int>",&T_hltL1Mu12EG7L3MuFiltered23);
  mytree_->Branch("T_hltL1Mu12EG7L3IsoMuFiltered23",				"std::vector<int>",&T_hltL1Mu12EG7L3IsoMuFiltered23);
  mytree_->Branch("T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter",				"std::vector<int>",&T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter);
  mytree_->Branch("T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter",				"std::vector<int>",&T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter);
  mytree_->Branch("T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter",				"std::vector<int>",&T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter);
  mytree_->Branch("T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",				"std::vector<int>",&T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter);
  mytree_->Branch("T_hltEle17Ele8GsfTrackIsoLeg1Filter",				"std::vector<int>",&T_hltEle17Ele8GsfTrackIsoLeg1Filter);
  mytree_->Branch("T_hltEle17Ele8GsfTrackIsoLeg2Filter",				"std::vector<int>",&T_hltEle17Ele8GsfTrackIsoLeg2Filter);
  mytree_->Branch("T_hltEle27WP80GsfTrackIsoFilter",				"std::vector<int>",&T_hltEle27WP80GsfTrackIsoFilter);
 



  //reco muon properties
  mytree_->Branch("T_mu_pt",	 "std::vector<float>",	 &T_mu_pt);
  mytree_->Branch("T_mu_eta",	 "std::vector<float>",	 &T_mu_eta);
  mytree_->Branch("T_mu_phi",	 "std::vector<float>",	 &T_mu_phi);
  mytree_->Branch("T_mu_iso",    "std::vector<float>",   &T_mu_iso);

  //reco electron properties
  mytree_->Branch("T_ele_pt",	 "std::vector<float>",	 &T_ele_pt);
  mytree_->Branch("T_ele_eta",	 "std::vector<float>",	 &T_ele_eta);
  mytree_->Branch("T_ele_phi",	 "std::vector<float>",	 &T_ele_phi);

  //
  ////////////////////////////////
}

// ------------ METHOD called once each job just after ending the event loop  ------------
void 
TriggerMatching::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  TriggerMatching::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  TriggerMatching::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  TriggerMatching::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  TriggerMatching::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerMatching::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


bool
TriggerMatching::hasWasMother(const reco::GenParticle  p)
{

  bool foundH = false;
  if (p.numberOfMothers()==0) return foundH;
  const reco::Candidate  *part = (p.mother());
  // loop on the mother particles to check if is has a W has mother
  //std::cout<<"======= MC matching ======"<<std::endl;
  //int i = 1;
  while ((part->numberOfMothers()>0)) {
    //std::cout<<"pdgid: "<<part->pdgId()<<std::endl;
    //std::cout<<"number of mothers: "<<part->numberOfMothers()<<std::endl;
    //std::cout<<i<<std::endl; i++;

    const reco::Candidate  *MomPart = part->mother();
    //std::cout<<"mother pdgid = "<<MomPart->pdgId()<<std::endl;

    if( fabs(MomPart->pdgId()) != 24 && fabs(MomPart->pdgId()) != 11 && fabs(MomPart->pdgId()) != 13 ) return foundH;

    if (fabs(MomPart->pdgId())==24 && fabs(MomPart->mother()->pdgId()) == 6){ //24 = W+, 25 = H, 6 = top      
      //std::cout<<"number of W mothers: "<<MomPart->numberOfMothers()<<std::endl;
      //std::cout<<"mother of W is "<<MomPart->mother()->pdgId()<<std::endl;
      foundH = true;
      //std::cout<<"matched"<<std::endl;
      break;
    }
    part = MomPart;
  }
  return foundH;
}



void
TriggerMatching::beginEvent()
{
  /////////////////////////////////////////////////////////////////////////////
  //   create the vectors for the tree branches
  ////////////////////////////////////////////////////////////////////////////

  //modules for trigger matching
  T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q = new std::vector<int>;
  T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02 = new std::vector<int>;
  T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02 = new std::vector<int>;
  T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8 = new std::vector<int>;
  T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17 = new std::vector<int>;
  T_hltDiMuonGlb17Glb8DzFiltered0p2 = new std::vector<int>;
  T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17 = new std::vector<int>;
  T_hltDiMuonGlbFiltered17TrkFiltered8 = new std::vector<int>;
  T_hltDiMuonGlb17Trk8DzFiltered0p2 = new std::vector<int>;
  //std::vector<int> *T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17 = new std::vector<int>;
  //std::vector<int> *T_hltDiMuonGlb17Glb8DzFiltered0p2 = new std::vector<int>;
  //T_hltL3MuonRelTrkIsolationVVL = new std::vector<int>;
  T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4 = new std::vector<int>;
  //std::vector<int> *T_hltDiMuonGlbFiltered17TrkFiltered8 = new std::vector<int>;
  //std::vector<int> *T_hltDiMuonGlb17Trk8DzFiltered0p2 = new std::vector<int>;
  //T_hltGlbTrkMuonRelTrkIsolationVVL = new std::vector<int>;
  T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4 = new std::vector<int>;
  T_hltL3NoFiltersfL1sMu12L3Filtered17 = new std::vector<int>;

  T_hltL3fL1sMu25erL1f0L2f25L3Filtered30  = new std::vector<int>;
  T_hltDiMuonGlbFiltered30TrkFiltered11  = new std::vector<int>;
  T_hltDiMuonGlb30Trk11DzFiltered0p2  = new std::vector<int>;
  T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8  = new std::vector<int>;
  T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8  = new std::vector<int>;
  T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter  = new std::vector<int>;
  T_hltL1Mu12EG7L3MuFiltered23  = new std::vector<int>;
  T_hltL1Mu12EG7L3IsoMuFiltered23  = new std::vector<int>;
  T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter  = new std::vector<int>;
  T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter  = new std::vector<int>;
  T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter  = new std::vector<int>;
  T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter  = new std::vector<int>;
  T_hltEle17Ele8GsfTrackIsoLeg1Filter  = new std::vector<int>;
  T_hltEle17Ele8GsfTrackIsoLeg2Filter  = new std::vector<int>;
  T_hltEle27WP80GsfTrackIsoFilter  = new std::vector<int>;

  //reco mu properties
  T_mu_pt	= new std::vector<float>;
  T_mu_eta	= new std::vector<float>;
  T_mu_phi	= new std::vector<float>;
  T_mu_iso      = new std::vector<float>;

  //reco electron properties
  T_ele_pt	= new std::vector<float>;
  T_ele_eta	= new std::vector<float>;
  T_ele_phi	= new std::vector<float>;

}



void
TriggerMatching::endEvent()
{
  /////////////////////////////////////////////////////////////////////////////
  //   delete the vectors for the tree branches
  ////////////////////////////////////////////////////////////////////////////

  //modules for trigger matching
  delete T_hltL3fL1sMu16L1f0L2f16QL3Filtered40Q;
  delete T_hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15IterTrk02;
  delete T_hltL3fL1sMu16L1f0TkFiltered24QL3crIsoRhoFiltered0p15IterTrk02;
  delete T_hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8;
  delete T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  delete T_hltDiMuonGlb17Glb8DzFiltered0p2;
  delete T_hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  delete T_hltDiMuonGlbFiltered17TrkFiltered8;
  delete T_hltDiMuonGlb17Trk8DzFiltered0p2;
  //std::vector<int> *T_hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17;
  //std::vector<int> *T_hltDiMuonGlb17Glb8DzFiltered0p2;
  //delete T_hltL3MuonRelTrkIsolationVVL;
  delete T_hltDiMuonGlb17Glb8DzFiltered0p2RelTrkIsoFiltered0p4;
  //std::vector<int> *T_hltDiMuonGlbFiltered17TrkFiltered8;
  //std::vector<int> *T_hltDiMuonGlb17Trk8DzFiltered0p2;
  //delete T_hltGlbTrkMuonRelTrkIsolationVVL;
  delete T_hltDiMuonGlb17Trk8DzFiltered0p2RelTrkIsoFiltered0p4;
  delete T_hltL3NoFiltersfL1sMu12L3Filtered17;

  delete T_hltL3fL1sMu25erL1f0L2f25L3Filtered30;
  delete T_hltDiMuonGlbFiltered30TrkFiltered11;
  delete T_hltDiMuonGlb30Trk11DzFiltered0p2;
  delete T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8;
  delete T_hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3IsoFiltered8;
  delete T_hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter;
  delete T_hltL1Mu12EG7L3MuFiltered23;
  delete T_hltL1Mu12EG7L3IsoMuFiltered23;
  delete T_hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter;
  delete T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg1Filter;
  delete T_hltEle23Ele12CaloIdTrackIdIsoTrackIsoLeg2Filter;
  delete T_hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter;
  delete T_hltEle17Ele8GsfTrackIsoLeg1Filter;
  delete T_hltEle17Ele8GsfTrackIsoLeg2Filter;
  delete T_hltEle27WP80GsfTrackIsoFilter;



  
  //reco mu properties
  delete T_mu_pt;
  delete T_mu_eta;
  delete T_mu_phi;
  delete T_mu_iso;

  //reco electron properties
  delete T_ele_pt;
  delete T_ele_eta;
  delete T_ele_phi;

}



//define this as a plug-in
DEFINE_FWK_MODULE(TriggerMatching);

