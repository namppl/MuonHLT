/** \class MuonTriggerEfficiencyAnalyzer
 *  Class to measure muon trigger efficiencies (very rough)
 *
 *  $Date: 2012/11/16 15:00:34 $
 *  $Revision: 1.2 $
 *  \authors D. Trocino - Northeastern University <daniele.trocino@cern.ch>
 */

/**
CB list of possible improvements
1. Implement isolation cuts in TAG and PROBE
2. Allow to run on more than one single path at a time
3. Compute turn-ons in |eta|: [0,0.9], [0.9,1.2], [1.2,2.1], [1.2,2.4]
4. Check L1 selection (does not seem to work now)
5. check vertex selection instructions
**/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTEventAnalyzerAOD.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"

#include <map>
#include <string>
#include <memory>
#include <iomanip>

#include "TH1F.h"
#include "TH2F.h"

class MuonTriggerEfficiencyAnalyzer : public edm::EDAnalyzer {

 public:
  /// default constructor
  MuonTriggerEfficiencyAnalyzer(const edm::ParameterSet& cfg);
  /// default destructor
  virtual ~MuonTriggerEfficiencyAnalyzer() {};
  /// everything that needs to be done before the event loop

 private:
  virtual void beginJob();
  /// everything that needs to be done after the event loop
  virtual void endJob();
  /// everything that needs to be done before each run
  virtual void beginRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done after each run
  virtual void endRun(const edm::Run & run, const edm::EventSetup & eventSetup);
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

  /// input tag for mouns
  edm::InputTag vertexes_;
  /// input tag for mouns
  edm::InputTag muons_;
  /// file service
  edm::Service<TFileService> outfile_;
  /// histograms
  std::map<std::string, TH1*> hists_;

  // Trigger process
  std::string triggerProcess_;
  // Trigger names
  std::string tagTriggerName_;
  std::string triggerName_;
  std::string probeFilterL1_;
  std::string probeFilterL2_;
  std::string probeFilterL3_;
  std::string probeFilterL3Iso_;
	
  // Trigger indexes
  int tagTriggerIndex_;
  int triggerIndex_;
  // HLTConfig
  HLTConfigProvider hltConfig_;
  // Rerun
  bool useRerun_;

  // Trigger acceptance
  double tagPtCut_;
  double probePtCut_;
  double tagEtaCut_;
  double probeEtaCut_;

  // ID
  std::string tagId_;
  std::string probeId_;
  bool useIso_;

  // Isolation
  bool isIsolated(reco::Muon muon, std::string isolationType, double isolationCut);
  std::string isolationType_;
  double isolationCut_;

  // Max number of offline muons allowed in the event
  unsigned int nMaxMuons_;

  // Services
  edm::ESHandle<MagneticField> magneticField_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<MuonDetLayerGeometry> detLayerGeometry_;
};

/// default constructor
MuonTriggerEfficiencyAnalyzer::MuonTriggerEfficiencyAnalyzer(const edm::ParameterSet& cfg): 
  vertexes_(cfg.getParameter<edm::InputTag>("vertexes")), 
  muons_(cfg.getParameter<edm::InputTag>("muons")), 
  triggerProcess_(cfg.getParameter<std::string>("triggerProcess")), 
  tagTriggerName_(cfg.getParameter<std::string>("tagTriggerName")), 
  triggerName_(cfg.getParameter<std::string>("triggerName")), 
  probeFilterL1_(cfg.getParameter<std::string>("probeFilterL1")),
  probeFilterL2_(cfg.getParameter<std::string>("probeFilterL2")),
  probeFilterL3_(cfg.getParameter<std::string>("probeFilterL3")),
  probeFilterL3Iso_(cfg.getParameter<std::string>("probeFilterL3Iso")),
  useRerun_(cfg.getParameter<bool>("useRerun")),
  tagPtCut_(cfg.getParameter<double>("tagPtCut")),
  probePtCut_(cfg.getParameter<double>("probePtCut")),
  tagEtaCut_(cfg.getParameter<double>("tagEtaCut")),
  probeEtaCut_(cfg.getParameter<double>("probeEtaCut")),
  tagId_(cfg.getParameter<std::string>("tagID")),
  probeId_(cfg.getParameter<std::string>("probeID")),
  useIso_(cfg.getParameter<bool>("useIso")),
  isolationType_(cfg.getParameter<std::string>("isolationType")),
  isolationCut_(cfg.getParameter<double>("isolationCut")),
  nMaxMuons_(cfg.getUntrackedParameter<unsigned int>("maxNumberMuons", 999999))
{}

void MuonTriggerEfficiencyAnalyzer::beginJob() {
  std::cout<<"beginJobs"<<std::endl;
  std::cout<<"tag pT cut = "<<tagPtCut_<<", probe pT cut = "<<probePtCut_<<std::endl;
  std::cout<<"tag eta cut = "<<tagEtaCut_<<", probe eta cut = "<<probeEtaCut_<<std::endl;
  double eta_bins[] = {-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4}; 
  int eta_bin_n = sizeof(eta_bins)/sizeof(double); 
  double absEta_bins[] = {.0, .9, 1.2, 2.1, 2.4};
  int absEta_bin_n = sizeof(absEta_bins)/sizeof(double);
  std::cout<<"1"<<std::endl;
  hists_["muonPt_tag"]   = outfile_->make<TH1F>("muonPt_tag"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_tag"]  = outfile_->make<TH1F>("muonEta_tag" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_tag"]  = outfile_->make<TH1F>("muonAbsEta_tag" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_tag"]  = outfile_->make<TH1F>("muonPhi_tag" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_tag"] = outfile_->make<TH1F>("nvtx_tag",     "nvtx", 50, 0., 50.);

  std::cout<<"2"<<std::endl;
  hists_["muonPt_tag"]->Sumw2();
  hists_["muonEta_tag"]->Sumw2();
  hists_["muonAbsEta_tag"]->Sumw2();
  hists_["muonPhi_tag"]->Sumw2();
  hists_["muonNvtx_tag"]->Sumw2();

  hists_["muonPt_probe_den"]   = outfile_->make<TH1F>("muonPt_probe_den"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_den"]  = outfile_->make<TH1F>("muonEta_probe_den" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_den"]  = outfile_->make<TH1F>("muonAbsEta_probe_den" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_den"]  = outfile_->make<TH1F>("muonPhi_probe_den" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_den"] = outfile_->make<TH1F>("muonNvtx_probe_den", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_den"]->Sumw2();
  hists_["muonEta_probe_den"]->Sumw2();
  hists_["muonAbsEta_probe_den"]->Sumw2();
  hists_["muonPhi_probe_den"]->Sumw2();
  hists_["muonNvtx_probe_den"]->Sumw2();

  std::cout<<"3"<<std::endl;
  hists_["muonPt_probe_L1"]   = outfile_->make<TH1F>("muonPt_probe_L1"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L1"]  = outfile_->make<TH1F>("muonEta_probe_L1" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L1"]  = outfile_->make<TH1F>("muonAbsEta_probe_L1" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L1"]  = outfile_->make<TH1F>("muonPhi_probe_L1" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L1"] = outfile_->make<TH1F>("muonNvtx_probe_L1", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L1"]->Sumw2();
  hists_["muonEta_probe_L1"]->Sumw2();
  hists_["muonAbsEta_probe_L1"]->Sumw2();
  hists_["muonPhi_probe_L1"]->Sumw2();
  hists_["muonNvtx_probe_L1"]->Sumw2();

  hists_["muonPt_probe_L2"]   = outfile_->make<TH1F>("muonPt_probe_L2"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L2"]  = outfile_->make<TH1F>("muonEta_probe_L2" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L2"]  = outfile_->make<TH1F>("muonAbsEta_probe_L2" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L2"]  = outfile_->make<TH1F>("muonPhi_probe_L2" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L2"] = outfile_->make<TH1F>("muonNvtx_probe_L2", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L2"]->Sumw2();
  hists_["muonEta_probe_L2"]->Sumw2();
  hists_["muonAbsEta_probe_L2"]->Sumw2();
  hists_["muonPhi_probe_L2"]->Sumw2();
  hists_["muonNvtx_probe_L2"]->Sumw2();
  
  hists_["muonPt_probe_L3"]   = outfile_->make<TH1F>("muonPt_probe_L3"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3"]  = outfile_->make<TH1F>("muonEta_probe_L3" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3"]  = outfile_->make<TH1F>("muonPhi_probe_L3" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3"] = outfile_->make<TH1F>("muonNvtx_probe_L3", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3"]->Sumw2();
  hists_["muonEta_probe_L3"]->Sumw2();
  hists_["muonAbsEta_probe_L3"]->Sumw2();
  hists_["muonPhi_probe_L3"]->Sumw2();
  hists_["muonNvtx_probe_L3"]->Sumw2();
  
  hists_["muonPt_probe_L3OIState"]   = outfile_->make<TH1F>("muonPt_probe_L3OIState"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIState"]  = outfile_->make<TH1F>("muonEta_probe_L3OIState" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIState"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIState" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIState"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIState" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIState"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIState", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIState"]->Sumw2();
  hists_["muonEta_probe_L3OIState"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIState"]->Sumw2();
  hists_["muonPhi_probe_L3OIState"]->Sumw2();
  hists_["muonNvtx_probe_L3OIState"]->Sumw2();
  
  hists_["muonPt_probe_L3OIHit"]   = outfile_->make<TH1F>("muonPt_probe_L3OIHit"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIHit"]  = outfile_->make<TH1F>("muonEta_probe_L3OIHit" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIHit"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIHit" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIHit"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIHit" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIHit"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIHit", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIHit"]->Sumw2();
  hists_["muonEta_probe_L3OIHit"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIHit"]->Sumw2();
  hists_["muonPhi_probe_L3OIHit"]->Sumw2();
  hists_["muonNvtx_probe_L3OIHit"]->Sumw2();
  
  hists_["muonPt_probe_L3IOHit"]   = outfile_->make<TH1F>("muonPt_probe_L3IOHit"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3IOHit"]  = outfile_->make<TH1F>("muonEta_probe_L3IOHit" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3IOHit"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3IOHit" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3IOHit"]  = outfile_->make<TH1F>("muonPhi_probe_L3IOHit" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3IOHit"] = outfile_->make<TH1F>("muonNvtx_probe_L3IOHit", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3IOHit"]->Sumw2();
  hists_["muonEta_probe_L3IOHit"]->Sumw2();
  hists_["muonAbsEta_probe_L3IOHit"]->Sumw2();
  hists_["muonPhi_probe_L3IOHit"]->Sumw2();
  hists_["muonNvtx_probe_L3IOHit"]->Sumw2();

  hists_["muonPt_probe_L3Global"]   = outfile_->make<TH1F>("muonPt_probe_L3Global"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3Global"]  = outfile_->make<TH1F>("muonEta_probe_L3Global" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3Global"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3Global" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3Global"]  = outfile_->make<TH1F>("muonPhi_probe_L3Global" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3Global"] = outfile_->make<TH1F>("muonNvtx_probe_L3Global", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3Global"]->Sumw2();
  hists_["muonEta_probe_L3Global"]->Sumw2();
  hists_["muonAbsEta_probe_L3Global"]->Sumw2();
  hists_["muonPhi_probe_L3Global"]->Sumw2();
  hists_["muonNvtx_probe_L3Global"]->Sumw2();

  hists_["muonPt_probe_L3Iso"]   = outfile_->make<TH1F>("muonPt_probe_L3Iso"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3Iso"]  = outfile_->make<TH1F>("muonEta_probe_L3Iso" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3Iso"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3Iso" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3Iso"]  = outfile_->make<TH1F>("muonPhi_probe_L3Iso" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3Iso"] = outfile_->make<TH1F>("muonNvtx_probe_L3Iso", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3Iso"]->Sumw2();
  hists_["muonEta_probe_L3Iso"]->Sumw2();
  hists_["muonAbsEta_probe_L3Iso"]->Sumw2();
  hists_["muonPhi_probe_L3Iso"]->Sumw2();
  hists_["muonNvtx_probe_L3Iso"]->Sumw2();
  
  hists_["muonPt_probe_L3OIStateGlobal"]   = outfile_->make<TH1F>("muonPt_probe_L3OIStateGlobal"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIStateGlobal"]  = outfile_->make<TH1F>("muonEta_probe_L3OIStateGlobal" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIStateGlobal"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIStateGlobal" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIStateGlobal"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIStateGlobal" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIStateGlobal"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIStateGlobal", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIStateGlobal"]->Sumw2();
  hists_["muonEta_probe_L3OIStateGlobal"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIStateGlobal"]->Sumw2();
  hists_["muonPhi_probe_L3OIStateGlobal"]->Sumw2();
  hists_["muonNvtx_probe_L3OIStateGlobal"]->Sumw2();
  
  hists_["muonPt_probe_L3OIHitGlobal"]   = outfile_->make<TH1F>("muonPt_probe_L3OIHitGlobal"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIHitGlobal"]  = outfile_->make<TH1F>("muonEta_probe_L3OIHitGlobal" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIHitGlobal"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIHitGlobal" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIHitGlobal"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIHitGlobal" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIHitGlobal"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIHitGlobal", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIHitGlobal"]->Sumw2();
  hists_["muonEta_probe_L3OIHitGlobal"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIHitGlobal"]->Sumw2();
  hists_["muonPhi_probe_L3OIHitGlobal"]->Sumw2();
  hists_["muonNvtx_probe_L3OIHitGlobal"]->Sumw2();
  
  hists_["muonPt_probe_L3IOHitGlobal"]   = outfile_->make<TH1F>("muonPt_probe_L3IOHitGlobal"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3IOHitGlobal"]  = outfile_->make<TH1F>("muonEta_probe_L3IOHitGlobal" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3IOHitGlobal"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3IOHitGlobal" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3IOHitGlobal"]  = outfile_->make<TH1F>("muonPhi_probe_L3IOHitGlobal" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3IOHitGlobal"] = outfile_->make<TH1F>("muonNvtx_probe_L3IOHitGlobal", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3IOHitGlobal"]->Sumw2();
  hists_["muonEta_probe_L3IOHitGlobal"]->Sumw2();
  hists_["muonAbsEta_probe_L3IOHitGlobal"]->Sumw2();
  hists_["muonPhi_probe_L3IOHitGlobal"]->Sumw2();
  hists_["muonNvtx_probe_L3IOHitGlobal"]->Sumw2();
  
  hists_["muonPt_probe_L3TkTrack"]   = outfile_->make<TH1F>("muonPt_probe_L3TkTrack"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3TkTrack"]  = outfile_->make<TH1F>("muonEta_probe_L3TkTrack" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3TkTrack"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3TkTrack" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3TkTrack"]  = outfile_->make<TH1F>("muonPhi_probe_L3TkTrack" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3TkTrack"] = outfile_->make<TH1F>("muonNvtx_probe_L3TkTrack", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3TkTrack"]->Sumw2();
  hists_["muonEta_probe_L3TkTrack"]->Sumw2();
  hists_["muonAbsEta_probe_L3TkTrack"]->Sumw2();
  hists_["muonPhi_probe_L3TkTrack"]->Sumw2();
  hists_["muonNvtx_probe_L3TkTrack"]->Sumw2();
  
  hists_["muonPt_probe_L3OIStateTkTrack"]   = outfile_->make<TH1F>("muonPt_probe_L3OIStateTkTrack"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIStateTkTrack"]  = outfile_->make<TH1F>("muonEta_probe_L3OIStateTkTrack" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIStateTkTrack"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIStateTkTrack" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIStateTkTrack"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIStateTkTrack" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIStateTkTrack"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIStateTkTrack", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIStateTkTrack"]->Sumw2();
  hists_["muonEta_probe_L3OIStateTkTrack"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIStateTkTrack"]->Sumw2();
  hists_["muonPhi_probe_L3OIStateTkTrack"]->Sumw2();
  hists_["muonNvtx_probe_L3OIStateTkTrack"]->Sumw2();
  
  hists_["muonPt_probe_L3OIHitTkTrack"]   = outfile_->make<TH1F>("muonPt_probe_L3OIHitTkTrack"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIHitTkTrack"]  = outfile_->make<TH1F>("muonEta_probe_L3OIHitTkTrack" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIHitTkTrack"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIHitTkTrack" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIHitTkTrack"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIHitTkTrack" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIHitTkTrack"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIHitTkTrack", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIHitTkTrack"]->Sumw2();
  hists_["muonEta_probe_L3OIHitTkTrack"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIHitTkTrack"]->Sumw2();
  hists_["muonPhi_probe_L3OIHitTkTrack"]->Sumw2();
  hists_["muonNvtx_probe_L3OIHitTkTrack"]->Sumw2();
  
  hists_["muonPt_probe_L3IOHitTkTrack"]   = outfile_->make<TH1F>("muonPt_probe_L3IOHitTkTrack"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3IOHitTkTrack"]  = outfile_->make<TH1F>("muonEta_probe_L3IOHitTkTrack" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3IOHitTkTrack"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3IOHitTkTrack" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3IOHitTkTrack"]  = outfile_->make<TH1F>("muonPhi_probe_L3IOHitTkTrack" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3IOHitTkTrack"] = outfile_->make<TH1F>("muonNvtx_probe_L3IOHitTkTrack", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3IOHitTkTrack"]->Sumw2();
  hists_["muonEta_probe_L3IOHitTkTrack"]->Sumw2();
  hists_["muonAbsEta_probe_L3IOHitTkTrack"]->Sumw2();
  hists_["muonPhi_probe_L3IOHitTkTrack"]->Sumw2();
  hists_["muonNvtx_probe_L3IOHitTkTrack"]->Sumw2();

  hists_["muonPt_probe_L3Seed"]   = outfile_->make<TH1F>("muonPt_probe_L3Seed"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3Seed"]  = outfile_->make<TH1F>("muonEta_probe_L3Seed" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3Seed"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3Seed" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3Seed"]  = outfile_->make<TH1F>("muonPhi_probe_L3Seed" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3Seed"] = outfile_->make<TH1F>("muonNvtx_probe_L3Seed", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3Seed"]->Sumw2();
  hists_["muonEta_probe_L3Seed"]->Sumw2();
  hists_["muonAbsEta_probe_L3Seed"]->Sumw2();
  hists_["muonPhi_probe_L3Seed"]->Sumw2();
  hists_["muonNvtx_probe_L3Seed"]->Sumw2();
  
  hists_["muonPt_probe_L3OIStateSeed"]   = outfile_->make<TH1F>("muonPt_probe_L3OIStateSeed"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIStateSeed"]  = outfile_->make<TH1F>("muonEta_probe_L3OIStateSeed" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIStateSeed"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIStateSeed" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIStateSeed"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIStateSeed" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIStateSeed"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIStateSeed", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIStateSeed"]->Sumw2();
  hists_["muonEta_probe_L3OIStateSeed"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIStateSeed"]->Sumw2();
  hists_["muonPhi_probe_L3OIStateSeed"]->Sumw2();
  hists_["muonNvtx_probe_L3OIStateSeed"]->Sumw2();
  
  hists_["muonPt_probe_L3OIHitSeed"]   = outfile_->make<TH1F>("muonPt_probe_L3OIHitSeed"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3OIHitSeed"]  = outfile_->make<TH1F>("muonEta_probe_L3OIHitSeed" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3OIHitSeed"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3OIHitSeed" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3OIHitSeed"]  = outfile_->make<TH1F>("muonPhi_probe_L3OIHitSeed" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3OIHitSeed"] = outfile_->make<TH1F>("muonNvtx_probe_L3OIHitSeed", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3OIHitSeed"]->Sumw2();
  hists_["muonEta_probe_L3OIHitSeed"]->Sumw2();
  hists_["muonAbsEta_probe_L3OIHitSeed"]->Sumw2();
  hists_["muonPhi_probe_L3OIHitSeed"]->Sumw2();
  hists_["muonNvtx_probe_L3OIHitSeed"]->Sumw2();
  
  hists_["muonPt_probe_L3IOHitSeed"]   = outfile_->make<TH1F>("muonPt_probe_L3IOHitSeed"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe_L3IOHitSeed"]  = outfile_->make<TH1F>("muonEta_probe_L3IOHitSeed" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe_L3IOHitSeed"]  = outfile_->make<TH1F>("muonAbsEta_probe_L3IOHitSeed" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe_L3IOHitSeed"]  = outfile_->make<TH1F>("muonPhi_probe_L3IOHitSeed" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe_L3IOHitSeed"] = outfile_->make<TH1F>("muonNvtx_probe_L3IOHitSeed", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe_L3IOHitSeed"]->Sumw2();
  hists_["muonEta_probe_L3IOHitSeed"]->Sumw2();
  hists_["muonAbsEta_probe_L3IOHitSeed"]->Sumw2();
  hists_["muonPhi_probe_L3IOHitSeed"]->Sumw2();
  hists_["muonNvtx_probe_L3IOHitSeed"]->Sumw2();


  std::cout<<"4"<<std::endl;
  hists_["muonPt_probe1"]   = outfile_->make<TH1F>("muonPt_probe1"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe1"]  = outfile_->make<TH1F>("muonEta_probe1" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe1"]  = outfile_->make<TH1F>("muonAbsEta_probe1" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe1"]  = outfile_->make<TH1F>("muonPhi_probe1" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe1"] = outfile_->make<TH1F>("muonNvtx_probe1", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe1"]->Sumw2();
  hists_["muonEta_probe1"]->Sumw2();
  hists_["muonAbsEta_probe1"]->Sumw2();
  hists_["muonPhi_probe1"]->Sumw2();
  hists_["muonNvtx_probe1"]->Sumw2();

  hists_["muonPt_probe2"]   = outfile_->make<TH1F>("muonPt_probe2"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe2"]  = outfile_->make<TH1F>("muonEta_probe2" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe2"]  = outfile_->make<TH1F>("muonAbsEta_probe2" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe2"]  = outfile_->make<TH1F>("muonPhi_probe2" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe2"] = outfile_->make<TH1F>("muonNvtx_probe2", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe2"]->Sumw2();
  hists_["muonEta_probe2"]->Sumw2();
  hists_["muonAbsEta_probe2"]->Sumw2();
  hists_["muonPhi_probe2"]->Sumw2();
  hists_["muonNvtx_probe2"]->Sumw2();

  hists_["muonPt_probe3"]   = outfile_->make<TH1F>("muonPt_probe3"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe3"]  = outfile_->make<TH1F>("muonEta_probe3" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe3"]  = outfile_->make<TH1F>("muonAbsEta_probe3" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe3"]  = outfile_->make<TH1F>("muonPhi_probe3" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe3"] = outfile_->make<TH1F>("muonNvtx_probe3", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe3"]->Sumw2();
  hists_["muonEta_probe3"]->Sumw2();
  hists_["muonAbsEta_probe3"]->Sumw2();
  hists_["muonPhi_probe3"]->Sumw2();
  hists_["muonNvtx_probe3"]->Sumw2();

  hists_["muonPt_probe4"]   = outfile_->make<TH1F>("muonPt_probe4"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe4"]  = outfile_->make<TH1F>("muonEta_probe4" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe4"]  = outfile_->make<TH1F>("muonAbsEta_probe4" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe4"]  = outfile_->make<TH1F>("muonPhi_probe4" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe4"] = outfile_->make<TH1F>("muonNvtx_probe4", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe4"]->Sumw2();
  hists_["muonEta_probe4"]->Sumw2();
  hists_["muonAbsEta_probe4"]->Sumw2();
  hists_["muonPhi_probe4"]->Sumw2();
  hists_["muonNvtx_probe4"]->Sumw2();

  hists_["muonPt_probe5"]   = outfile_->make<TH1F>("muonPt_probe5"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe5"]  = outfile_->make<TH1F>("muonEta_probe5" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe5"]  = outfile_->make<TH1F>("muonAbsEta_probe5" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe5"]  = outfile_->make<TH1F>("muonPhi_probe5" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe5"] = outfile_->make<TH1F>("muonNvtx_probe5", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe5"]->Sumw2();
  hists_["muonEta_probe5"]->Sumw2();
  hists_["muonAbsEta_probe5"]->Sumw2();
  hists_["muonPhi_probe5"]->Sumw2();
  hists_["muonNvtx_probe5"]->Sumw2();

  hists_["muonPt_probe6"]   = outfile_->make<TH1F>("muonPt_probe6"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe6"]  = outfile_->make<TH1F>("muonEta_probe6" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe6"]  = outfile_->make<TH1F>("muonAbsEta_probe6" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe6"]  = outfile_->make<TH1F>("muonPhi_probe6" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe6"] = outfile_->make<TH1F>("muonNvtx_probe6", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe6"]->Sumw2();
  hists_["muonEta_probe6"]->Sumw2();
  hists_["muonAbsEta_probe6"]->Sumw2();
  hists_["muonPhi_probe6"]->Sumw2();
  hists_["muonNvtx_probe6"]->Sumw2();

  hists_["muonPt_probe7"]   = outfile_->make<TH1F>("muonPt_probe7"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe7"]  = outfile_->make<TH1F>("muonEta_probe7" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe7"]  = outfile_->make<TH1F>("muonAbsEta_probe7" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe7"]  = outfile_->make<TH1F>("muonPhi_probe7" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe7"] = outfile_->make<TH1F>("muonNvtx_probe7", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe7"]->Sumw2();
  hists_["muonEta_probe7"]->Sumw2();
  hists_["muonAbsEta_probe7"]->Sumw2();
  hists_["muonPhi_probe7"]->Sumw2();
  hists_["muonNvtx_probe7"]->Sumw2();

  hists_["muonPt_probe8"]   = outfile_->make<TH1F>("muonPt_probe8"  , "pt"  ,  40,  0., 200.);
  hists_["muonEta_probe8"]  = outfile_->make<TH1F>("muonEta_probe8" , "eta" ,  eta_bin_n-1, eta_bins);
  hists_["muonAbsEta_probe8"]  = outfile_->make<TH1F>("muonAbsEta_probe8" , "|eta|" ,  absEta_bin_n-1, absEta_bins);
  hists_["muonPhi_probe8"]  = outfile_->make<TH1F>("muonPhi_probe8" , "phi" ,  100, -5.,   5.); 
  hists_["muonNvtx_probe8"] = outfile_->make<TH1F>("muonNvtx_probe8", "nvtx", 50, 0., 50.);

  hists_["muonPt_probe8"]->Sumw2();
  hists_["muonEta_probe8"]->Sumw2();
  hists_["muonAbsEta_probe8"]->Sumw2();
  hists_["muonPhi_probe8"]->Sumw2();
  hists_["muonNvtx_probe8"]->Sumw2();

  std::cout<<"5"<<std::endl;
  hists_["mumuMass_all"] = outfile_->make<TH1F>("mumuMass_all", "mass",   90, 30., 120.);
  hists_["mumuMass_den"] = outfile_->make<TH1F>("mumuMass_den", "mass",   90, 30., 120.);
  hists_["mumuMass_L1"] = outfile_->make<TH1F>("mumuMass_L1", "mass",   90, 30., 120.);
  hists_["mumuMass_L2"] = outfile_->make<TH1F>("mumuMass_L2", "mass",   90, 30., 120.);
  hists_["mumuMass_L3"] = outfile_->make<TH1F>("mumuMass_L3", "mass",   90, 30., 120.);
  hists_["mumuMass_L3Seed"] = outfile_->make<TH1F>("mumuMass_L3Seed", "mass",   90, 30., 120.);
  hists_["mumuMass_L3OIStateSeed"] = outfile_->make<TH1F>("mumuMass_L3OIStateSeed", "mass",   90, 30., 120.);
  hists_["mumuMass_L3OIHitSeed"] = outfile_->make<TH1F>("mumuMass_L3OIHitSeed", "mass",   90, 30., 120.);
  hists_["mumuMass_L3IOHitSeed"] = outfile_->make<TH1F>("mumuMass_L3IOHitSeed", "mass",   90, 30., 120.);
  hists_["mumuMass_L3TkTrack"] = outfile_->make<TH1F>("mumuMass_L3TkTrack", "mass",   90, 30., 120.);
  hists_["mumuMass_L3OIStateTkTrack"] = outfile_->make<TH1F>("mumuMass_L3OIStateTkTrack", "mass",   90, 30., 120.);
  hists_["mumuMass_L3OIHitTkTrack"] = outfile_->make<TH1F>("mumuMass_L3OIHitTkTrack", "mass",   90, 30., 120.);
  hists_["mumuMass_L3IOHitTkTrack"] = outfile_->make<TH1F>("mumuMass_L3IOHitTkTrack", "mass",   90, 30., 120.);
  hists_["mumuMass_L3Global"] = outfile_->make<TH1F>("mumuMass_L3Global", "mass",   90, 30., 120.);
  hists_["mumuMass_L3OIStateGlobal"] = outfile_->make<TH1F>("mumuMass_L3OIStateGlobal", "mass",   90, 30., 120.);
  hists_["mumuMass_L3OIHitGlobal"] = outfile_->make<TH1F>("mumuMass_L3OIHitGlobal", "mass",   90, 30., 120.);
  hists_["mumuMass_L3IOHitGlobal"] = outfile_->make<TH1F>("mumuMass_L3IOHitGlobal", "mass",   90, 30., 120.);
  hists_["mumuMass1"]     = outfile_->make<TH1F>("mumuMass1",     "mass",   90, 30., 120.);
  hists_["mumuMass2"]     = outfile_->make<TH1F>("mumuMass2",     "mass",   90, 30., 120.);
  hists_["mumuMass3"]     = outfile_->make<TH1F>("mumuMass3",     "mass",   90, 30., 120.);
  hists_["mumuMass4"]     = outfile_->make<TH1F>("mumuMass4",     "mass",   90, 30., 120.);

  hists_["mumuMass_all"]->Sumw2();
  hists_["mumuMass_den"]->Sumw2();
  hists_["mumuMass_L1"]->Sumw2();
  hists_["mumuMass_L2"]->Sumw2();
  hists_["mumuMass_L3"]->Sumw2();
  hists_["mumuMass_L3Seed"]->Sumw2();
  hists_["mumuMass_L3OIStateSeed"]->Sumw2();
  hists_["mumuMass_L3OIHitSeed"]->Sumw2();
  hists_["mumuMass_L3IOHitSeed"]->Sumw2();
  hists_["mumuMass_L3TkTrack"]->Sumw2();
  hists_["mumuMass_L3OIStateTkTrack"]->Sumw2();
  hists_["mumuMass_L3OIHitTkTrack"]->Sumw2();
  hists_["mumuMass_L3IOHitTkTrack"]->Sumw2();
  hists_["mumuMass_L3Global"]->Sumw2();
  hists_["mumuMass_L3OIStateGlobal"]->Sumw2();
  hists_["mumuMass_L3OIHitGlobal"]->Sumw2();
  hists_["mumuMass_L3IOHitGlobal"]->Sumw2();
  hists_["mumuMass1"]->Sumw2();
  hists_["mumuMass2"]->Sumw2();
  hists_["mumuMass3"]->Sumw2();
  hists_["mumuMass4"]->Sumw2();

  std::cout<<"6"<<std::endl;
  hists_["muonPt12_den"]  = outfile_->make<TH2F>("muonPt12_den"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_den"] = outfile_->make<TH2F>("muonEta12_den" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_den"] = outfile_->make<TH2F>("muonPhi12_den" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_den"]->Sumw2();
  hists_["muonEta12_den"]->Sumw2();
  hists_["muonPhi12_den"]->Sumw2();

  hists_["muonPt12_L1"]  = outfile_->make<TH2F>("muonPt12_L1"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_L1"] = outfile_->make<TH2F>("muonEta12_L1" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_L1"] = outfile_->make<TH2F>("muonPhi12_L1" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_L1"]->Sumw2();
  hists_["muonEta12_L1"]->Sumw2();
  hists_["muonPhi12_L1"]->Sumw2();

  hists_["muonPt12_L2"]  = outfile_->make<TH2F>("muonPt12_L2"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_L2"] = outfile_->make<TH2F>("muonEta12_L2" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_L2"] = outfile_->make<TH2F>("muonPhi12_L2" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_L2"]->Sumw2();
  hists_["muonEta12_L2"]->Sumw2();
  hists_["muonPhi12_L2"]->Sumw2();
  
  hists_["muonPt12_L3"]  = outfile_->make<TH2F>("muonPt12_L3"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_L3"] = outfile_->make<TH2F>("muonEta12_L3" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_L3"] = outfile_->make<TH2F>("muonPhi12_L3" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_L3"]->Sumw2();
  hists_["muonEta12_L3"]->Sumw2();
  hists_["muonPhi12_L3"]->Sumw2();
   
  hists_["muonPt12_L3Seed"]  = outfile_->make<TH2F>("muonPt12_L3Seed"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_L3Seed"] = outfile_->make<TH2F>("muonEta12_L3Seed" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_L3Seed"] = outfile_->make<TH2F>("muonPhi12_L3Seed" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_L3Seed"]->Sumw2();
  hists_["muonEta12_L3Seed"]->Sumw2();
  hists_["muonPhi12_L3Seed"]->Sumw2();
  
  hists_["muonPt12_L3TkTrack"]  = outfile_->make<TH2F>("muonPt12_L3TkTrack"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_L3TkTrack"] = outfile_->make<TH2F>("muonEta12_L3TkTrack" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_L3TkTrack"] = outfile_->make<TH2F>("muonPhi12_L3TkTrack" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_L3TkTrack"]->Sumw2();
  hists_["muonEta12_L3TkTrack"]->Sumw2();
  hists_["muonPhi12_L3TkTrack"]->Sumw2();
  
  hists_["muonPt12_L3Global"]  = outfile_->make<TH2F>("muonPt12_L3Global"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_L3Global"] = outfile_->make<TH2F>("muonEta12_L3Global" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_L3Global"] = outfile_->make<TH2F>("muonPhi12_L3Global" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_L3Global"]->Sumw2();
  hists_["muonEta12_L3Global"]->Sumw2();
  hists_["muonPhi12_L3Global"]->Sumw2();
  std::cout<<"7"<<std::endl;
  hists_["muonPt12_1"]  = outfile_->make<TH2F>("muonPt12_1"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_1"] = outfile_->make<TH2F>("muonEta12_1" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_1"] = outfile_->make<TH2F>("muonPhi12_1" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_1"]->Sumw2();
  hists_["muonEta12_1"]->Sumw2();
  hists_["muonPhi12_1"]->Sumw2();
  
  hists_["muonPt12_2"]  = outfile_->make<TH2F>("muonPt12_2"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_2"] = outfile_->make<TH2F>("muonEta12_2" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_2"] = outfile_->make<TH2F>("muonPhi12_2" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_2"]->Sumw2();
  hists_["muonEta12_2"]->Sumw2();
  hists_["muonPhi12_2"]->Sumw2();
    
  hists_["muonPt12_3"]  = outfile_->make<TH2F>("muonPt12_3"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_3"] = outfile_->make<TH2F>("muonEta12_3" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_3"] = outfile_->make<TH2F>("muonPhi12_3" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_3"]->Sumw2();
  hists_["muonEta12_3"]->Sumw2();
  hists_["muonPhi12_3"]->Sumw2();
      
  hists_["muonPt12_4"]  = outfile_->make<TH2F>("muonPt12_4"  , "p_{T,1} vs p_{T,2}",  80, 20., 100., 80, 20., 100.);
  hists_["muonEta12_4"] = outfile_->make<TH2F>("muonEta12_4" , "eta" , eta_bin_n-1, eta_bins, eta_bin_n-1, eta_bins);
  hists_["muonPhi12_4"] = outfile_->make<TH2F>("muonPhi12_4" , "phi" , 16, -3.2, 3.2, 16, -3.2, 3.2); 

  hists_["muonPt12_4"]->Sumw2();
  hists_["muonEta12_4"]->Sumw2();
  hists_["muonPhi12_4"]->Sumw2();
  
  hists_["deltaR_trobj_tag"]   = outfile_->make<TH1F>("deltaR_trobj_tag" ,   "#DeltaR(trig,#mu)" , 600, 0., 6.0); 
  hists_["deltaR_trobj_probe_L1"] = outfile_->make<TH1F>("deltaR_trobj_probe_L1" , "#DeltaR(trig,#mu)" , 600, 0., 6.0); 
  hists_["deltaR_trobj_probe_L2"] = outfile_->make<TH1F>("deltaR_trobj_probe_L2" , "#DeltaR(trig,#mu)" , 600, 0., 6.0);  
  hists_["deltaR_trobj_probe_L3"] = outfile_->make<TH1F>("deltaR_trobj_probe_L3" , "#DeltaR(trig,#mu)" , 600, 0., 6.0);  
}

void MuonTriggerEfficiencyAnalyzer::endJob() {
  hists_["muonPt_probe1"]->Divide(   hists_["muonPt_probe_L1"],   hists_["muonPt_probe_den"],   1.0, 1.0, "B" );
  hists_["muonEta_probe1"]->Divide(  hists_["muonEta_probe_L1"],  hists_["muonEta_probe_den"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe1"]->Divide(  hists_["muonAbsEta_probe_L1"],  hists_["muonAbsEta_probe_den"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe1"]->Divide(  hists_["muonPhi_probe_L1"],  hists_["muonPhi_probe_den"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe1"]->Divide( hists_["muonNvtx_probe_L1"], hists_["muonNvtx_probe_den"], 1.0, 1.0, "B" );
  
  hists_["muonPt_probe2"]->Divide(   hists_["muonPt_probe_L2"],   hists_["muonPt_probe_L1"],   1.0, 1.0, "B" );
  hists_["muonEta_probe2"]->Divide(  hists_["muonEta_probe_L2"],  hists_["muonEta_probe_L1"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe2"]->Divide(  hists_["muonAbsEta_probe_L2"],  hists_["muonAbsEta_probe_L1"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe2"]->Divide(  hists_["muonPhi_probe_L2"],  hists_["muonPhi_probe_L1"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe2"]->Divide( hists_["muonNvtx_probe_L2"], hists_["muonNvtx_probe_L1"], 1.0, 1.0, "B" );
  
  hists_["muonPt_probe3"]->Divide(   hists_["muonPt_probe_L3Seed"],   hists_["muonPt_probe_L2"],   1.0, 1.0, "B" );
  hists_["muonEta_probe3"]->Divide(  hists_["muonEta_probe_L3Seed"],  hists_["muonEta_probe_L2"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe3"]->Divide(  hists_["muonAbsEta_probe_L3Seed"],  hists_["muonAbsEta_probe_L2"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe3"]->Divide(  hists_["muonPhi_probe_L3Seed"],  hists_["muonPhi_probe_L2"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe3"]->Divide( hists_["muonNvtx_probe_L3Seed"], hists_["muonNvtx_probe_L2"], 1.0, 1.0, "B" ); 

  hists_["muonPt_probe4"]->Divide(   hists_["muonPt_probe_L3TkTrack"],   hists_["muonPt_probe_L3Seed"],   1.0, 1.0, "B" );
  hists_["muonEta_probe4"]->Divide(  hists_["muonEta_probe_L3TkTrack"],  hists_["muonEta_probe_L3Seed"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe4"]->Divide(  hists_["muonAbsEta_probe_L3TkTrack"],  hists_["muonAbsEta_probe_L3Seed"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe4"]->Divide(  hists_["muonPhi_probe_L3TkTrack"],  hists_["muonPhi_probe_L3Seed"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe4"]->Divide( hists_["muonNvtx_probe_L3TkTrack"], hists_["muonNvtx_probe_L3Seed"], 1.0, 1.0, "B" );

  hists_["muonPt_probe5"]->Divide(   hists_["muonPt_probe_L3Global"],   hists_["muonPt_probe_L3TkTrack"],   1.0, 1.0, "B" );
  hists_["muonEta_probe5"]->Divide(  hists_["muonEta_probe_L3Global"],  hists_["muonEta_probe_L3TkTrack"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe5"]->Divide(  hists_["muonAbsEta_probe_L3Global"],  hists_["muonAbsEta_probe_L3TkTrack"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe5"]->Divide(  hists_["muonPhi_probe_L3Global"],  hists_["muonPhi_probe_L3TkTrack"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe5"]->Divide( hists_["muonNvtx_probe_L3Global"], hists_["muonNvtx_probe_L3TkTrack"], 1.0, 1.0, "B" );

  hists_["muonPt_probe6"]->Divide(   hists_["muonPt_probe_L3"],   hists_["muonPt_probe_L3Global"],   1.0, 1.0, "B" );
  hists_["muonEta_probe6"]->Divide(  hists_["muonEta_probe_L3"],  hists_["muonEta_probe_L3Global"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe6"]->Divide(  hists_["muonAbsEta_probe_L3"],  hists_["muonAbsEta_probe_L3Global"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe6"]->Divide(  hists_["muonPhi_probe_L3"],  hists_["muonPhi_probe_L3Global"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe6"]->Divide( hists_["muonNvtx_probe_L3"], hists_["muonNvtx_probe_L3Global"], 1.0, 1.0, "B" );

  hists_["muonPt_probe7"]->Divide(   hists_["muonPt_probe_L3"],   hists_["muonPt_probe_L2"],   1.0, 1.0, "B" );
  hists_["muonEta_probe7"]->Divide(  hists_["muonEta_probe_L3"],  hists_["muonEta_probe_L2"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe7"]->Divide(  hists_["muonAbsEta_probe_L3"],  hists_["muonAbsEta_probe_L2"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe7"]->Divide(  hists_["muonPhi_probe_L3"],  hists_["muonPhi_probe_L2"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe7"]->Divide( hists_["muonNvtx_probe_L3"], hists_["muonNvtx_probe_L2"], 1.0, 1.0, "B" );

  hists_["muonPt_probe8"]->Divide(   hists_["muonPt_probe_L3"],   hists_["muonPt_probe_den"],   1.0, 1.0, "B" );
  hists_["muonEta_probe8"]->Divide(  hists_["muonEta_probe_L3"],  hists_["muonEta_probe_den"],  1.0, 1.0, "B" );
  hists_["muonAbsEta_probe8"]->Divide(  hists_["muonAbsEta_probe_L3"],  hists_["muonAbsEta_probe_den"],  1.0, 1.0, "B" );
  hists_["muonPhi_probe8"]->Divide(  hists_["muonPhi_probe_L3"],  hists_["muonPhi_probe_den"],  1.0, 1.0, "B" );
  hists_["muonNvtx_probe8"]->Divide( hists_["muonNvtx_probe_L3"], hists_["muonNvtx_probe_den"], 1.0, 1.0, "B" );
/*
  hists_["mumuMass1"]->Divide(  hists_["mumuMass_L1"],  hists_["mumuMass_den"],  1.0, 1.0, "B" ); 
  hists_["mumuMass2"]->Divide(  hists_["mumuMass_L2"],  hists_["mumuMass_L1"],  1.0, 1.0, "B" ); 
  hists_["mumuMass3"]->Divide(  hists_["mumuMass_L3"],  hists_["mumuMass_L2"],  1.0, 1.0, "B" ); 
  hists_["mumuMass4"]->Divide(  hists_["mumuMass_L3"],  hists_["mumuMass_den"],  1.0, 1.0, "B" ); 

  hists_["muonPt12_1"]->Divide(  hists_["muonPt12_L1"],  hists_["muonPt12_den"],  1.0, 1.0, "B" ); 
  hists_["muonEta12_1"]->Divide( hists_["muonEta12_L1"], hists_["muonEta12_den"], 1.0, 1.0, "B" ); 
  hists_["muonPhi12_1"]->Divide( hists_["muonPhi12_L1"], hists_["muonPhi12_den"], 1.0, 1.0, "B" );
  
  hists_["muonPt12_2"]->Divide(  hists_["muonPt12_L2"],  hists_["muonPt12_L1"],  1.0, 1.0, "B" ); 
  hists_["muonEta12_2"]->Divide( hists_["muonEta12_L2"], hists_["muonEta12_L1"], 1.0, 1.0, "B" ); 
  hists_["muonPhi12_2"]->Divide( hists_["muonPhi12_L2"], hists_["muonPhi12_L1"], 1.0, 1.0, "B" );
  
  hists_["muonPt12_3"]->Divide(  hists_["muonPt12_L3"],  hists_["muonPt12_L2"],  1.0, 1.0, "B" ); 
  hists_["muonEta12_3"]->Divide( hists_["muonEta12_L3"], hists_["muonEta12_L2"], 1.0, 1.0, "B" ); 
  hists_["muonPhi12_3"]->Divide( hists_["muonPhi12_L3"], hists_["muonPhi12_L2"], 1.0, 1.0, "B" );
  
  hists_["muonPt12_4"]->Divide(  hists_["muonPt12_L3"],  hists_["muonPt12_den"],  1.0, 1.0, "B" ); 
  hists_["muonEta12_4"]->Divide( hists_["muonEta12_L3"], hists_["muonEta12_den"], 1.0, 1.0, "B" ); 
  hists_["muonPhi12_4"]->Divide( hists_["muonPhi12_L3"], hists_["muonPhi12_den"], 1.0, 1.0, "B" ); 
*/
}

void MuonTriggerEfficiencyAnalyzer::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {
  std::cout<<"beginRun"<<std::endl;

  bool changed = true;
  if( hltConfig_.init(run, eventSetup, triggerProcess_, changed) ) {
  }
  else {
    std::cout << "Warning, didn't find process " << triggerProcess_.c_str() << std::endl;
    // Now crash
    assert(false);
  }

  triggerIndex_ = -1; 
  tagTriggerIndex_ = -1; 

  for(unsigned iHltPath=0; iHltPath<hltConfig_.size(); ++iHltPath) {
    std::string tempName = hltConfig_.triggerName(iHltPath);
    if(tempName.find(triggerName_) != std::string::npos) {
      std::cout<<tempName<<std::endl;
      triggerIndex_ = int(iHltPath);
    }
    if(tempName.find(tagTriggerName_) != std::string::npos) {
      tagTriggerIndex_ = int(iHltPath);
    }

    if( triggerIndex_>-1 && tagTriggerIndex_>-1 ) break; 
  } // end for each path

  if( triggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find trigger " <<  triggerName_.c_str() << std::endl;
    // Now crash
    assert(false);    
  }
  if( tagTriggerIndex_ == -1 ) {
    std::cout << "Warning, didn't find tag trigger " <<  tagTriggerName_.c_str() << std::endl;
    // Now crash
    assert(false);    
  }
}

void MuonTriggerEfficiencyAnalyzer::endRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void MuonTriggerEfficiencyAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &eventSetup) {
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  std::cout<<"Start event loop"<<std::endl;
  using reco::Muon;

  edm::Handle<reco::VertexCollection> pvHandle; 
  event.getByLabel(vertexes_, pvHandle);
  const reco::VertexCollection vertices = *pvHandle.product();
  unsigned int nGoodVtx = 0; 
  for(reco::VertexCollection::const_iterator it=vertices.begin(); it!=vertices.end(); ++it) {
    if( it->ndof()>4                     && 
	(fabs(it->z())<=24.)             && 
	(fabs(it->position().rho())<=2.)   ) 
      nGoodVtx++;
  }
  if( nGoodVtx==0 ) return;
  const reco::Vertex & pv = vertices[0];

  // Handle to the muon collection
  edm::Handle<std::vector<Muon> > muons;
  event.getByLabel(muons_, muons);

  if( nMaxMuons_>0 && muons->size()>nMaxMuons_ ) return; 

  // Get trigger results
  edm::Handle<edm::TriggerResults> triggerResults;
  event.getByLabel(edm::InputTag("TriggerResults", "", triggerProcess_), triggerResults);

  if(!triggerResults.isValid()) {
    std::cout << "Trigger results not valid" << std::endl;
    return;
  } 

  if( !triggerResults->accept(tagTriggerIndex_) ) return; // there are no tags

  // Get trigger summary 
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  event.getByLabel(edm::InputTag("hltTriggerSummaryAOD", "", triggerProcess_), triggerEvent);

  if(!triggerEvent.isValid()) { 
    std::cout << "TriggerEvent not valid" << std::endl;
    return;
  }

  // Sanity check
  assert(triggerResults->size()==hltConfig_.size());

  // Get trigger objects from trigger summary
  const trigger::TriggerObjectCollection & toc = triggerEvent->getObjects();

  // Modules in tag trigger path
  const std::vector<std::string>& tagModuleLabels(hltConfig_.moduleLabels(tagTriggerIndex_));
  assert( tagModuleLabels.size()==hltConfig_.size(tagTriggerIndex_) );
  const unsigned int tagModuleIndex( hltConfig_.size(tagTriggerIndex_)-2 ); // index of last filter (excluding HLTEndBool)
  const unsigned int tagFilterIndex( triggerEvent->filterIndex( edm::InputTag( tagModuleLabels[tagModuleIndex], "", triggerProcess_) ) );
  assert( tagFilterIndex < triggerEvent->sizeFilters() );
  const trigger::Vids & tagVids( triggerEvent->filterIds(tagFilterIndex) );
  const trigger::Keys & tagKeys( triggerEvent->filterKeys(tagFilterIndex) );
  assert( tagVids.size()==tagKeys.size() );
  const unsigned int nTagTrig(tagVids.size());

  // Loop muon collection and fill histograms
  for(std::vector<Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1) {
    if( tagId_ == "Tight" && !muon::isTightMuon( (*mu1), pv ) )
      continue; 
    else if( tagId_ == "Medium" && !muon::isMediumMuon( (*mu1) ) )
      continue;
    else if( tagId_ == "Loose" && !muon::isLooseMuon( (*mu1) ) )
      continue;
    if( (*mu1).pt()>tagPtCut_ && fabs((*mu1).eta())<tagEtaCut_ && isIsolated(*mu1, isolationType_, isolationCut_) ) {
      // Is this a suitable tag? 
      double maxTagDeltaR = 0.3; 
      double finTagDeltaR = 10.0; 
      bool isTagTrigMatch = false; 
      for(unsigned int i=0; i!=nTagTrig; ++i) {
	const trigger::TriggerObject & tagTo = toc[tagKeys[i]];
	double tmpTagDeltaR = deltaR( (*mu1), tagTo ); 
	if( tmpTagDeltaR<finTagDeltaR ) {
	  finTagDeltaR = tmpTagDeltaR;
	  if( tmpTagDeltaR<maxTagDeltaR ) {
	    isTagTrigMatch = true;
	    //break;
	  }
	}
      }

      hists_["deltaR_trobj_tag"]->Fill(finTagDeltaR); 

      if(isTagTrigMatch==false) continue; 

      // Go on and look for a probe
      for(std::vector<Muon>::const_iterator mu2=muons->begin(); mu2!=muons->end(); ++mu2) {
	if( mu2==mu1 ) continue; 

        if( probeId_ == "Tight" && !muon::isTightMuon( (*mu2), pv ) )
          continue; 
        else if( probeId_ == "Medium" && !muon::isMediumMuon( (*mu2) ) )
          continue;
        else if( probeId_ == "Loose" && !muon::isLooseMuon( (*mu2) ) )
          continue;
	if( /*(*mu2).pt()>probePtCut_ &&*/ fabs((*mu2).eta())<probeEtaCut_ && isIsolated(*mu2, isolationType_, isolationCut_) ) {
	  if( mu1->charge()*mu2->charge()<0 ) { // check only muon pairs of unequal charge 

	    double mumuMass = (mu1->p4()+mu2->p4()).mass();
	    hists_["mumuMass_all"]->Fill( mumuMass );

	    if( mumuMass>86. && mumuMass<96. ) { // check only muon pairs compatible with Z decay
	      // Ok, probe candidate found
	      // Check if more requirements apply (i.e. probe needs be matched to some trigger filter)

	      // Go on and fill denominator plots
	      hists_["muonPt_tag" ]->Fill( mu1->pt () );
	      hists_["muonEta_tag"]->Fill( mu1->eta() );
	      hists_["muonAbsEta_tag"]->Fill( fabs(mu1->eta()) );
	      hists_["muonPhi_tag"]->Fill( mu1->phi() );
	      hists_["muonNvtx_tag"]->Fill( nGoodVtx );

	      hists_["muonPt_probe_den" ]->Fill( mu2->pt () );
              if((*mu2).pt()>probePtCut_) {
	        hists_["muonEta_probe_den"]->Fill( mu2->eta() );
	        hists_["muonAbsEta_probe_den"]->Fill( fabs(mu2->eta()) );
	        hists_["muonPhi_probe_den"]->Fill( mu2->phi() );
	        hists_["muonNvtx_probe_den"]->Fill( nGoodVtx );
              }
	      hists_["mumuMass_den"]->Fill( mumuMass );

	      hists_["muonPt12_den"]->Fill( mu1->pt(), mu2->pt() );
	      hists_["muonEta12_den"]->Fill( mu1->eta(), mu2->eta() );
	      hists_["muonPhi12_den"]->Fill( mu1->phi(), mu2->phi() );



	      // Modules in probe trigger path
	      const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_));
	      const unsigned int m(hltConfig_.size(triggerIndex_));
	      assert( moduleLabels.size()==m );
	      const unsigned int lastModuleIndex(triggerResults->index(triggerIndex_));
	      assert(lastModuleIndex<m);

  
	      // Check if probe passes
	      bool isPassingProbeL1 = false; 
	      bool isPassingProbeL2 = false; 
	      bool isPassingProbeL3 = false; 
	      bool isPassingProbeL3OIStateSeed = false; 
	      bool isPassingProbeL3OIHitSeed = false; 
	      bool isPassingProbeL3IOHitSeed = false; 
	      bool isPassingProbeL3TkTrack = false; 
	      //bool isPassingProbeL3OIStateTkTrack = false; 
	      //bool isPassingProbeL3OIHitTkTrack = false; 
	      //bool isPassingProbeL3IOHitTkTrack = false; 
	      bool isPassingProbeL3Global = false; 
	      bool isPassingProbeL3Iso = false; 
	      //bool isPassingProbeL3OIStateGlobal = false; 
	      //bool isPassingProbeL3OIHitGlobal = false; 
	      //bool isPassingProbeL3IOHitGlobal = false; 
	      unsigned int filterIndex = triggerEvent->sizeFilters(); 
	      double finProbeDeltaR = 999999.; 
	      double maxProbeDeltaR = 999999.;
	      double muEta = 999999.;
	      double muPhi = 999999.;
              //double l2pt = -999999.;

	      /*if( probeFilterNum_.length()==0 ) { // full path efficiency
		if( triggerResults->accept(triggerIndex_) ) {
		  const unsigned int moduleIndex( hltConfig_.size(triggerIndex_)-2 ); // index of last filter (excluding HLTEndBool)
		  filterIndex = triggerEvent->filterIndex( edm::InputTag( moduleLabels[moduleIndex], "", triggerProcess_) ); 
		}
		}*/
	      // efficiency of an intermediate filter
	      // Results from TriggerEvent product - Attention: must look only for
	      // modules actually run in this path for this event!
	      for(unsigned int j=0; j<=lastModuleIndex; ++j) { 
		if( probeFilterL1_.compare(moduleLabels[j])!=0 ) continue;
                std::cout<<"L1 module: "<<moduleLabels[j]<<std::endl; 
		filterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterL1_, "", triggerProcess_));
		std::cout<<"Debug"<<std::endl;
                break; 
	      }
std::cout<<"L1"<<std::endl;
	      //L1 matching begins
	      if( filterIndex<triggerEvent->sizeFilters() ) {
		const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
		const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
		assert( vids.size()==keys.size() );
		const unsigned int nProbeTrig(vids.size());
		const std::string  moduleType(hltConfig_.moduleType(probeFilterL1_));
		std::cout<<"Debug2"<<std::endl;
		if( moduleType.compare("HLTLevel1GTSeed")==0 ) { // L1 matching
		  // Propagation to MB2/ME2 for L1 matching
		  eventSetup.get<IdealMagneticFieldRecord>().get(magneticField_);
		  eventSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", propagator_);
		  eventSetup.get<MuonRecoGeometryRecord>().get(detLayerGeometry_);
		  GlobalPoint pos( mu2->innerTrack()->outerPosition().x(), 
				   mu2->innerTrack()->outerPosition().y(), 
				   mu2->innerTrack()->outerPosition().z() );
		  GlobalVector mom( mu2->innerTrack()->outerMomentum().x(), 
				    mu2->innerTrack()->outerMomentum().y(), 
				    mu2->innerTrack()->outerMomentum().z() );
		  const FreeTrajectoryState state(pos, mom, (*mu2).charge(), &*magneticField_);
		  const DetLayer *detLayer = NULL;
		  if( fabs(mu2->eta())>0.9 ) {
		    if( mu2->eta()>0. )
		      detLayer = detLayerGeometry_->idToLayer(CSCDetId(1, 2, 0, 0, 0));
		    else 
		      detLayer = detLayerGeometry_->idToLayer(CSCDetId(2, 2, 0, 0, 0));
		    if(detLayer==NULL && fabs(mu2->eta())<1.2) 
		      detLayer = detLayerGeometry_->idToLayer(DTChamberId(0, 2, 0));
		  }
		  else 
		    detLayer = detLayerGeometry_->idToLayer(DTChamberId(0, 2, 0)); 
		  if(detLayer==NULL) {
		    std::cout << "WARNING: detLayer is NULL at eta=" << mu2->eta() << std::endl;
		    break; 
		  }
		  TrajectoryStateOnSurface tsos = propagator_->propagate(state, detLayer->surface());
		  if(!tsos.isValid()) {
		    std::cout << "WARNING: tsos is not valid at eta=" << mu2->eta() << std::endl;
		    break; 
		  }
                  std::cout<<"Debug3"<<std::endl;
		  muEta = tsos.globalPosition().eta(); 
		  muPhi = tsos.globalPosition().phi(); 
		  maxProbeDeltaR = 0.3;
		}
		else {
		  std::cout<<"Wrong L1 path name"<<std::endl;
		  break;
		}
		for(unsigned int i=0; i<nProbeTrig; ++i) {
		  const trigger::TriggerObject & to = toc[keys[i]];
		  double tmpProbeDeltaR = deltaR( muEta, muPhi, to.eta(), to.phi() ); 
		  if( tmpProbeDeltaR<finProbeDeltaR ) {
		    finProbeDeltaR = tmpProbeDeltaR;
		    if( tmpProbeDeltaR<maxProbeDeltaR ) {
		      isPassingProbeL1 = true;
		      //break;
		    }
		  }
		}
	      }//L1 matching ends
	      hists_["deltaR_trobj_probe_L1"]->Fill(finProbeDeltaR);
	      
	      if(isPassingProbeL1==false) continue;

	      
	      hists_["muonPt_probe_L1" ]->Fill( mu2->pt () );
              if((*mu2).pt()>probePtCut_) {
       	        hists_["muonEta_probe_L1"]->Fill( mu2->eta() );
	        hists_["muonAbsEta_probe_L1"]->Fill( fabs(mu2->eta()) );
	        hists_["muonPhi_probe_L1"]->Fill( mu2->phi() );
	        hists_["muonNvtx_probe_L1"]->Fill( nGoodVtx );
	      }

	      hists_["mumuMass_L1"]->Fill( mumuMass );

	      hists_["muonPt12_L1"]->Fill( mu1->pt(), mu2->pt() );
	      hists_["muonEta12_L1"]->Fill( mu1->eta(), mu2->eta() );
	      hists_["muonPhi12_L1"]->Fill( mu1->phi(), mu2->phi() );

	      
	      finProbeDeltaR = 999999.; 
	     

	      for(unsigned int j=0; j<=lastModuleIndex; ++j) { 
		if( probeFilterL2_.compare(moduleLabels[j])!=0 ) continue;
                std::cout<<"L2 module:"<<moduleLabels[j]<<std::endl; 
		filterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterL2_, "", triggerProcess_));
		break; 
              }
 
std::cout<<"L2"<<std::endl;
	      //L2 matching begins
	      if( filterIndex<triggerEvent->sizeFilters() ) {
		const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
		const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
		assert( vids.size()==keys.size() );
		const unsigned int nProbeTrig(vids.size());
		const std::string  moduleType(hltConfig_.moduleType(probeFilterL2_));
		if( moduleType.compare("HLTMuonL2PreFilter")==0 ) { // L2 matching
		  muEta = (*mu2).eta(); 
		  muPhi = (*mu2).phi(); 
		  maxProbeDeltaR = 0.3;
		}
		else {
		  std::cout<<"Wrong L2 path name:"<<moduleType<<std::endl;
		  break;
		}
		  
		for(unsigned int i=0; i<nProbeTrig; ++i) {
		  const trigger::TriggerObject & to = toc[keys[i]];
		  double tmpProbeDeltaR = deltaR( muEta, muPhi, to.eta(), to.phi() ); 
		  if( tmpProbeDeltaR<finProbeDeltaR ) {
		    finProbeDeltaR = tmpProbeDeltaR;
		    if( tmpProbeDeltaR<maxProbeDeltaR ) {
		      isPassingProbeL2 = true;
                      //l2pt = to.pt();
		      //break;
		    }
		  }
		}
	      }//L2 matching ends
	      hists_["deltaR_trobj_probe_L2"]->Fill(finProbeDeltaR);
	      
	      if(isPassingProbeL2==false) continue;

	      hists_["muonPt_probe_L2" ]->Fill( mu2->pt () );
              if((*mu2).pt()>probePtCut_) {
	        hists_["muonEta_probe_L2"]->Fill( mu2->eta() );
	        hists_["muonAbsEta_probe_L2"]->Fill( fabs(mu2->eta()) );
	        hists_["muonPhi_probe_L2"]->Fill( mu2->phi() );
	        hists_["muonNvtx_probe_L2"]->Fill( nGoodVtx );
              }
	      hists_["mumuMass_L2"]->Fill( mumuMass );

	      hists_["muonPt12_L2"]->Fill( mu1->pt(), mu2->pt() );
	      hists_["muonEta12_L2"]->Fill( mu1->eta(), mu2->eta() );
	      hists_["muonPhi12_L2"]->Fill( mu1->phi(), mu2->phi() );

	      finProbeDeltaR = 999999.;
 
              if(useRerun_){
std::cout<<"L3 Seed"<<std::endl;
 	      //L3 seed matching begins
              edm::Handle<reco::TrackCollection> L3TkTrackOIStateHandle;
              edm::Handle<reco::TrackCollection> L3TkTrackOIHitHandle;
              edm::Handle<reco::TrackCollection> L3TkTrackIOHitHandle;
              edm::Handle<L3MuonTrajectorySeedCollection> L3SeedOIStateHandle;
              edm::Handle<L3MuonTrajectorySeedCollection> L3SeedOIHitHandle;
              edm::Handle<L3MuonTrajectorySeedCollection> L3SeedIOHitHandle;
	      edm::Handle<reco::TrackCollection> L3GlobalHandle;
              event.getByLabel("hltL3TrajSeedOIState",L3SeedOIStateHandle);
              event.getByLabel("hltL3TrajSeedOIHit",L3SeedOIHitHandle);
              event.getByLabel("hltL3TrajSeedIOHit",L3SeedIOHitHandle);
	      event.getByLabel("hltL3Muons", L3GlobalHandle);
	
	      muEta = (*mu2).eta(); 
              muPhi = (*mu2).phi();
              maxProbeDeltaR = 0.3;

std::cout<<"L3 Seed oistate for"<<std::endl;
              if(!L3SeedOIStateHandle.isValid()) continue;
              for(unsigned i=0; i!=L3SeedOIStateHandle->size(); i++) {
                const L3MuonTrajectorySeed seed = L3SeedOIStateHandle->at(i);
                double tmpProbeDeltaR = deltaR( muEta, muPhi, seed.l2Track()->eta(), seed.l2Track()->phi() );
                if( tmpProbeDeltaR<finProbeDeltaR ) {
                  finProbeDeltaR = tmpProbeDeltaR;
                  if( tmpProbeDeltaR<maxProbeDeltaR ) {
                    isPassingProbeL3OIStateSeed = true;
                    //break;
                  }
                }
              }
	      finProbeDeltaR = 999999.;
	      
             
	      if(!L3SeedOIHitHandle.isValid()) continue;
	      for(unsigned i=0; i!=L3SeedOIHitHandle->size(); i++) {
		const L3MuonTrajectorySeed seed = L3SeedOIHitHandle->at(i);
		double tmpProbeDeltaR = deltaR( muEta, muPhi, seed.l2Track()->eta(), seed.l2Track()->phi() );
		if( tmpProbeDeltaR<finProbeDeltaR ) {
		  finProbeDeltaR = tmpProbeDeltaR;
		  if( tmpProbeDeltaR<maxProbeDeltaR ) {
		    isPassingProbeL3OIHitSeed = true;
		    //break;
		  }
		}
	      }	      
	      finProbeDeltaR = 999999.;
                
	      if(!L3SeedIOHitHandle.isValid()) continue;
	      for(unsigned i=0; i!=L3SeedIOHitHandle->size(); i++) {
		const L3MuonTrajectorySeed seed = L3SeedIOHitHandle->at(i);
		double tmpProbeDeltaR = deltaR( muEta, muPhi, seed.l2Track()->eta(), seed.l2Track()->phi() );
		if( tmpProbeDeltaR<finProbeDeltaR ) {
		  finProbeDeltaR = tmpProbeDeltaR;
		  if( tmpProbeDeltaR<maxProbeDeltaR ) {
		    isPassingProbeL3IOHitSeed = true;
		    //break;
		  }
		}
	      }

	      if(!isPassingProbeL3OIStateSeed&&!isPassingProbeL3OIHitSeed&&!isPassingProbeL3IOHitSeed) continue;
	    
	      hists_["muonPt_probe_L3Seed" ]->Fill( mu2->pt () );
              if( (*mu2).pt()>probePtCut_  ) {
	        hists_["muonEta_probe_L3Seed"]->Fill( mu2->eta() );
	        hists_["muonAbsEta_probe_L3Seed"]->Fill( fabs(mu2->eta()) );
	        hists_["muonPhi_probe_L3Seed"]->Fill( mu2->phi() );
	        hists_["muonNvtx_probe_L3Seed"]->Fill( nGoodVtx );
              }
	      finProbeDeltaR = 999999.;

	      std::cout<<"L3 track"<<std::endl;
              //L3 tk track matching begins
 
	      event.getByLabel("hltL3TkTracksFromL2OIState", L3TkTrackOIStateHandle);
	      event.getByLabel("hltL3TkTracksFromL2OIHit", L3TkTrackOIHitHandle);
	      event.getByLabel("hltL3TkTracksFromL2IOHit", L3TkTrackIOHitHandle);
              muEta = (*mu2).eta();
              muPhi = (*mu2).phi();
              maxProbeDeltaR = 0.3;
	      
              for(unsigned i=0; i!=L3TkTrackOIStateHandle->size(); i++) {
                const reco::Track tk=L3TkTrackOIStateHandle->at(i);
                double tmpProbeDeltaR = deltaR( muEta, muPhi, tk.eta(), tk.phi() );
                if( tmpProbeDeltaR<finProbeDeltaR ) {
                  finProbeDeltaR = tmpProbeDeltaR;
                  if( tmpProbeDeltaR<maxProbeDeltaR ) {
                    isPassingProbeL3TkTrack = true;
                    //break;
                  }
                }
              }
              if(!isPassingProbeL3TkTrack) {
	        for(unsigned i=0; i!=L3TkTrackOIHitHandle->size(); i++) {
                  const reco::Track tk=L3TkTrackOIHitHandle->at(i);
                  double tmpProbeDeltaR = deltaR( muEta, muPhi, tk.eta(), tk.phi() );
                  if( tmpProbeDeltaR<finProbeDeltaR ) {
                    finProbeDeltaR = tmpProbeDeltaR;
                    if( tmpProbeDeltaR<maxProbeDeltaR ) {
                      isPassingProbeL3TkTrack = true;
                      //break;
                    }
                  }
                }
              }
              if(!isPassingProbeL3TkTrack) {
	        for(unsigned i=0; i!=L3TkTrackIOHitHandle->size(); i++) {
                  const reco::Track tk=L3TkTrackIOHitHandle->at(i);
                  double tmpProbeDeltaR = deltaR( muEta, muPhi, tk.eta(), tk.phi() );
                  if( tmpProbeDeltaR<finProbeDeltaR ) {
                    finProbeDeltaR = tmpProbeDeltaR;
                    if( tmpProbeDeltaR<maxProbeDeltaR ) {
                      isPassingProbeL3TkTrack = true;
                      //break;
                    }
                  }
                }
	      }
              if(!isPassingProbeL3TkTrack) continue;
	      
              hists_["muonPt_probe_L3TkTrack" ]->Fill( mu2->pt () );
              if((*mu2).pt()>probePtCut_ ) {
	        hists_["muonEta_probe_L3TkTrack"]->Fill( mu2->eta() );
	        hists_["muonAbsEta_probe_L3TkTrack"]->Fill( fabs(mu2->eta()) );
	        hists_["muonPhi_probe_L3TkTrack"]->Fill( mu2->phi() );
	        hists_["muonNvtx_probe_L3TkTrack"]->Fill( nGoodVtx );
              }
	      //hists_["mumuMass_L3TkTrack"]->Fill( mumuMass );

	      //hists_["muonPt12_L3TkTrack"]->Fill( mu1->pt(), mu2->pt() );
	      //hists_["muonEta12_L3TkTrack"]->Fill( mu1->eta(), mu2->eta() );
	      //hists_["muonPhi12_L3TkTrack"]->Fill( mu1->phi(), mu2->phi() );

	      finProbeDeltaR = 999999.;

std::cout<<"L3 global matching"<<std::endl;
              //L3 global matching begins
              muEta = (*mu2).eta();
              muPhi = (*mu2).phi();
              maxProbeDeltaR = 0.3;
	      
              for(unsigned i=0; i!=L3GlobalHandle->size(); i++) {
                const reco::Track l3=L3GlobalHandle->at(i); 
                double tmpProbeDeltaR = deltaR( muEta, muPhi, l3.eta(), l3.phi() );
                if( tmpProbeDeltaR<finProbeDeltaR ) {
                  finProbeDeltaR = tmpProbeDeltaR;
                  if( tmpProbeDeltaR<maxProbeDeltaR ) {
                    isPassingProbeL3Global = true;
                    //break;
                  }
                }
              }
	      
              if(!isPassingProbeL3Global) continue;             
	      
              hists_["muonPt_probe_L3Global" ]->Fill( mu2->pt () );
              if((*mu2).pt()>probePtCut_ ) {
  	        hists_["muonEta_probe_L3Global"]->Fill( mu2->eta() );
	        hists_["muonAbsEta_probe_L3Global"]->Fill( fabs(mu2->eta()) );
	        hists_["muonPhi_probe_L3Global"]->Fill( mu2->phi() );
	        hists_["muonNvtx_probe_L3Global"]->Fill( nGoodVtx );
              }
	      finProbeDeltaR = 999999.;	      

              }
std::cout<<"L3"<<std::endl;
	      //L3 matching begins
	      for(unsigned int j=0; j<=lastModuleIndex; ++j) { 
		if( probeFilterL3_.compare(moduleLabels[j])!=0 ) continue; 
		std::cout<<"L3 module:"<<moduleLabels[j]<<"("<<probeFilterL3_<<")"<<std::endl;
                filterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterL3_, "", triggerProcess_));
		break; 
              }
	      if( filterIndex<triggerEvent->sizeFilters() ) {
		const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
		const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
		assert( vids.size()==keys.size() );
		const unsigned int nProbeTrig(vids.size());
		const std::string  moduleType(hltConfig_.moduleType(probeFilterL3_));

//		if( moduleType.compare("HLTMuonL3PreFilter")==0 ) { // L3 et al. matching
		  muEta = (*mu2).eta(); 
		  muPhi = (*mu2).phi(); 
		  maxProbeDeltaR = 0.3;
//		}
//		else {
//		  std::cout<<"Wrong L3 path name:"<<moduleType<<std::endl;
//		  break;
//		}

		for(unsigned int i=0; i<nProbeTrig; ++i) {
		  const trigger::TriggerObject & to = toc[keys[i]];
		  double tmpProbeDeltaR = deltaR( muEta, muPhi, to.eta(), to.phi() ); 
		  if( tmpProbeDeltaR<finProbeDeltaR ) {
		    finProbeDeltaR = tmpProbeDeltaR;
		    if( tmpProbeDeltaR<maxProbeDeltaR ) {
		      isPassingProbeL3 = true;
		      //break;
		    }
		  }
		}
	      } //L3 matching ends
	      
	      hists_["deltaR_trobj_probe_L3"]->Fill(finProbeDeltaR); 

	      if(!isPassingProbeL3) continue; 

	      // Ok, probe passed
	      // Go on and fill L3erator plots
	      hists_["muonPt_probe_L3" ]->Fill( mu2->pt () );
              if((*mu2).pt()>probePtCut_ ) {
	        hists_["muonEta_probe_L3"]->Fill( mu2->eta() );
	        hists_["muonAbsEta_probe_L3"]->Fill( fabs(mu2->eta()) );
	        hists_["muonPhi_probe_L3"]->Fill( mu2->phi() );
	        hists_["muonNvtx_probe_L3"]->Fill( nGoodVtx );
              }
	      hists_["mumuMass_L3"]->Fill( mumuMass );

	      hists_["muonPt12_L3"]->Fill( mu1->pt(), mu2->pt() );
	      hists_["muonEta12_L3"]->Fill( mu1->eta(), mu2->eta() );
	      hists_["muonPhi12_L3"]->Fill( mu1->phi(), mu2->phi() );

              if(useIso_) {
	        finProbeDeltaR = 999999.;	      
                filterIndex = triggerEvent->filterIndex(edm::InputTag(probeFilterL3Iso_, "", triggerProcess_));
	        if( filterIndex<triggerEvent->sizeFilters() ) {
                  const trigger::Vids & vids( triggerEvent->filterIds(filterIndex) );
                  const trigger::Keys & keys( triggerEvent->filterKeys(filterIndex) );
                  assert( vids.size()==keys.size() );
                  const unsigned int nProbeTrig(vids.size());
                  for(unsigned int i=0; i<nProbeTrig; ++i) {
                    const trigger::TriggerObject & to = toc[keys[i]];
                    double tmpProbeDeltaR = deltaR( muEta, muPhi, to.eta(), to.phi() );
                    if( tmpProbeDeltaR<finProbeDeltaR ) {
                      finProbeDeltaR = tmpProbeDeltaR;
                      if( tmpProbeDeltaR<maxProbeDeltaR ) {
                        isPassingProbeL3Iso = true;
                      }
                    }
                  }
                }
                if(!isPassingProbeL3Iso) continue;
                hists_["muonPt_probe_L3Iso" ]->Fill( mu2->pt () );
                if((*mu2).pt()>probePtCut_ ) {
                  hists_["muonEta_probe_L3Iso"]->Fill( mu2->eta() );
                  hists_["muonAbsEta_probe_L3Iso"]->Fill( fabs(mu2->eta()) );
                  hists_["muonPhi_probe_L3Iso"]->Fill( mu2->phi() );
                  hists_["muonNvtx_probe_L3Iso"]->Fill( nGoodVtx );
                }
                
              }
	    } 
	  }
	}
      }
    }
  }
}

bool MuonTriggerEfficiencyAnalyzer::isIsolated(reco::Muon muon, std::string isolationType, double isolationCut) {
  double isolationVariable=.0;
  if( isolationType=="TrkIso" )
    isolationVariable = muon.isolationR03().sumPt/muon.pt();
  else if( isolationType=="PFIso" )
    isolationVariable = (muon.pfIsolationR04().sumChargedHadronPt + muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt)/muon.pt();
  else //Wrong isolation type
    std::cout<<"Wrong isolation type"<<std::endl;

  if( isolationVariable<isolationCut )
    return true;
  else
    return false;
}


// define this as a plug-in
DEFINE_FWK_MODULE(MuonTriggerEfficiencyAnalyzer);
