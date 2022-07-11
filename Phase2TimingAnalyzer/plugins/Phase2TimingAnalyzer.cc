// -*- C++ -*-
//
// Package:    Phase2TimingAnalyzer/Phase2TimingAnalyzer
// Class:      Phase2TimingAnalyzer
//
/**\class Phase2TimingAnalyzer Phase2TimingAnalyzer.cc Phase2TimingAnalyzer/Phase2TimingAnalyzer/plugins/Phase2TimingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthew Citron <mcitron@ucsb.edu> 10/19/2017
//         Created:  Tue, 30 Nov 2021 21:39:01 GMT
//
//

// system include files
#include <memory>
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "Phase2Timing/Phase2TimingAnalyzer/plugins/JetTimingTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

//
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "TMath.h"
#include "TTree.h"

#include "TGeoPolygon.h" //may or may not need! for ctau addition
//include electron and photon data format



//
// class declaration
//
struct tree_struc_{ //structs group several related variables, unlike an array it
  int nrecojets;    //it doesnt have to the same data type
  int ngen;
  int ntrack;
  int nLLP;
  //ctau beginning might need some of the header files and other files to do calculation
  std::vector<float> e_ctau;
  std::vector<float> e_eta;
  std::vector<float> e_phi;
  std::vector<float> e_pt;
  std::vector<float> e_vx;
  std::vector<float> e_vy;
  std::vector<float> e_vz;
  std::vector<float> e_ebeta;
  std::vector<float> e_ebphi;
  std::vector<float> e_ebdelay;
  std::vector<float> e_hgeta;
  std::vector<float> e_hgphi;
  std::vector<float> e_hgdelay;

  std::vector<float> LLP_eta;
  std::vector<float> LLP_phi;
  std::vector<float> LLP_pt;
  std::vector<float> LLP_mass;

  std::vector<float> track_eta;
  std::vector<float> track_phi;
  std::vector<float> track_pt;
  std::vector<float> track_vx;
  std::vector<float> track_vy;
  std::vector<float> track_vz;

  std::vector<float> recojet_pt;
  std::vector<float> recojet_eta;
  std::vector<float> recojet_phi;
  std::vector<float> recojet_e;
  std::vector<float> recojet_ECALtime;
  std::vector<float> recojet_ECALenergy;
  std::vector<float> recojet_ECALnCells;
  std::vector<float> recojet_MTDtime;
  std::vector<float> recojet_MTDenergy;
  std::vector<float> recojet_MTDnCells;
  std::vector<float> recojet_MTDClutime;
  std::vector<float> recojet_MTDCluenergy;
  std::vector<float> recojet_MTDnClus;
  std::vector<float> recojet_HGCALtime;
  std::vector<float> recojet_HGCALenergy;
  std::vector<float> recojet_HGCALnTracksters;
  std::vector<float> recojet_closestGenIndex;
  std::vector<float> recojet_closestGenR;
  std::vector<float> recojet_closestEbGenIndex;
  std::vector<float> recojet_closestEbGenR;
  std::vector<float> recojet_closestHgGenIndex;
  std::vector<float> recojet_closestHgGenR;
};

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class Phase2TimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Phase2TimingAnalyzer(const edm::ParameterSet&); 
  ~Phase2TimingAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private: 
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void initTreeStructure();
  void clearVectors();
  JetTimingTools _jetTimingTools;

  // ---------- member data -------------------- //
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexCollectionToken_;
  edm::EDGetTokenT<std::vector<reco::Photon>> photonCollectionToken_; //here is token for Photon
  edm::EDGetTokenT<std::vector<reco::Track>> trackCollectionToken_;

  edm::Service<TFileService> fs;
  
  const edm::EDGetTokenT< edm::View<reco::GenParticle> > _genParticles; 
  edm::Handle< edm::View<reco::GenParticle> > _genParticlesH;
  const edm::EDGetTokenT< edm::View<reco::PFJet> > _recoak4PFJets; 
  edm::Handle< edm::View<reco::PFJet> > _recoak4PFJetsH;
  const edm::EDGetTokenT<edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>> ecalRecHitsEBToken_;
  edm::Handle< edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>> > _ecalRecHitsEBH;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterEMToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterMergeToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterHADToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterTrkEMToken;
  const edm::EDGetTokenT<edm::View<ticl::Trackster>> _tracksterTrkToken;
  const edm::EDGetTokenT<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>> mtdRecHitsBTLToken_;
  edm::Handle< edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>> > _mtdRecHitsBTLH;
  const edm::EDGetTokenT<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>> mtdRecHitsETLToken_;
  edm::Handle< edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>> > _mtdRecHitsETLH;
  const edm::EDGetTokenT<FTLClusterCollection> btlRecCluToken_;
  edm::Handle<FTLClusterCollection> _btlRecCluH;
  const edm::EDGetTokenT<FTLClusterCollection> etlRecCluToken_;
  edm::Handle<FTLClusterCollection> _etlRecCluH;  

  // setup tree;                                                                                                                                             
  TTree* tree;
  tree_struc_ tree_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//if I want to add anothe collection with info I need this is wehre I would add it
// constructors and destructor
//
Phase2TimingAnalyzer::Phase2TimingAnalyzer(const edm::ParameterSet& iConfig):

  _jetTimingTools(consumesCollector()),

  vertexCollectionToken_(consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"))),
  photonCollectionToken_(consumes<std::vector<reco::Photon>>(edm::InputTag("photons"))), //good one
  trackCollectionToken_(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"))),
  _genParticles(consumes< edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  _genParticlesH(),
  _recoak4PFJets(consumes< edm::View<reco::PFJet> >(iConfig.getParameter<edm::InputTag>("recoak4PFJets"))),
  _recoak4PFJetsH(),
  ecalRecHitsEBToken_{consumes<edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>>(												       iConfig.getParameter<edm::InputTag>("ebRecHitsColl"))},
  _ecalRecHitsEBH(),
  _tracksterEMToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersEM"))),
  _tracksterMergeToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersMerge"))),
  _tracksterHADToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersHAD"))),
  _tracksterTrkEMToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersTrkEM"))),
  _tracksterTrkToken(consumes<edm::View<ticl::Trackster>>(iConfig.getParameter<edm::InputTag>("ticlTrackstersTrk"))),
  mtdRecHitsBTLToken_{consumes<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>>(												       iConfig.getParameter<edm::InputTag>("mtdBTLRecHitsColl"))},
  _mtdRecHitsBTLH(),
  mtdRecHitsETLToken_{consumes<edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>>(												       iConfig.getParameter<edm::InputTag>("mtdETLRecHitsColl"))},
  _mtdRecHitsETLH(),
  btlRecCluToken_(consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recBTLCluTag"))),
  _btlRecCluH(),
  etlRecCluToken_(consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recETLCluTag"))),
  _etlRecCluH()
{
}

Phase2TimingAnalyzer::~Phase2TimingAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void Phase2TimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  _jetTimingTools.init(iSetup);
  Handle< std::vector<reco::Vertex> > vertexCollectionH;
  Handle< std::vector<reco::Track> > trackCollectionH;
  Handle<View<ticl::Trackster>> tracksterEMH;
  Handle<View<ticl::Trackster>> tracksterMergeH;
  Handle<View<ticl::Trackster>> tracksterHADH;
  Handle<View<ticl::Trackster>> tracksterTrkEMH;
  Handle<View<ticl::Trackster>> tracksterTrkH;

  Handle< std::vector<reco::Photon> > photonCollectionTokenH; //photon, good one

  iEvent.getByToken(vertexCollectionToken_,vertexCollectionH);
  iEvent.getByToken(photonCollectionToken_,photonCollectionTokenH); //good one
  iEvent.getByToken(trackCollectionToken_,trackCollectionH);
  iEvent.getByToken(_genParticles, _genParticlesH);
  iEvent.getByToken(_recoak4PFJets, _recoak4PFJetsH);
  iEvent.getByToken(ecalRecHitsEBToken_, _ecalRecHitsEBH);
  iEvent.getByToken(_tracksterEMToken, tracksterEMH);
  iEvent.getByToken(_tracksterMergeToken, tracksterMergeH);
  iEvent.getByToken(_tracksterHADToken, tracksterHADH);
  iEvent.getByToken(_tracksterTrkEMToken, tracksterTrkEMH);
  iEvent.getByToken(_tracksterTrkToken, tracksterTrkH);
  iEvent.getByToken(mtdRecHitsBTLToken_, _mtdRecHitsBTLH);
  iEvent.getByToken(mtdRecHitsETLToken_, _mtdRecHitsETLH);
  
  iEvent.getByToken(btlRecCluToken_,_btlRecCluH);
  iEvent.getByToken(etlRecCluToken_,_etlRecCluH);


  //  auto _btlRecCluH = makeValid(iEvent.getHandle(btlRecCluToken_));
  // auto _etlRecCluH = makeValid(iEvent.getHandle(etlRecCluToken_));

  //  auto const& ecalRecHitsEB = iEvent.get(ecalRecHitsEBToken_);

  //variable declaration
  int nelectron = 0;
  int nrecojets = 0;
  int ngen = 0;
  int ntrack = 0;
  int nLLP = 0; //counter for LLP
  //ctau spot 2
  std::vector<float> e_ctau;
  std::vector<float> e_eta;
  std::vector<float> e_phi;
  std::vector<float> e_pt;
  std::vector<float> e_vx;
  std::vector<float> e_vy;
  std::vector<float> e_vz;
  std::vector<float> e_ebeta;
  std::vector<float> e_ebphi;
  std::vector<float> e_ebdelay;
  std::vector<float> e_hgeta;
  std::vector<float> e_hgphi;
  std::vector<float> e_hgdelay;

  std::vector<float> LLP_eta; //writing out the LLP pt, eta and phi
  std::vector<float> LLP_phi;
  std::vector<float> LLP_pt;
  std::vector<float> LLP_mass;

  std::vector<float> track_eta;
  std::vector<float> track_phi;
  std::vector<float> track_pt;
  std::vector<float> track_vx;
  std::vector<float> track_vy;
  std::vector<float> track_vz;

  std::vector<float> recojet_pt;
  std::vector<float> recojet_eta;
  std::vector<float> recojet_phi;
  std::vector<float> recojet_e;
  std::vector<float> recojet_ECALtime;
  std::vector<float> recojet_ECALenergy;
  std::vector<float> recojet_ECALnCells;
  std::vector<float> recojet_MTDtime;
  std::vector<float> recojet_MTDenergy;
  std::vector<float> recojet_MTDnCells;
  std::vector<float> recojet_MTDClutime;
  std::vector<float> recojet_MTDCluenergy;
  std::vector<float> recojet_MTDnClus;
  std::vector<float> recojet_HGCALtime;
  std::vector<float> recojet_HGCALenergy;
  std::vector<float> recojet_HGCALnTracksters;
  std::vector<float> recojet_closestGenIndex;
  std::vector<float> recojet_closestGenR;
  std::vector<float> recojet_closestEbGenIndex;
  std::vector<float> recojet_closestEbGenR;
  std::vector<float> recojet_closestHgGenIndex;
  std::vector<float> recojet_closestHgGenR;

  
  bool debug=0;


  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON GENPARTICLES --------------------------------------"<<std::endl; 
  for (const auto & genpar_iter : *_genParticlesH){

    if (genpar_iter.mother(0) == NULL)continue;
    //genparticle is the ground truth infomation we put into the software, physics of the paticles

    reco::GenParticle * genParticleMother = (reco::GenParticle *) genpar_iter.mother();
    std::vector<double> ecalIntersection = _jetTimingTools.surfaceIntersection(genpar_iter,*genParticleMother,130);
    std::vector<double> hgcalIntersection = _jetTimingTools.endCapIntersection(genpar_iter,*genParticleMother,300,520);

    if(abs(genpar_iter.pdgId()) !=11  || genParticleMother->pdgId()!=6000113) continue; //pdgID for electron w mother LLP
    float vx = genpar_iter.vertex().x();
    float vy = genpar_iter.vertex().y();
    float vz = genpar_iter.vertex().z();
    nelectron++; //not used anywhere???
    ngen++;
    
    //for ctau calculation
    double displacement = TMath::Sqrt((genpar_iter.vertex()-genParticleMother->vertex()).Mag2());
    double genParticleBeta = genParticleMother->p()/genParticleMother->energy();
    double genParticleGamma = 1./TMath::Sqrt(1.-genParticleBeta*genParticleBeta);
    double ctau = displacement*10 / (genParticleBeta*genParticleGamma);
    e_ctau.push_back(ctau); //ctau distance from beam to decay
    
    e_pt.push_back(genpar_iter.pt());
    e_eta.push_back(genpar_iter.eta());
    e_phi.push_back(genpar_iter.phi());
    e_vx.push_back(vx);
    e_vy.push_back(vy);
    e_vz.push_back(vz);
    e_ebphi.push_back(ecalIntersection[1]);
    e_ebeta.push_back(ecalIntersection[0]);
    e_ebdelay.push_back(ecalIntersection[3]);
    e_hgphi.push_back(hgcalIntersection[1]);
    e_hgeta.push_back(hgcalIntersection[0]);
    e_hgdelay.push_back(hgcalIntersection[3]);
  }

  auto const& ecalRecHitsEB = iEvent.get(ecalRecHitsEBToken_);
  auto const& mtdRecHitsBTL = iEvent.get(mtdRecHitsBTLToken_);
  auto const& mtdRecHitsETL = iEvent.get(mtdRecHitsETLToken_);
  //  auto const& mtdClusBTL = iEvent.get(btlRecCluToken_);
  //  auto const& mtdClusETL = iEvent.get(btlRecCluToken_);

  // Loop over all tracks in a certain collection

  for (const auto & track_iter : *trackCollectionH){ //loops over each track
      track_pt.push_back(track_iter.pt());
      track_eta.push_back(track_iter.eta());
      track_phi.push_back(track_iter.phi()); //see if there is .vertex to calc vertex for tracks
      
      track_vx.push_back(track_iter.vertex().x());
      track_vy.push_back(track_iter.vertex().y());
      track_vz.push_back(track_iter.vertex().z());
      ntrack++;
  }
  
  for (const auto & genpar_iter : *_genParticlesH){ //to get LLP info 
    if (genpar_iter.mother(0) == NULL)continue;
    //genparticle is the ground truth infomation we put into the software, physics of the paticle
    
    if(abs(genpar_iter.pdgId())!=6000113) continue; //pdgID for LLP w mother Higgs, || genParticleMother->pdgId()!=25
    LLP_pt.push_back(genpar_iter.pt());
    LLP_eta.push_back(genpar_iter.eta());
    LLP_phi.push_back(genpar_iter.phi());
    LLP_mass.push_back(genpar_iter.mass());
    nLLP++;
  }


  if(debug)std::cout<<" [DEBUG MODE] --------------- LOOP ON RECO JETS --------------------------------------"<<std::endl; 
  for (const auto & recojet_iter : *_recoak4PFJetsH){

    if(recojet_iter.pt()<20)continue;
    if(fabs(recojet_iter.eta())>3)continue;
    
    nrecojets++;
    recojet_pt.push_back(recojet_iter.pt());
    recojet_eta.push_back(recojet_iter.eta());
    recojet_phi.push_back(recojet_iter.phi()); //recojet_phiPhase2TimingAnalyzer.cc.push_back(recojet_iter.phi());
    recojet_e.push_back(recojet_iter.energy());
    int closestGenIndex = -1;
    float closestGenR = 999;
    int closestEbGenIndex = -1;
    float closestEbGenR = 999;
    int closestHgGenIndex = -1;
    float closestHgGenR = 999;
    for (int ig = 0; ig < ngen; ig++){
	TLorentzVector jetVec;
	jetVec.SetPtEtaPhiM(recojet_iter.pt(),recojet_iter.eta(),recojet_iter.phi(),0);
	TLorentzVector pos;
	pos.SetPtEtaPhiM(1,e_eta[ig],e_phi[ig],0);
	TLorentzVector ebPos;
	ebPos.SetPtEtaPhiM(1,e_ebeta[ig],e_ebphi[ig],0);
	TLorentzVector hgPos;
	hgPos.SetPtEtaPhiM(1,e_hgeta[ig],e_hgphi[ig],0);
	
	float ebDeltaR = jetVec.DeltaR(ebPos);
	if (ebDeltaR < 0.4 and closestEbGenR > ebDeltaR){
	    closestEbGenR = ebDeltaR;
	    closestEbGenIndex = ig;
	}
	float hgDeltaR = jetVec.DeltaR(hgPos);
	if (hgDeltaR < 0.4 and closestHgGenR > hgDeltaR){
	    closestHgGenR = hgDeltaR;
	    closestHgGenIndex = ig;
	}
	float genDeltaR = jetVec.DeltaR(pos);
	if (genDeltaR < 0.4 and closestGenR > genDeltaR){
	    closestGenR = genDeltaR;
	    closestGenIndex = ig;
	}
    }
    recojet_closestGenIndex.push_back(closestGenIndex);
    recojet_closestGenR.push_back(closestGenR);
    recojet_closestEbGenIndex.push_back(closestEbGenIndex);
    recojet_closestEbGenR.push_back(closestEbGenR);
    recojet_closestHgGenIndex.push_back(closestHgGenIndex);
    recojet_closestHgGenR.push_back(closestHgGenR);
    

    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM ECAL --------------------------------------"<<std::endl; 
    float weightedECALTimeCell = 0;
    float totalECALEnergyCell = 0;
    unsigned int ECALnCells = 0;
    if(fabs(recojet_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromEcalCells(recojet_iter, ecalRecHitsEB, weightedECALTimeCell, totalECALEnergyCell, ECALnCells);
    recojet_ECALenergy.push_back(totalECALEnergyCell);
    recojet_ECALnCells.push_back(ECALnCells);
    if(ECALnCells>0)
      recojet_ECALtime.push_back(weightedECALTimeCell);
    else
      recojet_ECALtime.push_back(-50);



    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CELLS--------------------------------------"<<std::endl; 
    float weightedMTDTimeCell = 0;
    float totalMTDEnergyCell = 0;
    unsigned int MTDnCells = 0;
    if(fabs(recojet_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDCells(recojet_iter, mtdRecHitsBTL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,1);
    else if(fabs(recojet_iter.eta())>1.4442 && fabs(recojet_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDCells(recojet_iter, mtdRecHitsETL, weightedMTDTimeCell, totalMTDEnergyCell, MTDnCells,0);
    }

    recojet_MTDenergy.push_back(totalMTDEnergyCell);
    recojet_MTDnCells.push_back(MTDnCells);
    if(MTDnCells>0)
      recojet_MTDtime.push_back(weightedMTDTimeCell);
    else
      recojet_MTDtime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM MTD CLUSTERS--------------------------------------"<<std::endl; 
    float weightedMTDTimeClu = 0;
    float totalMTDEnergyClu = 0;
    unsigned int MTDnClus = 0;
    if(fabs(recojet_iter.eta())<1.4442)
      _jetTimingTools.jetTimeFromMTDClus(recojet_iter, _btlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,1);
    else if(fabs(recojet_iter.eta())>1.4442 && fabs(recojet_iter.eta()) <3.0){
      _jetTimingTools.jetTimeFromMTDClus(recojet_iter, _etlRecCluH, weightedMTDTimeClu, totalMTDEnergyClu, MTDnClus,0);
    }
    recojet_MTDCluenergy.push_back(totalMTDEnergyClu);
    recojet_MTDnClus.push_back(MTDnClus);
    if(MTDnClus>0)
      recojet_MTDClutime.push_back(weightedMTDTimeClu);
    else
      recojet_MTDClutime.push_back(-50);


    if(debug)std::cout<<" [DEBUG MODE] ------------- COMPUTE JET TIME FROM HGCAL --------------------------------------"<<std::endl; 
    float weightedHGCALTimeTrackster = 0;
    float totalHGCALEnergyTrackster = 0;
    unsigned int HGCALnTracksters = 0;
  // Handle<View<ticl::Trackster>> tracksterEMH;
  // Handle<View<ticl::Trackster>> tracksterMergeH;
  // Handle<View<ticl::Trackster>> tracksterHADH;
  // Handle<View<ticl::Trackster>> tracksterTrkEMH;
  // Handle<View<ticl::Trackster>> tracksterTrkH;
    std::vector<ticl::Trackster> tracksters;
    // tracksters.reserve(tracksters.size() + distance(tracksterEMH->begin(),tracksterEMH->end()));
    // tracksters.insert(tracksters.end(),tracksterEMH->begin(),tracksterEMH->end());
    tracksters.reserve(tracksters.size() + distance(tracksterMergeH->begin(),tracksterMergeH->end()));
    tracksters.insert(tracksters.end(),tracksterMergeH->begin(),tracksterMergeH->end());
    // tracksters.reserve(tracksters.size() + distance(tracksterHADH->begin(),tracksterHADH->end()));
    // tracksters.insert(tracksters.end(),tracksterHADH->begin(),tracksterHADH->end());
    // tracksters.reserve(tracksters.size() + distance(tracksterTrkEMH->begin(),tracksterTrkEMH->end()));
    // tracksters.insert(tracksters.end(),tracksterTrkEMH->begin(),tracksterTrkEMH->end());
    // tracksters.reserve(tracksters.size() + distance(tracksterTrkH->begin(),tracksterTrkH->end()));
    // tracksters.insert(tracksters.end(),tracksterTrkH->begin(),tracksterTrkH->end());
    _jetTimingTools.jetTimeFromHgcalTracksters(recojet_iter, tracksters, weightedHGCALTimeTrackster, totalHGCALEnergyTrackster, HGCALnTracksters);
    recojet_HGCALenergy.push_back(totalHGCALEnergyTrackster);
    recojet_HGCALnTracksters.push_back(HGCALnTracksters);
    if(HGCALnTracksters>0)
      recojet_HGCALtime.push_back(weightedHGCALTimeTrackster);
    else
      recojet_HGCALtime.push_back(-50);


  }

  // --- setup tree values                                                                                                                                   
  initTreeStructure();
  clearVectors();

  tree_.nLLP = nLLP;
  for (int iL = 0; iL < nLLP; iL++){
    tree_.LLP_pt.push_back(LLP_pt[iL]);
    tree_.LLP_eta.push_back(LLP_eta[iL]);
    tree_.LLP_phi.push_back(LLP_phi[iL]);
    tree_.LLP_mass.push_back(LLP_mass[iL]);
  }

  tree_.ngen        = ngen;
  for (int ig = 0; ig < ngen; ig++){
    tree_.e_ctau.push_back(e_ctau[ig]);
    tree_.e_pt.push_back(e_pt[ig]);
    tree_.e_eta.push_back(e_eta[ig]);
    tree_.e_phi.push_back(e_phi[ig]); 
    tree_.e_vx.push_back(e_vx[ig]);
    tree_.e_vy.push_back(e_vy[ig]);
    tree_.e_vz.push_back(e_vz[ig]); 
    tree_.e_ebphi.push_back(e_ebphi[ig]); 
    tree_.e_ebeta.push_back(e_ebeta[ig]); 
    tree_.e_ebdelay.push_back(e_ebdelay[ig]); 
    tree_.e_hgphi.push_back(e_hgphi[ig]); 
    tree_.e_hgeta.push_back(e_hgeta[ig]); 
    tree_.e_hgdelay.push_back(e_hgdelay[ig]); 
  }

  tree_.ntrack = ntrack;
  for (int it = 0; it < ntrack; it++){
    tree_.track_pt.push_back(track_pt[it]);
    tree_.track_eta.push_back(track_eta[it]);
    tree_.track_phi.push_back(track_phi[it]);

    tree_.track_vx.push_back(track_vx[it]);
    tree_.track_vy.push_back(track_vy[it]);
    tree_.track_vz.push_back(track_vz[it]);

  }

  tree_.nrecojets        = nrecojets;
  for (int ij = 0; ij < nrecojets; ij++){
    tree_.recojet_closestGenIndex.push_back(recojet_closestGenIndex[ij]);
    tree_.recojet_closestGenR.push_back(recojet_closestGenR[ij]);
    tree_.recojet_closestEbGenIndex.push_back(recojet_closestEbGenIndex[ij]);
    tree_.recojet_closestEbGenR.push_back(recojet_closestEbGenR[ij]);
    tree_.recojet_closestHgGenIndex.push_back(recojet_closestHgGenIndex[ij]);
    tree_.recojet_closestHgGenR.push_back(recojet_closestHgGenR[ij]);
    tree_.recojet_pt.push_back(recojet_pt[ij]);
    tree_.recojet_eta.push_back(recojet_eta[ij]);
    tree_.recojet_phi.push_back(recojet_phi[ij]); 
    tree_.recojet_e.push_back(recojet_e[ij]);
    tree_.recojet_ECALtime.push_back(recojet_ECALtime[ij]);
    tree_.recojet_ECALenergy.push_back(recojet_ECALenergy[ij]);
    tree_.recojet_ECALnCells.push_back(recojet_ECALnCells[ij]);
    tree_.recojet_MTDtime.push_back(recojet_MTDtime[ij]);
    tree_.recojet_MTDenergy.push_back(recojet_MTDenergy[ij]);
    tree_.recojet_MTDnCells.push_back(recojet_MTDnCells[ij]);
    tree_.recojet_MTDClutime.push_back(recojet_MTDClutime[ij]);
    tree_.recojet_MTDCluenergy.push_back(recojet_MTDCluenergy[ij]);
    tree_.recojet_MTDnClus.push_back(recojet_MTDnClus[ij]);
    tree_.recojet_HGCALtime.push_back(recojet_HGCALtime[ij]);
    tree_.recojet_HGCALenergy.push_back(recojet_HGCALenergy[ij]);
    tree_.recojet_HGCALnTracksters.push_back(recojet_HGCALnTracksters[ij]);
  }


  // --- fill tree
  tree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void Phase2TimingAnalyzer::beginJob() {

  bool verbose_ = true;
  if (verbose_) std::cout << "Starting job" << std::endl;
  // --- set up output tree                                                                                                                                  
  tree = fs->make<TTree>("tree","tree");
  tree->Branch("ngen",              &tree_.ngen,                "ngen/I");
  
  tree->Branch("LLP_eta", &tree_.LLP_eta);
  tree->Branch("LLP_phi", &tree_.LLP_phi);
  tree->Branch("LLP_pt", &tree_.LLP_pt);
  tree->Branch("LLP_mass", &tree_.LLP_mass);

  tree->Branch("e_eta", &tree_.e_eta);
  tree->Branch("e_phi", &tree_.e_phi);
  tree->Branch("e_pt", &tree_.e_pt);
  tree->Branch("e_vx", &tree_.e_vx);
  tree->Branch("e_ebeta", &tree_.e_ebeta);
  tree->Branch("e_ebphi", &tree_.e_ebphi);
  tree->Branch("e_ebdelay", &tree_.e_ebdelay);
  tree->Branch("e_hgeta", &tree_.e_hgeta);
  tree->Branch("e_hgphi", &tree_.e_hgphi);
  tree->Branch("e_hgdelay", &tree_.e_hgdelay);
  tree->Branch("e_vx", &tree_.e_vx);
  tree->Branch("e_vy", &tree_.e_vy);
  tree->Branch("e_vz", &tree_.e_vz);
  
  tree->Branch("e_ctau", &tree_.e_ctau);

  tree->Branch("track_eta", &tree_.track_eta);
  tree->Branch("track_phi", &tree_.track_phi);
  tree->Branch("track_pt", &tree_.track_pt);
  tree->Branch("track_vx", &tree_.track_vx);
  tree->Branch("track_vy", &tree_.track_vy);
  tree->Branch("track_vz", &tree_.track_vz);


  tree->Branch("nrecojets",              &tree_.nrecojets,                "nrecojets/I");
  tree->Branch("recoJet_pt",             &tree_.recojet_pt);
  tree->Branch("recoJet_eta",             &tree_.recojet_eta);
  tree->Branch("recoJet_phi",             &tree_.recojet_phi);
  tree->Branch("recoJet_e",             &tree_.recojet_e);
  tree->Branch("recoJet_ECALtime",             &tree_.recojet_ECALtime);
  tree->Branch("recoJet_ECALenergy",             &tree_.recojet_ECALenergy);
  tree->Branch("recoJet_ECALnCells",             &tree_.recojet_ECALnCells);
  tree->Branch("recoJet_MTDtime",             &tree_.recojet_MTDtime);
  tree->Branch("recoJet_MTDenergy",             &tree_.recojet_MTDenergy);
  tree->Branch("recoJet_MTDnCells",             &tree_.recojet_MTDnCells);
  tree->Branch("recoJet_MTDClutime",             &tree_.recojet_MTDClutime);
  tree->Branch("recoJet_MTDCluenergy",             &tree_.recojet_MTDCluenergy);
  tree->Branch("recoJet_MTDnClus",             &tree_.recojet_MTDnClus);
  tree->Branch("recoJet_HGCALtime",             &tree_.recojet_HGCALtime);
  tree->Branch("recoJet_HGCALenergy",             &tree_.recojet_HGCALenergy);
  tree->Branch("recoJet_HGCALnTracksters",             &tree_.recojet_HGCALnTracksters);
  tree->Branch("recoJet_closestGenIndex",&tree_.recojet_closestGenIndex);
  tree->Branch("recoJet_closestGenR",&tree_.recojet_closestGenR);
  tree->Branch("recoJet_closestEbGenIndex",&tree_.recojet_closestEbGenIndex);
  tree->Branch("recoJet_closestEbGenR",&tree_.recojet_closestEbGenR);
  tree->Branch("recoJet_closestHgGenIndex",&tree_.recojet_closestHgGenIndex);
  tree->Branch("recoJet_closestHgGenR",&tree_.recojet_closestHgGenR);
}


// ------------ initialize trees ------------                                                                                                                 
void Phase2TimingAnalyzer::initTreeStructure()
{
}


void Phase2TimingAnalyzer::clearVectors()
{

  tree_.recojet_pt.clear();
  tree_.recojet_eta.clear();
  tree_.recojet_phi.clear();
  tree_.recojet_e.clear();
  tree_.recojet_ECALtime.clear();
  tree_.recojet_ECALenergy.clear();
  tree_.recojet_ECALnCells.clear();
  tree_.recojet_MTDtime.clear();
  tree_.recojet_MTDenergy.clear();
  tree_.recojet_MTDnCells.clear();
  tree_.recojet_MTDClutime.clear();
  tree_.recojet_MTDCluenergy.clear();
  tree_.recojet_MTDnClus.clear();
  tree_.recojet_HGCALtime.clear();
  tree_.recojet_HGCALenergy.clear();
  tree_.recojet_HGCALnTracksters.clear();
  tree_.recojet_closestGenIndex.clear();
  tree_.recojet_closestGenR.clear();
  tree_.recojet_closestEbGenIndex.clear();
  tree_.recojet_closestEbGenR.clear();
  tree_.recojet_closestHgGenIndex.clear();
  tree_.recojet_closestHgGenR.clear();

  tree_.track_pt.clear();
  tree_.track_eta.clear();
  tree_.track_phi.clear();
  tree_.track_vx.clear();
  tree_.track_vy.clear();
  tree_.track_vz.clear();

  tree_.e_pt.clear();
  tree_.e_eta.clear();
  tree_.e_phi.clear();
  tree_.e_ctau.clear();

  tree_.LLP_pt.clear();
  tree_.LLP_eta.clear();
  tree_.LLP_phi.clear();
  tree_.LLP_mass.clear();

  tree_.e_vx.clear();
  tree_.e_vy.clear();
  tree_.e_vz.clear();
  tree_.e_ebphi.clear();
  tree_.e_ebeta.clear();
  tree_.e_ebdelay.clear();
  tree_.e_hgphi.clear();
  tree_.e_hgeta.clear();
  tree_.e_hgdelay.clear();
}

// ------------ method called once each job just after ending the event loop  ------------
void Phase2TimingAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2TimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2TimingAnalyzer);
