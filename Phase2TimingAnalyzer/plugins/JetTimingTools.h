#ifndef JetTimingTools_hh
#define JetTimingTools_hh
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/JetReco/interface/Jet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TMath.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDNumberingBuilder/interface/MTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"

#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
}  // namespace edm

class JetTimingTools {
public:
  explicit JetTimingTools(edm::ConsumesCollector &&);
  // virtual ~JetTimingTools(){};
  void init(const edm::EventSetup &es);
  std::vector<double> surfaceIntersection(const reco::GenParticle&, reco::GenParticle&, double);
  std::vector<double> endCapIntersection(const reco::GenParticle&, reco::GenParticle&, double,double);
  void jetTimeFromEcalCells(const reco::Jet&,
                            const edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>&,
                            float&,
                            float&,
                            uint&);
  void jetTimeFromHcalCells(const reco::Jet&,
                            const edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>&,
                            float&,
                            float&,
                            uint&,
			    bool);
  void jetTimeFromHgcalTracksters(const reco::Jet&,
                            const std::vector<ticl::Trackster>&,
                            float&,
                            float&,
                            uint&);
  void jetTimeFromMTDCells(const reco::Photon&,
                            const edm::SortedCollection<FTLRecHit, edm::StrictWeakOrdering<FTLRecHit>>&,
                            float&,
                            float&,
			   uint&,
			   bool);//const reco::Jet&
  void jetTimeFromMTDClus(const reco::Jet&,
			      const edm::Handle<FTLClusterCollection>&,
			      float&,
			      float&,
			      uint&,
			      bool);
  void setMatchingRadius(double );
  void setEcalCellEnergyThreshold(double);
  void setEcalCellTimeThreshold(double );
  void setEcalCellTimeErrorThreshold(double );
  void setHcalCellEnergyThreshold(double);
  void setHcalCellTimeThreshold(double );
  void setHgcalTracksterEnergyThreshold(double);
  void setHgcalTracksterTimeThreshold(double );
  void setHgcalTracksterTimeErrorThreshold(double );
  void setMTDCellEnergyThreshold(double);
  void setMTDCellTimeThreshold(double );
  void setMTDCellTimeErrorThreshold(double );

private:
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
    edm::ESHandle<CaloGeometry> caloGeometry_;
    
    edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdGeometryToken_;
    edm::ESHandle<MTDGeometry> mtdGeometry_;

    edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdTopologyToken_;
    edm::ESHandle<MTDTopology> mtdTopology_;

    
    // Configurables for timing definition
    double ecalCellEnergyThresh_;
    double ecalCellTimeThresh_;
    double ecalCellTimeErrorThresh_;
    double hcalCellEnergyThresh_;
    double hcalCellTimeThresh_;
    double hcalCellTimeErrorThresh_;
    double hgcalTracksterEnergyThresh_;
    double hgcalTracksterTimeThresh_;
    double hgcalTracksterTimeErrorThresh_;
    double mtdCellEnergyThresh_;
    double mtdCellTimeThresh_;
    double mtdCellTimeErrorThresh_;
    double matchingRadius2_;
};

#endif //JetTimingTools_hh
