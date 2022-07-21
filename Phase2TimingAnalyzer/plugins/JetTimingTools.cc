#include "Phase2Timing/Phase2TimingAnalyzer/plugins/JetTimingTools.h"


JetTimingTools::JetTimingTools(edm::ConsumesCollector && cc):
    caloGeometryToken_(cc.esConsumes()),
    mtdGeometryToken_(cc.esConsumes()),
    mtdTopologyToken_(cc.esConsumes()),
    ecalCellEnergyThresh_(0.5),
    ecalCellTimeThresh_(12.5),
    ecalCellTimeErrorThresh_(100),
    hcalCellEnergyThresh_(1.),
    hcalCellTimeThresh_(12.5),
    hgcalTracksterEnergyThresh_(0.5),
    hgcalTracksterTimeThresh_(12.5),
    hgcalTracksterTimeErrorThresh_(100),
    mtdCellEnergyThresh_(0.5),
    mtdCellTimeThresh_(25.),
    mtdCellTimeErrorThresh_(100),
    matchingRadius2_(0.16)
{

}

void  JetTimingTools::setMatchingRadius(double matchingRadius) {matchingRadius2_ = matchingRadius*matchingRadius; }

void  JetTimingTools::setHcalCellEnergyThreshold(double ecalCellEnergyThresh){ ecalCellEnergyThresh_ = ecalCellEnergyThresh;}

void  JetTimingTools::setHcalCellTimeThreshold(double ecalCellTimeThresh){ ecalCellTimeThresh_ = ecalCellTimeThresh; }

void  JetTimingTools::setEcalCellEnergyThreshold(double ecalCellEnergyThresh){ ecalCellEnergyThresh_ = ecalCellEnergyThresh;}

void  JetTimingTools::setEcalCellTimeThreshold(double ecalCellTimeThresh){ ecalCellTimeThresh_ = ecalCellTimeThresh; }

void  JetTimingTools::setEcalCellTimeErrorThreshold(double ecalCellTimeErrorThresh){ ecalCellTimeErrorThresh_ = ecalCellTimeErrorThresh;}

void  JetTimingTools::setHgcalTracksterEnergyThreshold(double hgcalTracksterEnergyThresh){ hgcalTracksterEnergyThresh_ = hgcalTracksterEnergyThresh;}

void  JetTimingTools::setHgcalTracksterTimeThreshold(double hgcalTracksterTimeThresh){ hgcalTracksterTimeThresh_ = hgcalTracksterTimeThresh; }

void  JetTimingTools::setHgcalTracksterTimeErrorThreshold(double hgcalTracksterTimeErrorThresh){ hgcalTracksterTimeErrorThresh_ = hgcalTracksterTimeErrorThresh;}

void  JetTimingTools::setMTDCellEnergyThreshold(double mtdCellEnergyThresh){ mtdCellEnergyThresh_ = mtdCellEnergyThresh;}

void  JetTimingTools::setMTDCellTimeThreshold(double mtdCellTimeThresh){ mtdCellTimeThresh_ = mtdCellTimeThresh; }

void  JetTimingTools::setMTDCellTimeErrorThreshold(double mtdCellTimeErrorThresh){ mtdCellTimeErrorThresh_ = mtdCellTimeErrorThresh;}

void JetTimingTools::init(const edm::EventSetup &es){
  caloGeometry_ = es.getHandle(caloGeometryToken_);
  mtdGeometry_ = es.getHandle(mtdGeometryToken_);
  mtdTopology_ = es.getHandle(mtdTopologyToken_);
}

//Calculate gen delay
std::vector<double> JetTimingTools::surfaceIntersection(const reco::GenParticle& genParticle, reco::GenParticle& genParticleMother, double radius){
    //all in cm
    // double radius = 130;
    double lightSpeed = 29979245800;

    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> particlePos = genParticle.vertex();
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> parentPos = genParticleMother.vertex();
    double genParticleBeta = genParticleMother.p()/genParticleMother.energy();
    double parentTime = TMath::Sqrt((particlePos-parentPos).Mag2())/(lightSpeed*genParticleBeta);
    double a = (genParticle.px()*genParticle.px() + genParticle.py()*genParticle.py())
    *lightSpeed*lightSpeed/(genParticle.energy()*genParticle.energy());
    double b = 2*(genParticle.vx()*genParticle.px() + genParticle.vy()*genParticle.py())*lightSpeed/genParticle.energy();
    double c = (genParticle.vx()*genParticle.vx() + genParticle.vy()*genParticle.vy()) - radius*radius;
    if (c > 0) return {1000.,1000.,1000.,1000.};
    double sqrt_disc = TMath::Sqrt(b*b - 4*a*c);
    double t1 = (-b + sqrt_disc)/(2*a);
    double xPos = t1*(genParticle.px()/genParticle.energy())*lightSpeed + genParticle.vx();
    double yPos = t1*(genParticle.py()/genParticle.energy())*lightSpeed + genParticle.vy();
    double zPos = t1*(genParticle.pz()/genParticle.energy())*lightSpeed + genParticle.vz();
    double phi = TMath::ATan2(yPos,xPos);
    double theta = TMath::ACos(zPos/TMath::Sqrt(xPos*xPos+yPos*yPos+zPos*zPos));
    double eta = -TMath::Log(TMath::Tan(theta/2.));
    // double t1Light = TMath::Sqrt(zPos*zPos+radius*radius)/lightSpeed;
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> intersectionXYZ(xPos,yPos,zPos);
    double t1Light = TMath::Sqrt((intersectionXYZ-parentPos).Mag2())/lightSpeed;
    std::vector<double> outVec = {eta,phi,t1Light,(t1+parentTime-t1Light),t1};
    return outVec;
}
std::vector<double> JetTimingTools::endCapIntersection(const reco::GenParticle& genParticle, reco::GenParticle& genParticleMother, double zMin, double zMax){
    //all in cm
    // double radius = 130;
    double lightSpeed = 29979245800;

    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> particlePos = genParticle.vertex();
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> parentPos = genParticleMother.vertex();
    double genParticleBeta = genParticleMother.p()/genParticleMother.energy();
    double parentTime = TMath::Sqrt((particlePos-parentPos).Mag2())/(lightSpeed*genParticleBeta);
    double xPos;
    double yPos;
    double zPos;
    double t1;
    if (fabs(parentPos.z()) > zMax) return {1000.,1000.,1000.,1000.};
    else if (fabs(particlePos.z()) > zMin and fabs(particlePos.z()) < zMax)
    {
        xPos = particlePos.x();
        yPos = particlePos.y();
        zPos = particlePos.z();
	t1 = 0;
    }
    else {
	if (genParticle.pz() > 0) t1 = (zMin-genParticle.vz())/((genParticle.pz()/genParticle.energy())*lightSpeed);
	else t1 = (-zMin-genParticle.vz())/((genParticle.pz()/genParticle.energy())*lightSpeed);
	xPos = t1*(genParticle.px()/genParticle.energy())*lightSpeed + genParticle.vx();
	yPos = t1*(genParticle.py()/genParticle.energy())*lightSpeed + genParticle.vy();
	zPos = t1*(genParticle.pz()/genParticle.energy())*lightSpeed + genParticle.vz();
    }
    double phi = TMath::ATan2(yPos,xPos);
    double theta = TMath::ACos(zPos/TMath::Sqrt(xPos*xPos+yPos*yPos+zPos*zPos));
    double eta = -TMath::Log(TMath::Tan(theta/2.));
    ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>> intersectionXYZ(xPos,yPos,zPos);
    double t1Light = TMath::Sqrt((intersectionXYZ-parentPos).Mag2())/lightSpeed;
    std::vector<double> outVec = {eta,phi,t1Light,(t1+parentTime-t1Light),t1};
    return outVec;
}
//calculate jet time from hcal cells
void JetTimingTools::jetTimeFromHcalCells(
    const reco::Jet& jet,
    const edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>& hcalRecHits,
    float& weightedTimeCell,
    float& totalEmEnergyCell,
    uint& nCells,
    bool useTDCHCAL){
    for (auto const& hcalRH : hcalRecHits) {
    if (hcalRH.energy() < hcalCellEnergyThresh_)
      continue;
    // if (hcalRH.timeError() <= 0. || hcalRH.timeError() > hcalCellTimeErrorThresh_)
    //   continue;
    if (useTDCHCAL and fabs(hcalRH.timeFalling()) > hcalCellTimeThresh_)
      continue;
    else if (!useTDCHCAL and fabs(hcalRH.time()) > hcalCellTimeThresh_)
	continue;
    auto const pos = caloGeometry_->getPosition(hcalRH.detid());
    if (reco::deltaR2(jet, pos) > matchingRadius2_)
      continue;
    if (useTDCHCAL) weightedTimeCell += hcalRH.timeFalling() * hcalRH.energy() * sin(pos.theta());
    else  weightedTimeCell += hcalRH.time() * hcalRH.energy() * sin(pos.theta());
    totalEmEnergyCell += hcalRH.energy() * sin(pos.theta());
    nCells++;
  }
  if (totalEmEnergyCell > 0) {
    weightedTimeCell /= totalEmEnergyCell;
    }
}
//calculate jet time from ecal cells
void JetTimingTools::jetTimeFromEcalCells(
    const reco::Jet& jet,
    const edm::SortedCollection<EcalRecHit, edm::StrictWeakOrdering<EcalRecHit>>& ecalRecHits,
    float& weightedTimeCell,
    float& totalEmEnergyCell,
    uint& nCells) {
    for (auto const& ecalRH : ecalRecHits) {
    if (ecalRH.checkFlag(EcalRecHit::kSaturated) || ecalRH.checkFlag(EcalRecHit::kLeadingEdgeRecovered) ||
        ecalRH.checkFlag(EcalRecHit::kPoorReco) || ecalRH.checkFlag(EcalRecHit::kWeird) ||
        ecalRH.checkFlag(EcalRecHit::kDiWeird))
      continue;
    if (ecalRH.energy() < ecalCellEnergyThresh_)
      continue;
    if (ecalRH.timeError() <= 0. || ecalRH.timeError() > ecalCellTimeErrorThresh_)
      continue;
    if (fabs(ecalRH.time()) > ecalCellTimeThresh_)
      continue;

    auto const pos = caloGeometry_->getPosition(ecalRH.detid());
    if (reco::deltaR2(jet, pos) > matchingRadius2_)
      continue;
    weightedTimeCell += ecalRH.time() * ecalRH.energy() * sin(pos.theta());
    totalEmEnergyCell += ecalRH.energy() * sin(pos.theta());
    nCells++;
  }
  if (totalEmEnergyCell > 0) {
    weightedTimeCell /= totalEmEnergyCell;
    }
}
//calculate jet time from hgcal tracksters
void JetTimingTools::jetTimeFromHgcalTracksters(
    const reco::Jet& jet,
    const std::vector<ticl::Trackster>& tracksters,
    float& weightedTimeTrackster,
    float& totalEnergyTrackster,
    uint& nTracksters) {
    for (auto const& trackster : tracksters) {
    if (trackster.regressed_energy() < hgcalTracksterEnergyThresh_)
      continue;
    if (trackster.timeError() <= 0. || trackster.timeError() > hgcalTracksterTimeErrorThresh_)
      continue;
    if (fabs(trackster.time()) > hgcalTracksterTimeThresh_)
      continue;
    auto const pos = trackster.barycenter();
    if (reco::deltaR2(jet, pos) > matchingRadius2_)
      continue;
    weightedTimeTrackster += trackster.time() * trackster.regressed_energy();
    totalEnergyTrackster += trackster.regressed_energy();
    nTracksters++;
  }
  if (totalEnergyTrackster > 0) {
    weightedTimeTrackster /= totalEnergyTrackster;
    }
}

//calculate jet time from mtd rechits
void JetTimingTools::jetTimeFromMTDCells(
    const reco::Photon& photon,
    const edm::SortedCollection<FTLRecHit,edm::StrictWeakOrdering<FTLRecHit> >& mtdRecHits,
    float& weightedTimeCell,
    float& totalEmEnergyCell,
    uint& nCells, bool isBTL) {//reco::Photon& photon
  
  for (auto const& mtdRH : mtdRecHits) {
    /* if (mtdRH.energy() < mtdCellEnergyThresh_)
            continue;
    if (mtdRH.timeError() <= 0. || mtdRH.timeError() > mtdCellTimeErrorThresh_)
      continue;
    if (fabs(mtdRH.time()) > mtdCellTimeThresh_)
      continue;
    */
    double dR;
    double gp_x;
    double gp_y;
    double gp_z;
    double gp_theta;

    if(isBTL){
      BTLDetId detId = mtdRH.id();
      DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(mtdTopology_->getMTDTopologyMode()));
      const MTDGeomDet* thedet = mtdGeometry_->idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());      
      Local3DPoint local_point(0., 0., 0.);
      local_point = topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
      const auto& global_point = thedet->toGlobal(local_point);
      dR=reco::deltaR2(photon, global_point);
      gp_x = global_point.x();
      gp_y = global_point.y();
      gp_z = global_point.z();
      gp_theta = global_point.theta();
    }else{
      ETLDetId detId = mtdRH.id();
      DetId geoId = detId.geographicalId();
      const MTDGeomDet* thedet = mtdGeometry_->idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      Local3DPoint local_point(topo.localX(mtdRH.row()), topo.localY(mtdRH.column()), 0.);
      const auto& global_point = thedet->toGlobal(local_point);
      dR=reco::deltaR2(photon, global_point);
      gp_x = global_point.x();
      gp_y = global_point.y();
      gp_z = global_point.z();
      gp_theta = global_point.theta();
    }

    if (dR > matchingRadius2_)
      continue;
    double tof = (TMath::Sqrt((gp_x)*(gp_x)+(gp_y)*(gp_y)+(gp_z)*(gp_z))/10)*1./(2.9979); 
    //      std::cout<< "MTD RECHIT: eta: "<<gp_eta<<" phi: "<<gp_phi<<" x: "<<gp_x<<" y: "<<gp_y<<" z: "<<gp_z<<" e: "<<mtdRH.energy()<<" time: "<<mtdRH.time()<<" tof: "<<tof<<" timeDiff: "<<(mtdRH.time()-tof)<<std::endl;
    weightedTimeCell += (mtdRH.time()-tof) * mtdRH.energy() * sin(gp_theta);
    totalEmEnergyCell += mtdRH.energy() * sin(gp_theta);
    nCells++;
  }
  if (totalEmEnergyCell > 0) {
    weightedTimeCell /= totalEmEnergyCell;
  }
}



//calculate jet time from mtd clusters
void JetTimingTools::jetTimeFromMTDClus(
    const reco::Jet& jet,
    const edm::Handle<FTLClusterCollection>& mtdRecClusters,
    float& weightedTimeClu,
    float& totalEnergyClu,
    uint& nClu, bool isBTL) {
  
  double dR;
  double gp_x;
  double gp_y;
  double gp_z;
  double gp_theta;
  for (auto const& mtdClu : *mtdRecClusters) {
    for(const auto& cluster : mtdClu){
      if(isBTL){ 
	if (cluster.energy() < 1.)
	  continue;    
	
	BTLDetId cluId = cluster.id();
	DetId detIdObject(cluId);
	const auto& genericDet = mtdGeometry_->idToDetUnit(detIdObject);
	if(genericDet == nullptr)
	  continue;
	const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(genericDet->topology());
	const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
	
	// --- Cluster position in the module reference frame
	Local3DPoint local_point(topo.localX(cluster.x()), topo.localY(cluster.y()), 0.);
	const auto& global_point = genericDet->toGlobal(local_point);
	dR=reco::deltaR2(jet, global_point);
	gp_x = global_point.x();
	gp_y = global_point.y();
	gp_z = global_point.z();
      gp_theta = global_point.theta();
      } else { 
	if (cluster.energy() < 0.001)
	  continue;    
	
	ETLDetId cluId = cluster.id();
	DetId detIdObject(cluId);
	const auto& genericDet = mtdGeometry_->idToDetUnit(detIdObject);
	if(genericDet == nullptr)
	  continue;

	const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(genericDet->topology());
	const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
	
	// --- Cluster position in the module reference frame
	Local3DPoint local_point(topo.localX(cluster.x()), topo.localY(cluster.y()), 0.);
	const auto& global_point = genericDet->toGlobal(local_point);
	dR=reco::deltaR2(jet, global_point);
	gp_x = global_point.x();
	gp_y = global_point.y();
	gp_z = global_point.z();
	gp_theta = global_point.theta();
      }


    if (dR > matchingRadius2_)
      continue;
    double tof = (TMath::Sqrt((gp_x)*(gp_x)+(gp_y)*(gp_y)+(gp_z)*(gp_z))/10)*1./(2.9979); 
    //    std::cout<< "MTD CLUSTER:: "<<" x: "<<gp_x<<" y: "<<gp_y<<" z: "<<gp_z<<" e: "<<cluster.energy()<<" time: "<<cluster.time()<<" tof: "<<tof<<" timeDiff: "<<(cluster.time()-tof)<<std::endl;
    weightedTimeClu += (cluster.time()) * cluster.energy() * sin(gp_theta);
    totalEnergyClu += cluster.energy() * sin(gp_theta);
    nClu++;

    }
  }
 
  if (totalEnergyClu > 0) {
    weightedTimeClu /= totalEnergyClu;
  }

}



