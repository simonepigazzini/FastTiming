#ifndef _FTL_DUMP_HITS_
#define _FTL_DUMP_HITS_

#include "TMath.h"

#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "RecoMTD/TransientTrackingRecHit/interface/MTDTransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderWithPropagator.h"
#include "RecoTracker/TransientTrackingRecHit/interface/Traj2TrackHits.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDRingForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetRing.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "PrecisionTiming/FTLAnalysis/interface/FTLHitsTree.h"

#include <memory>
using namespace std;



class FTLDumpHits : public edm::EDAnalyzer
{
public:
  explicit FTLDumpHits(const edm::ParameterSet& pSet);
  ~FTLDumpHits() {};
  
  //---utils
  
  //---methods
  virtual void beginJob() override {};
  virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
  virtual void endJob() override {};
  
  std::string PrintPosition(const GlobalPoint& gp);
  std::string PrintPosition(const LocalPoint& lp);
  
private:
  const MTDGeometry* mtdGeometry_;
  
  //---inputs
  edm::Handle<reco::GenParticleCollection> genParticlesHandle_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

  edm::Handle<std::vector<PSimHit> > simHitsBTLHandle_;
  edm::EDGetTokenT<std::vector<PSimHit> > simHitsBTLToken_;
  edm::Handle<FTLRecHitCollection> recHitsBTLHandle_;
  edm::EDGetTokenT<FTLRecHitCollection> recHitsBTLToken_;    
  edm::Handle<FTLClusterCollection> clustersBTLHandle_;
  edm::EDGetTokenT<FTLClusterCollection> clustersBTLToken_;    

  edm::Handle<std::vector<PSimHit> > simHitsETLHandle_;
  edm::EDGetTokenT<std::vector<PSimHit> > simHitsETLToken_;
  edm::Handle<FTLRecHitCollection> recHitsETLHandle_;
  edm::EDGetTokenT<FTLRecHitCollection> recHitsETLToken_;    
  edm::Handle<FTLClusterCollection> clustersETLHandle_;
  edm::EDGetTokenT<FTLClusterCollection> clustersETLToken_;    

  edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_;
  edm::Handle<edm::View<reco::Track> > tracksHandle_;
  edm::EDGetTokenT<vector<SimVertex> >                 genVtxToken_;
  edm::Handle<vector<SimVertex> >                      genVtxHandle_;    
  edm::ESHandle<TransientTrackBuilder> builder;

  //---options
  BTLDetId::CrysLayout crysLayout_;
  double track_hit_DRMax_;
  double track_hit_distMax_;
  bool verbosity_;
  
  //---outputs
  FTLHitsTree outTree_;
  edm::Service<TFileService> fs_;  
  
};



FTLDumpHits::FTLDumpHits(const edm::ParameterSet& pSet):
  genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
  simHitsBTLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsBTLTag"))),
  recHitsBTLToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("recHitsBTLTag"))),
  clustersBTLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersBTLTag"))),
  simHitsETLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsETLTag"))),
  recHitsETLToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("recHitsETLTag"))),
  clustersETLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersETLTag"))),
  tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag"))),
  genVtxToken_(consumes<vector<SimVertex> >(pSet.getUntrackedParameter<edm::InputTag>("genVtxTag"))),
  crysLayout_((BTLDetId::CrysLayout)(pSet.getUntrackedParameter<int>("crysLayout"))),
  track_hit_DRMax_(pSet.getParameter<double>("track_hit_DRMax")),
  track_hit_distMax_(pSet.getParameter<double>("track_hit_distMax")),
  verbosity_(pSet.getParameter<bool>("verbosity"))
{
  outTree_ = FTLHitsTree(pSet.getUntrackedParameter<string>("treeName").c_str(), "FTLHits tree for FTL studies");
}



void FTLDumpHits::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
  outTree_.Reset();
  
  //---get the MTD geometry
  edm::ESHandle<MTDGeometry> geoHandle;
  setup.get<MTDDigiGeometryRecord>().get(geoHandle);
  mtdGeometry_ = geoHandle.product();
  
  edm::ESHandle<MTDDetLayerGeometry> layerGeo;
  setup.get<MTDRecoGeometryRecord>().get(layerGeo);
  
  //--- get the B field
  edm::ESHandle<MagneticField> theField;
  setup.get<IdealMagneticFieldRecord>().get(theField);
  PropagationDirection dir(alongMomentum);
  SteppingHelixPropagator* propagator = new SteppingHelixPropagator(theField.product(),dir);
  propagator -> setMaterialMode(false);
  propagator -> setNoErrorPropagation(false);

  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  //--- load gen particles
  event.getByToken(genParticlesToken_, genParticlesHandle_);
  auto genParticles = *genParticlesHandle_.product();
  
  //---load sim hits
  event.getByToken(simHitsBTLToken_, simHitsBTLHandle_);
  auto simHitsBTL = *simHitsBTLHandle_.product();
  
  //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
  event.getByToken(recHitsBTLToken_, recHitsBTLHandle_);
  auto recHitsBTL = FTLRecHitCollection();
  if(recHitsBTLHandle_.isValid())
    recHitsBTL = *recHitsBTLHandle_.product();

  event.getByToken(clustersBTLToken_, clustersBTLHandle_);
  // auto clusters = FTLClusterCollection();
  // if(clustersBTLHandle_.isValid())
  auto clustersBTL = *clustersBTLHandle_.product();


  //---load sim hits
  event.getByToken(simHitsETLToken_, simHitsETLHandle_);
  auto simHitsETL = *simHitsETLHandle_.product();
  
  //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
  event.getByToken(recHitsETLToken_, recHitsETLHandle_);
  auto recHitsETL = FTLRecHitCollection();
  if(recHitsETLHandle_.isValid())
    recHitsETL = *recHitsETLHandle_.product();

  event.getByToken(clustersETLToken_, clustersETLHandle_);
  // auto clusters = FTLClusterCollection();
  // if(clustersETLHandle_.isValid())
  auto clustersETL = *clustersETLHandle_.product();
  
  //---load tracks
  event.getByToken(tracksToken_,tracksHandle_);
  auto tracks = *tracksHandle_.product();
  

  event.getByToken(genVtxToken_, genVtxHandle_);    
  const SimVertex* genPV = NULL;
  if(genVtxHandle_.isValid())
    genPV = &(genVtxHandle_.product()->at(0));
  
  
  //---fill the tree - simHits
  outTree_.BTLsimHits_n = 0;
  for(auto simHit : simHitsBTL)
  {
    BTLDetId id = simHit.detUnitId();
    DetId geoId = id.geographicalId( crysLayout_);
    const auto& det = mtdGeometry_ -> idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
    
    double energy = simHit.energyLoss()*1000.;
    double time   = simHit.tof();
    
    int RR = id.mtdRR();
    int module = id.module();
    int modType = id.modType();
    int crystal = id.crystal();
    int ieta = id.ieta(crysLayout_);
    int iphi = id.iphi(crysLayout_);
    
    LocalPoint lp_entry(simHit.entryPoint().x()/10., simHit.entryPoint().y()/10., simHit.entryPoint().z()/10.);
    LocalPoint lp_exit ( simHit.exitPoint().x()/10.,  simHit.exitPoint().y()/10.,  simHit.exitPoint().z()/10.);
    GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,id.row(topo.nrows()),id.column(topo.nrows())));
    GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit,id.row(topo.nrows()),id.column(topo.nrows())));
    
    outTree_.BTLsimHits_n += 1;
    
    outTree_.BTLsimHits_energy->push_back(energy);
    outTree_.BTLsimHits_time->push_back(time);
    outTree_.BTLsimHits_rr->push_back(RR);
    outTree_.BTLsimHits_module->push_back(module);
    outTree_.BTLsimHits_modType->push_back(modType);
    outTree_.BTLsimHits_crystal->push_back(crystal);
    outTree_.BTLsimHits_ieta->push_back(ieta);
    outTree_.BTLsimHits_iphi->push_back(iphi);
    outTree_.BTLsimHits_entry_local_x->push_back(lp_entry.x());
    outTree_.BTLsimHits_entry_local_y->push_back(lp_entry.y());
    outTree_.BTLsimHits_entry_local_z->push_back(lp_entry.z());
    outTree_.BTLsimHits_entry_global_R->push_back(sqrt(gp_entry.perp2()));
    outTree_.BTLsimHits_exit_local_x->push_back(lp_exit.x());
    outTree_.BTLsimHits_exit_local_y->push_back(lp_exit.y());
    outTree_.BTLsimHits_exit_local_z->push_back(lp_exit.z());
    outTree_.BTLsimHits_exit_global_R->push_back(sqrt(gp_exit.perp2()));
  }
  
  
  
  //---fill the tree - recHits
  outTree_.BTLrecHits_n = 0;
  for(auto recHit : recHitsBTL)
  {
    BTLDetId id = recHit.id();
    DetId geoId = id.geographicalId( crysLayout_);
    const auto& det = mtdGeometry_ -> idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
    
    double energy = recHit.energy();
    double time   = recHit.time();
    
    MeasurementPoint mp(recHit.row(),recHit.column());
    LocalPoint lp = topo.localPosition(mp);
    GlobalPoint gp = det->toGlobal(lp);
    
    int RR = id.mtdRR();
    int module = id.module();
    int modType = id.modType();
    int crystal = id.crystal();
    int ieta = id.ieta(crysLayout_);
    int iphi = id.iphi(crysLayout_);
    
    outTree_.BTLrecHits_n += 1;
    
    outTree_.BTLrecHits_energy->push_back(energy);
    outTree_.BTLrecHits_time->push_back(time);
    outTree_.BTLrecHits_rr->push_back(RR);
    outTree_.BTLrecHits_module->push_back(module);
    outTree_.BTLrecHits_modType->push_back(modType);
    outTree_.BTLrecHits_crystal->push_back(crystal);
    outTree_.BTLrecHits_ieta->push_back(ieta);
    outTree_.BTLrecHits_iphi->push_back(iphi);
    outTree_.BTLrecHits_local_x->push_back(lp.x());
    outTree_.BTLrecHits_local_y->push_back(lp.y());
    outTree_.BTLrecHits_local_z->push_back(lp.z());
    outTree_.BTLrecHits_global_R->push_back(sqrt(gp.perp2()));
  }


  //---fill the tree - recHits
  outTree_.BTLclusters_n = 0;
  
  for (auto clusIt : clustersBTL)
    {    
      DetId id = clusIt.detId();
      const auto& det = mtdGeometry_ -> idToDet(id);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      for ( auto cluster : clusIt)
	{
	  float energy = cluster.energy();
	  float time   = cluster.time();
	  float x=cluster.x();
	  float y=cluster.y();
	  int size=cluster.size();
	  int sizeX=cluster.sizeX();
	  int sizeY=cluster.sizeY();
	  float seed_energy=cluster.seed().energy();
	  float seed_time=cluster.seed().time();
	  int seed_x=cluster.seed().x();
	  int seed_y=cluster.seed().y();
	  
	  MeasurementPoint mp(cluster.x(),cluster.y());
	  LocalPoint lp = topo.localPosition(mp);
	  GlobalPoint gp = det->toGlobal(lp);
	  
	  MTDDetId mtdId(id);
	  int RR = 0;
	  int module = 0;
	  int modType = 0;
	  int crystal = 0;
	  int ieta = 0;
	  int iphi = 0;
	  
	  if ( mtdId.mtdSubDetector() == MTDDetId::BTL )
	    {
	      BTLDetId btlId(id);
	      RR = btlId.mtdRR();
	      module = btlId.module();
	      modType = btlId.modType();
	      crystal = btlId.crystal();
	      ieta = btlId.ieta(crysLayout_);
	      iphi = btlId.iphi(crysLayout_);
	    }
	  
	  outTree_.BTLclusters_n += 1;    
	  outTree_.BTLclusters_size->push_back(size);
	  outTree_.BTLclusters_size_x->push_back(sizeX);
	  outTree_.BTLclusters_size_y->push_back(sizeY);
	  outTree_.BTLclusters_energy->push_back(energy);
	  outTree_.BTLclusters_time->push_back(time);
	  outTree_.BTLclusters_rr->push_back(RR);
	  outTree_.BTLclusters_module->push_back(module);
	  outTree_.BTLclusters_modType->push_back(modType);
	  outTree_.BTLclusters_crystal->push_back(crystal);
	  outTree_.BTLclusters_ieta->push_back(ieta);
	  outTree_.BTLclusters_iphi->push_back(iphi);
	  outTree_.BTLclusters_x->push_back(x);
	  outTree_.BTLclusters_y->push_back(y);
	  outTree_.BTLclusters_seed_energy->push_back(seed_energy);
	  outTree_.BTLclusters_seed_time->push_back(seed_time);
	  outTree_.BTLclusters_seed_x->push_back(seed_x);
	  outTree_.BTLclusters_seed_y->push_back(seed_y);
	  outTree_.BTLclusters_local_x->push_back(lp.x());
	  outTree_.BTLclusters_local_y->push_back(lp.y());
	  outTree_.BTLclusters_local_z->push_back(lp.z());
	  outTree_.BTLclusters_global_R->push_back(sqrt(gp.perp2()));
	}
    }


  outTree_.ETLsimHits_n = 0;
  for(auto simHit : simHitsETL)
  {
    ETLDetId id = simHit.detUnitId();
    DetId geoId = id.geographicalId();
    const auto& det = mtdGeometry_ -> idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
    
    double energy = simHit.energyLoss()*1000.;
    double time   = simHit.tof();

    int RR = id.mtdRR();
    int module = id.module();
    int modType = id.modType();
    int crystal = 0;
    int ieta = 0;
    int iphi = 0;

    //need to add protection of oot
    // ETL is already in module-local coordinates so just scale to cm from mm
    Local3DPoint simscaled(0.1*simHit.entryPoint().x(),0.1*simHit.entryPoint().y(),0.1*simHit.entryPoint().z());
    const auto& thepixel = topo.pixel(simscaled); // mm -> cm here is the switch                               
    const uint8_t row(thepixel.first), col(thepixel.second);

    LocalPoint lp_entry(simHit.entryPoint().x()/10., simHit.entryPoint().y()/10., simHit.entryPoint().z()/10.);
    LocalPoint lp_exit ( simHit.exitPoint().x()/10.,  simHit.exitPoint().y()/10.,  simHit.exitPoint().z()/10.);
    GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,row,col));
    GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit,row,col));
    
    outTree_.ETLsimHits_n += 1;
    
    outTree_.ETLsimHits_energy->push_back(energy);
    outTree_.ETLsimHits_time->push_back(time);
    outTree_.ETLsimHits_rr->push_back(RR);
    outTree_.ETLsimHits_module->push_back(module);
    outTree_.ETLsimHits_modType->push_back(modType);
    outTree_.ETLsimHits_crystal->push_back(crystal);
    outTree_.ETLsimHits_ieta->push_back(ieta);
    outTree_.ETLsimHits_iphi->push_back(iphi);
    outTree_.ETLsimHits_entry_local_x->push_back(lp_entry.x());
    outTree_.ETLsimHits_entry_local_y->push_back(lp_entry.y());
    outTree_.ETLsimHits_entry_local_z->push_back(lp_entry.z());
    outTree_.ETLsimHits_entry_global_R->push_back(sqrt(gp_entry.perp2()));
    outTree_.ETLsimHits_exit_local_x->push_back(lp_exit.x());
    outTree_.ETLsimHits_exit_local_y->push_back(lp_exit.y());
    outTree_.ETLsimHits_exit_local_z->push_back(lp_exit.z());
    outTree_.ETLsimHits_exit_global_R->push_back(sqrt(gp_exit.perp2()));
  }
  
  
  
  //---fill the tree - recHits
  outTree_.ETLrecHits_n = 0;
  for(auto recHit : recHitsETL)
  {
    ETLDetId id = recHit.id();
    DetId geoId = id.geographicalId( );
    const auto& det = mtdGeometry_ -> idToDet(geoId);
    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
    
    double energy = recHit.energy();
    double time   = recHit.time();
    
    MeasurementPoint mp(recHit.row(),recHit.column());
    LocalPoint lp = topo.localPosition(mp);
    GlobalPoint gp = det->toGlobal(lp);
    
    int RR = id.mtdRR();
    int module = id.module();
    int modType = id.modType();
    int crystal = 0;
    int ieta = 0;
    int iphi = 0;
    
    outTree_.ETLrecHits_n += 1;
    
    outTree_.ETLrecHits_energy->push_back(energy);
    outTree_.ETLrecHits_time->push_back(time);
    outTree_.ETLrecHits_rr->push_back(RR);
    outTree_.ETLrecHits_module->push_back(module);
    outTree_.ETLrecHits_modType->push_back(modType);
    outTree_.ETLrecHits_crystal->push_back(crystal);
    outTree_.ETLrecHits_ieta->push_back(ieta);
    outTree_.ETLrecHits_iphi->push_back(iphi);
    outTree_.ETLrecHits_local_x->push_back(lp.x());
    outTree_.ETLrecHits_local_y->push_back(lp.y());
    outTree_.ETLrecHits_local_z->push_back(lp.z());
    outTree_.ETLrecHits_global_R->push_back(sqrt(gp.perp2()));
  }


  //---fill the tree - recHits
  outTree_.ETLclusters_n = 0;
  
  for (auto clusIt : clustersETL)
    {    
      DetId id = clusIt.detId();
      const auto& det = mtdGeometry_ -> idToDet(id);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      for ( auto cluster : clusIt)
	{
	  float energy = cluster.energy();
	  float time   = cluster.time();
	  float x=cluster.x();
	  float y=cluster.y();
	  int size=cluster.size();
	  int sizeX=cluster.sizeX();
	  int sizeY=cluster.sizeY();
	  float seed_energy=cluster.seed().energy();
	  float seed_time=cluster.seed().time();
	  int seed_x=cluster.seed().x();
	  int seed_y=cluster.seed().y();
	  
	  MeasurementPoint mp(cluster.x(),cluster.y());
	  LocalPoint lp = topo.localPosition(mp);
	  GlobalPoint gp = det->toGlobal(lp);
	  
	  MTDDetId mtdId(id);
	  int RR = 0;
	  int module = 0;
	  int modType = 0;
	  int crystal = 0;
	  int ieta = 0;
	  int iphi = 0;
	  
	  if ( mtdId.mtdSubDetector() == MTDDetId::ETL )
	    {
	      ETLDetId btlId(id);
	      RR = btlId.mtdRR();
	      module = btlId.module();
	      modType = btlId.modType();
	      // crystal = btlId.crystal();
	      // ieta = btlId.ieta();
	      // iphi = btlId.iphi();
	    }
	  
	  outTree_.ETLclusters_n += 1;    
	  outTree_.ETLclusters_size->push_back(size);
	  outTree_.ETLclusters_size_x->push_back(sizeX);
	  outTree_.ETLclusters_size_y->push_back(sizeY);
	  outTree_.ETLclusters_energy->push_back(energy);
	  outTree_.ETLclusters_time->push_back(time);
	  outTree_.ETLclusters_rr->push_back(RR);
	  outTree_.ETLclusters_module->push_back(module);
	  outTree_.ETLclusters_modType->push_back(modType);
	  outTree_.ETLclusters_crystal->push_back(crystal);
	  outTree_.ETLclusters_ieta->push_back(ieta);
	  outTree_.ETLclusters_iphi->push_back(iphi);
	  outTree_.ETLclusters_x->push_back(x);
	  outTree_.ETLclusters_y->push_back(y);
	  outTree_.ETLclusters_seed_energy->push_back(seed_energy);
	  outTree_.ETLclusters_seed_time->push_back(seed_time);
	  outTree_.ETLclusters_seed_x->push_back(seed_x);
	  outTree_.ETLclusters_seed_y->push_back(seed_y);
	  outTree_.ETLclusters_local_x->push_back(lp.x());
	  outTree_.ETLclusters_local_y->push_back(lp.y());
	  outTree_.ETLclusters_local_z->push_back(lp.z());
	  outTree_.ETLclusters_global_R->push_back(sqrt(gp.perp2()));
	}
    }

  //--- fill the tree - tracks
  int idx=0;
  for(unsigned iTrack = 0; iTrack < tracks.size(); ++iTrack)
  {
    auto track = tracks.at(iTrack);
    // skip neutrals
    if( track.charge() == 0 ) continue;
    //    if( fabs(track.eta()) > 1.5 ) continue;
    if( track.pt() < 0.7 ) continue;
    
    // match with gen particles
    float DRMin = 999999;
    int genPdgId = 0;
    float genEta = -999.;
    float genPhi = -999.;
    float genPt = -999.;
    for(auto& genPart : genParticles)
    {
      if( genPart.status() != 1 ) continue;
      if( genPart.charge() == 0 ) continue;
      
      float Deta = track.eta()-genPart.eta();
      float Dphi = deltaPhi(track.phi(),genPart.phi());
      float DR   = sqrt(Deta*Deta+Dphi*Dphi);
      
      if( DR < DRMin )
      {
        DRMin = DR;
        
        genPdgId = genPart.pdgId();
        genEta   = genPart.eta();
        genPhi   = genPart.phi();
        genPt    = genPart.pt();
      }
    }
    
    if (genPV)
      {
        outTree_.track_mcMatch_genVtx_x -> push_back(genPV->position().x());
        outTree_.track_mcMatch_genVtx_y -> push_back(genPV->position().y());
        outTree_.track_mcMatch_genVtx_z -> push_back(genPV->position().z());
        outTree_.track_mcMatch_genVtx_t -> push_back(genPV->position().t()*1E9); //ns
      }
    else
      {
        outTree_.track_mcMatch_genVtx_x -> push_back(-999.);
        outTree_.track_mcMatch_genVtx_y -> push_back(-999.);
        outTree_.track_mcMatch_genVtx_z -> push_back(-999.);
        outTree_.track_mcMatch_genVtx_t -> push_back(-999.);
      }
    
    outTree_.track_idx -> push_back(idx);
    outTree_.track_pt -> push_back(track.pt());
    outTree_.track_eta -> push_back(track.eta());
    outTree_.track_phi -> push_back(track.phi());
    outTree_.track_x -> push_back(track.vx());
    outTree_.track_y -> push_back(track.vy());
    outTree_.track_z -> push_back(track.vz());
    outTree_.track_t -> push_back(track.t0());
    outTree_.track_energy -> push_back(sqrt(track.momentum().mag2()));
    outTree_.track_normalizedChi2 -> push_back(track.normalizedChi2());
    outTree_.track_numberOfValidHits -> push_back(track.numberOfValidHits());
    outTree_.track_numberOfLostHits -> push_back(track.numberOfLostHits());
    outTree_.track_isHighPurity -> push_back(track.quality(reco::TrackBase::TrackQuality::highPurity));
    outTree_.track_hasMTD -> push_back(track.isTimeOk());
    outTree_.track_mcMatch_genPdgId -> push_back(genPdgId);
    outTree_.track_mcMatch_genPt -> push_back(genPt);
    outTree_.track_mcMatch_genEta -> push_back(genEta);
    outTree_.track_mcMatch_genPhi -> push_back(genPhi);
    outTree_.track_mcMatch_DR -> push_back(DRMin);
        
    if( verbosity_ ) std::cout << "*** track " << iTrack << " / " << tracks.size() << "   pt: " << track.pt() << "   eta: " << track.eta() << "   phi: " << track.phi() << std::endl;
    if( verbosity_ ) std::cout << "*** match with gen particle   DR: " << DRMin << "   gen pdgId: " << genPdgId << "   gen eta: " << genEta << "   gen phi: " << genPhi << "   genPt: " << genPt << std::endl;
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    outTree_.BTLmatchedSimHits_n->resize(idx+1);
    outTree_.BTLmatchedRecHits_n->resize(idx+1);
    outTree_.BTLmatchedClusters_n->resize(idx+1);
    outTree_.BTLmatchedSimHits_idx->resize(idx+1);
    outTree_.BTLmatchedSimHits_energy->resize(idx+1);
    outTree_.BTLmatchedSimHits_energyCorr->resize(idx+1);
    outTree_.BTLmatchedSimHits_time->resize(idx+1);
    outTree_.BTLmatchedSimHits_rr->resize(idx+1);
    outTree_.BTLmatchedSimHits_module->resize(idx+1);
    outTree_.BTLmatchedSimHits_modType->resize(idx+1);
    outTree_.BTLmatchedSimHits_crystal->resize(idx+1);
    outTree_.BTLmatchedSimHits_ieta->resize(idx+1);
    outTree_.BTLmatchedSimHits_iphi->resize(idx+1);
    outTree_.BTLmatchedSimHits_entry_local_x->resize(idx+1);
    outTree_.BTLmatchedSimHits_entry_local_y->resize(idx+1);
    outTree_.BTLmatchedSimHits_entry_local_z->resize(idx+1);
    outTree_.BTLmatchedSimHits_entry_global_R->resize(idx+1);
    outTree_.BTLmatchedSimHits_exit_local_x->resize(idx+1);
    outTree_.BTLmatchedSimHits_exit_local_y->resize(idx+1);
    outTree_.BTLmatchedSimHits_exit_local_z->resize(idx+1);
    outTree_.BTLmatchedSimHits_exit_global_R->resize(idx+1);
    outTree_.BTLmatchedSimHits_track_Deta->resize(idx+1);
    outTree_.BTLmatchedSimHits_track_Dphi->resize(idx+1);
    outTree_.BTLmatchedSimHits_track_DR->resize(idx+1);
    outTree_.BTLmatchedSimHits_track_Dz->resize(idx+1);
    outTree_.BTLmatchedSimHits_track_RDphi->resize(idx+1);
    outTree_.BTLmatchedSimHits_track_dist->resize(idx+1);
    outTree_.BTLmatchedRecHits_idx->resize(idx+1);
    outTree_.BTLmatchedRecHits_energy->resize(idx+1);
    outTree_.BTLmatchedRecHits_energyCorr->resize(idx+1);
    outTree_.BTLmatchedRecHits_time->resize(idx+1);
    outTree_.BTLmatchedRecHits_rr->resize(idx+1);
    outTree_.BTLmatchedRecHits_module->resize(idx+1);
    outTree_.BTLmatchedRecHits_modType->resize(idx+1);
    outTree_.BTLmatchedRecHits_crystal->resize(idx+1);
    outTree_.BTLmatchedRecHits_ieta->resize(idx+1);
    outTree_.BTLmatchedRecHits_iphi->resize(idx+1);
    outTree_.BTLmatchedRecHits_local_x->resize(idx+1);
    outTree_.BTLmatchedRecHits_local_y->resize(idx+1);
    outTree_.BTLmatchedRecHits_local_z->resize(idx+1);
    outTree_.BTLmatchedRecHits_global_R->resize(idx+1);
    outTree_.BTLmatchedRecHits_track_Deta->resize(idx+1);
    outTree_.BTLmatchedRecHits_track_Dphi->resize(idx+1);
    outTree_.BTLmatchedRecHits_track_DR->resize(idx+1);
    outTree_.BTLmatchedRecHits_track_Dz->resize(idx+1);
    outTree_.BTLmatchedRecHits_track_RDphi->resize(idx+1);
    outTree_.BTLmatchedRecHits_track_dist->resize(idx+1);
    outTree_.BTLmatchedRecHits_sietaieta->resize(idx+1);
    outTree_.BTLmatchedRecHits_siphiiphi->resize(idx+1);
    outTree_.BTLmatchedClusters_idx->resize(idx+1);
    outTree_.BTLmatchedClusters_energy->resize(idx+1);
    outTree_.BTLmatchedClusters_energyCorr->resize(idx+1);
    outTree_.BTLmatchedClusters_time->resize(idx+1);
    outTree_.BTLmatchedClusters_rr->resize(idx+1);
    outTree_.BTLmatchedClusters_module->resize(idx+1);
    outTree_.BTLmatchedClusters_modType->resize(idx+1);
    outTree_.BTLmatchedClusters_crystal->resize(idx+1);
    outTree_.BTLmatchedClusters_ieta->resize(idx+1);
    outTree_.BTLmatchedClusters_iphi->resize(idx+1);
    outTree_.BTLmatchedClusters_size->resize(idx+1);
    outTree_.BTLmatchedClusters_size_x->resize(idx+1);
    outTree_.BTLmatchedClusters_size_y->resize(idx+1);
    outTree_.BTLmatchedClusters_local_x->resize(idx+1);
    outTree_.BTLmatchedClusters_local_y->resize(idx+1);
    outTree_.BTLmatchedClusters_local_z->resize(idx+1);
    outTree_.BTLmatchedClusters_global_R->resize(idx+1);
    outTree_.BTLmatchedClusters_track_Deta->resize(idx+1);
    outTree_.BTLmatchedClusters_track_Dphi->resize(idx+1);
    outTree_.BTLmatchedClusters_track_DR->resize(idx+1);
    outTree_.BTLmatchedClusters_track_Dz->resize(idx+1);
    outTree_.BTLmatchedClusters_track_RDphi->resize(idx+1);
    outTree_.BTLmatchedClusters_track_dist->resize(idx+1);

    outTree_.ETLmatchedSimHits_n->resize(idx+1);
    outTree_.ETLmatchedRecHits_n->resize(idx+1);
    outTree_.ETLmatchedClusters_n->resize(idx+1);
    outTree_.ETLmatchedSimHits_idx->resize(idx+1);
    outTree_.ETLmatchedSimHits_energy->resize(idx+1);
    outTree_.ETLmatchedSimHits_energyCorr->resize(idx+1);
    outTree_.ETLmatchedSimHits_time->resize(idx+1);
    outTree_.ETLmatchedSimHits_rr->resize(idx+1);
    outTree_.ETLmatchedSimHits_module->resize(idx+1);
    outTree_.ETLmatchedSimHits_modType->resize(idx+1);
    outTree_.ETLmatchedSimHits_crystal->resize(idx+1);
    outTree_.ETLmatchedSimHits_ieta->resize(idx+1);
    outTree_.ETLmatchedSimHits_iphi->resize(idx+1);
    outTree_.ETLmatchedSimHits_entry_local_x->resize(idx+1);
    outTree_.ETLmatchedSimHits_entry_local_y->resize(idx+1);
    outTree_.ETLmatchedSimHits_entry_local_z->resize(idx+1);
    outTree_.ETLmatchedSimHits_entry_global_R->resize(idx+1);
    outTree_.ETLmatchedSimHits_exit_local_x->resize(idx+1);
    outTree_.ETLmatchedSimHits_exit_local_y->resize(idx+1);
    outTree_.ETLmatchedSimHits_exit_local_z->resize(idx+1);
    outTree_.ETLmatchedSimHits_exit_global_R->resize(idx+1);
    outTree_.ETLmatchedSimHits_track_Deta->resize(idx+1);
    outTree_.ETLmatchedSimHits_track_Dphi->resize(idx+1);
    outTree_.ETLmatchedSimHits_track_DR->resize(idx+1);
    outTree_.ETLmatchedSimHits_track_Dz->resize(idx+1);
    outTree_.ETLmatchedSimHits_track_RDphi->resize(idx+1);
    outTree_.ETLmatchedSimHits_track_dist->resize(idx+1);
    outTree_.ETLmatchedRecHits_idx->resize(idx+1);
    outTree_.ETLmatchedRecHits_energy->resize(idx+1);
    outTree_.ETLmatchedRecHits_energyCorr->resize(idx+1);
    outTree_.ETLmatchedRecHits_time->resize(idx+1);
    outTree_.ETLmatchedRecHits_rr->resize(idx+1);
    outTree_.ETLmatchedRecHits_module->resize(idx+1);
    outTree_.ETLmatchedRecHits_modType->resize(idx+1);
    outTree_.ETLmatchedRecHits_crystal->resize(idx+1);
    outTree_.ETLmatchedRecHits_ieta->resize(idx+1);
    outTree_.ETLmatchedRecHits_iphi->resize(idx+1);
    outTree_.ETLmatchedRecHits_local_x->resize(idx+1);
    outTree_.ETLmatchedRecHits_local_y->resize(idx+1);
    outTree_.ETLmatchedRecHits_local_z->resize(idx+1);
    outTree_.ETLmatchedRecHits_global_R->resize(idx+1);
    outTree_.ETLmatchedRecHits_track_Deta->resize(idx+1);
    outTree_.ETLmatchedRecHits_track_Dphi->resize(idx+1);
    outTree_.ETLmatchedRecHits_track_DR->resize(idx+1);
    outTree_.ETLmatchedRecHits_track_Dz->resize(idx+1);
    outTree_.ETLmatchedRecHits_track_RDphi->resize(idx+1);
    outTree_.ETLmatchedRecHits_track_dist->resize(idx+1);
    outTree_.ETLmatchedRecHits_sietaieta->resize(idx+1);
    outTree_.ETLmatchedRecHits_siphiiphi->resize(idx+1);
    outTree_.ETLmatchedClusters_idx->resize(idx+1);
    outTree_.ETLmatchedClusters_energy->resize(idx+1);
    outTree_.ETLmatchedClusters_energyCorr->resize(idx+1);
    outTree_.ETLmatchedClusters_time->resize(idx+1);
    outTree_.ETLmatchedClusters_rr->resize(idx+1);
    outTree_.ETLmatchedClusters_module->resize(idx+1);
    outTree_.ETLmatchedClusters_modType->resize(idx+1);
    outTree_.ETLmatchedClusters_crystal->resize(idx+1);
    outTree_.ETLmatchedClusters_ieta->resize(idx+1);
    outTree_.ETLmatchedClusters_iphi->resize(idx+1);
    outTree_.ETLmatchedClusters_size->resize(idx+1);
    outTree_.ETLmatchedClusters_size_x->resize(idx+1);
    outTree_.ETLmatchedClusters_size_y->resize(idx+1);
    outTree_.ETLmatchedClusters_local_x->resize(idx+1);
    outTree_.ETLmatchedClusters_local_y->resize(idx+1);
    outTree_.ETLmatchedClusters_local_z->resize(idx+1);
    outTree_.ETLmatchedClusters_global_R->resize(idx+1);
    outTree_.ETLmatchedClusters_track_Deta->resize(idx+1);
    outTree_.ETLmatchedClusters_track_Dphi->resize(idx+1);
    outTree_.ETLmatchedClusters_track_DR->resize(idx+1);
    outTree_.ETLmatchedClusters_track_Dz->resize(idx+1);
    outTree_.ETLmatchedClusters_track_RDphi->resize(idx+1);
    outTree_.ETLmatchedClusters_track_dist->resize(idx+1);
        
    //---get compatible layers/Dets
    std::vector<GlobalPoint> gp_ext;
    auto tTrack = builder->build(track);
    TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
    float theMaxChi2 = 25.;
    float theNSigma = 5.;
    std::unique_ptr<MeasurementEstimator> theEstimator = std::make_unique<Chi2MeasurementEstimator>(theMaxChi2,theNSigma);
    SteppingHelixPropagator prop(theField.product(),anyDirection);

    //try BTL
    bool inBTL=false;
    const vector<const DetLayer*>& layersBTL = layerGeo->allBTLLayers();
    for (const DetLayer* ilay : layersBTL) 
      {
	pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos,prop,*theEstimator);	
	if (!comp.first) continue;
	vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos,prop,*theEstimator);
	for( const auto& detWithState : compDets ) 
	  {
	    gp_ext.push_back(detWithState.second.globalPosition());
	    if (!inBTL)
	      {
		outTree_.track_eta_atBTL -> push_back(gp_ext.back().eta());
		outTree_.track_phi_atBTL -> push_back(gp_ext.back().phi());
		inBTL=true;
	      }
	  }
      }

    //try ETL
    bool inETL=false;
    const vector<const DetLayer*>& layersETL = layerGeo->allETLLayers();
    for (const DetLayer* ilay : layersETL) 
      {
	const BoundDisk& disk = static_cast<const MTDRingForwardDoubleLayer*>(ilay)->specificSurface();
	const double diskZ = disk.position().z();
	if( tsos.globalPosition().z() * diskZ < 0 ) continue; // only propagate to the disk that's on the same side
	pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos,prop,*theEstimator);	
	if (!comp.first) continue;
	vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos,prop,*theEstimator);
	for( const auto& detWithState : compDets ) 
	  {
	    gp_ext.push_back(detWithState.second.globalPosition());
	    if (!inETL)
	      {
		outTree_.track_eta_atETL -> push_back(gp_ext.back().eta());
		outTree_.track_phi_atETL -> push_back(gp_ext.back().phi());
		inETL=true;
	      }
	  }
      }

    if( !inBTL )
      {
	outTree_.track_eta_atBTL -> push_back(-999.);
	outTree_.track_phi_atBTL -> push_back(-999.);
      }

    if( !inETL )
      {
	outTree_.track_eta_atETL -> push_back(-999.);
	outTree_.track_phi_atETL -> push_back(-999.);
      }

    {
    //---get associated simHits
    if( verbosity_ ) std::cout << "*** simHits - n tot: " << simHitsBTL.size() << std::endl;
    int simHitIt = 0;
    for(auto simHit : simHitsBTL)
    {
      BTLDetId id = simHit.detUnitId();
      DetId geoId = id.geographicalId( crysLayout_ );
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      
      double energy = simHit.energyLoss()*1000.;
      double time   = simHit.tof();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = id.crystal();
      int ieta = id.ieta(crysLayout_);
      int iphi = id.iphi(crysLayout_);
      
      LocalPoint lp_entry(   simHit.entryPoint().x()/10.,   simHit.entryPoint().y()/10.,   simHit.entryPoint().z()/10.);
      LocalPoint lp_mid  (simHit.localPosition().x()/10.,simHit.localPosition().y()/10.,simHit.localPosition().z()/10.);
      LocalPoint lp_exit (    simHit.exitPoint().x()/10.,    simHit.exitPoint().y()/10.,    simHit.exitPoint().z()/10.);
      GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,id.row(topo.nrows()),id.column(topo.nrows())));
      GlobalPoint gp_mid   = det->toGlobal(topo.pixelToModuleLocalPoint(lp_mid,id.row(topo.nrows()),id.column(topo.nrows())));
      GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit,id.row(topo.nrows()),id.column(topo.nrows())));
      
      float eta = gp_mid.eta();
      float phi = gp_mid.phi();

      int closestPoint=-1;
      float minDist=999;
      for (unsigned int ic=0; ic<gp_ext.size();++ic)
	{
	  GlobalPoint diff(gp_ext[ic].x()-gp_mid.x(),gp_ext[ic].y()-gp_mid.y(),gp_ext[ic].z()-gp_mid.z());
	  if (diff.mag()<minDist)
	    {
	      closestPoint=ic;
	      minDist=diff.mag();
	    }
	}

      if (closestPoint == -1)
	continue;

      GlobalPoint gp_track = gp_ext[closestPoint];      
      //GlobalPoint gp_track = gp_ext[abs(id.row(topo.nrows())-int())];

      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp_mid.z()      : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp_mid-gp_track).mag()      : -999.;
      if( DR < 0.05 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
                                   << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << id.row(topo.nrows()) << "   column: " << id.column(topo.nrows()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   entryPoint:   local: " << PrintPosition(lp_entry)        << "   global: " << PrintPosition(gp_entry) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":     midPoint:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid)   << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":    exitPoint:   local: " << PrintPosition(lp_exit)         << "   global: " << PrintPosition(gp_exit)  << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ": hit position:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid) << "   DR: " << DR << "   dist: " << (gp_track-gp_mid).mag() << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.BTLmatchedSimHits_n->at(idx) += 1;
      
      outTree_.BTLmatchedSimHits_idx->at(idx).push_back(idx);
      outTree_.BTLmatchedSimHits_energy->at(idx).push_back(energy);
      outTree_.BTLmatchedSimHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.BTLmatchedSimHits_time->at(idx).push_back(time);
      outTree_.BTLmatchedSimHits_rr->at(idx).push_back(RR);
      outTree_.BTLmatchedSimHits_module->at(idx).push_back(module);
      outTree_.BTLmatchedSimHits_modType->at(idx).push_back(modType);
      outTree_.BTLmatchedSimHits_crystal->at(idx).push_back(crystal);
      outTree_.BTLmatchedSimHits_ieta->at(idx).push_back(ieta);
      outTree_.BTLmatchedSimHits_iphi->at(idx).push_back(iphi);
      outTree_.BTLmatchedSimHits_entry_local_x->at(idx).push_back(lp_entry.x());
      outTree_.BTLmatchedSimHits_entry_local_y->at(idx).push_back(lp_entry.y());
      outTree_.BTLmatchedSimHits_entry_local_z->at(idx).push_back(lp_entry.z());
      outTree_.BTLmatchedSimHits_entry_global_R->at(idx).push_back(sqrt(gp_entry.perp2()));
      outTree_.BTLmatchedSimHits_exit_local_x->at(idx).push_back(lp_exit.x());
      outTree_.BTLmatchedSimHits_exit_local_y->at(idx).push_back(lp_exit.y());
      outTree_.BTLmatchedSimHits_exit_local_z->at(idx).push_back(lp_exit.z());
      outTree_.BTLmatchedSimHits_exit_global_R->at(idx).push_back(sqrt(gp_exit.perp2()));
      outTree_.BTLmatchedSimHits_track_Deta->at(idx).push_back(fabs(Deta));
      outTree_.BTLmatchedSimHits_track_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.BTLmatchedSimHits_track_DR->at(idx).push_back(DR);
      outTree_.BTLmatchedSimHits_track_Dz->at(idx).push_back(Dz);
      outTree_.BTLmatchedSimHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.BTLmatchedSimHits_track_dist->at(idx).push_back(dist);
      
      ++simHitIt;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    
    //---find associated recHits
    float sieie=0, sipip=0;
    float ss_hit_count=0;
    int recHitIt = 0;
    if( verbosity_ ) std::cout << "*** recHits - tot: " << recHitsBTL.size() << std::endl;
    for(auto recHit : recHitsBTL)
    {
      BTLDetId id = recHit.id();
      DetId geoId = id.geographicalId( crysLayout_ );
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      
      double energy = recHit.energy();
      double time   = recHit.time();
      
      MeasurementPoint mp(recHit.row(),recHit.column());
      LocalPoint lp = topo.localPosition(mp);
      GlobalPoint gp = det->toGlobal(lp);

      float eta = gp.eta();
      float phi = gp.phi();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = id.crystal();
      int ieta = id.ieta(crysLayout_);
      int iphi = id.iphi(crysLayout_);
      
      int closestPoint=-1;
      float minDist=999;
      for (unsigned int ic=0; ic<gp_ext.size();++ic)
	{
	  GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
	  if (diff.mag()<minDist)
	    {
	      closestPoint=ic;
	      minDist=diff.mag();
	    }
	}
      
      if (closestPoint == -1)
	continue;

      GlobalPoint gp_track = gp_ext[closestPoint];      
      //      GlobalPoint gp_track = gp_ext[abs(id.row(topo.nrows())-int())];
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
      
      // if( DR < 0.2 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
                                   << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << recHit.row() << "- " << id.row(topo.nrows()) << "   column: " << recHit.column() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << "   DR: " << DR << "   dist: " << (gp_track-det->position()).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": hit position:   local: " << PrintPosition(lp) << "   global: " << PrintPosition(gp) << "DR: " << DR << "   dist: " << (gp_track-gp).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": extrapolated track " << abs(recHit.row()-int()) << " : " << PrintPosition(gp_track) << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.BTLmatchedRecHits_n->at(idx) += 1;
      
      outTree_.BTLmatchedRecHits_idx->at(idx).push_back(idx);
      outTree_.BTLmatchedRecHits_energy->at(idx).push_back(energy);
      outTree_.BTLmatchedRecHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.BTLmatchedRecHits_time->at(idx).push_back(time);
      outTree_.BTLmatchedRecHits_rr->at(idx).push_back(RR);
      outTree_.BTLmatchedRecHits_module->at(idx).push_back(module);
      outTree_.BTLmatchedRecHits_modType->at(idx).push_back(modType);
      outTree_.BTLmatchedRecHits_crystal->at(idx).push_back(crystal);
      outTree_.BTLmatchedRecHits_ieta->at(idx).push_back(ieta);
      outTree_.BTLmatchedRecHits_iphi->at(idx).push_back(iphi);
      outTree_.BTLmatchedRecHits_local_x->at(idx).push_back(lp.x());
      outTree_.BTLmatchedRecHits_local_y->at(idx).push_back(lp.y());
      outTree_.BTLmatchedRecHits_local_z->at(idx).push_back(lp.z());
      outTree_.BTLmatchedRecHits_global_R->at(idx).push_back(sqrt(gp.perp2()));
      outTree_.BTLmatchedRecHits_track_Deta->at(idx).push_back(fabs(Deta));
      outTree_.BTLmatchedRecHits_track_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.BTLmatchedRecHits_track_Dz->at(idx).push_back(Dz);
      outTree_.BTLmatchedRecHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.BTLmatchedRecHits_track_DR->at(idx).push_back(DR);
      outTree_.BTLmatchedRecHits_track_dist->at(idx).push_back(dist);
      
      if(recHit.energy() > 0.5)
      {
        sieie += energy*pow(eta-gp_track.eta(),2);
        sipip += energy*pow(phi-gp_track.phi(),2);
        ss_hit_count += energy;
      }
      
      ++recHitIt;
    }
  
    //---find associated recHits
    int clusterIt = 0;

    for (auto clusIt : clustersBTL)
      {    
	DetId id = clusIt.detId();
	const auto& det = mtdGeometry_ -> idToDet(id);
	const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
	const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
	for ( auto cluster : clusIt)
	  {
	    
	    MTDDetId mtdId(id);
	    int RR = 0;
	    int module = 0;
	    int modType = 0;
	    int crystal = 0;
	    int ieta = 0;
	    int iphi = 0;
	    
	    if ( mtdId.mtdSubDetector() == MTDDetId::BTL )
	      {
		BTLDetId btlId(id);
		RR = btlId.mtdRR();
		module = btlId.module();
		modType = btlId.modType();
		crystal = btlId.crystal();
		ieta = btlId.ieta(crysLayout_);
		iphi = btlId.iphi(crysLayout_);
	      }
	    
	    double energy = cluster.energy();
	    double time   = cluster.time();
	    int size=cluster.size();
	    int sizeX=cluster.sizeX();
	    int sizeY=cluster.sizeY();
	    
	    MeasurementPoint mp(cluster.x(),cluster.y());
	    LocalPoint lp = topo.localPosition(mp);
	    GlobalPoint gp = det->toGlobal(lp);
	    
	    float eta = gp.eta();
	    float phi = gp.phi();
	    
	    int closestPoint=-1;
	    float minDist=999;
	    for (unsigned int ic=0; ic<gp_ext.size();++ic)
	      {
		GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
		if (diff.mag()<minDist)
		  {
		    closestPoint=ic;
		    minDist=diff.mag();
		  }
	      }
	    
	    if (closestPoint == -1)
	      continue;

	    GlobalPoint gp_track = gp_ext[closestPoint];      
	    //	    GlobalPoint gp_track = gp_ext[abs(cluster.seed().x-int())];
	    
	    float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
	    float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
	    float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
	    float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
	    float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
	    float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
	    	    
	    
	    if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
	    if( gp_track.mag() <= 0. ) continue;
	    
	    outTree_.BTLmatchedClusters_n->at(idx) += 1;
	    
	    outTree_.BTLmatchedClusters_idx->at(idx).push_back(idx);
	    outTree_.BTLmatchedClusters_energy->at(idx).push_back(energy);
	    outTree_.BTLmatchedClusters_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
	    outTree_.BTLmatchedClusters_time->at(idx).push_back(time);
	    outTree_.BTLmatchedClusters_rr->at(idx).push_back(RR);
	    outTree_.BTLmatchedClusters_module->at(idx).push_back(module);
	    outTree_.BTLmatchedClusters_modType->at(idx).push_back(modType);
	    outTree_.BTLmatchedClusters_crystal->at(idx).push_back(crystal);
	    outTree_.BTLmatchedClusters_ieta->at(idx).push_back(ieta);
	    outTree_.BTLmatchedClusters_iphi->at(idx).push_back(iphi);
	    outTree_.BTLmatchedClusters_size->at(idx).push_back(size);
	    outTree_.BTLmatchedClusters_size_x->at(idx).push_back(sizeX);
	    outTree_.BTLmatchedClusters_size_y->at(idx).push_back(sizeY);
	    outTree_.BTLmatchedClusters_local_x->at(idx).push_back(lp.x());
	    outTree_.BTLmatchedClusters_local_y->at(idx).push_back(lp.y());
	    outTree_.BTLmatchedClusters_local_z->at(idx).push_back(lp.z());
	    outTree_.BTLmatchedClusters_global_R->at(idx).push_back(sqrt(gp.perp2()));
	    outTree_.BTLmatchedClusters_track_Deta->at(idx).push_back(fabs(Deta));
	    outTree_.BTLmatchedClusters_track_Dphi->at(idx).push_back(fabs(Dphi));
	    outTree_.BTLmatchedClusters_track_Dz->at(idx).push_back(Dz);
	    outTree_.BTLmatchedClusters_track_RDphi->at(idx).push_back(RDphi);
	    outTree_.BTLmatchedClusters_track_DR->at(idx).push_back(DR);
	    outTree_.BTLmatchedClusters_track_dist->at(idx).push_back(dist);
	    
	    ++clusterIt;
	  }
      }
  
  
    if( verbosity_ ) std::cout << "---\n\n\n" << std::endl;
    
    outTree_.BTLmatchedRecHits_sietaieta->at(idx).push_back( ss_hit_count>0 ? sqrt(sieie)/ss_hit_count : -999. );
    outTree_.BTLmatchedRecHits_siphiiphi->at(idx).push_back( ss_hit_count>0 ? sqrt(sipip)/ss_hit_count : -999. );

    }
    ///ETL
    {
    //---get associated simHits
    if( verbosity_ ) std::cout << "*** simHits - n tot: " << simHitsETL.size() << std::endl;
    int simHitIt = 0;
    for(auto simHit : simHitsETL)
    {
      ETLDetId id = simHit.detUnitId();
      DetId geoId = id.geographicalId(  );
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      
      double energy = simHit.energyLoss()*1000.;
      double time   = simHit.tof();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = 0;
      int ieta = 0;
      int iphi = 0;

      // ETL is already in module-local coordinates so just scale to cm from mm
      Local3DPoint simscaled(0.1*simHit.entryPoint().x(),0.1*simHit.entryPoint().y(),0.1*simHit.entryPoint().z());
      const auto& thepixel = topo.pixel(simscaled); // mm -> cm here is the switch                               
      const uint8_t row(thepixel.first), col(thepixel.second);
      
      LocalPoint lp_entry(   simHit.entryPoint().x()/10.,   simHit.entryPoint().y()/10.,   simHit.entryPoint().z()/10.);
      LocalPoint lp_mid  (simHit.localPosition().x()/10.,simHit.localPosition().y()/10.,simHit.localPosition().z()/10.);
      LocalPoint lp_exit (    simHit.exitPoint().x()/10.,    simHit.exitPoint().y()/10.,    simHit.exitPoint().z()/10.);
      GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry ,row,col));
      GlobalPoint gp_mid   = det->toGlobal(topo.pixelToModuleLocalPoint(lp_mid   ,row,col));
      GlobalPoint gp_exit  = det->toGlobal(topo.pixelToModuleLocalPoint(lp_exit  ,row,col));
      
      float eta = gp_mid.eta();
      float phi = gp_mid.phi();

      int closestPoint=-1;
      float minDist=999;
      for (unsigned int ic=0; ic<gp_ext.size();++ic)
	{
	  GlobalPoint diff(gp_ext[ic].x()-gp_mid.x(),gp_ext[ic].y()-gp_mid.y(),gp_ext[ic].z()-gp_mid.z());
	  if (diff.mag()<minDist)
	    {
	      closestPoint=ic;
	      minDist=diff.mag();
	    }
	}

      if (closestPoint == -1)
	continue;

      GlobalPoint gp_track = gp_ext[closestPoint];      

      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp_mid.z()      : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp_mid-gp_track).mag()      : -999.;
      if( DR < 0.05 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":   entryPoint:   local: " << PrintPosition(lp_entry)        << "   global: " << PrintPosition(gp_entry) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":     midPoint:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid)   << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":    exitPoint:   local: " << PrintPosition(lp_exit)         << "   global: " << PrintPosition(gp_exit)  << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << std::endl;
        if( verbosity_ ) std::cout << ">>> " << simHitIt << ": hit position:   local: " << PrintPosition(lp_mid)          << "   global: " << PrintPosition(gp_mid) << "   DR: " << DR << "   dist: " << (gp_track-gp_mid).mag() << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.ETLmatchedSimHits_n->at(idx) += 1;
      
      outTree_.ETLmatchedSimHits_idx->at(idx).push_back(idx);
      outTree_.ETLmatchedSimHits_energy->at(idx).push_back(energy);
      outTree_.ETLmatchedSimHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.ETLmatchedSimHits_time->at(idx).push_back(time);
      outTree_.ETLmatchedSimHits_rr->at(idx).push_back(RR);
      outTree_.ETLmatchedSimHits_module->at(idx).push_back(module);
      outTree_.ETLmatchedSimHits_modType->at(idx).push_back(modType);
      outTree_.ETLmatchedSimHits_crystal->at(idx).push_back(crystal);
      outTree_.ETLmatchedSimHits_ieta->at(idx).push_back(ieta);
      outTree_.ETLmatchedSimHits_iphi->at(idx).push_back(iphi);
      outTree_.ETLmatchedSimHits_entry_local_x->at(idx).push_back(lp_entry.x());
      outTree_.ETLmatchedSimHits_entry_local_y->at(idx).push_back(lp_entry.y());
      outTree_.ETLmatchedSimHits_entry_local_z->at(idx).push_back(lp_entry.z());
      outTree_.ETLmatchedSimHits_entry_global_R->at(idx).push_back(sqrt(gp_entry.perp2()));
      outTree_.ETLmatchedSimHits_exit_local_x->at(idx).push_back(lp_exit.x());
      outTree_.ETLmatchedSimHits_exit_local_y->at(idx).push_back(lp_exit.y());
      outTree_.ETLmatchedSimHits_exit_local_z->at(idx).push_back(lp_exit.z());
      outTree_.ETLmatchedSimHits_exit_global_R->at(idx).push_back(sqrt(gp_exit.perp2()));
      outTree_.ETLmatchedSimHits_track_Deta->at(idx).push_back(fabs(Deta));
      outTree_.ETLmatchedSimHits_track_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.ETLmatchedSimHits_track_DR->at(idx).push_back(DR);
      outTree_.ETLmatchedSimHits_track_Dz->at(idx).push_back(Dz);
      outTree_.ETLmatchedSimHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.ETLmatchedSimHits_track_dist->at(idx).push_back(dist);
      
      ++simHitIt;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    
    //---find associated recHits
    float sieie=0, sipip=0;
    float ss_hit_count=0;
    int recHitIt = 0;
    if( verbosity_ ) std::cout << "*** recHits - tot: " << recHitsETL.size() << std::endl;
    for(auto recHit : recHitsETL)
    {
      ETLDetId id = recHit.id();
      DetId geoId = id.geographicalId(  );
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      
      double energy = recHit.energy();
      double time   = recHit.time();
      
      MeasurementPoint mp(recHit.row(),recHit.column());
      LocalPoint lp = topo.localPosition(mp);
      GlobalPoint gp = det->toGlobal(lp);

      float eta = gp.eta();
      float phi = gp.phi();
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = 0;
      int ieta = 0;
      int iphi = 0;
      
      int closestPoint=-1;
      float minDist=999;
      for (unsigned int ic=0; ic<gp_ext.size();++ic)
	{
	  GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
	  if (diff.mag()<minDist)
	    {
	      closestPoint=ic;
	      minDist=diff.mag();
	    }
	}
      
      if (closestPoint == -1)
	continue;

      GlobalPoint gp_track = gp_ext[closestPoint];      
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
      
      // if( DR < 0.2 && DR > 0. )
      {
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << "   DR: " << DR << "   dist: " << (gp_track-det->position()).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": hit position:   local: " << PrintPosition(lp) << "   global: " << PrintPosition(gp) << "DR: " << DR << "   dist: " << (gp_track-gp).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": extrapolated track " << abs(recHit.row()-int()) << " : " << PrintPosition(gp_track) << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.ETLmatchedRecHits_n->at(idx) += 1;
      
      outTree_.ETLmatchedRecHits_idx->at(idx).push_back(idx);
      outTree_.ETLmatchedRecHits_energy->at(idx).push_back(energy);
      outTree_.ETLmatchedRecHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.ETLmatchedRecHits_time->at(idx).push_back(time);
      outTree_.ETLmatchedRecHits_rr->at(idx).push_back(RR);
      outTree_.ETLmatchedRecHits_module->at(idx).push_back(module);
      outTree_.ETLmatchedRecHits_modType->at(idx).push_back(modType);
      outTree_.ETLmatchedRecHits_crystal->at(idx).push_back(crystal);
      outTree_.ETLmatchedRecHits_ieta->at(idx).push_back(ieta);
      outTree_.ETLmatchedRecHits_iphi->at(idx).push_back(iphi);
      outTree_.ETLmatchedRecHits_local_x->at(idx).push_back(lp.x());
      outTree_.ETLmatchedRecHits_local_y->at(idx).push_back(lp.y());
      outTree_.ETLmatchedRecHits_local_z->at(idx).push_back(lp.z());
      outTree_.ETLmatchedRecHits_global_R->at(idx).push_back(sqrt(gp.perp2()));
      outTree_.ETLmatchedRecHits_track_Deta->at(idx).push_back(fabs(Deta));
      outTree_.ETLmatchedRecHits_track_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.ETLmatchedRecHits_track_Dz->at(idx).push_back(Dz);
      outTree_.ETLmatchedRecHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.ETLmatchedRecHits_track_DR->at(idx).push_back(DR);
      outTree_.ETLmatchedRecHits_track_dist->at(idx).push_back(dist);
      
      if(recHit.energy() > 0.5)
      {
        sieie += energy*pow(eta-gp_track.eta(),2);
        sipip += energy*pow(phi-gp_track.phi(),2);
        ss_hit_count += energy;
      }
      
      ++recHitIt;
    }
  
    //---find associated recHits
    int clusterIt = 0;

    for (auto clusIt : clustersETL)
      {    
	DetId id = clusIt.detId();
	const auto& det = mtdGeometry_ -> idToDet(id);
	const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
	const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
	for ( auto cluster : clusIt)
	  {
	    
	    MTDDetId mtdId(id);
	    int RR = 0;
	    int module = 0;
	    int modType = 0;
	    int crystal = 0;
	    int ieta = 0;
	    int iphi = 0;
	    
	    if ( mtdId.mtdSubDetector() == MTDDetId::ETL )
	      {
		ETLDetId btlId(id);
		RR = btlId.mtdRR();
		module = btlId.module();
		modType = btlId.modType();
		// crystal = btlId.crystal();
		// ieta = btlId.ieta();
		// iphi = btlId.iphi();
	      }
	    
	    double energy = cluster.energy();
	    double time   = cluster.time();
	    int size=cluster.size();
	    int sizeX=cluster.sizeX();
	    int sizeY=cluster.sizeY();
	    
	    MeasurementPoint mp(cluster.x(),cluster.y());
	    LocalPoint lp = topo.localPosition(mp);
	    GlobalPoint gp = det->toGlobal(lp);
	    
	    float eta = gp.eta();
	    float phi = gp.phi();
	    
	    int closestPoint=-1;
	    float minDist=999;
	    for (unsigned int ic=0; ic<gp_ext.size();++ic)
	      {
		GlobalPoint diff(gp_ext[ic].x()-gp.x(),gp_ext[ic].y()-gp.y(),gp_ext[ic].z()-gp.z());
		if (diff.mag()<minDist)
		  {
		    closestPoint=ic;
		    minDist=diff.mag();
		  }
	      }
	    
	    if (closestPoint == -1)
	      continue;

	    GlobalPoint gp_track = gp_ext[closestPoint];      
	    
	    float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
	    float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
	    float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
	    float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
	    float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
	    float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
	    
	    if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
	    if( gp_track.mag() <= 0. ) continue;
	    
	    outTree_.ETLmatchedClusters_n->at(idx) += 1;
	    
	    outTree_.ETLmatchedClusters_idx->at(idx).push_back(idx);
	    outTree_.ETLmatchedClusters_energy->at(idx).push_back(energy);
	    outTree_.ETLmatchedClusters_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
	    outTree_.ETLmatchedClusters_time->at(idx).push_back(time);
	    outTree_.ETLmatchedClusters_rr->at(idx).push_back(RR);
	    outTree_.ETLmatchedClusters_module->at(idx).push_back(module);
	    outTree_.ETLmatchedClusters_modType->at(idx).push_back(modType);
	    outTree_.ETLmatchedClusters_crystal->at(idx).push_back(crystal);
	    outTree_.ETLmatchedClusters_ieta->at(idx).push_back(ieta);
	    outTree_.ETLmatchedClusters_iphi->at(idx).push_back(iphi);
	    outTree_.ETLmatchedClusters_size->at(idx).push_back(size);
	    outTree_.ETLmatchedClusters_size_x->at(idx).push_back(sizeX);
	    outTree_.ETLmatchedClusters_size_y->at(idx).push_back(sizeY);
	    outTree_.ETLmatchedClusters_local_x->at(idx).push_back(lp.x());
	    outTree_.ETLmatchedClusters_local_y->at(idx).push_back(lp.y());
	    outTree_.ETLmatchedClusters_local_z->at(idx).push_back(lp.z());
	    outTree_.ETLmatchedClusters_global_R->at(idx).push_back(sqrt(gp.perp2()));
	    outTree_.ETLmatchedClusters_track_Deta->at(idx).push_back(fabs(Deta));
	    outTree_.ETLmatchedClusters_track_Dphi->at(idx).push_back(fabs(Dphi));
	    outTree_.ETLmatchedClusters_track_Dz->at(idx).push_back(Dz);
	    outTree_.ETLmatchedClusters_track_RDphi->at(idx).push_back(RDphi);
	    outTree_.ETLmatchedClusters_track_DR->at(idx).push_back(DR);
	    outTree_.ETLmatchedClusters_track_dist->at(idx).push_back(dist);
	    
	    ++clusterIt;
	  }
      }
  
  
    if( verbosity_ ) std::cout << "---\n\n\n" << std::endl;
    
    outTree_.ETLmatchedRecHits_sietaieta->at(idx).push_back( ss_hit_count>0 ? sqrt(sieie)/ss_hit_count : -999. );
    outTree_.ETLmatchedRecHits_siphiiphi->at(idx).push_back( ss_hit_count>0 ? sqrt(sipip)/ss_hit_count : -999. );
    }
    ++idx;
  }
  
  outTree_.GetTTreePtr()->Fill();
}


std::string FTLDumpHits::PrintPosition(const GlobalPoint& gp)
{
  std::stringstream output;
  
  output << "(";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.x() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.y() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << gp.z();
  output << ") cm";
  
  output << "   R: " << std::fixed << std::setprecision(3) << std::setw(7) << gp.perp();
  output << " cm";
  
  output << "   eta: " << std::setprecision(3) << std::setw(6) << gp.eta(); 
  output << "   phi: " << std::setprecision(3) << std::setw(6) << gp.phi();
  
  return output.str();
}

std::string FTLDumpHits::PrintPosition(const LocalPoint& lp)
{
  std::stringstream output;
  
  output << "(";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.x() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.y() << ",";
  output << std::fixed << std::setprecision(3) << std::setw(8) << lp.z();
  output << ") cm";
  
  return output.str();
}
DEFINE_FWK_MODULE(FTLDumpHits);

#endif
