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
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
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
  edm::Handle<std::vector<PSimHit> > simHitsHandle_;
  edm::EDGetTokenT<std::vector<PSimHit> > simHitsToken_;
  edm::Handle<FTLRecHitCollection> recHitsHandle_;
  edm::EDGetTokenT<FTLRecHitCollection> recHitsToken_;    
  edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_;
  edm::Handle<edm::View<reco::Track> > tracksHandle_;
  
  //---options
  BTLDetId::CrysLayout crysLayout_;
  double track_hit_DRMax_;
  double track_hit_distMax_;
  bool verbosity_;
  
  //---outputs
  FTLHitsTree outTree_;
  edm::Service<TFileService> fs_;  
  
  //---other
  MeasurementEstimator* theEstimator;
};



FTLDumpHits::FTLDumpHits(const edm::ParameterSet& pSet):
  genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
  simHitsToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsTag"))),
  recHitsToken_(consumes<FTLRecHitCollection>(pSet.getUntrackedParameter<edm::InputTag>("recHitsTag"))),
  tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("tracksTag"))),
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
  
  //--- load gen particles
  event.getByToken(genParticlesToken_, genParticlesHandle_);
  auto genParticles = *genParticlesHandle_.product();
  
  //---load sim hits
  event.getByToken(simHitsToken_, simHitsHandle_);
  auto simHits = *simHitsHandle_.product();
  
  //---load the FTL collection if present in the EventContent (avoid crash with standard geometry)
  event.getByToken(recHitsToken_, recHitsHandle_);
  auto recHits = FTLRecHitCollection();
  if(recHitsHandle_.isValid())
    recHits = *recHitsHandle_.product();
  
  //---load tracks
  event.getByToken(tracksToken_,tracksHandle_);
  auto tracks = *tracksHandle_.product();
  
  
  
  //---fill the tree - simHits
  outTree_.simHits_n = 0;
  for(auto simHit : simHits)
  {
    BTLDetId id = simHit.detUnitId();
    DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+14*(id.modType()-1),0,1);
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
    
    outTree_.simHits_n += 1;
    
    outTree_.simHits_energy->push_back(energy);
    outTree_.simHits_time->push_back(time);
    outTree_.simHits_rr->push_back(RR);
    outTree_.simHits_module->push_back(module);
    outTree_.simHits_modType->push_back(modType);
    outTree_.simHits_crystal->push_back(crystal);
    outTree_.simHits_ieta->push_back(ieta);
    outTree_.simHits_iphi->push_back(iphi);
    outTree_.simHits_entry_local_x->push_back(lp_entry.x());
    outTree_.simHits_entry_local_y->push_back(lp_entry.y());
    outTree_.simHits_entry_local_z->push_back(lp_entry.z());
    outTree_.simHits_entry_global_R->push_back(sqrt(gp_entry.perp2()));
    outTree_.simHits_exit_local_x->push_back(lp_exit.x());
    outTree_.simHits_exit_local_y->push_back(lp_exit.y());
    outTree_.simHits_exit_local_z->push_back(lp_exit.z());
    outTree_.simHits_exit_global_R->push_back(sqrt(gp_exit.perp2()));
  }
  
  
  
  //---fill the tree - recHits
  outTree_.recHits_n = 0;
  for(auto recHit : recHits)
  {
    BTLDetId id = recHit.id();
    DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+14*(id.modType()-1),0,1);
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
    
    outTree_.recHits_n += 1;
    
    outTree_.recHits_energy->push_back(energy);
    outTree_.recHits_time->push_back(time);
    outTree_.recHits_rr->push_back(RR);
    outTree_.recHits_module->push_back(module);
    outTree_.recHits_modType->push_back(modType);
    outTree_.recHits_crystal->push_back(crystal);
    outTree_.recHits_ieta->push_back(ieta);
    outTree_.recHits_iphi->push_back(iphi);
    outTree_.recHits_local_x->push_back(lp.x());
    outTree_.recHits_local_y->push_back(lp.y());
    outTree_.recHits_local_z->push_back(lp.z());
    outTree_.recHits_global_R->push_back(sqrt(gp.perp2()));
  }
  
  
  
  //--- fill the tree - tracks
  int idx=0;
  for(unsigned iTrack = 0; iTrack < tracks.size(); ++iTrack)
  {
    auto track = tracks.at(iTrack);
    // skip neutrals
    if( track.charge() == 0 ) continue;
    if( fabs(track.eta()) > 1.5 ) continue;
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
    
    
    outTree_.track_idx -> push_back(idx);
    outTree_.track_pt -> push_back(track.pt());
    outTree_.track_eta -> push_back(track.eta());
    outTree_.track_phi -> push_back(track.phi());
    outTree_.track_energy -> push_back(sqrt(track.momentum().mag2()));
    outTree_.track_normalizedChi2 -> push_back(track.normalizedChi2());
    outTree_.track_numberOfValidHits -> push_back(track.numberOfValidHits());
    outTree_.track_numberOfLostHits -> push_back(track.numberOfLostHits());
    outTree_.track_isHighPurity -> push_back(track.quality(reco::TrackBase::TrackQuality::highPurity));
    outTree_.track_mcMatch_genPdgId -> push_back(genPdgId);
    outTree_.track_mcMatch_genPt -> push_back(genPt);
    outTree_.track_mcMatch_genEta -> push_back(genEta);
    outTree_.track_mcMatch_genPhi -> push_back(genPhi);
    outTree_.track_mcMatch_DR -> push_back(DRMin);
        
    if( verbosity_ ) std::cout << "*** track " << iTrack << " / " << tracks.size() << "   pt: " << track.pt() << "   eta: " << track.eta() << "   phi: " << track.phi() << std::endl;
    if( verbosity_ ) std::cout << "*** match with gen particle   DR: " << DRMin << "   gen pdgId: " << genPdgId << "   gen eta: " << genEta << "   gen phi: " << genPhi << "   genPt: " << genPt << std::endl;
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    outTree_.matchedSimHits_n->resize(idx+1);
    outTree_.matchedRecHits_n->resize(idx+1);
    
    outTree_.matchedSimHits_idx->resize(idx+1);
    outTree_.matchedSimHits_energy->resize(idx+1);
    outTree_.matchedSimHits_energyCorr->resize(idx+1);
    outTree_.matchedSimHits_time->resize(idx+1);
    outTree_.matchedSimHits_rr->resize(idx+1);
    outTree_.matchedSimHits_module->resize(idx+1);
    outTree_.matchedSimHits_modType->resize(idx+1);
    outTree_.matchedSimHits_crystal->resize(idx+1);
    outTree_.matchedSimHits_ieta->resize(idx+1);
    outTree_.matchedSimHits_iphi->resize(idx+1);
    outTree_.matchedSimHits_entry_local_x->resize(idx+1);
    outTree_.matchedSimHits_entry_local_y->resize(idx+1);
    outTree_.matchedSimHits_entry_local_z->resize(idx+1);
    outTree_.matchedSimHits_entry_global_R->resize(idx+1);
    outTree_.matchedSimHits_exit_local_x->resize(idx+1);
    outTree_.matchedSimHits_exit_local_y->resize(idx+1);
    outTree_.matchedSimHits_exit_local_z->resize(idx+1);
    outTree_.matchedSimHits_exit_global_R->resize(idx+1);
    outTree_.matchedSimHits_track_Deta->resize(idx+1);
    outTree_.matchedSimHits_track_Dphi->resize(idx+1);
    outTree_.matchedSimHits_track_DR->resize(idx+1);
    outTree_.matchedSimHits_track_Dz->resize(idx+1);
    outTree_.matchedSimHits_track_RDphi->resize(idx+1);
    outTree_.matchedSimHits_track_dist->resize(idx+1);
    
    outTree_.matchedRecHits_idx->resize(idx+1);
    outTree_.matchedRecHits_energy->resize(idx+1);
    outTree_.matchedRecHits_energyCorr->resize(idx+1);
    outTree_.matchedRecHits_time->resize(idx+1);
    outTree_.matchedRecHits_rr->resize(idx+1);
    outTree_.matchedRecHits_module->resize(idx+1);
    outTree_.matchedRecHits_modType->resize(idx+1);
    outTree_.matchedRecHits_crystal->resize(idx+1);
    outTree_.matchedRecHits_ieta->resize(idx+1);
    outTree_.matchedRecHits_iphi->resize(idx+1);
    outTree_.matchedRecHits_local_x->resize(idx+1);
    outTree_.matchedRecHits_local_y->resize(idx+1);
    outTree_.matchedRecHits_local_z->resize(idx+1);
    outTree_.matchedRecHits_global_R->resize(idx+1);
    outTree_.matchedRecHits_track_Deta->resize(idx+1);
    outTree_.matchedRecHits_track_Dphi->resize(idx+1);
    outTree_.matchedRecHits_track_DR->resize(idx+1);
    outTree_.matchedRecHits_track_Dz->resize(idx+1);
    outTree_.matchedRecHits_track_RDphi->resize(idx+1);
    outTree_.matchedRecHits_track_dist->resize(idx+1);
    outTree_.matchedRecHits_sietaieta->resize(idx+1);
    outTree_.matchedRecHits_siphiiphi->resize(idx+1);
    
    
    // track extrapolation
    const GlobalPoint vtx_outer(track.outerPosition().x(),track.outerPosition().y(),track.outerPosition().z());
    GlobalVector vec_outer(track.outerMomentum().x(),track.outerMomentum().y(),track.outerMomentum().z());
    CurvilinearTrajectoryError err_outer(track.outerStateCovariance());
    const FreeTrajectoryState fts_outer(GlobalTrajectoryParameters(vtx_outer,vec_outer,track.charge(),theField.product()),err_outer);
      
    const Surface::RotationType dummyRot;
    
    std::vector<float> cyl_R;
    if( (crysLayout_ == BTLDetId::CrysLayout::tile) )
    {
      cyl_R.push_back(117.450);
      cyl_R.push_back(117.456);
      cyl_R.push_back(117.473);
      cyl_R.push_back(117.501);
      cyl_R.push_back(117.540);
      cyl_R.push_back(117.590);
      cyl_R.push_back(117.653);
      cyl_R.push_back(117.723);
      cyl_R.push_back(117.810);
    }
    if( crysLayout_ == BTLDetId::CrysLayout::bar )
    {
      cyl_R.push_back(117.450);
      cyl_R.push_back(117.541);
      cyl_R.push_back(117.812);
    }
    if( crysLayout_ == BTLDetId::CrysLayout::barzflat )
    {
      cyl_R.push_back(117.450);
      cyl_R.push_back(117.450); cyl_R.push_back(117.451); cyl_R.push_back(117.453); cyl_R.push_back(117.456); cyl_R.push_back(117.459); cyl_R.push_back(117.462); cyl_R.push_back(117.467); cyl_R.push_back(117.473);
      cyl_R.push_back(117.479); cyl_R.push_back(117.485); cyl_R.push_back(117.493); cyl_R.push_back(117.500); cyl_R.push_back(117.510); cyl_R.push_back(117.519); cyl_R.push_back(117.530); cyl_R.push_back(117.540);
      cyl_R.push_back(117.552); cyl_R.push_back(117.565); cyl_R.push_back(117.578); cyl_R.push_back(117.591); cyl_R.push_back(117.606); cyl_R.push_back(117.621); cyl_R.push_back(117.637); cyl_R.push_back(117.654);
      cyl_R.push_back(117.671); cyl_R.push_back(117.689); cyl_R.push_back(117.708); cyl_R.push_back(117.727); cyl_R.push_back(117.747); cyl_R.push_back(117.768); cyl_R.push_back(117.789); cyl_R.push_back(117.812);
    }
    
    int cylIt = 0;
    std::map<int,GlobalPoint> gp_ext;
    for(auto val : cyl_R)
    {
      Cylinder::ConstCylinderPointer theTargetCylinder = Cylinder::build(val,Surface::PositionType(0.,0.,0.),dummyRot);
      
      std::pair<TrajectoryStateOnSurface,double> aTsosPath_outer(propagator->propagateWithPath(fts_outer,*theTargetCylinder));
      gp_ext[cylIt] = GlobalPoint(0.,0.,0.);
      if( aTsosPath_outer.first.isValid() )
      {
        GlobalPoint temp(aTsosPath_outer.first.globalPosition().x(),aTsosPath_outer.first.globalPosition().y(),aTsosPath_outer.first.globalPosition().z());
        gp_ext[cylIt] = GlobalPoint(temp);
        if( verbosity_ ) std::cout << "*** track extrapolation: " << PrintPosition(temp) << std::endl;
        
        if( cylIt == 0 )
        {
          outTree_.track_eta_atBTL -> push_back(gp_ext[cylIt].eta());
          outTree_.track_phi_atBTL -> push_back(gp_ext[cylIt].phi());
        }
      }
      else
      {
        if( cylIt == 0 )
        {
          outTree_.track_eta_atBTL -> push_back(-999.);
          outTree_.track_phi_atBTL -> push_back(-999.);
        }
      }
      
      ++cylIt;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    
    /*
    //---get compatible layers
    const vector<const DetLayer*>& layers = layerGeo -> allBTLLayers();
    
    GlobalTrajectoryParameters gtp(vtx_inner,vec_inner,track.charge(),theField.product());
    SteppingHelixPropagator prop(theField.product(),alongMomentum);
    
    float theMaxChi2 = 25.;
    float theNSigma = 3.;
    theEstimator = new Chi2MeasurementEstimator(theMaxChi2,theNSigma);
    
    if( verbosity_ ) std::cout << "*** get compatible layers" << std::endl;
    int it = 0;
    for(auto ilay = layers.begin(); ilay!=layers.end(); ++ilay)
    {
      const MTDTrayBarrelLayer* layer = (const MTDTrayBarrelLayer*)(*ilay);
      const BoundCylinder& cyl = layer -> specificSurface();
      TrajectoryStateOnSurface tsos(gtp,cyl);
      
      std::pair<bool,TrajectoryStateOnSurface> comp = layer -> compatible(tsos,prop,*theEstimator);
      
      if( verbosity_ ) std::cout << ">>> layer " << it << "   cylinder R: " << cyl.radius() << " cm   cylinder half length: " << cyl.bounds().length()/2. << " is compatible: " << comp.first << std::endl;
      
      std::vector<DetLayer::DetWithState> compDets = layer->compatibleDets(tsos,prop,*theEstimator);
      
      if( verbosity_ )
      {
        std::cout << ">>>>>> number of compatibleDets: " << compDets.size() << std::endl;
        
        int it2 = 0;
        for(auto compIt : compDets)
        {
          std::cout << ">>>>>>>>> compatibleDet " << it2 << ":   final state pos: " << PrintPosition(compIt.second.globalPosition())          << std::endl;
          std::cout << ">>>>>>>>> compatibleDet " << it2 << ":           det pos: " << PrintPosition(compIt.first->position())                << std::endl;
          std::cout << ">>>>>>>>> compatibleDet " << it2 << ":          distance: " << (tsos.globalPosition()-compIt.first->position()).mag() << std::endl;
        }
      }
      
      ++it;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;    
    */
    
    
    
    //---get associated simHits
    if( verbosity_ ) std::cout << "*** simHits - n tot: " << simHits.size() << std::endl;
    int simHitIt = 0;
    for(auto simHit : simHits)
    {
      BTLDetId id = simHit.detUnitId();
      DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+14*(id.modType()-1),0,1);
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
      
      GlobalPoint gp_track = gp_ext[abs(id.row(topo.nrows())-int(topo.nrows()/2))];
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp_mid.z()      : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp_mid-gp_track).mag()      : -999.;
      
      if( DR < 0.05 && DR > 0. )
      {
        if( verbosity_ )  std::cout << ">>> topology:   nRows: " << topo.nrows() << "   nColumns: " << topo.ncolumns() << "   pitchx: " << topo.pitch().first << "   pitchy: " << topo.pitch().second << std::endl;
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
      
      outTree_.matchedSimHits_n->at(idx) += 1;
      
      outTree_.matchedSimHits_idx->at(idx).push_back(idx);
      outTree_.matchedSimHits_energy->at(idx).push_back(energy);
      outTree_.matchedSimHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.matchedSimHits_time->at(idx).push_back(time);
      outTree_.matchedSimHits_rr->at(idx).push_back(RR);
      outTree_.matchedSimHits_module->at(idx).push_back(module);
      outTree_.matchedSimHits_modType->at(idx).push_back(modType);
      outTree_.matchedSimHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedSimHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedSimHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedSimHits_entry_local_x->at(idx).push_back(lp_entry.x());
      outTree_.matchedSimHits_entry_local_y->at(idx).push_back(lp_entry.y());
      outTree_.matchedSimHits_entry_local_z->at(idx).push_back(lp_entry.z());
      outTree_.matchedSimHits_entry_global_R->at(idx).push_back(sqrt(gp_entry.perp2()));
      outTree_.matchedSimHits_exit_local_x->at(idx).push_back(lp_exit.x());
      outTree_.matchedSimHits_exit_local_y->at(idx).push_back(lp_exit.y());
      outTree_.matchedSimHits_exit_local_z->at(idx).push_back(lp_exit.z());
      outTree_.matchedSimHits_exit_global_R->at(idx).push_back(sqrt(gp_exit.perp2()));
      outTree_.matchedSimHits_track_Deta->at(idx).push_back(fabs(Deta));
      outTree_.matchedSimHits_track_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.matchedSimHits_track_DR->at(idx).push_back(DR);
      outTree_.matchedSimHits_track_Dz->at(idx).push_back(Dz);
      outTree_.matchedSimHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedSimHits_track_dist->at(idx).push_back(dist);
      
      ++simHitIt;
    }
    if( verbosity_ ) std::cout << "---" << std::endl;
    
    
    //---find associated recHits
    float sieie=0, sipip=0;
    float ss_hit_count=0;
    int recHitIt = 0;
    if( verbosity_ ) std::cout << "*** recHits - tot: " << recHits.size() << std::endl;
    for(auto recHit : recHits)
    {
      BTLDetId id = recHit.id();
      DetId geoId = BTLDetId(id.mtdSide(),id.mtdRR(),id.module()+14*(id.modType()-1),0,1);
      const auto& det = mtdGeometry_ -> idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
      
      double energy = recHit.energy();
      double time   = recHit.time();
      
      MeasurementPoint mp(id.row(topo.nrows()),id.column(topo.nrows()));
      LocalPoint lp = topo.localPosition(mp);
      GlobalPoint gp = det->toGlobal(lp);
      
      int RR = id.mtdRR();
      int module = id.module();
      int modType = id.modType();
      int crystal = id.crystal();
      int ieta = id.ieta(crysLayout_);
      int iphi = id.iphi(crysLayout_);
      float eta = gp.eta();
      float phi = gp.phi();
      
      GlobalPoint gp_track = gp_ext[abs(id.row(topo.nrows())-int(topo.nrows()/2))];
      
      float Deta  = gp_track.mag() > 0. ? eta-gp_track.eta()           : -999.;
      float Dphi  = gp_track.mag() > 0. ? deltaPhi(phi,gp_track.phi()) : -999.;
      float DR    = gp_track.mag() > 0. ? sqrt(Deta*Deta+Dphi*Dphi)    : -999.;
      float Dz    = gp_track.mag() > 0. ? gp_track.z()-gp.z()          : -999.;
      float RDphi = gp_track.mag() > 0. ? sqrt(gp_track.perp2())*Dphi  : -999.;
      float dist  = gp_track.mag() > 0. ? (gp-gp_track).mag()          : -999.;
      
      // if( DR < 0.2 && DR > 0. )
      {
        if( verbosity_ )  std::cout << ">>> topology:   nRows: " << topo.nrows() << "   nColumns: " << topo.ncolumns() << "   pitchx: " << topo.pitch().first << "   pitchy: " << topo.pitch().second << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":   energy: " << energy << " MeV   time: " << time << " ns"
                                   << "   RR: " << RR << "   module: " << module << "   modType: " << modType << "   crystal: " << crystal
                                   << "   ieta: " << ieta << "   iphi: " << iphi << "   row: " << recHit.row() << "- " << id.row(topo.nrows()) << "   column: " << recHit.column() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ":  detPosition:  global: " << PrintPosition(det->position()) << "   DR: " << DR << "   dist: " << (gp_track-det->position()).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": hit position:   local: " << PrintPosition(lp) << "   global: " << PrintPosition(gp) << "DR: " << DR << "   dist: " << (gp_track-gp).mag() << std::endl;
        if( verbosity_ ) std::cout << ">>> " << recHitIt << ": extrapolated track " << abs(recHit.row()-int(topo.nrows()/2)) << " : " << PrintPosition(gp_track) << std::endl;
      }
      
      if( DR > track_hit_DRMax_ || dist > track_hit_distMax_ ) continue;
      if( gp_track.mag() <= 0. ) continue;
      
      outTree_.matchedRecHits_n->at(idx) += 1;
      
      outTree_.matchedRecHits_idx->at(idx).push_back(idx);
      outTree_.matchedRecHits_energy->at(idx).push_back(energy);
      outTree_.matchedRecHits_energyCorr->at(idx).push_back(energy*fabs(sin(track.theta())));
      outTree_.matchedRecHits_time->at(idx).push_back(time);
      outTree_.matchedRecHits_rr->at(idx).push_back(RR);
      outTree_.matchedRecHits_module->at(idx).push_back(module);
      outTree_.matchedRecHits_modType->at(idx).push_back(modType);
      outTree_.matchedRecHits_crystal->at(idx).push_back(crystal);
      outTree_.matchedRecHits_ieta->at(idx).push_back(ieta);
      outTree_.matchedRecHits_iphi->at(idx).push_back(iphi);
      outTree_.matchedRecHits_local_x->at(idx).push_back(lp.x());
      outTree_.matchedRecHits_local_y->at(idx).push_back(lp.y());
      outTree_.matchedRecHits_local_z->at(idx).push_back(lp.z());
      outTree_.matchedRecHits_global_R->at(idx).push_back(sqrt(gp.perp2()));
      outTree_.matchedRecHits_track_Deta->at(idx).push_back(fabs(Deta));
      outTree_.matchedRecHits_track_Dphi->at(idx).push_back(fabs(Dphi));
      outTree_.matchedRecHits_track_Dz->at(idx).push_back(Dz);
      outTree_.matchedRecHits_track_RDphi->at(idx).push_back(RDphi);
      outTree_.matchedRecHits_track_DR->at(idx).push_back(DR);
      outTree_.matchedRecHits_track_dist->at(idx).push_back(dist);
      
      if(recHit.energy() > 0.5)
      {
        sieie += energy*pow(eta-gp_track.eta(),2);
        sipip += energy*pow(phi-gp_track.phi(),2);
        ss_hit_count += energy;
      }
      
      ++recHitIt;
    }
    if( verbosity_ ) std::cout << "---\n\n\n" << std::endl;
    
    outTree_.matchedRecHits_sietaieta->at(idx).push_back( ss_hit_count>0 ? sqrt(sieie)/ss_hit_count : -999. );
    outTree_.matchedRecHits_siphiiphi->at(idx).push_back( ss_hit_count>0 ? sqrt(sipip)/ss_hit_count : -999. );
    
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
