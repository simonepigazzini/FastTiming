#ifndef _MTD_NEUTRALS_ANALIZER_
#define _MTD_NEUTRALS_ANALIZER_

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
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDRingForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetRing.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "PrecisionTiming/FTLAnalysis/interface/MTDNeutralsTree.h"

using namespace std;                             

class MTDNeutralsAnalyzer : public edm::EDAnalyzer
{
public:                             

    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;
    
    explicit MTDNeutralsAnalyzer(const edm::ParameterSet& pSet);
    ~MTDNeutralsAnalyzer() {};

    //---methods
    virtual void beginJob() override {};
    virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
    virtual void endJob() override {};
    
private:
    //---inputs
    const MTDGeometry* mtdGeometry_;

    //---gen particles
    edm::Handle<reco::GenParticleCollection> genParticlesHandle_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

    //---sim MTD hits
    edm::Handle<std::vector<PSimHit> > simHitsBTLHandle_;
    edm::EDGetTokenT<std::vector<PSimHit> > simHitsBTLToken_;    
    edm::Handle<std::vector<PSimHit> > simHitsETLHandle_;
    edm::EDGetTokenT<std::vector<PSimHit> > simHitsETLToken_;    
    
    //---BTL clusters
    edm::Handle<FTLRecHitCollection> recHitsBTLHandle_;
    edm::EDGetTokenT<FTLRecHitCollection> recHitsBTLToken_;    
    edm::Handle<FTLClusterCollection> clustersBTLHandle_;
    edm::EDGetTokenT<FTLClusterCollection> clustersBTLToken_;    
    //---ETL clusters
    edm::Handle<FTLRecHitCollection> recHitsETLHandle_;
    edm::EDGetTokenT<FTLRecHitCollection> recHitsETLToken_;    
    edm::Handle<FTLClusterCollection> clustersETLHandle_;
    edm::EDGetTokenT<FTLClusterCollection> clustersETLToken_;    
    
    //---PF
    edm::Handle<reco::PFCandidateCollection> pfCandidatesHandle_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;

    //---vertices
    edm::EDGetTokenT<genXYZ>                  genXYZToken_;
    edm::Handle<genXYZ>                       genXYZHandle_;
    edm::EDGetTokenT<float>                   genT0Token_;
    edm::Handle<float>                        genT0Handle_;
    edm::EDGetTokenT<vector<reco::Vertex> >   vtx3DToken_;
    edm::Handle<vector<reco::Vertex> >        vtx3DHandle_;
    edm::EDGetTokenT<vector<reco::Vertex> >   vtx4DToken_;
    edm::Handle<vector<reco::Vertex> >        vtx4DHandle_;    
    
    //---options
    BTLDetId::CrysLayout crysLayout_;
    
    //---outputs
    MTDNeutralsTree outTree_;
    edm::Service<TFileService> fs_;  
  
};

MTDNeutralsAnalyzer::MTDNeutralsAnalyzer(const edm::ParameterSet& pSet):
    genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
    simHitsBTLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsBTLTag"))), 
    simHitsETLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsETLTag"))),   
    clustersBTLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersBTLTag"))),    
    clustersETLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersETLTag"))),
    pfCandidatesToken_(consumes<reco::PFCandidateCollection>(pSet.getUntrackedParameter<edm::InputTag>("pfCandidatesTag"))),    
    genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
    genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
    vtx3DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx3DTag"))),
    vtx4DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DTag"))),
    crysLayout_((BTLDetId::CrysLayout)(pSet.getUntrackedParameter<int>("crysLayout")))
{
    outTree_ = MTDNeutralsTree(pSet.getUntrackedParameter<string>("outTreeName").c_str(), "4D TOFPID studies");
}

void MTDNeutralsAnalyzer::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    //---load gen particles
    event.getByToken(genParticlesToken_, genParticlesHandle_);
    auto genParticles = *genParticlesHandle_.product();

    //---load sim hits
    event.getByToken(simHitsBTLToken_, simHitsBTLHandle_);
    auto simHitsBTL = *simHitsBTLHandle_.product();
    event.getByToken(simHitsETLToken_, simHitsETLHandle_);
    auto simHitsETL = *simHitsETLHandle_.product();

    
    //---get the MTD geometry
    edm::ESHandle<MTDGeometry> geoHandle;
    setup.get<MTDDigiGeometryRecord>().get(geoHandle);
    mtdGeometry_ = geoHandle.product();
  
    edm::ESHandle<MTDDetLayerGeometry> layerGeo;
    setup.get<MTDRecoGeometryRecord>().get(layerGeo);
    
    //---MTD clusters
    event.getByToken(clustersBTLToken_, clustersBTLHandle_);
    auto clustersBTL = *clustersBTLHandle_.product();
    event.getByToken(clustersETLToken_, clustersETLHandle_);
    auto clustersETL = *clustersETLHandle_.product();
    
    //---load PFCandidates
    event.getByToken(pfCandidatesToken_, pfCandidatesHandle_);
    auto pfCandidates = *pfCandidatesHandle_.product();    
    
    //---load gen, sim and reco vertices
    // GEN
    event.getByToken(genXYZToken_, genXYZHandle_);
    event.getByToken(genT0Token_, genT0Handle_);
    auto xyz = genXYZHandle_.product();
    auto t = *genT0Handle_.product();
    auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
    auto genPV = SimVertex(v, t).position();
    // 3D
    event.getByToken(vtx3DToken_, vtx3DHandle_);
    auto vtxs3D = *vtx3DHandle_.product();
    // Full 4D
    event.getByToken(vtx4DToken_, vtx4DHandle_);
    auto vtxs4D = *vtx4DHandle_.product();

    for(auto& cand : pfCandidates)
    {
        outTree_.Reset();
        //---fill global info
        outTree_.event = event.id().event();
        outTree_.lumi = event.id().luminosityBlock();
        outTree_.run = event.id().run();    

        if(cand.pt()>0.5)
        {
            outTree_.genpv_t = genPV.t();
            outTree_.genpv_z = genPV.z();            
            outTree_.pt = cand.pt();
            outTree_.eta = cand.eta();
            outTree_.phi = cand.phi();
            outTree_.ecalEnergy = cand.ecalEnergy();
            outTree_.hcalEnergy = cand.hcalEnergy();
            outTree_.particleId = cand.particleId();

            float dr_min=1e6;
            reco::GenParticle* gen_match = NULL;
            for(auto& gen_part : genParticles)
            {
                if(gen_part.pdgId() == 22 || gen_part.pdgId() == 111)
                {
                    float dr = deltaR(cand.eta(), cand.phi(), gen_part.eta(), gen_part.phi());
                    if(dr<dr_min)
                    {
                        dr_min = dr;
                        gen_match = &gen_part;
                    }
                }
            }
            if(gen_match)
            {
                outTree_.gen_pt = gen_match->pt();
                outTree_.gen_eta = gen_match->eta();                
                outTree_.gen_phi = gen_match->phi();            
                outTree_.gen_DR = dr_min;
                outTree_.gen_pdgId = gen_match->pdgId();
            }
            else
            {
                outTree_.gen_pt = -1;
                outTree_.gen_eta = -1;
                outTree_.gen_phi = -1;
                outTree_.gen_DR = -1;
                outTree_.gen_pdgId = -1;
            }
          
            if(std::abs(cand.eta()) < 1.5)
            {
                //---BTL clusters loop
                for(auto& clusIt : clustersBTL)
                {    
                    DetId id = clusIt.detId();
                    const auto& det = mtdGeometry_ -> idToDet(id);
                    const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
                    const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
                    for(auto& cluster : clusIt)
                    {
                        float energy = cluster.energy();
                        float time   = cluster.time();
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
	  
                        if ( mtdId.mtdSubDetector() == MTDDetId::BTL )
                        {
                            BTLDetId btlId(id);
                            RR = btlId.mtdRR();
                            module = btlId.module();
                            modType = btlId.modType();
                        }
	  
                        outTree_.clus_n += 1;    
                        outTree_.clus_det->push_back(1);
                        outTree_.clus_size->push_back(size);
                        outTree_.clus_size_x->push_back(sizeX);
                        outTree_.clus_size_y->push_back(sizeY);
                        outTree_.clus_energy->push_back(energy);
                        outTree_.clus_time->push_back(time);
                        outTree_.clus_rr->push_back(RR);
                        outTree_.clus_module->push_back(module);
                        outTree_.clus_modType->push_back(modType);
                        outTree_.clus_eta->push_back(gp.eta());
                        outTree_.clus_phi->push_back(gp.phi());
                        outTree_.clus_seed_energy->push_back(seed_energy);
                        outTree_.clus_seed_time->push_back(seed_time);
                        outTree_.clus_seed_x->push_back(seed_x);
                        outTree_.clus_seed_y->push_back(seed_y);
                        outTree_.clus_local_x->push_back(lp.x());
                        outTree_.clus_local_y->push_back(lp.y());
                        outTree_.clus_local_z->push_back(lp.z());
                        outTree_.clus_global_R->push_back(sqrt(gp.perp2()));
                        outTree_.clus_global_dist->push_back(sqrt(pow(gp.x()-genPV.x(), 2)+pow(gp.y()-genPV.y(), 2)+pow(gp.z()-genPV.z(), 2)));
                        outTree_.clus_neu_DEta->push_back(gp.eta()-cand.eta());
                        outTree_.clus_neu_DPhi->push_back(deltaPhi(cand.phi(), gp.phi()));  
                        outTree_.clus_neu_DR->push_back(deltaR(cand.eta(), cand.phi(), gp.eta(), gp.phi()));
                    }
                }
                unsigned int min_dr_pos = std::min_element(outTree_.clus_neu_DR->begin(), outTree_.clus_neu_DR->end())-outTree_.clus_neu_DR->begin();
                if(min_dr_pos < outTree_.clus_neu_DR->size())
                {
                    outTree_.minDR = outTree_.clus_neu_DR->at(min_dr_pos);
                    outTree_.minDEta = outTree_.clus_neu_DEta->at(min_dr_pos);                
                    outTree_.minDPhi = outTree_.clus_neu_DPhi->at(min_dr_pos);
                    outTree_.tof = outTree_.clus_global_dist->at(min_dr_pos)/2.99792458e1;
                    outTree_.mtdTime = outTree_.clus_time->at(min_dr_pos);
                    outTree_.mtdEnergy = outTree_.clus_energy->at(min_dr_pos);
                }

                //---sim hits
                if(gen_match)
                {
                    float minDR_simh = 1e6;
                    for(auto& simHit : simHitsBTL)
                    {
                        if(simHit.tof()<0 || simHit.tof()>25 || simHit.energyLoss()*1000<1.5)
                            continue;

                        BTLDetId id = simHit.detUnitId();
                        DetId geoId = id.geographicalId( crysLayout_ );
                        const auto& det = mtdGeometry_ -> idToDet(geoId);
                        const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
                        const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
                                
                        LocalPoint lp_entry(simHit.entryPoint().x()/10., simHit.entryPoint().y()/10., simHit.entryPoint().z()/10.);
                        GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,id.row(topo.nrows()),id.column(topo.nrows())));
                        auto dr = deltaR(gp_entry.eta(), gp_entry.phi(), gen_match->eta(), gen_match->phi());
                        if(dr < minDR_simh)
                        {
                            minDR_simh = dr;
                            outTree_.simh_energy = simHit.energyLoss()*1000.;
                            outTree_.simh_time = simHit.tof();
                            outTree_.simh_tof = sqrt(pow(gp_entry.x()-genPV.x(), 2)+pow(gp_entry.y()-genPV.y(), 2)+pow(gp_entry.z()-genPV.z(), 2))/2.99792458e1;
                            outTree_.simh_DR = dr;
                            if(min_dr_pos < outTree_.clus_neu_DR->size())
                            {
                                auto r_eta = outTree_.clus_eta->at(min_dr_pos);
                                auto r_phi = outTree_.clus_phi->at(min_dr_pos);
                                outTree_.simh_recoh_DR = deltaR(gp_entry.eta(), gp_entry.phi(), r_eta, r_phi);
                                outTree_.simh_recoh_DPhi = deltaPhi(gp_entry.phi(), r_phi);                                
                            }
                        }
                    }
                }
            }
            // else if(std::abs(cand.eta()) < 3.)
            // {
            //     for (auto clusIt : clustersETL)
            //     {    
            //         DetId id = clusIt.detId();
            //         const auto& det = mtdGeometry_ -> idToDet(id);
            //         const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
            //         const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
            //         for ( auto cluster : clusIt)
            //         {
            //             float energy = cluster.energy();
            //             float time   = cluster.time();
            //             int size=cluster.size();
            //             int sizeX=cluster.sizeX();
            //             int sizeY=cluster.sizeY();
            //             float seed_energy=cluster.seed().energy();
            //             float seed_time=cluster.seed().time();
            //             int seed_x=cluster.seed().x();
            //             int seed_y=cluster.seed().y();
	  
            //             MeasurementPoint mp(cluster.x(),cluster.y());
            //             LocalPoint lp = topo.localPosition(mp);
            //             GlobalPoint gp = det->toGlobal(lp);
	  
            //             MTDDetId mtdId(id);
            //             int RR = 0;
            //             int module = 0;
            //             int modType = 0;
	  
            //             if ( mtdId.mtdSubDetector() == MTDDetId::ETL )
            //             {
            //                 ETLDetId btlId(id);
            //                 RR = btlId.mtdRR();
            //                 module = btlId.module();
            //                 modType = btlId.modType();
            //             }

            //             outTree_.clus_n += 1;    
            //             outTree_.clus_det->push_back(1);
            //             outTree_.clus_size->push_back(size);
            //             outTree_.clus_size_x->push_back(sizeX);
            //             outTree_.clus_size_y->push_back(sizeY);
            //             outTree_.clus_energy->push_back(energy);
            //             outTree_.clus_time->push_back(time);
            //             outTree_.clus_rr->push_back(RR);
            //             outTree_.clus_module->push_back(module);
            //             outTree_.clus_modType->push_back(modType);
            //             outTree_.clus_eta->push_back(gp.eta());
            //             outTree_.clus_phi->push_back(gp.phi());
            //             outTree_.clus_seed_energy->push_back(seed_energy);
            //             outTree_.clus_seed_time->push_back(seed_time);
            //             outTree_.clus_seed_x->push_back(seed_x);
            //             outTree_.clus_seed_y->push_back(seed_y);
            //             outTree_.clus_local_x->push_back(lp.x());
            //             outTree_.clus_local_y->push_back(lp.y());
            //             outTree_.clus_local_z->push_back(lp.z());
            //             outTree_.clus_global_R->push_back(sqrt(gp.perp2()));
            //             outTree_.clus_global_dist->push_back(sqrt(pow(gp.x()-genPV.x(), 2)+pow(gp.y()-genPV.y(), 2)+pow(gp.z()-genPV.z(), 2)));                        
            //             outTree_.clus_neu_DEta->push_back(gp.eta()-cand.eta());                        
            //             outTree_.clus_neu_DPhi->push_back(deltaPhi(cand.phi(), gp.phi()));  
            //             outTree_.clus_neu_DR->push_back(deltaR(cand.eta(), cand.phi(), gp.eta(), gp.phi()));                        
            //         }                    
            //     }
            //     unsigned int min_dr_pos = std::min_element(outTree_.clus_neu_DR->begin(), outTree_.clus_neu_DR->end())-outTree_.clus_neu_DR->begin();
            //     if(min_dr_pos < outTree_.clus_neu_DR->size())
            //     {
            //         outTree_.minDR = outTree_.clus_neu_DR->at(min_dr_pos);
            //         outTree_.minDEta = outTree_.clus_neu_DEta->at(min_dr_pos);
            //         outTree_.minDPhi = outTree_.clus_neu_DPhi->at(min_dr_pos);
            //         outTree_.tof = outTree_.clus_global_dist->at(min_dr_pos)/2.99792458e1;                    
            //         outTree_.mtdTime = outTree_.clus_time->at(min_dr_pos);
            //         outTree_.mtdEnergy = outTree_.clus_energy->at(min_dr_pos);
            //     }
            // }
            
            //---Fill trees
            outTree_.GetTTreePtr()->Fill();                    
        }
    }
}


DEFINE_FWK_MODULE(MTDNeutralsAnalyzer);

#endif
