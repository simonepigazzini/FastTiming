#ifndef _MTD_4D_VERTEXING_ANALIZER_
#define _MTD_4D_VERTEXING_ANALIZER_

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "PrecisionTiming/FTLAnalysis/interface/MTD4DTree.h"
#include "PrecisionTiming/FTLAnalysis/interface/MTDTOFPIDTree.h"

using namespace std;                             

class MTD4DVertexingAnalyzer : public edm::EDAnalyzer
{
public:                             

    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;
    
    explicit MTD4DVertexingAnalyzer(const edm::ParameterSet& pSet);
    ~MTD4DVertexingAnalyzer() {};

    //---methods
    virtual void beginJob() override {};
    virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
    virtual void endJob() override {};
    
private:
    //---inputs
    //---tracks
    edm::Handle<reco::GenParticleCollection> genParticlesHandle_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

    edm::Handle<TrackingParticleCollection> simParticlesHandle_;
    edm::EDGetTokenT<TrackingParticleCollection> simParticlesToken_;
    edm::EDGetTokenT<reco::RecoToSimCollection> trkRecoToSimMapToken_;
    edm::Handle<reco::RecoToSimCollection> trkRecoToSimMapHandle_;    
    edm::EDGetTokenT<reco::SimToRecoCollection> trkSimToRecoMapToken_;
    edm::Handle<reco::SimToRecoCollection> trkSimToRecoMapHandle_;    
    
    edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_;
    edm::Handle<edm::View<reco::Track> > tracksHandle_;
    edm::EDGetTokenT<edm::View<reco::Track> > extTracksToken_;
    edm::Handle<edm::View<reco::Track> > extTracksHandle_;

    //---TOFPID
    edm::EDGetTokenT<edm::ValueMap<float> > t0PIDToken_;
    edm::Handle<edm::ValueMap<float> >      t0PIDHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > sigmat0PIDToken_;
    edm::Handle<edm::ValueMap<float> >      sigmat0PIDHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > probPiPIDToken_;
    edm::Handle<edm::ValueMap<float> >      probPiPIDHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > probPPIDToken_;
    edm::Handle<edm::ValueMap<float> >      probPPIDHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > probKPIDToken_;
    edm::Handle<edm::ValueMap<float> >      probKPIDHandle_;
    
    //---vertices
    edm::EDGetTokenT<genXYZ>                genXYZToken_;
    edm::Handle<genXYZ>                     genXYZHandle_;
    edm::EDGetTokenT<float>                 genT0Token_;
    edm::Handle<float>                      genT0Handle_;
    edm::EDGetTokenT<vector<reco::Vertex> > vtx4DToken_;
    edm::Handle<vector<reco::Vertex> >      vtx4DHandle_;    
    edm::EDGetTokenT<vector<reco::Vertex> > vtx4DNoPIDToken_;
    edm::Handle<vector<reco::Vertex> >      vtx4DNoPIDHandle_;    

    //---options

    //---outputs
    MTD4DTree vtxsTree_;
    MTDTOFPIDTree trksTree_;
    edm::Service<TFileService> fs_;  
  
};

MTD4DVertexingAnalyzer::MTD4DVertexingAnalyzer(const edm::ParameterSet& pSet):
    genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
    trkRecoToSimMapToken_(consumes<reco::RecoToSimCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackAndTrackingParticlesAssociatorMapTag"))),
    trkSimToRecoMapToken_(consumes<reco::SimToRecoCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackAndTrackingParticlesAssociatorMapTag"))),
    tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("generalTracksTag"))),
    extTracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("extendedTracksTag"))),
    t0PIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("t0TOFPIDTag"))),
    sigmat0PIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("sigmat0TOFPIDTag"))),
    probPiPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probPiTOFPIDTag"))),
    probPPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probPTOFPIDTag"))),
    probKPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probKTOFPIDTag"))),    
    genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
    genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
    vtx4DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DTag"))),
    vtx4DNoPIDToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DNoPIDTag")))        
{
    vtxsTree_ = MTD4DTree(pSet.getUntrackedParameter<string>("vtxsTreeName").c_str(), "4D vertexing studies");
    trksTree_ = MTDTOFPIDTree(pSet.getUntrackedParameter<string>("trksTreeName").c_str(), "4D TOFPID studies");
}

void MTD4DVertexingAnalyzer::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    vtxsTree_.Reset();
    trksTree_.Reset();
  
    //--- load gen particles
    event.getByToken(genParticlesToken_, genParticlesHandle_);
    auto genParticles = *genParticlesHandle_.product();

    //--- load sim particles (aka trackingParticles)
    // event.getByToken(trkRecoToSimMapToken_, trkRecoToSimMapHandle_);
    // auto trkRecoToSimMap = *trkRecoToSimMapHandle_.product();
    // event.getByToken(trkSimToRecoMapToken_, trkSimToRecoMapHandle_);
    // auto trkSimToRecoMap = *trkSimToRecoMapHandle_.product();  
  
    //---load general tracks
    event.getByToken(tracksToken_,tracksHandle_);
    auto tracks = *tracksHandle_.product();

    //---load extended tracks
    event.getByToken(extTracksToken_, extTracksHandle_);
    auto extTracks = *extTracksHandle_.product();

    //---load TOFPID information
    event.getByToken(t0PIDToken_, t0PIDHandle_);
    auto t0PID = *t0PIDHandle_.product();    
    event.getByToken(sigmat0PIDToken_, sigmat0PIDHandle_);
    auto sigmat0PID = *sigmat0PIDHandle_.product();    
    event.getByToken(probPiPIDToken_, probPiPIDHandle_);
    auto probPiPID = *probPiPIDHandle_.product();    
    event.getByToken(probPPIDToken_, probPPIDHandle_);
    auto probPPID = *probPPIDHandle_.product();    
    event.getByToken(probKPIDToken_, probKPIDHandle_);
    auto probKPID = *probKPIDHandle_.product();    
    
    //---load sim and reco vertices
    // SIM
    event.getByToken(genXYZToken_, genXYZHandle_);
    event.getByToken(genT0Token_, genT0Handle_);
    auto xyz = genXYZHandle_.product();
    auto t = *genT0Handle_.product();
    auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
    auto genPV = SimVertex(v, t).position();
    // Full 4D
    event.getByToken(vtx4DToken_, vtx4DHandle_);
    auto vtxs4D = *vtx4DHandle_.product();
    // 4D no PID
    event.getByToken(vtx4DNoPIDToken_, vtx4DNoPIDHandle_);
    auto vtxs4DNoPID = *vtx4DNoPIDHandle_.product();

    float min_dz=1e6, min_dzt=1e6;
    int best4D_dz_idx=-1, best4D_dzt_idx=-1;
    for(unsigned int iVtx=0; iVtx<vtxs4D.size(); ++iVtx)
    {
        auto& vtx4D = vtxs4D[iVtx];
        if(vtx4D.isValid() && !vtx4D.isFake())
        {
            auto dz = std::abs(vtx4D.z()-genPV.z());
            auto dzt = sqrt(dz*dz + pow((vtx4D.t()-genPV.t())*30, 2));
            if(dz < min_dz)
            {
                min_dz = dz;
                best4D_dz_idx = iVtx;
            }
            if(dzt < min_dzt)
            {
                min_dzt = dzt;
                best4D_dzt_idx = iVtx;
            }
        }
    }
    auto& best_dz_4D_vtx = vtxs4D[best4D_dz_idx];
    auto& best_dzt_4D_vtx = vtxs4D[best4D_dzt_idx];    

    //---Fill output tree
    // 0th 4D vtx    
    vtxsTree_.vtx4D_0_valid = vtxs4D[0].isValid() && !vtxs4D[0].isFake();
    vtxsTree_.vtx4D_0_dz = vtxs4D[0].z()-genPV.z();
    vtxsTree_.vtx4D_0_dt = vtxs4D[0].t()-genPV.t();
    vtxsTree_.vtx4D_0_chi2 = vtxs4D[0].chi2()/vtxs4D[0].ndof();
    vtxsTree_.vtx4D_0_ntrks = vtxs4D[0].nTracks();
    // 0th 4D w/o PID
    vtxsTree_.vtx4DNoPID_Nnofake = 0;
    for(auto& vtx : vtxs4DNoPID)
    {
        if(vtx.isValid() && !vtx.isFake())
            vtxsTree_.vtx4DNoPID_Nnofake++;
    }
    vtxsTree_.vtx4DNoPID_N = vtxs4DNoPID.size();
    vtxsTree_.vtx4DNoPID_0_valid = vtxs4DNoPID[0].isValid() && !vtxs4DNoPID[0].isFake();
    vtxsTree_.vtx4DNoPID_0_dz = vtxs4DNoPID[0].z()-genPV.z();
    vtxsTree_.vtx4DNoPID_0_dt = vtxs4DNoPID[0].t()-genPV.t();
    vtxsTree_.vtx4DNoPID_0_chi2 = vtxs4DNoPID[0].chi2()/vtxs4DNoPID[0].ndof();
    vtxsTree_.vtx4DNoPID_0_ntrks = vtxs4DNoPID[0].nTracks();
    // dz matching 4D std
    vtxsTree_.vtx4D_best_idx = best4D_dz_idx;
    vtxsTree_.vtx4D_best_dz = best_dz_4D_vtx.z()-genPV.z();
    vtxsTree_.vtx4D_best_dt = best_dz_4D_vtx.t()-genPV.t();
    vtxsTree_.vtx4D_best_chi2 = best_dz_4D_vtx.chi2()/best_dz_4D_vtx.ndof();
    vtxsTree_.vtx4D_best_ntrks = best_dz_4D_vtx.nTracks();
    // dzt matching
    vtxsTree_.vtx4D_best_dzt_idx = best4D_dzt_idx;
    vtxsTree_.vtx4D_best_dzt_dz = best_dzt_4D_vtx.z()-genPV.z();
    vtxsTree_.vtx4D_best_dzt_dt = best_dzt_4D_vtx.t()-genPV.t();

    vtxsTree_.GetTTreePtr()->Fill();

    //---tracks loop
    //  - Check PID hypothesis, and vertex association comparing 4D and 4DnoPID
    for(unsigned int itrack=0; itrack<extTracks.size(); ++itrack)
    {
        auto& track = extTracks[itrack];
        reco::TrackBaseRef track_ref(tracksHandle_, itrack);
        //---get TOFPIDProducer recomputed t0,sigmat0 and particle probs
        float t_t0=-1, t_sigmat0=-1, t_probPi=-1, t_probP=-1, t_probK=-1;
        if(t0PID.contains(track_ref.id()))
        {
            t_t0 = t0PID[track_ref];
            t_sigmat0 = sigmat0PID[track_ref];
            t_probPi = probPiPID[track_ref];
            t_probP = probPPID[track_ref];
            t_probK = probKPID[track_ref];
        }
        
        // match with gen particles
        float DRMin = 1e9;
        int genPdgId = 0;
        float genEta = -999.;
        float genPhi = -999.;
        float genPt = -999.;
        for(auto& genPart : genParticles)
	{
            if( genPart.status() != 1 ) continue;
            if( genPart.charge() == 0 ) continue;

            float DR   = deltaR(track.eta(), track.phi(), genPart.eta(), genPart.phi());
      
            if( DR < DRMin )
	    {
                DRMin = DR;
        
                genPdgId = genPart.pdgId();
                genEta   = genPart.eta();
                genPhi   = genPart.phi();
                genPt    = genPart.pt();
	    }
	}
        
        trksTree_.trk_idx -> push_back(itrack);
        trksTree_.trk_pt -> push_back(track.pt());
        trksTree_.trk_eta -> push_back(track.eta());
        trksTree_.trk_phi -> push_back(track.phi());
        trksTree_.trk_x -> push_back(track.vx());
        trksTree_.trk_y -> push_back(track.vy());
        trksTree_.trk_z -> push_back(track.vz());
        trksTree_.trk_sigmaz -> push_back(track.dzError());
        trksTree_.trk_t -> push_back(track.t0());
        trksTree_.trk_sigmat -> push_back(track.t0Error());
        trksTree_.trk_PID_t -> push_back(t_t0);
        trksTree_.trk_PID_sigmat -> push_back(t_sigmat0);
        trksTree_.trk_PID_probPi -> push_back(t_probPi);
        trksTree_.trk_PID_probP -> push_back(t_probP);
        trksTree_.trk_PID_probK -> push_back(t_probK);        
        trksTree_.trk_energy -> push_back(sqrt(track.momentum().mag2()));
        trksTree_.trk_normalizedChi2 -> push_back(track.normalizedChi2());
        trksTree_.trk_numberOfValidHits -> push_back(track.numberOfValidHits());
        trksTree_.trk_numberOfLostHits -> push_back(track.numberOfLostHits());
        trksTree_.trk_isHighPurity -> push_back(track.quality(reco::TrackBase::TrackQuality::highPurity));
        trksTree_.trk_hasMTD -> push_back(track.isTimeOk());
        trksTree_.trk_genPdgId -> push_back(genPdgId);
        trksTree_.trk_genPt -> push_back(genPt);
        trksTree_.trk_genEta -> push_back(genEta);
        trksTree_.trk_genPhi -> push_back(genPhi);
        trksTree_.trk_genDR -> push_back(DRMin);
        trksTree_.trk_genVtx_x -> push_back(genPV.x());
        trksTree_.trk_genVtx_y -> push_back(genPV.y());        
        trksTree_.trk_genVtx_z -> push_back(genPV.z());        
        trksTree_.trk_genVtx_t -> push_back(genPV.t());        
    }

    trksTree_.GetTTreePtr()->Fill();
}


DEFINE_FWK_MODULE(MTD4DVertexingAnalyzer);

#endif
