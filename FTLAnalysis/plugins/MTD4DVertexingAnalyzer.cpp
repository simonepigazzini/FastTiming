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

#include "DataFormats/RecoCandidate/interface/RecoChargedRefCandidate.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "PrecisionTiming/FTLAnalysis/interface/TrackPUIDMVA.h"
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
    edm::EDGetTokenT<edm::ValueMap<int> > generalToExtendedTrkMapToken_;
    edm::Handle<edm::ValueMap<int> > generalToExtendedTrkMapHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > extTracksMTDtimeToken_;
    edm::Handle<edm::ValueMap<float> > extTracksMTDtimeHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > pathLengthToken_;
    edm::Handle<edm::ValueMap<float> > pathLengthHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > btlMatchChi2Token_;
    edm::Handle<edm::ValueMap<float> > btlMatchChi2Handle_;
    edm::EDGetTokenT<edm::ValueMap<float> > btlMatchTimeChi2Token_;
    edm::Handle<edm::ValueMap<float> > btlMatchTimeChi2Handle_;
    edm::EDGetTokenT<edm::ValueMap<float> > etlMatchChi2Token_;
    edm::Handle<edm::ValueMap<float> > etlMatchChi2Handle_;
    edm::EDGetTokenT<edm::ValueMap<float> > etlMatchTimeChi2Token_;
    edm::Handle<edm::ValueMap<float> > etlMatchTimeChi2Handle_;
    
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
    edm::EDGetTokenT<genXYZ>                  genXYZToken_;
    edm::Handle<genXYZ>                       genXYZHandle_;
    edm::EDGetTokenT<float>                   genT0Token_;
    edm::Handle<float>                        genT0Handle_;
    edm::EDGetTokenT<vector<TrackingVertex> > simVtxToken_;
    edm::Handle<vector<TrackingVertex> >      simVtxHandle_;    
    edm::EDGetTokenT<vector<reco::Vertex> >   vtx3DToken_;
    edm::Handle<vector<reco::Vertex> >        vtx3DHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> >   vtxScores3DToken_;
    edm::Handle<edm::ValueMap<float> >        vtxScores3DHandle_;
    edm::EDGetTokenT<vector<reco::Vertex> >   vtx4DToken_;
    edm::Handle<vector<reco::Vertex> >        vtx4DHandle_;    
    edm::EDGetTokenT<edm::ValueMap<float> >   vtxScores4DToken_;
    edm::Handle<edm::ValueMap<float> >        vtxScores4DHandle_;
    edm::EDGetTokenT<vector<reco::Vertex> >   vtx4DNoPIDToken_;
    edm::Handle<vector<reco::Vertex> >        vtx4DNoPIDHandle_;    
    
    //---options

    //---outputs
    TrackPUIDMVA mva3D_;
    TrackPUIDMVA mva4D_;    
    MTD4DTree vtxsTree_;
    MTDTOFPIDTree trksTree_;
    edm::Service<TFileService> fs_;  
  
};

MTD4DVertexingAnalyzer::MTD4DVertexingAnalyzer(const edm::ParameterSet& pSet):
    genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
    simParticlesToken_(consumes<TrackingParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackingParticlesTag"))),
    trkRecoToSimMapToken_(consumes<reco::RecoToSimCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackAndTrackingParticlesAssociatorMapTag"))),
    trkSimToRecoMapToken_(consumes<reco::SimToRecoCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackAndTrackingParticlesAssociatorMapTag"))),
    tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("generalTracksTag"))),
    extTracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("extendedTracksTag"))),
    generalToExtendedTrkMapToken_(consumes<edm::ValueMap<int> >(pSet.getUntrackedParameter<edm::InputTag>("generalToExtendedTrkMapTag"))),    
    extTracksMTDtimeToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("extTracksMTDtimeTag"))),
    pathLengthToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("extTracksPathLengthTag"))),
    btlMatchChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("btlMatchChi2Tag"))),
    btlMatchTimeChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("btlMatchTimeChi2Tag"))),
    etlMatchChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("etlMatchChi2Tag"))),
    etlMatchTimeChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("etlMatchTimeChi2Tag"))),
    t0PIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("t0TOFPIDTag"))),
    sigmat0PIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("sigmat0TOFPIDTag"))),
    probPiPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probPiTOFPIDTag"))),
    probPPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probPTOFPIDTag"))),
    probKPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probKTOFPIDTag"))),    
    genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
    genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
    simVtxToken_(consumes<vector<TrackingVertex> >(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),    
    vtx3DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx3DTag"))),
    vtxScores3DToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("vtx3DTag"))),
    vtx4DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DTag"))),
    vtxScores4DToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DTag"))),
    vtx4DNoPIDToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DNoPIDTag"))),
    mva3D_(pSet.getParameter<edm::FileInPath>("trackPUID_3DBDT_weights_file").fullPath(), false),
    mva4D_(pSet.getParameter<edm::FileInPath>("trackPUID_4DBDT_weights_file").fullPath(), true)
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
    
    //--- load sim particles
    event.getByToken(simParticlesToken_, simParticlesHandle_);
    auto simParticles = *simParticlesHandle_.product();
    
    //--- load sim particles (aka trackingParticles)
    event.getByToken(trkRecoToSimMapToken_, trkRecoToSimMapHandle_);
    auto trkRecoToSimMap = *trkRecoToSimMapHandle_.product();
    event.getByToken(trkSimToRecoMapToken_, trkSimToRecoMapHandle_);
    auto trkSimToRecoMap = *trkSimToRecoMapHandle_.product();  
  
    //---load general tracks
    event.getByToken(tracksToken_,tracksHandle_);
    auto tracks = *tracksHandle_.product();

    //---load extended tracks
    event.getByToken(extTracksToken_, extTracksHandle_);
    auto extTracks = *extTracksHandle_.product();
    event.getByToken(generalToExtendedTrkMapToken_, generalToExtendedTrkMapHandle_);
    auto generalToExtendedTrkMap = *generalToExtendedTrkMapHandle_.product();

    //---extended tracks path length
    event.getByToken(pathLengthToken_, pathLengthHandle_);
    auto extPathLenght = *pathLengthHandle_.product();
    event.getByToken(extTracksMTDtimeToken_, extTracksMTDtimeHandle_);
    auto extTracksMTDtime = *extTracksMTDtimeHandle_.product();
    event.getByToken(btlMatchChi2Token_, btlMatchChi2Handle_);
    auto btlMatchChi2 = *btlMatchChi2Handle_.product();
    event.getByToken(btlMatchTimeChi2Token_, btlMatchTimeChi2Handle_);
    auto btlMatchTimeChi2 = *btlMatchTimeChi2Handle_.product();
    event.getByToken(etlMatchChi2Token_, etlMatchChi2Handle_);
    auto etlMatchChi2 = *etlMatchChi2Handle_.product();
    event.getByToken(etlMatchTimeChi2Token_, etlMatchTimeChi2Handle_);
    auto etlMatchTimeChi2 = *etlMatchTimeChi2Handle_.product();
    
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
    
    //---load gen, sim and reco vertices
    // GEN
    event.getByToken(genXYZToken_, genXYZHandle_);
    event.getByToken(genT0Token_, genT0Handle_);
    auto xyz = genXYZHandle_.product();
    auto t = *genT0Handle_.product();
    auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
    auto genPV = SimVertex(v, t).position();
    // SIM
    // event.getByToken(simVtxToken_, simVtxHandle_);
    // auto vtxsSim = *simVtxHandle_.product();    
    // 3D
    event.getByToken(vtx3DToken_, vtx3DHandle_);
    auto vtxs3D = *vtx3DHandle_.product();
    event.getByToken(vtxScores3DToken_, vtxScores3DHandle_);
    auto vtxScores3D = *vtxScores3DHandle_.product();
    // Full 4D
    event.getByToken(vtx4DToken_, vtx4DHandle_);
    auto vtxs4D = *vtx4DHandle_.product();
    event.getByToken(vtxScores4DToken_, vtxScores4DHandle_);
    auto vtxScores4D = *vtxScores4DHandle_.product();
    // 4D no PID
    event.getByToken(vtx4DNoPIDToken_, vtx4DNoPIDHandle_);
    auto vtxs4DNoPID = *vtx4DNoPIDHandle_.product();

    //---gen sumpt of charged particle within |eta| acceptance
    int n_gen_charged=0;
    for(auto& part : genParticles)
    {
        if(part.status() == 1  && part.charge() != 0 && std::abs(part.eta())<4)
        {
            vtxsTree_.gen_sumpt += part.pt();
            if(part.pt()>0.7)
                ++n_gen_charged;
        }        
    }
    
    //---tracks loop
    //  - Check PID hypothesis, and vertex association comparing 4D and 4DnoPID
    std::map<int, reco::GenParticleRef> trackToGenPartMap;
    float vtx3D_0_dz_sumpt=0, vtx4D_0_dzdt_sumpt=0;
    float vtx3D_0_dz_genm_sumpt=0, vtx4D_0_dzdt_genm_sumpt=0;
    float vtx3D_0_3sigma_sumpt=0, vtx4D_0_3sigma_sumpt=0, vtx4D_0_circle3s_sumpt=0;
    float vtx3D_0_3sigma_genm_sumpt=0, vtx4D_0_3sigma_genm_sumpt=0, vtx4D_0_circle3s_genm_sumpt=0;
    for(unsigned int itrack=0; itrack<tracks.size(); ++itrack)
    {
        
        auto& track = tracks[itrack];
        reco::TrackBaseRef track_ref(tracksHandle_, itrack);
        reco::TrackBaseRef ext_track_ref;
        if(generalToExtendedTrkMap[track_ref] != -1)
            ext_track_ref = reco::TrackBaseRef(extTracksHandle_, generalToExtendedTrkMap[track_ref]);

        if(track.pt() < 0.7)
            continue;
        
        //---get TOFPIDProducer recomputed t0,sigmat0 and particle probs
        float t_t0=-1, t_sigmat0=-1, t_probPi=-1, t_probP=-1, t_probK=-1;
        if(t0PID.contains(track_ref.id()))
        {
            t_t0 = t0PID[track_ref];
            t_sigmat0 = sigmat0PID[track_ref]>0 ? 0.035 : -1; // FIXME once the correct value will be available in the upstream samples
            t_probPi = probPiPID[track_ref];
            t_probP = probPPID[track_ref];
            t_probK = probKPID[track_ref];
        }

        //---Match with tracking particles
        if(trkRecoToSimMap.find(track_ref) != trkRecoToSimMap.end())
        {
            auto& simTrack = trkRecoToSimMap[track_ref].begin()->first;
            trksTree_.trk_simIsFromPV->push_back(true);
            trksTree_.trk_simPt->push_back(simTrack->pt());
            trksTree_.trk_simEta->push_back(simTrack->eta());
            trksTree_.trk_simPhi->push_back(simTrack->phi());
            trksTree_.trk_simZ->push_back(simTrack->vz());
	}
        else
        {
            bool sim_found=false;
            float minDR = 1e6;
            TrackingParticle matched_sim_trk;
            for(auto& simtrk : simParticles)
            {
                auto dr = deltaR(track.eta(), track.phi(), simtrk.eta(), simtrk.phi());
                if(dr < minDR)
                {
                    minDR = dr;
                    matched_sim_trk = simtrk;
                    sim_found = true;
                }                
            }
            if(!sim_found)
            {
                trksTree_.trk_simIsFromPV->push_back(false);
                trksTree_.trk_simPt->push_back(-1);
                trksTree_.trk_simEta->push_back(0);
                trksTree_.trk_simPhi->push_back(0);
                trksTree_.trk_simZ->push_back(0);
            }
            else
            {
                trksTree_.trk_simIsFromPV->push_back(false);
                trksTree_.trk_simPt->push_back(matched_sim_trk.pt());
                trksTree_.trk_simEta->push_back(matched_sim_trk.eta());
                trksTree_.trk_simPhi->push_back(matched_sim_trk.phi());
                trksTree_.trk_simZ->push_back(matched_sim_trk.vz());
                trksTree_.trk_simDR->push_back(minDR);
            }                
        }

        // match with gen particles
        float DRMin = 1e9;
        int genPdgId = 0;
        float genEta = -999.;
        float genPhi = -999.;
        float genPt = -999.;        
        for(unsigned int iPart=0; iPart<genParticles.size(); ++iPart)
        {
            auto& genPart = genParticles[iPart];
            
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

        //---get MTD hits from pattern
        const auto& pattern = track.hitPattern();
        
        trksTree_.trk_idx -> push_back(itrack);
        trksTree_.trk_pt -> push_back(track.pt());
        trksTree_.trk_eta -> push_back(track.eta());
        trksTree_.trk_phi -> push_back(track.phi());
        trksTree_.trk_dxy_3D -> push_back(track.dxy(vtxs3D[0].position()));
        trksTree_.trk_dxy_4D -> push_back(track.dxy(vtxs4D[0].position()));        
        trksTree_.trk_dz_3D -> push_back(track.dz(vtxs3D[0].position()));
        trksTree_.trk_dz_4D -> push_back(track.dz(vtxs4D[0].position()));        
        trksTree_.trk_z -> push_back(track.vz());
        trksTree_.trk_sigmaxy -> push_back(track.dxyError());        
        trksTree_.trk_sigmaz -> push_back(track.dzError());
        trksTree_.trk_mtdt -> push_back(!ext_track_ref.isNull() ? extTracksMTDtime[ext_track_ref] : -1);
        trksTree_.trk_path_len -> push_back(!ext_track_ref.isNull() ? extPathLenght[ext_track_ref] : -1);
        trksTree_.trk_PID_t -> push_back(t_t0);
        trksTree_.trk_PID_sigmat -> push_back(t_sigmat0);
        trksTree_.trk_PID_probPi -> push_back(t_probPi);
        trksTree_.trk_PID_probP -> push_back(t_probP);
        trksTree_.trk_PID_probK -> push_back(t_probK);        
        trksTree_.trk_chi2 -> push_back(track.chi2());
        trksTree_.trk_ndof -> push_back(track.ndof());
        trksTree_.trk_isHighPurity -> push_back(track.quality(reco::TrackBase::TrackQuality::highPurity));
        trksTree_.trk_numberOfValidHits -> push_back(track.numberOfValidHits());
        trksTree_.trk_numberOfLostHits -> push_back(track.numberOfLostHits());
        trksTree_.trk_numberOfValidPixelBarrelHits -> push_back(pattern.numberOfValidPixelBarrelHits());
        trksTree_.trk_numberOfValidPixelEndcapHits -> push_back(pattern.numberOfValidPixelEndcapHits());        
        trksTree_.trk_numberOfValidHitsBTL -> push_back(pattern.numberOfValidTimingBTLHits());
        trksTree_.trk_numberOfValidHitsETL -> push_back(pattern.numberOfValidTimingETLHits());
        trksTree_.trk_hasMTD -> push_back(track.isTimeOk());
        trksTree_.trk_puid_3D -> push_back(mva3D_(track_ref, vtxs3D[0]));
        trksTree_.trk_puid_4D -> push_back(mva4D_(track_ref, ext_track_ref, vtxs4D[0],
                                                  t0PID, sigmat0PID, btlMatchChi2, btlMatchTimeChi2, etlMatchChi2, etlMatchTimeChi2,
                                                  extTracksMTDtime, extPathLenght));
        trksTree_.trk_genPdgId -> push_back(genPdgId);
        trksTree_.trk_genPt -> push_back(genPt);
        trksTree_.trk_genEta -> push_back(genEta);
        trksTree_.trk_genPhi -> push_back(genPhi);
        trksTree_.trk_genDR -> push_back(DRMin);
        trksTree_.trk_genVtx_x -> push_back(genPV.x());
        trksTree_.trk_genVtx_y -> push_back(genPV.y());        
        trksTree_.trk_genVtx_z -> push_back(genPV.z());        
        trksTree_.trk_genVtx_t -> push_back(genPV.t());

        trksTree_.pv4d_valid->push_back(vtxs4D[0].isValid() && !vtxs4D[0].isFake());
        trksTree_.pv4d_z->push_back(vtxs4D[0].z());
        trksTree_.pv4d_t->push_back(vtxs4D[0].t());
        trksTree_.pv4d_chi2->push_back(vtxs4D[0].chi2()/vtxs4D[0].ndof());
        trksTree_.pv4d_ntrks->push_back(vtxs4D[0].nTracks());
        trksTree_.pv3d_valid->push_back(vtxs3D[0].isValid() && !vtxs3D[0].isFake());
        trksTree_.pv3d_z->push_back(vtxs3D[0].z());
        trksTree_.pv3d_t->push_back(vtxs3D[0].t());
        trksTree_.pv3d_chi2->push_back(vtxs3D[0].chi2()/vtxs3D[0].ndof());
        trksTree_.pv3d_ntrks->push_back(vtxs3D[0].nTracks());
        
        //---compute sumpt for vtx0 of 3D and 4D collections (sigmat is wrong in CMSSW < 10_4_0_mtd5 so set i manually to 35ps)        
        if(std::abs(track.vz()-vtxs4D[0].z()) < 0.1 && (std::abs(t_t0-vtxs4D[0].t())<(3*0.035) || t_sigmat0 == -1))
        {
            vtx4D_0_dzdt_sumpt += track.pt();
            if(std::abs(genPt/track.pt()-1) < 0.05 && DRMin < 0.03)
                vtx4D_0_dzdt_genm_sumpt += track.pt();
        }
        if(std::abs(track.vz()-vtxs3D[0].z()) < 0.1)
        {
            vtx3D_0_dz_sumpt += track.pt();
            if(std::abs(genPt/track.pt()-1) < 0.05 && DRMin < 0.03)
                vtx3D_0_dz_genm_sumpt += track.pt();
        }

        //---compute dz,dt significance sumpt
        if(std::abs(track.vz()-vtxs4D[0].z()) < 3*sqrt(pow(track.dzError(), 2) + pow(vtxs4D[0].zError(),2)) &&
           (std::abs(t_t0-vtxs4D[0].t())<(3*0.035) || t_sigmat0 == -1))
        {
            vtx4D_0_3sigma_sumpt += track.pt();
            if(std::abs(genPt/track.pt()-1) < 0.05 && DRMin < 0.03)
                vtx4D_0_3sigma_genm_sumpt += track.pt();
        }
        if(std::abs(track.vz()-vtxs3D[0].z()) < 3*sqrt(pow(track.dzError(), 2) + pow(vtxs3D[0].zError(),2)))
        {
            vtx3D_0_3sigma_sumpt += track.pt();
            if(std::abs(genPt/track.pt()-1) < 0.05 && DRMin < 0.03)
                vtx3D_0_3sigma_genm_sumpt += track.pt();
        }
        //---3 sigma cut on (dz/sigma_z)^2+(dt/sigma_t)^2
        float dzs2 = pow(track.vz()-vtxs4D[0].z(), 2)/(pow(track.dzError(), 2) + pow(vtxs4D[0].zError(),2));
        float dts2 = pow((t_t0-vtxs4D[0].t())/0.035, 2);
        if((dzs2 + dts2)<9 || (dzs2<9 && t_sigmat0 == -1))
        {
            vtx4D_0_circle3s_sumpt += track.pt();
            if(std::abs(genPt/track.pt()-1) < 0.05 && DRMin < 0.03)
                vtx4D_0_circle3s_genm_sumpt += track.pt();
        }
    }

    trksTree_.n_gen_charged = n_gen_charged;        
    
    vtxsTree_.vtx3D_0_dz_sumpt = vtx3D_0_dz_sumpt;
    vtxsTree_.vtx3D_0_dz_genm_sumpt = vtx3D_0_dz_genm_sumpt;
    vtxsTree_.vtx4D_0_dzdt_sumpt = vtx4D_0_dzdt_sumpt;
    vtxsTree_.vtx4D_0_dzdt_genm_sumpt = vtx4D_0_dzdt_genm_sumpt;
    vtxsTree_.vtx3D_0_3sigma_sumpt = vtx3D_0_3sigma_sumpt;
    vtxsTree_.vtx3D_0_3sigma_genm_sumpt = vtx3D_0_3sigma_genm_sumpt;
    vtxsTree_.vtx4D_0_3sigma_sumpt = vtx4D_0_3sigma_sumpt;
    vtxsTree_.vtx4D_0_3sigma_genm_sumpt = vtx4D_0_3sigma_genm_sumpt;
    vtxsTree_.vtx4D_0_circle3s_sumpt = vtx4D_0_circle3s_sumpt;
    vtxsTree_.vtx4D_0_circle3s_genm_sumpt = vtx4D_0_circle3s_genm_sumpt;
        
    //---sim vertices
    //vtxsTree_.sim_n_vtxs = vtxsSim.size();
    
    //---3D vertices
    float min_dz=1e6, min_dzt=1e6;
    int best3D_idx=-1;
    int n_good_vtxs=0;
    for(unsigned int iVtx=0; iVtx<vtxs3D.size(); ++iVtx)
    {
        auto& vtx3D = vtxs3D[iVtx];
        reco::VertexRef vtx3D_ref(vtx3DHandle_, iVtx);
        if(vtx3D.isValid() && !vtx3D.isFake())
        {
            auto dz = std::abs(vtx3D.z()-genPV.z());
            if(dz < min_dz)
            {
                min_dz = dz;
                best3D_idx = iVtx;
            }
            ++n_good_vtxs;
        }

        //---Fill output tree
        ++vtxsTree_.vtx3D_n_vtxs;
        vtxsTree_.vtx3D_valid->push_back(vtxs3D[iVtx].isValid() && !vtxs3D[iVtx].isFake());
        vtxsTree_.vtx3D_z->push_back(vtxs3D[iVtx].z());
        vtxsTree_.vtx3D_t->push_back(-999);
        vtxsTree_.vtx3D_chi2->push_back(vtxs3D[iVtx].chi2()/vtxs3D[iVtx].ndof());
        vtxsTree_.vtx3D_ntrks->push_back(vtxs3D[iVtx].nTracks());

        //---count number of gen matched tracks associated to this vertex (trackWeight > 0.5)
        int assoc_trks=0;
        float sumpt=0, sumpt_genmatch=0;
        for(auto trkRef=vtx3D.tracks_begin(); trkRef!=vtx3D.tracks_end(); ++trkRef)
        {
            if(vtx3D.trackWeight(*trkRef) > 0.5)
            {
                if(trackToGenPartMap.find(trkRef->key()) != trackToGenPartMap.end())
                {
                    auto& genPart = trackToGenPartMap[trkRef->key()];                
                    if(std::abs(genPart->pt()/trkRef->get()->pt()-1) < 0.05 &&
                       deltaR(trkRef->get()->eta(), trkRef->get()->phi(), genPart->eta(), genPart->phi())<0.03)
                    {
                        ++assoc_trks;
                        sumpt_genmatch += trkRef->get()->pt();
                    }
                }
                sumpt += trkRef->get()->pt();
            }
        }
        vtxsTree_.vtx3D_ntrks_genmatch->push_back(assoc_trks);
        vtxsTree_.vtx3D_sumpt_genmatch->push_back(sumpt_genmatch);        
        vtxsTree_.vtx3D_sumpt->push_back(sumpt);
        vtxsTree_.vtx3D_score_sumpt2->push_back(vtxScores3D[vtx3D_ref]);
    }
    vtxsTree_.vtx3D_n_good = n_good_vtxs;
    vtxsTree_.vtx3D_best_dz = best3D_idx;
    
    //---4D vertices
    min_dz=1e6, min_dzt=1e6;
    int best4D_dz_idx=-1, best4D_dzt_idx=-1;
    n_good_vtxs=0;
    for(unsigned int iVtx=0; iVtx<vtxs4D.size(); ++iVtx)
    {
        auto& vtx4D = vtxs4D[iVtx];
        reco::VertexRef vtx4D_ref(vtx4DHandle_, iVtx);
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
            ++n_good_vtxs;
        }
        
        //---Fill output tree
        ++vtxsTree_.vtx4D_n_vtxs;
        vtxsTree_.vtx4D_valid->push_back(vtxs4D[iVtx].isValid() && !vtxs4D[iVtx].isFake());
        vtxsTree_.vtx4D_z->push_back(vtxs4D[iVtx].z());
        vtxsTree_.vtx4D_t->push_back(vtxs4D[iVtx].t());
        vtxsTree_.vtx4D_chi2->push_back(vtxs4D[iVtx].chi2()/vtxs4D[iVtx].ndof());
        vtxsTree_.vtx4D_ntrks->push_back(vtxs4D[iVtx].nTracks());

        //---count number of gen matched tracks associated to this vertex (trackWeight > 0.5)
        int assoc_trks=0;
        float sumpt=0, sumpt_genmatch=0;
        for(auto trkRef=vtx4D.tracks_begin(); trkRef!=vtx4D.tracks_end(); ++trkRef)
        {
            if(vtx4D.trackWeight(*trkRef) > 0.5)
            {                    
                if(trackToGenPartMap.find(trkRef->key()) != trackToGenPartMap.end())
                {
                    auto& genPart = trackToGenPartMap[trkRef->key()];                
                    if(std::abs(genPart->pt()/trkRef->get()->pt()-1) < 0.05 &&
                       deltaR(trkRef->get()->eta(), trkRef->get()->phi(), genPart->eta(), genPart->phi())<0.03)
                    {
                        ++assoc_trks;
                        sumpt_genmatch += trkRef->get()->pt();
                    }
                }
                sumpt += trkRef->get()->pt();
            }
        }
        vtxsTree_.vtx4D_ntrks_genmatch->push_back(assoc_trks);
        vtxsTree_.vtx4D_sumpt_genmatch->push_back(sumpt_genmatch);        
        vtxsTree_.vtx4D_sumpt->push_back(sumpt);
        vtxsTree_.vtx4D_score_sumpt2->push_back(vtxScores4D[vtx4D_ref]);
    }
    vtxsTree_.vtx4D_n_good = n_good_vtxs;
    vtxsTree_.vtx4D_best_dz = best4D_dz_idx;
    vtxsTree_.vtx4D_best_dzt = best4D_dzt_idx;

    //---4DNoPID vertices
    min_dz=1e6, min_dzt=1e6;
    int best4DNoPID_dz_idx=-1, best4DNoPID_dzt_idx=-1;
    n_good_vtxs=0;
    for(unsigned int iVtx=0; iVtx<vtxs4DNoPID.size(); ++iVtx)
    {
        auto& vtx4DNoPID = vtxs4DNoPID[iVtx];
        if(vtx4DNoPID.isValid() && !vtx4DNoPID.isFake())
        {
            auto dz = std::abs(vtx4DNoPID.z()-genPV.z());
            auto dzt = sqrt(dz*dz + pow((vtx4DNoPID.t()-genPV.t())*30, 2));
            if(dz < min_dz)
            {
                min_dz = dz;
                best4DNoPID_dz_idx = iVtx;
            }
            if(dzt < min_dzt)
            {
                min_dzt = dzt;
                best4DNoPID_dzt_idx = iVtx;
            }
            ++n_good_vtxs;
        }
        
        //---Fill output tree
        ++vtxsTree_.vtx4DNoPID_n_vtxs;
        vtxsTree_.vtx4DNoPID_valid->push_back(vtxs4DNoPID[iVtx].isValid() && !vtxs4DNoPID[iVtx].isFake());
        vtxsTree_.vtx4DNoPID_z->push_back(vtxs4DNoPID[iVtx].z());
        vtxsTree_.vtx4DNoPID_t->push_back(vtxs4DNoPID[iVtx].t());
        vtxsTree_.vtx4DNoPID_chi2->push_back(vtxs4DNoPID[iVtx].chi2()/vtxs4DNoPID[iVtx].ndof());
        vtxsTree_.vtx4DNoPID_ntrks->push_back(vtxs4DNoPID[iVtx].nTracks());

        //---count number of gen matched tracks associated to this vertex (trackWeight > 0.5)
        int assoc_trks=0;
        float sumpt=0, sumpt_genmatch=0;
        for(auto trkRef=vtx4DNoPID.tracks_begin(); trkRef!=vtx4DNoPID.tracks_end(); ++trkRef)
        {
            if(vtx4DNoPID.trackWeight(*trkRef) > 0.5)
            {                    
                if(trackToGenPartMap.find(trkRef->key()) != trackToGenPartMap.end())
                {
                    auto& genPart = trackToGenPartMap[trkRef->key()];                
                    if(std::abs(genPart->pt()/trkRef->get()->pt()-1) < 0.05 &&
                       deltaR(trkRef->get()->eta(), trkRef->get()->phi(), genPart->eta(), genPart->phi())<0.03)
                    {
                        ++assoc_trks;
                        sumpt_genmatch += trkRef->get()->pt();
                    }
                }
                sumpt += trkRef->get()->pt();
            }
        }
        vtxsTree_.vtx4DNoPID_ntrks_genmatch->push_back(assoc_trks);
        vtxsTree_.vtx4DNoPID_sumpt_genmatch->push_back(sumpt_genmatch);        
        vtxsTree_.vtx4DNoPID_sumpt->push_back(sumpt);
    }
    vtxsTree_.vtx4DNoPID_n_good = n_good_vtxs;
    vtxsTree_.vtx4DNoPID_best_dz = best4DNoPID_dz_idx;
    vtxsTree_.vtx4DNoPID_best_dzt = best4DNoPID_dzt_idx;

    //---Gen PV 
    vtxsTree_.gen_z = genPV.z();
    vtxsTree_.gen_t = genPV.t();        

    //---Fill trees
    vtxsTree_.GetTTreePtr()->Fill();
    trksTree_.GetTTreePtr()->Fill();
}


DEFINE_FWK_MODULE(MTD4DVertexingAnalyzer);

#endif
