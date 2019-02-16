#ifndef __TRACK_PU_ID_DUMPER__
#define __TRACK_PU_ID_DUMPER__

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
#include "PrecisionTiming/FTLAnalysis/interface/MTDTrackPUIDTree.h"

using namespace std;                             

class TrackPUIDDumper : public edm::EDAnalyzer
{
public:                             

    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;
    
    explicit TrackPUIDDumper(const edm::ParameterSet& pSet);
    ~TrackPUIDDumper() {};

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
    edm::EDGetTokenT<edm::ValueMap<float> > btlMatchChi2Token_;
    edm::Handle<edm::ValueMap<float> > btlMatchChi2Handle_;
    edm::EDGetTokenT<edm::ValueMap<float> > btlMatchTimeChi2Token_;
    edm::Handle<edm::ValueMap<float> > btlMatchTimeChi2Handle_;
    edm::EDGetTokenT<edm::ValueMap<float> > etlMatchChi2Token_;
    edm::Handle<edm::ValueMap<float> > etlMatchChi2Handle_;
    edm::EDGetTokenT<edm::ValueMap<float> > etlMatchTimeChi2Token_;
    edm::Handle<edm::ValueMap<float> > etlMatchTimeChi2Handle_;
    edm::EDGetTokenT<edm::ValueMap<float> > extTracksMTDtimeToken_;
    edm::Handle<edm::ValueMap<float> > extTracksMTDtimeHandle_;
    edm::EDGetTokenT<edm::ValueMap<float> > pathLengthToken_;
    edm::Handle<edm::ValueMap<float> > pathLengthHandle_;
    
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
    edm::EDGetTokenT<vector<reco::Vertex> >   vtx3DToken_;
    edm::Handle<vector<reco::Vertex> >        vtx3DHandle_;
    edm::EDGetTokenT<vector<reco::Vertex> >   vtx4DToken_;
    edm::Handle<vector<reco::Vertex> >        vtx4DHandle_;    
    
    //---options

    //---outputs
    TrackPUIDMVA mva3D_;
    TrackPUIDMVA mva4D_;    
    MTDTrackPUIDTree trks3DTree_;
    MTDTrackPUIDTree trks4DTree_;    
    edm::Service<TFileService> fs_;  
  
};

TrackPUIDDumper::TrackPUIDDumper(const edm::ParameterSet& pSet):
    genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
    simParticlesToken_(consumes<TrackingParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackingParticlesTag"))),
    trkRecoToSimMapToken_(consumes<reco::RecoToSimCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackAndTrackingParticlesAssociatorMapTag"))),
    trkSimToRecoMapToken_(consumes<reco::SimToRecoCollection>(pSet.getUntrackedParameter<edm::InputTag>("trackAndTrackingParticlesAssociatorMapTag"))),
    tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("generalTracksTag"))),
    extTracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("extendedTracksTag"))),
    btlMatchChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("btlMatchChi2Tag"))),
    btlMatchTimeChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("btlMatchTimeChi2Tag"))),
    etlMatchChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("etlMatchChi2Tag"))),
    etlMatchTimeChi2Token_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("etlMatchTimeChi2Tag"))) ,
    extTracksMTDtimeToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("extTracksMTDtimeTag"))),
    pathLengthToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("extTracksPathLengthTag"))),
    t0PIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("t0TOFPIDTag"))),
    sigmat0PIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("sigmat0TOFPIDTag"))),
    probPiPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probPiTOFPIDTag"))),
    probPPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probPTOFPIDTag"))),
    probKPIDToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("probKTOFPIDTag"))),    
    genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
    genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
    vtx3DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx3DTag"))),
    vtx4DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DTag"))),
    mva3D_(pSet.getParameter<edm::FileInPath>("trackPUID_3DBDT_weights_file").fullPath(), false),
    mva4D_(pSet.getParameter<edm::FileInPath>("trackPUID_4DBDT_weights_file").fullPath(), true)
{
    trks3DTree_ = MTDTrackPUIDTree((pSet.getUntrackedParameter<string>("trksTreeName")+"3D").c_str(), "Track PU-ID training tree");
    trks4DTree_ = MTDTrackPUIDTree((pSet.getUntrackedParameter<string>("trksTreeName")+"4D").c_str(), "Track PU-ID training tree");    
}

void TrackPUIDDumper::analyze(edm::Event const& event, edm::EventSetup const& setup)
{ 
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

    //---extended tracks path length
    event.getByToken(btlMatchChi2Token_, btlMatchChi2Handle_);
    auto btlMatchChi2 = *btlMatchChi2Handle_.product();
    event.getByToken(btlMatchTimeChi2Token_, btlMatchTimeChi2Handle_);
    auto btlMatchTimeChi2 = *btlMatchTimeChi2Handle_.product();
    event.getByToken(etlMatchChi2Token_, etlMatchChi2Handle_);
    auto etlMatchChi2 = *etlMatchChi2Handle_.product();
    event.getByToken(etlMatchTimeChi2Token_, etlMatchTimeChi2Handle_);
    auto etlMatchTimeChi2 = *etlMatchTimeChi2Handle_.product();
    event.getByToken(pathLengthToken_, pathLengthHandle_);
    auto extPathLenght = *pathLengthHandle_.product();
    event.getByToken(extTracksMTDtimeToken_, extTracksMTDtimeHandle_);
    auto extTracksMTDtime = *extTracksMTDtimeHandle_.product();
    
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
    // 3D
    event.getByToken(vtx3DToken_, vtx3DHandle_);
    auto vtxs3D = *vtx3DHandle_.product();
    // Full 4D
    event.getByToken(vtx4DToken_, vtx4DHandle_);
    auto vtxs4D = *vtx4DHandle_.product();

    for(unsigned int itrack=0; itrack<extTracks.size(); ++itrack)
    {
        trks3DTree_.Reset();
        trks4DTree_.Reset();

        auto& track = tracks[itrack];
        auto& ext_track = extTracks[itrack];
        reco::TrackBaseRef track_ref(tracksHandle_, itrack);
        reco::TrackBaseRef ext_track_ref(extTracksHandle_, itrack);
        
        //---Remove tracks with very low pt
        if(track.pt() < 0.4)
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

        // match with gen particles
        float DRMin = 1e9;
        int genPdgId = 0;
        float genEta = -999.;
        float genPhi = -999.;
        float genPt = -999.;        
        if(trkRecoToSimMap.find(track_ref) != trkRecoToSimMap.end())
        {
            auto& simTrack = trkRecoToSimMap[track_ref].begin()->first;
            trks3DTree_.simIsFromPV = true;
            trks3DTree_.simPt = simTrack->pt();
            trks3DTree_.simEta = simTrack->eta();
            trks3DTree_.simPhi = simTrack->phi();
            trks3DTree_.simZ = simTrack->vz(); 
            trks4DTree_.simIsFromPV = true;
            trks4DTree_.simPt = simTrack->pt();
            trks4DTree_.simEta = simTrack->eta();
            trks4DTree_.simPhi = simTrack->phi();
            trks4DTree_.simZ = simTrack->vz();

            for(unsigned int iPart=0; iPart<genParticles.size(); ++iPart)
            {
                auto& genPart = genParticles[iPart];
            
                if( genPart.status() != 1 ) continue;
                if( genPart.charge() == 0 ) continue;

                float DR   = deltaR(ext_track.eta(), ext_track.phi(), genPart.eta(), genPart.phi());
                
                if( DR < DRMin )
                {
                    DRMin = DR;
                    
                    genPdgId = genPart.pdgId();
                    genEta   = genPart.eta();
                    genPhi   = genPart.phi();
                    genPt    = genPart.pt();
                }
            }
	}
        else
        {
            trks3DTree_.simIsFromPV = false;
            trks3DTree_.simPt = -1;
            trks3DTree_.simEta = 0;
            trks3DTree_.simPhi = 0;
            trks3DTree_.simZ = 0;
            trks4DTree_.simIsFromPV = false;
            trks4DTree_.simPt = -1;
            trks4DTree_.simEta = 0;
            trks4DTree_.simPhi = 0;
            trks4DTree_.simZ = 0;
        }
        
        //---fill 3D tree
        const auto& pattern3D = track.hitPattern();        

        trks3DTree_.idx = itrack;
        trks3DTree_.pt = track.pt();
        trks3DTree_.eta = track.eta();
        trks3DTree_.phi = track.phi();
        trks3DTree_.x = track.vx();
        trks3DTree_.y = track.vy();
        trks3DTree_.z = track.vz();
        trks3DTree_.dz = track.dz(vtxs3D[0].position());
        trks3DTree_.dxy = track.dxy(vtxs3D[0].position());
        trks3DTree_.dzErr = track.dzError();
        trks3DTree_.dxyErr = track.dxyError();
        trks3DTree_.chi2 = track.chi2();
        trks3DTree_.ndof = track.ndof();
        trks3DTree_.t0 = -1;
        trks3DTree_.sigmat0 = -1;
        trks3DTree_.mtdt = -1;
        trks3DTree_.path_len = -1;
        trks3DTree_.probPi = -1;
        trks3DTree_.probP = -1;
        trks3DTree_.probK = -1;
        trks3DTree_.btlMatchChi2 = -1;
        trks3DTree_.btlMatchTimeChi2 = -1;
        trks3DTree_.etlMatchChi2 = -1;
        trks3DTree_.etlMatchTimeChi2 = -1;
        trks3DTree_.normalizedChi2 = track.normalizedChi2();
        trks3DTree_.numberOfValidHits = track.numberOfValidHits();
        trks3DTree_.numberOfLostHits = track.numberOfLostHits();
        trks3DTree_.numberOfValidPixelBarrelHits = pattern3D.numberOfValidPixelBarrelHits();
        trks3DTree_.numberOfValidPixelEndcapHits = pattern3D.numberOfValidPixelEndcapHits();
        trks3DTree_.numberOfValidHitsBTL = pattern3D.numberOfValidTimingBTLHits();
        trks3DTree_.numberOfValidHitsETL = pattern3D.numberOfValidTimingETLHits();
        trks3DTree_.hasMTD = track.isTimeOk();
        trks3DTree_.genPdgId = genPdgId;
        trks3DTree_.genPt = genPt;
        trks3DTree_.genEta = genEta;
        trks3DTree_.genPhi = genPhi;
        trks3DTree_.genDR = DRMin;
        trks3DTree_.genVtx_x = genPV.x();
        trks3DTree_.genVtx_y = genPV.y();        
        trks3DTree_.genVtx_z = genPV.z();        
        trks3DTree_.genVtx_t = genPV.t();
        trks3DTree_.pv_valid = vtxs3D[0].isValid() && !vtxs3D[0].isFake();
        trks3DTree_.pv_ntrks = vtxs3D[0].nTracks();
        trks3DTree_.pv_chi2 = vtxs3D[0].chi2()/vtxs3D[0].ndof();
        trks3DTree_.pv_x = vtxs3D[0].x();
        trks3DTree_.pv_y = vtxs3D[0].y();
        trks3DTree_.pv_z = vtxs3D[0].z();
        trks3DTree_.pv_t = -1.;
        trks3DTree_.puid = mva3D_(track_ref, vtxs3D[0]);
       
        //---fill 4D tree
        const auto& pattern4D = ext_track.hitPattern();        

        trks4DTree_.idx = itrack;
        trks4DTree_.pt = track.pt();
        trks4DTree_.eta = track.eta();
        trks4DTree_.phi = track.phi();
        trks4DTree_.x = track.vx();
        trks4DTree_.y = track.vy();
        trks4DTree_.z = track.vz();
        trks4DTree_.dz = track.dz(vtxs4D[0].position());
        trks4DTree_.dxy = track.dxy(vtxs4D[0].position());
        trks4DTree_.dzErr = track.dzError();
        trks4DTree_.dxyErr = track.dxyError();
        trks4DTree_.chi2 = track.chi2();
        trks4DTree_.ndof = track.ndof();
        trks4DTree_.t0 = t_t0;
        trks4DTree_.sigmat0 = t_sigmat0;
        trks4DTree_.mtdt = extTracksMTDtime.contains(ext_track_ref.id()) ? extTracksMTDtime[ext_track_ref] : -1;
        trks4DTree_.path_len = extPathLenght.contains(ext_track_ref.id()) ? extPathLenght[ext_track_ref] : -1;
        trks4DTree_.probPi = t_probPi;
        trks4DTree_.probP = t_probP;
        trks4DTree_.probK = t_probK;        
        trks4DTree_.btlMatchChi2 = btlMatchChi2.contains(ext_track_ref.id()) ? btlMatchChi2[ext_track_ref] : -1;
        trks4DTree_.btlMatchTimeChi2 = btlMatchTimeChi2.contains(ext_track_ref.id()) ? btlMatchTimeChi2[ext_track_ref] : -1;
        trks4DTree_.etlMatchChi2 = etlMatchChi2.contains(ext_track_ref.id()) ? etlMatchChi2[ext_track_ref] : -1;
        trks4DTree_.etlMatchTimeChi2 = etlMatchTimeChi2.contains(ext_track_ref.id()) ? etlMatchTimeChi2[ext_track_ref] : -1;        
        trks4DTree_.normalizedChi2 = ext_track.normalizedChi2();
        trks4DTree_.numberOfValidHits = track.numberOfValidHits();
        trks4DTree_.numberOfLostHits = track.numberOfLostHits();
        trks4DTree_.numberOfValidPixelBarrelHits = pattern4D.numberOfValidPixelBarrelHits();
        trks4DTree_.numberOfValidPixelEndcapHits = pattern4D.numberOfValidPixelEndcapHits();
        trks4DTree_.numberOfValidHitsBTL = pattern4D.numberOfValidTimingBTLHits();
        trks4DTree_.numberOfValidHitsETL = pattern4D.numberOfValidTimingETLHits();
        trks4DTree_.hasMTD = ext_track.isTimeOk();
        trks4DTree_.genPdgId = genPdgId;
        trks4DTree_.genPt = genPt;
        trks4DTree_.genEta = genEta;
        trks4DTree_.genPhi = genPhi;
        trks4DTree_.genDR = DRMin;
        trks4DTree_.genVtx_x = genPV.x();
        trks4DTree_.genVtx_y = genPV.y();        
        trks4DTree_.genVtx_z = genPV.z();        
        trks4DTree_.genVtx_t = genPV.t();
        trks4DTree_.pv_valid = vtxs4D[0].isValid() && !vtxs4D[0].isFake();
        trks4DTree_.pv_ntrks = vtxs4D[0].nTracks();
        trks4DTree_.pv_chi2 = vtxs4D[0].chi2()/vtxs4D[0].ndof();
        trks4DTree_.pv_x = vtxs4D[0].x();
        trks4DTree_.pv_y = vtxs4D[0].y();
        trks4DTree_.pv_z = vtxs4D[0].z();
        trks4DTree_.pv_t = vtxs4D[0].t();
        trks4DTree_.puid = mva4D_(track_ref, ext_track_ref, vtxs4D[0],
                                  t0PID, sigmat0PID, btlMatchChi2, btlMatchTimeChi2, etlMatchChi2, etlMatchTimeChi2,
                                  extTracksMTDtime, extPathLenght);
        
        //---Fill tree
        trks3DTree_.GetTTreePtr()->Fill();
        trks4DTree_.GetTTreePtr()->Fill();
    }
    
}


DEFINE_FWK_MODULE(TrackPUIDDumper);

#endif
