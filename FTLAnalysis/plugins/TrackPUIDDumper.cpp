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
    MTDTrackPUIDTree trksTree_;
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
    vtx4DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DTag")))
{
    trksTree_ = MTDTrackPUIDTree(pSet.getUntrackedParameter<string>("trksTreeName").c_str(), "4D TOFPID studies");
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
        trksTree_.Reset();
        
        auto& track = extTracks[itrack];
        reco::TrackBaseRef track_ref(tracksHandle_, itrack);
        reco::TrackBaseRef ext_track_ref(extTracksHandle_, itrack);

        //---Remove tracks with very low pt
        if(track.pt() < 0.5)
            continue;
        
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
        if(trkRecoToSimMap.find(track_ref) != trkRecoToSimMap.end())
        {
            auto& simTrack = trkRecoToSimMap[track_ref].begin()->first;
            trksTree_.simIsFromPV = true;
            trksTree_.simPt = simTrack->pt();
            trksTree_.simEta = simTrack->eta();
            trksTree_.simPhi = simTrack->phi();
            trksTree_.simZ = simTrack->vz();
        
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
	}
        else
        {
            trksTree_.simIsFromPV = false;
            trksTree_.simPt = -1;
            trksTree_.simEta = 0;
            trksTree_.simPhi = 0;
            trksTree_.simZ = 0;
        }
        
        //---get MTD hits from pattern
        const auto& pattern = track.hitPattern();        

        trksTree_.idx = itrack;
        trksTree_.pt = track.pt();
        trksTree_.eta = track.eta();
        trksTree_.phi = track.phi();
        trksTree_.x = track.vx();
        trksTree_.y = track.vy();
        trksTree_.z = track.vz();
        trksTree_.dzErr = track.dzError();
        trksTree_.dxyErr = track.dxyError();
        trksTree_.chi2 = track.chi2();
        trksTree_.ndof = track.ndof();
        trksTree_.t0 = t_t0;
        trksTree_.sigmat0 = t_sigmat0;
        trksTree_.mtdt = extTracksMTDtime[ext_track_ref];
        trksTree_.path_len = extPathLenght[ext_track_ref];
        trksTree_.probPi = t_probPi;
        trksTree_.probP = t_probP;
        trksTree_.probK = t_probK;        
        trksTree_.energy = sqrt(track.momentum().mag2());
        trksTree_.btlMatchChi2 = btlMatchChi2.contains(ext_track_ref.id()) ? btlMatchChi2[ext_track_ref] : -1;
        trksTree_.btlMatchTimeChi2 = btlMatchTimeChi2.contains(ext_track_ref.id()) ? btlMatchTimeChi2[ext_track_ref] : -1;
        trksTree_.etlMatchChi2 = etlMatchChi2.contains(ext_track_ref.id()) ? etlMatchChi2[ext_track_ref] : -1;
        trksTree_.etlMatchTimeChi2 = etlMatchTimeChi2.contains(ext_track_ref.id()) ? etlMatchTimeChi2[ext_track_ref] : -1;        
        trksTree_.normalizedChi2 = track.normalizedChi2();
        trksTree_.numberOfValidHits = track.numberOfValidHits();
        trksTree_.numberOfLostHits = track.numberOfLostHits();
        trksTree_.numberOfValidPixelBarrelHits = pattern.numberOfValidPixelBarrelHits();
        trksTree_.numberOfValidPixelEndcapHits = pattern.numberOfValidPixelEndcapHits();
        trksTree_.numberOfValidStripTIBHits = pattern.numberOfValidStripTIBHits();
        trksTree_.numberOfValidStripTIDHits = pattern.numberOfValidStripTIDHits();
        trksTree_.numberOfValidStripTOBHits = pattern.numberOfValidStripTOBHits();
        trksTree_.numberOfValidStripTECHits = pattern.numberOfValidStripTECHits();
        trksTree_.numberOfValidHitsBTL = pattern.numberOfValidTimingBTLHits();
        trksTree_.numberOfValidHitsETL = pattern.numberOfValidTimingETLHits();
        trksTree_.isHighPurity = track.quality(reco::TrackBase::TrackQuality::highPurity);
        trksTree_.hasMTD = track.isTimeOk();
        trksTree_.genPdgId = genPdgId;
        trksTree_.genPt = genPt;
        trksTree_.genEta = genEta;
        trksTree_.genPhi = genPhi;
        trksTree_.genDR = DRMin;
        trksTree_.genVtx_x = genPV.x();
        trksTree_.genVtx_y = genPV.y();        
        trksTree_.genVtx_z = genPV.z();        
        trksTree_.genVtx_t = genPV.t();
        trksTree_.pv3d_valid = vtxs3D[0].isValid() && !vtxs4D[0].isFake();
        trksTree_.pv3d_ntrks = vtxs3D[0].nTracks();
        trksTree_.pv3d_chi2 = vtxs3D[0].chi2()/vtxs3D[0].ndof();
        trksTree_.pv3d_x = vtxs3D[0].x();
        trksTree_.pv3d_y = vtxs3D[0].y();
        trksTree_.pv3d_z = vtxs3D[0].z();
        trksTree_.pv4d_valid = vtxs4D[0].isValid() && !vtxs4D[0].isFake();
        trksTree_.pv4d_ntrks = vtxs4D[0].nTracks();
        trksTree_.pv4d_chi2 = vtxs4D[0].chi2()/vtxs4D[0].ndof();
        trksTree_.pv4d_x = vtxs4D[0].x();
        trksTree_.pv4d_y = vtxs4D[0].y();
        trksTree_.pv4d_z = vtxs4D[0].z();
        trksTree_.pv4d_t = vtxs4D[0].t();
        // trksTree_.puid_3D = mva3D_(track_ref, vtxs3D[0]);
        // trksTree_.puid_4D = mva4D_(track_ref, ext_track_ref, vtxs4D[0],
        //                            t0PID, sigmat0PID, btlMatchChi2, btlMatchTimeChi2, etlMatchChi2, etlMatchTimeChi2,
        //                            extTracksMTDtime, extPathLenght);
        
        //---Fill tree
        trksTree_.GetTTreePtr()->Fill();
    }
    
}


DEFINE_FWK_MODULE(TrackPUIDDumper);

#endif
