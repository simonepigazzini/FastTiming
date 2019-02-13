#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "PrecisionTiming/FTLAnalysis/interface/MVAComputer.h"

class TrackPUIDMVA
{
public:
    //---ctors---
    TrackPUIDMVA(std::string weights_file, bool is4D=false);

    //---dtor---
    ~TrackPUIDMVA() {};

    //---getters---
    // 4D
    float operator() (
        reco::TrackBaseRef& trk, reco::TrackBaseRef& ext_trk, reco::Vertex& vtx,
        edm::ValueMap<float>& t0s,
        edm::ValueMap<float>& sigma_t0s,
        edm::ValueMap<float>& btl_chi2s,
        edm::ValueMap<float>& btl_time_chi2s,                      
        edm::ValueMap<float>& etl_chi2s,
        edm::ValueMap<float>& etl_time_chi2s,
        edm::ValueMap<float>& tmtds,
        edm::ValueMap<float>& trk_lengths
        );

    // 3D
    float operator() (reco::TrackBaseRef& trk, reco::Vertex& vtx);

private:
    bool                       is4D_;
    MVAComputer::mva_variables vars_;
    MVAComputer                mva_;
};
    
