#ifndef _MTD_NEUTRALS_SUMET_TOY_
#define _MTD_NEUTRALS_SUMET_TOY_

#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"

#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "PrecisionTiming/FTLAnalysis/interface/MTDNeutralsTree.h"

using namespace std;

class MTDNeutralSumEtToy : public edm::EDAnalyzer
{
public:                             

    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;
    
    explicit MTDNeutralSumEtToy(const edm::ParameterSet& pSet);
    ~MTDNeutralSumEtToy() {};

    //---methods
    virtual void beginJob() override {};
    virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
    virtual void endJob() override {};
    
private:
    //---inputs
    
    //---PF
    edm::Handle<reco::PFCandidateCollection> pfCandidatesHandle_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;

    //---vertices
    edm::EDGetTokenT<genXYZ>                  genXYZToken_;
    edm::Handle<genXYZ>                       genXYZHandle_;
    edm::EDGetTokenT<float>                   genT0Token_;
    edm::Handle<float>                        genT0Handle_;
    
    //---options
    vector<double> timeResolutions_;
    
    //---outputs
    MTDNeutralsToyTree outTree_;
    edm::Service<TFileService> fs_;    
};

MTDNeutralSumEtToy::MTDNeutralSumEtToy(const edm::ParameterSet& pSet):
    pfCandidatesToken_(consumes<reco::PFCandidateCollection>(pSet.getUntrackedParameter<edm::InputTag>("pfCandidatesTag"))),    
    genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
    genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
    timeResolutions_(pSet.getUntrackedParameter<vector<double> >("timeResolutions"))
{
    outTree_ = MTDNeutralsToyTree(pSet.getUntrackedParameter<string>("outTreeName").c_str(), "Neutrals toy MC");
}

void MTDNeutralSumEtToy::analyze(edm::Event const& event, edm::EventSetup const& setup)
{
    //---load gen vertex
    event.getByToken(genXYZToken_, genXYZHandle_);
    event.getByToken(genT0Token_, genT0Handle_);
    auto xyz = genXYZHandle_.product();
    auto t = *genT0Handle_.product();
    auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
    auto genPV = SimVertex(v, t).position();

    //---NuGun only!
    // auto v = math::XYZVectorD(0, 0, gRandom->Gaus(0, 4.2));
    // auto genPV = SimVertex(v, gRandom->Gaus(0, 0.18)).position();
    
    //---load PFCandidates
    event.getByToken(pfCandidatesToken_, pfCandidatesHandle_);
    auto pfCandidates = *pfCandidatesHandle_.product();    

    //---efficiency parametrization from:
    //   https://spigazzi.web.cern.ch/spigazzi/precision_timing/TDR/neutrals/SingleGamma/0PU_chi2best/neutrals_efficiency_vs_eta.png
    TF1 param_neu_eff_vs_eta("param_neu_eff_vs_eta", "pol2", 0, 1.5);
    param_neu_eff_vs_eta.SetParameters(0.306374, -0.0155606, 0.121197);

    //---clean output tree
    outTree_.Reset();
    outTree_.vtx_z = genPV.z();
    outTree_.vtx_t = genPV.t();
    //outTree_.t_res = &timeResolutions_;
    outTree_.barrel_mtd_noeff_npho->resize(timeResolutions_.size(), 0);
    outTree_.barrel_mtd_noeff_sumet->resize(timeResolutions_.size(), 0.);
    outTree_.barrel_mtd_npho->resize(timeResolutions_.size(), 0);
    outTree_.barrel_mtd_sumet->resize(timeResolutions_.size(), 0.);
    
    for(auto& cand : pfCandidates)
    {
        //---looking at neutral in ECAL only
        if(cand.pt()>1. && cand.particleId()==4)
        {
            if(std::abs(cand.eta())<1.5)
            {                
                int ires=0;
                for(auto& t_res : timeResolutions_)
                {
                    //---eff + 3 sigma dt cut, for different time resolutions
                    auto eff = gRandom->Uniform();
                    auto dt = std::abs(gRandom->Gaus(0, 0.18)-genPV.t());
                    if((eff<=param_neu_eff_vs_eta.Eval(cand.eta()) && dt<3*t_res) ||
                       eff>param_neu_eff_vs_eta.Eval(cand.eta()))
                    {
                        outTree_.barrel_mtd_npho->at(ires) += 1;
                        outTree_.barrel_mtd_sumet->at(ires) += cand.pt();
                    }

                    //---100% eff, 3 sigma cut only
                    if(dt<3*t_res)
                    {
                        outTree_.barrel_mtd_noeff_npho->at(ires) += 1;
                        outTree_.barrel_mtd_noeff_sumet->at(ires) += cand.pt();
                    }
                    
                    ++ires;
                }

                //---sumEt without MTD
                ++outTree_.barrel_npho;
                outTree_.barrel_sumet += cand.pt();
            }
        }
    }

    outTree_.GetTTreePtr()->Fill();
}

DEFINE_FWK_MODULE(MTDNeutralSumEtToy);

#endif
