#ifndef _MTD_TOFPID_TREE_
#define _MTD_TOFPID_TREE_

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME MTDTOFPIDTree

#define DATA_TABLE                              \
    DATA(int, event)                            \
    DATA(int, lumi)                             \
    DATA(int, run)                              \
    DATA(float, vtx_0_valid)                    \
    DATA(float, vtx_0_z)                        \
    DATA(float, vtx_0_t)                        \
    DATA(float, vtx_0_chi2)                     \
    DATA(float, vtx_0_ntrks)                    \
    DATA(int, n_gen_charged)

#define DATA_CLASS_TABLE                                \
    DATA(vector<int>,   trk_idx)                        \
    DATA(vector<float>, trk_pt)                         \
    DATA(vector<float>, trk_eta)                        \
    DATA(vector<float>, trk_phi)                        \
    DATA(vector<float>, trk_x)                          \
    DATA(vector<float>, trk_y)                          \
    DATA(vector<float>, trk_z)                          \
    DATA(vector<float>, trk_sigmaz)                     \
    DATA(vector<float>, trk_t)                          \
    DATA(vector<float>, trk_sigmat)                     \
    DATA(vector<float>, trk_mtdt)                       \
    DATA(vector<float>, trk_path_len)                   \
    DATA(vector<float>, trk_PID_t)                      \
    DATA(vector<float>, trk_PID_sigmat)                 \
    DATA(vector<float>, trk_PID_probPi)                 \
    DATA(vector<float>, trk_PID_probP)                  \
    DATA(vector<float>, trk_PID_probK)                  \
    DATA(vector<float>, trk_energy)                     \
    DATA(vector<float>, trk_normalizedChi2)             \
    DATA(vector<int>,   trk_numberOfValidHits)          \
    DATA(vector<int>,   trk_numberOfLostHits)           \
    DATA(vector<int>,   trk_numberOfValidHitsBTL)       \
    DATA(vector<int>,   trk_numberOfValidHitsETL)       \
    DATA(vector<int>,   trk_isHighPurity)               \
    DATA(vector<int>,   trk_hasMTD)                     \
    DATA(vector<float>, trk_genPdgId)                   \
    DATA(vector<float>, trk_genPt)                      \
    DATA(vector<float>, trk_genEta)                     \
    DATA(vector<float>, trk_genPhi)                     \
    DATA(vector<float>, trk_genDR)                      \
    DATA(vector<float>, trk_genVtx_x)                   \
    DATA(vector<float>, trk_genVtx_y)                   \
    DATA(vector<float>, trk_genVtx_z)                   \
    DATA(vector<float>, trk_genVtx_t)          


#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
