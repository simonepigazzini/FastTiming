#ifndef _MTD_4D_TREE_
#define _MTD_4D_TREE_

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME MTD4DTree

#define DATA_TABLE                              \
    DATA(int, event)                            \
    DATA(int, lumi)                             \
    DATA(int, run)                              \
    DATA(int, vtx3D_n_vtxs)                     \
    DATA(int, vtx3D_n_good)                     \
    DATA(int, vtx3D_best_dz)                    \
    DATA(int, vtx4D_n_vtxs)                     \
    DATA(int, vtx4D_n_good)                     \
    DATA(int, vtx4D_best_dz)                    \
        DATA(int, vtx4D_best_dzt)               \
    DATA(int, vtx4DNoPID_n_vtxs)                \
    DATA(int, vtx4DNoPID_n_good)                \
    DATA(int, vtx4DNoPID_best_dz)               \
    DATA(int, vtx4DNoPID_best_dzt)              \
    DATA(float, gen_sumpt)                      \
    DATA(float, gen_z)                          \
    DATA(float, gen_t)                          \
    DATA(float, sim_z)                          \
    DATA(float, sim_t)                          \
    DATA(int, sim_n_vtxs)                       \
    DATA(int, sim_n_trks)                       

#define DATA_CLASS_TABLE                                \
    DATA(vector<int>, vtx3D_valid)                      \
    DATA(vector<float>, vtx3D_z)                        \
    DATA(vector<float>, vtx3D_t)                        \
    DATA(vector<float>, vtx3D_chi2)                     \
    DATA(vector<float>, vtx3D_ntrks)                    \
    DATA(vector<float>, vtx3D_ntrks_genmatch)           \
    DATA(vector<float>, vtx3D_sumpt_genmatch)           \
    DATA(vector<float>, vtx3D_sumpt)                    \
    DATA(vector<int>, vtx4D_valid)                      \
    DATA(vector<float>, vtx4D_z)                        \
    DATA(vector<float>, vtx4D_t)                        \
    DATA(vector<float>, vtx4D_chi2)                     \
    DATA(vector<float>, vtx4D_ntrks)                    \
    DATA(vector<float>, vtx4D_ntrks_genmatch)           \
    DATA(vector<float>, vtx4D_sumpt_genmatch)           \
    DATA(vector<float>, vtx4D_sumpt)                    \
    DATA(vector<int>, vtx4DNoPID_valid)                 \
    DATA(vector<float>, vtx4DNoPID_z)                   \
    DATA(vector<float>, vtx4DNoPID_t)                   \
    DATA(vector<float>, vtx4DNoPID_chi2)                \
    DATA(vector<float>, vtx4DNoPID_ntrks)               \
    DATA(vector<float>, vtx4DNoPID_ntrks_genmatch)      \
    DATA(vector<float>, vtx4DNoPID_sumpt_genmatch)      \
    DATA(vector<float>, vtx4DNoPID_sumpt)                    

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
