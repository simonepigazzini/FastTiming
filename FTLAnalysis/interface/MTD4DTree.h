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
    DATA(float, vtx4DNoPID_N)                   \
    DATA(float, vtx4DNoPID_Nnofake)             \
    DATA(float, vtx4D_0_valid)                  \
    DATA(float, vtx4D_0_dz)                     \
    DATA(float, vtx4D_0_dt)                     \
    DATA(float, vtx4D_0_chi2)                   \
    DATA(float, vtx4D_0_ntrks)                  \
    DATA(float, vtx4DNoPID_0_valid)             \
    DATA(float, vtx4DNoPID_0_dz)                \
    DATA(float, vtx4DNoPID_0_dt)                \
    DATA(float, vtx4DNoPID_0_chi2)              \
    DATA(float, vtx4DNoPID_0_ntrks)             \
    DATA(int, vtx4D_best_idx)                   \
    DATA(float, vtx4D_best_dz)                  \
    DATA(float, vtx4D_best_dt)                  \
    DATA(float, vtx4D_best_chi2)                \
    DATA(float, vtx4D_best_ntrks)               \
    DATA(int, vtx4D_best_dzt_idx)               \
    DATA(float, vtx4D_best_dzt_dz)              \
    DATA(float, vtx4D_best_dzt_dt)

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
