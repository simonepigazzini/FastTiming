#ifndef _MTD_TRACK_PU_ID_TREE_
#define _MTD_TRACK_PU_ID_TREE_

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME MTDTrackPUIDTree

#define DATA_TABLE                              \
    DATA(float, idx)                            \
    DATA(float, pt)                             \
    DATA(float, eta)                            \
    DATA(float, phi)                            \
    DATA(float, x)                              \
    DATA(float, y)                              \
    DATA(float, z)                              \
    DATA(float, dz)                             \
    DATA(float, dxy)                            \
    DATA(float, dzErr)                          \
    DATA(float, dxyErr)                         \
    DATA(float, chi2)                           \
    DATA(float, ndof)                           \
    DATA(float, t0)                             \
    DATA(float, sigmat0)                        \
    DATA(float, mtdt)                           \
    DATA(float, path_len)                       \
    DATA(float, probPi)                         \
    DATA(float, probP)                          \
    DATA(float, probK)                          \
    DATA(float, btlMatchChi2)                   \
    DATA(float, btlMatchTimeChi2)               \
    DATA(float, etlMatchChi2)                   \
    DATA(float, etlMatchTimeChi2)               \
    DATA(float, normalizedChi2)                 \
    DATA(float, numberOfValidHits)              \
    DATA(float, numberOfLostHits)               \
    DATA(float, numberOfValidPixelBarrelHits)   \
    DATA(float, numberOfValidPixelEndcapHits)   \
    DATA(float, numberOfValidHitsBTL)           \
    DATA(float, numberOfValidHitsETL)           \
    DATA(float, hasMTD)                         \
    DATA(float, genPdgId)                       \
    DATA(float, genPt)                          \
    DATA(float, genEta)                         \
    DATA(float, genPhi)                         \
    DATA(float, genDR)                          \
    DATA(bool, simIsFromPV)                     \
    DATA(float, simPt)                          \
    DATA(float, simEta)                         \
    DATA(float, simPhi)                         \
    DATA(float, simZ)                           \
    DATA(float, genVtx_x)                       \
    DATA(float, genVtx_y)                       \
    DATA(float, genVtx_z)                       \
    DATA(float, genVtx_t)                       \
    DATA(float, pv_valid)                       \
    DATA(float, pv_ntrks)                       \
    DATA(float, pv_chi2)                        \
    DATA(float, pv_x)                           \
    DATA(float, pv_y)                           \
    DATA(float, pv_z)                           \
    DATA(float, pv_t)                           \
    DATA(float, puid)

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
