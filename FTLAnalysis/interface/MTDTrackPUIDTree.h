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
    DATA(float, dzErr)                          \
    DATA(float, dxyErr)                         \
        DATA(float, chi2)                       \
    DATA(float, ndof)                           \
    DATA(float, t0)                             \
    DATA(float, sigmat0)                        \
    DATA(float, mtdt)                           \
    DATA(float, path_len)                       \
    DATA(float, probPi)                         \
    DATA(float, probP)                          \
    DATA(float, probK)                          \
    DATA(float, energy)                         \
    DATA(float, btlMatchChi2)                   \
    DATA(float, btlMatchTimeChi2)               \
    DATA(float, etlMatchChi2)                   \
    DATA(float, etlMatchTimeChi2)               \
    DATA(float, normalizedChi2)                 \
    DATA(int, numberOfValidHits)                \
    DATA(int, numberOfLostHits)                 \
    DATA(int, numberOfValidPixelBarrelHits)     \
    DATA(int, numberOfValidPixelEndcapHits)     \
    DATA(int, numberOfValidStripTIBHits)        \
    DATA(int, numberOfValidStripTIDHits)        \
    DATA(int, numberOfValidStripTOBHits)        \
    DATA(int, numberOfValidStripTECHits)        \
    DATA(int, numberOfValidHitsBTL)             \
    DATA(int, numberOfValidHitsETL)             \
    DATA(int, isHighPurity)                     \
    DATA(bool, hasMTD)                          \
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
    DATA(float, pv3d_valid)                     \
    DATA(float, pv3d_ntrks)                     \
    DATA(float, pv3d_chi2)                      \
    DATA(float, pv3d_x)                         \
    DATA(float, pv3d_y)                         \
    DATA(float, pv3d_z)                         \
    DATA(float, pv4d_valid)                     \
    DATA(float, pv4d_ntrks)                     \
    DATA(float, pv4d_chi2)                      \
    DATA(float, pv4d_x)                         \
    DATA(float, pv4d_y)                         \
    DATA(float, pv4d_z)                         \
    DATA(float, pv4d_t)                         

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
