#ifndef _MTD_NEUTRALS__TREE_
#define _MTD_NEUTRALS__TREE_

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME MTDNeutralsTree

#define DATA_TABLE                              \
    DATA(int, run)                              \
    DATA(int, lumi)                             \
    DATA(int, event)                            \
    DATA(float, genpv_t)                        \
    DATA(float, genpv_z)                        \
    DATA(float, pt)                             \
    DATA(float, eta)                            \
    DATA(float, phi)                            \
    DATA(float, ecalEnergy)                     \
        DATA(float, hcalEnergy)                 \
    DATA(float, particleId)                     \
    DATA(int, clus_n)                           \
    DATA(float, gen_pt)                         \
    DATA(float, gen_eta)                        \
    DATA(float, gen_phi)                        \
    DATA(float, gen_DR)                         \
    DATA(int, gen_pdgId)                        \
    DATA(float, minDR)                          \
    DATA(float, minDEta)                        \
    DATA(float, minDPhi)                        \
    DATA(float, tof)                            \
    DATA(float, mtdTime)                        \
    DATA(float, mtdEnergy)                      \
    DATA(float, simh_DR)                        \
    DATA(float, simh_recoh_DR)                  \
    DATA(float, simh_recoh_DPhi)                \
    DATA(float, simh_energy)                    \
    DATA(float, simh_time)                      \
    DATA(float, simh_tof)

#define DATA_CLASS_TABLE                        \
    DATA(vector<float>, clus_det)               \
    DATA(vector<float>, clus_size)              \
    DATA(vector<float>, clus_size_x)            \
    DATA(vector<float>, clus_size_y)            \
    DATA(vector<float>, clus_energy)            \
    DATA(vector<float>, clus_time)              \
    DATA(vector<float>, clus_rr)                \
    DATA(vector<float>, clus_module)            \
    DATA(vector<float>, clus_modType)           \
    DATA(vector<float>, clus_eta)               \
    DATA(vector<float>, clus_phi)               \
    DATA(vector<float>, clus_seed_energy)       \
    DATA(vector<float>, clus_seed_time)         \
    DATA(vector<float>, clus_seed_x)            \
    DATA(vector<float>, clus_seed_y)            \
    DATA(vector<float>, clus_local_x)           \
    DATA(vector<float>, clus_local_y)           \
    DATA(vector<float>, clus_local_z)           \
    DATA(vector<float>, clus_global_R)          \
    DATA(vector<float>, clus_global_dist)       \
    DATA(vector<float>, clus_neu_DPhi)          \
    DATA(vector<float>, clus_neu_DEta)          \
    DATA(vector<float>, clus_neu_DR)
    
#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
