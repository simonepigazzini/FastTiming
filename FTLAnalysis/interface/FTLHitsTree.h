#ifndef FTL_HITS_TREE
#define FTL_HITS_TREE

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

//---Define the TTree branches
#define DYNAMIC_TREE_NAME FTLHitsTree

#define DATA_TABLE                       \
  DATA(int, event)                       \
  DATA(int, lumi)			    \
  DATA(int, run)			    \
  DATA(int, BTLsimHits_n)                   \
  DATA(int, BTLrecHits_n)                   \
  DATA(int, BTLclusters_n)                   \
  DATA(int, ETLsimHits_n)                   \
  DATA(int, ETLrecHits_n)                   \
  DATA(int, ETLclusters_n)                   

#define DATA_CLASS_TABLE                                                \
  DATA(vector<int>,   track_idx)                                        \
  DATA(vector<float>, track_pt)                                         \
  DATA(vector<float>, track_eta)                                        \
  DATA(vector<float>, track_phi)                                        \
  DATA(vector<float>, track_eta_atBTL)                                  \
  DATA(vector<float>, track_phi_atBTL)                                  \
  DATA(vector<float>, track_eta_atETL)                                  \
  DATA(vector<float>, track_phi_atETL)                                  \
  DATA(vector<float>, track_x)                                          \
  DATA(vector<float>, track_y)                                          \
  DATA(vector<float>, track_z)                                          \
  DATA(vector<float>, track_t)                                          \
  DATA(vector<float>, track_energy)                                     \
  DATA(vector<float>, track_normalizedChi2)                             \
  DATA(vector<int>,   track_numberOfValidHits)                          \
  DATA(vector<int>,   track_numberOfLostHits)                           \
  DATA(vector<int>,   track_isHighPurity)                               \
  DATA(vector<int>,   track_hasMTD)                                     \
  DATA(vector<float>, track_mcMatch_genPdgId)                           \
  DATA(vector<float>, track_mcMatch_genPt)                              \
  DATA(vector<float>, track_mcMatch_genEta)                             \
  DATA(vector<float>, track_mcMatch_genPhi)                             \
  DATA(vector<float>, track_mcMatch_genVtx_x)                           \
  DATA(vector<float>, track_mcMatch_genVtx_y)                           \
  DATA(vector<float>, track_mcMatch_genVtx_z)                           \
  DATA(vector<float>, track_mcMatch_genVtx_t)                           \
  DATA(vector<float>, track_mcMatch_DR)                                 \
  DATA(vector<float>, BTLsimHits_energy)                                   \
  DATA(vector<float>, BTLsimHits_time)                                     \
  DATA(vector<int>,   BTLsimHits_rr)                                       \
  DATA(vector<int>,   BTLsimHits_module)                                   \
  DATA(vector<int>,   BTLsimHits_modType)                                  \
  DATA(vector<int>,   BTLsimHits_crystal)                                  \
  DATA(vector<int>,   BTLsimHits_ieta)                                     \
  DATA(vector<int>,   BTLsimHits_iphi)                                     \
  DATA(vector<float>, BTLsimHits_entry_local_x)                            \
  DATA(vector<float>, BTLsimHits_entry_local_y)                            \
  DATA(vector<float>, BTLsimHits_entry_local_z)                            \
  DATA(vector<float>, BTLsimHits_entry_global_R)                           \
  DATA(vector<float>, BTLsimHits_exit_local_x)                             \
  DATA(vector<float>, BTLsimHits_exit_local_y)                             \
  DATA(vector<float>, BTLsimHits_exit_local_z)                             \
  DATA(vector<float>, BTLsimHits_exit_global_R)                            \
  DATA(vector<float>, BTLrecHits_energy)                                   \
  DATA(vector<float>, BTLrecHits_time)                                     \
  DATA(vector<int>,   BTLrecHits_rr)                                       \
  DATA(vector<int>,   BTLrecHits_module)                                   \
  DATA(vector<int>,   BTLrecHits_modType)                                  \
  DATA(vector<int>,   BTLrecHits_crystal)                                  \
  DATA(vector<int>,   BTLrecHits_ieta)                                     \
  DATA(vector<int>,   BTLrecHits_iphi)                                     \
  DATA(vector<float>, BTLrecHits_local_x)                                  \
  DATA(vector<float>, BTLrecHits_local_y)                                  \
  DATA(vector<float>, BTLrecHits_local_z)                                  \
  DATA(vector<float>, BTLrecHits_global_R)                                 \
  DATA(vector<float>, BTLclusters_energy)                                   \
  DATA(vector<float>, BTLclusters_time)                                     \
  DATA(vector<float>, BTLclusters_x)                                   \
  DATA(vector<float>, BTLclusters_y)                                     \
  DATA(vector<float>, BTLclusters_seed_energy)                                   \
  DATA(vector<float>, BTLclusters_seed_time)                                     \
  DATA(vector<int>,   BTLclusters_seed_x)                                   \
  DATA(vector<int>,   BTLclusters_seed_y)                                     \
  DATA(vector<int>,   BTLclusters_size)                                       \
  DATA(vector<int>,   BTLclusters_size_x)                                       \
  DATA(vector<int>,   BTLclusters_size_y)                                       \
  DATA(vector<int>,   BTLclusters_rr)                                       \
  DATA(vector<int>,   BTLclusters_module)                                   \
  DATA(vector<int>,   BTLclusters_modType)                                  \
  DATA(vector<int>,   BTLclusters_crystal)                                  \
  DATA(vector<int>,   BTLclusters_ieta)                                     \
  DATA(vector<int>,   BTLclusters_iphi)                                     \
  DATA(vector<float>, BTLclusters_local_x)                                  \
  DATA(vector<float>, BTLclusters_local_y)                                  \
  DATA(vector<float>, BTLclusters_local_z)                                  \
  DATA(vector<float>, BTLclusters_global_R)                                 \
  DATA(vector<int>,   BTLmatchedSimHits_n)                                 \
  DATA(vector<int>,   BTLmatchedRecHits_n)                                 \
  DATA(vector<int>,   BTLmatchedClusters_n)                                 \
  DATA(vector, BTLmatchedSimHits_idx,           <vector<int> >)            \
  DATA(vector, BTLmatchedSimHits_energy,        <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_energyCorr,    <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_time,          <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_rr,            <vector<int> >)            \
  DATA(vector, BTLmatchedSimHits_module,        <vector<int> >)            \
  DATA(vector, BTLmatchedSimHits_modType,       <vector<int> >)            \
  DATA(vector, BTLmatchedSimHits_crystal,       <vector<int> >)            \
  DATA(vector, BTLmatchedSimHits_ieta,          <vector<int> >)            \
  DATA(vector, BTLmatchedSimHits_iphi,          <vector<int> >)            \
  DATA(vector, BTLmatchedSimHits_entry_local_x, <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_entry_local_y, <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_entry_local_z, <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_entry_global_R,<vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_exit_local_x,  <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_exit_local_y,  <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_exit_local_z,  <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_exit_global_R, <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_track_Deta,    <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_track_Dphi,    <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_track_DR,      <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_track_Dz,      <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_track_RDphi,   <vector<float> >)          \
  DATA(vector, BTLmatchedSimHits_track_dist,    <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_idx,           <vector<int> >)            \
  DATA(vector, BTLmatchedRecHits_energy,        <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_energyCorr,    <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_time,          <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_rr,            <vector<int> >)            \
  DATA(vector, BTLmatchedRecHits_module,        <vector<int> >)            \
  DATA(vector, BTLmatchedRecHits_modType,       <vector<int> >)            \
  DATA(vector, BTLmatchedRecHits_crystal,       <vector<int> >)            \
  DATA(vector, BTLmatchedRecHits_ieta,          <vector<int> >)            \
  DATA(vector, BTLmatchedRecHits_iphi,          <vector<int> >)            \
  DATA(vector, BTLmatchedRecHits_local_x,       <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_local_y,       <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_local_z,       <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_global_R,      <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_track_Deta,    <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_track_Dphi,    <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_track_DR,      <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_track_Dz,      <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_track_RDphi,   <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_track_dist,    <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_sietaieta,     <vector<float> >)          \
  DATA(vector, BTLmatchedRecHits_siphiiphi,     <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_idx,           <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_energy,        <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_energyCorr,    <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_time,          <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_rr,            <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_module,        <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_modType,       <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_crystal,       <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_ieta,          <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_iphi,          <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_size,          <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_size_x,        <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_size_y,        <vector<int> >)            \
  DATA(vector, BTLmatchedClusters_local_x,       <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_local_y,       <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_local_z,       <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_global_R,      <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_track_Deta,    <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_track_Dphi,    <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_track_DR,      <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_track_Dz,      <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_track_RDphi,   <vector<float> >)          \
  DATA(vector, BTLmatchedClusters_track_dist,    <vector<float> >)          \
  DATA(vector<float>, ETLsimHits_energy)                                   \
  DATA(vector<float>, ETLsimHits_time)                                     \
  DATA(vector<int>,   ETLsimHits_rr)                                       \
  DATA(vector<int>,   ETLsimHits_module)                                   \
  DATA(vector<int>,   ETLsimHits_modType)                                  \
  DATA(vector<int>,   ETLsimHits_crystal)                                  \
  DATA(vector<int>,   ETLsimHits_ieta)                                     \
  DATA(vector<int>,   ETLsimHits_iphi)                                     \
  DATA(vector<float>, ETLsimHits_entry_local_x)                            \
  DATA(vector<float>, ETLsimHits_entry_local_y)                            \
  DATA(vector<float>, ETLsimHits_entry_local_z)                            \
  DATA(vector<float>, ETLsimHits_entry_global_R)                           \
  DATA(vector<float>, ETLsimHits_exit_local_x)                             \
  DATA(vector<float>, ETLsimHits_exit_local_y)                             \
  DATA(vector<float>, ETLsimHits_exit_local_z)                             \
  DATA(vector<float>, ETLsimHits_exit_global_R)                            \
  DATA(vector<float>, ETLrecHits_energy)                                   \
  DATA(vector<float>, ETLrecHits_time)                                     \
  DATA(vector<int>,   ETLrecHits_rr)                                       \
  DATA(vector<int>,   ETLrecHits_module)                                   \
  DATA(vector<int>,   ETLrecHits_modType)                                  \
  DATA(vector<int>,   ETLrecHits_crystal)                                  \
  DATA(vector<int>,   ETLrecHits_ieta)                                     \
  DATA(vector<int>,   ETLrecHits_iphi)                                     \
  DATA(vector<float>, ETLrecHits_local_x)                                  \
  DATA(vector<float>, ETLrecHits_local_y)                                  \
  DATA(vector<float>, ETLrecHits_local_z)                                  \
  DATA(vector<float>, ETLrecHits_global_R)                                 \
  DATA(vector<float>, ETLclusters_energy)                                   \
  DATA(vector<float>, ETLclusters_time)                                     \
  DATA(vector<float>, ETLclusters_x)                                   \
  DATA(vector<float>, ETLclusters_y)                                     \
  DATA(vector<float>, ETLclusters_seed_energy)                                   \
  DATA(vector<float>, ETLclusters_seed_time)                                     \
  DATA(vector<int>,   ETLclusters_seed_x)                                   \
  DATA(vector<int>,   ETLclusters_seed_y)                                     \
  DATA(vector<int>,   ETLclusters_size)                                       \
  DATA(vector<int>,   ETLclusters_size_x)                                       \
  DATA(vector<int>,   ETLclusters_size_y)                                       \
  DATA(vector<int>,   ETLclusters_rr)                                       \
  DATA(vector<int>,   ETLclusters_module)                                   \
  DATA(vector<int>,   ETLclusters_modType)                                  \
  DATA(vector<int>,   ETLclusters_crystal)                                  \
  DATA(vector<int>,   ETLclusters_ieta)                                     \
  DATA(vector<int>,   ETLclusters_iphi)                                     \
  DATA(vector<float>, ETLclusters_local_x)                                  \
  DATA(vector<float>, ETLclusters_local_y)                                  \
  DATA(vector<float>, ETLclusters_local_z)                                  \
  DATA(vector<float>, ETLclusters_global_R)                                 \
  DATA(vector<int>,   ETLmatchedSimHits_n)                                 \
  DATA(vector<int>,   ETLmatchedRecHits_n)                                 \
  DATA(vector<int>,   ETLmatchedClusters_n)                                 \
  DATA(vector, ETLmatchedSimHits_idx,           <vector<int> >)            \
  DATA(vector, ETLmatchedSimHits_energy,        <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_energyCorr,    <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_time,          <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_rr,            <vector<int> >)            \
  DATA(vector, ETLmatchedSimHits_module,        <vector<int> >)            \
  DATA(vector, ETLmatchedSimHits_modType,       <vector<int> >)            \
  DATA(vector, ETLmatchedSimHits_crystal,       <vector<int> >)            \
  DATA(vector, ETLmatchedSimHits_ieta,          <vector<int> >)            \
  DATA(vector, ETLmatchedSimHits_iphi,          <vector<int> >)            \
  DATA(vector, ETLmatchedSimHits_entry_local_x, <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_entry_local_y, <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_entry_local_z, <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_entry_global_R,<vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_exit_local_x,  <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_exit_local_y,  <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_exit_local_z,  <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_exit_global_R, <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_track_Deta,    <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_track_Dphi,    <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_track_DR,      <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_track_Dz,      <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_track_RDphi,   <vector<float> >)          \
  DATA(vector, ETLmatchedSimHits_track_dist,    <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_idx,           <vector<int> >)            \
  DATA(vector, ETLmatchedRecHits_energy,        <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_energyCorr,    <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_time,          <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_rr,            <vector<int> >)            \
  DATA(vector, ETLmatchedRecHits_module,        <vector<int> >)            \
  DATA(vector, ETLmatchedRecHits_modType,       <vector<int> >)            \
  DATA(vector, ETLmatchedRecHits_crystal,       <vector<int> >)            \
  DATA(vector, ETLmatchedRecHits_ieta,          <vector<int> >)            \
  DATA(vector, ETLmatchedRecHits_iphi,          <vector<int> >)            \
  DATA(vector, ETLmatchedRecHits_local_x,       <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_local_y,       <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_local_z,       <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_global_R,      <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_track_Deta,    <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_track_Dphi,    <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_track_DR,      <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_track_Dz,      <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_track_RDphi,   <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_track_dist,    <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_sietaieta,     <vector<float> >)          \
  DATA(vector, ETLmatchedRecHits_siphiiphi,     <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_idx,           <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_energy,        <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_energyCorr,    <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_time,          <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_rr,            <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_module,        <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_modType,       <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_crystal,       <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_ieta,          <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_iphi,          <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_size,          <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_size_x,        <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_size_y,        <vector<int> >)            \
  DATA(vector, ETLmatchedClusters_local_x,       <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_local_y,       <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_local_z,       <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_global_R,      <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_track_Deta,    <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_track_Dphi,    <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_track_DR,      <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_track_Dz,      <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_track_RDphi,   <vector<float> >)          \
  DATA(vector, ETLmatchedClusters_track_dist,    <vector<float> >)          

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
