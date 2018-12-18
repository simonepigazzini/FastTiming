import ROOT as R
import math as M
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input',dest='input')
parser.add_argument('--output',dest='output')
parser.add_argument('--layout',dest='layout')
args = parser.parse_args()

f=R.TFile(args.input)

dh=f.Get("DumpHits")

histos = {}

# time offset introduced in mtd2 due to light collection time (pitch.z()/2 * 7.5ps for bars)
if args.layout == "barzflat":
    t_offset=0.210
elif args.layout == "barphi":
    t_offset=0.160
elif args.layout == "tile":
    t_offset=0.200

histos["track_pt"]=R.TH1F("track_pt","track_pt",100,0.,10.)
histos["track_eta"]=R.TH1F("track_eta","track_eta",100,0.,1.5)
histos["track_phi"]=R.TH1F("track_phi","track_phi",100,-M.pi,M.pi)

histos["matchedTrack_pt"]=R.TH1F("matchedTrack_pt","matchedTrack_pt",100,0.,10.)
histos["matchedTrack_eta"]=R.TH1F("matchedTrack_eta","matchedTrack_eta",100,0.,1.5)
histos["matchedTrack_phi"]=R.TH1F("matchedTrack_phi","matchedTrack_phi",100,-M.pi,M.pi)
histos["matchedTrack_nCluster"]=R.TH1F("matchedTrack_nCluster","matchedTrack_nCluster",10,0.,10.)
histos["matchedTrack_nHits"]=R.TH1F("matchedTrack_nHits","matchedTrack_nHits",10,0.,10.)

histos["mtdTrack_pt"]=R.TH1F("mtdTrack_pt","mtdTrack_pt",100,0.,10.)
histos["mtdTrack_eta"]=R.TH1F("mtdTrack_eta","mtdTrack_eta",100,0.,1.5)
histos["mtdTrack_phi"]=R.TH1F("mtdTrack_phi","mtdTrack_phi",100,-M.pi,M.pi)
histos["mtdTrack_dt"]=R.TH1F("mtdTrack_dt","mtdTrack_dt",100,-0.2,0.2)
histos["mtdTrack_dt_vs_pt"]=R.TH2F("mtdTrack_dt_vs_pt","mtdTrack_dt_vs_pt",20,0,4,100,-0.2,0.2)
histos["mtdTrack_dt_vs_eta"]=R.TH2F("mtdTrack_dt_vs_eta","mtdTrack_dt_vs_eta",40,0,1.5,100,-0.2,0.2)
histos["mtdTrack_dz"]=R.TH1F("mtdTrack_dz","mtdTrack_dz",200,-0.05,0.05)
histos["mtdTrack_dz_vs_pt"]=R.TH2F("mtdTrack_dz_vs_pt","mtdTrack_dz_vs_pt",20,0,4,200,-0.05,0.05)
histos["mtdTrack_dz_vs_eta"]=R.TH2F("mtdTrack_dz_vs_eta","mtdTrack_dz_vs_eta",40,0,1.5,200,-0.05,0.05)
histos["mtdTrack_ptRes"]=R.TH1F("mtdTrack_ptRes","mtdTrack_ptRes",200,-0.05,0.05)
histos["mtdTrack_ptRes_vs_pt"]=R.TH2F("mtdTrack_ptRes_vs_pt","mtdTrack_ptRes_vs_pt",20,0,4,200,-0.05,0.05)
histos["mtdTrack_ptRes_vs_eta"]=R.TH2F("mtdTrack_ptRes_vs_eta","mtdTrack_ptRes_vs_eta",40,0,1.5,200,-0.05,0.05)

histos["cluster_energy"]=R.TH1F("cluster_energy","cluster_energy",100,0.,20.)
histos["cluster_time"]=R.TH1F("cluster_time","cluster_time",100,0.,25.)
histos["cluster_size"]=R.TH1F("cluster_size","cluster_size",20,0.,20.)
histos["cluster_sizeX"]=R.TH1F("cluster_sizeX","cluster_sizeX",20,0.,20.)
histos["cluster_sizeY"]=R.TH1F("cluster_sizeY","cluster_sizeY",20,0.,20.)
histos["cluster_seedEnergyRatio"]=R.TH1F("cluster_seedEnergyRatio","cluster_seedEnergyRatio",110,0.,1.1)

histos["matchedCluster_energy"]=R.TH1F("matchedCluster_energy","matchedCluster_energy",100,0.,20.)
histos["matchedCluster_time"]=R.TH1F("matchedCluster_time","matchedCluster_time",100,0.,25.)
histos["matchedCluster_DR"]=R.TH1F("matchedCluster_DR","matchedCluster_DR",100,0.,0.05)
histos["matchedCluster_size"]=R.TH1F("matchedCluster_size","matchedCluster_size",20,0.,20.)
histos["matchedCluster_size_vs_pt"]=R.TH2F("matchedCluster_size_vs_pt","matchedCluster_size_vs_pt",100,0.,4.,20,-0.5,19.5)
histos["matchedCluster_size_vs_eta"]=R.TH2F("matchedCluster_size_vs_eta","matchedCluster_size_vs_eta",100,0.,1.5,20,-0.5,19.5)
histos["matchedCluster_sizeX"]=R.TH1F("matchedCluster_sizeX","matchedCluster_sizeX",20,0.,20.)
histos["matchedCluster_sizeY"]=R.TH1F("matchedCluster_sizeY","matchedCluster_sizeY",20,0.,20.)

ievent=0
for event in dh:
    ievent=ievent+1
    if (ievent%100 ==0):
        print "Analysing event %d"%ievent
        
    for iclus in range(0,event.clusters_n):
        histos["cluster_energy"].Fill(event.clusters_energy[iclus])
        histos["cluster_time"].Fill(event.clusters_time[iclus])
        histos["cluster_size"].Fill(event.clusters_size[iclus])
        histos["cluster_sizeX"].Fill(event.clusters_size_x[iclus])
        histos["cluster_sizeY"].Fill(event.clusters_size_y[iclus])
        histos["cluster_seedEnergyRatio"].Fill(event.clusters_seed_energy[iclus]/event.clusters_energy[iclus])
    for itrack in range(0,len(event.track_idx)):
        histos["track_pt"].Fill(event.track_pt[itrack])
        histos["track_eta"].Fill(abs(event.track_eta_atBTL[itrack]))
        histos["track_phi"].Fill(event.track_phi_atBTL[itrack])
        histos["matchedTrack_nCluster"].Fill(event.matchedClusters_n[itrack])

        if (event.matchedClusters_n[itrack]>0):
            histos["matchedTrack_pt"].Fill(event.track_pt[itrack])
            histos["matchedTrack_eta"].Fill(abs(event.track_eta_atBTL[itrack]))
            histos["matchedTrack_phi"].Fill(event.track_phi_atBTL[itrack])

        if (event.track_hasMTD[itrack]>0):
            histos["mtdTrack_pt"].Fill(event.track_pt[itrack])
            histos["mtdTrack_eta"].Fill(abs(event.track_eta_atBTL[itrack]))
            histos["mtdTrack_phi"].Fill(event.track_phi_atBTL[itrack])
            histos["mtdTrack_dt"].Fill(event.track_t[itrack]-event.track_mcMatch_genVtx_t[itrack]-t_offset)
            histos["mtdTrack_dt_vs_pt"].Fill(event.track_pt[itrack],event.track_t[itrack]-event.track_mcMatch_genVtx_t[itrack]-t_offset)
            histos["mtdTrack_dt_vs_eta"].Fill(abs(event.track_eta_atBTL[itrack]),event.track_t[itrack]-event.track_mcMatch_genVtx_t[itrack]-t_offset)
            histos["mtdTrack_dz"].Fill(event.track_z[itrack]-event.track_mcMatch_genVtx_z[itrack])
            histos["mtdTrack_dz_vs_pt"].Fill(event.track_pt[itrack],event.track_z[itrack]-event.track_mcMatch_genVtx_z[itrack])
            histos["mtdTrack_dz_vs_eta"].Fill(abs(event.track_eta_atBTL[itrack]),event.track_z[itrack]-event.track_mcMatch_genVtx_z[itrack])
            histos["mtdTrack_ptRes"].Fill(event.track_pt[itrack]/event.track_mcMatch_genPt[itrack]-1.)
            histos["mtdTrack_ptRes_vs_pt"].Fill(event.track_pt[itrack],event.track_pt[itrack]/event.track_mcMatch_genPt[itrack]-1.)
            histos["mtdTrack_ptRes_vs_eta"].Fill(abs(event.track_eta_atBTL[itrack]),event.track_pt[itrack]/event.track_mcMatch_genPt[itrack]-1.)

        for iclus in range(0,event.matchedClusters_n[itrack]):
            histos["matchedCluster_energy"].Fill(event.matchedClusters_energy[itrack][iclus])
            histos["matchedCluster_time"].Fill(event.matchedClusters_time[itrack][iclus])
            histos["matchedCluster_DR"].Fill(event.matchedClusters_track_DR[itrack][iclus])
            histos["matchedCluster_size"].Fill(event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_size_vs_pt"].Fill(event.track_pt[itrack],event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_size_vs_eta"].Fill(event.track_eta_atBTL[itrack],event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_sizeX"].Fill(event.matchedClusters_size_x[itrack][iclus])
            histos["matchedCluster_sizeY"].Fill(event.matchedClusters_size_y[itrack][iclus])

histos["effCluster_pt"]=R.TGraphAsymmErrors(histos["matchedTrack_pt"],histos["track_pt"])
histos["effCluster_eta"]=R.TGraphAsymmErrors(histos["matchedTrack_eta"],histos["track_eta"])
histos["effCluster_phi"]=R.TGraphAsymmErrors(histos["matchedTrack_phi"],histos["track_phi"])

histos["effMtd_pt"]= R.TGraphAsymmErrors(histos["mtdTrack_pt"],histos["track_pt"])
histos["effMtd_eta"]=R.TGraphAsymmErrors(histos["mtdTrack_eta"],histos["track_eta"])
histos["effMtd_phi"]=R.TGraphAsymmErrors(histos["mtdTrack_phi"],histos["track_phi"])

fOut=R.TFile(args.output,"RECREATE")
for hn, histo in histos.iteritems():
    if isinstance(histo,R.TH1F):
        histo.SetMinimum(0.)
    if isinstance(histo,R.TGraphAsymmErrors):
        histo.SetMinimum(0.)
        histo.SetMaximum(1.1)
    histo.Write()
fOut.Close()
print "Saved histos in "+args.output
