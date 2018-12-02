import ROOT as R
import math as M

f=R.TFile("/tmp/dumpHits_mb.root")

dh=f.Get("DumpHits")

histos = {}
histos["track_pt"]=R.TH1F("track_pt","track_pt",100,0.,10.)
histos["track_eta"]=R.TH1F("track_eta","track_eta",100,0.,1.5)
histos["track_phi"]=R.TH1F("track_phi","track_phi",100,-M.pi,M.pi)
histos["matchedTrack_pt"]=R.TH1F("matchedTrack_pt","matchedTrack_pt",100,0.,10.)
histos["matchedTrack_eta"]=R.TH1F("matchedTrack_eta","matchedTrack_eta",100,0.,1.5)
histos["matchedTrack_phi"]=R.TH1F("matchedTrack_phi","matchedTrack_phi",100,-M.pi,M.pi)
histos["matchedTrack_nCluster"]=R.TH1F("matchedTrack_nCluster","matchedTrack_nCluster",10,0.,10.)
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
histos["matchedCluster_sizeX"]=R.TH1F("matchedCluster_sizeX","matchedCluster_sizeX",20,0.,20.)
histos["matchedCluster_sizeY"]=R.TH1F("matchedCluster_sizeY","matchedCluster_sizeY",20,0.,20.)
histos["matchedCluster_size_vs_pt"]=R.TH2F("matchedCluster_size_vs_pt","matchedCluster_size_vs_pt",100,0.,10.,20,-0.5,19.5)
histos["matchedCluster_size_vs_eta"]=R.TH2F("matchedCluster_size_vs_eta","matchedCluster_size_vs_eta",100,0.,1.5,20,-0.5,19.5)

for event in dh:
    for iclus in range(0,event.clusters_n):
        histos["cluster_energy"].Fill(event.clusters_energy[iclus])
        histos["cluster_time"].Fill(event.clusters_time[iclus])
        histos["cluster_size"].Fill(event.clusters_size[iclus])
        histos["cluster_sizeX"].Fill(event.clusters_size_x[iclus])
        histos["cluster_sizeY"].Fill(event.clusters_size_y[iclus])
        histos["cluster_seedEnergyRatio"].Fill(event.clusters_seed_energy[iclus]/event.clusters_energy[iclus])
    for itrack in range(0,len(event.track_idx)):
        histos["track_pt"].Fill(event.track_pt[itrack])
        histos["track_eta"].Fill(abs(event.track_eta[itrack]))
        histos["track_phi"].Fill(event.track_phi[itrack])
        histos["matchedTrack_nCluster"].Fill(event.matchedClusters_n[itrack])
        if (event.matchedClusters_n[itrack]>0):
            histos["matchedTrack_pt"].Fill(event.track_pt[itrack])
            histos["matchedTrack_eta"].Fill(abs(event.track_eta[itrack]))
            histos["matchedTrack_phi"].Fill(event.track_phi[itrack])
        for iclus in range(0,event.matchedClusters_n[itrack]):
            histos["matchedCluster_energy"].Fill(event.matchedClusters_energy[itrack][iclus])
            histos["matchedCluster_time"].Fill(event.matchedClusters_time[itrack][iclus])
            histos["matchedCluster_DR"].Fill(event.matchedClusters_track_DR[itrack][iclus])
            histos["matchedCluster_size"].Fill(event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_size_vs_pt"].Fill(event.track_pt[itrack],event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_size_vs_eta"].Fill(event.track_eta[itrack],event.matchedClusters_size[itrack][iclus])
            histos["matchedCluster_sizeX"].Fill(event.matchedClusters_size_x[itrack][iclus])
            histos["matchedCluster_sizeY"].Fill(event.matchedClusters_size_y[itrack][iclus])

histos["effCluster_pt"]=R.TGraphAsymmErrors(histos["matchedTrack_pt"],histos["track_pt"])
histos["effCluster_eta"]=R.TGraphAsymmErrors(histos["matchedTrack_eta"],histos["track_eta"])
histos["effCluster_phi"]=R.TGraphAsymmErrors(histos["matchedTrack_phi"],histos["track_phi"])

fOut=R.TFile("clusterPlots.root","RECREATE")
for hn, histo in histos.iteritems():
    if isinstance(histo,R.TH1F):
        histo.SetMinimum(0.)
    histo.Write()
fOut.Close()
