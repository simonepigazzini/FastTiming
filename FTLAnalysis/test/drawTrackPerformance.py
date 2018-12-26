import ROOT as R
import math as M
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input',dest='input')
parser.add_argument('--output',dest='output')
parser.add_argument('--layout',dest='layout')
parser.add_argument('--events',dest='events',default='-1')
args = parser.parse_args()

f=R.TFile(args.input)

dh=f.Get("DumpHits")

histos = {}

# time offset introduced in mtd2 due to light collection time (pitch.z()/2 * 7.5ps for bars)
t_offset={}
t_offset['ETL']=0.
if args.layout == "barzflat":
    t_offset['BTL']=0.210
elif args.layout == "barphi":
    t_offset['BTL']=0.160
elif args.layout == "tile":
    t_offset['BTL']=0.200

histos["track_pt"]=R.TH1F("track_pt","track_pt",100,0.,10.)
histos["track_eta"]=R.TH1F("track_eta","track_eta",100,0.,3,)
histos["track_phi"]=R.TH1F("track_phi","track_phi",100,-M.pi,M.pi)

histos["mtdTrack_pt"]=R.TH1F("mtdTrack_pt","mtdTrack_pt",100,0.,10.)
histos["mtdTrack_eta"]=R.TH1F("mtdTrack_eta","mtdTrack_eta",100,0.,3,)
histos["mtdTrack_phi"]=R.TH1F("mtdTrack_phi","mtdTrack_phi",100,-M.pi,M.pi)
histos["mtdTrack_dt"]=R.TH1F("mtdTrack_dt","mtdTrack_dt",100,-0.2,0.2)
histos["mtdTrack_dt_vs_pt"]=R.TH2F("mtdTrack_dt_vs_pt","mtdTrack_dt_vs_pt",50,0,10,100,-0.2,0.2)
histos["mtdTrack_dt_vs_eta"]=R.TH2F("mtdTrack_dt_vs_eta","mtdTrack_dt_vs_eta",80,0,3,100,-0.2,0.2)
histos["mtdTrack_dz"]=R.TH1F("mtdTrack_dz","mtdTrack_dz",200,-0.05,0.05)
histos["mtdTrack_dz_vs_pt"]=R.TH2F("mtdTrack_dz_vs_pt","mtdTrack_dz_vs_pt",50,0,10,200,-0.05,0.05)
histos["mtdTrack_dz_vs_eta"]=R.TH2F("mtdTrack_dz_vs_eta","mtdTrack_dz_vs_eta",80,0,3,200,-0.05,0.05)
histos["mtdTrack_ptRes"]=R.TH1F("mtdTrack_ptRes","mtdTrack_ptRes",200,-0.05,0.05)
histos["mtdTrack_ptRes_vs_pt"]=R.TH2F("mtdTrack_ptRes_vs_pt","mtdTrack_ptRes_vs_pt",50,0,10,200,-0.05,0.05)
histos["mtdTrack_ptRes_vs_eta"]=R.TH2F("mtdTrack_ptRes_vs_eta","mtdTrack_ptRes_vs_eta",80,0,3,200,-0.05,0.05)


for det in ["BTL","ETL"]:
    histos[det+"track_pt"]= R.TH1F (det+"track_pt" ,det+"track_pt",100,0.,10.)
    histos[det+"track_eta"]=R.TH1F(det+"track_eta",det+"track_eta",100,0.,3,)
    histos[det+"track_phi"]=R.TH1F(det+"track_phi",det+"track_phi",100,-M.pi,M.pi)

    histos[det+"matchedTrack_nCluster"]=R.TH1F(det+"matchedTrack_nCluster",det+"matchedTrack_nCluster",10,0.,10.)
    histos[det+"matchedTrack_nRecHits"]=R.TH1F(det+"matchedTrack_nRecHits",det+"matchedTrack_nRecHits",10,0.,10.)

    histos[det+"cluster_energy"]=R.TH1F(det+"cluster_energy",det+"cluster_energy",100,0.,20.)
    histos[det+"cluster_time"]=R.TH1F(det+"cluster_time",det+"cluster_time",100,0.,25.)
    histos[det+"cluster_size"]=R.TH1F(det+"cluster_size",det+"cluster_size",20,0.,20.)
    histos[det+"cluster_sizeX"]=R.TH1F(det+"cluster_sizeX",det+"cluster_sizeX",20,0.,20.)
    histos[det+"cluster_sizeY"]=R.TH1F(det+"cluster_sizeY",det+"cluster_sizeY",20,0.,20.)
    histos[det+"cluster_seedEnergyRatio"]=R.TH1F(det+"cluster_seedEnergyRatio",det+"cluster_seedEnergyRatio",110,0.,1.1)
    
    histos[det+"recHit_energy"]=R.TH1F(det+"recHit_energy",det+"recHit_energy",100,0.,20.)
    histos[det+"recHit_time"]=R.TH1F(det+"recHit_time",det+"recHit_time",100,0.,25.)

    histos[det+"matchedClusterTrack_pt"]=R.TH1F(det+"matchedClusterTrack_pt",det+"matchedClusterTrack_pt",100,0.,10.)
    histos[det+"matchedClusterTrack_eta"]=R.TH1F(det+"matchedClusterTrack_eta",det+"matchedClusterTrack_eta",100,0.,3,)
    histos[det+"matchedClusterTrack_phi"]=R.TH1F(det+"matchedClusterTrack_phi",det+"matchedClusterTrack_phi",100,-M.pi,M.pi)
    histos[det+"matchedRecHitTrack_pt"]=R.TH1F(det+"matchedRecHitTrack_pt",det+"matchedRecHitTrack_pt",100,0.,10.)
    histos[det+"matchedRecHitTrack_eta"]=R.TH1F(det+"matchedRecHitTrack_eta",det+"matchedRecHitTrack_eta",100,0.,3,)
    histos[det+"matchedRecHitTrack_phi"]=R.TH1F(det+"matchedRecHitTrack_phi",det+"matchedRecHitTrack_phi",100,-M.pi,M.pi)
    histos[det+"matchedBestRecHitTrack_pt"]=R.TH1F(det+"matchedBestRecHitTrack_pt",det+"matchedBestRecHitTrack_pt",100,0.,10.)
    histos[det+"matchedBestRecHitTrack_eta"]=R.TH1F(det+"matchedBestRecHitTrack_eta",det+"matchedBestRecHitTrack_eta",100,0.,3,)
    histos[det+"matchedBestRecHitTrack_phi"]=R.TH1F(det+"matchedBestRecHitTrack_phi",det+"matchedBestRecHitTrack_phi",100,-M.pi,M.pi)

    histos[det+"matchedCluster_energy"]=R.TH1F(det+"matchedCluster_energy",det+"matchedCluster_energy",100,0.,20.)
    histos[det+"matchedCluster_time"]=R.TH1F(det+"matchedCluster_time",det+"matchedCluster_time",100,0.,25.)
    histos[det+"matchedCluster_DR"]=R.TH1F(det+"matchedCluster_DR",det+"matchedCluster_DR",100,0.,0.05)
    histos[det+"matchedCluster_size"]=R.TH1F(det+"matchedCluster_size",det+"matchedCluster_size",20,0.,20.)
    histos[det+"matchedCluster_size_vs_pt"]=R.TH2F(det+"matchedCluster_size_vs_pt",det+"matchedCluster_size_vs_pt",100,0.,4.,20,-0.5,19.5)
    histos[det+"matchedCluster_size_vs_eta"]=R.TH2F(det+"matchedCluster_size_vs_eta",det+"matchedCluster_size_vs_eta",100,0.,3,20,-0.5,19.5)
    histos[det+"matchedCluster_sizeX"]=R.TH1F(det+"matchedCluster_sizeX",det+"matchedCluster_sizeX",20,0.,20.)
    histos[det+"matchedCluster_sizeY"]=R.TH1F(det+"matchedCluster_sizeY",det+"matchedCluster_sizeY",20,0.,20.)

    histos[det+"bestCluster_energy"]=R.TH1F(det+"bestCluster_energy",det+"bestCluster_energy",100,0.,20.)
    histos[det+"bestCluster_time"]=R.TH1F(det+"bestCluster_time",det+"bestCluster_time",100,0.,25.)
    histos[det+"bestCluster_DR"]=R.TH1F(det+"bestCluster_DR",det+"bestCluster_DR",100,0.,0.05)
    histos[det+"bestCluster_size"]=R.TH1F(det+"bestCluster_size",det+"bestCluster_size",20,0.,20.)
    histos[det+"bestCluster_size_vs_pt"]=R.TH2F(det+"bestCluster_size_vs_pt",det+"bestCluster_size_vs_pt",100,0.,4.,20,-0.5,19.5)
    histos[det+"bestCluster_size_vs_eta"]=R.TH2F(det+"bestCluster_size_vs_eta",det+"bestCluster_size_vs_eta",100,0.,3,20,-0.5,19.5)
    histos[det+"bestCluster_sizeX"]=R.TH1F(det+"bestCluster_sizeX",det+"bestCluster_sizeX",20,0.,20.)
    histos[det+"bestCluster_sizeY"]=R.TH1F(det+"bestCluster_sizeY",det+"bestCluster_sizeY",20,0.,20.)

    histos[det+"matchedRecHit_energy"]=R.TH1F(det+"matchedRecHit_energy",det+"matchedRecHit_energy",100,0.,20.)
    histos[det+"matchedRecHit_time"]=R.TH1F(det+"matchedRecHit_time",det+"matchedRecHit_time",100,0.,25.)
    histos[det+"matchedRecHit_DR"]=R.TH1F(det+"matchedRecHit_DR",det+"matchedRecHit_DR",100,0.,0.05)

    histos[det+"bestRecHit_energy"]=R.TH1F(det+"bestRecHit_energy",det+"bestRecHit_energy",100,0.,20.)
    histos[det+"bestRecHit_time"]=R.TH1F(det+"bestRecHit_time",det+"bestRecHit_time",100,0.,25.)
    histos[det+"bestRecHit_DR"]=R.TH1F(det+"bestRecHit_DR",det+"bestRecHit_DR",100,0.,0.05)
    histos[det+"bestRecHit_time_vs_pt"]=R.TH2F(det+"bestRecHit_time_vs_pt",det+"bestRecHit_time_vs_pt",50,0,10,100,0,25)
    histos[det+"bestRecHit_time_vs_eta"]=R.TH2F(det+"bestRecHit_time_vs_eta",det+"bestRecHit_time_vs_eta",80,0,3,100,0,25)

det_id = { 'BTL':1  , 'ETL':2 }
etaCut = { 'BTL':[0,1.5]  , 'ETL':[1.5,3] }

for ievent,event in enumerate(dh):
    if (int(args.events) != -1 and ievent>int(args.events)):
        break
    if (ievent%100 ==0):
        print "Analysing event %d"%ievent

    for det in ["BTL","ETL"]:
        for iclus in range(0,event.clusters_n):
            if ( event.clusters_det[iclus] !=  det_id[det] ):
                continue 
            histos[det+"cluster_energy"].Fill(event.clusters_energy[iclus])
            histos[det+"cluster_time"].Fill(event.clusters_time[iclus])
            histos[det+"cluster_size"].Fill(event.clusters_size[iclus])
            histos[det+"cluster_sizeX"].Fill(event.clusters_size_x[iclus])
            histos[det+"cluster_sizeY"].Fill(event.clusters_size_y[iclus])
            histos[det+"cluster_seedEnergyRatio"].Fill(event.clusters_seed_energy[iclus]/event.clusters_energy[iclus])
            
            for ihit in range(0,event.recHits_n):
                if ( event.recHits_det[ihit] !=  det_id[det] ):
                    continue 
            histos[det+"recHit_energy"].Fill(event.recHits_energy[ihit])
            histos[det+"recHit_time"].Fill(event.recHits_time[ihit])

    for itrack in range(0,len(event.track_idx)):
        histos["track_pt"].Fill(event.track_pt[itrack])
        histos["track_eta"].Fill(abs(event.track_eta[itrack]))
        histos["track_phi"].Fill(event.track_phi[itrack])
        
        if (event.track_hasMTD[itrack]>0):
            t_off = t_offset['BTL'] if event.track_eta_atBTL>-100. else t_offset['ETL']
            histos["mtdTrack_pt"].Fill(event.track_pt[itrack])
            histos["mtdTrack_eta"].Fill(abs(event.track_eta[itrack]))
            histos["mtdTrack_phi"].Fill(event.track_phi[itrack])
            histos["mtdTrack_dt"].Fill(event.track_t[itrack]-event.track_mcMatch_genVtx_t[itrack]-t_off)
            histos["mtdTrack_dt_vs_pt"].Fill(event.track_pt[itrack],event.track_t[itrack]-event.track_mcMatch_genVtx_t[itrack]-t_off)
            histos["mtdTrack_dt_vs_eta"].Fill(abs(event.track_eta[itrack]),event.track_t[itrack]-event.track_mcMatch_genVtx_t[itrack]-t_off)
            histos["mtdTrack_dz"].Fill(event.track_z[itrack]-event.track_mcMatch_genVtx_z[itrack])
            histos["mtdTrack_dz_vs_pt"].Fill(event.track_pt[itrack],event.track_z[itrack]-event.track_mcMatch_genVtx_z[itrack])
            histos["mtdTrack_dz_vs_eta"].Fill(abs(event.track_eta[itrack]),event.track_z[itrack]-event.track_mcMatch_genVtx_z[itrack])
            histos["mtdTrack_ptRes"].Fill(event.track_pt[itrack]/event.track_mcMatch_genPt[itrack]-1.)
            histos["mtdTrack_ptRes_vs_pt"].Fill(event.track_pt[itrack],event.track_pt[itrack]/event.track_mcMatch_genPt[itrack]-1.)
            histos["mtdTrack_ptRes_vs_eta"].Fill(abs(event.track_eta[itrack]),event.track_pt[itrack]/event.track_mcMatch_genPt[itrack]-1.)

        for det in ["BTL","ETL"]:
            if ( 
                abs(event.track_eta[itrack]) < etaCut[det][0] or
                abs(event.track_eta[itrack]) > etaCut[det][1]
                ):
                continue
            
            histos[det+"track_pt"].Fill(event.track_pt[itrack])
            histos[det+"track_eta"].Fill(abs(event.track_eta[itrack]))
            histos[det+"track_phi"].Fill(event.track_phi[itrack])

            goodClusters=0
            bestClus=-1
            bestDR=9999
            for iclus in range(0,event.matchedClusters_n[itrack]):
                if (event.matchedClusters_time[itrack][iclus]>20):
                    continue
                if (event.matchedClusters_det[itrack][iclus]!=det_id[det]):
                    continue
                goodClusters=goodClusters+1
                if (event.matchedClusters_track_DR[itrack][iclus]<bestDR):
                    bestDR=event.matchedClusters_track_DR[itrack][iclus]
                    bestClus=iclus
                histos[det+"matchedCluster_energy"].Fill(event.matchedClusters_energy[itrack][iclus])
                histos[det+"matchedCluster_time"].Fill(event.matchedClusters_time[itrack][iclus])
                histos[det+"matchedCluster_DR"].Fill(event.matchedClusters_track_DR[itrack][iclus])
                histos[det+"matchedCluster_size"].Fill(event.matchedClusters_size[itrack][iclus])
                histos[det+"matchedCluster_size_vs_pt"].Fill(event.track_pt[itrack],event.matchedClusters_size[itrack][iclus])
                histos[det+"matchedCluster_size_vs_eta"].Fill(event.track_eta[itrack],event.matchedClusters_size[itrack][iclus])
                histos[det+"matchedCluster_sizeX"].Fill(event.matchedClusters_size_x[itrack][iclus])
                histos[det+"matchedCluster_sizeY"].Fill(event.matchedClusters_size_y[itrack][iclus])

            if (goodClusters>0):
                histos[det+"matchedClusterTrack_pt"].Fill(event.track_pt[itrack])
                histos[det+"matchedClusterTrack_eta"].Fill(abs(event.track_eta[itrack]))
                histos[det+"matchedClusterTrack_phi"].Fill(event.track_phi[itrack])

            if (bestClus>=0):
                histos[det+"bestCluster_energy"].Fill(event.matchedClusters_energy[itrack][bestClus])
                histos[det+"bestCluster_time"].Fill(event.matchedClusters_time[itrack][bestClus])
                histos[det+"bestCluster_DR"].Fill(event.matchedClusters_track_DR[itrack][bestClus])
                histos[det+"bestCluster_size"].Fill(event.matchedClusters_size[itrack][bestClus])
                histos[det+"bestCluster_size_vs_pt"].Fill(event.track_pt[itrack],event.matchedClusters_size[itrack][bestClus])
                histos[det+"bestCluster_size_vs_eta"].Fill(event.track_eta[itrack],event.matchedClusters_size[itrack][bestClus])
                histos[det+"bestCluster_sizeX"].Fill(event.matchedClusters_size_x[itrack][bestClus])
                histos[det+"bestCluster_sizeY"].Fill(event.matchedClusters_size_y[itrack][bestClus])
            
            goodRecHits=0
            bestHit=-1
            bestDR=9999
            for ihit in range(0,event.matchedRecHits_n[itrack]):
                if (event.matchedRecHits_time[itrack][ihit]>20):
                    continue
                if (event.matchedRecHits_det[itrack][ihit]!=det_id[det]):
                    continue
                goodRecHits=goodRecHits+1
                if (event.matchedRecHits_track_DR[itrack][ihit]<bestDR):
                    bestDR=event.matchedRecHits_track_DR[itrack][ihit]
                    bestHit=ihit
                histos[det+"matchedRecHit_energy"].Fill(event.matchedRecHits_energy[itrack][ihit])
                histos[det+"matchedRecHit_time"].Fill(event.matchedRecHits_time[itrack][ihit])
                histos[det+"matchedRecHit_DR"].Fill(event.matchedRecHits_track_DR[itrack][ihit])
                
            if (goodRecHits>0):
                histos[det+"matchedRecHitTrack_pt"].Fill(event.track_pt[itrack])
                histos[det+"matchedRecHitTrack_eta"].Fill(abs(event.track_eta[itrack]))
                histos[det+"matchedRecHitTrack_phi"].Fill(event.track_phi[itrack])
                    
            if (bestHit>=0):
                histos[det+"bestRecHit_energy"].Fill(event.matchedRecHits_energy[itrack][bestHit])
                histos[det+"bestRecHit_time"].Fill(event.matchedRecHits_time[itrack][bestHit])
                histos[det+"bestRecHit_DR"].Fill(event.matchedRecHits_track_DR[itrack][bestHit])
                histos[det+"bestRecHit_time_vs_pt"].Fill(event.track_pt[itrack],event.matchedRecHits_time[itrack][bestHit])
                histos[det+"bestRecHit_time_vs_eta"].Fill(abs(event.track_eta[itrack]),event.matchedRecHits_time[itrack][bestHit])
                histos[det+"matchedBestRecHitTrack_pt"].Fill(event.track_pt[itrack])
                histos[det+"matchedBestRecHitTrack_eta"].Fill(abs(event.track_eta[itrack]))
                histos[det+"matchedBestRecHitTrack_phi"].Fill(event.track_phi[itrack])
                    
            histos[det+"matchedTrack_nCluster"].Fill(goodClusters)
            histos[det+"matchedTrack_nRecHits"].Fill(goodRecHits)

for det in ["BTL","ETL"]:
    histos[det+"effCluster_pt"]=R.TGraphAsymmErrors( histos[det+"matchedClusterTrack_pt"], histos[det+"track_pt"])
    histos[det+"effCluster_eta"]=R.TGraphAsymmErrors(histos[det+"matchedClusterTrack_eta"],histos[det+"track_eta"])
    histos[det+"effCluster_phi"]=R.TGraphAsymmErrors(histos[det+"matchedClusterTrack_phi"],histos[det+"track_phi"])

    histos[det+"effRecHit_pt"]=R.TGraphAsymmErrors( histos[det+"matchedRecHitTrack_pt"],histos [det+"track_pt"])
    histos[det+"effRecHit_eta"]=R.TGraphAsymmErrors(histos[det+"matchedRecHitTrack_eta"],histos[det+"track_eta"])
    histos[det+"effRecHit_phi"]=R.TGraphAsymmErrors(histos[det+"matchedRecHitTrack_phi"],histos[det+"track_phi"])

    histos[det+"effBestRecHit_pt"]=R.TGraphAsymmErrors( histos[det+"matchedBestRecHitTrack_pt"], histos[det+"track_pt"])
    histos[det+"effBestRecHit_eta"]=R.TGraphAsymmErrors(histos[det+"matchedBestRecHitTrack_eta"],histos[det+"track_eta"])
    histos[det+"effBestRecHit_phi"]=R.TGraphAsymmErrors(histos[det+"matchedBestRecHitTrack_phi"],histos[det+"track_phi"])

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
