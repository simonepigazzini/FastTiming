#!/bin/tcsh


set temp=(`getopt -s tcsh --long ntuples,plots -- $argv:q`)
eval set argv=\($temp:q\)

set ntuples=0
set plots=1

#while (1)
#    switch($1:q)
#    case --ntuples:
#	echo "Running ntuplizer" ; set ntuples=1; shift; 
#	breaksw;
#    case --plots:
#	echo "Running plots" ; set plots=1; shift;
#	breaksw;
#    case --:
#	 shift;
#	 break;
#    default:
#	echo "Internal error!" ; exit 1;
#    endsw
#end

if ( $ntuples == 1) then
#TILE
cmsRun runHits_cfg.py eosdirs=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_tile_v2/ pattern='step3*.root' output=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_tile_v2/DumpHits.root crysLayout=tile > & ! /tmp/dumpHits_tile_singleMu.log &
cmsRun runHits_cfg.py eosdirs=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_tile_v2/ pattern='step3*.root' output=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_tile_v2/DumpHits.root crysLayout=tile > & ! /tmp/dumpHits_tile_singlePi.log &
#BARZ
cmsRun runHits_cfg.py eosdirs=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_barz_v2/ pattern='step3*.root' output=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_barz_v2/DumpHits.root crysLayout=barzflat > & ! /tmp/dumpHits_barz_singleMu.log &
cmsRun runHits_cfg.py eosdirs=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_barz_v2/ pattern='step3*.root' output=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_barz_v2/DumpHits.root crysLayout=barzflat > & ! /tmp/dumpHits_barz_singlePi.log &
# BARPHI
cmsRun runHits_cfg.py eosdirs=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_barphi_v2/ pattern='step3*.root' output=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_barphi_v2/DumpHits.root crysLayout=barphi > & ! /tmp/dumpHits_barphi_singleMu.log &
cmsRun runHits_cfg.py eosdirs=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_barphi_v2/ pattern='step3*.root' output=/eos/cms/store/group/cafcern/meridian/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_barphi_v2/DumpHits.root crysLayout=barphi > & ! /tmp/dumpHits_barphi_singlePi.log &
endif

if ( $plots == 1) then
python drawTrackPerformance.py --input=~/eoscaf/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_tile_v2/DumpHits.root --output=tile_SingleMu.root --layout=tile
python drawTrackPerformance.py --input=~/eoscaf/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_barz_v2/DumpHits.root --output=barzflat_SingleMu.root --layout=barzflat
python drawTrackPerformance.py --input=~/eoscaf/MTD/10_4_0_mtd2_test/SingleMu_FlatPt_BTL_barphi_v2/DumpHits.root --output=barphi_SingleMu.root --layout=barphi

python drawTrackPerformance.py --input=~/eoscaf/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_tile_v2/DumpHits.root --output=tile_singlePi.root --layout=tile
python drawTrackPerformance.py --input=~/eoscaf/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_barz_v2/DumpHits.root --output=barzflat_singlePi.root --layout=barzflat
python drawTrackPerformance.py --input=~/eoscaf/MTD/10_4_0_mtd2_test/SinglePi_FlatPt_BTL_barphi_v2/DumpHits.root --output=barphi_singlePi.root --layout=barphi
endif
