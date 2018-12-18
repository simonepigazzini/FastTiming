#!/bin/tcsh

python plot.py --input barzflat_SingleMu.root --output /eos/home-m/meridian/www/plots/MTD/10_4_0_mtd2/barzflat_singleMu --title 'Bar Z (3x56 mm^{2}) Flat (4.5mm gap) Single #mu p_{T} [0.7-4] GeV'
python plot.py --input barphi_SingleMu.root --output /eos/home-m/meridian/www/plots/MTD/10_4_0_mtd2/barphi_singleMu --title 'Bar#phi (3x50 mm^{2} )stag Single#mu p_{T} [0.7-4] GeV'
python plot.py --input tile_SingleMu.root --output /eos/home-m/meridian/www/plots/MTD/10_4_0_mtd2/tile_singleMu --title 'Tile (11.5x11.5 mm^{2}) Single#mu p_{T} [0.7-4] GeV'

python plot.py --input barzflat_singlePi.root --output /eos/home-m/meridian/www/plots/MTD/10_4_0_mtd2/barzflat_singlePi --title 'Bar Z (3x56 mm^{2}) Flat (4.5mm gap) Single #pi p_{T} [0.7-4] GeV'
python plot.py --input barphi_singlePi.root --output /eos/home-m/meridian/www/plots/MTD/10_4_0_mtd2/barphi_singlePi --title 'Bar#phi (3x50 mm^{2} )stag Single#pi p_{T} [0.7-4] GeV'
python plot.py --input tile_singlePi.root --output /eos/home-m/meridian/www/plots/MTD/10_4_0_mtd2/tile_singlePi --title 'Tile (11.5x11.5 mm^{2}) Single#pi p_{T} [0.7-4] GeV'
