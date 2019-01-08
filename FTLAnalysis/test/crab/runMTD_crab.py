#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException


def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = '',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()
    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.requestName = 'MTDRECO'
        config.General.workArea = options.workArea
        config.General.transferOutputs = True
        config.General.transferLogs = False

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'runMTD_cfg.py'
        config.JobType.pyCfgParams = ['runMTDReco=True','crysLayout=barzflat','output=mtdReco.root']

        config.Data.inputDataset = None
        config.Data.inputDBS = 'global'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 1
        config.Data.outLFNDirBase = '/store/group/dpg_mtd/comm_mtd/meridian/10_4_0_mtd3'
        config.Data.publication = True
        config.Data.outputDatasetTag = None
        config.Data.allowNonValidInputDataset = True

        config.Site.storageSite = 'T2_CH_CERN'
#        config.Site.whitelist = [ 'T2_CH_CERN','T2_US_Nebraska','T2_US_Wisconsin' ]

        config.User.voRole = 'priorityuser'

        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [
#                         '/RelValDYToLL_M_50_14TeV/CMSSW_10_4_0_mtd2_patch1-103X_upgrade2023_realistic_v2_2023D35noPU-v1/GEN-SIM-RECO',
#                          '/RelValDYToLL_M_50_14TeV/CMSSW_10_4_0_mtd2_patch1-PU25ns_103X_upgrade2023_realistic_v2_2023D35PU200-v1/GEN-SIM-RECO',
                          '/RelValSingleMuFlatPt_0p7to10/CMSSW_10_4_0_mtd2_patch1-PU25ns_103X_upgrade2023_realistic_v2_2023D35PU200-v1/GEN-SIM-RECO',
                          '/RelValNuGun/CMSSW_10_4_0_mtd2_patch1-PU25ns_103X_upgrade2023_realistic_v2_2023D35PU200-v1/GEN-SIM-RECO',
                          '/RelValMinBias_14TeV/CMSSW_10_4_0_mtd2_patch1-103X_upgrade2023_realistic_v2_2023D35noPU-v1/GEN-SIM-RECO',
                          '/RelValSingleMuFlatPt_0p7to10_pythia8/CMSSW_10_4_0_mtd2_patch1-103X_upgrade2023_realistic_v2_2023D35noPU-v2/GEN-SIM-RECO',
                          '/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/CMSSW_10_4_0_mtd2_patch1-103X_upgrade2023_realistic_v2_2023D35noPU-v2/GEN-SIM-RECO',
                          '/RelValSingleKaonFlatPt_0p7to10/CMSSW_10_4_0_mtd2_patch1-103X_upgrade2023_realistic_v2_2023D35noPU-v1/GEN-SIM-RECO',
                          '/RelValSingleProtonFlatPt_0p7to10/CMSSW_10_4_0_mtd2_patch1-103X_upgrade2023_realistic_v2_2023D35noPU-v1/GEN-SIM-RECO',
                          '/DYToLL_M-50_14TeV_pythia8/PhaseIIMTDTDRAutumn18DR-PU200_pilot_103X_upgrade2023_realistic_v2_ext2-v2/GEN-SIM-RECO'
                        ]

        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            config.General.requestName = '%s_%s_V4' % (inDS.split('/')[1],inDS.split('/')[2])
            config.General.requestName = config.General.requestName.translate(None, '_')
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s' % (config.General.requestName)
            # Submit.
            try:
                print "Submitting for input dataset %s" % (inDS)
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Submission for input dataset %s failed: %s" % (inDS, hte.headers)
            except ClientException as cle:
                print "Submission for input dataset %s failed: %s" % (inDS, cle)

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print "-"*len(msg)
            print msg
            print "-"*len(msg)
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
            except ClientException as cle:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)


if __name__ == '__main__':
    main()
