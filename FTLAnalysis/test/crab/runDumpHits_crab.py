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

        config.General.requestName = 'runHits_analysis'
        config.General.workArea = options.workArea
        config.General.transferOutputs = True
        config.General.transferLogs = False
        
        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'runHits_cfg.py'
        config.JobType.pyCfgParams = ['useMTDTrack=True','crysLayout=barzflat','output=DumpHits.root']

        config.Data.inputDataset = None
        config.Data.inputDBS = 'phys03'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 1
        config.Data.outLFNDirBase = '/store/user/meridian/MTD'
        config.Data.publication = False
        config.Data.outputDatasetTag = '10_4_0_mtd3_runHits_analysis_v5'
        config.Data.allowNonValidInputDataset = True
        config.Data.useParent = True

        config.Site.storageSite = 'T2_CH_CERN'
        config.User.voRole = 'priorityuser'

        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [
#            '/RelValDYToLL_M_50_14TeV/meridian-CMSSW_10_4_0_mtd2_patch1-103X_upgrade2023_realistic_v2_2023D35noPU-v1-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValSingleMuFlatPt_0p7to10_pythia8/meridian-RelValSingleMuFlatPt0p7to10pythia8CMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v2-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValMinBias_14TeV/meridian-RelValMinBias14TeVCMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v1-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/meridian-RelValSinglePiFlatPt0p7to10pythia8cfiCMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v2-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValSingleKaonFlatPt_0p7to10/meridian-RelValSingleKaonFlatPt0p7to10CMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v1-479d09f3e9ff3659dd49d9006e28a0a3/USER'
#### V2 ### chi2 cut @ 50
#            '/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/meridian-RelValSinglePiFlatPt0p7to10pythia8cfiCMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v2V2-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#### V3 ### chi2 cut @ 1000
#            '/RelValSinglePiFlatPt_0p7to10_pythia8_cfi/meridian-RelValSinglePiFlatPt0p7to10pythia8cfiCMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v2V3-479d09f3e9ff3659dd49d9006e28a0a3/USER'
#            '/RelValDYToLL_M_50_14TeV/meridian-RelValDYToLLM5014TeVCMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v1V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValMinBias_14TeV/meridian-RelValMinBias14TeVCMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v1V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValSingleMuFlatPt_0p7to10_pythia8/meridian-RelValSingleMuFlatPt0p7to10pythia8CMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v2V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValSingleKaonFlatPt_0p7to10/meridian-RelValSingleKaonFlatPt0p7to10CMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v1V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValSingleProtonFlatPt_0p7to10/meridian-RelValSingleProtonFlatPt0p7to10CMSSW1040mtd2patch1-103Xupgrade2023realisticv22023D35noPU-v1V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
            '/RelValDYToLL_M_50_14TeV/meridian-RelValDYToLLM5014TeVCMSSW1040mtd2patch1-PU25ns103Xupgrade2023realisticv22023D35PU200-v1V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/DYToLL_M-50_14TeV_pythia8/meridian-DYToLLM-5014TeVpythia8PhaseIIMTDTDRAutumn18DR-PU200pilot103Xupgrade2023realisticv2ext2-v2V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValNuGun/meridian-RelValNuGunCMSSW1040mtd2patch1-PU25ns103Xupgrade2023realisticv22023D35PU200-v1V3-479d09f3e9ff3659dd49d9006e28a0a3/USER',
#            '/RelValSingleMuFlatPt_0p7to10/meridian-RelValSingleMuFlatPt0p7to10CMSSW1040mtd2patch1-PU25ns103Xupgrade2023realisticv22023D35PU200-v1V3-479d09f3e9ff3659dd49d9006e28a0a3/USER'
            ]

        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            config.General.requestName = 'runHits_%s' % (inDS.split('/')[1])
            config.General.requestName = config.General.requestName.translate(None, '_')
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s_v3' % (config.General.requestName)
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
