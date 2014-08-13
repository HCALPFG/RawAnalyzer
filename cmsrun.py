import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing()

options.register('file',
                 '', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "file path for storing output")

options.parseArguments()

write_to_file = bool ( options.file != "" ) 
file_path = options.file

process = cms.Process("Analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
# -1 means run on all events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#default is HcalTBSource but you can change to PoolSource if you like
process.source = cms.Source("HcalTBSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/user/d/drew/USC_223708.root'
#        'file:/afs/cern.ch/user/d/drew/USC_223495.root'
#        '/store/group/comm_hcal/LS1/USC_223495.root'
         # '/store/group/comm_hcal/LS1/USC_222759.root'
         '/store/group/comm_hcal/LS1/USC_224625.root'
#         '/store/group/comm_hcal/LS1/USC_223775.root'
    )
)

process.analyzer = cms.EDAnalyzer('RawAnalyzer',
   writeToFile = cms.untracked.bool ( write_to_file ),
   filePath = cms.untracked.string ( file_path )                                   
)
#only type out every 2000 events, up to you
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.p = cms.Path(process.analyzer)
