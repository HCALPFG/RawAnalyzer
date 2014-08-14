import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register(
	'file','',VarParsing.VarParsing.multiplicity.singleton,
	VarParsing.VarParsing.varType.string,
	'File path for storing output')
options.parseArguments()
file_path = options.file
#print file_path

process = cms.Process("RawAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) # -1 means run on all events

#default is HcalTBSource but you can change to PoolSource if you like
process.source = cms.Source("HcalTBSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/user/d/drew/USC_223708.root'
        'file:/afs/cern.ch/user/d/drew/USC_223495.root'  #HO pedestal
#        '/store/group/comm_hcal/LS1/USC_223495.root'      #HO pedestal, local
#         '/store/group/comm_hcal/LS1/USC_222759.root'
#         '/store/group/comm_hcal/LS1/USC_223775.root'
#	  '/store/group/comm_hcal/LS1/USC_224285.root'    #not sure, takend 31/7/2014
#         '/store/group/comm_hcal/LS1/USC_224625.root'
    )
)

process.analyzer = cms.EDAnalyzer('RawAnalyzer',
	outputFile = cms.untracked.string(file_path)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 2000  #type out ever <n> events
process.p = cms.Path(process.analyzer)
