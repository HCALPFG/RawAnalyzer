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
#process.source = cms.Source("HcalTBSource",
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/user/d/drew/USC_223708.root'
#        'file:/afs/cern.ch/user/d/drew/USC_223495.root'  #HO pedestal
#        '/store/group/comm_hcal/LS1/USC_223495.root'      #HO pedestal, local
#         '/store/group/comm_hcal/LS1/USC_222759.root'
#         '/store/group/comm_hcal/LS1/USC_223775.root'
#	  '/store/group/comm_hcal/LS1/USC_224285.root'    #not sure, takend 31/7/2014
#         '/store/group/comm_hcal/LS1/USC_224625.root'
'root://eoscms//eos/cms/store/data/Run2015B/HcalNZS/RAW/v1/000/251/244/00000/46080143-C025-E511-9CB7-02163E014166.root'
#'root://eoscms//eos/cms/store/data/Run2015B/HcalNZS/RAW/v1/000/251/883/00000/20C23681-852B-E511-9FBC-02163E01413E.root'
#'root://eoscms//eos/cms/store/data/Run2015B/HcalNZS/RAW/v1/000/251/883/00000/369E8A59-802B-E511-B85E-02163E01259F.root'
#'root://eoscms//eos/cms/store/data/Run2015B/HcalNZS/RAW/v1/000/251/883/00000/488F97C1-8F2B-E511-86B8-02163E0144D2.root'
#'root://eoscms//eos/cms/store/data/Run2015B/HcalNZS/RAW/v1/000/251/883/00000/FAE69354-7E2B-E511-80D7-02163E0125C8.root'
    )
)

process.analyzer = cms.EDAnalyzer('RawAnalyzer',
	debugit = cms.untracked.bool(False),
	outputFile = cms.untracked.string(file_path),
	badevlist = cms.vint32(
	153647285,	152905909,	153143477,	153217205,	151718625,	153024693,	150641153,	151460577,
	152364043,	152889525,	153151669,	151148928,	153471157,	149944833,	151407329,	152529024,
	150403585,	151124352,	152368139,	152451200,	152950965,	153135285,	154125042,	154268402,
	152261643,	150718977,	152737973,	153409717,	153800866,	151321313,	152910005,	153348277,
	154002162,	149846529,	150489601,	150526465,	151370465,	152959157,	153262261,	153916146,
	150202881,	152750261,  153004213),
	modval = cms.untracked.int32(103)
)
process.TFileService = cms.Service("TFileService",fileName = cms.string("RawAnalyzer.root") )
process.MessageLogger.cerr.FwkReport.reportEvery = 2000  #type out ever <n> events
process.p = cms.Path(process.analyzer)
