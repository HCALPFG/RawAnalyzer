// -*- C++ -*-
//
// Package:    RawAnalyzer
// Class:      RawAnalyzer
// 
/**\class RawAnalyzer RawAnalyzer.cc subsystem/RawAnalyzer/src/RawAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Drew Baden
//         Created:  Fri Dec 11 21:26:21 CET 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>
#include <DataFormats/FEDRawData/interface/FEDHeader.h>
#include <DataFormats/FEDRawData/interface/FEDTrailer.h>
#include <DataFormats/FEDRawData/interface/FEDNumbering.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include "TH1.h"
#include "TH2.h"


using namespace edm;
using namespace std;

#define bold  "\33[0;1m"
#define red   "\033[0;31m"        /* 0 -> normal ;  31 -> red */
#define cyan  "\033[1;36m"        /* 1 -> bold ;  36 -> cyan */
#define green "\033[4;32m"        /* 4 -> underline ;  32 -> green */
#define blue  "\033[9;34m"        /* 9 -> strike ;  34 -> blue */
 
#define black  "\033[0;30m"
#define brown  "\033[0;33m"
#define magenta  "\033[0;35m"
#define gray  "\033[0;37m"
 
#define none   "\033[0m"        /* to flush the previous property */

char ESC=27;
//
// HTR types, starting from 0:  HO,HBHE1+,HBHE2+,HBHE3+,HBHE4+,HBHE5+,HBHE6+,HBHE7+,
//                              HF,HBHE1-,HBHE2-,HBHE3-,HBHE4-,HBHE5-,HBHE6-,HBHE7-
const uint32_t ntptype[16] = {0,24,16,16,16,16,16,16,2,24,16,16,16,16,16,16};
const char* htr_types[16] = {"HO","HBHE1+","HBHE2+","HBHE3+","HBHE4+","HBHE5+","HBHE6+","HBHE7+",
                        "HF","HBHE1-","HBHE2-","HBHE3-","HBHE4-","HBHE5-","HBHE6-","HBHE7-"};
const char* err_types[15] = {"CT","HM","TM","FK","CE","LK","BE","CK","OD","LW","FE","RL","EE","BZ","OW"};
bool searching = false;
bool isFEDopen = false;
bool printbegin = true;
int nloopse = 0;
int ncountse = 0;
int find_evn = 0;
int evn, sr, evn1, evn2, orn, htrn, bcn, formatVer, npresamp, ndll, ntps, fw, subV, pipeline, htype;
int tpsamps, qiesamps, ndccw, nwc, reserved, us, cm;
int ow, bz, ee, rl, le, lw, od, ck, be, tmb, hm, ct, odle, odce, odfe, isone, ttcready, dll;
int hspigot[30];
int wspigot[30];
int ospigot[30];
int spigots = 0;
bool redirect = false;
FILE* fout = NULL;
uint32_t* hdata;  //here is the 16-bit HTR payload data
uint32_t n16;     //number of 16-bit words for this HTR payload
uint32_t* tpgs;
uint32_t nqies;
uint32_t* qies;
uint32_t mant,range,capid,dv,er,addr,fib;
int nFEDs = 0;
int iFEDs[48];
const uint32_t* payload;
bool debugit = false;
int whichFED = -1;
size_t size = 0;
int criteria = 0;
int fed1, fed2;
bool not12_13 = true;
Handle<FEDRawDataCollection> rawdata;

void payload_return(uint32_t s,
 	int* nwords,int* htrerr,int* lrberr,int* ee,int* ep,int* eb,int* ev,int* et);
void htr_data_from_payload(const uint32_t *payload, int spigot);
void htr_fill_header();
void htr_data_delete();
void print_htr_payload_headers(bool header, int spigot, bool extra);
void tpgs_from_htr();
void tpgs_delete();
void print_htr_tpgs();
void qies_from_htr();
void qies_delete();
void print_htr_qies();
void qies_unpack(int which);
void qies_unpack_word(uint32_t data);
void htr_data_print();
void setup_spigots(bool printout);
bool find_nonzero_tpgs(bool printit);
bool find_self_similar_qies(bool printit);
bool find_htr_header_errors();
bool check_event_numbers(int irun, int iev);
void getFeds(int* nfeds, int*ifeds, bool printout);
bool checkFedBcN(int nfeds, int* ifeds, int* fedBcN);
bool checkFedEvN(int nfeds, int* ifeds, int* fedEvN);
bool checkFedOrN(int nfeds, int* ifeds, int* fedOrN);
int get_spigot_bcn(int ispigot);
int get_spigot_evn(int ispigot);
int get_spigot_orn(int ispigot);
bool checkFedBcNIdle(int nfeds, int* ifeds, int* fedBcNIdle);
bool check_htr_header_errors(int irun,int iev);
//
// class decleration
//

class RawAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RawAnalyzer(const edm::ParameterSet&);
      ~RawAnalyzer();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      std::set<int> FEDids_;

      // ----------member data ---------------------------
};

void payload_return(uint32_t s,
 	int* nwords,int* htrerr,int* lrberr,int* ee,int* ep,int* eb,int* ev,int* et) {
	*nwords = (s & 0x3FF);
	*htrerr = (s>>24) & 0xFF;
	*lrberr = (s>>16) & 0xFF;
	*ee = (s>>15) & 0x1;
	*ep = (s>>14) & 0x1;
	*eb = (s>>13) & 0x1;
	*ev = (s>>12) & 0x1;
	*et = (s>>11) & 0x1;
}
void htr_fill_header() {
	if (debugit) {
	  cout << "htr_fill_header...hdata 0-7: ";
	  for (int i=0; i<8; i++) printf("0x%X ",hdata[i]);
	  cout << endl;
	}
	evn1 = (hdata[0] & 0xFF);
	sr = hdata[0] >> 15;
	evn = evn1 + (hdata[1] * 256);
	if (debugit) cout << dec << "EVN " << evn << " HTRn " << htrn << endl;
	ow = (hdata[2] & 0x1);		//overflow warning
	bz = (hdata[2]>>1) & 0x1;		//internal buffers busy, fast report
	ee = (hdata[2]>>2) & 0x1;		//empty event
	rl = (hdata[2]>>3) & 0x1;		//rejected previous L1A
	le = (hdata[2]>>4) & 0x1;		//latency error
	lw = (hdata[2]>>5) & 0x1;		//latency warning
	od = (hdata[2]>>6) & 0x1;		//optical data error
	ck = (hdata[2]>>7) & 0x1;		//clocking problems
	be = (hdata[2]>>8) & 0x1;		//bunch error
	odle = (hdata[2]>>9) & 0x1;               //link er if OD set
	odce = (hdata[2]>>10) & 0x1;               //capid er if OD set
	odfe = (hdata[2]>>11) & 0x1;               //fe er if OD set
	tmb = (hdata[2]>>12) & 0x1;	//test mode (1=patterns)
	hm = (hdata[2]>>13) & 0x1;	//histo mode
	ct = (hdata[2]>>14) & 0x1;	//calibration trigger
	isone = (hdata[2]>>15) & 0x1;	//should be 1
	if (debugit) cout << hex << "hdata[2]=0x" << hdata[2] << dec << endl;
	orn = hdata[3] >> 11;
	htrn = hdata[3] & 0x7FF;	
	bcn = hdata[4] & 0xFFF;
	formatVer = hdata[4] >> 12;
	ntps = (hdata[5]>>8);
	npresamp = (hdata[5]>>3)&0x1F;
	ndll = (hdata[5]>>1)&0x3;
	us = (hdata[6] >> 15);
	cm = (hdata[6] >> 14) & 0x1;
	reserved = (hdata[6]>>12)&0x3;
	subV = (hdata[6]>>12)&0xF;
	fw = hdata[6] & 0xFF;
	htype = (hdata[7]>>8) & 0xFF;
	pipeline = hdata[7]&0xFF;
//	if (htype > 15) cout << "+++++++++ hcalType is " << htype << " illegal!!!!" << endl;
	if (debugit) cout << "htype " << htype;
	if (htype == 0) {
	  tpsamps = 0;
	}
	else {
	  tpsamps = ntps/ntptype[htype];
	}
	qiesamps = hdata[n16-4]>>11;
	ndccw = hdata[n16-2];
	nwc = hdata[n16-3];
	ttcready = hdata[5] & 0x1;
	dll = (hdata[5]>>1) & 0x3;
	evn2 = hdata[n16-1] >> 8;
  	if (debugit) cout << "htr_fill_header done" << endl;
}

void htr_data_print() {
	//
	// printout in 16 bit words
	//
	//
	htr_fill_header();
	//
	// first 9 words are the header
	//
	printf("Header (Ver 3, 0x59 and onwards)\n");
	printf("Header 1: 0x%4.4X   {SR=%d,0000000,EVN[7:0]=0x%X}\n",hdata[0],sr,evn1);
	printf("Header 2: 0x%4.4X   {EVN[23:8}=%d, full EVN=%d}\n",hdata[1],hdata[1],evn);
	printf("Header 3: 0x%4.4X   {1,CT,HM,TM,FK,CE,LK,BE,CK,OD,LW,LE,RL,EE,BZ,OW}\n",hdata[2]);
	printf("Header 4: 0x%4.4X   {ORN[4:0]=%d,HTRmodnumber[10:0]=%d}\n",hdata[3],orn,htrn);
	printf("Header 5: 0x%4.4X   {Ver[3:0]=%d,BCN[11:0]=%d}\n",hdata[4],formatVer,bcn);
	printf("Header 6: 0x%4.4X   {#TPs[7:0]=%d,#Presamples[4:0]=%d,DLL[1:0]=%d,TTCready=%d}\n",
		hdata[5],ntps,npresamp,ndll,ttcready);
	printf("Header 7: 0x%4.4X   {US=%d,CM=%d,Reserved=%d,SubV[3:0]=%d,fw[7:0]=0x%X}\n",
			hdata[6],us,cm,reserved,subV,fw);
	printf("Header 8: 0x%4.4X   {HCALtype[7:0]=%d (is %s), Pipelength[7:0]=%d}\n",
		hdata[7],htype,htr_types[htype],pipeline);
	if (redirect) {
		fprintf(fout,"Header (Ver 3, 0x59 and onwards)\n");
		fprintf(fout,"Header 1: 0x%4.4X   {SR=%d,0000000,EVN[7:0]=0x%X}\n",hdata[0],sr,evn1);
		fprintf(fout,"Header 2: 0x%4.4X   {EVN[23:8}=%d, full EVN=%d}\n",hdata[1],hdata[1],evn);
		fprintf(fout,"Header 3: 0x%4.4X   {1,CT,HM,TM,FK,CE,LK,BE,CK,OD,LW,LE,RL,EE,BZ,OW}\n",hdata[2]);
		fprintf(fout,"Header 4: 0x%4.4X   {ORN[4:0]=%d,HTRmodnumber[10:0]=%d}\n",hdata[3],orn,htrn);
		fprintf(fout,"Header 5: 0x%4.4X   {Ver[3:0]=%d,BCN[11:0]=%d}\n",hdata[4],formatVer,bcn);
		fprintf(fout,"Header 6: 0x%4.4X   {#TPs[7:0]=%d,#Presamples[4:0]=%d,DLL[1:0]=%d,TTCready=%d}\n",
		hdata[5],ntps,npresamp,ndll,ttcready);
		fprintf(fout,"Header 7: 0x%4.4X   {US=%d,CM=%d,Reserved=%d,SubV[3:0]=%d,fw[7:0]=0x%X}\n",
			hdata[6],us,cm,reserved,subV,fw);
		fprintf(fout,"Header 8: 0x%4.4X   {HCALtype[7:0]=%d (is %s), Pipelength[7:0]=%d}\n",
		hdata[7],htype,htr_types[htype],pipeline);
	}
	int ipt = 8;
	//
	// next come the TPGs
	//
	int format_version = hdata[4] >> 12;
	if (format_version == 5) {
		printf(
		"%3d TPG channels, %2d time samples per channel, format {SLBid[2:0],SLBch[1:0],Z,SOI,TP[8:0]}\n",
			ntptype[htype],tpsamps);
		if (redirect) fprintf(fout,
		  "%3d TPG channels, %2d time samples per channel, format {SLBid[2:0],SLBch[1:0],Z,SOI,TP[8:0]}\n",
			ntptype[htype],tpsamps);
		int jpt = 0;
		for (int i=0; i<(int) ntptype[htype]; i++) {
		  printf("%2d: ",jpt+1);
		  if (redirect) fprintf(fout,"%2d: ",jpt+1);
		  jpt++;
		  for (int j=0; j<tpsamps; j++) {
		    printf("0x%4.4X={%d,%d,%d,%d,0x%3.3X} ",hdata[ipt],hdata[ipt]>>13,
		        (hdata[ipt]>>11)&3,(hdata[ipt]>>10)&1,(hdata[ipt]>>9)&1,hdata[ipt]&0x1FF); 
		    if (redirect) fprintf(fout,"0x%4.4X={%d,%d,%d,%d,0x%3.3X} ",hdata[ipt],hdata[ipt]>>13,
		        (hdata[ipt]>>11)&3,(hdata[ipt]>>10)&1,(hdata[ipt]>>9)&1,hdata[ipt]&0x1FF); 
		    ipt++;
		  }
		  printf("\n");
		  if (redirect) fprintf(fout,"\n");
		}
	}
	else if (format_version == 6) {
		//
		// version 6 means we can do zero suppression, so it's not trivial to know what has been 
		// suppressed.  it is worth going into it a bit.
		//
		// each HTR type has a different number of TPs that it sends to the RCT, and this is as in the
		// array ntptype above:  0 for HO, 2 for HF, 24 for HBHE1+,HBHE1-, and 16 otherwise
		// this is because HO has no TPs, HF sends regions (2 per board), HBHE1 has no longitudinal
		// segmentation, and HBHEn (n>1) has some segmentation.  voila.
		// so, here are the relevant numbers:
		// NTPS is the number of TP words in the data (from header word 6), and the HTR type is from
		// header word 8.   The number of TPs sent to the RCT for that type comes from ntptype[htype].
		// Call that number TP.   Then the number of TP samples (call it TPS) kept is easy to calculate:
		//
		// NTPS = TPS * TP
		//
		// however, if there's zero suppression, all TPS of them are not put into the data here.
		// Unfortunately, there's no way of knowing TPS from looking at this data, you have to know it
		// some other way (via one of the many HCAL data bases, and good luck to you)
		//
		if (htype == 0) {
			//
			// HO is different than HBHE and HF as far as TPG go
			//
			printf("%3d TPG channels present, format {00000ZS0,bits[8:1]}\n",ntps);
			if (redirect) fprintf(fout,"%3d TPG channels present, format {00000ZS0,bits[8:1]}\n",ntps);
			int jpt = 8;
			for (int i=0; i<ntps; i++) {
				int tword = hdata[jpt+i];
				int high5 = (tword>>11)&0x1F;
				int z = (tword>>10)&0x1;
				int soi = (tword>>9)&0x1;
				int mbits = tword&0xFF;
				printf("0x%4.4X={high5=0x%2.2X,z=%d,soi=%d,muon_bits=0x%2.2X}\n",
					tword,high5,z,soi,mbits);				
				if (redirect) fprintf(fout,"0x%4.4X={high5=0x%2.2X,z=%d,soi=%d,muon_bits=0x%2.2X}\n",
					tword,high5,z,soi,mbits);				
			}
		}
		else {
			printf(
			"%3d TPG channels present, format {SLBid[2:0],SLBch[1:0],Z,SOI,TP[8:0]}\n",	ntps);
			if (redirect) fprintf(fout,
			"%3d TPG channels present, format {SLBid[2:0],SLBch[1:0],Z,SOI,TP[8:0]}\n",	ntps);
			int jpt = 8;
			for (int i=0; i<ntps; i++) {
				int tword = hdata[jpt+i];
				int tp_tag = (tword>>11)&0x1F;
				int slb_id = (tp_tag >> 2)&0x7;
				int slb_ch = tp_tag & 0x3;
				int z = (tword>>10)&0x1;
				int soi = (tword>>9)&0x1;
				int tp = (tword>>8)&0x1FF;
				printf("0x%4.4X={slb_id=%d,slb_ch=%d,z=%d,soi=%d,tp=%d} ",
					tword,slb_id,slb_ch,z,soi,tp);
				if (redirect) fprintf(fout,"0x%4.4X={slb_id=%d,slb_ch=%d,z=%d,soi=%d,tp=%d} ",
					tword,slb_id,slb_ch,z,soi,tp);
			}
		}
	}
	else {
		//
		// complain like crazy
		//
		cout << "HEY, I HAVE NO IDEA WHAT FORMAT VERSION " << format_version << " MEANS!!!" << endl;
	}
	//
	// next come the QIE data, time to figure out which flavor.  5=old, 6="new"
	//
	int nqiewords = 0;
	if (format_version == 5) {
		int jpt = 0;
		nqies = n16 - ntps - 20;
		int nqiech = nqies/qiesamps;
		printf("%3d QIE channels, %2d time samples per channel\n ",nqiech,qiesamps);
//		printf("format {Fiber[2:0],QIE[1:0],ER,DV,Cap[1:0],Range[1:0],Mant[4:0]}:\n");
		printf("   fi    QIE_CH       CAPID       ER          DV     Mant/Range\n");
		if (redirect) {
		  fprintf(fout,"%3d QIE channels, %2d time samples per channel\n",nqiech,qiesamps);
//		  fprintf(fout,"format {Fiber[2:0],QIE[1:0],ER,DV,Cap[1:0],Range[1:0],Mant[4:0]}:\n");
		  fprintf(fout,"   fi    QIE_CH       CAPID       ER          DV     Mant/Range\n");
		}
		for (int i=0; i<nqiech; i++) {
		  printf("%2d: ",jpt+1);
		  if (redirect) fprintf(fout,"%2d: ",jpt+1);
		  jpt++;
		  int* qie_ch = new int[qiesamps];
		  int* capids = new int[qiesamps];
		  int* ers = new int[qiesamps];
		  int* dvs = new int[qiesamps];
		  int* mants = new int[qiesamps];
		  int* ranges = new int[qiesamps];
		  int fibs = -1;
		  for (int j=0; j<qiesamps; j++) {
		    qies_unpack_word(hdata[ipt++]);
		    mants[j] = mant;
		    ranges[j] = range;
		    capids[j] = capid;
		    dvs[j] = dv;
		    ers[j] = er;
		    qie_ch[j] = addr;
		    fibs = fib;
		  }
		  printf("%d ",fibs+1);
		  for (int j=0; j<qiesamps; j++) printf("%d",qie_ch[j]);
		  printf(" ");
		  for (int j=0; j<qiesamps; j++) printf("%d",capids[j]);
		  printf(" ");
		  for (int j=0; j<qiesamps; j++) printf("%d",ers[j]);
		  printf(" ");
		  for (int j=0; j<qiesamps; j++) printf("%d",dvs[j]);
		  printf(" ");
		  for (int j=0; j<qiesamps; j++) printf("%2d/%d ",mants[j],ranges[j]);
		  printf("\n");
		  if (redirect) {
		    fprintf(fout,"%d ",fibs);
		    for (int j=0; j<qiesamps; j++) fprintf(fout,"%d",qie_ch[j]);
		    fprintf(fout," ");
		    for (int j=0; j<qiesamps; j++) fprintf(fout,"%d",capids[j]);
		    fprintf(fout," ");
		    for (int j=0; j<qiesamps; j++) fprintf(fout,"%d",ers[j]);
		    fprintf(fout," ");
		    for (int j=0; j<qiesamps; j++) fprintf(fout,"%d",dvs[j]);
		    fprintf(fout," ");
		    for (int j=0; j<qiesamps; j++) fprintf(fout,"%2d/%d",mants[j],ranges[j]);
		    fprintf(fout,"\n");
		  }
		}
	}
	else if (format_version == 6) {
		cout << "Version 6!" << endl;
		int flavor = 0;
		int jpt = 8 + ntps;
		nqiewords = n16 - jpt - 12;  // 12 for the extra(8), pre-trailer(3), and trailer(1) 
					// the data is packed so this is the number of words, NOT #qies!
		for (int i=0; i<nqiewords; i++) {
		  int iword = hdata[jpt+i];
		  int bit31 = (iword>>15)&1;
		  if (bit31 == 1) {
		    //
		    // channel header
		    //
		    flavor = (iword>>12)&0x7;
		    int errf0 = (iword>>10)&0x1;
		    int errf1 = (iword>>11)&0x1;
		    int capid = (iword>>8)&0x3;
		    int channelid = iword&0xFF;
		    int thisqie = channelid&0x3;
		    int thisfiber = (channelid>>2)&0x7;
		    printf(
		       "0x%4.4X Header: Flavor=%d CapIdErr=%d  LinkErr=%d  1st Capid=%d Qie=%d Fiber=%d\n",
			iword,flavor,errf0,errf1,capid,thisqie,thisfiber);
		    if (redirect) fprintf(fout,
		       "0x%4.4X Header: Flavor=%d CapIdErr=%d  LinkErr=%d  1st Capid=%d Qie=%d Fiber=%d\n",
			iword,flavor,errf0,errf1,capid,thisqie,thisfiber);
		  }
		  else {
		    if (flavor == 5) {
			int qie1 = (iword>>8)&0x7F;
			int r1 = (qie1>>6)&0x3;
			int m1 = qie1&0x1F;
			int qie0 = iword&0x7F;
			int r0 = (qie0>>6)&0x3;
			int m0 = qie0&0x1F;
			printf("0x%4.4X  Qie1(m/r)=%2d/%d  Qie0(m/r)=%2d/%d\n",iword,m1,r1,m0,r0);
			if (redirect) fprintf(fout,"0x%4.4X  Qie1(m/r)=%2d/%d  Qie0(m/r)=%2d/%d\n",
				iword,m1,r1,m0,r0);
		    }
		    else if (flavor == 6) {
			int qie = (iword>>8)&0x7F;
			int r = (qie>>6)&0x3;
			int m = qie&0x1F;
			int capid = (iword>>8)&0x3;
			int dv = (iword>>10)&0x1;
			int er = (iword>>11)&0x1;
			printf("0x%4.4X  Qie(m/r)=%2d/%d  DV=%d  ER=%d  Capid=%d\n",iword,m,r,dv,er,capid);
			if (redirect) fprintf(fout,
				"0x%4.4X  Qie(m/r)=%2d/%d  DV=%d  ER=%d  Capid=%d\n",iword,m,r,dv,er,capid);
		    }
		    else cout << "FLAVOR " << flavor << " IS UNKNOWN AND CONFUSES ME TERRIBLY!!!" << endl;
		  }
		}
	}
	//
	// then the parity word (if 0xFFFF)
	//
	int jpt = 8 + ntps + nqiewords;
	if (hdata[jpt] == 0xFFFF) {
	  printf("Parity word 0xFFFF\n");
	  if (redirect) fprintf(fout,"Parity word 0xFFFF\n");
	  jpt++;
	}
	//
	// and then the 8 "extra" words
	//
	if (cm == 0) {
	  if (us == 0) {
		printf("CM=US=0 Normal mode data\n");
		if (redirect) fprintf(fout,"CM=US=0 Normal mode data\n");
		for (int i=0; i<8; i++) {
		  int iword = hdata[jpt+i];
		  printf("Extra %d: 0x%4.4X  {Empty=%d,Full=%d,LatCnt[1:0]=%d,BCNtime[11:0]=%d\n",
		  	i+1,iword,iword>>15,iword>>14,(iword>>13)&0x3,iword&0xFFF);
		  if (redirect) fprintf(fout,
			"Extra %d: 0x%4.4X  {Empty=%d,Full=%d,LatCnt[1:0]=%d,BCNtime[11:0]=%d\n",
		  	i+1,iword,iword>>15,iword>>14,(iword>>13)&0x3,iword&0xFFF);
		}
	  }
	  else {
	   	printf("CM=0 US=1 Unsuppressed data detected\n");
	   	if (redirect) fprintf(fout,"CM=0 US=1 Unsuppressed data detected\n");
		int iword = hdata[jpt];
		int maskDigi8to1 = iword&0xFF;
		int maskTP8to1 = (iword>>8)&0xFF;
		printf("Extra 1: 0x%4.4X  Mask TP  8- 1: 0x%3.3X  Mask Digi  8- 1: 0x%3.3X\n",
			iword,maskTP8to1,maskDigi8to1);
		if (redirect) fprintf(fout,"Extra 1: 0x%4.4X  Mask TP  8- 1: 0x%3.3X  Mask Digi  8- 1: 0x%3.3X\n",
			iword,maskTP8to1,maskDigi8to1);
		iword = hdata[jpt+1];
		int maskDigi16to9 = iword&0xFF;
		int maskTP16to9 = (iword>>8)&0xFF;
		printf("Extra 2: 0x%4.4X  Mask TP 16- 9: 0x%3.3X  Mask Digi 16- 9: 0x%3.3X\n",
			iword,maskTP16to9,maskDigi16to9);
		if (redirect) fprintf(fout,"Extra 1: 0x%4.4X  Mask TP 16- 9: 0x%3.3X  Mask Digi 16- 9: 0x%3.3X\n",
			iword,maskTP16to9,maskDigi16to9);
		iword = hdata[jpt+2];
		int maskDigi24to17 = iword&0xFF;
		int maskTP24to17 = (iword>>8)&0xFF;
		printf("Extra 3: 0x%4.4X  Mask TP 24-17: 0x%3.3X  Mask Digi  24-17: 0x%3.3X\n",
			iword,maskTP24to17,maskDigi24to17);
		if (redirect) fprintf(fout,"Extra 1: 0x%4.4X  Mask TP 24-17: 0x%3.3X  Mask Digi 24-17: 0x%3.3X\n",
			iword,maskTP24to17,maskDigi24to17);
		iword = hdata[jpt+3];
		int threshDigi1 = iword&0xFF;
		int threshDigi24 = (iword>>8)&0xFF;
		printf("Extra 4: 0x%4.4X  Threshold Digi24: 0x%3.3X  Threshold Digi1: 0x%3.3X\n",
			iword,threshDigi24,threshDigi1);
		if (redirect) fprintf(fout,"Extra 4: 0x%4.4X  Threshold Digi24: 0x%3.3X  Threshold Digi1: 0x%3.3X\n",
			iword,threshDigi24,threshDigi1);
		iword = hdata[jpt+4];
		int syncFcnt1 = iword&0x3;
		int syncFcnt2 = (iword>>2)&0x3;
		int syncFcnt3 = (iword>>4)&0x3;
		int syncFcnt4 = (iword>>6)&0x3;
		int threshTP3to0 = (iword>>12)&0xF;
		printf("Extra 5: 0x%4.4X  Threshold TP[3:0]: 0x%3.3X  Sycn Count Fibers: 4=%d 3=%d 2=%d 1=%d\n",
			iword,threshTP3to0,syncFcnt4,syncFcnt3,syncFcnt2,syncFcnt1);
		if (redirect) fprintf(fout,
			"Extra 5: 0x%4.4X  Threshold TP[3:0]: 0x%3.3X  Sycn Count Fibers: 4=%d 3=%d 2=%d 1=%d\n",
			iword,threshTP3to0,syncFcnt4,syncFcnt3,syncFcnt2,syncFcnt1);
		iword = hdata[jpt+5];
		int syncFcnt5 = iword&0x3;
		int syncFcnt6 = (iword>>2)&0x3;
		int syncFcnt7 = (iword>>4)&0x3;
		int syncFcnt8 = (iword>>6)&0x3;
		int threshTP7to4 = (iword>>12)&0xF;
		printf("Extra 6: 0x%4.4X  Threshold TP[7:4]: 0x%3.3X  Sycn Count Fibers: 8=%d 7=%d 6=%d 5=%d\n",
			iword,threshTP7to4,syncFcnt8,syncFcnt7,syncFcnt6,syncFcnt5);
		if (redirect) fprintf(fout,
			"Extra 6: 0x%4.4X  Threshold TP[7:4]: 0x%3.3X  Sycn Count Fibers: 8=%d 7=%d 6=%d 5=%d\n",
			iword,threshTP7to4,syncFcnt8,syncFcnt7,syncFcnt6,syncFcnt5);
		iword = hdata[jpt+6];
		int bcnofrxbc0 = iword&0xFFF;
		int ZSmask18to16 = (iword>>12)&0x7;
		printf("Extra 7: 0x%4.4X  ZS_Mask[18:16]: 0x%2.2X  BCN-of-RxBC0 = %d\n",
			iword,ZSmask18to16,bcnofrxbc0);
		if (redirect) fprintf(fout,"Extra 7: 0x%4.4X  ZS_Mask[18:16]: 0x%2.2X  BCN-of-RxBC0 = %d\n",
			iword,ZSmask18to16,bcnofrxbc0);
		iword = hdata[jpt+7];
		int ZSmask15to0 = iword;
		printf("Extra 8: 0x%4.4X  ZS_Mask[15:0]: 0x%4.4X\n",iword,ZSmask15to0);
		if (redirect) fprintf(fout,"Extra 8: 0x%4.4X  ZS_Mask[15:0]: 0x%4.4X\n",iword,ZSmask15to0);
	  }
	}
	else {
		printf("CM=1 Compact mode, no extra words here\n");
		if (redirect) fprintf(fout,"CM=1 Compact mode, no extra words here\n");
	}
	//
	// then 4 trailer words
	//
	jpt = jpt + 8;
	printf("Pre-Trailer 3: 0x%4.4X  {DAQsamples[4:0]=%d,DAQwords[10:0]=%d}\n",
		hdata[jpt],hdata[jpt]>>11,hdata[jpt]&0x3FF);
        printf("Pre-Trailer 2: 0x%4.4X  {CRC}\n",hdata[jpt+1]);
//	printf("Pre-Trailer 2: 0x%4.4X  {0000,WordCount[11:0]=%d}\n",hdata[ipt+1],hdata[ipt+1]&0xFFF);
	printf("Pre-Trailer 1: 0x%4.4X  {all zeros}\n",hdata[jpt+2]);
	printf("Trailer:       0x%4.4X  {EVN[7:0]=0x%X,00000000}\n",hdata[jpt+3],hdata[jpt+3]>>8);
	if (redirect) {
	  fprintf(fout,"Pre-Trailer 3: 0x%4.4X  {DAQsamples[4:0]=%d,DAQwords[10:0]=%d}\n",
		hdata[jpt],hdata[jpt]>>11,hdata[jpt]&0x3FF);
          fprintf(fout,"Pre-Trailer 2: 0x%4.4X  {CRC}\n",hdata[jpt+1]);
//        fprintf(fout,"Pre-Trailer 2: 0x%4.4X  {0000,WordCount[11:0]=%d}\n",hdata[ipt+1],hdata[ipt+1]&0xFFF);
	  fprintf(fout,"Pre-Trailer 1: 0x%4.4X  {all zeros}\n",hdata[jpt+2]);
	  fprintf(fout,"Trailer:       0x%4.4X  {EVN[7:0]=0x%X,00000000}\n",hdata[jpt+3],hdata[jpt+3]>>8);
	}
//	for (int i=0; i<4; i++) {
//	  printf("Trailer %d: 0x%4.4X\n",i+1,hdata[ipt]);
//	  if (redirect) fprintf(fout,"Extra %d: 0x%4.4X\n",i+1,hdata[ipt]);
//	  ipt++;
//	}

}
void htr_data_from_payload(const uint32_t *payload, int spigot) {
	//
	// sheesh, now go to 16 bit words. so we can follow the documentation
	//
	if (debugit) cout << "htr_data_from_payload: ";
	int iptr = ospigot[spigot];
	int nwords = wspigot[spigot];
	n16 = 2*nwords;
	if (debugit) cout << " ospigot " << iptr << " nwords " << nwords;
	hdata = new uint32_t [n16];
	for (int j=0; j<nwords; j++) {
	    hdata[2*j] = payload[iptr+j] & 0xFFFF;
//		    printf("%d %d 0x%x\n",j,2*j,hdata[2*j]);
	    hdata[2*j+1] = (payload[iptr+j] >> 16) & 0xFFFF;
//		    printf("%d %d 0x%x\n",j,2*j+1,hdata[2*j+1]);
	}
	if (debugit) cout << " leaving htr_data_from_payload" << endl;
}

void htr_data_delete() {
	delete [] hdata;
}

void print_htr_payload_headers(bool header, int spigot, bool extra) {
	if (header) {
          printf("                                      ");
	  printf("                 Samples      CHTFCLBCOLFREBO\n");
          printf(" Spigot    EvN    BcN  OrN  HTR   fw  #TP");
	  printf("     Htype  TP  QIE DCCw   TMMKEKEKDWELEZW TTC DLL\n");
        }
	if (debugit) cout << "calling htr_fill_header" << endl;
	htr_fill_header();
	//
	// the following are all errors if == 1 except for the last one
	//
	if (debugit) cout << "**** nwc *** " << nwc << "  "  << hex << nwc << dec << endl;
	printf("   %2d    %4d   %5d  %4d %3d 0x%2.2X   %2d     0x%2.2x   %2d   %2d  %3d   ",
	    spigot,evn,bcn,orn,htrn,fw,ntps,htype,tpsamps,qiesamps,ndccw);
	// print out the bits that tell about errors, but only if it's set
	int nerror = 0;
	int ierror[15];
	if (ct == 0) cout << " "; 	else {cout << "1";ierror[nerror++]=0;}
	if (hm == 0) cout << " ";	else {cout << "1";ierror[nerror++]=1;}
	if (tmb == 0) cout << " ";	else {cout << "1";ierror[nerror++]=2;}
	if (odfe == 0) cout << " ";	else {cout << "1";ierror[nerror++]=3;}
	if (odce == 0) cout << " ";	else {cout << "1";ierror[nerror++]=4;}
	if (odle == 0) cout << " ";	else {cout << "1";ierror[nerror++]=5;}
	if (be == 0) cout << " ";	else {cout << "1";ierror[nerror++]=6;}
	if (ck == 0) cout << " ";	else {cout << "1";ierror[nerror++]=7;}
	if (od == 0) cout << " ";	else {cout << "1";ierror[nerror++]=8;}
	if (lw == 0) cout << " ";	else {cout << "1";ierror[nerror++]=9;}
	if (le == 0) cout << " ";	else {cout << "1";ierror[nerror++]=10;}
	if (rl == 0) cout << " ";	else {cout << "1";ierror[nerror++]=11;}
	if (ee == 0) cout << " ";	else {cout << "1";ierror[nerror++]=12;}
	if (bz == 0) cout << " ";	else {cout << "1";ierror[nerror++]=13;}
	if (ow == 0) cout << " ";	else {cout << "1";ierror[nerror++]=14;}
	printf("  %2d  %2d",ttcready,dll);
	if (nerror > 0) {
	  cout << " Error bits: ";
	  for (int i=0; i<nerror; i++) cout << err_types[ierror[i]] << " ";
	}
//	printf(" %d",isone);
	if (redirect) {
	  fprintf(fout,"   %2d    %4d   %5d  %4d  0x%2.2X   %2d     0x%2.2x   %2d   %2d  %3d   ",
	    spigot,evn,bcn,htrn,fw,ntps,htype,tpsamps,qiesamps,ndccw);
	  fprintf(fout," %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
		ct,hm,tmb,odfe,odce,odle,be,ck,od,lw,le,rl,ee,bz,ow);
	  fprintf(fout,"  %2d %2d",ttcready,dll);
//	  fprintf(fout," %d",isone);
	}
	if (evn2 == evn1) cout << endl;
	else cout << " ***" << endl;
	//
	// extra stuff?  \033[1m set bold, \033[m resets it  ("\033[7m" is reverse video)
	//
	if (extra) {
		if (ow) printf("\033[1m OW bit set: Overflow warning!\033[m\n");
		if (bz) printf("\033[1m BZ bit set: Internal buffers busy (fast report)\033[m\n");
		if (ee) printf("\033[1m EE bit set: Empty event!\033[m\n");
		if (rl) printf("\033[1m RL bit set: Previous L1A rejected!\033[m\n");
		if (le) printf("\033[1m LE bit set: Latency error detected!\033[m\n");
		if (lw) printf("\033[1m LW bit set: Latency warning!\033[m\n");
		if (od) {
		        printf("\033[1m OD bit set: Optical Data error detected!\033[m");
		        printf("\033[1m   (could be turned on by orbit gap stuff)\033[m\n");
		}
		if (ck)   printf("\033[1m CK bit set: Clock problems, TTC or DLL\033[m\n");
		if (be)   printf("\033[1m BE bit set: Bunch error, internal to HTR counter\033[m\n");
		if (tmb)  printf("\033[1m TM bit set: Must be in pattern mode\033[m\n");
		if (ct)   printf("\033[1m CT bit set: Must be in calibration mode\033[m\n");
		if (!isone) printf("\033[1m Error in header!!!  Payload word 3 should have MSB=1\033[m\n");		
	     if (redirect) {
		if (ow) fprintf(fout,"\033[1m OW bit set: Overflow warning!\033[m\n");
		if (bz) fprintf(fout,"\033[1m BZ bit set: Internal buffers busy (fast report)\033[m\n");
		if (ee) fprintf(fout,"\033[1m EE bit set: Empty event!\033[m\n");
		if (rl) fprintf(fout,"\033[1m RL bit set: Previous L1A rejected!\033[m\n");
		if (le) fprintf(fout,"\033[1m LE bit set: Latency error detected!\033[m\n");
		if (lw) fprintf(fout,"\033[1m LW bit set: Latency warning!\033[m\n");
		if (od) {
		        fprintf(fout,"\033[1m OD bit set: Optical Data error detected!\033[m");
		        fprintf(fout,"\033[1m   (could be turned on by orbit gap stuff)\033[m\n");
		}
		if (ck)   fprintf(fout,"\033[1m CK bit set: Clock problems, TTC or DLL\033[m\n");
		if (be)   fprintf(fout,"\033[1m BE bit set: Bunch error, internal to HTR counter\033[m\n");
		if (tmb)  fprintf(fout,"\033[1m TM bit set: Must be in pattern mode\033[m\n");
		if (ct)   fprintf(fout,"\033[1m CT bit set: Must be in calibration mode\033[m\n");
		if (!isone) fprintf(fout,"\033[1m Error in header!!!  Payload word 3 should have MSB=1\033[m\n");		
	     }
	}
}

void tpgs_from_htr() {
	//
	// collect all of the TPs into arrays
	//
	tpgs = new uint32_t [ntps];
	for (int j=0; j<ntps; j++) tpgs[j] = hdata[j+8];

}

void tpgs_delete() {
	delete [] tpgs;
}

void qies_from_htr() {
	//
	// collect all of the qies into arrays
	//
	nqies = n16 - ntps - 20;
        if (debugit) cout << "qies_from_htr: nqies is " << nqies << endl;
	qies = new uint32_t [nqies];
	for (int j=0; j<(int) nqies; j++) qies[j] = hdata[j+ntps+8];

}

void qies_delete() {
	delete [] qies;
}

void qies_unpack(int which) {
	//
	// "which" goes from 0 to number of DAQ data words from trailer word -4
	//
	qies_unpack_word(qies[which]);
}

void qies_unpack_word(uint32_t data) {

        if (debugit) cout << "qies_unpack: qies[]=0x" << hex << data << dec << endl;
	mant = data & 0x1F;
	range = (data>>5) & 0x3;
	capid = (data>>7) & 0x3;
	dv = (data>>9) & 0x1;
	er = (data>>10) & 0x1;
	addr = (data>>11) & 0x3;
	fib = (data>>13) & 0x7;
}

void print_htr_qies() {
	//
	// qiesamps is the variable that you use to know how many time samples
	//
	int nqiech = nqies/qiesamps;
	printf("There are %d QIE channels present, %d time samples per channel.  mant/range values: \n",
	 	nqiech,qiesamps);
	if (redirect) fprintf(fout,"There are %d QIE channels present, %d time samples per channel.  mant/range values: \n",
		nqiech,qiesamps);
	int ipt = 0;
	for (int i=0; i<nqiech; i++) {
	  printf("%3d: ",i+1);
	  if (redirect) fprintf(fout,"%3d: ",i+1);
	  for (int j=0; j<qiesamps; j++) {
	    if (debugit) cout << "qies_unpack...";
	    qies_unpack(ipt++);
	    if (debugit) cout << "done" << endl;
	    printf("%2d/%d ",mant,range);
	    if (redirect) fprintf(fout,"%2d/%d ",mant,range);
	  }
	  printf("\n");
	  if (redirect) fprintf(fout,"\n");
	}	

}

void print_htr_tpgs() {
	//
	// format this dump 1 line per channel, using number of time samples
	// note: NO ZERO SUPPRESSION!!!!
	//
	uint32_t jptr = 0;
	printf("\n Ch  (Z means 0x0 sent to RCT,  * means SOI)\n");
	for (uint32_t j=0; j<ntptype[htype]; j++) {
	  printf("%3d ",j+1);
	  for (int k=0; k<tpsamps; k++) {
//		    uint32_t tp_tag = (tpgs[jptr] >> 11) & 0x1F;
//		    uint32_t slb_id = tp_tag >> 2;
//		    uint32_t slb_ch = tp_tag & 0x3;
	    uint32_t tp_z = (tpgs[jptr] >> 10) & 0x1;
	    uint32_t tp_soi = (tpgs[jptr] >> 9) & 0x1;
	    uint32_t tp_data = (tpgs[jptr] & 0x1FF);
	    if (tp_z == 1) printf("Z");
	    else printf(" ");
	    printf("0x%3.3X",tp_data);
	    if (tp_soi == 1) printf("* ");
	    else printf("  ");
	    if (redirect) {
	      if (tp_z == 1) fprintf(fout,"Z");
	      else fprintf(fout," ");
	      fprintf(fout,"0x%3.3X",tp_data);
	      if (tp_soi == 1) fprintf(fout,"* ");
	      else fprintf(fout,"  ");
	    }
	    jptr++;
	  }
	  printf("\n");
	  if (redirect) fprintf(fout,"\n");
	}
}

void setup_spigots(bool printout) {

	    int nwords,htrerr,lrberr,ee,ep,eb,ev,et;
	    spigots = 0;
	    if (printout) cout << "  Spigot  #Words  HTRerr  LRBerr  E-P-B-V-T " << endl;
	    for (int i=0; i<18; i++) {
	       uint32_t s = payload[i+6] & 0xFFFFFFFF;
	       payload_return(s,&nwords,&htrerr,&lrberr,&ee,&ep,&eb,&ev,&et);
	       if (nwords > 0) {
	         hspigot[spigots] = s;
		 wspigot[spigots] = nwords;
	         if (printout) printf("  %4d    %5d    0x%2.2x    0x%2.2x   %d-%d-%d-%d-%d\n",
	       	     spigots,nwords,htrerr,lrberr,ee,ep,eb,ev,et);
		 spigots++;
	       }
	    }
	    if (printout) {
		cout << " There are " << spigots << " spigots" << endl;
	        printf(" What does E-P-B-V-T mean you ask?\n");
		printf(
		"  %sE%s = HTR input enabled in DCC         %sP%s = Data received from this HTR for this event\n",
			bold,none,bold,none);
		printf(
		"  %sB%s = BX and ORN from HTR matches DCC  %sV%s = HTR event number matches TTC event number\n",
			bold,none,bold,none);
		printf("  %sT%s = LRB data was truncated on reading\n",bold,none);
	    }
	    ospigot[0] = 24;
	    for (int i=1; i<spigots; i++) ospigot[i] = ospigot[i-1] + wspigot[i-1];
	    if (debugit) for (int i=0; i<spigots; i++) 
	    	cout << dec << "Spigot " << i << " offset " << ospigot[i] << endl;

}

bool find_htr_header_errors() {
	return false;
}

bool find_nonzero_tpgs(bool printit) {

	    bool global = false;

	    for (int i=fed1; i<fed2; i++) {
	      const FEDRawData& data = rawdata->FEDData(i);
	      size=data.size();
              FEDHeader header(data.data());
	      FEDTrailer trailer(data.data()+size-8);
	      payload=(uint32_t*)(data.data());
	      //
	      // setup for all the spigots
	      //
//	      cout << "spigots " << spigots << endl;
	      if (printit) cout << "FED " << i << endl;
	      setup_spigots(false);
	      if (printit) cout << "spigots " << spigots << endl;
	      //
	      // loop over spigots
	      //
	      for (int j=0; j<spigots; j++) {
	        if (printit) cout << "  Spigot " << j << endl;
	        htr_data_from_payload(payload,j);
		htr_fill_header();
		tpgs_from_htr();
//		print_htr_tpgs();
		//
		// find nonzero tpgs
		//
		int jptr = 0;
		int kptr = 0;
//		cout << "looping over TPs for spigot " << j << ", number TPs ntptype[ " << 
//		   htype << "]=" << ntptype[htype] << " tpsamps= " << tpsamps << endl;
		for (int k=0; k<(int) ntptype[htype]; k++) {
		  bool foundit = false;
		  kptr = jptr;
		  for (int l=0; l<tpsamps; l++) {
		    int tpdata = tpgs[jptr] & 0x1FF;
		    if (tpdata > 0) foundit = true;
		    jptr++;
		  }
		  if (foundit) {
		    global = true;
		    cout << "Nonzero TPG FED " << i << " spigot " << j;
		    if (printit) {
		      jptr = kptr;
		      for (int l=0; l<tpsamps; l++) {
		        int tpdata = tpgs[jptr] & 0x1FF;
		        int soi = (tpgs[jptr]>>9) & 0x1;
		        if (tpdata > 0) cout << "    TP=" << k << " TS=" << l << " TPG=0x" <<
		           tpdata << " SOI=" << soi << endl;
		        jptr++;
		      }
		    }
		  }
		}
	      }
	      cout << endl;
	    }
	    return global;
}

bool find_self_similar_qies(bool printit) {
	//
	// loop over all FEDs
	//
	bool global = false;
	for (int i=fed1; i<fed2; i++) {
	      const FEDRawData& data = rawdata->FEDData(i);
	      size=data.size();
              FEDHeader header(data.data());
	      FEDTrailer trailer(data.data()+size-8);
	      payload=(uint32_t*)(data.data());
	      //
	      // setup for all the spigots
	      //
//	      cout << "spigots " << spigots << endl;
	      if (printit) cout << "FED " << i << endl;
	      setup_spigots(false);
	      if (printit) cout << "spigots " << spigots << endl;
	      //
	      // loop over spigots
	      //
	      for (int j=0; j<spigots; j++) {
	        if (printit) cout << "  Spigot " << j << endl;
	        htr_data_from_payload(payload,j);
		htr_fill_header();
		tpgs_from_htr();
//		print_htr_tpgs();
		//
		// find nonzero tpgs
		//
		int jptr = 0;
		int kptr = 0;
//		cout << "looping over TPs for spigot " << j << ", number TPs ntptype[ " << 
//		   htype << "]=" << ntptype[htype] << " tpsamps= " << tpsamps << endl;
		for (int k=0; k<(int) ntptype[htype]; k++) {
		  bool foundit = false;
		  kptr = jptr;
		  for (int l=0; l<tpsamps; l++) {
		    int tpdata = tpgs[jptr] & 0x1FF;
		    if (tpdata > 0) foundit = true;
		    jptr++;
		  }
		  if (foundit) {
		    global = true;
		    cout << "Nonzero TPG FED " << i << " spigot " << j;
		    if (printit) {
		      jptr = kptr;
		      for (int l=0; l<tpsamps; l++) {
		        int tpdata = tpgs[jptr] & 0x1FF;
		        int soi = (tpgs[jptr]>>9) & 0x1;
		        if (tpdata > 0) cout << "    TP=" << k << " TS=" << l << " TPG=0x" <<
		           tpdata << " SOI=" << soi << endl;
		        jptr++;
		      }
		    }
		  }
		}
	      }
	      cout << endl;
	    }
	    return global;
}

void getFeds(int* nfeds, int* ifeds, bool printout) {
	int kfeds = 0;
	for (int i=700; i<732; i++) {
	  const FEDRawData& data = rawdata->FEDData(i);
	  size=data.size();
	  if (size > 0) {
//		cout << "Found FED " << i << endl;
		FEDHeader header(data.data());
		FEDTrailer trailer(data.data()+size-8);
		payload=(uint32_t*)(data.data());
		ifeds[kfeds++] = i;
	  }
	}
	*nfeds = kfeds;
	if (kfeds > 0 && printout) {
		cout << "Feds: ";
		for (int i=0; i<kfeds; i++) cout << " " << ifeds[i];
		cout << endl;
	}
	if (kfeds == 0) cout << " No feds for this event????  " << endl;
}	

bool checkFedBcN(int nfeds, int* ifeds, int* fedBcN) {
	//
	// see if all BcN are the same for all FEDs (should be)
	//
	bool first = true;
	int theBcN = -1;
	int numsame = 1;
	bool okok1 = true;
	for (int i=0; i<nfeds; i++) {
	  const FEDRawData& data = rawdata->FEDData(ifeds[i]);
	  FEDHeader header(data.data());
	  FEDTrailer trailer(data.data()+size-8);
	  payload=(uint32_t*)(data.data());
	  int bcn =  (payload[0] >> 20) & 0xFFF;  // 12 bits
	  if (first) {theBcN = bcn; first = false;}
	  else {if (bcn == theBcN) numsame++;}	  
	  //
	  // what about the spigots?
	  //
	  setup_spigots(false);
	  bool okok = true;
	  for (int j=0; j<spigots; j++) {
	    int sbcn = get_spigot_bcn(j);
	    if (sbcn != bcn) okok = false;
	  }
	  okok1 = okok1 && okok;
	  if (!okok) cout << "  ===> FED " << ifeds[i] << 
		" has BcN mismatch between FED and one of the spigots " << endl;
	}
	*fedBcN = theBcN;
	if (numsame == nfeds && okok1) return true;
	else return false;
}
bool checkFedEvN(int nfeds, int* ifeds, int* fedEvN) {
	//
	// see if all EvN are the same for all FEDs (should be)
	//
	bool first = true;
	int theEvN = -1;
	int numsame = 1;
	bool okok1 = true;
	for (int i=0; i<nfeds; i++) {
	  const FEDRawData& data = rawdata->FEDData(ifeds[i]);
	  FEDHeader header(data.data());
	  FEDTrailer trailer(data.data()+size-8);
	  payload=(uint32_t*)(data.data());
	  int evn =  payload[1] & 0xFFFFFF; // 24 bits
	  if (first) {theEvN = evn; first = false;}
	  else {if (evn == theEvN) numsame++;}	  
	  //
	  // what about the spigots?
	  //
	  setup_spigots(false);
	  bool okok = true;
	  for (int j=0; j<spigots; j++) {
	    int sevn = get_spigot_evn(j);
	    if (sevn != evn) okok = false;
	  }
	  okok1 = okok1 && okok;
	  if (!okok) cout << "  ===> FED " << ifeds[i] << 
		" has EvN mismatch between FED and one of the spigots " << endl;
	}
	*fedEvN = theEvN;
	if (numsame == nfeds && okok1) return true;
	else return false;
}

bool checkFedOrN(int nfeds, int* ifeds, int* fedOrN) {
	//
	// see if all OrN are the same for all FEDs (should be)
	//
	bool first = true;
	int theOrN = -1;
	int numsame = 1;
	bool okok1 = true;
	for (int i=0; i<nfeds; i++) {
	  const FEDRawData& data = rawdata->FEDData(ifeds[i]);
	  FEDHeader header(data.data());
	  FEDTrailer trailer(data.data()+size-8);
	  payload=(uint32_t*)(data.data());
	  int orn = (payload[2] >> 4) & 0xFFFFFFFF;
	  if (first) {theOrN = orn; first = false;}
	  else {if (orn == theOrN) numsame++;}
	  //
	  // what about the spigots?
	  //
	  setup_spigots(false);
	  bool okok = true;
	  for (int j=0; j<spigots; j++) {
	    int sorn = get_spigot_orn(j);
	    if (sorn != (orn & 0x1F) ) okok = false;
	  }
	  okok1 = okok1 && okok;
	  if (!okok) cout << "  ===> FED " << ifeds[i] << 
		" has OrN mismatch between FED and one of the spigots " << endl;
	}
	*fedOrN = theOrN;
	if (numsame == nfeds && okok1) return true;
	else return false;
}

bool checkFedBcNIdle(int nfeds, int* ifeds, int* fedBcNIdle) {
	//
	// see if all BcN Idle words of all spigots are the same for all FEDs (should be)
	//
//	int numsame = 1;
	bool okok1 = true;
	int theidle = -1;
	for (int i=0; i<nfeds; i++) {
	  const FEDRawData& data = rawdata->FEDData(ifeds[i]);
	  FEDHeader header(data.data());
	  FEDTrailer trailer(data.data()+size-8);
	  payload=(uint32_t*)(data.data());
	  //
	  // loop over the spigots?
	  //
	  setup_spigots(false);
	  bool okok = true;
	  for (int j=0; j<spigots; j++) {
	    bool okok2 = true;
	    bool dont = not12_13 && ( (j == 12) || (j == 13) );
	    htr_data_from_payload(payload,j);
	    htr_fill_header();
	    if (!dont) {
	      if ( (cm == 0) && (us == 0) && (formatVer == 6)) {
		// 
		// should be format version 6, and normal data, the only place where the bcn 
		// idle values are stored in the "extra" block
		//
	        int jpt = n16 - 12;  // pointer to extra,
		bool first = true;
		theidle = 0;
		for (int k=0; k<8; k++) {
		  int bcnidle = hdata[jpt+k] & 0xFFF;
		  if (first) {
			theidle = bcnidle;
			first=false;
		  }
		  if (theidle != bcnidle) okok2 = false;
		}
	      }
	    }
	    if (!okok2) {
		cout << "  ===> FED " << ifeds[i] << " spigot " << j << 
			" does not have all BcN IDLEs the same! " << endl;
	    }
	    htr_data_delete();
//	    cout << "after htr_data_delete"<< endl;
	    okok = okok && okok2;
	  }
	  okok1 = okok1 && okok;
	  if (!okok) cout << "  ===> FED " << ifeds[i] << 
		" has BCN Idle mismatch between FED and one of the spigots " << endl;
	}
	*fedBcNIdle = theidle;
	return okok1;
}

int get_spigot_bcn(int ispigot) {
        htr_data_from_payload(payload,ispigot);
	htr_fill_header();
	int bcn = hdata[4] & 0xFFF;
        htr_data_delete();
	return bcn;	
}

int get_spigot_evn(int ispigot) {
        htr_data_from_payload(payload,ispigot);
	htr_fill_header();
	int evn1 = (hdata[0] & 0xFF);
	int evn = evn1 + (hdata[1] * 256);
        htr_data_delete();
	return evn;
}

int get_spigot_orn(int ispigot) {
        htr_data_from_payload(payload,ispigot);
	htr_fill_header();
	int orn = hdata[3] >> 11;  // only 5 bits!
        htr_data_delete();
	return orn;	
}

int fedOrN_previous = 0;
int OrN_delta;
bool check_event_numbers(int irun, int iev) {
	int nfeds = 0;
	int ifeds[50];  // ususally 32 max
	getFeds(&nfeds, ifeds, false);
	if (nfeds < 1) cout << "NO FEDS FOR THIS EVENT!!!!!" << endl;
	int fedBcN, fedEvN, fedOrN, fedBcnIdle;
	//
	// first check if the FEDs and their spigots are all consistent with the same BcN, EvN, and OrN
	//
	bool bcnok = checkFedBcN(nfeds, ifeds, &fedBcN);
	bool evnok = checkFedEvN(nfeds, ifeds, &fedEvN);
	bool ornok = checkFedOrN(nfeds, ifeds, &fedOrN);
	OrN_delta = fedOrN - fedOrN_previous;
	fedOrN_previous = fedOrN;
	printf("Run %7d Ev %6d       Orn %10d : delta %10d\n",irun,iev,fedOrN,OrN_delta);	
	bool bcnidle = checkFedBcNIdle(nfeds, ifeds, &fedBcnIdle);
	bool result = bcnok && evnok && ornok && bcnidle;
	if (!result) cout << "  ===> Problem with BCN/EVN/ORN!!!" << endl;
//	else cout << " Run " << irun << " " << iev << " looks good..." << endl;
//	else cout << "." << endl;
//	return result && (iev < 4279);
//	return result;
	return true;
}

bool check_htr_header_errors(int irun,int iev) {
	int nfeds = 0;
	int ifeds[50];  // ususally 32 max
	getFeds(&nfeds, ifeds, false);
	if (nfeds < 1) cout << "NO FEDS FOR THIS EVENT!!!!!" << endl;
	for (int i=0; i<nfeds; i++) {
	  const FEDRawData& data = rawdata->FEDData(ifeds[i]);
	  FEDHeader header(data.data());
	  FEDTrailer trailer(data.data()+size-8);
	  payload=(uint32_t*)(data.data());
	  setup_spigots(false);
	  for (int j=0; j<spigots; j++) {
		htr_data_from_payload(payload,j);
		htr_fill_header();
		int nerror = 0;
		int ierror[15];
		if (ct == 1) ierror[nerror++]=0;
		if (hm == 1) ierror[nerror++]=1;
		if (tmb == 1) ierror[nerror++]=2;
		if (odfe == 1) ierror[nerror++]=3;
		if (odce == 1) ierror[nerror++]=4;
		if (odle == 1) ierror[nerror++]=5;
		if (be == 1) ierror[nerror++]=6;
		if (ck == 1) ierror[nerror++]=7;
		if (od == 1) ierror[nerror++]=8;
		if (lw == 1) ierror[nerror++]=9;
		if (le == 1) ierror[nerror++]=10;
		if (rl == 1) ierror[nerror++]=11;
		if (ee == 1) ierror[nerror++]=12;
		if (bz == 1) ierror[nerror++]=13;
		if (ow == 1) ierror[nerror++]=14;
		if (nerror > 0) {
		  printf("Run %7d Event %5d FED %3d Spigot %2d Errors:",irun,iev,ifeds[i],j);
		  for (int k=0; k<nerror; k++) cout << err_types[ierror[k]] << " ";
		  cout << endl;
		}
	  }
	}
	ncountse++;
//	if ( (ncountse < nloopse) && (iev < 4279))  return true;
	if (ncountse < nloopse) return true;
	else return false;
}
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RawAnalyzer::RawAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


RawAnalyzer::~RawAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void RawAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& iSetup) {

   int this_evn = e.id().event();
   isFEDopen = false;
   if (printbegin) cout << "---> Run: " << e.id().run() << " Event: " << e.id().event() << endl;
//
//   e.getByType(rawdata);   <=== this is old.  Seth Cooper changed it!  (7/2014)
//
   if (!e.getByLabel("source",rawdata)) {
	cout << "Hmm, getByLabel returns FALSE for run " << e.id().run() << " event " << 
		e.id().event() << endl;
	return;
   }
   if (searching) {
     if (criteria == 0) {
       if (this_evn != find_evn) return;
       searching = false;
     }
     else if (criteria == 1) {
       if (!find_nonzero_tpgs(false)) return;
       searching = false;
     }
     else if (criteria == 2) {
       //
       // search for QIEs that have the same value for all (or most) TS
       //
       if (!find_self_similar_qies(true)) return;
       searching = false;
     }
     else if (criteria == 3) {
       //
       // stop on events where there are any errors in the HTR header
       //
       if (!find_htr_header_errors()) return;
       searching = false;
     }
     else if (criteria == 4) {
       //
       // cycle through and printout stuff about ORN, BCN, EVN, BCN idle, etc
       //
      if (check_event_numbers(int (e.id().run()),int (e.id().event()))) {
//fflush(stdout);
	return;
      }
      printbegin = true;
      searching = false;
     }
     else if (criteria == 5) {
	//
	// report any HTR that has any errors set in it's header.  This will produce a lot of output!
	//
	if (check_htr_header_errors(int (e.id().run()),int (e.id().event()))) return;
        printbegin = true;
	searching = false;
		printf("   %sbit0/OW%s=overflow warning goes away if L1A rate reduced\n",bold,none);
		printf("   %sbit1/BZ%s=internal buffers busy\n",bold,none);
		printf("   %sbit2/EE%s=empty event (because of previous BZ, includes only 1st 5 and last 3 words)\n",
			bold,none);
		printf("   %sbit3/RL%s=trig rule violation (rules 1 adn 2 of TDR 16.4.3)\n",bold,none);
		printf("   %sbit4/FE%s=detects FE idle words or CCA corrupted words that should be suppressed\n",
			bold,none);
		printf("   %sbit5/LW%s=latency warning, Jose fifo has tuned latency since previous Tx\n",
			bold,none);
		printf("   %sbit6/OD%s=so called fiber data error, is OR of 9,10,11 below\n",bold,none);
		printf("   %sbit7/CK%s=clocking error, ~TTCready || DLL_unlock\n",bold,none);
		printf("   %sbit8/BE%s=Bunch count error if BCC does not wrap nicely when BC0 received\n",bold,none);
		printf("   %sbit9/LK%s=OR of 8 fibers, each is ER || ~DV || Frame_error\n",bold,none);
		printf("   %sbit10/CE%s=OR of 8 fibers, each is capid rotation error\n",bold,none);
		printf("   %sbit11/FK%s=OR of 8 fibers, each is if FE flags are not right\n",bold,none);
		printf("   %sbit12/TM%s=test mode (0=data, 1=pattern or counter etc)\n",bold,none);
		printf("   %sbit13/HM%s=histo mode (if 1)\n",bold,none);
		printf("   %sbit13/CT%s=calib trig (if 1, otherwise L1A event)\n",bold,none);
		printf("   %sDLL%s  counts number of times DLL unlock since last reset (trust not)\n",bold,none);
		printf("   %sTTC%s  TTCready, should always be asserted\n",bold,none);
     }
   }
//
//   Handle<FEDRawDataCollection> rawdata;
   nFEDs = 0;
   cout << "+++++++++ Looking for FEDs ++++++++++++++" << endl;
//   for (int i = 700; i<FEDNumbering::lastFEDId(); i++) {
   for (int i = 700; i<732; i++) {
	const FEDRawData& data = rawdata->FEDData(i);
	size=data.size();

	if (size>0 && (FEDids_.empty() || FEDids_.find(i)!=FEDids_.end())) {
	  cout << " Found FED# " << setw(4) << i << " " << setw(8) << size << " bytes " << endl;
	  iFEDs[nFEDs++] = i;
	}
   }
   cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
   //
   // this might be old fashioned but what the heck...
   //	    
   int done = 0;
   while (done == 0) {	      
	  cout << "===========================================";
	  cout << "===========================================" << endl;
	  cout << "  0  Next (^C then 0 to quit)              ";
	  cout << "  1  Select FED and report header info     " << endl;
	  if (redirect)
	  cout << "  2  Close output file, direct ONLY to term";
	  else
	  cout << "  2  Direct output to term AND file        ";
	  cout << "  3  Find event number and stop            " << endl;
	  cout << "  4  DCC hex dump (unformatted)            ";
	  cout << "  5  DCC Event header HTR info             " << endl;
	  cout << "  6  All HTR payload headers               ";
	  cout << "  7  HTR payload (formatted)               " << endl;
	  cout << "  8  HTR payload (unformatted)             ";
	  cout << "  9  Search for anomalies (there are many) " << endl;
	  cout << "What is your pleasure?> ";
	  int iw;
	  scanf("%d",&iw);
	  if (iw == 0) done = 1;
	  else if (iw == 1) {
	    //
	    // select which FED you want to see
	    //
	    cout << "There are " << nFEDs << " here: ";
	    for (int i=0; i<nFEDs; i++) cout << iFEDs[i] << " ";
	    cout << endl;
	    cout << "Select FED (0=all)> ";
	    cin >> whichFED;
	    isFEDopen = true;
	    if (whichFED == 0) cout << "Be careful, you still need to choose 1 FED for other menu items!!!!"<<endl;
	    for (int ifed=0; ifed<nFEDs; ifed++) {
		if (whichFED == 0 || whichFED == iFEDs[ifed]) {
		    const FEDRawData& data = rawdata->FEDData(iFEDs[ifed]);
		    size=data.size();
	            FEDHeader header(data.data());
		    FEDTrailer trailer(data.data()+size-8);
		    cout << " L1Id: " << setw(8) << header.lvl1ID();
		    cout << " BXId: " << setw(4) << header.bxID();
		    cout << endl;
		    payload=(uint32_t*)(data.data());
		    //
		    // very first 4 words...	
		    //
		    cout << "==============================================================================" << endl;
		    cout << " DCC Header (3 64-bit words):" << endl;
//		    printf(" Data = 0x%8.8X%8.8X means:\n",payload[1],payload[0]);
 	 	     int ifirst = payload[0] & 0xf; // should be 0x8
		     int fov = (payload[0] >> 4) & 0xF;  // what is FOV?
		     int fedn = (payload[0] >> 8) & 0xFFF;  // 12 bits
		     int bcn =  (payload[0] >> 20) & 0xFFF;  // 12 bits
		     int evn =  payload[1] & 0xFFFFFF; // 24 bits
		     int evtTy = (payload[1] >> 24) & 0xF;
		      int isfive = (payload[1] >> 28) & 0xF;
		     printf("   %sFirst nibble%s=0x%X (8) ",bold,none,ifirst);
		     printf("%sFOV%s=0x%X %sFED%s=%d %sBcN%s=%d (0x%X) %sEvN%s=%d %sEvty%s=%d (1=phys, 2=calib)",
			bold,none,fov,bold,none,fedn,bold,none,bcn,bcn,bold,none,evn,bold,none,evtTy);
		     printf("%sLast nibble%s=0x%X (5)\n",bold,none,isfive);
//		    printf(" Data = 0x%8.8X%8.8X means:\n",payload[3],payload[2]);
		     int if2 = payload[2] & 0xF;
		     int orn = (payload[2] >> 4) & 0xFFFFFFFF;
		     int reserved = (payload[3] >> 4) & 0xFFFFF;
		     int calTy = (payload[2] >> 24) & 0xF;
		     int last2 = (payload[2] >> 28) & 0xF;
		     printf("   %sFirst nibble%s=0x%X (0) ",bold,none,if2);
		     printf("%sOrN%s=%d %sreserved%s=0x%X %scalTy%s=%d (laser setting calib trig)",
			bold,none,orn,bold,none,reserved,bold,none,calTy);
		     printf("%sLast nibble%s=0x%X (0)\n",bold,none,last2);
		    //
		    // next comes the dcc header/status info
		    //
//		    printf(" Data = 0x%8.8X%8.8X means:\n",payload[5],payload[4]);
		     int formver = payload[4] & 0xFF;
		     int TTS = (payload[4] >> 8) & 0xF;
		     int HTRstatus = (payload[4] >> 14) & 0x7FFF;
		     int DCCstatus = payload[5] & 0x3FF;
		     int DCCrev = (payload[5] >> 16) & 0xFF;
		     printf("   %sFormat ver%s=0x%X %sTTS%s=0x%X %sHTRstatus%s=0x%X %sDCCstatus%s=0x%X %sDCCrev%s=0x%x\n",
			bold,none,formver,bold,none,TTS,bold,none,HTRstatus,bold,none,DCCstatus,bold,none,DCCrev);
		    //
		    // setup for looking at spigots for this FED
		    //
		    printf("   Note: HTRs status: ([14:0] 1 bit per HTR from E*P*B*V*not(T), see below)\n");
		    cout << "  Now for the HTRs (spigots)..." << endl;
		    setup_spigots(true);
		}
	     }
	  }
	  else if (iw == 2) {
	    bool temp;
	    if (redirect) {
	      //
	      // close file
	      //
	      fclose(fout);
	      temp = false;
	    }
	    else {
	      //
	      // open file
	      //
	      cout << "Filename: (will overwrite any existing file) ";
	      char fname[128];
	      cin >> fname;
	      fout = fopen(fname,"w");
	      cout << endl;
	      temp = true;
	    }
	    redirect = temp;
	  }
	  else if (iw == 3) {
	    searching = true;
	    criteria = 0;
	    cout << "Event number: ";
	    scanf("%d",&find_evn);
	    return;
	  }
	  else if (iw == 4) {
	    //
	    // unformatted hex dump of the entire DCC payload (all spigots)
	    //
	    if (isFEDopen) {
		    int ndcc32 = size/4;   //size is in bytes
		    if (debugit) cout << "# DCC words is " << ndcc32 << endl;
		    for (int i=0; i<ndcc32; i++) {
			if (i == 0) {
			    printf("DCC Header:\n");
			    if (redirect) fprintf(fout,"DCC Header:\n");
			}
			else if (i == 6) {
			    printf("HTR payload summaries:\n");
			    if (redirect) fprintf(fout,"HTR payload summaries:\n");
			}
			for (int j=0; j<spigots; j++) if (i == ospigot[j]) { 
			    printf("HTR payload %d:\n",j+1);
			    if (redirect) fprintf(fout,"HTR payload %d:\n",j+1);
			}
			printf("  %4d 0x%8.8X\n",i,payload[i]);
			if (redirect) fprintf(fout,"  %4d 0x%8.8X\n",i,payload[i]);
		    }
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (iw == 5) {
	    //
	    // dump out the header information for this event
	    //
	    if (isFEDopen) {
	        cout << " There are " << spigots << " spigots" << endl;
	        cout << "  Spigot  #DCCWords  HTRerr  LRBerr  E-P-B-V-T " << endl;
		for (int i=0; i<spigots; i++) {
	          uint32_t s = hspigot[i];
	          int nwords,htrerr,lrberr,ee,ep,eb,ev,et;
	          payload_return(s,&nwords,&htrerr,&lrberr,&ee,&ep,&eb,&ev,&et);
	          if (nwords > 0) {
		    printf("  %4d       %5d    0x%2.2x    0x%2.2x   %d-%d-%d-%d-%d\n",
	       	       i,nwords,htrerr,lrberr,ee,ep,eb,ev,et);
		    if (redirect) fprintf(fout,"  %4d       %5d    0x%2.2x    0x%2.2x   %d-%d-%d-%d-%d\n",
	       	       i,nwords,htrerr,lrberr,ee,ep,eb,ev,et);
		  }
	        }
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (iw == 6) {
	    //
	    // formatted dump of the headers for all HTRs (spigots)
	    //
	    cout << "This FED (1) or all FEDs (0) :";
	    int kfed;
	    cin >> kfed;
	    if (kfed == 0 || (kfed > 0 && isFEDopen)) {
		if (kfed == 0) {
	    	  for (int ifed=0; ifed<nFEDs; ifed++) {
		      const FEDRawData& data = rawdata->FEDData(iFEDs[ifed]);
		      size=data.size();
	              FEDHeader header(data.data());
		      FEDTrailer trailer(data.data()+size-8);
		      cout << "====================================================================" << endl;
		      cout << "FED " << iFEDs[ifed] << " L1Id: " << setw(8) << header.lvl1ID();
		      cout << " BXId: " << setw(4) << header.bxID() << endl;
		      payload=(uint32_t*)(data.data());
		      setup_spigots(true);
		      for (int j=0; j<spigots; j++) {
		        htr_data_from_payload(payload,j);
		        if (j == 0) print_htr_payload_headers(true,j,false);
		        else print_htr_payload_headers(false,j,false);
		        htr_data_delete();
		      }
		  }
		}
		else {
		  //
		  // just do the FED that has been opened
		  //
		  cout << "Printout for " << spigots << " spigots of FED " << whichFED << endl;
	          for (int i=0; i<spigots; i++) {
		    //
		    // extract the header for this spigot
		    //
		    htr_data_from_payload(payload,i);
		    //
		    // print it out
		    //
		    if (debugit) cout << "ok, printout for spigot " << i << endl;
		    if (i == 0) print_htr_payload_headers(true,i,false);
		    else print_htr_payload_headers(false,i,false);
		    //
		    // now delete the spigot
		    //
		    htr_data_delete();
		  }
		}
	        //
		// just to be helpful, printout info about what all those things mean
		//
		printf("   %sbit0/OW%s=overflow warning goes away if L1A rate reduced\n",bold,none);
		printf("   %sbit1/BZ%s=internal buffers busy\n",bold,none);
		printf("   %sbit2/EE%s=empty event (because of previous BZ, includes only 1st 5 and last 3 words)\n",
			bold,none);
		printf("   %sbit3/RL%s=trig rule violation (rules 1 adn 2 of TDR 16.4.3)\n",bold,none);
		printf("   %sbit4/FE%s=detects FE idle words or CCA corrupted words that should be suppressed\n",
			bold,none);
		printf("   %sbit5/LW%s=latency warning, Jose fifo has tuned latency since previous Tx\n",
			bold,none);
		printf("   %sbit6/OD%s=so called fiber data error, is OR of 9,10,11 below\n",bold,none);
		printf("   %sbit7/CK%s=clocking error, ~TTCready || DLL_unlock\n",bold,none);
		printf("   %sbit8/BE%s=Bunch count error if BCC does not wrap nicely when BC0 received\n",bold,none);
		printf("   %sbit9/LK%s=OR of 8 fibers, each is ER || ~DV || Frame_error\n",bold,none);
		printf("   %sbit10/CE%s=OR of 8 fibers, each is capid rotation error\n",bold,none);
		printf("   %sbit11/FK%s=OR of 8 fibers, each is if FE flags are not right\n",bold,none);
		printf("   %sbit12/TM%s=test mode (0=data, 1=pattern or counter etc)\n",bold,none);
		printf("   %sbit13/HM%s=histo mode (if 1)\n",bold,none);
		printf("   %sbit13/CT%s=calib trig (if 1, otherwise L1A event)\n",bold,none);
		printf("   %sDLL%s  counts number of times DLL unlock since last reset (trust not)\n",bold,none);
		printf("   %sTTC%s  TTCready, should always be asserted\n",bold,none);
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (iw == 7) {
	    //
	    // dump payload for this spigot, but prompt for which one first
	    //
	    if (isFEDopen) {
		int spign;
		cout << "Which spigot (0-" << spigots-1 << "): ";
		scanf("%d",&spign);
		//
		// extract the data for this spigot
		//
		if (debugit) cout << "htr_data_from_payload...";
		htr_data_from_payload(payload,spign);
		if (debugit) cout << "done" << endl;
		//
		// printout the header info
		//
		if (debugit) cout << "printing htr payload headers...";
		print_htr_payload_headers(true,spign,true);
		if (debugit) cout << "done" << endl;
		//
		// extract the TPGs
		//
		if (debugit) cout << "tpgs_from_htr...";
		tpgs_from_htr();
		if (debugit) cout << "done" << endl;
		//
		// printout the TPs
		//
		if (debugit) cout << "print_htr_tpgs...";
		print_htr_tpgs();
		if (debugit) cout << "done" << endl;
		//
		// now delete the tpgs
		//
		if (debugit) cout << "tpgs_delete...";
		tpgs_delete();
		if (debugit) cout << "done" << endl;
		//
		// now look at the raw data
		//
		if (debugit) cout << "qies_from_htr...";
		qies_from_htr();
		if (debugit) cout << "done" << endl;
		//
		// now do a formatted dump
		//
		if (debugit) cout << "print_htr_qies...";
		print_htr_qies();
		if (debugit) cout << "done" << endl;
		//
		// now delete the qies
		//
		if (debugit) cout << "qies_delete...";
		qies_delete();
		if (debugit) cout << "done" << endl;
		//
		// and delete the entire payload for this HTR
		//
		if (debugit) cout << "htr_data_delete...";
		htr_data_delete();
		if (debugit) cout << "done" << endl;
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (iw == 8) {
	    //
	    // dump payload for this spigot, but prompt for which one first
	    //
	    if (isFEDopen) {
		int spign;
		cout << "Which spigot (0-" << spigots-1 << "): ";
		scanf("%d",&spign);
		//
		// extract the data for this spigot
		//
		htr_data_from_payload(payload,spign);
		//
		// printout in unformatted (or close to it anyway)
		//
		htr_data_print();
		//
		// and delete the entire payload for this HTR
		//
		htr_data_delete();
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (iw == 9) {
	    cout << "===========================================";
	    cout << "===========================================" << endl;
	    cout << "  1  Search for non-zero TPG" << endl;
	    cout << "  2  Search for self-similar QIEs" << endl;
	    cout << "  3  Search for event with error(s) in HTR header" << endl;
	    cout << "  4  Printout BCN, ORN, Evn, BCNidle, etc" << endl;
	    cout << "  5  Report on any HTR that has any error bits set in the header" << endl;
	    cout << "What is your pleasure?> ";
	    scanf("%d",&criteria);
	    if (criteria > 0 && criteria < 6) {
	      searching = true;
	      printbegin = false;
	      if (criteria == 4) {
		int iw;
		cout << "But not 12 and 13? (1=right, 0=use both): ";
		cin >> iw;
		if (iw == 0) not12_13 = false;
		else not12_13 = true;
	      }
	      else if (criteria == 5) {
		cout << "How many events to loop over? :";
		cin >> nloopse;
		ncountse = 0;
	      }
	      return;
	    }
	    else cout << "Hey, follow instructions, read the menu!" << endl;
	  }
	  else cout << "Hmm, can't follow instructions?  Try again" << endl;
	    
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
RawAnalyzer::beginJob()
{
  cout << "In the beginning....." << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RawAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(RawAnalyzer);
