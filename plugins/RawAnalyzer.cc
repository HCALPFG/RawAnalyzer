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
#include <bitset>
#include <fstream>
#include <stdarg.h>
#include <readline/readline.h>
#include <readline/history.h>
#include "TH1.h"
#include "TH2.h"


using namespace edm;
using namespace std;

#define bold  "\33[0;1m"
#define red   "\033[0;31m"        // 0 -> normal ;  31 -> red 
#define cyan  "\033[1;36m"        // 1 -> bold ;  36 -> cyan
#define green "\033[4;32m"        // 4 -> underline ;  32 -> green
#define blue  "\033[9;34m"        // 9 -> strike ;  34 -> blue 
#define black  "\033[0;30m"
#define brown  "\033[0;33m"
#define magenta  "\033[0;35m"
#define gray  "\033[0;37m"
#define none   "\033[0m"        // to flush the previous property

const char* CodeVersion="08AUG2014-v0";

int nargs = 0;
char *args[100]; //max of 10 arguments on the stack
//
// HTR types, starting from 0:  HO,HBHE1+,HBHE2+,HBHE3+,HBHE4+,HBHE5+,HBHE6+,HBHE7+,
//                              HF,HBHE1-,HBHE2-,HBHE3-,HBHE4-,HBHE5-,HBHE6-,HBHE7-
const uint32_t ntptype[19] = {0,24,16,16,16,16,16,16,2,24,16,16,16,16,16,16,0,0,0};
const char* htr_types[19] = {"HO    ","HBHE1+","HBHE2+","HBHE3+","HBHE4+","HBHE5+","HBHE6+","HBHE7+",
                             "HF    ","HBHE1-","HBHE2-","HBHE3-","HBHE4-","HBHE5-","HBHE6-","HBHE7-",
			     "CALIB ","DEBUG ","ZDC   "};
const uint32_t nslbtype[19] = {0,6,4,4,4,4,4,4,1,6,4,4,4,4,4,4,0,0,0}; // ntptype/4
const char* err_types[15] = {"CT","HM","TM","FK","CE","LK","BE","CK","OD","LW","FE","RL","EE","BZ","OW"};
bool searching = false;
bool isFEDopen = false;
bool printbegin = true;
int nloopse = 0;
int ncountse = 0;
int qieMANTcut = 0;
int qieEXPcut = 0;
int find_evn = 0;
int evn, sr, evn1, evn2, orn, htrn, bcn, formatVer, npresamp, ndll, ntps, fw, subV, pipeline, htype;
int tpsamps, qiesamps, ndccw, nwc, reserved, us, cm;
int ow, bz, ee, rl, le, lw, od, ck, be, tmb, hm, ct, odle, odce, odfe, isone, ttcready, dll;
int hspigot[30];
int wspigot[30];
int ospigot[30];
int spigots = 0;
FILE* fout = NULL;
uint32_t* hdata;  //here is the 16-bit HTR payload data
uint32_t n16;     //number of 16-bit words for this HTR payload
uint32_t* tpgs;
//
// nominal qie to fC....close enough I think
int qie2fC[4][32]={
{-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,  //0-15
17,19,21,23,25,27,29,32,35,38,42,46,50,54,59,64}, //16-31
{59,64,69,74,79,84,89,94,99,104,109,114,119,124,129,137,  //32-47
147,157,167,177,187,197,209,224,239,254,272,292,312,334,359,384}, //48-63
{359,384,409,434,459,484,509,534,559,584,609,634,659,684,709,747, //64-79
797,847,897,947,997,1047,1109,1184,1259,1334,1422,1522,1622,1734,1859,1984}, //80-95
{1859,1984,2109,2234,2359,2484,2609,2734,2859,2984,3109,3234,3359,3484,3609,3797, //96-101
4047,4297,4547,4797,5047,5297,5609,5984,6359,6734,7172,7672,8172,8734,9359,9984}}; //102-127
// qies, 240 max (24 x 10 time samples), with 0=mant,1=range,2=fiber,3=qie,4=capid,5=error
#define MANT 0
#define RANGE 1
#define FIBER 2
#define QIE 3
#define CAPID 4
#define ERROR 5
int nqies;
int qies[240][6];
// qiestats[FED-700(100)][SPIGOT(16)][FIBER(8)][QIE(3)][CAPID(4)][0=number,1=sum(x),2=sum(x^2)]
int qiestats[100][16][8][3][4][3];
int qiestatsAny[100][16];
int qieFlavor = 0;
uint32_t mant,range,capid,dv,er,addr,fib;
int nFEDs = 0;
int iFEDs[100];  // there should not be that many, fewer than 20 probably but you never know
int iFEDsize[100];
const uint32_t* payload;
bool debugit = false;
int whichFED = -1;
size_t size = 0;
int criteria = 0;
int fed1, fed2;
bool not12_13 = true;
Handle<FEDRawDataCollection> rawdata;

void payload_spigot_header_return(uint32_t s,
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
void print_htr_qies(int spigot);
void qies_unpack(int which);
void qies_unpack_word(uint32_t data);
void htr_data_print();
void htr_data_print_formatted(int spign);
void setup_spigots(bool printout);
void find_htr_header_errors();
void check_event_numbers(int irun, int iev);
void getFeds(int* nfeds, int* ifeds, bool printout, int irun, int iev);
bool checkFedBcN(int nfeds, int* ifeds, int* fedBcN);
bool checkFedEvN(int nfeds, int* ifeds, int* fedEvN);
bool checkFedOrN(int nfeds, int* ifeds, int* fedOrN);
int get_spigot_bcn(int ispigot);
int get_spigot_evn(int ispigot);
int get_spigot_orn(int ispigot);
bool checkFedBcNIdle(int nfeds, int* ifeds, int* fedBcNIdle);
void check_htr_header_errors(int irun,int iev);
char* vparse_input(const char* prompt, int num, ...);
int getINT(const char* prompt);
void moveup();
int getHEX(const char* prompt);
void print_error_warnings();
void check_qie(int type, int irun, int iev);
void qies_stats_print();

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

void setup_spigots(bool printout) {
	//
	// this takes the DCC (FED) payload and sets up pointers to the individual spigots
	// note: there are only 12 spigots max for the DCC2 (as of July 2014)
	//
	// this has to be called as soon as you decide on which FED you want to look at
	//
	int nwords,htrerr,lrberr,ee,ep,eb,ev,et;
	spigots = 0;
	if (printout) cout << "  Spigot  #Words  HTRerr  LRBerr  E-P-B-V-T " << endl;
	for (int i=0; i<18; i++) {
	       uint32_t s = payload[i+6] & 0xFFFFFFFF;
	       payload_spigot_header_return(s,&nwords,&htrerr,&lrberr,&ee,&ep,&eb,&ev,&et);
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

void htr_data_from_payload(const uint32_t *payload, int spigot) {
	//
	// start with a particular spigot, then take 32-bit data from "payload"
	// and create 16-bit "hdata" words that holds the data for that spigot
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
	//
	// delete the data collected in htr_data_from_payload
	delete [] hdata;
}

void payload_spigot_header_return(uint32_t s,
 	int* nwords,int* htrerr,int* lrberr,int* ee,int* ep,int* eb,int* ev,int* et) {
	//
	// input should be a pointer to the n'th spigot of the DCC header, it tells you the
	// number of words for the spigot, and some error bits
	//
	// this is used by setup_spigots
	//
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
	//
	// this function digs out various spigot header words and puts it into local variables
	// for ease of access
	//
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
//	if (htype > 15) cout << "  HTYPE is " << htype << " which is a surprise and probably an error!!!!" << endl;
//	else cout << "  This payload appears to be from the " << ntptype[htype]	<< " subdetector...." << endl;
	pipeline = hdata[7]&0xFF;
	if (debugit) cout << "htype " << htype;
	if ( (htype == 0) || (htype > 15) ) {
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

void htr_data_print_formatted(int spign) {
	//
	// first the headers
	//
	htr_fill_header();
	print_htr_payload_headers(true,spign,true);
	//
	// next the TPs
	//
	tpgs_from_htr();
	print_htr_tpgs();
	tpgs_delete();
	//
	// then the QIEs
	//
	//
	// next come the QIE data, time to figure out which flavor.  5=old, 6="new"
	//
	int nqiewords = 0;
	cout << "QIE values:" << endl;
	if (formatVer == 5) {
	  cout << "    Format version=5, no packing, also obsolete as of long ago (it is currently 2014 at least)" << endl;
	}
	else if (formatVer == 6) {
	  cout << "Format version=6, QIE packing in blocks, depends on flavor which is found in 1st word of each block" << endl;
	  int flavor = 0;
	  int jpt = 8 + ntps;
	  nqiewords = n16 - jpt - 12;  // 12 for the extra(8), pre-trailer(3), and trailer(1) 
				// the data is packed so this is the number of words, NOT #qies!
	  //
	  // need to know which flavor.  peek into 1st word, it should tell us
	  //
	  int if1 = hdata[jpt];
	  int ib1 = (if1>>15)&0x1;
	  int fl1 = (if1>>12)&0x7;
	  if (fl1 == 5) {
	    cout << "   Flavor  Errf0/1 Capid Fib QIE ---- QIE samples, mant/range ---- (SOI in bold)" << endl;
	  }
	  else if (fl1 == 6) {
	    cout << "   Flavor  Errf0/1 Capid Fib QIE (SOI in bold)" << endl;
	  }
	  else {
		cout << "What?  Flavor " << if1 << " is illegal, something is seriously wrong" << endl;
		return;
	  }
	  if (ib1 == 0) {
		cout << "What?  MSB is not 1, something is seriously wrong" << endl;
		return;
	  }
	  bool errf0OK = true;
	  bool errf1OK = true;
	  bool capidOK = true;
	  bool allfibers = true;
	  int capid1 = -1;
	  int fib1 = -1;
	  int isamp = 0;
	  for (int i=0; i<nqiewords; i++) {
		int iword = hdata[jpt+i];
		int bit31 = (iword>>15)&1;
		if (bit31 == 1) {
		    isamp = 0;
		    if (i>0) cout << endl;
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
		    printf("      %d        %d/%d   %d    %d   %d",flavor,errf0,errf1,capid,thisfiber,thisqie);
		    //
		    // perform some basic checks
		    //
		    if (errf0 == 1) errf0OK = false;
		    if (errf1 == 1) errf1OK = false;
		    // check that all CAPIDs are the same (they should be)
		    if (i == 0) capid1 = capid;
		    else {
			if (capid != capid1) capidOK = false;
			capid1 = capid;
		    }
		    // check if all fibers are there, assume incrementing
		    if ( (thisqie == 0) && (thisfiber != (fib1+1)) ) allfibers = false;
		    fib1 = thisfiber;
		}
		else {
		    if (flavor == 5) {
			int qie1 = (iword>>8)&0x7F;
			int r1 = (qie1>>6)&0x3;
			int m1 = qie1&0x1F;
			int qie0 = iword&0x7F;
			int r0 = (qie0>>6)&0x3;
			int m0 = qie0&0x1F;
			if (isamp == npresamp) printf("%s %2d/%d%s",bold,m0,r0,none);
			else printf(" %2d/%d",m0,r0);
			isamp++;
			if (isamp == npresamp) printf("%s %2d/%d%s",bold,m1,r1,none);
			else printf(" %2d/%d",m1,r1);
			isamp++;
		    }
		    else if (flavor == 6) {
			int qie = (iword>>8)&0x7F;
			int r = (qie>>6)&0x3;
			int m = qie&0x1F;
			int capid = (iword>>8)&0x3;
			int dv = (iword>>10)&0x1;
			int er = (iword>>11)&0x1;
			if (isamp == npresamp) printf(
			"%s   Qie(mant/range)=%2d/%d  DV=%d  ER=%d  Capid=%d%s\n",
				bold,m,r,dv,er,capid,none);
			else printf("   Qie(mant/range)=%2d/%d  DV=%d  ER=%d  Capid=%d\n",m,r,dv,er,capid);
			isamp++;
		    }
		    else cout << "FLAVOR " << flavor << " IS UNKNOWN AND CONFUSES ME TERRIBLY!!!" << endl;
		}
	  }
	  //
	  // report on the consistency checks (always done for the flavor header)
	  //
	  cout << endl << "  ===> Error checking...";
	  if (errf0OK) cout << "no CAPID errors...";
	  else cout << "CAPID ERRORS!!!! (errf0)...";
	  if (errf1OK) cout << "no LINK errors...";
	  else cout << "LINK ERRORS!!!! (errf1)...";
	  if (capidOK) cout << "all CAPIDs match...";
	  else cout << "CAPID mismatch!!!...";
	  if (allfibers) cout << "all fibers present" << endl;
	  else cout << "NOT all fibers present!!!!" << endl;
	}
	//
	// then the parity word (if 0xFFFF)
	//
	int jpt = 8 + ntps + nqiewords;
	if (hdata[jpt] == 0xFFFF) {
	  printf("Parity word 0xFFFF detected (put in to make an even number of payload words sent to the DAQ)\n");
	  jpt++;
	}
	//
	// and then the 8 "extra" words
	//
	cout << "Next come 8 'Extra-info' words" << endl;
	if (cm == 0) {
	  if (us == 0) {
		printf("   CM=US=0 Normal mode data\n");
		printf("   Fiber Empty Full Count    BCN idle\n");
		bool bcnidleOK = true;
		int bcn1 = -1;
		for (int i=0; i<8; i++) {
		  int iword = hdata[jpt+i];
		  int empty = (iword>>15)&1;
		  int full = (iword>>14)&1;
		  int latcnt = (iword>>12)&0x3;
		  int bcnidle = iword&0xFFF;
		  printf("      %d    %d     %d     %d   %d (0x%X)\n",i+1,empty,full,latcnt,bcnidle,bcnidle);
		  if (i == 0) bcn1 = bcnidle;
		  else {
		    if (bcnidle != bcn1) bcnidleOK = false;
		    bcn1 = bcnidle;
		  }
		}
		if (bcnidleOK) cout << "  ===> All BCN idle counters match..." << endl;
		else cout << "  ===> Problem with BCN idles - at least one of them does not match!!!!" << endl;
	  }
	  else {
	   	printf("   CM=0 US=1 Unsuppressed data detected\n");
		// 1
		int iword = hdata[jpt];
		int maskDigi8to1 = iword&0xFF;
		int maskTP8to1 = (iword>>8)&0xFF;
		printf("   Mask TP  8- 1: 0x%3.3X  Mask Digi  8- 1: 0x%3.3X\n",maskTP8to1,maskDigi8to1);
		// 2
		iword = hdata[jpt+1];
		int maskDigi16to9 = iword&0xFF;
		int maskTP16to9 = (iword>>8)&0xFF;
		printf("   Mask TP 16- 9: 0x%3.3X  Mask Digi 16- 9: 0x%3.3X\n",maskTP16to9,maskDigi16to9);
		// 3
		iword = hdata[jpt+2];
		int maskDigi24to17 = iword&0xFF;
		int maskTP24to17 = (iword>>8)&0xFF;
		printf("   Mask TP 24-17: 0x%3.3X  Mask Digi  24-17: 0x%3.3X\n",maskTP24to17,maskDigi24to17);
		// 4
		iword = hdata[jpt+3];
		int threshDigi1 = iword&0xFF;
		int threshDigi24 = (iword>>8)&0xFF;
		printf("   Threshold Digi24: 0x%3.3X  Threshold Digi1: 0x%3.3X\n",threshDigi24,threshDigi1);
		// 5
		iword = hdata[jpt+4];
		int syncFcnt1 = iword&0x3;
		int syncFcnt2 = (iword>>2)&0x3;
		int syncFcnt3 = (iword>>4)&0x3;
		int syncFcnt4 = (iword>>6)&0x3;
		int threshTP3to0 = (iword>>12)&0xF;
		printf("   Threshold TP[3:0]: 0x%3.3X  Sycn Count Fibers: 4=%d 3=%d 2=%d 1=%d\n",
			threshTP3to0,syncFcnt4,syncFcnt3,syncFcnt2,syncFcnt1);
		// 6
		iword = hdata[jpt+5];
		int syncFcnt5 = iword&0x3;
		int syncFcnt6 = (iword>>2)&0x3;
		int syncFcnt7 = (iword>>4)&0x3;
		int syncFcnt8 = (iword>>6)&0x3;
		int threshTP7to4 = (iword>>12)&0xF;
		printf("   Threshold TP[7:4]: 0x%3.3X  Sycn Count Fibers: 8=%d 7=%d 6=%d 5=%d\n",
			threshTP7to4,syncFcnt8,syncFcnt7,syncFcnt6,syncFcnt5);
		// 7
		iword = hdata[jpt+6];
		int bcnofrxbc0 = iword&0xFFF;
		int ZSmask18to16 = (iword>>12)&0x7;
		printf("   ZS_Mask[18:16]: 0x%2.2X  BCN-of-RxBC0 = %d\n",ZSmask18to16,bcnofrxbc0);
		// 8
		iword = hdata[jpt+7];
		int ZSmask15to0 = iword;
		printf("   ZS_Mask[15:0]: 0x%4.4X\n",ZSmask15to0);
	  }
	}
	else {
		printf("   CM=1 Compact mode, this means that there are no extra words here\n");
	}
	//
	// then 4 trailer words
	//
	cout << "Next come 3 'Pre-trailer' words followed by the 'Trailer'" << endl;
	jpt = jpt + 8;
	printf("   Pre-Trailer 3: DAQsamples[4:0]=%d   DAQwords[10:0]=%d\n",
		(hdata[jpt]>>11)&0x1F,hdata[jpt]&0x3FF);
        printf("   Pre-Trailer 2: CRC (or perhaps replaced with word count?) = 0x%X\n",hdata[jpt+1]);
	printf("   Pre-Trailer 1: Overwritten to 0x%X by the DCC\n",hdata[jpt+2]);
	printf("   Trailer:       EVN[7:0]=0x%X (low byte is 0x%X, possibly overwritten by DCC)\n",(hdata[jpt+3]>>8)&0xFF,hdata[jpt+3]&0xFF);
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
	printf("Header 5: 0x%4.4X   {Version[3:0]=%d,BCN[11:0]=%d}\n",hdata[4],formatVer,bcn);
	printf("Header 6: 0x%4.4X   {#TPs[7:0]=%d,#Presamples[4:0]=%d,DLL[1:0]=%d,TTCready=%d}\n",
		hdata[5],ntps,npresamp,ndll,ttcready);
	printf("Header 7: 0x%4.4X   {US=%d,CM=%d,Reserved=%d,SubV[3:0]=%d,fw[7:0]=0x%X}\n",
			hdata[6],us,cm,reserved,subV,fw);
	printf("Header 8: 0x%4.4X   {HCALtype[7:0]=%d (is %s), Pipelength[7:0]=%d}\n",
		hdata[7],htype,htr_types[htype],pipeline);
	int ipt = 8;
	//
	// next come the TPGs
	//
	int format_version = hdata[4] >> 12;
	if (format_version == 5) {
		printf(
		"%3d TPG channels, %2d time samples per channel, format {SLBid[2:0],SLBch[1:0],Z,SOI,TP[8:0]}\n",
			ntptype[htype],tpsamps);
		int jpt = 0;
		for (int i=0; i<(int) ntptype[htype]; i++) {
		  printf("%2d: ",jpt+1);
		  jpt++;
		  for (int j=0; j<tpsamps; j++) {
		    printf("0x%4.4X={%d,%d,%d,%d,0x%3.3X} ",hdata[ipt],hdata[ipt]>>13,
		        (hdata[ipt]>>11)&3,(hdata[ipt]>>10)&1,(hdata[ipt]>>9)&1,hdata[ipt]&0x1FF); 
		    ipt++;
		  }
		  printf("\n");
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
			int jpt = 8;
			for (int i=0; i<ntps; i++) {
				int tword = hdata[jpt+i];
				int high5 = (tword>>11)&0x1F;
				int z = (tword>>10)&0x1;
				int soi = (tword>>9)&0x1;
				int mbits = tword&0xFF;
				printf("0x%4.4X={high5=0x%2.2X,z=%d,soi=%d,muon_bits=0x%2.2X}\n",
					tword,high5,z,soi,mbits);				
			}
		}
		else {
			printf(
			"%3d TPG channels present, format {SLBid[2:0],SLBch[1:0],Z,SOI,TP[8:0]}\n",	ntps);
			int jpt = 8;
			for (int i=0; i<ntps; i++) {
				int tword = hdata[jpt+i];
				int tp_tag = (tword>>11)&0x1F;
				int slb_id = (tp_tag >> 2)&0x7;
				int slb_ch = tp_tag & 0x3;
				int z = (tword>>10)&0x1;
				int soi = (tword>>9)&0x1;
				int tp = tword&0x1FF;
				printf("0x%4.4X={slb_id=%d,slb_ch=%d,z=%d,soi=%d,tp=%d}\n",
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
/*		int jpt = 0;
		nqies = n16 - ntps - 20;
		int nqiech = nqies/qiesamps;
		printf("%3d QIE channels, %2d time samples per channel\n ",nqiech,qiesamps);
//		printf("format {Fiber[2:0],QIE[1:0],ER,DV,Cap[1:0],Range[1:0],Mant[4:0]}:\n");
		printf("   fi    QIE_CH       CAPID       ER          DV     Mant/Range\n");
		for (int i=0; i<nqiech; i++) {
		  printf("%2d: ",jpt+1);
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
		}*/
		cout << "HTR has format version 5 which is obsolete!!!!" << endl;
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
		  }
		  else {
		    if (flavor == 5) {
			int qie1 = (iword>>8)&0x7F;
			int r1 = (qie1>>6)&0x3;
			int m1 = qie1&0x1F;
			int qie0 = iword&0x7F;
			int r0 = (qie0>>6)&0x3;
			int m0 = qie0&0x1F;
			printf("0x%4.4X  Qie1(mant/range)=%2d/%d  Qie0(mant/range)=%2d/%d\n",iword,m1,r1,m0,r0);
		    }
		    else if (flavor == 6) {
			int qie = (iword>>8)&0x7F;
			int r = (qie>>6)&0x3;
			int m = qie&0x1F;
			int capid = (iword>>8)&0x3;
			int dv = (iword>>10)&0x1;
			int er = (iword>>11)&0x1;
			printf("0x%4.4X  Qie(mant/range)=%2d/%d  DV=%d  ER=%d  Capid=%d\n",iword,m,r,dv,er,capid);
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
	  jpt++;
	}
	//
	// and then the 8 "extra" words
	//
	if (cm == 0) {
	  if (us == 0) {
		printf("CM=US=0 Normal mode data\n");
		for (int i=0; i<8; i++) {
		  int iword = hdata[jpt+i];
		  printf("Extra %d: 0x%4.4X  {Empty=%d,Full=%d,LatCnt[1:0]=%d,BCNtime[11:0]=%d\n",
		  	i+1,iword,(iword>>15)&1,(iword>>14)&1,(iword>>12)&0x3,iword&0xFFF);
		}
	  }
	  else {
	   	printf("CM=0 US=1 Unsuppressed data detected\n");
		int iword = hdata[jpt];
		int maskDigi8to1 = iword&0xFF;
		int maskTP8to1 = (iword>>8)&0xFF;
		printf("Extra 1: 0x%4.4X  Mask TP  8- 1: 0x%3.3X  Mask Digi  8- 1: 0x%3.3X\n",
			iword,maskTP8to1,maskDigi8to1);
		iword = hdata[jpt+1];
		int maskDigi16to9 = iword&0xFF;
		int maskTP16to9 = (iword>>8)&0xFF;
		printf("Extra 2: 0x%4.4X  Mask TP 16- 9: 0x%3.3X  Mask Digi 16- 9: 0x%3.3X\n",
			iword,maskTP16to9,maskDigi16to9);
		iword = hdata[jpt+2];
		int maskDigi24to17 = iword&0xFF;
		int maskTP24to17 = (iword>>8)&0xFF;
		printf("Extra 3: 0x%4.4X  Mask TP 24-17: 0x%3.3X  Mask Digi  24-17: 0x%3.3X\n",
			iword,maskTP24to17,maskDigi24to17);
		iword = hdata[jpt+3];
		int threshDigi1 = iword&0xFF;
		int threshDigi24 = (iword>>8)&0xFF;
		printf("Extra 4: 0x%4.4X  Threshold Digi24: 0x%3.3X  Threshold Digi1: 0x%3.3X\n",
			iword,threshDigi24,threshDigi1);
		iword = hdata[jpt+4];
		int syncFcnt1 = iword&0x3;
		int syncFcnt2 = (iword>>2)&0x3;
		int syncFcnt3 = (iword>>4)&0x3;
		int syncFcnt4 = (iword>>6)&0x3;
		int threshTP3to0 = (iword>>12)&0xF;
		printf("Extra 5: 0x%4.4X  Threshold TP[3:0]: 0x%3.3X  Sycn Count Fibers: 4=%d 3=%d 2=%d 1=%d\n",
			iword,threshTP3to0,syncFcnt4,syncFcnt3,syncFcnt2,syncFcnt1);
		iword = hdata[jpt+5];
		int syncFcnt5 = iword&0x3;
		int syncFcnt6 = (iword>>2)&0x3;
		int syncFcnt7 = (iword>>4)&0x3;
		int syncFcnt8 = (iword>>6)&0x3;
		int threshTP7to4 = (iword>>12)&0xF;
		printf("Extra 6: 0x%4.4X  Threshold TP[7:4]: 0x%3.3X  Sycn Count Fibers: 8=%d 7=%d 6=%d 5=%d\n",
			iword,threshTP7to4,syncFcnt8,syncFcnt7,syncFcnt6,syncFcnt5);
		iword = hdata[jpt+6];
		int bcnofrxbc0 = iword&0xFFF;
		int ZSmask18to16 = (iword>>12)&0x7;
		printf("Extra 7: 0x%4.4X  ZS_Mask[18:16]: 0x%2.2X  BCN-of-RxBC0 = %d\n",
			iword,ZSmask18to16,bcnofrxbc0);
		iword = hdata[jpt+7];
		int ZSmask15to0 = iword;
		printf("Extra 8: 0x%4.4X  ZS_Mask[15:0]: 0x%4.4X\n",iword,ZSmask15to0);
	  }
	}
	else {
		printf("CM=1 Compact mode, no extra words here\n");
	}
	//
	// then 4 trailer words
	//
	jpt = jpt + 8;
	printf("Pre-Trailer 3: 0x%4.4X  {DAQsamples[4:0]=%d,DAQwords[10:0]=%d}\n",
		hdata[jpt],hdata[jpt]>>11,hdata[jpt]&0x3FF);
        printf("Pre-Trailer 2: 0x%4.4X  {CRC}\n",hdata[jpt+1]);
	printf("Pre-Trailer 1: 0x%4.4X  {all zeros or possibly overwritten by DCC}\n",hdata[jpt+2]);
	printf("Trailer:       0x%4.4X  {EVN[7:0]=0x%X,low byte=0x%X}\n",hdata[jpt+3],(hdata[jpt+3]>>8)&0xFF,hdata[jpt+3]&0xFF);

}

void print_htr_payload_headers(bool header, int spigot, bool extra) {
	if (header) {
          printf("                               ");
	  printf("Samples               CHTFCLBCOLFREBO\n");
          printf(" Spigot    EvN    BcN   OrN HTR  fw   #TP");
	  printf("      Htype      TP  QIE Presamp DCCw  TMMKEKEKDWELEZW TTC DLL\n");
        }
	if (debugit) cout << "calling htr_fill_header" << endl;
	htr_fill_header();
	//
	// the following are all errors if == 1 except for the last one
	//
	if (debugit) cout << "**** nwc *** " << nwc << "  "  << hex << nwc << dec << endl;
	printf("   %2d    %4d   %5d  %4d %3d 0x%2.2X   %2d   0x%2.2x=%6s   %2d   %2d   %2d     %3d  ",
	    spigot,evn,bcn,orn,htrn,fw,ntps,htype,htr_types[htype],tpsamps,qiesamps,npresamp,ndccw);
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
	}
}

void tpgs_from_htr() {
	//
	// collect all of the TPs into arrays
	//
	tpgs = new uint32_t [ntps];
	for (int j=0; j<ntps; j++) {
		tpgs[j] = hdata[j+8];
		if (debugit) cout << hex << j << " 0x" << tpgs[j] << endl;
	}
}

void tpgs_delete() {
	delete [] tpgs;
}

void qies_from_htr() {
	//
	// collect all of the qies into arrays.  depends on formatVer, and flavor for formatVer=6
	// since it's now Aug 2014, I'm not going to bother making this work for formatVer=5.  
	//
	if (formatVer == 5) {
		cout << "Sorry, but HTR format version 5 is no longer supported (at least not by me)" << endl;
		return;
	}
	//
	// ok, now we know which spigot, and hdata contains the words for that spigot providing
	// htr_data_from_payload has already been called to fill hdata
	//
	// for the VME system, the most we can have is 24 channels times 10 time samples
	//
	nqies = 0;
	int jpt = ntps + 8;	// 8 header words, plus number of TPs
	int nqiewords = n16 - jpt - 12;  // 12 for the extra(8), pre-trailer(3), and trailer(1) 
	int flavor=0,errf0=0,errf1=0,capid=0,channelid=0,thisqie=0,thisfiber=0;
//	cout << "n16= " << n16 << " jpt= " << jpt << " nqiewords=" << nqiewords << endl;
	for (int i=0; i<nqiewords; i++) {
	  int iword = hdata[i+jpt];
//	  cout << "  i=" << i << " hdata = 0x" << hex << iword << dec << endl;	  
	  int bit31 = (iword>>15)&1;
	  if (bit31 == 1) {
	    //
	    // channel header
	    //
	    flavor = (iword>>12)&0x7;
	    qieFlavor = flavor;
	    errf0 = (iword>>10)&0x1;
	    errf1 = (iword>>11)&0x1;
	    capid = (iword>>8)&0x3;
	    channelid = iword&0xFF;
	    thisqie = channelid&0x3;
	    thisfiber = (channelid>>2)&0x7;
//	    printf(
//	       "0x%4.4X Header: Flavor=%d CapIdErr=%d  LinkErr=%d  1st Capid=%d Qie=%d Fiber=%d\n",
//			iword,flavor,errf0,errf1,capid,thisqie,thisfiber);
	  }
	  else {
	    if (flavor == 5) {
		//
		// compact mode.   increment capid by hand.  check ERROR word to be sure it's ok
		//
		int qie0 = iword&0x7F;
		qies[nqies][RANGE] = (qie0>>6)&0x3;
		qies[nqies][MANT] = qie0&0x1F;
		qies[nqies][FIBER] = thisfiber;
		qies[nqies][QIE] = thisqie;
		qies[nqies][CAPID] = capid++;	// increment!
		if (capid == 4) capid = 0;
		if (errf0 || errf1) {
			if (errf0) qies[nqies][ERROR] = 1;
			if (errf1) qies[nqies][ERROR] = qies[nqies][ERROR]&0x2;
		}
		else qies[nqies][ERROR] = 0;
		nqies++;
		// now get the upper half, which is the next qie
		int qie1 = (iword>>8)&0x7F;
		qies[nqies][RANGE] = (qie1>>6)&0x3;
		qies[nqies][MANT] = qie1&0x1F;
		qies[nqies][FIBER] = thisfiber;
		qies[nqies][QIE] = thisqie;
		qies[nqies][CAPID] = capid++;	// increment!
		if (capid == 4) capid = 0;
		if (errf0 || errf1) {
			if (errf0) qies[nqies][ERROR] = 1;
			if (errf1) qies[nqies][ERROR] = qies[nqies][ERROR]&0x2;
		}
		else qies[nqies][ERROR] = 0;
		nqies++;
	    }
	    else if (flavor == 6) {
		int qie = (iword>>8)&0x7F;
		int r = (qie>>6)&0x3;
		int m = qie&0x1F;
		int capid = (iword>>8)&0x3;
		int dv = (iword>>10)&0x1;
		int er = (iword>>11)&0x1;
		qies[nqies][RANGE] = r;
		qies[nqies][MANT] = m;
		qies[nqies][FIBER] = thisfiber;
		qies[nqies][QIE] = thisqie;
		qies[nqies][CAPID] = capid;	// increment!
		if (errf0 || errf1 || (er==1) || (dv=0) ) {
			if (errf0) qies[nqies][ERROR] = 1;
			if (errf1) qies[nqies][ERROR] = qies[nqies][ERROR]&0x2;
			if (er==1) qies[nqies][ERROR] = qies[nqies][ERROR]&0x4;
			if (dv==0) qies[nqies][ERROR] = qies[nqies][ERROR]&0x8;
		}
		else qies[nqies][ERROR] = 0;
		nqies++;
	    }
	    else {
		cout << "FLAVOR " << flavor << " IS UNKNOWN AND CONFUSES ME TERRIBLY!!!" << endl;
		return;
	    }
	  }
	}
}


/*
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
*/

void print_htr_qies(int spigot) {
	//
	// qiesamps is the variable that you use to know how many time samples
	//
	printf("\n\nBehold info for %d QIEs from %s spigot %d: %d time samples/QIE, %d HCAL towers. QIE block flavor is %d\n",
	 	nqies,htr_types[htype],spigot,qiesamps,nqies/qiesamps,qieFlavor);
	cout << " Fib QIE  ===> Samples, format is 'capid:mant/range(error)'  (SOI=bold)  See below about 'errors'." << endl;
	int old_qie = -1;
	int isamp = 0;
	for (int i=0; i<nqies; i++) {
	  int ifiber = qies[i][FIBER];
	  int iqie = qies[i][QIE];
	  int icap = qies[i][CAPID];
	  int mant = qies[i][MANT];
	  int range = qies[i][RANGE];
	  int ierr = qies[i][ERROR];
	  if (old_qie == iqie) {
	  	//
		// print out info for qie here
		//
		if (isamp == npresamp) {
		  if (ierr == 0) printf("%s%1d:%2.2d/%1d       %s",bold,icap,mant,range,none);
		  else printf("%s%1d:%2.2d/%1d(0x%X) %s",bold,icap,mant,range,ierr,none);
		}
		else {
		  if (ierr == 0) printf("%1d:%2.2d/%1d       ",icap,mant,range);
		  else printf("%1d:%2.2d/%1d(0x%X) ",icap,mant,range,ierr);
		}
		isamp++;
	  }
	  else {
		//
		// print fiber and qie number (SOI will never be the 1st one....)
		// 
		if (i == 0) printf("  %d   %d:  ",ifiber,iqie);
		else printf("\n  %d   %d:  ",ifiber,iqie);
		if (ierr == 0) printf("%1d:%2.2d/%1d       ",icap,mant,range);
		else printf("%1d:%2.2d/%1d(0x%X) ",icap,mant,range,ierr);
		old_qie = iqie;
		isamp = 1;
	  }
	}
	cout << endl;
	cout << "   [Note: if error=0 then not shown, otherwise is packed like this: {DV,ER,Errf1,Errf0}]"<< endl;
}

void check_qie(int type, int irun, int iev) {
	//
	// cycle through the QIEs, print info according to type.
	// type=0 just do statistics
	// type=1 statistics AND apply mantissa and range cuts
	//
	int nfeds = 0;
	int ifeds[50];  // ususally 32 max
	getFeds(&nfeds, ifeds, false, irun, iev);
	if (nfeds < 1) cout << "NO FEDS FOR THIS EVENT!!!!!" << endl;
	for (int i=0; i<nfeds; i++) {
	  const FEDRawData& data = rawdata->FEDData(ifeds[i]);
	  int f = ifeds[i] - 700;
	  FEDHeader header(data.data());
	  FEDTrailer trailer(data.data()+size-8);
	  payload=(uint32_t*)(data.data());
	  //
	  // what about the spigots?
	  //
	  setup_spigots(false);
	  for (int j=0; j<spigots; j++) {
	    htr_data_from_payload(payload,j);
	    htr_fill_header();
	    qies_from_htr();
	    for (int k=0; k<nqies; k++) {
		int m = qies[k][MANT];
		int r = qies[k][RANGE];
		int ifib = qies[k][FIBER];
		int iqie = qies[k][QIE];
		int icap = qies[k][CAPID];
		int fc = qie2fC[r][m];
		qiestats[f][j][ifib][iqie][icap][0]++;
		qiestatsAny[f][j] = 1;
		qiestats[f][j][ifib][iqie][icap][1] += fc;
		qiestats[f][j][ifib][iqie][icap][2] += fc*fc;
		bool rangeGT = r > qieEXPcut;
		bool rangeEQmantGT = (r == qieEXPcut) && (m >= qieMANTcut);
		if ( (type == 1) && (rangeGT || rangeEQmantGT) )
		  printf("===>Run=%6d Ev=%7d FED=%3d Spigot=%2d Fiber=%1d QIE=%1d CAPID=%1d Mant=%2d Exp=%1d\n",
		  irun,iev,ifeds[i],j,ifib,iqie,icap,m,r);
	    }
	    htr_data_delete();
	  }
	}
}

void qies_stats_print() {
	//
	// loop over qiestats, check if there were any increments [word 0]
	//
	cout << "Linearized QIE average and SD: " << endl;
	for (int i=0; i<100; i++) {
	  int ifed = i + 700;	// FED number
	  for (int j=0; j<16; j++) {
	    // spigot
	    if (qiestatsAny[i][j] > 0) {
	      cout << "FED number " << ifed << " Spigot " << j << ":" << endl;
	      printf(" Fiber QIE Capid        Av         SD     Num\n");
	      for (int k=0; k<8; k++) {
	        // fiber
	        for (int l=0; l<3; l++) {
		  // qie
		  for (int m=0; m<4; m++) {
		    // capid
		    int num = qiestats[i][j][k][l][m][0];
		    if (num > 0) {
		      double xnum = (double) num;
		      double av = qiestats[i][j][k][l][m][1]/xnum;
		      double av2 = qiestats[i][j][k][l][m][2]/xnum;
		      double sd = sqrt(av2 - av*av);
		      printf("   %1d    %1d    %1d   %9.2f  %9.2f  %6d\n",
			k,l,m,av,sd,num);
		    }
		  }
	        }
	      }
	    }
	  }
	}
}

void print_htr_tpgs() {
	cout << "TPs are next: " << endl;
	if (htype == 0) {
		//
		// HO is different than HBHE and HF as far as TPG go
		//
		//
		// first get the NumSamples
		//
		int ntpsamp = 0;
		for (int i=0; i<ntps; i++) {
			int tword = tpgs[i];
//			cout << "0x" << hex << tword << dec << endl;
			int high5 = (tword>>11)&0x1F;
			if (high5 == 0) ntpsamp++; // should be the same number for all 3 TPs (each one is a different set of 8 muon bits)
		}
		cout << "  This is HO - TPs are very different, and not all that well used (caveat emptor holds)" << endl;
		cout << "  I see " << ntps << " TP words here and " << ntpsamp << " samples for each muon bit:" << endl;
		printf("   =>Muon bits:222221111111111         \n");
		printf("    Sample SOI 432109876543210987654321\n");
		for (int i=0; i<ntpsamp; i++) {
			int t1 = tpgs[i];		// 1-8
			int t2 = tpgs[i+4];		// 9-16
			int t3 = tpgs[i+8];		// 17-24
			int soi1 = (t1>>9)&0x1;  // should be the same for all 3 words
			int soi2 = (t2>>9)&0x1;  // should be the same for all 3 words
			int soi3 = (t3>>9)&0x1;  // should be the same for all 3 words
			int t18 = t1&0xFF;
			int t28 = t2&0xFF;
			int t38 = t3&0xFF;
			if ( (soi1 != soi2) || (soi1 != soi3) || (soi2 != soi3) ) cout << "Problem!!! SOIs do not agree on TP " << i << endl;
			cout << "       " << i << "    " << soi1 << "  " << 
				(bitset<8>) t38 << (bitset<8>) t28 << (bitset<8>) t18 << endl;
		}
		cout << "   SOI is the 'sample of interest', 1 corresponds to BX of the L1A." << endl;
		cout << "   There are 24 HO channles per HTR and so 24 possible muon bits could be set, see above." << endl;
	}
	else {
		//
		// HBHE or HF, standard TPs.   assumes TPGS_FROM_HTR already called
		//
		int ntpsamp = ntps/ntptype[htype]; // should never be dividing by 0!
		cout << "  This is " << htr_types[htype] << ": ";
		cout << nslbtype[htype] << " SLB/HTR, " << ntptype[htype] << " TPs per HTR, " << 
			ntps << " trigger primitives in the data here " << endl;
		cout << "  There should be " << ntpsamp << " samples per channel and SOI should be at the same sample for all TPs" << endl;
		cout << "  SLB ID = 1-6 specifies which SLB" << endl;
		cout << "  SLB CH = 0-3 means A0,A1,C0,C1 for TOP FPGA and B0,B1,D0,D1 for BOT FPGA" << endl;
		cout << "  TP is the 9-bit number that gets sent to RCT" << endl;
//		cout << "    Sample SLB_ID  SLB_CH Z  SOI   TP" << endl;
		cout << "  SLB  CH:  ===> Samples (format n=TP(Z), SOI in bold)" << endl;
		int isamp = 1;
		int old_id = -1;
		int old_ch = -1;
		for (int i=0; i<ntps; i++) {
			int tword = tpgs[i];
			int tp_tag = (tword>>11)&0x1F;
			int slb_id = (tp_tag >> 2)&0x7;
			int slb_ch = tp_tag & 0x3;
			int z = (tword>>10)&0x1;
			int soi = (tword>>9)&0x1;
			int tp = tword&0x1FF;
			if ( (slb_id == old_id) && (slb_ch == old_ch) ) {
			  //
			  // just print the Z and TP, bold if SOI
			  //
			  if (soi == 1) printf("%s%d=0x%2.2X(%1d) %s",bold,isamp,tp,z,none);
			  else printf("%d=0x%2.2X(%1d) ",isamp,tp,z);
			}
			else {
			  //
			  // print SLB and channel number
			  //
			  if (slb_id == old_id) {
				if (i == 0) printf("        %1d:  ",slb_ch);
				else printf("\n        %1d:  ",slb_ch);
			  }
			  else {
				if (i == 0) printf("    %1d   %1d:  ",slb_id,slb_ch);
				else printf("\n    %1d   %1d:  ",slb_id,slb_ch);
			  }
			  if (soi == 1) printf("%s%d=0x%2.2X(%1d) %s",bold,isamp,tp,z,none);
			  else printf("%d=0x%2.2X(%1d) ",isamp,tp,z);
			  old_id = slb_id;
			  old_ch = slb_ch;
			}
//			printf("        %d     %d       %d   %d   %d  0x%2.2X\n",isamp,slb_id,slb_ch,z,soi,tp);
			isamp++;
			if (isamp > ntpsamp) isamp = 1;
		}
	}
}

void find_htr_header_errors() {
}

void print_error_warnings() {
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

void getFeds(int* nfeds, int* ifeds, bool printout, int irun, int iev) {
	int kfeds = 0;
	for (int i=700; i<800; i++) {
//	for (int i=700; i<FEDNumbering::lastFEDId(); i++) {
	  const FEDRawData& data = rawdata->FEDData(i);
//	  FEDHeader header(data.data());
//	  FEDTrailer trailer(data.data());
	  size=data.size();
//	    cout << "i=" << i << " size=" << size << " ok=" <<  rawdata.isValid() << endl;
	  if (size > 0) {
//		cout << " header check " << header.check() << " trailer check " << trailer.check() <<
//		" trailer status " << trailer.evtStatus() << endl;
//		cout << "Found FED " << i << endl;
		FEDHeader header(data.data());
		FEDTrailer trailer(data.data()+size-8);
		payload=(uint32_t*)(data.data());
		ifeds[kfeds++] = i;
	  }
	}
	*nfeds = kfeds;
	if (kfeds > 0 && printout) {
		printf("Run %6d Event %6d Feds: ",irun,iev);
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
	// (assume getFed already called)
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
void check_event_numbers(int irun, int iev) {
	int nfeds = 0;
	int ifeds[50];  // ususally 32 max
	getFeds(&nfeds, ifeds, false, irun, iev);
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
//	printf("==> Run %7d Ev %6d       Orn %10d : delta %10d\n",irun,iev,fedOrN,OrN_delta);	
	bool bcnidle = checkFedBcNIdle(nfeds, ifeds, &fedBcnIdle);
	bool result = bcnok && evnok && ornok && bcnidle;
	if (!result) cout << "  ===> Problem with BCN/EVN/ORN!!!" << endl;
//	else cout << " Run " << irun << " " << iev << " looks good..." << endl;
//	else cout << "." << endl;
//	return result && (iev < 4279);
//	return result;
	fflush(stdout);
	return;
}

void check_htr_header_errors(int irun,int iev) {
	int nfeds = 0;
	int ifeds[50];  // ususally 32 max
	getFeds(&nfeds, ifeds, false, irun, iev);
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
	return;
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
RawAnalyzer::RawAnalyzer(const edm::ParameterSet& iConfig) {
   //now do what ever initialization is needed

}


RawAnalyzer::~RawAnalyzer() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void RawAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& iSetup) {

   int this_run = e.id().run();
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
       //
       // stop on events where there are any errors in the HTR header
       //
       find_htr_header_errors();
       fflush(stdout);
       ncountse++;
       if (ncountse < nloopse) return;
       fflush(stdout);
       printbegin = true;
       searching = false;
     }
     else if (criteria == 2) {
       //
       // cycle through and printout stuff about ORN, BCN, EVN, BCN idle, etc
       //
       check_event_numbers(this_run,this_evn);
       fflush(stdout);
       ncountse++;
       if (ncountse < nloopse) return;
       fflush(stdout);
       printbegin = true;
       searching = false;
     }
     else if (criteria == 3) {
	//
	// report any HTR that has any errors set in it's header.  This will produce a lot of output!
	//
	check_htr_header_errors(this_run,this_evn);
	fflush(stdout);
	ncountse++;
	if (ncountse < nloopse) return;
	void print_error_warnings();
	fflush(stdout);
        printbegin = true;
	searching = false;
     }
     else if (criteria == 4) {
       //
       // cycle through and printout stuff if any QIE passes cuts
       //
       check_qie(1, this_run,this_evn);
       fflush(stdout);
       ncountse++;
       if (ncountse < nloopse) return;
       fflush(stdout);
//       qies_stats_print();
       printbegin = true;
       searching = false;
     }
     else if (criteria == 5) {
       //
       // cycle through and printout qie stats at the end
       //
       check_qie(0, this_run,this_evn);
       fflush(stdout);
       ncountse++;
       if (ncountse < nloopse) return;
       fflush(stdout);
       qies_stats_print();
       printbegin = true;
       searching = false;
     }
     else if (criteria == 6) {
       //
       // cycle through and printout list of FEDS so we can see if events are consistent
       //
       int nfeds = 0;
       int ifeds[50];  // ususally 32 max
       getFeds(&nfeds, ifeds, true, this_run, this_evn);
       fflush(stdout);
       ncountse++;
       if (ncountse < nloopse) return;
       fflush(stdout);
       printbegin = true;
       searching = false;
     }
   }
//
//   Handle<FEDRawDataCollection> rawdata;
   nFEDs = 0;
//   cout << "+++++++++ Looking for FEDs ++++++++++++++" << endl;
//   for (int i = 700; i<FEDNumbering::lastFEDId(); i++) {
   for (int i = 700; i<900; i++) {
	const FEDRawData& data = rawdata->FEDData(i);
	size=data.size();

	if (size>0 && (FEDids_.empty() || FEDids_.find(i)!=FEDids_.end())) {
//	  cout << " Found FED# " << setw(4) << i << " " << setw(8) << size << " bytes " << endl;
	  iFEDsize[nFEDs] = size;
	  iFEDs[nFEDs++] = i;
	}
   }
   cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
   if (nFEDs>0) {
      cout << "FED(#bytes): ";
      for (int i=0; i<nFEDs; i++) {
	printf("%d(%d) ",iFEDs[i],iFEDsize[i]);
        if ( (i+1)%10 == 0) cout << endl << "             ";
      }
      cout << endl;
   }
   else cout << "No FEDs found here.   Maybe this is not an HCAL file?" << endl;
   cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
   //
   // this might be old fashioned but what the heck...
   //	    
   char* creturn;
   int done = 0;
   while (done == 0) {	      
	  cout << "=========== Run " << this_run << " Event " << this_evn << " ===========";
	  cout << "===========================================" << endl;
	  cout << "  NEXT       Next                                  " << endl;
	  cout << "  FED        Select FED and report header info     " << endl;
	  cout << "  EVENT      Find event number and stop            " << endl;
	  cout << "  DCCHEX     DCC hex dump (unformatted)            " << endl;
	  cout << "  SHEADERS   All HTR payload headers               " << endl;
	  cout << "  SPEEK      HTR payload dump (formatted or not)   " << endl;
	  cout << "  QIE        Dump QIE data                         " << endl;
	  cout << "  LOOP       Loop and search for specifics (to uncover anomolies, and there are many) " << endl;
	  cout << "  QUIT       try to quit (throws timeout exception)" << endl;
	  creturn = vparse_input("Main>",9,"QUIT","NEXT","FED","EVENT","DCCHEX","SHEADERS","SPEEK","QIE","LOOP");
	  if (!strcmp(creturn,"UNKNOWN")) cout << "Hmm, not a legal command.  Try again? " << endl << endl;
	  else if (!strcmp(creturn,"QUIT")) throw cms::Exception("Timeout");
	  else if (!strcmp(creturn,"NEXT")) return;
	  else if (!strcmp(creturn,"FED")) {
	    //
	    // select which FED you want to see
	    //
	    cout << "There are " << nFEDs << " here: ";
	    for (int i=0; i<nFEDs; i++) cout << iFEDs[i] << " ";
	    cout << endl;
	    whichFED = getINT("Select FED (0=all)> ");
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
	  else if (!strcmp(creturn,"EVENT")) {
	    searching = true;
	    criteria = 0;
	    find_evn = getINT("Enter event number to look for> ");
	    return;
	  }
	  else if (!strcmp(creturn,"DCCHEX")) {
	    //
	    // unformatted hex dump of the entire DCC payload (all spigots)
	    //
	    if (isFEDopen) {
		    int ndcc32 = size/4;   //size is in bytes
		    if (debugit) cout << "# DCC words is " << ndcc32 << endl;
		    for (int i=0; i<ndcc32; i++) {
			if (i == 0) {
			    printf("DCC Header:\n");
			}
			else if (i == 6) {
			    printf("HTR payload summaries:\n");
			}
			for (int j=0; j<spigots; j++) if (i == ospigot[j]) { 
			    printf("HTR payload %d:\n",j+1);
			}
			printf("  %4d 0x%8.8X\n",i,payload[i]);
		    }
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (!strcmp(creturn,"SHEADERS")) {
	    //
	    // formatted dump of the headers for all HTRs (spigots)
	    //
	    int kfed = getINT("This FED (1) or all FEDs (0) :");
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
		print_error_warnings();
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (!strcmp(creturn,"SPEEK")) {
	    //
	    // dump payload for this spigot, but prompt for whether formatted or not, and which one first
	    //
	    if (isFEDopen) {
		cout << "There are " << spigots << " spigots here (starting at 0). ";
		int spign = getINT("Which spigot do you want? ");
		int doform = getINT("enter 1=formatted (easier to see), 0=unformatted (a bit more detailed): ");
		//
		// extract the data for this spigot
		//
		htr_data_from_payload(payload,spign);
		if (doform == 1) {
			//
			// printout in unformatted (or close to it anyway)
			//
			htr_data_print_formatted(spign);
		}
		else {
			//
			// printout in unformatted (or close to it anyway)
			//
			htr_data_print();
		}
		//
		// and delete the entire payload for this HTR
		//
		htr_data_delete();
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (!strcmp(creturn,"QIE")) {
	    //
	    // take a look at the QIE raw data, en mass
	    //
	    if (isFEDopen) {
		cout << "There are " << spigots << " spigots here (starting at 0). ";
		int spign = getINT("Which spigot do you want? ");
		htr_data_from_payload(payload,spign);
		htr_fill_header();
		qies_from_htr();
		print_htr_qies(spign);
		htr_data_delete();
	    }
	    else {
		cout << "PLEASE OPEN A FED BEFORE DOING ANYTHING!!!!!" << endl;
	    }
	  }
	  else if (!strcmp(creturn,"LOOP")) {
	    int done2 = 0;
	    criteria = 0;
	    cout << "===========================================";
	    cout << "===========================================" << endl;
	    cout << "  HERR         Search for event with error(s) in HTR header" << endl;
	    cout << "  EVNUM        Printout inconsistencies in BCN, ORN, Evn, BCNidle, etc" << endl;
	    cout << "  BERR         Report on any HTR that has any error bits set in the header" << endl;
	    cout << "  QIECUT       Set a cut and stop when you find any QIE that is greater or equal" << endl;
	    cout << "  QIESTATS     Run through events and calculate QIE Mean and SD per CAPID (and printout)" << endl;
	    cout << "  FEDLIST      Run through events and printout list of FEDs" << endl;
	    cout << "  QUIT         Return to MAIN" << endl;
	    while (done2 == 0) {
	      char* creturn2 = vparse_input("LOOP>",7,"QUIT","HERR","EVNUM","BERR","QIECUT","QIESTATS","FEDLIST");
	      if (!strcmp(creturn2,"UNKNOWN")) cout << "Hmm, not a legal command.  Try again? " << endl << endl;
	      else if (!strcmp(creturn2,"QUIT")) done2 = 1;
	      else if (!strcmp(creturn2,"HERR")) {
		criteria = 1;
		nloopse = getINT("How many events to loop over? :");
		ncountse = 0;
		printbegin = false;
		searching = true;
		return;
	      }
	      else if (!strcmp(creturn2,"EVNUM")) {
		criteria = 2;
		int iw = getINT("1=do not consider spigots 12 and 13 (calibration), 0=use both: ");		
		if (iw == 0) not12_13 = false;
		else not12_13 = true;
		nloopse = getINT("How many events to loop over? :");
		ncountse = 0;
		printbegin = false;
		searching = true;
		return;
	      }
	      else if (!strcmp(creturn2,"BERR")) {
		criteria = 3;
		nloopse = getINT("How many events to loop over? :");
		ncountse = 0;
		printbegin = false;
		searching = true;
		return;
	      }
	      else if (!strcmp(creturn2,"QIECUT")) {
		criteria = 4;
		nloopse = getINT("How many events to loop over? :");
		ncountse = 0;
		qieMANTcut = getINT("QIE cut (mantissa, 0-31): ");
		qieEXPcut = getINT("QIE cut (range, 0-3): ");
		cout << "OK, starting loop. There will be no reporting, please be patient...." << endl;
		printbegin = false;
		searching = true;
		//
		// clear some arrays
		//
	        for (int i=0; i<100; i++) {
		 for (int j=0; j<16; j++) {
		  qiestatsAny[i][j] = 0;
		  for (int k=0; k<8; k++) {
		   for (int l=0; l<3; l++) {
		    for (int m=0; m<4; m++) {
		     for (int n=0; n<3; n++) qiestats[i][j][k][l][m][n] = 0;
		    }//m
		   }//l
		  }//k
	 	 }//j
	        }//i
	        return;
	      }
	      else if (!strcmp(creturn2,"QIESTATS")) {
		criteria = 5;
		nloopse = getINT("How many events to loop over? :");
		cout << "OK, starting loop. There will be no reporting, please be patient...." << endl;
		ncountse = 0;
		printbegin = false;
		searching = true;
		//
		// clear some arrays
		//
	        for (int i=0; i<100; i++) {
		 for (int j=0; j<16; j++) {
		  qiestatsAny[i][j] = 0;
		  for (int k=0; k<8; k++) {
		   for (int l=0; l<3; l++) {
		    for (int m=0; m<4; m++) {
		     for (int n=0; n<3; n++) qiestats[i][j][k][l][m][n] = 0;
		    }//m
		   }//l
		  }//k
	 	 }//j
	        }//i
	        return;
	      }
	      else if (!strcmp(creturn2,"FEDLIST")) {
		criteria = 6;
		nloopse = getINT("How many events to loop over? :");
		ncountse = 0;
		printbegin = false;
		searching = true;
		return;
	      }
	    }
	  }
	    
   }
}


// ------------ method called once each job just before starting event loop  ------------
void RawAnalyzer::beginJob() {
  cout << "******************************************************************************************************************" << endl;
  cout << "* This is a sort of expert system for looking at raw HCAL data.   Much of it assumes an understanding of what's  *" << endl;
  cout << "* in that data.  To do that, you will need to pay attention to 2 documents, both of which might be evolving.     *" << endl;
  cout << "* Tullio's HTR documentation:  http://cmsdoc.cern.ch/cms/HCAL/document/CountingHouse/HTR/design/HTR_MainFPGA.pdf *" << endl;
  cout << "* Eric's DCC documentation:  http://cmsdoc.cern.ch/cms/HCAL/document/CountingHouse/DCC/FormatGuide.pdf           *" << endl; 
   printf("* Best of luck.   This version is %s%12s%s                                       (Drew Baden, drew@umd.edu)  *\n",
		bold,CodeVersion,none);
  cout << "* (Give someone a fish, they eat for a day.  Teach them to fish, they eat for a lifetime!)                       *" << endl;
  cout << "******************************************************************************************************************" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RawAnalyzer::endJob() {
}

char* vparse_input(const char* prompt, int num, ...)
{
	char* pch;
	char* pch2;
	char* cmd;
	char* cwhich;
	va_list arguments;
	int i;
	char cans[500];
	char bcans[500];
	//
	// get the input
	//
	do {
	   cmd=readline(prompt);
	   for (i=0; i<500 && cmd[i]!=0; i++) {
	      cans[i]=toupper(cmd[i]);
	      bcans[i] = cmd[i];
	   }
	   bcans[i] = 0;
	   cans[i++]=0;
	   if (i>0) add_history(bcans);
	   free(cmd);
	   printf("\n\n");
	   pch = strtok(cans," ");
	} while (pch == NULL);
	int alen = strlen(pch);
	//
	// put any other arguments onto the stack
	//
	nargs = 0;
	pch2 = strtok(bcans," ");
	pch2 = strtok(NULL," ");
	while (pch2 != NULL) {
	  args[nargs] = new char[100];
	  strcpy(args[nargs],pch2);
	  pch2 = strtok(NULL," ");
	  nargs++;
	}
	//
	// init arg list for variable number of args
	//
	va_start(arguments, num);
	//
	// loop over possible arguments, look for a match
	//
	for (i=0; i<num; i++) {
	  cwhich = va_arg(arguments, char*);
	  int icmp = strncmp(cwhich,pch,alen);
	  if (icmp == 0) return cwhich;
	}
	return (char*) "UNKNOWN";
}
int getINT(const char* prompt) {

	int temp;
	//
	// anything on the stack?
	//
	if (nargs == 0) {
	  printf("\n%s",prompt);
	  scanf("%d",&temp);
	  return temp;
	}
	else sscanf(args[0],"%d",&temp);
	moveup();
	return temp;
}
void moveup() {
	//
	// nargs-1 points to the top of the stack
	//
	if (nargs == 1) {
	  nargs = 0;
	  return;
	}
	//
	// more than 1 arg on the stack - move them all up
	//
	for (int i=0; i<nargs; i++) args[i] = args[i+1];
	nargs--;
}

//define this as a plug-in
DEFINE_FWK_MODULE(RawAnalyzer);
