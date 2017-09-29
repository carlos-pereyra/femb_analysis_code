#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>
#include "TH1.h"
#include "TSpectrum.h"
#include "TROOT.h"
using namespace std;
TTree *tr_rawdata1;                  //ROOT TREE tr_rawdata variables
TFile *inputFile1;
unsigned short subrunIn1, chanIn1;    //output tree and variable
std::vector<unsigned short> *wfIn1, *wfIn2;
const int const_numSubrun = 64;
const int const_numChan = 128;
const char file1[100] = "20170816T145947_fembTest_gainenc_test_g3_s2_extpulse.root";
const char path1[100] = "Root_Binary_Files/";
std::vector<unsigned short> wf_root1[64][128]; //store waveforms


Double_t fitf(Double_t reltime, Double_t gain)
{
    Double_t fitval = 4.31054*exp(-2.94809*reltime)*gain
    -2.6202*exp(-2.82833*reltime)*cos(1.19361*reltime)*gain
    -2.6202*exp(-2.82833*reltime)*cos(1.19361*reltime)*cos(2.38722*reltime)*gain
    +0.464924*exp(-2.40318*reltime)*cos(2.5928*reltime)*gain
    +0.464924*exp(-2.40318*reltime)*cos(2.5928*reltime)*cos(5.18561*reltime)*gain
    +0.762456*exp(-2.82833*reltime)*sin(1.19361*reltime)*gain
    -0.762456*exp(-2.82833*reltime)*cos(2.38722*reltime)*sin(1.19361*reltime)*gain
    +0.762456*exp(-2.82833*reltime)*cos(1.19361*reltime)*sin(2.38722*reltime)*gain
    -2.620200*exp(-2.82833*reltime)*sin(1.19361*reltime)*sin(2.38722*reltime)*gain
    -0.327684*exp(-2.40318*reltime)*sin(2.5928*reltime)*gain +
    +0.327684*exp(-2.40318*reltime)*cos(5.18561*reltime)*sin(2.5928*reltime)*gain
    -0.327684*exp(-2.40318*reltime)*cos(2.5928*reltime)*sin(5.18561*reltime)*gain
    +0.464924*exp(-2.40318*reltime)*sin(2.5928*reltime)*sin(5.18561*reltime)*gain;
    
    return fitval;
}

void testfit(){
    //-------------------------------
    //----   1. read root file   ----
    //-------------------------------
    char location1[200];
    sprintf(location1,"%s%s",path1,file1);
    TFile *inputFile1 = new TFile(location1, "READ"); //read input root file
    //file1
    tr_rawdata1 = (TTree*) inputFile1->Get("femb_wfdata"); //initialize tr_rawdata branches
    if( !tr_rawdata1 ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
        }

        //file1
        tr_rawdata1->SetBranchAddress("subrun", &subrunIn1);  // initialize subrun branch
        tr_rawdata1->SetBranchAddress("chan", &chanIn1);      // initialize channel branch
        tr_rawdata1->SetBranchAddress("wf", &wfIn1);          // initialize waveform branch
        Long64_t nEntries1(tr_rawdata1->GetEntries());        //tr_rawdata branch-row length
        tr_rawdata1->GetEntry(0);                            // initialize tr_rawdata (pointer) to null
        for(Long64_t entry(0); entry<nEntries1; ++entry) {   // loop over input waveforms, group waveforms by subrun
            tr_rawdata1->GetEntry(entry);
            //make sure channels and subrun values are ok
            if( subrunIn1 < 0 || subrunIn1 >= const_numSubrun ) continue;
            if( chanIn1 < 0 || chanIn1 >= const_numChan ) continue;
            for( unsigned int s = 0 ; s < wfIn1->size() ; s++ ){        //store waveform vector in array for quick access
                wf_root1[subrunIn1][chanIn1].push_back( wfIn1->at(s) );
            }
        }
        cout << "max subrun1: " << subrunIn1 << " max channel2: " << chanIn1 << endl;
        inputFile1->Close();
    //---------------------------------
    
    TCanvas *c0 = new TCanvas("c0", "c0", 100, 10, 1000, 800);
    TH1F *hist = new TH1F("","",100,-20,20);
    hist->GetXaxis()->SetRange(-20,20);
    for (int i = -20; i<21; i++) {
        hist->SetBinContent(i,fitf(i,25));
    }
    hist->Draw();
    c0->Update();

};
