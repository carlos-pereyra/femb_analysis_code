/*
 Purpose: To read from "gainMeasurement_femb_1-parseBinaryFile.root" demonstrate the waveforms and identify the peaks using a fit.
 
 1. Purpose: Acquire and analyze ADC Signal:
    a. Read Binary File with TFile::Open
    b. Read Tree from file
    c. Loop over entries in tree
    d. Store wf per subrun per channel in an object
    e. Plot the wf per subrun per channel (not done yet)
    f. Write each plot into a root file (not done yet)
 
 Identifying Peaks - Process:
 2. Purpose: create histogram of waveform amplitudes
    a. define histogram object
    b. pipe object to peak finding algorithm
    c. get height of positive and negative peaks
    d. subtract background noise to obtain true positive and negative peaks
 */

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

using namespace std;

unsigned short subrunIn1, chanIn1, subrunIn2, chanIn2;    //output tree and variable
std::vector<unsigned short> *wfIn1, *wfIn2;
TTree *tr_rawdata1, *tr_rawdata2;                  //ROOT TREE tr_rawdata variables
TFile *inputFile1, *inputFile2;

const int DBG = 0;
const int const_numSubrun = 64;
const int const_numChan = 128;

//------------ Control ------------
int subrunLow=2, subrunHigh=10, subrunPoint = 10; //controls peak finding limits (subrun and channel) - need baseline info (mean and rms) - found from subrun 1
int channelLow=127, channelHigh=127, channelPoint = 127;
int subrunBaseline=1;
const char file1[100] = "20170816T145947_fembTest_gainenc_test_g3_s3_extpulse";
const char file2[100] = "20170816T145947_fembTest_gainenc_test_g3_s3_extpulse";
const char path1[100] = "Root_Binary_Files/";
const char path2[100] = "Root_Binary_Files/";

Double_t baseline1, baseline2;
std::vector<unsigned short> wf_root1[64][128]; //store waveforms
std::vector<unsigned short> wf_root2[64][128]; //store waveforms

class Analyze {
    //store waveform information
    //class must be set in sequential order
    //int pos_npeaks,neg_npeaks; //n number of positive peaks
    
    public:
    
    struct waveform_struct{
        int subrun, channel, pos_npeaks, neg_npeaks, wf_elements;
        Double_t pos_peak_rms, pos_peak_mean, neg_peak_rms, neg_peak_mean;
        std::vector<Int_t> pos_peak_x;
        std::vector<Double_t> pos_peak_y;
        std::vector<Int_t> neg_peak_x;
        std::vector<Double_t> neg_peak_y;
    };
    waveform_struct wave;
    std::vector<waveform_struct> vector_struct[64][128]; //storage of structures
    
    //-----------------
    //----   set   ----
    //-----------------
    void set_run_channel(int, int);
    void set_peaks(int,int,std::vector<unsigned short> wf[64][128], Double_t, Int_t); // determine pos & neg peaks
    
    // set peak mean & rms values
    void set_pos_peak_mean(int, int);
    void set_neg_peak_mean(int,int);
    void set_pos_peak_rms(int,int);
    void set_neg_peak_rms(int,int);
    void set_run_info();
    
    //-----------------
    //----   get   ----
    //-----------------
    int get_pos_npeaks(int,int); //full number of peaks (pos & neg)?
    int get_neg_npeaks(int,int); //full number of peaks (pos & neg)?

    // peak locations & amplitudes
    Int_t get_pos_peak_x(int,int,int);
    Int_t get_neg_peak_x(int,int,int);
    Double_t get_pos_peak_y(int,int,int);
    Double_t get_neg_peak_y(int,int,int);
    
    // set peak mean & rms values
    Double_t get_pos_peak_mean(int,int);
    Double_t get_neg_peak_mean(int,int);
    Double_t get_pos_peak_rms(int,int);
    Double_t get_neg_peak_rms(int,int);
};
//----------------------------
//----   set definitions  ----
//----------------------------
void Analyze::set_run_channel(int run, int channel){//1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}
void Analyze::set_peaks(int subrun, int channel, std::vector<unsigned short> wf[64][128], Double_t baseline, Int_t interval){//2. get peak locations
    //initialize high resolution peak search
    wave.wf_elements = wf[subrun][channel].size();
    TH1F *wfH = new TH1F("","ADC Signal",wf_root1[subrunPoint][channelPoint].size(),0,wf_root1[subrunPoint][channelPoint].size());
    for (int j=0; j < wf[subrun][channel].size()-2; j++) wfH -> SetBinContent(j,wf[subrun][channel].at(j+1));
    
    TSpectrum *search = new TSpectrum(45);
    Int_t nfound = search->Search(wfH, 8, "goff", .25); //get peaks to get interval...
    Double_t *xpeaks = search->GetPositionX();  //get x-position
    double a[nfound];
    for (int i = 0; i<nfound; i++) {
        a[i] = xpeaks[i];
    }
    Long64_t n1 = nfound;
    Int_t idx[nfound];
    TMath::Sort(nfound,a,idx,0);
    //int interval = a[idx[1]]-a[idx[0]];
    cout << "subrun: " << subrun << " channel: " << channel << "\tinterval: " << interval << endl;
    int stop = (wf[subrun][channel].size()/interval);
    
    wave.pos_peak_x.clear();
    wave.pos_peak_y.clear();
    wave.neg_peak_x.clear();
    wave.neg_peak_y.clear();
    
    if(DBG) cout << "\n================================================================" << endl;
    if(DBG) cout << "\tSetting Low & High Peaks (Subrun: " << subrun << " ,Channel: " << channel << ")\n" << endl;
    if(DBG) cout << "\tX_High\tHigh_Peak\tX_Low\tLow_Peak"<< endl;
    if(DBG) cout << "\t======\t=========\t=====\t========"<< endl;
    for (int i = 0; i < stop; i++) {
        wfH->GetXaxis()->SetRange(interval*(i),interval*(i+1)); // problem with wfH
        
        Double_t binMin = wfH->GetMinimumBin();
        Double_t currentMin = wfH->GetMinimum();
        Double_t binMax = wfH->GetMaximumBin();
        Double_t currentMax = wfH->GetMaximum();
        
        wave.neg_peak_x.push_back(binMin);      // store low x-positions
        wave.neg_peak_y.push_back(currentMin);  // store low y-positions

        if (currentMax>(baseline+60)) {
            wave.pos_peak_x.push_back(binMax);      // store high x-positions
            wave.pos_peak_y.push_back(currentMax);  // store high y-positions
        }
        if(DBG) cout << "\t" << binMax << "\t   " << currentMax << "\t\t" << binMin << "\t"  << currentMin << endl;
    }
    if(DBG) cout << "\n" << endl;
    wfH->GetXaxis()->SetRange(1,wf[subrun][channel].size());
    wave.pos_npeaks = wave.pos_peak_y.size();
    wave.neg_npeaks = wave.neg_peak_y.size();
}
void Analyze::set_pos_peak_mean(int subrun, int channel){
    int sum = 0;
    TH1F *hpos = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).pos_npeaks; i++) hpos->Fill(vector_struct[subrun][channel].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak mean: " << hpos->GetMean() << endl;
    vector_struct[subrun][channel].at(0).pos_peak_mean = hpos->GetMean();
}
void Analyze::set_neg_peak_mean(int subrun, int channel){
    int sum = 0;
    TH1F *hneg = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).neg_npeaks; i++) hneg->Fill(vector_struct[subrun][channel].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak mean: " << hneg->GetMean() << endl;
    vector_struct[subrun][channel].at(0).neg_peak_mean = hneg->GetMean();
}
void Analyze::set_pos_peak_rms(int subrun, int channel){
    int sum = 0;
    TH1F *pos_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).pos_npeaks; i++) pos_rms->Fill(vector_struct[subrun][channel].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak rms: " << pos_rms->GetRMS() << endl;
    vector_struct[subrun][channel].at(0).pos_peak_rms = pos_rms->GetRMS();
}
void Analyze::set_neg_peak_rms(int subrun, int channel){
    int sum = 0;
    TH1F *neg_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).neg_npeaks; i++) neg_rms->Fill(vector_struct[subrun][channel].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak rms: " << neg_rms->GetRMS() << endl;
    vector_struct[subrun][channel].at(0).neg_peak_rms = neg_rms->GetRMS();
}
void Analyze::set_run_info() {
    vector_struct[wave.subrun][wave.channel].push_back(wave);
}
//----------------------------
//----   get definitions  ----
//----------------------------
int Analyze::get_pos_npeaks(int subrun,int channel){ // number of positive peaks only
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.pos_npeaks;
}
int Analyze::get_neg_npeaks(int subrun,int channel){ // number of positive peaks only
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.neg_npeaks;
}

Int_t Analyze::get_pos_peak_x(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.pos_peak_x.at(element);
}
Int_t Analyze::get_neg_peak_x(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.neg_peak_x.at(element);
}
Double_t Analyze::get_pos_peak_y(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    //if(DBG) cout << "\tGetting positive peak: " << w.pos_peak_y.at(element) << "\t at bin: " << w.pos_peak_x.at(element) << endl;
    if (element < w.pos_peak_x.size()) return w.pos_peak_y.at(element);
    else {
        cout << "requested element is outside positive peak matrix size (Analyze::get_pos_peak_y)" << endl;
        return 0;
    }
}
Double_t Analyze::get_neg_peak_y(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    //if(DBG) cout << "\tGetting negative peak: " << w.neg_peak_y.at(element) << "\t at bin: " << w.neg_peak_x.at(element) << endl;
    if (element < w.neg_peak_x.size()) return w.neg_peak_y.at(element);
    else {
        cout << "requested element is outside positive peak matrix size (Analyze::get_neg_peak_y)" << endl;
        return 0;
    }
}

Double_t Analyze::get_pos_peak_mean(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting positive peak mean: " << w.pos_peak_mean << endl;
    return w.pos_peak_mean;
}
Double_t Analyze::get_neg_peak_mean(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting negative peak mean: " << w.neg_peak_rms << endl;
    return w.neg_peak_mean;
}
Double_t Analyze::get_pos_peak_rms(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting positive peak rms: " << w.pos_peak_rms << endl;
    return w.pos_peak_rms;
}
Double_t Analyze::get_neg_peak_rms(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting negative peak rms: " << w.neg_peak_rms << endl;
    return w.neg_peak_rms;
}
//----------------------------------------------------------------------------
//----                                main                                ----
//----------------------------------------------------------------------------
void test(){
    //----------------------------
    //----   read root file   ----
    //----------------------------
    char location1[200],location2[200];
    sprintf(location1,"%s%s",path1,file1);
    sprintf(location2,"%s%s",path2,file2);
    
    TFile *inputFile1 = new TFile(location1, "READ"); //read input root file
    TFile *inputFile2 = new TFile(location2, "READ"); //read input root file
    //file1
    tr_rawdata1 = (TTree*) inputFile1->Get("femb_wfdata"); //initialize tr_rawdata branches
    if( !tr_rawdata1 ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    //file2
    tr_rawdata2 = (TTree*) inputFile2->Get("femb_wfdata"); //initialize tr_rawdata branches
    if( !tr_rawdata2 ){
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
    //file2
    tr_rawdata2->SetBranchAddress("subrun", &subrunIn2);  // initialize subrun branch
    tr_rawdata2->SetBranchAddress("chan", &chanIn2);      // initialize channel branch
    tr_rawdata2->SetBranchAddress("wf", &wfIn2);          // initialize waveform branch
    Long64_t nEntries2(tr_rawdata2->GetEntries());        //tr_rawdata branch-row length
    tr_rawdata2->GetEntry(0);                            // initialize tr_rawdata (pointer) to null
    for(Long64_t entry(0); entry<nEntries2; ++entry) {   // loop over input waveforms, group waveforms by subrun
        tr_rawdata2->GetEntry(entry);
        //make sure channels and subrun values are ok
        if( subrunIn2 < 0 || subrunIn2 >= const_numSubrun ) continue;
        if( chanIn2 < 0 || chanIn2 >= const_numChan ) continue;
        for( unsigned int s = 0 ; s < wfIn2->size() ; s++ ){        //store waveform vector in array for quick access
            wf_root2[subrunIn2][chanIn2].push_back( wfIn2->at(s) );
        }
    }
    
    cout << "max subrun1: " << subrunIn1 << " max channel2: " << chanIn1 << endl;
    cout << "max subrun2: " << subrunIn1 << " max channel2: " << chanIn1 << endl;
    inputFile1->Close();
    inputFile2->Close();
    
    //--------------------------
    //----   Show Results   ----
    //--------------------------
    TCanvas *c0 = new TCanvas("c0", "c0", 100, 10, 1600, 800);
    
    TH1F *h1_1 = new TH1F("h1_1","",2000,1,2000); //histogram subrun 1 wf #1
    TH1F *h1_2 = new TH1F("h1_2","",2000,1,2000); //histogram subrun 1 wf #2
    TH1F *hb_1 = new TH1F("hb_1","",2000,1,5000); //histogram subrun 1 baseline #1
    TH1F *hb_2 = new TH1F("hb_2","",2000,1,5000); //histogram subrun 1 baseline #2

    TH1F *h2_1 = new TH1F("h2_1","",2000,1,2000); //histogram subrun 2 wf #1
    TH1F *h2_2 = new TH1F("h2_2","",2000,1,2000); //histogram subrun 2 wf #2

    TH1F *h3_1 = new TH1F("h3_1","",2000,1,2000); //histogram wf
    TH1F *h3_2 = new TH1F("h3_2","",2000,1,2000); //histogram wf

    TH1F *h4_1 = new TH1F("h4_1","",2000,1,2000); //histogram wf
    TH1F *h4_2 = new TH1F("h4_2","",2000,1,2000); //histogram wf

    TH1F *h5_1 = new TH1F("h5_1","",2000,1,2000); //histogram wf
    TH1F *h5_2 = new TH1F("h5_2","",2000,1,2000); //histogram wf

    TH1F *h6_1 = new TH1F("h6_1","",2000,1,2000); //histogram wf
    TH1F *h6_2 = new TH1F("h6_2","",2000,1,2000); //histogram wf

    TH1F *h7_1 = new TH1F("h7_1","",2000,1,2000); //histogram wf
    TH1F *h7_2 = new TH1F("h7_2","",2000,1,2000); //histogram wf

    TH1F *h8_1 = new TH1F("h8_1","",2000,1,2000); //histogram wf
    TH1F *h8_2 = new TH1F("h8_2","",2000,1,2000); //histogram wf

    TH1F *h9_1 = new TH1F("h9_1","",2000,1,2000); //histogram wf
    TH1F *h9_2 = new TH1F("h9_2","",2000,1,2000); //histogram wf

    TH1F *h10_1 = new TH1F("h10_1","",2000,1,2000); //histogram wf
    TH1F *h10_2 = new TH1F("h10_2","",2000,1,2000); //histogram wf
    
    char titleHist[40],titleGraph[40],titleBaseline[40];
    c0->Divide(5,3);
    
    //---------------------------
    //----   1. ADC Signal   ----
    //---------------------------
    c0->cd(1);
    subrunPoint = 1;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h1_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h1_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h1_1->SetTitle(titleGraph);
    h1_1->SetMaximum(4100);
    h1_1->SetMinimum(0);
    h1_1->Draw();
    h1_2->SetLineColor(kRed);
    h1_2->Draw("same");
    
    c0->cd(2);
    subrunPoint = 2;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h2_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h2_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h2_1->SetTitle(titleGraph);
    h2_1->SetMaximum(4100);
    h2_1->SetMinimum(0);
    h2_1->Draw();
    h2_2->SetLineColor(kRed);
    h2_2->Draw("same");
    
    c0->cd(3);
    subrunPoint = 3;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h3_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h3_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h3_1->SetTitle(titleGraph);
    h3_1->SetMaximum(4100);
    h3_1->SetMinimum(0);
    h3_1->Draw();
    h3_2->SetLineColor(kRed);
    h3_2->Draw("same");
    
    c0->cd(4);
    subrunPoint = 4;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h4_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h4_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h4_1->SetTitle(titleGraph);
    h4_1->SetMaximum(4100);
    h4_1->SetMinimum(0);
    h4_1->Draw();
    h4_2->SetLineColor(kRed);
    h4_2->Draw("same");
    
    c0->cd(5);
    subrunPoint = 5;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h5_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h5_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h5_1->SetTitle(titleGraph);
    h5_1->SetMaximum(4100);
    h5_1->SetMinimum(0);
    h5_1->Draw();
    h5_2->SetLineColor(kRed);
    h5_2->Draw("same");
    
    c0->cd(6);
    subrunPoint = 6;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h6_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h6_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h6_1->SetTitle(titleGraph);
    h6_1->SetMaximum(4100);
    h6_1->SetMinimum(0);
    h6_1->Draw();
    h6_2->SetLineColor(kRed);
    h6_2->Draw("same");
    
    c0->cd(7);
    subrunPoint = 7;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h7_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h7_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h7_1->SetTitle(titleGraph);
    h7_1->SetMaximum(4100);
    h7_1->SetMinimum(0);
    h7_1->Draw();
    h7_2->SetLineColor(kRed);
    h7_2->Draw("same");
    
    c0->cd(8);
    subrunPoint = 8;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h8_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h8_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h8_1->SetTitle(titleGraph);
    h8_1->SetMaximum(4100);
    h8_1->SetMinimum(0);
    h8_1->Draw();
    h8_2->SetLineColor(kRed);
    h8_2->Draw("same");
    
    c0->cd(9);
    subrunPoint = 9;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h9_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h9_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h9_1->SetTitle(titleGraph);
    h9_1->SetMaximum(4100);
    h9_1->SetMinimum(0);
    h9_1->Draw();
    h9_2->SetLineColor(kRed);
    h9_2->Draw("same");
    
    c0->cd(10);
    subrunPoint = 10;
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) h10_1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunPoint][channelPoint].size(); s++) h10_2->SetBinContent(s , wf_root2[subrunPoint][channelPoint].at(s)); //draw input waveform
    h10_1->SetTitle(titleGraph);
    h10_1->SetMaximum(4100);
    h10_1->SetMinimum(0);
    h10_1->Draw();
    h10_2->SetLineColor(kRed);
    h10_2->Draw("same");

    //-------------------------------------------------------------------------------------
    //----   1. ADC Counts for wf#1 & wf#2 Baseline Comparison
    //-------------------------------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "c1", 100, 10, 1400, 700);
    for(int s = 0; s < wf_root1[subrunBaseline][channelPoint].size(); s++) hb_1->Fill(wf_root1[subrunBaseline][channelPoint].at(s)); //draw input waveform
    for(int s = 0; s < wf_root2[subrunBaseline][channelPoint].size(); s++) hb_2->Fill(wf_root2[subrunBaseline][channelPoint].at(s)); //draw input waveform
    baseline1 = hb_1->GetMean(); //Get Baseline Mean & RMS
    baseline2 = hb_2->GetMean();
    sprintf(titleGraph,"Baseline (Subrun: %d, Channel %d)",subrunBaseline,channelPoint);
    
    hb_1->GetXaxis()->SetTitle("ADC Counts");
    hb_1->GetYaxis()->SetTitle("Bin Count");
    hb_1->SetTitle(titleGraph);
    hb_1->Draw();
    hb_2->SetLineColor(kRed);
    hb_2->Draw("same");
    c1->SetLogy();
    c1->Update();
    
    //-----------------------------------
    //----   find peak information   ---- wf#1
    //-----------------------------------
    Analyze foo;
    for (int s=subrunLow; s <(subrunHigh+1); s++) { //+1 to reach max iterator
        for (int c=channelLow; c <(channelHigh+1); c++) { //+1 to reach max iterator
            foo.set_run_channel(s,c);
            foo.set_peaks(s, c, wf_root1, baseline1, 497); // store x-positions - 497 interval
            foo.set_run_info(); // set all attributes and save to vector - need this to access all information
        }
    }
    //--------------------------------------------------------------------------------
    //----   find peak statistics for variable subrun settings - single channel   ---- wf#1
    //--------------------------------------------------------------------------------
    for (int s=subrunLow; s <(subrunHigh+1); s++) { // = maximum subrun
        for (int c=channelLow; c <(channelHigh+1); c++) { // = maximum channel
            cout << "iterations [s++, c++] " << s << " , " << c << endl;
            foo.set_pos_peak_mean(s,c);
            foo.set_neg_peak_mean(s,c);
            foo.set_pos_peak_rms(s,c);
            foo.set_neg_peak_rms(s,c);
        }
    }
    
    c0->cd(11);
    TH1F *pp = new TH1F("pp","",10,0,10);
    TH1F *np = new TH1F("np","",10,0,10);
    sprintf(titleHist,"%s (Channel %d)",file1,channelPoint);
    for (int j = subrunLow; j < (subrunHigh+1); j++) pp->SetBinContent(j+1,foo.get_pos_peak_mean(j,channelPoint));
    for (int j = subrunLow; j < (subrunHigh+1); j++) np->SetBinContent(j+1,foo.get_neg_peak_mean(j,channelPoint));
    pp->SetFillStyle(3003);
    pp->SetFillColor(kBlack);
    pp->SetTitle(titleHist);
    pp->GetXaxis()->SetTitle("Subrun");
    pp->GetYaxis()->SetTitle("Mean Peaks");
    pp->SetMaximum(4100);
    pp->SetMinimum(0);
    pp->Draw();
    np->SetFillStyle(3003);
    np->SetFillColor(kBlue);
    np->Draw("same");
    
    // fit pos peak
    TF1 *g1 = new TF1("m1","pol1",subrunLow,10);
    Double_t par1[9];
    pp->Fit("m1","pol1","N",subrunLow,10);
    //g1->GetParameters(&par1[0]);
    //g1->SetLineColorAlpha(kBlue,0.35);
    //g1->SetLineStyle(2);
    g1->Draw("same");
    
    // fit neg peak
    TF1 *g2 = new TF1("m2","pol1",subrunLow,10);
    Double_t par2[9];
    np->Fit("m2","pol1","N",subrunLow,10);
    //g2->GetParameters(&par2[0]);
    //g2->SetLineColorAlpha(kMagenta,0.35);
    //g2->SetLineStyle(9);
    g2->Draw("same");
    
    c0->Update();
    sprintf(titleHist,"ADC_Performance_Plots/20170816T145947/%s.pdf",file1);
    c0->Print(titleHist);
    
    
    //-------------------------------------------------------------------------------------
    //----   2. ADC Signal Overlay with Located Peaks (Visual Peak Inspection)  ---- wf#1 -
    //-------------------------------------------------------------------------------------
    TCanvas *c2 = new TCanvas("c2", "c2", 100, 10, 1400, 700);
    TH1F *s1 = new TH1F("s1","",wf_root1[subrunPoint][channelPoint].size(),1,wf_root1[subrunPoint][channelPoint].size()); //histogram wf
    for (int j = 1; j < foo.get_pos_npeaks(subrunPoint, channelPoint); j++) s1->SetBinContent(foo.get_pos_peak_x(subrunPoint,channelPoint,j), foo.get_pos_peak_y(subrunPoint,channelPoint,j)); //draw positive peaks
    for (int j = 1; j < foo.get_neg_npeaks(subrunPoint, channelPoint); j++) s1->SetBinContent(foo.get_neg_peak_x(subrunPoint,channelPoint,j), foo.get_neg_peak_y(subrunPoint,channelPoint,j)); //draw negative peaks
    
    TH1F *wf1 = new TH1F("wf1","",wf_root1[subrunPoint][channelPoint].size(),1,wf_root1[subrunPoint][channelPoint].size()); //histogram wf
    for(int s = 0; s < wf_root1[subrunPoint][channelPoint].size(); s++) wf1->SetBinContent(s , wf_root1[subrunPoint][channelPoint].at(s)); //draw input waveform
    sprintf(titleGraph,"WF#1 & WF#2 Overlay (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    wf1->SetTitle(titleGraph);
    wf1->SetMaximum(4150);
    wf1->SetMinimum(0);
    wf1->Draw();
    s1->SetLineColor(kRed);
    s1->Draw("same");
    
    c2->Update();
    c2->Print("ADC_Performance_Plots/20170816T145947_fembTest_gainenc_test_g3_s2_extpulse/h2.pdf");

    //----------------------------------------------------------
    //----   3. Pos Peaks for Various Subruns (Histogram)   ---- wf #1
    //----------------------------------------------------------
    TCanvas *c3 = new TCanvas("c3", "c3", 100, 10, 1400, 700);
    c3->Divide(2,1);
    c3->cd(1);
    TH1F *ppp = new TH1F("pp","",10,0,10);
    TH1F *npp = new TH1F("np","",10,0,10);
    sprintf(titleHist,"%s (Channel %d)",file1,channelPoint);
    for (int j = subrunLow; j < (subrunHigh+1); j++) ppp->SetBinContent(j+1,foo.get_pos_peak_mean(j,channelPoint));
    for (int j = subrunLow; j < (subrunHigh+1); j++) npp->SetBinContent(j+1,foo.get_neg_peak_mean(j,channelPoint));
    ppp->SetFillStyle(3008);
    ppp->SetFillColor(kBlack);
    ppp->SetTitle(titleHist);
    ppp->GetXaxis()->SetTitle("Subrun");
    ppp->GetYaxis()->SetTitle("Mean Peaks");
    ppp->Draw();
    npp->SetFillStyle(3003);
    npp->SetFillColor(kBlue);
    npp->Draw("same");
    
    // fit pos peak
    TF1 *g11 = new TF1("m1","pol1",subrunLow,10);
    Double_t par11[9];
    ppp->Fit("m1","pol1","N",subrunLow,10);
    //g1->GetParameters(&par1[0]);
    g11->SetLineColorAlpha(kBlue,0.35);
    g11->SetLineStyle(2);
    g11->Draw("same");
    
    // fit neg peak
    TF1 *g22 = new TF1("m2","pol1",subrunLow,10);
    Double_t par22[9];
    npp->Fit("m2","pol1","N",subrunLow,10);
    //g2->GetParameters(&par2[0]);
    //g2->SetLineColorAlpha(kMagenta,0.35);
    //g2->SetLineStyle(9);
    g22->Draw("same");

    c3->cd(2);
    c3->Update();

}


