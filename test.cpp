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
#include "TROOT.h"
using namespace std;

unsigned short subrunIn1, chanIn1;    //output tree and variable
std::vector<unsigned short> *wfIn1;
TTree *tr_rawdata1;                  //ROOT TREE tr_rawdata variables
TFile *inputFile1;

const int DBG = 0;
const int const_numSubrun = 64;
const int const_numChan = 128;

//------------ Control ------------
int subrunLow=2, subrunHigh=10, subrunPoint = 8; //controls peak finding limits (subrun and channel) - need baseline info (mean and rms) - found from subrun 1
int channelLow=125, channelHigh=127, channelPoint = 127;
int subrunBaseline=1;
const char file1[100] = "20170816T145947_fembTest_gainenc_test_g3_s3_extpulse.root";
const char path1[100] = "Root_Binary_Files/";

Double_t baseline1;
std::vector<unsigned short> wf_root1[64][128]; //store waveforms

void pretty_plot(TH1F *h, double line_color,
                 const char *histName, const char *axis_label_x, const char *axis_label_y,
                 int o, int xmin, int ymin, int xmax, int ymax){
    
    h->SetLineColor(line_color);
    h->SetTitle(histName);
    h->GetXaxis()->SetTitle(axis_label_x); //min x, max x
    h->GetYaxis()->SetTitle(axis_label_y); //min x, max x
    h->SetMaximum(ymax);//max y
    h->SetMinimum(ymin);//min y
    h->GetXaxis()->SetRange(xmin,xmax); //min x, max x
    h->SetStats(0);
    h->Draw();
    
    if (o==1) {
        TPad *pad = new TPad("food", histName, 0.6, 0.4, 0.85, 0.7);
        pad->Draw();
        pad->cd();
        h->SetStats(0);
        h->Draw();
    }
    
}

class Analyze {
    //store waveform information
    //class must be set in sequential order
    //int pos_npeaks,neg_npeaks; //n number of positive peaks
    
    public:
    
    struct waveform_struct{
        int subrun, channel, pos_npeaks, neg_npeaks, wf_elements;
        Double_t pos_peak_rms, pos_peak_mean, neg_peak_rms, neg_peak_mean, rise_time, width; //peak stats
        Double_t baseline;
        std::vector<Int_t> pos_peak_x;
        std::vector<Int_t> pos_peak_y;
        std::vector<Int_t> neg_peak_x;
        std::vector<Int_t> neg_peak_y;
    };
    waveform_struct wave;
    std::vector<waveform_struct> vector_struct[64][128]; //storage of structures
    
    //-----------------
    //----   set   ----
    //-----------------
    void set_run_channel(int, int);
    void set_run_info(int,int); //initializes vector of structures for specific subruns and channels

    void set_baseline(int,int,std::vector<unsigned short> wf[64][128]);
    
    // set peak mean & rms values
    void set_peaks(int,int,std::vector<unsigned short> wf[64][128], Int_t); // determine pos & neg peaks

    void set_pos_peak_mean(int, int);
    void set_neg_peak_mean(int,int);
    void set_pos_peak_rms(int,int);
    void set_neg_peak_rms(int,int);
    
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
    
    Double_t get_baseline(int,int);

    // set peak mean & rms values
    Double_t get_pos_peak_mean(int,int);
    Double_t get_neg_peak_mean(int,int);
    Double_t get_pos_peak_rms(int,int);
    Double_t get_neg_peak_rms(int,int);
    
    // peak quality
    Int_t get_risetime(int,int,std::vector<unsigned short> wf[64][128]);
    Int_t get_width(int,int,std::vector<unsigned short> wf[64][128]);
};

//----------------------------------------------------------------------------
//----                          set definitions                           ----
//----------------------------------------------------------------------------
void Analyze::set_run_channel(int run, int channel){//1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}

void Analyze::set_baseline(int subrun, int channel, std::vector<unsigned short> wf[64][128]){
    TString histName;
    histName.Form("Baseline Data (Subrun: %d, Channel: %d)", subrun, channelPoint); // rename!!! -cp
    
    TH1F *base_hist = new TH1F("","",5000,1,5000); //histogram subrun 1 baseline #1
    
    for(int s = 0; s < wf[subrun][channel].size(); s++) base_hist->Fill(wf[subrun][channel].at(s));
    
    Double_t baseline = base_hist->GetMean(); //Get Baseline Mean
    vector_struct[subrun][channel].at(0).baseline = baseline;
    
}

void Analyze::set_peaks(int subrun, int channel, std::vector<unsigned short> wf[64][128], Int_t interval){//2. get peak locations
    TH1F *wfH = new TH1F("","ADC Signal",wf[subrun][channelPoint].size(),0,wf[subrun][channelPoint].size());
    
    for(int s = 0; s < wf[subrunPoint][channelPoint].size(); s++) wfH->SetBinContent(s , wf[subrun][channel].at(s)); //draw input waveform

    int stop = (wf[subrun][channel].size()/interval);
    
    std::vector<Int_t> p_peak_x;
    std::vector<Int_t> p_peak_y;
    std::vector<Int_t> n_peak_x;
    std::vector<Int_t> n_peak_y;
    p_peak_x.clear(); //setup vectors instead inside the definition
    p_peak_y.clear();
    n_peak_x.clear();
    n_peak_y.clear();
    
    if(DBG) cout << "\n================================================================" << endl;
    if(DBG) cout << "\tSetting Low & High Peaks (Subrun: " << subrun << " ,Channel: " << channel << ")\n" << endl;
    if(DBG) cout << "\tX_High\tHigh_Peak\tX_Low\tLow_Peak"<< endl;
    if(DBG) cout << "\t======\t=========\t=====\t========"<< endl;
    
    Double_t baseline = vector_struct[subrun][channel].at(0).baseline;
    for (int i = 0; i < stop; i++) {
        wfH->GetXaxis()->SetRange(interval*(i),interval*(i+1)); // problem with wfH
        
        Double_t binMin = wfH->GetMinimumBin();
        Double_t currentMin = wfH->GetMinimum();
        Double_t binMax = wfH->GetMaximumBin();
        Double_t currentMax = wfH->GetMaximum();
    
        if (currentMax>(baseline+60)) {
            p_peak_x.push_back(binMax);
            vector_struct[subrun][channel].at(0).pos_peak_x.push_back(binMax);
            
            p_peak_y.push_back(currentMax);
            vector_struct[subrun][channel].at(0).pos_peak_y.push_back(currentMax);
        }
        
        if (currentMin<(baseline-60)) {
            n_peak_x.push_back(binMin);
            vector_struct[subrun][channel].at(0).neg_peak_x.push_back(binMin);
            
            n_peak_y.push_back(currentMin);
            vector_struct[subrun][channel].at(0).neg_peak_y.push_back(currentMin);
        }
    }
    if(DBG) for (int i = 0; i < 5; i++) cout << "\t" << vector_struct[subrun][channel].at(0).pos_peak_x.at(i) << "\t   " << vector_struct[subrun][channel].at(0).pos_peak_y.at(i) << "\t\t" << vector_struct[subrun][channel].at(0).neg_peak_x.at(i) << "\t"  << vector_struct[subrun][channel].at(0).neg_peak_y.at(i) << endl;
    if(DBG) cout << "\n" << endl;
    
    wfH->GetXaxis()->SetRange(1,wf[subrun][channel].size());
    
    vector_struct[subrun][channel].at(0).pos_npeaks = vector_struct[subrun][channel].at(0).pos_peak_y.size(); // not finished after equals sign
    vector_struct[subrun][channel].at(0).neg_npeaks = vector_struct[subrun][channel].at(0).neg_peak_y.size(); // not finished after equals sign
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

void Analyze::set_run_info(int subrun, int channel) {
    vector_struct[subrun][channel].push_back(wave);
}

//----------------------------------------------------------------------------
//----                          get definitions                           ----
//----------------------------------------------------------------------------
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

Double_t Analyze::get_baseline(int subrun, int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.baseline;
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

Int_t Analyze::get_risetime(int subrun, int channel, std::vector<unsigned short> wf[64][128]){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    float ratio;
    Double_t y1 = (w.pos_peak_y.at(0)-w.baseline), y2;
    Int_t risetime;
    
    for (int i = 0; i < 10; i++) {
        y2 = (wf[subrun][channel].at(w.pos_peak_x.at(0)-i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.01) break;
        risetime = i+1;
        
        if(DBG) cout << "ratio: " << ratio << endl;
        if(DBG) cout <<"x1: "<< w.pos_peak_x.at(0) <<" y1: "<<y1<< " x2: " << risetime << " y2: " << y2 << endl;
    }
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting risetime: " << risetime << endl;
    vector_struct[subrun][channel].at(0).rise_time = risetime;
    return risetime;
}

Int_t Analyze::get_width(int subrun, int channel, std::vector<unsigned short> wf[64][128]){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    float ratio;
    Double_t y1 = (w.pos_peak_y.at(0)-w.baseline), y2;
    Int_t width, L, R;
    
    for (int i = 0; i < 10; i++) {
        y2 = (wf[subrun][channel].at(w.pos_peak_x.at(0)-i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.02) break;
        L = i+1;
        
        if(DBG) cout << "L ratio: " << ratio << endl;
        if(DBG) cout <<"x1: "<< w.pos_peak_x.at(0) <<" y1: "<<y1<< " x2: " << L << " y2: " << y2 << endl;
    }
    
    for (int i = 0; i < 10; i++) {
        y2 = (wf[subrun][channel].at(w.pos_peak_x.at(0)+i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.02) break;
        R = i+1;
        
        if(DBG) cout << "R ratio: " << ratio << endl;
        if(DBG) cout <<"x1: "<< w.pos_peak_x.at(0) <<" y1: "<<y1<< " x2: " << R << " y2: " << y2 << endl;
    }
    
    width = L + R + 1;
    
    //if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting risetime: " << risetime << endl;
    vector_struct[subrun][channel].at(0).width = width;
    return width;
}

double fudge_sundae(double *xv, double *par){
    
    double x= 10./1000.*xv[0]; //try converting x to 10^-3
    
    //if (x <=  0) return 0;
    //if (x >= 10) return 0;
    
    x /= par[1];
    
    double phi_1 = x * 1.19361;
    double phi_2 = x * 2.59280;
    
    double f = 4.31054 * TMath::Exp(-0.119760*x) -5.240400*TMath::Cos(phi_1) + 1.524912*TMath::Sin(phi_1);
    
    f *= TMath::Exp(-0.425150*x);
    
    f += 0.929848*TMath::Cos(phi_2) -0.655368*TMath::Sin(phi_2);
    
    f *= TMath::Exp(-2.40318*x) * par[0];
    
    f -= 4.31054*exp(-2.94809*xv[0]/par[1])*par[0]-2.6202*exp(-2.82833*xv[0]/par[1])*cos(1.19361*xv[0]/par[1])*par[0]
    -2.6202*exp(-2.82833*xv[0]/par[1])*cos(1.19361*xv[0]/par[1])*cos(2.38722*xv[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*xv[0]/par[1])*cos(2.5928*xv[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*xv[0]/par[1])*cos(2.5928*xv[0]/par[1])*cos(5.18561*xv[0]/par[1])*par[0]
    +0.762456*exp(-2.82833*xv[0]/par[1])*sin(1.19361*xv[0]/par[1])*par[0]
    -0.762456*exp(-2.82833*xv[0]/par[1])*cos(2.38722*xv[0]/par[1])*sin(1.19361*xv[0]/par[1])*par[0]
    +0.762456*exp(-2.82833*xv[0]/par[1])*cos(1.19361*xv[0]/par[1])*sin(2.38722*xv[0]/par[1])*par[0]
    -2.6202*exp(-2.82833*xv[0]/par[1])*sin(1.19361*xv[0]/par[1])*sin(2.38722*xv[0]/par[1])*par[0]
    -0.327684*exp(-2.40318*xv[0]/par[1])*sin(2.5928*xv[0]/par[1])*par[0] +
    +0.327684*exp(-2.40318*xv[0]/par[1])*cos(5.18561*xv[0]/par[1])*sin(2.5928*xv[0]/par[1])*par[0]
    -0.327684*exp(-2.40318*xv[0]/par[1])*cos(2.5928*xv[0]/par[1])*sin(5.18561*xv[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*xv[0]/par[1])*sin(2.5928*xv[0]/par[1])*sin(5.18561*xv[0]/par[1])*par[0];
    
    return f;
}

//----------------------------------------------------------------------------
//----                                main                                ----
//----------------------------------------------------------------------------
void test(){
    //-------------------------------
    //----   1. read root file   ----
    //-------------------------------
    char location1[200];
    sprintf(location1,"%s%s",path1,file1);
    TFile *inputFile1 = new TFile(location1, "READ");   //read input root file
    tr_rawdata1 = (TTree*) inputFile1->Get("femb_wfdata");  //initialize tr_rawdata branches
    if( !tr_rawdata1 ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    tr_rawdata1->SetBranchAddress("subrun", &subrunIn1);    // initialize subrun branch
    tr_rawdata1->SetBranchAddress("chan", &chanIn1);        // initialize channel branch
    tr_rawdata1->SetBranchAddress("wf", &wfIn1);            // initialize waveform branch
    Long64_t nEntries1(tr_rawdata1->GetEntries());          //tr_rawdata branch-row length
    tr_rawdata1->GetEntry(0);                               // initialize tr_rawdata (pointer) to null
    for(Long64_t entry(0); entry<nEntries1; ++entry) {      // loop over input waveforms, group waveforms by subrun
        tr_rawdata1->GetEntry(entry);
        if( subrunIn1 < 0 || subrunIn1 >= const_numSubrun ) continue;
        if( chanIn1 < 0 || chanIn1 >= const_numChan ) continue;
        for( unsigned int s = 0 ; s < wfIn1->size() ; s++ ){        //store waveform vector in array for quick access
            wf_root1[subrunIn1][chanIn1].push_back( wfIn1->at(s) );
        }
    }
    cout << "max subrun1: " << subrunIn1 << " max channel2: " << chanIn1 << endl;
    inputFile1->Close();
    
    //-------------------------------------
    //----   3. set peak information   ----
    //-------------------------------------
    Analyze foo;
    for (int s=1; s <(subrunHigh+1); s++) { //+1 to reach max iterator
        for (int c=channelLow; c <(channelHigh+1); c++) { //+1 to reach max iterator
            foo.set_run_info(s,c); // set all attributes and save to vector - need this to access all information
        }
    }
    for (int s=subrunLow; s <(subrunHigh+1); s++) { //+1 to reach max iterator
        for (int c=channelLow; c <(channelHigh+1); c++) { //+1 to reach max iterator
            foo.set_run_channel(s,c);
            foo.set_baseline(s,c,wf_root1);
            foo.set_peaks(s,c,wf_root1, 400); // store x-positions - 497 interval
        }
    }
    //----------------------------------------------------------------------------------
    //----  4. find peak statistics for variable subrun settings - single channel   ----
    //----------------------------------------------------------------------------------
    for (int s=subrunLow; s <(subrunHigh+1); s++) { // = maximum subrun
        for (int c=channelLow; c <(channelHigh+1); c++) { // = maximum channel
            if(DBG) cout << "peak finding [subrun: " << s << " ,channel: " << c << "]" << endl;
            foo.set_pos_peak_mean(s,c);
            foo.set_neg_peak_mean(s,c);
            foo.set_pos_peak_rms(s,c);
            foo.set_neg_peak_rms(s,c);
        }
    }
    
    //------------------------------------
    //----   5. ADC Signal & Summary  ----
    //------------------------------------
    TCanvas *c0 = new TCanvas("c0", "c0", 100, 10, 1000, 800);
    c0->Divide(3,3);
    char titleHist[40],titleGraph[40],titleBaseline[40];
    for (int sub = 1; sub<9; sub++) {
        TString histName;
        histName.Form("(Subrun: %d, Channel: %d)", sub, channelPoint);
        TH1F *h1 = new TH1F(histName,"",4000,1,4000); //1000 - x axis
        for(int j = 0; j < wf_root1[sub][channelPoint].size(); j++) h1->SetBinContent(j , wf_root1[sub][channelPoint].at(j)); //draw input waveform
        c0->cd(sub);
        pretty_plot(h1, 1, histName, "time (no units)","ADC",0,0,0,4000,4096);
    }
    c0->Modified();
    c0->Update();
    
    c0->cd(9);
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
    g1->SetLineColorAlpha(kBlue,0.35);
    g1->SetLineStyle(2);
    g1->Draw("same");
    
    g1->GetParameters(&par1[0]);
    
    // fit neg peak
    TF1 *g2 = new TF1("m2","pol1",subrunLow,10);
    Double_t par2[9];
    np->Fit("m2","pol1","N",subrunLow,10);
    g2->SetLineColorAlpha(kMagenta,0.35);
    g2->SetLineStyle(9);
    g2->Draw("same");
    
    g2->GetParameters(&par2[0]);
    
    sprintf(titleHist,"ADC_Performance_Plots/20170816T145947/%s.png",file1);
    c0->Print(titleHist);
    
    // 1. focus on signal (0-2000). 2. focus on signal (0-500)(show risetime) 3. focus on baseline (2200-2400) (fit guassian)
    TCanvas *c1 = new TCanvas("c1", "c1", 100, 10, 1000, 800);
    
    c1->Divide(3,3);
    for (int i = 0; i < 3; i++) {
        
        TString histName;
        int pad1 = i*3+1;
        int pad2 = i*3+2;
        int pad3 = i*3+3;
        int subrun = (i+1)*3;
        
        //waveform & risetime ------------------------------------------------------------
        c1->cd(pad1);
        
        histName.Form("First Peak R. Fit (Subrun: %d, Channel: %d)", subrun, channelPoint);
        TH1F *h1 = new TH1F(histName,"",4000,1,4000); //histogram subrun 1 wf #1
        
        for(int j=0; j<500; j++) h1->SetBinContent(j , (wf_root1[subrun][channelPoint].at(foo.get_pos_peak_x(subrun,channelPoint,1)-foo.get_risetime(subrun, channelPoint, wf_root1)+j)-foo.get_baseline(subrun,channelPoint))/(1) ); // somethings wrong here!!!
        
        pretty_plot(h1, 1, histName, "time (no units)", "ADC", 0,0,0,25,2000);
        
        TF1 *f3 = new TF1("func3",fudge_sundae,0,500,2);
        
        f3->SetParameters(140.*1.012,1.0); // 14.0mV/fC for now with 2.0 us shaping time
        
        h1->Fit("func3","","",0,foo.get_width(9,channelPoint, wf_root1)); //fit range* muy importante

        auto legend = new TLegend(0.5,0.7,0.85,0.9);
        legend->AddEntry(f3,"Histogram filled with random numbers");
        legend->Draw();
        
        //baseline ------------------------------------------------------------
        c1->cd(pad2);
        
        histName.Form("Baseline (Subrun: %d, Channel: %d)", subrun, channelPoint);
        TH1F *h2 = new TH1F(histName,"",4096,1,4096); //histogram subrun 1 wf #1
        
        for(int j = 0; j < 10000; j++) h2->Fill(wf_root1[subrun][channelPoint].at(j)); //draw input waveform
        
        pretty_plot(h2, 1, histName, "ADC Bins", "Count",0,foo.get_baseline(subrun,channelPoint)-75,0,foo.get_baseline(subrun,channelPoint)+75,1000);
        
        //show fitted function to one peak ------------------------------------------------------------
        c1->cd(pad3);
        
        histName.Form("Response (Subrun: %d, Channel: %d)", subrun, channelPoint);
        TH1F *h3 = new TH1F(histName,"",1000,0,1000); //histogram subrun 1 wf #1
        
        //for (Int_t i=0;i!=999;i++) h3->SetBinContent(i, f3->Eval((i+1)) );

        pretty_plot(h3, 1, histName, "time (?? unknown units)", "ADC",0,0,0,25,2000);
        
    }
    
    c1->Modified();
    c1->Update();
    
    //-------------------------------------------------------------------------------------
    //----   6. ADC Signal Overlay with Located Peaks (Visual Peak Inspection)  ---- wf#1 -
    //-------------------------------------------------------------------------------------
    if (DBG) {
        TCanvas *c2 = new TCanvas("c2", "c2", 100, 10, 1400, 700);
        TH1F *s1 = new TH1F("s1","",wf_root1[subrunPoint][channelPoint].size(),1,wf_root1[subrunPoint][channelPoint].size()); //histogram wf
        for (int j = 1; j < foo.get_pos_npeaks(subrunPoint, channelPoint); j++) s1->SetBinContent(foo.get_pos_peak_x(subrunPoint,channelPoint,j), foo.get_pos_peak_y(subrunPoint,channelPoint,j)); //draw positive peaks
        for (int j = 1; j < foo.get_neg_npeaks(subrunPoint, channelPoint); j++) s1->SetBinContent(foo.get_neg_peak_x(subrunPoint,channelPoint,j), foo.get_neg_peak_y(subrunPoint,channelPoint,j)); //draw negative peaks
        
        TH1F *wf1 = new TH1F("wf1","",wf_root1[subrunPoint][channelPoint].size()+1,0,wf_root1[subrunPoint][channelPoint].size()); //histogram wf
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
    }
}


