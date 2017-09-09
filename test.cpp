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
#include <vector>       // std::vector

#include "TH1.h"
#include "TSpectrum.h"

using namespace std;

unsigned short subrunIn, chanIn;    //output tree and variable
std::vector<unsigned short> *wfIn;

TTree *tr_rawdata;                  //ROOT TREE tr_rawdata variables
TFile* inputFile;

const int DBG = 1;                  //Debug Mode
const int const_numSubrun = 64;
const int const_numChan = 128;

std::vector<unsigned short> wf_root[64][128]; //store waveforms

class Analyze {
    //store waveform information
    //class must be set in sequential order
    int nPeaks;
    
    public:
    
    struct waveform_struct{
        int subrun, channel, nPeaks;
        std::vector<int> pos_peak_x;
        std::vector<int> pos_peak_y;
        std::vector<int> neg_peak_x;
        std::vector<int> neg_peak_y;
    };
    waveform_struct wave;
    std::vector<waveform_struct> vector_struct[64][128]; //storage of structures
    void set_run_channel(int, int);
    void set_peaks(int,int,std::vector<unsigned short> wf[64][128], TH1F *wfH);
    void set_run_info();
    
    int get_npeaks(int,int);
    int get_pos_peak_x(int,int,int);
    int get_pos_peak_y(int,int,int);
    int get_neg_peak_x(int,int,int);
    int get_neg_peak_y(int,int,int);
    void get_run_info(int, int);
};

void Analyze::set_run_channel(int run, int channel){//1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}

void Analyze::set_peaks(int subrun, int channel, std::vector<unsigned short> wf[64][128], TH1F *wfH){//2. get peak locations
    //initialize high resolution peak search
    for (int j=0; j < wf[subrun][channel].size()-2; j++) wfH -> SetBinContent(j,wf_root[subrun][channel].at(j+1));
    
    TSpectrum *search = new TSpectrum(40);
    Int_t nfound = search->Search(wfH, 4, "goff", .1); //get peaks
    Double_t *xpeaks = search->GetPositionX();  //get x-position
    double a[nfound];
    for (int i = 0; i<nfound; i++) {
        a[i] = xpeaks[i];
    }
    Long64_t n1 = nfound;
    Int_t idx[nfound];
    TMath::Sort(nfound,a,idx,0);
    int interval = a[idx[1]]-a[idx[0]];
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
    
        wave.pos_peak_x.push_back(binMax);      // store high x-positions
        wave.pos_peak_y.push_back(currentMax);  // store high y-positions
        wave.neg_peak_x.push_back(binMin);      // store low x-positions
        wave.neg_peak_y.push_back(currentMin);  // store low y-positions
        
        if(DBG) cout << "\t" << binMax << "\t   " << currentMax << "\t\t" << binMin << "\t"  << currentMin << endl;
    }
    if(DBG) cout << "\n" << endl;
    wfH->GetXaxis()->SetRange(1,wf[subrun][channel].size());
    wave.nPeaks = wave.pos_peak_y.size();
}
void Analyze::set_run_info() { // fix this
    //4
    vector_struct[wave.subrun][wave.channel].push_back(wave);
}
int Analyze::get_npeaks(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.nPeaks;
}
int Analyze::get_pos_peak_x(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.pos_peak_x.at(element);
}
int Analyze::get_pos_peak_y(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "\tGetting positive peak: " << w.pos_peak_y.at(element) << "\t at bin: " << w.pos_peak_x.at(element) << endl;
    if (element < w.pos_peak_x.size()) return w.pos_peak_y.at(element);
    else {
        cout << "requested element is outside positive peak matrix size (Analyze::get_pos_peak_y)" << endl;
        return 0;
    }
}
int Analyze::get_neg_peak_x(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.neg_peak_x.at(element);
}
int Analyze::get_neg_peak_y(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "\tGetting negative peak: " << w.neg_peak_y.at(element) << "\t at bin: " << w.neg_peak_x.at(element) << endl;
    if (element < w.neg_peak_x.size()) return w.neg_peak_y.at(element);
    else {
        cout << "requested element is outside positive peak matrix size (Analyze::get_neg_peak_y)" << endl;
        return 0;
    }
}
void Analyze::get_run_info(int subrun, int channel) {
    waveform_struct food = vector_struct[subrun][channel].at(0);
    cout << "channel: " << food.channel << endl;
    cout << "subrun: " << food.subrun << endl;
    cout << "positive peak locations (x): " << food.pos_peak_x.size() << endl;
    cout << "positive peak heights (y): " << food.pos_peak_y.size() << endl;
}
//


void test(){
    TFile *inputFile = new TFile("gainMeasurement_femb_1-parseBinaryFile.root", "READ"); //read input root file
    tr_rawdata = (TTree*) inputFile->Get("femb_wfdata"); //initialize tr_rawdata branches
    if( !tr_rawdata ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    tr_rawdata->SetBranchAddress("subrun", &subrunIn); // initialize subrun branch
    tr_rawdata->SetBranchAddress("chan", &chanIn); // initialize channel branch
    tr_rawdata->SetBranchAddress("wf", &wfIn); // initialize waveform branch
    Long64_t nEntries(tr_rawdata->GetEntries()); //tr_rawdata branch-row length
    tr_rawdata->GetEntry(0); // initialize tr_rawdata (pointer) to null
    for(Long64_t entry(0); entry<nEntries; ++entry) {        // loop over input waveforms, group waveforms by subrun
        tr_rawdata->GetEntry(entry);
        //make sure channels and subrun values are ok
        if( subrunIn < 0 || subrunIn >= const_numSubrun ) continue;
        if( chanIn < 0 || chanIn >= const_numChan ) continue;
        for( unsigned int s = 0 ; s < wfIn->size() ; s++ ){        //store waveform vector in array for quick access
            wf_root[subrunIn][chanIn].push_back( wfIn->at(s) );
        }
    }
    cout << "max subrun: " << subrunIn << " max channel " << chanIn << endl;
    inputFile->Close();
    
    //////
    
    //find peak heights
    Analyze foo;
    int subrunLow=2, subrunHigh=subrunIn, subrunPoint = 2;
    int channelLow=127, channelHigh=chanIn, channelPoint = 127; //specific range of parameters *fix this -cp
    TH1F *wfH = new TH1F("wfH","ADC Signal",wf_root[subrunPoint][channelPoint].size(),0,wf_root[subrunPoint][channelPoint].size());
    for (int s=subrunLow; s <(subrunHigh+1); s++) {
        for (int c=channelLow; c <(channelHigh+1); c++) {
            cout << "iterations [s++, c++] " << s << " , " << c << endl;
            foo.set_run_channel(s,c);
            foo.set_peaks(s, c, wf_root, wfH); // store x-positions
            foo.set_run_info(); // set all attributes and save to vector
        }
    }

    //show results
    TCanvas *c0 = new TCanvas("c0", "c0", 200, 10, 800, 600);
    TGraph *gCh = new TGraph();
    TGraph *gCh1 = new TGraph();
    char titleHist[40],titleGraph[40];
    
    c0->Divide(1,2); //canvas divide into two windows
    
    c0->cd(1);
    for(int s = 0; s < wf_root[subrunPoint][channelPoint].size(); s++){ //draw input waveform
        gCh->SetPoint(gCh->GetN() , s , wf_root[subrunPoint][channelPoint].at(s) );
    }
    for (int j = 0; j < foo.get_npeaks(subrunPoint,channelPoint); j++) { //draw positive peaks
        gCh1->SetPoint(gCh1->GetN(), foo.get_pos_peak_x(subrunPoint,channelPoint,j), foo.get_pos_peak_y(subrunPoint,channelPoint,j) );
        gCh1->SetPoint(gCh1->GetN(), foo.get_neg_peak_x(subrunPoint,channelPoint,j), foo.get_neg_peak_y(subrunPoint,channelPoint,j) );
    }
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    gCh->SetTitle(titleGraph);
    gCh->Draw("ALP");
    gCh1->Draw("*");

    c0->cd(2); //histogram of peaks of subrun and channel
    TH1F *h2 = new TH1F("h2","high & low peak counts",6000,0,6000);
    sprintf(titleHist,"Peaks (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for (int j = 0; j < foo.get_npeaks(subrunPoint,channelPoint); j++) {
        h2->Fill(foo.get_pos_peak_y(subrunPoint,channelPoint,j));
    }
    h2->SetTitle(titleHist);
    h2->Draw();
    c0->Update();
    
}


