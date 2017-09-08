/*
 Purpose: To read from "gainMeasurement_femb_1-parseBinaryFile.root" demonstrate the waveforms and identify the peaks using a fit.
 
 Storing Waveform - Process:
 1. Read Binary File with TFile::Open
 2. Read Tree from file
 3. Loop over entries in tree
 4. Store wf per subrun per channel in an object
 5. Plot the wf per subrun per channel (not done yet)
 6. Write each plot into a root file (not done yet)
 
 Identifying Peaks - Process:
 1. create histogram of waveform amplitudes
 2. define _
 3. _
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
        std::vector<int> peaks_x_1d;
        std::vector<int> peaks_y_1d;
    };
    waveform_struct wave;
    std::vector<waveform_struct> vector_struct[64][128]; //storage of structures
    void set_run_channel(int, int);
    void set_peaks_x(int,int,std::vector<unsigned short> wf[64][128], TH1F *wfH);
    void set_peaks_y(int,int,std::vector<unsigned short> wf[64][128]);
    void set_run_info();
    
    int get_npeaks(int,int);
    int get_peaks_x(int,int,int);
    int get_peaks_y(int,int,int);
    void get_run_info(int, int);
};

void Analyze::set_run_channel(int run, int channel){//1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}

void Analyze::set_peaks_x(int subrun, int channel, std::vector<unsigned short> wf[64][128], TH1F *wfH){//2. get peak locations
    //initialize high resolution peak search
    for (int j=0; j<wf_root[subrun][channel].size(); j++) {
        wfH -> SetBinContent(j,wf_root[subrun][channel].at(j));
    }
    TSpectrum *search = new TSpectrum(40);   // get peaks
    Int_t nfoundHD = search->Search(wfH,8, "", .1); //identify peak locations(@sigma = 8, threshold = 30%)
    Double_t *xpeaksHD = search->GetPositionX(); //get x-position
    int elements = wf[wave.subrun][wave.channel].size();
    double_t wf_mat[elements];
    double_t empty_mat[elements]; //dummy 1d array
    
    for(int s = 0; s < elements; s++) wf_mat[s]=wf[wave.subrun][wave.channel].at(s); // convert vector to double_t array
    wave.peaks_x_1d.clear();
    for (int i = 0; i < nfoundHD; i++){
        wave.peaks_x_1d.push_back(xpeaksHD[i]); // store x-positions
    }
    sort(wave.peaks_x_1d.begin(), wave.peaks_x_1d.end()); //sort peaks lowest to highest
    wave.nPeaks = nfoundHD;
}

void Analyze::set_peaks_y(int subrun, int channel,std::vector<unsigned short> wf[64][128]){
    //3
    if(DBG) cout << "\n\ti\txmin\txmax\tlow peak\thigh peak" << endl;
    if(DBG) cout << "\t====\t====\t===\t=========\t========" << endl;
    
    wave.peaks_y_1d.clear();
    for (int i = 0; i < wave.nPeaks-1; i++){
        int begin = wave.peaks_x_1d.at(i);
        int end = wave.peaks_x_1d.at(i+1);
        int max = *std::max_element(wf[wave.subrun][wave.channel].begin()+begin,wf[wave.subrun][wave.channel].begin()+end);
        int min = *std::min_element(wf[wave.subrun][wave.channel].begin()+begin,wf[wave.subrun][wave.channel].begin()+end);
        if(DBG) cout<<"\t"<<i<<"\t"<<begin<<"\t"<<end<<"\t"<< min << "\t\t" << max << endl;
        wave.peaks_y_1d.push_back(max);
    }
    //int m = *std::max_element(wf[wave.subrun][wave.channel].begin()+wave.peaks_x_1d.at(wave.nPeaks-1),wf[wave.subrun][wave.channel].end()); //store last max peak
    //wave.peaks_y_1d.push_back(m);
    for (int i = 0; i < wave.nPeaks-1; i++){ // needs to be fixed *Carlos
        int begin = wave.peaks_x_1d.at(i);
        int end = wave.peaks_x_1d.at(i+1);
        int min = *std::min_element(wf[wave.subrun][wave.channel].begin()+begin,wf[wave.subrun][wave.channel].begin()+end);
        wave.peaks_y_1d.push_back(min);
    }
    if(DBG) cout << "\n\tsubrun:\t" << subrun << "\tchannel:\t" << channel << "\tsaved peaks:\t" << wave.nPeaks << endl;
    //int n = *std::min_element(wf[wave.subrun][wave.channel].begin()+wave.peaks_x_1d.at(wave.nPeaks-1),wf[wave.subrun][wave.channel].end());
    //wave.peaks_y_1d.push_back(n);
    //if(DBG) cout<<"\t"<<wave.nPeaks-1<<"\t"<<wave.peaks_x_1d.at(wave.nPeaks-1)<<"\t"<< 0 <<"\t"<< n << "\t\t" << m << endl; //?
}
void Analyze::set_run_info() { // fix this
    //4
    vector_struct[wave.subrun][wave.channel].push_back(wave);
}
int Analyze::get_npeaks(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.nPeaks;
}
int Analyze::get_peaks_x(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    //cout << "size of x_:" << w.peaks_x_1d.size() << endl;
    return w.peaks_x_1d.at(element);
}
int Analyze::get_peaks_y(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    //cout << "size of y_:" << w.peaks_y_1d.size() << endl;
    return w.peaks_y_1d.at(element);
}
void Analyze::get_run_info(int subrun, int channel) {
    waveform_struct food = vector_struct[subrun][channel].at(0);
    cout << "channel: " << food.channel << endl;
    cout << "subrun: " << food.subrun << endl;
    cout << "peak locations (x): " << food.peaks_x_1d.size() << endl;
    cout << "peak heights (y): " << food.peaks_y_1d.size() << endl;
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
    int subrunLow=2, subrunHigh=subrunIn, subrunPoint = 10;
    int channelLow=127, channelHigh=chanIn, channelPoint = 127; //specific range of parameters *fix this -cp
    /*for (int s=subrunLow; s <(subrunHigh+1); s++) {
        for (int c=channelLow; c <(channelHigh+1); c++) {
            TH1F *wfH = new TH1F("wfH","ADC Signal",wf_root[s][c].size(),0,wf_root[s][c].size());
            foo.set_run_channel(s,c);
            foo.set_peaks_x(s, c, wf_root, wfH); // store x-positions
            foo.set_peaks_y(s, c, wf_root);
            foo.set_run_info(); // set all attributes and save to vector
        }
    }

    //show results
    TCanvas *c0 = new TCanvas("c0", "c0", 200, 10, 800, 600);
    TGraph *gCh = new TGraph();
    TGraph *gCh1 = new TGraph();
    char titleHist[40],titleGraph[40];
    
    c0->Divide(1,2); //canvas divide into two windows
    
    c0->cd(1); //draw waveform
    for(int s = 0; s < wf_root[subrunPoint][channelPoint].size(); s++){
        gCh->SetPoint(gCh->GetN() , s , wf_root[subrunPoint][channelPoint].at(s) );
    }
    for (int j = 0; j < (foo.get_npeaks(subrunPoint,channelPoint)-1); j++) {
        gCh1->SetPoint(gCh1->GetN(), foo.get_peaks_x(subrunPoint,channelPoint,j), foo.get_peaks_y(subrunPoint,channelPoint,j) );
    }
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    gCh->SetTitle(titleGraph);
    gCh->Draw("ALP");
    gCh1->Draw("*");

    c0->cd(2); //histogram of peaks of subrun and channel
    TH1F *h2 = new TH1F("h2","high & low peak counts",6000,0,6000);
    sprintf(titleHist,"Peaks (Subrun: %d, Channel %d)",subrunPoint,channelPoint);
    for (int j = 0; j<(foo.get_npeaks(subrunPoint,channelPoint)-1); j++) {
        h2->Fill(foo.get_peaks_y(subrunPoint,channelPoint,j));
    }
    h2->SetTitle(titleHist);
    h2->Draw();
    c0->Update();
    //*/
    
    
    // signal histogram plot & peaks
    TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 800, 600);
    TSpectrum *search = new TSpectrum(39);   // get peaks
    Double_t source[wf_root[subrunPoint][channelPoint].size()];
    TH1F *wfH = new TH1F("wfH","ADC Signal",wf_root[subrunPoint][channelPoint].size(),0,wf_root[subrunPoint][channelPoint].size());
    TH1F *nH = new TH1F("nH","ADC noise",wf_root[subrunPoint][channelPoint].size(),0,wf_root[subrunPoint][channelPoint].size());
    TH1F *spikes = new TH1F("spikes","spikes",wf_root[subrunPoint][channelPoint].size(),0,wf_root[subrunPoint][channelPoint].size());
    
    for (int j=0; j < wf_root[subrunPoint][channelPoint].size()-2; j++) wfH -> SetBinContent(j,wf_root[subrunPoint][channelPoint].at(j+1));
    
    for (int i = 0; i < wf_root[subrunPoint][channelPoint].size()-2; i++) source[i] = wfH->GetBinContent(i + 1); // for background subtraction
    
    // Estimate the background
    search->Background(source,wf_root[subrunPoint][channelPoint].size()-2,5,TSpectrum::kBackIncreasingWindow,
                  TSpectrum::kBackOrder2,kFALSE,
                  TSpectrum::kBackSmoothing3,kFALSE);
    
    // Draw the estimated background
    for (int i = 0; i < wf_root[subrunPoint][channelPoint].size()-2; i++) nH->SetBinContent(i + 1,source[i]);
    
    c1->Clear(0);
    Int_t nfound = search->Search(wfH, 4, "goff", .1);
    Double_t *xpeaks = search->GetPositionX(); //get x-position
    double a[nfound];
    for (int i = 0; i<nfound; i++) {
        a[i] = xpeaks[i];
    }
    Long64_t n1 = nfound;
    Int_t idx[nfound];
    TMath::Sort(nfound,a,idx,0);
    int interval = a[idx[1]]-a[idx[0]];
    for (int i = 0; i < (wf_root[subrunPoint][channelPoint].size()/interval); i++) {
        wfH->GetXaxis()->SetRange(interval*(i),interval*(i+1));
        Double_t currentMin = wfH->GetMinimum();
        Double_t binMin = wfH->GetMinimumBin();
        Double_t currentMax = wfH->GetMaximum();
        Double_t binMax = wfH->GetMaximumBin();
        
        spikes->SetBinContent(binMax,currentMax);
        spikes->SetBinContent(binMin,currentMin);
    }
    
    //1. sort
    //2. get indices i, i+1, i+2
    //3. get binContents @ i,i+1,i+2
    //4. cout << binContents verify that those bins are correctly peaks.
    //5. get average period (peak1+peak2/2)
    //6. get maximum and minimum peak height between ranges i & i+2
    
    
    wfH->GetXaxis()->SetRange(1,wf_root[subrunPoint][channelPoint].size());
    wfH->Draw();
    nH->SetLineColor(kGreen);
    nH->Draw("SAME L");
    spikes->SetLineColor(kRed);
    spikes->Draw("SAME L");

    c1->Update();
}


