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

unsigned short subrunIn, chanIn; //output tree and variable
std::vector<unsigned short> *wfIn;

TTree *tr_rawdata; //ROOT TREE tr_rawdata variables
TFile* inputFile;

const int DBG = 0; //Debug Mode
const int const_numSubrun = 64;
const int const_numChan = 128;

unsigned short subrun=8; // specific parameters
unsigned short chan=127;

Double_t xp; //x-peaks
Double_t yp; //y-peaks

std::vector<unsigned short> wf_root[64][128]; //store waveforms

class Analyze {
    //store waveform information
    //class must be set in sequential order
    std::vector<int> peak_locations;
    int nPeaks;
    
    public:
    
    struct waveform_struct{
        int subrun, channel, nPeaks;
        std::vector<int> wave_1d; //amplitudes
        std::vector<int> peaks_x_1d;
        std::vector<int> peaks_y_1d;
    };
    waveform_struct wave;
    std::vector<waveform_struct> vector_struct[64][128]; //storage of structures
    void set_run_channel(int, int);
    void set_peaks_x(std::vector<unsigned short> wf[64][128]);
    void set_peaks_y(std::vector<unsigned short> wf[64][128]);
    void set_run_info();
    
    int get_peaks(int,int,int);
    void get_run_info(int, int);
};

void Analyze::set_run_channel(int run, int channel){
    //1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}

void Analyze::set_peaks_x(std::vector<unsigned short> wf[64][128]){
    //2. get peak locations
    
    //initialize high resolution peak search
    TSpectrum *searchHD = new TSpectrum(80);   // get peaks
    int elements = wf[wave.subrun][wave.channel].size();
    double_t wf_mat[elements];
    //----------------------------------------------------------------------
    for(int s = 0; s < elements; s++){
        wf_mat[s]=wf[wave.subrun][wave.channel].at(s); // convert vector to double_t array
    }
    double_t empty_mat[elements]; //dummy 1d array
    //identify peak locations(@sigma = 8, threshold = 30%)
    Int_t nfoundHD = searchHD->SearchHighRes(wf_mat,empty_mat, wf[wave.subrun][wave.channel].size(), 8, 35, kTRUE, 3, kTRUE, 3); //number of peaks
    cout << "nfoundHD: " << nfoundHD << endl;
    Double_t *xpeaksHD = searchHD->GetPositionX(); //get x-position
    //----------------------------------------------------------------------
    const Int_t nbins = wf[wave.subrun][wave.channel].size(); //number of elements in wf
    
    cout << "\n\ti\txmin\txmax" << endl;
    cout << "\t====\t====\t===" << endl;
    for (int i = 0; i < nfoundHD; i++){
        peak_locations.push_back(xpeaksHD[i]); // store x-positions
        cout<<"\t"<<i<<"\t"<< xpeaksHD[i] << endl;
    }
    sort(peak_locations.begin(), peak_locations.end()); //sort peaks lowest to highest
    wave.peaks_x_1d = peak_locations;
    wave.nPeaks = nfoundHD;
}

void Analyze::set_peaks_y(std::vector<unsigned short> wf[64][128]){
    //3
    cout << "\n\ti\txmin\txmax\tlow peak\thigh peak" << endl;
    cout << "\t====\t====\t===\t=========\t========" << endl;
    for (int i = 0; i < wave.nPeaks-1; i++){
        int begin = peak_locations[i];
        int end = peak_locations[i+1];
        int max = *std::max_element(wf[wave.subrun][wave.channel].begin()+begin,wf[wave.subrun][wave.channel].begin()+end);
        int min = *std::min_element(wf[wave.subrun][wave.channel].begin()+begin,wf[wave.subrun][wave.channel].begin()+end);
        cout<<"\t"<<i<<"\t"<<begin<<"\t"<<end<<"\t"<< min << "\t\t" << max << endl;
        wave.peaks_y_1d.push_back(max);
    }
    int m = *std::max_element(wf[wave.subrun][wave.channel].begin()+peak_locations[wave.nPeaks-1],wf[wave.subrun][wave.channel].end()); //store last max peak
    wave.peaks_y_1d.push_back(m);
    for (int i = 0; i < wave.nPeaks-1; i++){ // needs to be fixed *Carlos
        int begin = peak_locations[i];
        int end = peak_locations[i+1];
        int min = *std::min_element(wf[wave.subrun][wave.channel].begin()+begin,wf[wave.subrun][wave.channel].begin()+end);
        //vector<int>::iterator it;
        //auto it = find(wf[wave.subrun][wave.channel].begin(),wf[wave.subrun][wave.channel].end(),min); //- wf[wave.subrun][wave.channel].begin();
        //cout << "it: " << it << endl;
        //int pos = std::distance(wf[wave.subrun][wave.channel].begin(), it);
        //int location = std::distance(,min)
        //wave.peaks_x_1d.push_back(it);
        wave.peaks_y_1d.push_back(min);
    }
    int n = *std::min_element(wf[wave.subrun][wave.channel].begin()+peak_locations[wave.nPeaks-1],wf[wave.subrun][wave.channel].end());
    cout<<"\t"<<wave.nPeaks-1<<"\t"<<peak_locations[wave.nPeaks-1]<<"\t"<< 0 <<"\t"<< n << "\t\t" << m << endl; //?
    wave.peaks_y_1d.push_back(n);
    for (int j = 0; j < wf[wave.subrun][wave.channel].size(); j++) {
        wave.wave_1d.push_back(wf[wave.subrun][wave.channel].at(j));
    }
}

void Analyze::set_run_info() { // fix this
    //4
    vector_struct[wave.subrun][wave.channel].push_back(wave);
}

int Analyze::get_peaks(int subrun,int channel,int element){ // make this work
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.peaks_y_1d.at(element);
}

void Analyze::get_run_info(int subrun, int channel) {// make this work
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
        if(DBG)std::cout << wfIn->size() << " for entry : " << entry << std::endl;
        for( unsigned int s = 0 ; s < wfIn->size() ; s++ ){        //store waveform vector in array for quick access
            wf_root[subrunIn][chanIn].push_back( wfIn->at(s) );
        }
    }
    cout << "max subrun: " << subrunIn << " max channel " << chanIn << endl;
    inputFile->Close();
    
    //////
    
    //find peak heights
    Analyze foo;
    foo.set_run_channel(subrun,chan);
    foo.set_peaks_x(wf_root); // store x-positions
    foo.set_peaks_y(wf_root);
    foo.set_run_info(); // set all attributes and save to vector

    //show results
    TCanvas *c0 = new TCanvas("c0", "c0", 200, 10, 800, 600);
    TGraph *gCh = new TGraph();
    TGraph *gCh1 = new TGraph();
    char titleHist[40],titleGraph[40];
    
    c0->Divide(1,2); //canvas divide into two windows
    //s->Draw();
    for (int j; j<(foo.wave.peaks_y_1d.size()/2); j++) {
        gCh1->SetPoint(gCh1->GetN() , foo.wave.peaks_x_1d.at(j), foo.wave.peaks_y_1d.at(j) );
    }
    
    for(int s = 0; s < wf_root[subrun][chan].size(); s++){
        gCh->SetPoint(gCh->GetN() , s , wf_root[subrun][chan].at(s) );
    }
    c0->cd(1);
    gCh->Draw("ALP");
    gCh1->Draw("*");
    sprintf(titleGraph,"Waveform (Subrun: %d, Channel %d)",subrun,chan);
    gCh->SetTitle(titleGraph);
    
    c0->cd(2);
    TH1F *h2 = new TH1F("h2","high & low peak counts",6000,0,6000);
    sprintf(titleHist,"Peaks (Subrun: %d, Channel %d)",subrun,chan);
    h2->SetTitle(titleHist);
    for (int j; j<foo.wave.peaks_y_1d.size(); j++) {
        cout << j << "\tfood: " << foo.get_peaks(subrun,chan,j) << endl;
        h2->Fill(foo.get_peaks(subrun,chan,j));
    }
    h2->Draw();
    
    c0->Update();
}


