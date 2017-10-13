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

#include "analyze.h"
#include "management.h"

using namespace std;

unsigned short subrunIn1, chanIn1;    //output tree and variable
std::vector<unsigned short> *wfIn1;
TTree *tr_rawdata1;                  //ROOT TREE tr_rawdata variables
TFile *inputFile1;

const int const_numSubrun = 64;
const int const_numChan = 128;

const char file1[100] = "20170816T145947_fembTest_gainenc_test_g3_s2_extpulse.root";
const char path1[100] = "Root_Binary_Files/";

std::vector<unsigned short> wf_root1[64][128]; //store waveforms

//----------------------------------------------------------------------------
//----                                main                                ----
//----------------------------------------------------------------------------
void test(){
    
    vector<double> par = getParameters();
    int subrunLow = par.at(0);
    int subrunHigh = par.at(1);
    
    int channelLow = par.at(2);
    int channelHigh = par.at(3);
    
    int subrunPoint = par.at(4);
    int channelPoint = par.at(5);
    
    vector<string> files = getFiles();
    string file = files.at(0);
    
    //
    // read root file
    //
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
    cout << "max subrun: " << subrunIn1 << " max channel: " << chanIn1 << endl;
    inputFile1->Close();
    
    //
    // initialize peak information
    //
    
    Analyze foo;
    
    for (int s=subrunLow; s <=subrunHigh; s++) { //+1 to reach max iterator
        for (int c=channelLow; c <=channelHigh; c++) { //+1 to reach max iterator
            foo.set_run_info(s,c); // set all attributes and save to vector - need this to access all information

            foo.set_run_channel(s,c);
            foo.set_baseline(s,c,wf_root1);
            foo.set_peaks(s,c,wf_root1, 400); // store x-positions - 497 interval
            
            foo.set_pos_peak_mean(s,c);
            foo.set_neg_peak_mean(s,c);
            foo.set_pos_peak_rms(s,c);
            foo.set_neg_peak_rms(s,c);
        }
    }
    
    //
    // wf display
    //
    
    TCanvas *canv_wf = new TCanvas("canv_wf", "canv_wf", 100, 10, 1500, 250);
    canv_wf->Divide(subrunHigh-subrunLow);
    for (int i = subrunLow; i <= subrunHigh; i++) {
        TString histName;
        
        canv_wf->cd(i+1);
        
        Int_t s = subrunLow+i;
        Int_t c = channelPoint;
        
        histName.Form("(Subrun: %d, Channel: %d)",s,c);
        TH1F *hist_wf = new TH1F(histName,"",4000,1,4000); //1000 - x axis
        
        for(int j = 0; j < wf_root1[s][c].size(); j++) hist_wf->SetBinContent(j , wf_root1[s][c].at(j)); //draw input waveform
        
        pretty_plot(hist_wf, 1, histName, "time (no units)","ADC","","",0,0,800,4096);
    }
    canv_wf->Modified();
    canv_wf->Update();
    
    //
    // subrun vs peak amplitude (for a specific channel)
    //
    
    TCanvas *canv_subrun_summary = new TCanvas("canv_subrun", "canv_subrun", 100, 10, 1000, 800);
    for (int c = channelLow ; c <= channelHigh; c++) {
        TString histName;
        
        Int_t xmin = 0;
        Int_t xmax = 10;
        
        histName.Form("Unique Name: %d", c);
        TH1F *peaks_hist = new TH1F(histName,"",(xmax-xmin),xmin,xmax);
        
        for (int s = subrunLow; s <= subrunHigh; s++) peaks_hist->SetBinContent(s,foo.get_pos_peak_mean(s,c));
        
        pretty_plot(peaks_hist, 1,histName, "food","food","add","", 0, 0, xmax, 4100);
        
        auto legend = new TLegend(0.5,0.7,0.85,0.9);
        legend->AddEntry(peaks_hist,histName);
        legend->Draw("same");
        
    }
    
    canv_subrun_summary->Modified();
    canv_subrun_summary->Update();
    
    //
    // response function & peak
    //
    
    TCanvas *canv_resp = new TCanvas("canv_resp", "canv_resp", 100, 10, 1500, 750);
    canv_resp->Divide(5,3); //!!!!
    for (int i = subrunLow; i <= subrunHigh; i++) { //!!!!
        
        // note to self!! ->loop over the peaks & save time and gain parameters cp
        
        TString histName;

        canv_resp->cd(i+1);
        
        Int_t s = subrunLow+i;
        Int_t c = channelPoint;
        Int_t peak = 2; //1st peak, 2nd peak, etc...
        
        Int_t x = foo.get_pos_peak_x(s,c,peak);
        Int_t xShift = foo.get_risetime(s,c,peak,wf_root1);
        Int_t Baseline = foo.get_baseline(s,c);
        
        histName.Form("Peak R. Fit (Subrun: %d, Channel: %d)", s, c);
        TH1F *hist_peak = new TH1F(histName,"",4000,1,4000); //histogram subrun 1 wf #1
        if (x) xShift = xShift+1;
        for(int j=0; j<50; j++) hist_peak->SetBinContent(j,(wf_root1[s][c].at(x - xShift + j) - Baseline) );
        
        TF1 *resp = new TF1("response",fudge_sundae,0,500,2);
        //resp->SetParameters(140.*1.012,foo.get_risetime(s,c,peak,wf_root1)); // 14.0mV/fC for now with 2.0 us shaping time
        resp->SetParameters(2000,foo.get_risetime(s,c,peak,wf_root1)); // 14.0mV/fC for now with 2.0 us shaping time
        hist_peak->Fit("response","q","",0,foo.get_width(s,c,peak,wf_root1)); //fit range es muy importante
        //hist_peak->Fit("response","q","",0,foo.get_width(s,c,wf_root1)); //fit range es muy importante // silenced fit results
        
        pretty_plot(hist_peak, 1, histName, "time (no units)", "ADC","add","",0,0,foo.get_width(s,c,peak,wf_root1),2100);
        
        Double_t gain = resp->GetParameter(0);
        Double_t risetime = resp->GetParameter(1);

        foo.set_fit_gain(s,c,peak,gain);
        foo.set_fit_risetime(s,c,peak,risetime);
    }
    canv_resp->Modified();
    canv_resp->Update();
    
    //
    // baseline histogram
    //
    
    TCanvas *canv_base = new TCanvas("canv_base", "canv_base", 100, 10, 1500, 250);
    canv_base->Divide(5); //!!!!
    for (int i = 0; i < 5; i++) { //!!!!
        TString histName;

        canv_base->cd(i+1);

        Int_t s = subrunLow+i;
        Int_t c = channelPoint;
        
        histName.Form("Baseline (Subrun: %d, Channel: %d)",s,c);
        TH1F *hist_baseL = new TH1F(histName,"",4096,1,4096); //histogram subrun 1 wf #1
        
        for(int j = 0; j < 10000; j++) hist_baseL->Fill(wf_root1[s][c].at(j)); //draw input waveform
        
        pretty_plot(hist_baseL, 1, histName, "ADC Bins", "Count","","",foo.get_baseline(s,c)-75,0,foo.get_baseline(s,c)+75,1000);
        
    }
    canv_base->Modified();
    canv_base->Update();
    
    //
    // Terminal Print Out Summary
    //
    
    std::cout << std::setprecision(2) << std::fixed;
    
    cout << "\n\tScript Parameters:" << endl;
    cout << "  \t==================\n" << endl;
    
    cout << "\t\t Subrun Low: " << subrunLow << endl;
    cout << "\t\t Subrun High: " << subrunHigh << endl;
    cout << "\t\t Subrun Point: " << subrunPoint << "\n" << endl;
    
    cout << "\t\t Channel Low: " << channelLow << endl;
    cout << "\t\t Channel High: " << channelHigh << endl;
    cout << "\t\t Channel Point: " << channelPoint << "\n" << endl;
    
    cout << "\n\tAnalysis Summary:" << endl;
    cout << "  \t================:\n" << endl;

    cout << "\t\tSubrun\tChannel\t\tBaseline\tPos. Peaks\tPos. Peak MEAN\tPos. Peak RMS\tNeg. Peaks\tNeg. Peak MEAN\tNeg. Peak RMS\tPeak Risetime\tPeak Width" << endl;
    cout << "\t\t======\t=======\t\t========\t==========\t==============\t=============\t==========\t==============\t=============\t=============\t==========" << endl;
    for (int s = subrunLow; s <= subrunHigh; s++) {
        for (int c = channelLow; c<= channelHigh; c++) {
            cout << "\t\t" << s << "\t" << c << "\t\t" << foo.get_baseline(s,c) << "\t\t" << foo.get_pos_npeaks(s,c) << "\t\t" << foo.get_pos_peak_mean(s,c) << "\t\t" << foo.get_pos_peak_rms(s,c) << "\t\t"
            << foo.get_neg_npeaks(s,c) << "\t\t" << foo.get_neg_peak_mean(s,c) << "\t\t" << foo.get_neg_peak_rms(s,c) << "\t\t" << foo.get_risetime(s,c,2,wf_root1) << "\t\t" << foo.get_width(s,c,2,wf_root1) << endl;
        }
    }
    
    std::cout << std::setprecision(2) << std::fixed;

    cout << "\n\tResp. Fit Summary:" << endl;
    cout << "  \t=================:\n" << endl;
    
    cout << "\t\tSubrun\tChannel\t\tRisetime\tGain" << endl;
    cout << "\t\t======\t=======\t\t========\t====" << endl;
    for (int s = subrunLow; s <= subrunHigh; s++) {
        for (int c = channelLow; c<= channelHigh; c++) {
            cout << "\t\t" << s << "\t" << c << "\t\t" << foo.get_fit_risetime(s,c,2) << "\t\t" << foo.get_fit_gain(s,c,2) << endl;
        }
    }
    
    //sprintf(titleHist,"ADC_Performance_Plots/20170816T145947/%s.png",file1);
    //c1->Print(titleHist); //save to file
    
    //
    // Show Peak Location
    //
    
    if (DBG) {
        char titleGraph[40];
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


