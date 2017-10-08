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
using namespace std;

unsigned short subrunIn1, chanIn1;    //output tree and variable
std::vector<unsigned short> *wfIn1;
TTree *tr_rawdata1;                  //ROOT TREE tr_rawdata variables
TFile *inputFile1;

const int const_numSubrun = 64;
const int const_numChan = 128;

//------------ Control ------------
int subrunLow=2, subrunHigh=10, subrunPoint = 8; //controls peak finding limits (subrun and channel) - need baseline info (mean and rms) - found from subrun 1
int channelLow=125, channelHigh=127, channelPoint = 127;
int subrunBaseline=1;
const char file1[100] = "20170816T145947_fembTest_gainenc_test_g3_s0_extpulse.root";
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
    cout << "max subrun: " << subrunIn1 << " max channel: " << chanIn1 << endl;
    inputFile1->Close();
    
    //-------------------------------------
    //----   3. set peak information   ----
    //-------------------------------------
    Analyze foo;
    for (int s=1; s <=subrunHigh; s++) { //+1 to reach max iterator
        for (int c=channelLow; c <=channelHigh; c++) { //+1 to reach max iterator
            foo.set_run_info(s,c); // set all attributes and save to vector - need this to access all information
        }
    }
    for (int s=subrunLow; s <=subrunHigh; s++) { //+1 to reach max iterator
        for (int c=channelLow; c <=channelHigh; c++) { //+1 to reach max iterator
            foo.set_run_channel(s,c);
            foo.set_baseline(s,c,wf_root1);
            foo.set_peaks(s,c,wf_root1, 400); // store x-positions - 497 interval
            
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
    pp->Fit("m1","qpol1","N",subrunLow,10);
    g1->SetLineColorAlpha(kBlue,0.35);
    g1->SetLineStyle(2);
    g1->Draw("same");
    
    g1->GetParameters(&par1[0]);
    
    // fit neg peak
    TF1 *g2 = new TF1("m2","pol1",subrunLow,10);
    Double_t par2[9];
    np->Fit("m2","qpol1","N",subrunLow,10);
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
        
        h1->Fit("func3","q","",0,foo.get_width(9,channelPoint, wf_root1)); //fit range* muy importante

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
        TH1F *h3 = new TH1F(histName,"q",1000,0,1000); //histogram subrun 1 wf #1
        
        for (Int_t i=0;i!=999;i++) h3->SetBinContent(i, f3->Eval((i+1)));

        pretty_plot(h3, 1, histName, "time (?? unknown units)", "ADC",0,0,0,25,2000);
        
    }
    
    c1->Modified();
    c1->Update();
    
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

    cout << "\t\tSubrun\tChannel\tBaseline\tPos. Peaks\tPos. Peak MEAN\tPos. Peak RMS\tNeg. Peaks\tNeg. Peak MEAN\tNeg. Peak RMS\tPeak Risetime\tPeak Width" << endl;
    cout << "\t\t======\t=======\t========\t==========\t==============\t=============\t==========\t==============\t=============\t=============\t==========" << endl;
    for (int s = subrunLow; s <= subrunHigh; s++) {
        for (int c = channelLow; c<= channelHigh; c++) {
            cout << "\t\t" << s << "\t" << c << "\t" << foo.get_baseline(s,c) << "\t\t" << foo.get_pos_npeaks(s,c) << "\t\t" << foo.get_pos_peak_mean(s,c) << "\t\t" << foo.get_pos_peak_rms(s,c) << "\t\t"
            << foo.get_neg_npeaks(s,c) << "\t\t" << foo.get_neg_peak_mean(s,c) << "\t\t" << foo.get_neg_peak_rms(s,c) << "\t\t" << foo.get_risetime(s,c,wf_root1) << endl;//<< "\t\t" << foo.get_width(s,c,wf_root1) << endl;
        }
    }

    
    
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


