#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>
#include "TTree.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TROOT.h"

#include "Headers/analyze.h"
#include "Headers/management.h"

using namespace std;

unsigned short subrunIn, chanIn;    //output tree and variable
std::vector<unsigned short> *wfIn;
//std::vector<unsigned short> wfOut;

const int const_numSubrun = 64;
const int const_numChan = 128;

std::vector<unsigned short> wf[64][128][10]; //store waveforms [subrun][channel][file]

Double_t signalSizes_fpga[64] = {0.606,0.625,0.644,0.663,0.682,0.701,0.720,0.739,0.758,0.777,0.796,0.815,0.834,
    0.853,0.872,0.891,0.909,0.928,0.947,0.966,0.985,1.004,1.023,1.042,1.061,1.080,1.099,1.118,1.137,
    1.156,1.175,1.194,1.213,1.232,1.251,1.269,1.288,1.307,1.326,1.345,1.364,1.383,1.402,1.421,1.440,
    1.459,1.478, 1.497,1.516,1.535,1.554,1.573,1.592,1.611,1.629,1.648,1.667,1.686,1.705,1.724,1.743,
    1.762,1.781,1.800};


void Ntuple_dat(const char* param_file ,const char* input_file){
    
    Analyze foo;

    // hardset color
    vector<Int_t> color;
    color.push_back(1); //kBlack 1
    color.push_back(432); //kCyan 2
    color.push_back(800); //kBlue 3
    color.push_back(612); //kMagenta 4
    color.push_back(840); //kTeal 5
    color.push_back(632); //kRed 6
    color.push_back(800); //kOrange 7
    color.push_back(416); //kSpring 8
    color.push_back(920); //Grey 9
    color.push_back(880); //Violet 10
    
    // hard set - expected shaping time
    vector<Double_t> pulse_shape; //[gain][shape]
    pulse_shape.push_back(1); // g2 s0 ( s = 1, F = 0 )
    pulse_shape.push_back(2); // g2 s1 ( s = 2, F = 1  )
    pulse_shape.push_back(4); // g2 s2 ( s = 4, F = 2  )
    pulse_shape.push_back(6); // g2 s3 ( s = 6, F = 3  )
    pulse_shape.push_back(1); // g2 s0 ( s = 1, F = 4  )
    pulse_shape.push_back(2); // g2 s1 ( s = 2, F = 5  )
    pulse_shape.push_back(4); // g2 s2 ( s = 4, F = 6  )
    pulse_shape.push_back(6); // g2 s3 ( s = 6, F = 7  )
    
    //
    // read parameters
    //
    vector<string> words  = getTimeStamp(param_file);
    string time_stamp = words.at(0);
    
    vector<double> par = getParameters(param_file);
    int s_l = par.at(1);
    int s_h = par.at(2);
    int s_p = par.at(3);
    
    int c_l = par.at(4);
    int c_h = par.at(5);
    int c_p = par.at(6);
    
    int p_l = par.at(7);
    int p_h = par.at(8);
    int p_p = par.at(9);
    
    vector<string> file_vector = getFiles(input_file);
    int f_l = 0;
    int f_h = file_vector.size();
    int f_p = par.at(10);

    //
    // read root file
    //
    TTree *tr_rawdata;
    //TFile *inputFile;
    string file_name = file_vector.at(f_p); // F_P => Reads only one file
    
    TFile *inputFile = new TFile(file_name.c_str(), "READ");
    tr_rawdata = (TTree*) inputFile->Get("femb_wfdata");
    if( !tr_rawdata ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    
    tr_rawdata->SetBranchAddress("subrun", &subrunIn);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("chan", &chanIn);        // initialize channel branch
    tr_rawdata->SetBranchAddress("wf", &wfIn);            // initialize waveform branch
    
    Long64_t nEntries(tr_rawdata->GetEntries());          // 11 subrun * 128 channels = 1408
    tr_rawdata->GetEntry(0);
    
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        tr_rawdata->GetEntry(entry);
        if( subrunIn < 0 || subrunIn >= const_numSubrun ) continue;
        if( chanIn < 0 || chanIn >= const_numChan ) continue;
        for( unsigned int s = 0 ; s < wfIn->size() ; s++ ){//store waveform vector in array for quick access
            wf[subrunIn][chanIn][f_p].push_back( wfIn->at(s) );
        }
    }
    
    inputFile->Close();
    
    //
    // initialize peak information
    //
    
    for (int s=s_l; s<=s_h; s++) { // subrun loop
        for (int c=c_l; c<=c_h; c++) { // channel loop
            foo.set_run_info(s,c,f_p);    // initialize vector of structures size
            
            foo.set_run_channel(s,c);   // initialize vector of structures size
            foo.set_baseline(s,c,f_p,wf);
            foo.set_baseline_rms(s,c,f_p,wf);
            foo.set_peaks(s,c,f_p,wf, 400);
            
            foo.set_pos_peak_mean(s,c,f_p);
            foo.set_neg_peak_mean(s,c,f_p);
            foo.set_pos_peak_rms(s,c,f_p);
            foo.set_neg_peak_rms(s,c,f_p);
        }
    }
    
    //
    // fit to peaks
    //
    
    auto legend = new TLegend(0.1,0.7,0.95,0.9);
    TString legend_title;
    
    Double_t RT_f = 0; //f = fit
    Double_t RT_k = 0; //k = known
    Double_t GN_f = 0;
    Double_t GN_m = 0; //m = measured
    Double_t GN_r = 0; //r = ratio
    Double_t fit_range = 0;
    Int_t s = 0;
    Int_t c = 0;
    Int_t p = 0;
    Int_t x = 0;
    Int_t x_s = 0;
    Int_t bl = 0;
    Int_t charge = 0;
    Int_t file = f_p;
    Double_t bl_rms = 0;
    
    TFile f(Form("Results/%s/Tree_%s_F%d.root",time_stamp.c_str(),time_stamp.c_str(),f_p),"recreate");
    TTree roast_beef("roast_beef","a big tree with stuff");
    roast_beef.Branch("RT_f",&RT_f,"RT_f/D");
    roast_beef.Branch("RT_k",&RT_k,"RT_k/D");
    roast_beef.Branch("GN_f",&GN_f,"GN_f/D");
    roast_beef.Branch("GN_m",&GN_m,"GN_m/D");
    roast_beef.Branch("GN_r",&GN_r,"GN_r/D");
    roast_beef.Branch("charge",&charge,"charge/I");
    roast_beef.Branch("file",&file,"file/I");
    roast_beef.Branch("s",&s,"s/I");
    roast_beef.Branch("c",&c,"c/I");
    roast_beef.Branch("p",&p,"p/I");
    roast_beef.Branch("bl",&bl,"bl/I");
    roast_beef.Branch("bl_rms",&bl_rms,"bl_rms/D");
    //roast_beef.Branch("wfOut",&wfOut);
    for (c = c_l; c <= c_h; c++) {
        for (s = s_l; s <= s_h; s++) {
            //wfOut.clear();
            TH1F *hist_peak = new TH1F("","",4000,1,4000);
            for (p = p_l; p < foo.get_pos_npeaks(s,c,f_p); p++) {
                RT_k = pulse_shape.at(f_p);
                RT_f = pulse_shape.at(f_p)-1;
                bl = foo.get_baseline(s,c,f_p);
                bl_rms = foo.get_baseline_rms(s,c,f_p);
                charge = ((signalSizes_fpga[s]-signalSizes_fpga[0])*184*10000/1.6);
                Int_t idx = 0;
                while (std::abs(RT_f-RT_k)>0.4) {
                    x = foo.get_pos_peak_x(s,c,f_p,p);
                    x_s = idx;
                    
                    //if(x_s > x) x_s = 0;
                    for(int j = 0; j < 100; j++) hist_peak->SetBinContent(j,(wf[s][c][f_p].at(x - x_s + j) - bl));
                    
                    // set fit range
                    fit_range = 30;//pulse_width.at(f_p);
                    
                    // define response function
                    TF1 *resp = new TF1("response",fudge_sundae,0,500,2);
                    
                    // set parameters
                    resp->SetParameters(3000,x_s);
                    
                    // set fit
                    hist_peak->Fit("response","oq","",0, fit_range);
                    
                    // set fit data
                    RT_f = std::abs(resp->GetParameter(1));
                    GN_f = std::abs(resp->GetParameter(0)/10);
                    
                    idx++;
                    if (idx>12) break;
                }
                //for (int idx = 0; idx < wf[s][c][f_p].size(); idx++) wfOut.push_back(wf[s][c][f_p].at(idx));
                GN_m = foo.get_pos_peak_y(s,c,f_p,p);
                GN_r = GN_f/GN_m;
                //cout << " b_rms: " << b_rms << " RT_f: " << RT_f;
                //cout << " RT: " << RT_k << " GN_f: " << GN_f << " GN_m: " << GN_m;
                //cout << "\t| s: " << s << " c: " << c << " p: " << p;
                //cout << " f: " << f_p << endl;
                roast_beef.Fill();
            }
            Int_t l_s = 1;
            Int_t l_t = 1;
        }
        
    }
    roast_beef.Write();
}

