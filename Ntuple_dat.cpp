#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <array>
#include <iostream>
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

//std::vector<unsigned short> wf[64][128][10]; //store waveforms [subrun][channel][file]
std::vector<unsigned short> wf;

Double_t signalSizes_fpga[64] = {0.606,0.625,0.644,0.663,0.682,0.701,0.720,0.739,0.758,0.777,0.796,0.815,0.834,
    0.853,0.872,0.891,0.909,0.928,0.947,0.966,0.985,1.004,1.023,1.042,1.061,1.080,1.099,1.118,1.137,
    1.156,1.175,1.194,1.213,1.232,1.251,1.269,1.288,1.307,1.326,1.345,1.364,1.383,1.402,1.421,1.440,
    1.459,1.478, 1.497,1.516,1.535,1.554,1.573,1.592,1.611,1.629,1.648,1.667,1.686,1.705,1.724,1.743,
    1.762,1.781,1.800}; //V

void Ntuple_dat(const char* param_file ,const char* input_file){
    Analyze foo;
    //
    // read parameters
    vector<string> words  = getTimeStamp(param_file);
    std::string time_stamp = words.at(0);
    Char_t *ts = (char *) time_stamp.c_str() ;
    
    std::cout << "ts: " << ts << std::endl;
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
    Int_t file = f_p;
    Int_t temp = par.at(11);
    
    // hard set - expected shaping time
    vector<Double_t> pulse_shape; //[gain][shape]
    pulse_shape.push_back(1); // g2 s0 ( s = 1, F = 0  )
    pulse_shape.push_back(2); // g2 s1 ( s = 2, F = 1  )
    pulse_shape.push_back(4); // g2 s2 ( s = 4, F = 2  )
    pulse_shape.push_back(6); // g2 s3 ( s = 6, F = 3  )
    pulse_shape.push_back(1); // g2 s0 ( s = 1, F = 4  )
    pulse_shape.push_back(2); // g2 s1 ( s = 2, F = 5  )
    pulse_shape.push_back(4); // g2 s2 ( s = 4, F = 6  )
    pulse_shape.push_back(6); // g2 s3 ( s = 6, F = 7  )
    Int_t gain, gain_v, shape;
    Double_t shape_v;
    if (f_p==0){ gain = 2; gain_v = 14; shape = 0; shape_v = 0.5;}
    if (f_p==1){ gain = 2; gain_v = 14; shape = 1; shape_v = 1;}
    if (f_p==2){ gain = 2; gain_v = 14; shape = 2; shape_v = 2;}
    if (f_p==3){ gain = 2; gain_v = 14; shape = 3; shape_v = 4;}
    if (f_p==4){ gain = 3; gain_v = 25; shape = 0; shape_v = 0.5;}
    if (f_p==5){ gain = 3; gain_v = 25; shape = 1; shape_v = 1;}
    if (f_p==6){ gain = 3; gain_v = 25; shape = 2; shape_v = 2;}
    if (f_p==7){ gain = 3; gain_v = 25; shape = 3; shape_v = 4;}
    Double_t RT_f_fited = 0, RT_f = 0, RT_f_temp = 0, RT_k = pulse_shape.at(f_p), fpga_voltage = 0;
    Double_t GN_f_fited = 0, GN_f = 0, GN_f_temp = 0, GN_m = 0, chi_v = 1, chi_t = 0, x_s = 0;
    Int_t s = 0, c = 0, p = 0, x = 0, y = 0, bl = 0, bl_s = 0, charge = 0, bl_rms = 0, fit_range = 30;
    //
    // read root file
    TTree *tr_rawdata;
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
    Long64_t nEntries(tr_rawdata->GetEntries());          // 11 subrun * 128 channels = 1408 entries (wf)
    tr_rawdata->GetEntry(0);
    
    //
    // write root file
    cout << "" << endl;
    TFile f(Form("Data/%s/g%d_s%d_extpulse/Results.root", time_stamp.c_str(), gain, shape),"recreate");
    TTree roast_beef("roast_beef","a big tree with stuff");
    roast_beef.Branch("RT_f",&RT_f,"RT_f/D");
    roast_beef.Branch("GN_f",&GN_f,"GN_f/D");
    roast_beef.Branch("GN_m",&GN_m,"GN_m/D");
    roast_beef.Branch("charge",&charge,"charge/I");
    roast_beef.Branch("fpga_voltage",&fpga_voltage,"fpga_voltage/D");
    roast_beef.Branch("s",&s,"s/I");
    roast_beef.Branch("c",&c,"c/I");
    roast_beef.Branch("p",&p,"p/I");
    roast_beef.Branch("bl",&bl,"bl/I");
    roast_beef.Branch("bl_s",&bl_s,"bl_s/I");
    roast_beef.Branch("gain",&gain,"gain/I");
    roast_beef.Branch("shape",&shape,"shape/I");
    roast_beef.Branch("gain_v",&gain_v,"gain_v/I");
    roast_beef.Branch("shape_v",&shape_v,"shape_v/D");
    roast_beef.Branch("timestamp",&words.at(0));
    //roast_beef.Branch("temp",&temp,"temp/I");       //x
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        tr_rawdata->GetEntry(entry);
        if( subrunIn < 0 || subrunIn >= const_numSubrun ) continue;
        if( chanIn < 0 || chanIn >= const_numChan ) continue;
        wf.clear();
        for(unsigned int idx = 0 ; idx < wfIn->size() ; idx++) wf.push_back(wfIn->at(idx));
        
        foo.set_run_info(subrunIn,chanIn,f_p);    // initialize vector of structures size
        foo.set_run_channel(subrunIn,chanIn);
        foo.set_peaks(subrunIn,chanIn,f_p,shape,wf);
        foo.set_baseline(subrunIn,chanIn,f_p,wf);

        TH1F *hist_peak = new TH1F("","",50,0,49);
        if (subrunIn > 1) {
            bl = foo.get_baseline(1,chanIn,f_p);

            for (p = 0; p < foo.get_pos_npeaks(subrunIn,chanIn,f_p)-2; p++) {
                
                x = foo.get_pos_peak_x(subrunIn,chanIn,f_p,p);
                y = foo.get_pos_peak_y(subrunIn,chanIn,f_p,p);

                TF1 *resp_fit = new TF1("resp_fit",fudge_sundae,0,49,2);    // define response function
                
                // fit
                //--------------------------------------------------------------------
                //std::cout << " --- --- --- --- --- " << p;
                int i_a = 1;
                while (i_a < 20) {
                    x_s = i_a;
                    if (x_s>x) break;
                    for(int i_b = 0; i_b < 49; i_b++) {
                        if (x+i_b>wf.size()) break;
                        hist_peak->SetBinContent(i_b,(wf.at(x - x_s + i_b) - bl));
                    }
                    resp_fit->SetParameters(2500,x_s);                              // set parameters
                    hist_peak->Fit("resp_fit","oq","",0, fit_range); // set fit
                    RT_f_fited = std::abs(resp_fit->GetParameter(1));               // set fit data
                    GN_f_fited = std::abs(resp_fit->GetParameter(0)/10);
                    chi_v = resp_fit->GetChisquare();

                    if( (std::abs(chi_t)/std::abs(chi_v))>1 || i_a==1 ) {
                        RT_f_temp = RT_f_fited;
                        GN_f_temp = GN_f_fited;
                        chi_t = chi_v;
                    }
                    
                    i_a++;
                }
                //--------------------------------------------------------------------
                std::cout << " --- --- --- --- --- final results" << std::endl;
                std::cout << " s: " << subrunIn << " c: " << chanIn << " p: " << p << " x: " << x;
                std::cout << " GN_f_temp: " << GN_f_temp << " y: " << y-bl << std::endl;
                
                RT_f = RT_f_temp;
                GN_f = GN_f_temp;
                GN_m = y-bl;
                charge = 0; if(subrunIn>1) charge = ((signalSizes_fpga[subrunIn-1]-signalSizes_fpga[0])*184*10000/1.6); // assumes s>=2
                fpga_voltage = 0; if(subrunIn>1) fpga_voltage = (signalSizes_fpga[subrunIn-1]-signalSizes_fpga[0]);
                s = subrunIn;
                c = chanIn;
                p = p;
                gain = gain;
                shape = shape;
                gain_v = gain_v;
                shape_v = shape_v;
                ts = ts;
                bl = bl;
                bl_s = bl_s;
                roast_beef.Fill();
            }
        }
        //cout << "subrun: " << subrunIn << " channel: " << chanIn << endl;
    }
    inputFile->Close();
    roast_beef.Write();
}

