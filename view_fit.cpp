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
TTree *tr_rawdata;                  //ROOT TREE tr_rawdata variables
TFile *inputFile;

const int const_numSubrun = 64;
const int const_numChan = 128;

//std::vector<unsigned short> wf[64][128][10]; //store waveforms [subrun][channel][file]
std::vector<unsigned short> wf; //store waveforms [subrun][channel][file]

//----------------------------------------------------------------------------
//----                                main                                ----
//----------------------------------------------------------------------------
void view_fit(const char* param_file ,const char* input_file){
    
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
    Int_t temp = par.at(11);
    
    //
    // read root file
    //
    
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
        if( chanIn < 0 || chanIn >= const_numChan ) continue; wf.clear();
        for( unsigned int s = 0 ; s < wfIn->size() ; s++ ){//store waveform vector in array for quick access
            wf.push_back( wfIn->at(s) );
        }
        foo.set_run_info(subrunIn,chanIn,f_p);    // initialize vector of structures size
        foo.set_run_channel(subrunIn,chanIn);   // initialize vector of structures size
        foo.set_baseline(subrunIn,chanIn,f_p,wf);
        foo.set_baseline_rms(subrunIn,chanIn,f_p,wf);
        foo.set_peaks(subrunIn,chanIn,f_p,wf, 400);
    }
    inputFile->Close();
    // fit method 1.
    TCanvas *canv_resp_fix = new TCanvas("canv_resp", "canv_resp", 100, 10, 600, 750);
    canv_resp_fix->Divide(2,4);
    for (s = s_l; s <= s_h; s++) {
        inc++;
        canv_resp_fix->cd(inc);
        
        bl = foo.get_baseline(s,c_p,f_p);
        bl_rms = foo.get_baseline_rms(s,c_p,f_p);
        
        x = foo.get_pos_peak_x(s,c_p,f_p,p_p);
        y = foo.get_pos_peak_y(s,c_p,f_p,p_p);
        
        //------------------------------------------------------
        x_s = pulse_shape.at(f_p);
        TH1F *hist_peak_fix = new TH1F("","",50,0,49);
        for(int j = 0; j < 100; j++) hist_peak_fix->SetBinContent(j,(wf[s][c_p][f_p].at(x - x_s + j) - bl));
        // set fit range
        fit_range = (Int_t) 30;//pulse_width.at(f_p);
        // define response function
        TF1 *resp_fix = new TF1("response_fixed",fudge_sundae,0,49,2);
        // set parameters
        resp_fix->SetParameters(2500,x_s);
        resp_fix->FixParameter(1,pulse_shape.at(f_p));
        // set fit
        hist_peak_fix->Fit("response_fixed","q","",0, fit_range);
        // set fit data
        RT_f_fixed = std::abs(resp_fix->GetParameter(1));
        GN_f_fixed = std::abs(resp_fix->GetParameter(0))/10;
        chi_fixed = resp_fix->GetChisquare();
        //------------------------------------------------------
        GN_m = foo.get_pos_peak_y(s,c_p,f_p,p_p);
        cout << " x: " << x << " y: " << y << " x_s: " << x_s;
        cout << "\t\t chi2 fix: " << chi_fixed << " RT_fix: " << RT_f_fixed << " GN_fixed " << GN_f_fixed;
        cout << "\t | GN_m: "     << GN_m      << " s: "      << s << " c: " << c_p << " p: " << p_p << endl;
    }
    canv_resp_fix->Modified();
    canv_resp_fix->Update();
    
    
    // fit method 2.
    TCanvas *canv_resp_fit = new TCanvas("canv_resp_fit", "canv_resp_fit", 100, 10, 600, 750);
    canv_resp_fit->Divide(2,4);
    inc = 0;
    for (s = s_l; s <= s_h; s++) {
        inc++;
        canv_resp_fit->cd(inc);

        RT_f_fited = 0;
        RT_f = 0;
        Int_t idx = 0;
        
        //------------------------------------------------------
        // fit process
        TH1F *hist_peak = new TH1F("","",50,0,49);
        while (idx < 10) {
            x = foo.get_pos_peak_x(s,c_p,f_p,p_p);
            x_s = idx;
            TH1F *hist_peak_fit = new TH1F("","",50,0,49);
            for(int j = 0; j < 100; j++) hist_peak_fit->SetBinContent(j,(wf[s][c_p][f_p].at(x - x_s + j) - bl));
            // set fit range
            fit_range = (Int_t) 30;//pulse_width.at(f_p);
            // define response function
            TF1 *resp_fit = new TF1("response_fit",fudge_sundae,0,49,2);
            // set parameters
            resp_fit->SetParameters(2500,x_s);
            // set fit
            hist_peak_fit->Fit("response_fit","oq","",0, fit_range);
            // set fit data
            RT_f_fited = std::abs(resp_fit->GetParameter(1));
            GN_f_fited = std::abs(resp_fit->GetParameter(0)/10);
            if ( std::abs(RT_f_fited-RT_k)<std::abs(RT_f-RT_k) ){
                RT_f = RT_f_fited;
                GN_f = GN_f_fited;
                chi_fited = resp_fit->GetChisquare();
                for(int j = 0; j < 100; j++) hist_peak->SetBinContent(j,(wf[s][c_p][f_p].at(x - x_s + j) - bl));;
            }
            idx++;
        }
        //------------------------------------------------------
        GN_m = foo.get_pos_peak_y(s,c_p,f_p,p_p);
        cout << " x: " << x << " y: " << y << " x_s: " << x_s;
        cout << "\t\t chi2 fit: " << chi_fited << " RT_fit: " << RT_f << " GN_fited " << GN_f;
        cout << "\t | GN_m: "     << GN_m      << " s: "      << s << " c: " << c_p << " p: " << p_p << endl;
        
        pretty_plot(hist_peak, color.at(2), 1, 1, Form("c %d",c_p), "Bin", "ADC","","",0,0,20,2800);
    }
    canv_resp_fit->Modified();
    canv_resp_fit->Update();

    // image out to directory
    TImage *img_resp = TImage::Create();
    img_resp->FromPad(canv_resp_fit);
    TString output_file_name;
    output_file_name.Form("Results/%s/%s_S%d_C%d_F%d_P%d_INL_fit.png",time_stamp.c_str(),time_stamp.c_str(),s_p,c_p,f_p,p_p);
    img_resp->WriteImage(output_file_name);
    
}

