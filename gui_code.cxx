#include <TApplication.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TGWindow.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include <TSpectrum.h>
#include "TROOT.h"
#include "TMath.h"
#include "TLegend.h"
#include "dec_def.h"
#include "Headers/analyze.h"

std::vector<unsigned short> *wfIn;
std::vector<unsigned short> wf;

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

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h){
    // create a main frame
    fMain = new TGMainFrame(p,w,h);

    //--------------------------------------//
    //          horizontal level 0          //
    //--------------------------------------//
    hfm_0 = new TGHorizontalFrame(fMain, w, h);
    
    tbuf1 = new TGTextBuffer(50);
    tent1 = new TGTextEntry(hfm_0, tbuf1, Tentry1);
    tent1->Resize(w*2, 20);
    tent1->SetText("Data/20170808T143522/g2_s0_extpulse/Results.root");
    
    hfm_0->AddFrame(tent1, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 3, 4));
    
    //--------------------------------------//
    //          horizontal level 0_0        //
    //--------------------------------------//
    hfm_0_0 = new TGHorizontalFrame(fMain, w, h);

    tbuf2 = new TGTextBuffer(50);
    tent2 = new TGTextEntry(hfm_0_0, tbuf2, Tentry2);
    tent2->Resize(w*2, 20);
    tent2->SetText("Data/20170808T143522/g2_s0_extpulse/gainMeasurement_femb_1-parseBinaryFile.root");

    hfm_0_0->AddFrame(tent2, new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 5, 5, 3, 4));
    //--------------------------------------//
    //          horizontal level 1          //
    //--------------------------------------//
    hfm_1 = new TGHorizontalFrame(fMain, w, h);
    //top = new TGCompositeFrame(hfm_1,10,10, kSunkenFrame);
    
    //          left side: frame            //
    vf_1 = new TGVerticalFrame(hfm_1, 10, 10);
    hf_1 = new TGHorizontalFrame(vf_1, 10, 10); //internal to left vertical
    l_frame = new TGCompositeFrame(vf_1, 10, 10, kSunkenFrame);
    
    vf_1->AddFrame(hf_1,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,1,2));
    vf_1->AddFrame(l_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandY | kLHintsExpandX,0,0,1,2));

    hfm_1->AddFrame(vf_1, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    
    // left side: create canvas //* try moving this into DoDrawFEMB() -Cp
    canvas_femb = new TRootEmbeddedCanvas("canvas_femb", l_frame, w, 0.5*h);     // create canvas widget
    l_frame->AddFrame(canvas_femb, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1) );
    
    // left side: draw button
    TGTextButton *draw = new TGTextButton(hf_1, "&Draw Summary");
    draw->Connect("Clicked()","MyMainFrame",this,"DoDrawFEMB()");
    hf_1->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    //          right side: frame           //
    vf_2 = new TGVerticalFrame(hfm_1, 10, 10);
    hf_2 = new TGHorizontalFrame(vf_2, 10, 10); //internal to right vertical
    r_frame = new TGCompositeFrame(vf_2, 10, 10, kSunkenFrame);

    vf_2->AddFrame(hf_2,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,1,2));
    vf_2->AddFrame(r_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandY | kLHintsExpandX,0,0,1,2));

    hfm_1->AddFrame(vf_2, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    
    // right side: create canvas
    canvas_gain = new TRootEmbeddedCanvas("canvas_gain", r_frame, w, 0.5*h);     // create canvas widget
    r_frame->AddFrame(canvas_gain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1) );
    
    // right side: slider control
    hslider = new TGHSlider(hf_2,150,kSlider1 | kScaleDownRight,HSId1);
    //hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoSlider(Int_t)");
    hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDrawGain(Int_t)");
    hslider->SetRange(0,127);
    hslider->SetPosition(0);
    hf_2->AddFrame(hslider, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

    //--------------------------------------//
    //          horizontal level 2          //
    //--------------------------------------//
    hfm_2 = new TGHorizontalFrame(fMain, w, h);
    
    // create canvas
    canvas_wf = new TRootEmbeddedCanvas("canvas_wf", hfm_2, w*2, 0.5*h);     // create canvas widget
    hfm_2->AddFrame(canvas_wf, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX));
    
    //--------------------------------------//
    //          horizontal level 3          //
    //--------------------------------------//
    hfm_3 = new TGHorizontalFrame(fMain,w,h);
    
    // number entry
    Nent_chan = new TGNumberEntry(hfm_3, 1, 6, kENTRY1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,127);
    Nent_chan->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawFit(Long_t)"); //??
    
    Nent_sub = new TGNumberEntry(hfm_3, 1, 6, kENTRY2, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,10);
    Nent_sub->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawFit(Long_t)");
    
    Nent_peak = new TGNumberEntry(hfm_3, 1, 6, kENTRY3, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,10);
    Nent_peak->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawFit(Long_t)");
    
    // add all widgets to frame
    hfm_3->AddFrame(Nent_chan, new TGLayoutHints(kLHintsBottom | kLHintsLeft, 0, 0, 0, 0));
    hfm_3->AddFrame(Nent_sub, new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 0, 0));
    hfm_3->AddFrame(Nent_peak, new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 0, 0));
    
    //--------------------------------------//
    //          horizontal level 4          //
    //--------------------------------------//
    hfm_4 = new TGHorizontalFrame(fMain,w,h);

    // exit button
    TGTextButton *exit = new TGTextButton(hfm_4, "&Exit", "gApplication->Terminate(0)");
    hfm_4->AddFrame(exit, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 3, 4));
    
    //          fmain stuff         //

    // add children to parent frame
    fMain->AddFrame(hfm_0, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    fMain->AddFrame(hfm_0_0, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    fMain->AddFrame(hfm_1, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    fMain->AddFrame(hfm_2, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    fMain->AddFrame(hfm_3, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    fMain->AddFrame(hfm_4, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    
    // set a name to the main frame
    fMain->SetWindowName("Scroll and roll");
    
    // map all subwindows of main frame
    fMain->MapSubwindows();
    
    // initialize the layout algorithm
    fMain->Resize(fMain->GetDefaultSize());
    
    // map main frame
    fMain->MapWindow();
    
}
void MyMainFrame::DoDrawGain(Int_t pos){
    Double_t RT_f, GN_f, GN_m;
    Int_t charge, s, c, p, gain, shape;
    
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id = te->WidgetId();
    //std::cout << " DoDrawGain::text entry1: " << tbuf1->GetString() << std::endl;

    TTree *tr_rawdata;
    std::string file_name = tbuf1->GetString(); // F_P => Reads only one file
    TFile *inputFile = new TFile(file_name.c_str(), "READ");
    tr_rawdata = (TTree*) inputFile->Get("roast_beef");
    if( !tr_rawdata ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    tr_rawdata->SetBranchAddress("RT_f", &RT_f);        // initialize subrun branch
    tr_rawdata->SetBranchAddress("GN_f", &GN_f);        // initialize subrun branch
    tr_rawdata->SetBranchAddress("GN_m", &GN_m);        // initialize subrun branch
    tr_rawdata->SetBranchAddress("charge", &charge);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("s", &s);              // initialize subrun branch
    tr_rawdata->SetBranchAddress("c", &c);              // initialize subrun branch
    tr_rawdata->SetBranchAddress("p", &p);              // initialize subrun branch
    tr_rawdata->SetBranchAddress("gain", &gain);        // initialize subrun branch
    tr_rawdata->SetBranchAddress("shape", &shape);      // initialize subrun branch
    
    Long64_t nEntries(tr_rawdata->GetEntries());
    tr_rawdata->GetEntry(0);
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        if (c == pos) { gain = gain; }
    }
    
    hprof = new TProfile("hprof",Form("gain: [%d] shape: [%d] channel: [%d]", gain, shape, pos),10,0,218200);
    if(gain == 2) hprof->SetMaximum(25);
    if(gain == 3) hprof->SetMaximum(40);
    tr_rawdata->Draw("(0.33*6241)*(GN_f/charge):charge>>hprof",Form("c==%d",pos),"*");
    hprof->SetLineColor(1);
    hprof->SetXTitle("Channel");
    hprof->SetYTitle("Channel");
    
    hprof_m = new TProfile("hprof_m",Form("gain: [%d] shape: [%d] channel: [%d]", gain, shape, pos),10,0,218200);
    tr_rawdata->Draw("(0.33*6241)*(GN_m/charge):charge>>hprof_m",Form("c==%d",pos),"same");
    hprof_m->SetLineColor(4);

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry("hprof"  ," Fit results ","l");
    legend->AddEntry("hprof_m"," Measurement results ","lep");
    legend->Draw();
    
    TCanvas *fCanvas = canvas_gain->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
}
void MyMainFrame::DoDrawFEMB(){
    Double_t RT_f, GN_f, GN_m;
    Int_t charge, s, c, p;
    Char_t *pesto_pasta;
    
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id = te->WidgetId();
    
    TTree *tr_rawdata;
    std::string file_name = tbuf1->GetString(); // F_P => Reads only one file
    TFile *inputFile = new TFile(file_name.c_str(), "READ");
    tr_rawdata = (TTree*) inputFile->Get("roast_beef");
    if( !tr_rawdata ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    tr_rawdata->SetBranchAddress("RT_f", &RT_f);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("GN_f", &GN_f);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("GN_m", &GN_m);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("charge", &charge);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("s", &s);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("c", &c);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("p", &p);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("timestamp", &pesto_pasta);
    
    h2 = new TH2F("h2",Form("FEMB Summary %s", "20170808T143522"),128,0,128,10,0,218200);
    tr_rawdata->Draw("charge:c>>h2","(0.33*6241)*GN_f/charge*(p==0 && GN_f<1e6)","colz");
    h2->SetStats(0);
    h2->SetXTitle("Channel");
    h2->SetYTitle("Injected Charge (e-)");
    
    TCanvas *fCanvas_x = canvas_femb->GetCanvas();
    fCanvas_x->cd();
    fCanvas_x->Update();
}
void MyMainFrame::DoDrawFit(Long_t pos){
    TGNumberEntry *sender = (TGNumberEntry *) gTQSender;
    Int_t id = sender->WidgetId();
    
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id_t = te->WidgetId();
    
    wf.clear();
    const int const_numSubrun = 64, const_numChan = 128;
    unsigned short subrunIn, chanIn;    //output tree and variable
    
    TString food = tbuf2->GetString();
    TFile *inputFile = new TFile(food, "READ");
    
    TTree *tr_rawdata;
    tr_rawdata = (TTree*) inputFile->Get("femb_wfdata");
    if( !tr_rawdata ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    tr_rawdata->SetBranchAddress("subrun", &subrunIn);      // initialize subrun branch
    tr_rawdata->SetBranchAddress("chan", &chanIn);          // initialize channel branch
    tr_rawdata->SetBranchAddress("wf", &wfIn);              // initialize waveform branch
    
    Long64_t nEntries(tr_rawdata->GetEntries());            // 11 subrun * 128 channels = 1408 entries (wf)
    tr_rawdata->GetEntry(0);
    //
    Analyze foo;
    Int_t f_p = 3, x_s = 0, x = 0, y = 0, bl = 0, fit_range = 30, scale = 2;
    Double_t RT_f = 0, RT_f_fited = 0, GN_f_fited = 0, RT_f_temp = 0, GN_f_temp = 0, chi_v, chi_t = 0, i_T = 0;
    TH1F *hist_peak = new TH1F("","",50,0,49);
    TF1 *resp_fit = new TF1("resp_fit",fudge_sundae,0,49,2); // define response function
    
    Int_t s = Nent_sub->GetIntNumber();
    Int_t c = Nent_chan->GetIntNumber();
    Int_t p = Nent_peak->GetIntNumber();
    
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        tr_rawdata->GetEntry(entry);
        if( subrunIn < 0 || subrunIn >= const_numSubrun ) continue;
        if( chanIn < 0 || chanIn >= const_numChan ) continue;
        
        if( chanIn == c && subrunIn == 1){
            for( unsigned int s = 0 ; s < wfIn->size() ; s++) wf.push_back( wfIn->at(s) );
            foo.set_run_info(subrunIn,chanIn,f_p);          // initialize vector of structures size
            foo.set_peaks(subrunIn,chanIn,f_p,wf, 400);
            foo.set_baseline(subrunIn,chanIn,f_p,wf);
            bl = foo.get_baseline(1,chanIn,f_p);
        }
        
        if( chanIn == c && subrunIn == s){
            for( unsigned int s = 0 ; s < wfIn->size() ; s++) wf.push_back( wfIn->at(s) );
            foo.set_run_info(subrunIn,chanIn,f_p);          // initialize vector of structures size
            foo.set_run_channel(subrunIn,chanIn);
            foo.set_peaks(subrunIn,chanIn,f_p,wf, 400);
            
            x = foo.get_pos_peak_x(subrunIn,chanIn,f_p,p);
            y = foo.get_pos_peak_y(subrunIn,chanIn,f_p,p)-bl;
            
            //--------------------------------------------------------------------
            int i_a = 1;
            while (i_a < 20) {                 // fit process
                x_s = i_a;
                if (x_s>x) break;
                for(int i_b = 0; i_b < 70; i_b++) hist_peak->SetBinContent(i_b,(wf.at(x - x_s + i_b) - bl));
                resp_fit->SetParameters(2500,x_s);                              // set parameters
                hist_peak->Fit("resp_fit","oq","",0, fit_range); // set fit
                RT_f_fited = std::abs(resp_fit->GetParameter(1));               // set fit data
                GN_f_fited = std::abs(resp_fit->GetParameter(0)/10);
                chi_v = resp_fit->GetChisquare();
                
                std::cout << " ch_v: " << chi_v << " ch_t: " << chi_t << " RT_f: " << RT_f_temp;
                std::cout << " GN_f: " << GN_f_fited << " GN_m: " << y << std::endl;
                if((std::abs(chi_t)/std::abs(chi_v))>1 || i_a==1) {
                    RT_f_temp = RT_f_fited;
                    GN_f_temp = GN_f_fited;
                    chi_t = chi_v;
                    i_T = i_a;
                }
                i_a++;
            }
            //--------------------------------------------------------------------
        }
    }
    
    std::cout << "-----------------final results" << std::endl;
    std::cout << " RT_temp: " << RT_f_temp << " GN_temp: " << GN_f_temp << " GN_m: " << y << std::endl;
    
    TH1F *hist = new TH1F(Form("subrun: [%d] channel: [%d] peak: [%d]", s, c, p),"",50,0,49);
    for (Int_t idx = 0; idx < 50; idx++) {
        if (i_T>x) i_T = 0;
        hist->SetBinContent(idx,wf.at(x-i_T+idx)-bl);
    }
    hist->Draw();

    resp_fit->SetParameters(3000,RT_f_temp);                              // set parameters
    hist->Fit("resp_fit","oq","",0, fit_range); // set fit
    resp_fit->Draw("same");
    
    TCanvas *fCanvas_y = canvas_wf->GetCanvas();
    //fCanvas_y->SetTitle (Form("subrun: [%d] channel: [%d] peak: [%d]", s, c, p));
    fCanvas_y->cd();
    fCanvas_y->Update();
    
    
}
void MyMainFrame::DoDrawWf(Long_t pos){
    TGNumberEntry *sender = (TGNumberEntry *) gTQSender;
    Int_t id = sender->WidgetId();
    
    wf.clear();
    const int const_numSubrun = 64, const_numChan = 128;
    unsigned short subrunIn, chanIn;    //output tree and variable
    
    TString food = Form("Data/20170808T143522/g2_s0_extpulse/gainMeasurement_femb_1-parseBinaryFile.root");
    TFile *inputFile = new TFile(food, "READ");
    
    TTree *tr_rawdata;
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
    
    TH1F *h1 = new TH1F("h1","",19498,0,19498);
    
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        tr_rawdata->GetEntry(entry);
        if( subrunIn < 0 || subrunIn >= const_numSubrun ) continue;
        if( chanIn < 0 || chanIn >= const_numChan ) continue;
        if( chanIn == Nent_chan->GetIntNumber() && subrunIn == Nent_sub->GetIntNumber()){
            
            //--------------------------------------------------------------------
            for( unsigned int s = 0 ; s < wfIn->size()-1 ; s++ ) wf.push_back(wfIn->at(s));
            TH1F *h1 = new TH1F("h1",Form("Waveform Plot Channel: %d Subrun: %d", chanIn, subrunIn),wf.size(),1,wf.size());
            for (int idx = 0; idx<wf.size(); idx++) h1->SetBinContent(idx, wf.at(idx));
            
            TH1F *d = new TH1F("d","",wf.size(),0,wf.size());
            TH1F *e = new TH1F("e","",wf.size(),0,wf.size());
            TH1F *f = new TH1F("f","",wf.size(),0,wf.size());
            
            Double_t *source = new Double_t[wf.size()];
            TSpectrum *soup = new TSpectrum(80);
            
            for (Long_t idx = 0; idx < wf.size(); idx++) source[idx]=wf.at(idx);
            soup->Background(source,wf.size(),32,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kTRUE, TSpectrum::kBackSmoothing9,kTRUE);
            for (Long_t idx = 0; idx < wf.size()-1; idx++) d->SetBinContent(idx,source[idx]);
            for (Long_t idx = 0; idx < wf.size()-1; idx++) e->Fill(d->GetBinContent(idx));
            
            Int_t nfound = soup->Search(h1,8,"",0.001);
            Double_t *xpeaks = soup->GetPositionX();
            Double_t *ypeaks = soup->GetPositionY();
            
            for (int idx = 0; idx < nfound; idx++) {
                if (ypeaks[idx] > (e->GetMean() + e->GetRMS() + 130) ) {
                    std::cout << " i: " << idx << " x: " << xpeaks[idx] << " y: " << ypeaks[idx] << std::endl;
                    f->SetBinContent(xpeaks[idx],ypeaks[idx]);
                }
            }
            //--------------------------------------------------------------------
            
            h1->Draw();
            h1->SetMaximum(4096);
            
            d->SetLineColor(kRed);
            d->Draw("same");
            
            f->SetLineColor(kPink);
            f->Draw("same");
        }
    }
    
    TCanvas *fCanvas_y = canvas_wf->GetCanvas();
    fCanvas_y->cd();
    fCanvas_y->Update();
}
void MyMainFrame::DoText(const char */*text*/){
    // Handle text entry widgets.
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id = te->WidgetId();
    hslider->SetPosition(atoi(fTbh1->GetString()));
}
void MyMainFrame::DoSlider(Int_t pos){
    // Handle slider widgets.
    Int_t id;
    TGFrame *frm = (TGFrame *) gTQSender;
    if (frm->IsA()->InheritsFrom(TGSlider::Class())) {
        TGSlider *sl = (TGSlider*) frm;
        id = sl->WidgetId();
    } else {
        TGHSlider *sd = (TGHSlider *) frm;
        id = sd->WidgetId();
    }
    char buf[32];
    sprintf(buf, "%d", pos);
    // send to text box
    fTbh1->Clear();
    fTbh1->AddText(0, buf);
    // Re-align the cursor with the characters.
    //fTeh1->SetCursorPosition(fTeh1->GetCursorPosition());
    //fTeh1->Deselect();
    //gClient->NeedRedraw(fTeh1);

}
MyMainFrame::~MyMainFrame(){
    // clean up used widgets: frames, buttons, layout hints
    fMain->Cleanup();
    delete fMain;
    
    canvas_gain->Clear();
    delete canvas_gain;
    
    canvas_femb->Clear();
    delete canvas_femb;
    
    canvas_wf->Clear();
    delete canvas_wf;
}
void execute(){
    // popup the GUI
    new MyMainFrame(gClient->GetRoot(),400,600);
}

int main(int argc, char **argv){
    TApplication theApp("",&argc, argv);
    execute();
    theApp.Run();
    return 0;
}

// file browser
//pBrowser = new TGFileBrowser(p);
//p->SetEditable(kFALSE);
//pBrowser->AddFSDirectory("/", "/");
//pBrowser->GotoDir(gSystem->pwd());

/*tr_rawdata->Draw("wf:Iteration$ >> h1(19498,0,19498)",Form("(chan==%ld && subrun==%ld)",Nent_chan->GetIntNumber(),Nent_sub->GetIntNumber()),"L");
 
 TH1F *g1 = (TH1F*)gDirectory->Get("h1");
 Double_t *source = new Double_t[19498];
 TSpectrum *soup = new TSpectrum(100);
 for (Long_t idx = 0; idx < 19498; idx++) source[idx] = g1->GetBinContent(idx);
 soup->Background(source,wf.size(),12,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kTRUE, TSpectrum::kBackSmoothing9,kTRUE);
 for (Long_t idx = 0; idx < wf.size(); idx++) h2->SetBinContent(idx,source[idx]);
 Int_t nfound = soup->Search(h2,4,"",0.10);
 printf("Found %d candidate peaks, %f baseline\n",nfound, g1->GetMean());
 h2->Draw("same");
 
 g1->Draw("same");
 g1->SetMaximum(4096);
 
 h1->SetLineColor(kRed);
 h1->Draw("same");*/
