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
    
    // frame            //
    vf_1 = new TGVerticalFrame(hfm_1, 10, 10);
    hf_1 = new TGHorizontalFrame(vf_1, 10, 10); //internal to left vertical
    l_frame = new TGCompositeFrame(vf_1, 10, 10, kSunkenFrame);
    
    vf_1->AddFrame(hf_1,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,1,2));
    vf_1->AddFrame(l_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandY | kLHintsExpandX,0,0,1,2));

    hfm_1->AddFrame(vf_1, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    
    // create canvas //* try moving this into DoDrawFEMB() -Cp
    canvas_femb = new TRootEmbeddedCanvas("canvas_femb", l_frame, w, 0.5*h);     // create canvas widget
    l_frame->AddFrame(canvas_femb, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1) );
    
    // draw button
    TGTextButton *draw = new TGTextButton(hf_1, "&Draw Summary");
    draw->Connect("Clicked()","MyMainFrame",this,"DoDrawFEMB()");
    hf_1->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    vf_1->AddFrame(hf_1, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    // right //
    vf_2 = new TGVerticalFrame(hfm_1, 10, 10);
    hf_2 = new TGHorizontalFrame(vf_2, 10, 10); //internal to right vertical
    r_frame = new TGCompositeFrame(vf_2, 10, 10, kSunkenFrame);

    vf_2->AddFrame(hf_2,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,1,2));
    vf_2->AddFrame(r_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandY | kLHintsExpandX,0,0,1,2));

    hfm_1->AddFrame(vf_2, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    
    // create canvas
    canvas_gain = new TRootEmbeddedCanvas("canvas_gain", r_frame, w, 0.5*h);     // create canvas widget
    r_frame->AddFrame(canvas_gain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1) );
    
    // selection
    fCombo = new TGComboBox(hf_2, kComboID1);
    fCombo->Resize(w/2, 20);
    fCombo->AddEntry("Gain", kComboID1);
    fCombo->AddEntry("Residuals", kComboID2);
    //fCombo->Connect("Selected(Int_t)", "MyMainFrame", this, "Echo(Int_t)");
    
    // slider control
    hslider = new TGHSlider(hf_2,150,kSlider1 | kScaleDownRight,HSId1);
    hslider->SetRange(0,127);
    hslider->SetPosition(0);
    hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDraw(Int_t)");
    
    hf_2->AddFrame(hslider, new TGLayoutHints(kLHintsRight, 5, 5, 3, 4));
    hf_2->AddFrame(fCombo,  new TGLayoutHints(kLHintsLeft, 5, 5, 3, 4));

    vf_2->AddFrame(hf_2, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    //--------------------------------------//
    //          horizontal level 2          //
    //--------------------------------------//
    hfm_2 = new TGHorizontalFrame(fMain, w, h);
    
    // create canvas
    canvas_wf = new TRootEmbeddedCanvas("canvas_wf", hfm_2, w*2, 0.5*h, kSunkenFrame);     // create canvas widget
    hfm_2->AddFrame(canvas_wf, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX));
    
    //--------------------------------------//
    //          horizontal level 3          //
    //--------------------------------------//
    hfm_3 = new TGHorizontalFrame(fMain,w,h);
    
    // number entry
    Nent_chan = new TGNumberEntry(hfm_3, 1, 6, kENTRY1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,127);
    Nent_chan->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawFit()"); //??
    
    Nent_sub = new TGNumberEntry(hfm_3, 1, 6, kENTRY2, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,10);
    Nent_sub->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawFit()");
    
    Nent_peak = new TGNumberEntry(hfm_3, 1, 6, kENTRY3, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,10);
    Nent_peak->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawFit()");
    
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
    
    TGTextButton *wf = new TGTextButton(hfm_4, "&Draw Waveform");
    wf->Connect("Clicked()","MyMainFrame",this,"DoDrawWf()");
    hfm_4->AddFrame(wf, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 3, 4));
    
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

void MyMainFrame::Echo(Int_t iden){
    Int_t pos = hslider->GetPosition();
    switch (iden) {
        case kComboID1:
            //hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDrawGain(Int_t)");
            DoDrawGain(pos);
            std::cout << " kCom1 - id: " << iden << std::endl;
            break;
        case kComboID2:
            //hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDrawRes(Int_t)");
            DoDrawRes(pos);
            std::cout << " kCom2 - id: " << iden << std::endl;
            break;
        //default:
        //    break;
    }
}

void MyMainFrame::DoDraw(Int_t pos){
    Int_t food = fCombo->WidgetId();
    
    switch (fCombo->GetSelected()) {
        case kComboID1:
            //hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDrawGain(Int_t)");
            DoDrawGain(pos);
            //std::cout << " kCom1 - id: " << fCombo->GetSelected() << std::endl;
            break;
            
        case kComboID2:
            //hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDrawRes(Int_t)");
            DoDrawRes(pos);
            //std::cout << " kCom2 - id: " << fCombo->GetSelected() << std::endl;
            break;
            
        default:
            break;
    }
}

void MyMainFrame::DoDrawGain(Int_t pos){
    Int_t id;
    TGFrame *frm = (TGFrame *) gTQSender;
    if (frm->IsA()->InheritsFrom(TGSlider::Class())) {
        TGSlider *sl = (TGSlider*) frm;
        id = sl->WidgetId();
    } else {
        TGHSlider *sd = (TGHSlider *) frm;
        id = sd->WidgetId();
    }
    
    Double_t RT_f, GN_f, GN_m, shape_v;
    Int_t charge, s, c, p, gain, shape, gain_v;

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
    tr_rawdata->SetBranchAddress("gain_v", &gain_v);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("shape_v", &shape_v);  // initialize subrun branch
    
    Long64_t nEntries(tr_rawdata->GetEntries());
    tr_rawdata->GetEntry(0);
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        if (c == pos) { gain = gain; }
    }
    
    hprof_f = new TProfile("hprof_f",Form("Gain: %d (mV/fC), Shape: %1.1f (#mu s), Channel: %d", gain_v, shape_v, pos),10,0,218200);
    if(gain == 2) hprof_f->SetMaximum(25);
    if(gain == 3) hprof_f->SetMaximum(40);
    tr_rawdata->Draw("(0.33*6241)*(GN_f/charge):charge>>hprof_f",Form("c==%d",pos),"*");
    hprof_f->SetLineColor(1);
    hprof_f->SetXTitle("Injected Charge (e-)");
    hprof_f->SetYTitle("Gain (mV/fC)");
    hprof_f->SetStats(0);

    hprof_m = new TProfile("hprof_m","",10,0,218200);
    tr_rawdata->Draw("(0.33*6241)*(GN_m/charge):charge>>hprof_m",Form("c==%d",pos),"same");
    hprof_m->SetLineColor(4);
    
    hprof_gs = new TProfile("hprof_gs","",10,0,218200);
    tr_rawdata->Draw("gain_v:charge>>hprof_gs",Form("c==%d && s>1",pos),"same");
    hprof_gs->SetLineColor(3);
    
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry("hprof_f"," Response Func. Fit: results ","l");
    legend->AddEntry("hprof_m"," Peak Measurement: results ","lep");
    legend->Draw();
    
    TCanvas *fCanvas = canvas_gain->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
    fCanvas->Clear();
}


void MyMainFrame::DoDrawRes(Int_t pos){
    
    Double_t RT_f, GN_f, GN_m, shape_v, fpga_voltage;
    Int_t charge, s, c, p, gain, shape, gain_v;
    
    TTree *tr_rawdata;
    TFile *inputFile = new TFile(tbuf1->GetString(), "READ");
    
    tr_rawdata = (TTree*) inputFile->Get("roast_beef");
    if( !tr_rawdata )
    {   std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0); }
    tr_rawdata->SetBranchAddress("RT_f", &RT_f);        // units [bins]
    tr_rawdata->SetBranchAddress("GN_f", &GN_f);        // units [ADC]
    tr_rawdata->SetBranchAddress("GN_m", &GN_m);        // units [ADC]
    tr_rawdata->SetBranchAddress("charge", &charge);    // units [e-]
    tr_rawdata->SetBranchAddress("fpga_voltage", &fpga_voltage);    // units [V]
    tr_rawdata->SetBranchAddress("s", &s);              //
    tr_rawdata->SetBranchAddress("c", &c);              // initialize subrun branch
    tr_rawdata->SetBranchAddress("p", &p);              // initialize subrun branch
    tr_rawdata->SetBranchAddress("gain", &gain);        // initialize subrun branch
    tr_rawdata->SetBranchAddress("shape", &shape);      // initialize subrun branch
    tr_rawdata->SetBranchAddress("gain_v", &gain_v);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("shape_v", &shape_v);  // initialize subrun branch
    
    Long64_t nEntries(tr_rawdata->GetEntries());
    tr_rawdata->GetEntry(0);
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        tr_rawdata->GetEntry(entry);
        if (c == pos) { gain = gain; }
    }
    residual_f = new TProfile("residual_f",Form("Gain: %d (mV/fC), Shape: %1.1f (#mu s), Channel: %d", gain_v, shape_v, pos),10,0,218200);
    tr_rawdata->Draw("(((6241*GN_f/3)/gain_v)-charge):charge>>residual_f",Form("c==%d",pos),"*");
    residual_f->SetLineColor(1);
    residual_f->SetXTitle("Injected Charge (e-)");
    residual_f->SetYTitle("Residual (e-)");
    residual_f->SetStats(0);
    
    residual_m = new TProfile("residual_m","",10,0,218200);
    tr_rawdata->Draw("(((6241*GN_m/3)/gain_v)-charge):charge>>residual_m",Form("c==%d",pos),"*");
    residual_m->SetLineColor(4);

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry("residual_f","residual fit: results","l");
    legend->AddEntry("residual_m","residual meas.: results","lep");
    legend->Draw();
    
    TCanvas *fCanvas = canvas_gain->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
    fCanvas->Clear();
}


void MyMainFrame::DoDrawFEMB(){
    Double_t RT_f, GN_f, GN_m;
    Int_t charge, s, c, p;
    
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
    
    h2 = new TH2F("h2",Form("FEMB Summary %s", "20170808T143522"),128,0,128,10,0,218200);
    tr_rawdata->Draw("charge:c>>h2","(0.33*6241)*GN_m/charge*(p==0 && GN_f<1e6)","colz");
    h2->SetStats(0);
    h2->SetXTitle("Channel");
    h2->SetYTitle("Injected Charge (e-)");
    
    TCanvas *fCanvas_x = canvas_femb->GetCanvas();
    fCanvas_x->cd();
    fCanvas_x->Update();
}


void MyMainFrame::DoDrawFit(){
    TGNumberEntry *sender = (TGNumberEntry *) gTQSender;
    Int_t id = sender->WidgetId();
    
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id_t = te->WidgetId();
    
    wf.clear();
    const int const_numSubrun = 64, const_numChan = 128;
    unsigned short subrunIn, chanIn;    //output tree and variable
    
    TTree *tr_rawdata;
    TFile *inputFile = new TFile(tbuf2->GetString(), "READ");
    
    // get branches
    TTree *rawdata; Int_t shape;
    TFile *inputF = new TFile(tbuf1->GetString(), "READ");
    rawdata = (TTree*) inputF->Get("roast_beef");
    if( !rawdata )
    {   std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0); }
    rawdata->SetBranchAddress("shape", &shape);      // initialize subrun branch
    rawdata->GetEntry(0);
    
    // get wf branches
    tr_rawdata = (TTree*) inputFile->Get("femb_wfdata");
    if( !tr_rawdata )
    {   std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0); }
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
            foo.set_peaks(subrunIn,chanIn,f_p,shape,wf);
            foo.set_baseline(subrunIn,chanIn,f_p,wf);
            bl = foo.get_baseline(1,chanIn,f_p);
        }
        
        if( chanIn == c && subrunIn == s){
            for( unsigned int s = 0 ; s < wfIn->size() ; s++) wf.push_back( wfIn->at(s) );
            foo.set_run_info(subrunIn,chanIn,f_p);          // initialize vector of structures size
            foo.set_run_channel(subrunIn,chanIn);
            foo.set_peaks(subrunIn,chanIn,f_p,shape,wf);
            
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
void MyMainFrame::DoDrawWf(){
    wf.clear();

    TGNumberEntry *sender = (TGNumberEntry *) gTQSender;
    Int_t id = sender->WidgetId();
    
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id_t = te->WidgetId();
    
    const int const_numSubrun = 64, const_numChan = 128;
    unsigned short subrunIn, chanIn;    //output tree and variable
    
    // get branches
    TTree *rawdata; Int_t shape;
    TFile *inputF = new TFile(tbuf1->GetString(), "READ");
    rawdata = (TTree*) inputF->Get("roast_beef");
    if( !rawdata )
    {   std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0); }
    rawdata->SetBranchAddress("shape", &shape);      // initialize subrun branch
    rawdata->GetEntry(0);

    // get wf branches
    TTree *tr_rawdata;
    TFile *inputFile = new TFile(tbuf2->GetString(), "READ");
    tr_rawdata = (TTree*) inputFile->Get("femb_wfdata");
    if( !tr_rawdata )
    {   std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0); }
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
            for( unsigned int s = 0 ; s < wfIn->size() ; s++ ) wf.push_back(wfIn->at(s));
            TH1F *h1 = new TH1F("h1",Form("Waveform Plot Channel: %d Subrun: %d", chanIn, subrunIn),wfIn->size(),1,wfIn->size());
            
            // fill hist
            for (unsigned long idx = 0; idx<wf.size(); idx++) h1->SetBinContent(idx, wf.at(idx));
            TH1F *d = new TH1F("d","",wf.size(),0,wf.size());
            TH1F *e = new TH1F("e","",wf.size(),0,wf.size());
            TH1F *f = new TH1F("f","",wf.size(),0,wf.size());
            Double_t *source = new Double_t[wf.size()];
            TSpectrum *soup = new TSpectrum(80);
            
            // fill hist
            for (unsigned short idx = 0; idx < wf.size(); idx++) source[idx]=wf.at(idx);
            soup->Background(source,wf.size(),29,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kTRUE, TSpectrum::kBackSmoothing9,kTRUE);
            for (unsigned short idx = 0; idx < wf.size(); idx++) d->SetBinContent(idx,source[idx]);
            for (unsigned short idx = 0; idx < wf.size(); idx++) e->Fill(d->GetBinContent(idx));
            
            // search for peaks
            Int_t sigma = 0;
            if (shape == 0) sigma = 4;
            if (shape == 1) sigma = 4;
            if (shape == 2) sigma = 8;
            if (shape == 3) sigma = 8;
            Int_t nfound = soup->Search(h1,sigma,"",0.01);
            Double_t *xpeaks = soup->GetPositionX();
            Double_t *ypeaks = soup->GetPositionY();
            
            for (int idx = 0; idx < nfound; idx++) {
                if (ypeaks[idx] > (e->GetMean() + e->GetRMS() + 130) ) {
                    std::cout << " i: " << idx << " x: " << xpeaks[idx] << " y: " << ypeaks[idx] << std::endl;
                    //f->SetBinContent(xpeaks[idx]+1,ypeaks[idx]);
                }
            }
            //--------------------------------------------------------------------
            
            h1->Draw();
            h1->SetMaximum(4096);
            //h1->SetMaximum(4096);
            
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
