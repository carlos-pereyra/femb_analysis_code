#include <TApplication.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include "TSpectrum.h"
#include "TROOT.h"
#include "dec_def.h"

std::vector<unsigned short> *wfIn;
std::vector<unsigned short> wf;

Double_t myfunction(Double_t *x, Double_t *par)
{
    Float_t xx = x[0];
    Double_t f = TMath::Abs(par[0]*xx+par[1]);
    return f;
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h){
    // create a main frame
    fMain = new TGMainFrame(p,w,h);
    
    //--------------------------------------//
    //          horizontal level 1          //
    //--------------------------------------//
    hfm_1 = new TGHorizontalFrame(fMain, w, h);
    top = new TGCompositeFrame(hfm_1,10,10, kSunkenFrame);
    
    //          left side: frame            //
    vf_1 = new TGVerticalFrame(hfm_1, 10, 10);
    hf_1 = new TGHorizontalFrame(vf_1, 10, 10); //internal to left vertical
    l_frame = new TGCompositeFrame(vf_1, 10, 10, kSunkenFrame);
    
    vf_1->AddFrame(hf_1,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,1,2));
    vf_1->AddFrame(l_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandY | kLHintsExpandX,0,0,1,2));

    hfm_1->AddFrame(vf_1, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    
    // left side: create canvas
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
    
    // right side: text box
    //fTeh1 = new TGTextEntry(hf_2, fTbh1 = new TGTextBuffer(10), HId1);
    //fTbh1->AddText(0, "0");
    //fBly = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 3, 0);
    //fTeh1->Connect("TextChanged(char*)", "MyMainFrame", this, "DoText(char*)");
    //hf_2->AddFrame(fTeh1, fBly);
    
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
    hfm_3 = new TGHorizontalFrame(fMain, w, h);
    //bot = new TGCompositeFrame(hfm_3, 800, 10, kSunkenFrame);

    // exit button
    TGTextButton *exit = new TGTextButton(hfm_3, "&Exit", "gApplication->Terminate(0)");
    
    // number entry
    Nent_chan = new TGNumberEntry(hfm_3, 1, 6, kENTRY1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,127);
    Nent_chan->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawWf(Long_t)");

    Nent_sub = new TGNumberEntry(hfm_3, 1, 6, kENTRY2, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0,10);
    Nent_sub->Connect("ValueSet(Long_t)", "MyMainFrame", this, "DoDrawWf(Long_t)");
    
    std::cout << "stuff: " << kENTRY1 << std::endl;
    // add all widgets to frame
    hfm_3->AddFrame(exit, new TGLayoutHints(kLHintsLeft, 5, 5, 3, 4));
    hfm_3->AddFrame(Nent_chan, new TGLayoutHints(kLHintsCenterX, 0, 0, 0, 0));
    hfm_3->AddFrame(Nent_sub, new TGLayoutHints(kLHintsRight, 0, 0, 0, 0));
    
    //          fmain stuff         //
    
    // add children to parent frame
    fMain->AddFrame(hfm_1, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    fMain->AddFrame(hfm_2, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    fMain->AddFrame(hfm_3, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    
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
    
    TTree *tr_rawdata;
    std::string file_name = "Data/20170808T143522/g2_s2_extpulse/Results.root"; // F_P => Reads only one file
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
    if(gain == 2) hprof->SetMaximum(40);
    if(gain == 3) hprof->SetMaximum(40);
    tr_rawdata->Draw("(0.33*6241)*(GN_m/charge):charge>>hprof",Form("c==%d", pos),"*");
    
    TCanvas *fCanvas = canvas_gain->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
}
void MyMainFrame::DoDrawFEMB(){
    Double_t RT_f, GN_f, GN_m;
    Int_t charge, s, c, p;
    Char_t *pesto_pasta;
    
    TTree *tr_rawdata;
    std::string file_name = "Data/20170808T143522/g2_s2_extpulse/Results.root"; // F_P => Reads only one file
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
    
    h2 = new TH2F("h2",Form("FEMB Summary: %s", pesto_pasta),128,0,128,10,0,218200);
    tr_rawdata->Draw("charge:c>>h2","(0.33*6241)*GN_m/charge*(p==11)","colz");
    h2->SetStats(0);
    
    TCanvas *fCanvas_x = canvas_femb->GetCanvas();
    fCanvas_x->cd();
    fCanvas_x->Update();
}
void MyMainFrame::DoDrawWf(Long_t pos){
    TGNumberEntry *sender = (TGNumberEntry *) gTQSender;
    Int_t id = sender->WidgetId();
    
    wf.clear();
    const int const_numSubrun = 64, const_numChan = 128;
    unsigned short subrunIn, chanIn;    //output tree and variable
    
    TString food = Form("Data/20170808T143522/g2_s2_extpulse/gainMeasurement_femb_1-parseBinaryFile.root");
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
    
    //tr_rawdata->Draw("wf;1",Form("(chan==%d && subrun==%d)",Nent_chan->GetIntNumber(),Nent_sub->GetIntNumber()),"L");
    
    Long64_t nEntries(tr_rawdata->GetEntries());          // 11 subrun * 128 channels = 1408 entries (wf)
    tr_rawdata->GetEntry(0);
    
    for(Long64_t entry(0); entry < nEntries; ++entry) {
        tr_rawdata->GetEntry(entry);
        if( subrunIn < 0 || subrunIn >= const_numSubrun ) continue;
        if( chanIn < 0 || chanIn >= const_numChan ) continue;
        if( chanIn == Nent_chan->GetIntNumber() && subrunIn == Nent_sub->GetIntNumber()){
            for( unsigned int s = 0 ; s < wfIn->size() ; s++ ) wf.push_back(wfIn->at(s));
            TH1F *h1 = new TH1F("h1",Form("Waveform Plot Channel: %d Subrun: %d", chanIn, subrunIn),wf.size(),1,wf.size());
            for (int idx = 0; idx<wf.size(); idx++) h1->SetBinContent(idx, wf.at(idx));
            h1->Draw();
            h1->SetMaximum(4096);
            
            TH1F *d = new TH1F("d","",wf.size(),0,wf.size());
            Double_t *source = new Double_t[wf.size()];
            TSpectrum *s = new TSpectrum();
            /*for (int idx = 0; idx < wf.size(); idx++) source[idx]=wf.at(idx);
            s->Background(source,wf.size(),6,kBackIncreasingWindow,kBackOrder2,kFALSE, kBackSmoothing3,kFALSE);
            for (int idx = 0; idx < wf.size(); idx++) d->SetBinContent(idx,source[idx]);
            d->SetLineColor(kRed);
            d->Draw("same L");*/
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
