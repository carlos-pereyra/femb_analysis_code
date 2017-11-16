#include <TApplication.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include "dec_def.h"

Double_t myfunction(Double_t *x, Double_t *par)
{
    Float_t xx = x[0];
    Double_t f = TMath::Abs(par[0]*xx+par[1]);
    return f;
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h){
    // create a main frame
    fMain = new TGMainFrame(p,w,h);
    
    // create a horizontal frame widget with buttons
    hf_main = new TGHorizontalFrame(fMain, w, h);
    
    vf_1 = new TGVerticalFrame(hf_main, 10, 10);
    hf_1 = new TGHorizontalFrame(vf_1, 10, 10); //internal to left
    vf_1->AddFrame(hf_1,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,1,2));
    l_frame = new TGCompositeFrame(vf_1, 10, 10, kSunkenFrame);
    vf_1->AddFrame(l_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandY | kLHintsExpandX,0,0,1,2));

    vf_2 = new TGVerticalFrame(hf_main, 10, 10);
    hf_2 = new TGHorizontalFrame(vf_2, 10, 10); //internal to right
    vf_2->AddFrame(hf_2,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,1,2));
    r_frame = new TGCompositeFrame(vf_2, 10, 10, kSunkenFrame);
    vf_2->AddFrame(r_frame, new TGLayoutHints(kLHintsTop | kLHintsExpandY | kLHintsExpandX,0,0,1,2));

    hf_main->AddFrame(vf_1, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    hf_main->AddFrame(vf_2, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    
    // left frame
    //
    canvas_femb = new TRootEmbeddedCanvas("canvas_femb", l_frame, w, 0.5*h);     // create canvas widget
    l_frame->AddFrame(canvas_femb, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1) );
    
    // draw button
    TGTextButton *draw = new TGTextButton(hf_1, "&Draw Summary","DoDrawFEMB()");
    //draw->Connect("Clicked()","MyMainFrame",this,"DoDrawFEMB()");
    hf_1->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    // right frame
    //
    
    // create canvas
    canvas_gain = new TRootEmbeddedCanvas("canvas_gain", r_frame, w, 0.5*h);     // create canvas widget
    r_frame->AddFrame(canvas_gain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1) );
    
    // text box
    fTeh1 = new TGTextEntry(hf_2, fTbh1 = new TGTextBuffer(10), HId1);
    fTbh1->AddText(0, "0");
    fBly = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 3, 0);
    fTeh1->Connect("TextChanged(char*)", "MyMainFrame", this, "DoText(char*)");
    hf_2->AddFrame(fTeh1, fBly);
    
    // slider control
    TGHSlider *hslider = new TGHSlider(hf_2,150,kSlider1 | kScaleDownRight,HSId1);
    hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoSlider(Int_t)");
    hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDrawGain(Int_t)");
    hslider->SetRange(0,127);
    hslider->SetPosition(0);
    hf_2->AddFrame(hslider, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    // exit button
    TGTextButton *exit = new TGTextButton(hf_2, "&Exit", "gApplication->Terminate(0)");
    hf_2->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    // add children to parent frame
    fMain->AddFrame(hf_main, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

    // set a name to the main frame
    fMain->SetWindowName("Scroll and roll");
    
    // map all subwindows of main frame
    fMain->MapSubwindows();
    
    // initialize the layout algorithm
    fMain->Resize(fMain->GetDefaultSize());
    
    // map main frame
    fMain->MapWindow();
}
MyMainFrame::~MyMainFrame(){
    canvas_gain->Clear();
    delete canvas_gain;
    // clean up used widgets: frames, buttons, layout hints
    fMain->Cleanup();
    delete fMain;
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
    tr_rawdata->SetBranchAddress("RT_f", &RT_f);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("GN_f", &GN_f);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("GN_m", &GN_m);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("charge", &charge);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("s", &s);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("c", &c);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("p", &p);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("gain", &gain);    // initialize subrun branch
    tr_rawdata->SetBranchAddress("shape", &shape);    // initialize subrun branch
    
    hprof = new TProfile("hprof",Form("2017929 g3 s2 channel %d", pos),10,0,218200);
    hprof->SetMaximum(30);
    if(gain == 2) hprof->SetMaximum(20);
    if(gain == 3) hprof->SetMaximum(30);
    tr_rawdata->Draw("(0.33*6241)*(GN_m/charge):charge>>hprof",Form("c==%d", pos),"*");
    
    TCanvas *fCanvas = canvas_gain->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
}
void DoDrawFEMB(){
    Double_t RT_f1, GN_f1, GN_m1;
    Int_t charge1, s1, c1, p1;
    
    TTree *tr_rawdata1;
    std::string file_name1 = "Data/20170808T143522/g2_s2_extpulse/Results.root"; // F_P => Reads only one file
    TFile *inputFile1 = new TFile(file_name1.c_str(), "READ");
    tr_rawdata1 = (TTree*) inputFile1->Get("roast_beef");
    if( !tr_rawdata1 ){
        std::cout << "Error opening input file tree, exiting" << std::endl;
        gSystem->Exit(0);
    }
    tr_rawdata1->SetBranchAddress("RT_f", &RT_f1);    // initialize subrun branch
    tr_rawdata1->SetBranchAddress("GN_f", &GN_f1);    // initialize subrun branch
    tr_rawdata1->SetBranchAddress("GN_m", &GN_m1);    // initialize subrun branch
    tr_rawdata1->SetBranchAddress("charge", &charge1);    // initialize subrun branch
    tr_rawdata1->SetBranchAddress("s", &s1);    // initialize subrun branch
    tr_rawdata1->SetBranchAddress("c", &c1);    // initialize subrun branch
    tr_rawdata1->SetBranchAddress("p", &p1);    // initialize subrun branch
    
    //h2 = new TH2F("h2","a hist",128,0,128,10,0,218200);
    //tr_rawdata1->Draw("charge1:c1>>h2","(0.33*6241)*GN_f/charge1*(p1==20 && GN_f1 < 1e9)","colz");
    
    //TCanvas *fCanvas = canvas_femb->GetCanvas();
    //fCanvas->cd();
    //fCanvas->Update();
}
void MyMainFrame::DoText(const char */*text*/){
    // Handle text entry widgets.
    TGTextEntry *te = (TGTextEntry *) gTQSender;
    Int_t id = te->WidgetId();
    //hslider->SetPosition(atoi(fTbh1->GetString()));
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
    fTeh1->SetCursorPosition(fTeh1->GetCursorPosition());
    fTeh1->Deselect();
    gClient->NeedRedraw(fTeh1);

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
