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
    
    // create canvas widget
    fEcanvas = new TRootEmbeddedCanvas("Ecanvas", fMain,w,h);
    fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 10, 1) );
    
    // create a horizontal frame widget with buttons
    hframe = new TGHorizontalFrame(fMain, w, 40);
    
    // draw button
    //TGTextButton *draw = new TGTextButton(hframe,"&Draw");
    //draw->Connect("Clicked()","MyMainFrame", this, "DoDraw()");
    //hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4) );
    
    // text box
    fTeh1 = new TGTextEntry(hframe, fTbh1 = new TGTextBuffer(10), HId1);
    fTbh1->AddText(0, "0");
    fBly = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 3, 0);
    fTeh1->Connect("TextChanged(char*)", "MyMainFrame", this, "DoText(char*)");
    hframe->AddFrame(fTeh1, fBly);
    
    // slider control
    TGHSlider *hslider = new TGHSlider(hframe,150,kSlider1 | kScaleDownRight,HSId1);
    hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoSlider(Int_t)");
    hslider->Connect("PositionChanged(Int_t)","MyMainFrame",this,"DoDraw(Int_t)");
    hslider->SetRange(0,127);
    hslider->SetPosition(0);
    hframe->AddFrame(hslider, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    // exit button
    TGTextButton *exit = new TGTextButton(hframe, "&Exit", "gApplication->Terminate(0)");
    hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
    
    fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
    
    // set a name to the main frame
    fMain->SetWindowName("Simple Example");
    
    // map all subwindows of main frame
    fMain->MapSubwindows();
    
    // initialize the layout algorithm
    fMain->Resize(fMain->GetDefaultSize());
    
    // map main frame
    fMain->MapWindow();
}
MyMainFrame::~MyMainFrame(){
    //fEcanvas->Clear();
    //delete fEcanvas;
    // clean up used widgets: frames, buttons, layout hints
    fMain->Cleanup();
    delete fMain;
}
void MyMainFrame::DoDraw(Int_t pos){
    Double_t RT_f, GN_f, GN_m;
    Int_t charge, s, c, p;
    
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
    
    hprof = new TProfile("hprof",Form("2017929 g3 s2 channel %d", pos),10,0,218200);
    hprof->SetMaximum(30);
    //if(gain = 2) hprof->SetMaximum(20);
    //if(gain = 3) hprof->SetMaximum(30);
    tr_rawdata->Draw("(0.33*6241)*(GN_m/charge):charge>>hprof",Form("c==%d", pos),"*");
    
    TCanvas *fCanvas = fEcanvas->GetCanvas();
    fCanvas->cd();
    fCanvas->Update();
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
    new MyMainFrame(gClient->GetRoot(),600,400);
}

int main(int argc, char **argv){
    TApplication theApp("",&argc, argv);
    execute();
    theApp.Run();
    return 0;
}
