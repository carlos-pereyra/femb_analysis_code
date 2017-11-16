#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TGClient.h>

#include "TCanvas.h"

#include "TGSlider.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLayout.h"
#include "TGTextBuffer.h"
#include "TGSplitter.h"

#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"
#include "TH2.h"

class TGWindow;
class TGMainFrame;
class TRootEmbeddedCanvas;

class MyMainFrame {
    RQ_OBJECT("MyMainFrame")
private:
    TGMainFrame         *fMain;
    TGHorizontalFrame   *hf_main, *hf_1, *hf_2;
    TGVerticalFrame     *vf_1, *vf_2;
    TGCompositeFrame    *l_frame, *r_frame;
    TRootEmbeddedCanvas *canvas_gain;
    TRootEmbeddedCanvas *canvas_femb;
    
    TGTextEntry         *fTeh1;
    TGTextBuffer        *fTbh1;
    TGLayoutHints       *fBly;
    
    TH2F                *h2;
    TProfile            *hprof;
    
public:
    MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
    virtual ~MyMainFrame();
    void DoDrawGain(Int_t pos);
    void DoDrawFEMB();
    void DoText(const char *text);
    void DoSlider(Int_t pos);
};

enum ETestCommandIdentifiers {
    HId1,
    HSId1
};
