#include <TQObject.h>
#include <RQ_OBJECT.h>

#include "TCanvas.h"

#include "TGSlider.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLayout.h"
#include "TGTextBuffer.h"

#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"

class TGWindow;
class TGMainFrame;
class TRootEmbeddedCanvas;

class MyMainFrame {
    RQ_OBJECT("MyMainFrame")
private:
    TGMainFrame         *fMain;
    TGHorizontalFrame   *hframe;
    TRootEmbeddedCanvas *fEcanvas;
    
    TGTextEntry         *fTeh1;
    TGTextBuffer        *fTbh1;
    TGLayoutHints       *fBly;
    
    TProfile            *hprof;
    
public:
    MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
    virtual ~MyMainFrame();
    void DoDraw(Int_t pos);
    void DoText(const char *text);
    void DoSlider(Int_t pos);
};

enum ETestCommandIdentifiers {
    HId1,
    HSId1
};
