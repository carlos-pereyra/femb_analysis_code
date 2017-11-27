#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TGClient.h>

#include "TCanvas.h"

#include "TGFrame.h"
#include "TGSlider.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGLayout.h"
#include "TGTextBuffer.h"
#include "TGSplitter.h"
#include "TGFileBrowser.h"
#include "TGNumberEntry.h"

#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TSystem.h"
#include "TH2.h"
#include "TObject.h"
#include "TROOT.h"

#include <TSpectrum.h>
#include <TROOT.h>
#include <TClass.h>
#include <TApplication.h>
#include <TVirtualX.h>
#include <TVirtualPadEditor.h>
#include <TGResourcePool.h>
#include <TGListBox.h>
#include <TGListTree.h>
#include <TGFSContainer.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGMsgBox.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TGComboBox.h>
#include <TGTab.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>
#include <TGFileDialog.h>
#include <TGTextEdit.h>
#include <TGShutter.h>
#include <TGProgressBar.h>
#include <TGColorSelect.h>
#include <RQ_OBJECT.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TEnv.h>
#include <TFile.h>
#include <TKey.h>
#include <TGDockableFrame.h>
#include <TGFontDialog.h>

class TGWindow;
class TGMainFrame;
class TRootEmbeddedCanvas;

class MyMainFrame {
    RQ_OBJECT("MyMainFrame")

private:
    TGMainFrame         *fMain;
    TGHorizontalFrame   *hfm_1, *hfm_2, *hfm_3, *hf_1, *hf_2;
    TGVerticalFrame     *vf_1, *vf_2;
    TGCompositeFrame    *l_frame, *r_frame, *top, *mid, *bot;
    TRootEmbeddedCanvas *canvas_gain, *canvas_femb, *canvas_wf;
    
    TGTextEntry         *fTeh1;
    TGTextBuffer        *fTbh1;
    TGLayoutHints       *fBly;
    TGFileBrowser       *pBrowser;
    TGNumberEntry       *Nent_chan, *Nent_sub;
    TGHSlider           *hslider;
    
    TH2F                *h2;
    TProfile            *hprof;
    
public:
    MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
    virtual ~MyMainFrame();
    void DoDrawGain(Int_t pos);
    void DoDrawFEMB();
    void DoDrawWf(Long_t);
    void DoDrawFit(Long_t);
    void DoText(const char *text);
    void DoSlider(Int_t pos);
};

enum ETestCommandIdentifiers {
    HId1,
    HSId1,
    kENTRY1,
    kENTRY2,
    kNESInteger,
    kNELLimitMinMax,
    
    /*kBackOrder2 =0,
    kBackOrder4 =1,
    kBackOrder6 =2,
    kBackOrder8 =3,
    kBackIncreasingWindow =0,
    kBackDecreasingWindow =1,
    kBackSmoothing3 =3,
    kBackSmoothing5 =5,
    kBackSmoothing7 =7,
    kBackSmoothing9 =9,
    kBackSmoothing11 =11,
    kBackSmoothing13 =13,
    kBackSmoothing15 =15*/
};

