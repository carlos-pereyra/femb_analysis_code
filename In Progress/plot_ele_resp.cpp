#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TAxis.h"

//======================================================================
// useful for command line debugging access..
//----------------------------------------------------------------------
Double_t x[1000],y[5][1000];
TCanvas *c;

//======================================================================
double response(double *x, double *par){
    
    
    double f = 4.31054*exp(-2.94809*x[0]/par[1])*par[0]-2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*par[0]
    -2.6202*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*cos(2.38722*x[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*cos(5.18561*x[0]/par[1])*par[0]
    +0.762456*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0]
    -0.762456*exp(-2.82833*x[0]/par[1])*cos(2.38722*x[0]/par[1])*sin(1.19361*x[0]/par[1])*par[0]
    +0.762456*exp(-2.82833*x[0]/par[1])*cos(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0]
    -2.6202*exp(-2.82833*x[0]/par[1])*sin(1.19361*x[0]/par[1])*sin(2.38722*x[0]/par[1])*par[0]
    -0.327684*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0] +
    +0.327684*exp(-2.40318*x[0]/par[1])*cos(5.18561*x[0]/par[1])*sin(2.5928*x[0]/par[1])*par[0]
    -0.327684*exp(-2.40318*x[0]/par[1])*cos(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0]
    +0.464924*exp(-2.40318*x[0]/par[1])*sin(2.5928*x[0]/par[1])*sin(5.18561*x[0]/par[1])*par[0];
    
    if (x[0] >0&&x[0] < 10){
        return f;
    }else{
        return 0;
    }
}

//======================================================================
double compact_response(double *xv, double *par){
    
    double x= xv[0];
    
    if (x <=  0) return 0;
    if (x >= 10) return 0;
    
    x /= par[1];
    
    double phi_1 = x * 1.19361;
    double phi_2 = x * 2.59280;
    
    double f = 4.31054 * TMath::Exp(-0.119760*x) -5.240400*TMath::Cos(phi_1) + 1.524912*TMath::Sin(phi_1);
    
    f *= TMath::Exp(-0.425150*x);
    
    f += 0.929848*TMath::Cos(phi_2) -0.655368*TMath::Sin(phi_2);
    
    f *= TMath::Exp(-2.40318*x) * par[0];
    
    return f;
    
}


//======================================================================
void pretty_plot(TGraph *g, double line_width, double line_color,
                 const char *axis_label, const char *title_text){
    
    g -> SetLineWidth (line_width);
    g -> SetLineColor (line_color);
    g -> GetXaxis()->SetLimits(0.,10.);
    g -> GetXaxis()->SetTitle("Time (#mus)");
    g -> GetYaxis()->SetTitle(axis_label);
    g -> SetTitle(title_text);
    
    c = new TCanvas();
    c -> cd();
    
    g -> Draw("AL");

}



//======================================================================

void plot_ele_resp(){
    TF1 *f1 = new TF1("func1",response,0,10,2);
    TF1 *f2 = new TF1("func2",compact_response,0,10,2);
    
    //four gains 4.7, 7.8, 14, and 25 mV/fC
    
    f1->SetParameters(47.*1.012,0.5); // 4.7mV/fC for now with 0.5 us shaping time
    f2->SetParameters(47.*1.012,0.5); // 4.7mV/fC for now with 0.5 us shaping time
    
    for (Int_t i=0;i!=999;i++){
        x[i] = 10./1000.*(i+1);
        y[0][i] = f2 -> Eval(10./1000.*(i+1)) - f1->Eval(10./1000.*(i+1));
        y[1][i] = 100 * (f2 -> Eval(x[i]) - f1->Eval(x[i]))/f1 -> Eval(x[i]);
        y[2][i] = 100 * (f2 -> Eval(x[i]) - f1->Eval(x[i]))/f2 -> Eval(x[i]);
        y[3][i] = 200 *(f2 -> Eval(x[i]) - f1->Eval(x[i]))/ (f2 -> Eval(x[i])+f1->Eval(x[i]));
        y[4][i] = 1.0 / f1->Eval(x[i]);
        
    }
    
    x[999]=1000;
    
    for (int i = 0; i < 5; i++){
        y[i][999]=0;
    }
    
    
    TGraph *g1 = new TGraph(1000,x,y[0]);
    TGraph *g2 = new TGraph(1000,x,y[1]);
    TGraph *g3 = new TGraph(1000,x,y[2]);
    TGraph *g4 = new TGraph(1000,x,y[3]);
    TGraph *g5 = new TGraph(1000,x,y[4]);
    
    pretty_plot (g5, 2.5, 6, "1.0 / f_{1},   [fC/mV]","Inverse Response Function");
    
    delete c;
    
    c = new TCanvas();
    c -> cd();
    
    g5 -> GetYaxis() -> SetRangeUser(-100,100);
    g5 -> Draw("AL");
    
    pretty_plot (g4, 2.5, 6, "2 #times (f_{2} - f_{1}) / (f_{1} + f_{2})  [%]","Perc. Diff, Avg Func Denominator");
    pretty_plot (g3, 2.5, 4, "(f_{2} - f_{1}) / f_{2}  [%]","Perc. Diff, New Func Denominator");
    pretty_plot (g2, 2.5, 2, "(f_{2} - f_{1}) / f_{1}  [%]","Perc. Diff, Old Func Denominator");
    pretty_plot (g1, 2.5, 1, "#Delta Amplitude  [mV/fC]","Absoulte Difference in Evaluation");
    
}
