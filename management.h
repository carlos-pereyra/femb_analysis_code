
#include<cstring>

//vector<double> getParameters();
//vector<string> getFiles();

vector<double> getParameters()
{/* -----------------------------------------------------------------------------
  |  Function: getQuestions()                                                   |
  |  Bio: Gets and stores questions from a text file.                           |
  ------------------------------------------------------------------------------*/
    
    /*------------------------------
     |       Opening File...        |
     |                              |
     ------------------------------*/
    ifstream dataIn("Parameters.txt");
    if (!dataIn.is_open()) {
        cout << "\'Cannot Read Control File\' - Error in: getParameters()" << endl;
        cout << "Must Name Question File: \'Parameters.txt\' " << endl;
        abort();
    }
    
    /*------------------------------
     |    Place Questions into ...  |
     |       ...Vector Array        |
     ------------------------------*/
    string oneline;
    vector<double> myVector;
    
    while (getline(dataIn,oneline)) {
        //cout << "Read in from file: " << oneline << endl;
        if ( oneline.empty() ) {
            //cout << "Theres a empty space in teh questions file..." << endl;
            continue;
        }
        double input = (double)atof(oneline.c_str());
        
        myVector.push_back(input);
    }
    
    dataIn.close();
    return myVector;
    
}

vector<string> getFiles()
{/* -----------------------------------------------------------------------------
  |  Function: getfiles()                                                   |
  |  Bio: Gets and stores questions from a text file.                           |
  ------------------------------------------------------------------------------*/
    
    /*------------------------------
     |       Opening File...        |
     |                              |
     ------------------------------*/
    ifstream dataIn("files.txt");
    if (!dataIn.is_open()) {
        cout << "\'Cannot Read Control File\' - Error in: getParameters()" << endl;
        cout << "Must Name Question File: \'Parameters.txt\' " << endl;
        abort();
    }
    
    /*------------------------------
     |    Place Questions into ...  |
     |       ...Vector Array        |
     ------------------------------*/
    string oneline;
    vector<string> myVector;
    
    while (getline(dataIn,oneline)) {
        //cout << "Read in from file: " << oneline << endl;
        if ( oneline.empty() ) {
            //cout << "Theres a empty space in teh questions file..." << endl;
            continue;
        }
        //double input = (double)atof(oneline.c_str());
        
        myVector.push_back(oneline.c_str());
    }
    
    dataIn.close();
    return myVector;
    
}

void pretty_plot(TH1F *h,
                 double line_color,
                 const char *histName, const char *axis_label_x, const char *axis_label_y,
                 const char *addSummary, const char* smallWindow,
                 int xmin, int ymin, int xmax, int ymax){
    
    const char option[10] = "add";
    
    h->SetLineColor(line_color);
    h->SetTitle(histName);
    h->GetXaxis()->SetTitle(axis_label_x); //min x, max x
    h->GetYaxis()->SetTitle(axis_label_y); //min x, max x
    h->SetMaximum(ymax);//max y
    h->SetMinimum(ymin);//min y
    h->GetXaxis()->SetRange(xmin,xmax); //min x, max x
    h->SetStats(0);
    h->Draw("same");
    
    if (addSummary == "add") {
        TPaveText *pave = new TPaveText(.65*xmax,0.3*ymax,.95*xmax,0.45*ymax);
        pave->SetFillColor(42);
        //TText *t1=pave->AddText("You can move");
        //t1->SetTextColor(4);
        //t1->SetTextSize(0.1);
        pave->AddText("Title and Stats pads");
        pave->AddText("X and Y axis");
        pave->AddText("You can modify bin contents");
        pave->Draw();
    }
    
    if (smallWindow == option) {
        TPad *pad = new TPad("food", histName, 0.6, 0.4, 0.85, 0.7);
        pad->Draw();
        pad->cd();
        h->SetStats(0);
        h->Draw();
    }
    
}

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
