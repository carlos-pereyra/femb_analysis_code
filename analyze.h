
const int DBG = 0;

class Analyze {
    
    public:
    
    struct waveform_struct{
        int subrun, channel, pos_npeaks, neg_npeaks, wf_elements;
        Double_t pos_peak_mean, pos_peak_rms, neg_peak_mean, neg_peak_rms, rise_time, width; //peak stats
        Double_t baseline;
        std::vector<Int_t> pos_peak_x;
        std::vector<Int_t> pos_peak_y;
        std::vector<Int_t> neg_peak_x;
        std::vector<Int_t> neg_peak_y;
    };
    waveform_struct wave;
    std::vector<waveform_struct> vector_struct[64][128]; //storage of structures
    
    //-----------------
    //----   set   ----
    //-----------------
    void set_run_channel(int, int);
    void set_run_info(int,int); //initializes vector of structures for specific subruns and channels

    void set_baseline(int,int,std::vector<unsigned short> wf[64][128]);
    
    // set peak mean & rms values
    void set_peaks(int,int,std::vector<unsigned short> wf[64][128], Int_t); // determine pos & neg peaks

    void set_pos_peak_mean(int, int);
    void set_neg_peak_mean(int,int);
    void set_pos_peak_rms(int,int);
    void set_neg_peak_rms(int,int);
    
    //-----------------
    //----   get   ----
    //-----------------
    int get_pos_npeaks(int,int); //full number of peaks (pos & neg)?
    int get_neg_npeaks(int,int); //full number of peaks (pos & neg)?

    // peak locations & amplitudes
    Int_t get_pos_peak_x(int,int,int);
    Int_t get_neg_peak_x(int,int,int);
    Double_t get_pos_peak_y(int,int,int);
    Double_t get_neg_peak_y(int,int,int);
    
    Double_t get_baseline(int,int);

    // set peak mean & rms values
    Double_t get_pos_peak_mean(int,int);
    Double_t get_neg_peak_mean(int,int);
    Double_t get_pos_peak_rms(int,int);
    Double_t get_neg_peak_rms(int,int);
    
    // peak quality
    Int_t get_risetime(int,int,std::vector<unsigned short> wf[64][128]);
    Int_t get_width(int,int,std::vector<unsigned short> wf[64][128]);
};

//----------------------------------------------------------------------------
//----                          set definitions                           ----
//----------------------------------------------------------------------------
void Analyze::set_run_channel(int run, int channel){//1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}

void Analyze::set_baseline(int subrun, int channel, std::vector<unsigned short> wf[64][128]){
    TString histName;
    histName.Form("Baseline Data (Subrun: %d, Channel: %d)", subrun, channel); // rename!!! -cp
    
    TH1F *base_hist = new TH1F("","",5000,1,5000); //histogram subrun 1 baseline #1
    
    for(int s = 0; s < wf[subrun][channel].size(); s++) base_hist->Fill(wf[subrun][channel].at(s));
    
    Double_t baseline = base_hist->GetMean(); //Get Baseline Mean
    vector_struct[subrun][channel].at(0).baseline = baseline;
    
}

void Analyze::set_peaks(int subrun, int channel, std::vector<unsigned short> wf[64][128], Int_t interval){//2. get peak locations
    TH1F *wfH = new TH1F("","ADC Signal",wf[subrun][channel].size(),0,wf[subrun][channel].size());
    
    for(int s = 0; s < wf[subrun][channel].size(); s++) wfH->SetBinContent(s , wf[subrun][channel].at(s)); //draw input waveform

    int stop = (wf[subrun][channel].size()/interval);
    
    std::vector<Int_t> p_peak_x;
    std::vector<Int_t> p_peak_y;
    std::vector<Int_t> n_peak_x;
    std::vector<Int_t> n_peak_y;
    p_peak_x.clear(); //setup vectors instead inside the definition
    p_peak_y.clear();
    n_peak_x.clear();
    n_peak_y.clear();
    
    if(DBG) cout << "\n================================================================" << endl;
    if(DBG) cout << "\tSetting Low & High Peaks (Subrun: " << subrun << " ,Channel: " << channel << ")\n" << endl;
    if(DBG) cout << "\tX_High\tHigh_Peak\tX_Low\tLow_Peak"<< endl;
    if(DBG) cout << "\t======\t=========\t=====\t========"<< endl;
    
    Double_t baseline = vector_struct[subrun][channel].at(0).baseline;
    for (int i = 0; i < stop; i++) {
        wfH->GetXaxis()->SetRange(interval*(i),interval*(i+1)); // problem with wfH
        
        Double_t binMin = wfH->GetMinimumBin();
        Double_t currentMin = wfH->GetMinimum();
        Double_t binMax = wfH->GetMaximumBin();
        Double_t currentMax = wfH->GetMaximum();
    
        if (currentMax>(baseline+60)) {
            p_peak_x.push_back(binMax);
            vector_struct[subrun][channel].at(0).pos_peak_x.push_back(binMax);
            
            p_peak_y.push_back(currentMax);
            vector_struct[subrun][channel].at(0).pos_peak_y.push_back(currentMax);
        }
        
        if (currentMin<(baseline-60)) {
            n_peak_x.push_back(binMin);
            vector_struct[subrun][channel].at(0).neg_peak_x.push_back(binMin);
            
            n_peak_y.push_back(currentMin);
            vector_struct[subrun][channel].at(0).neg_peak_y.push_back(currentMin);
        }
    }
    if(DBG) for (int i = 0; i < 5; i++) cout << "\t" << vector_struct[subrun][channel].at(0).pos_peak_x.at(i) << "\t   " << vector_struct[subrun][channel].at(0).pos_peak_y.at(i) << "\t\t" << vector_struct[subrun][channel].at(0).neg_peak_x.at(i) << "\t"  << vector_struct[subrun][channel].at(0).neg_peak_y.at(i) << endl;
    if(DBG) cout << "\n" << endl;
    
    wfH->GetXaxis()->SetRange(1,wf[subrun][channel].size());
    
    vector_struct[subrun][channel].at(0).pos_npeaks = vector_struct[subrun][channel].at(0).pos_peak_y.size(); // not finished after equals sign
    vector_struct[subrun][channel].at(0).neg_npeaks = vector_struct[subrun][channel].at(0).neg_peak_y.size(); // not finished after equals sign
}

void Analyze::set_pos_peak_mean(int subrun, int channel){
    int sum = 0;
    TH1F *hpos = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).pos_npeaks; i++) hpos->Fill(vector_struct[subrun][channel].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak mean: " << hpos->GetMean() << endl;
    vector_struct[subrun][channel].at(0).pos_peak_mean = hpos->GetMean();
}

void Analyze::set_neg_peak_mean(int subrun, int channel){
    int sum = 0;
    TH1F *hneg = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).neg_npeaks; i++) hneg->Fill(vector_struct[subrun][channel].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak mean: " << hneg->GetMean() << endl;
    vector_struct[subrun][channel].at(0).neg_peak_mean = hneg->GetMean();
}

void Analyze::set_pos_peak_rms(int subrun, int channel){
    int sum = 0;
    TH1F *pos_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).pos_npeaks; i++) pos_rms->Fill(vector_struct[subrun][channel].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak rms: " << pos_rms->GetRMS() << endl;
    vector_struct[subrun][channel].at(0).pos_peak_rms = pos_rms->GetRMS();
}

void Analyze::set_neg_peak_rms(int subrun, int channel){
    int sum = 0;
    TH1F *neg_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).neg_npeaks; i++) neg_rms->Fill(vector_struct[subrun][channel].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak rms: " << neg_rms->GetRMS() << endl;
    vector_struct[subrun][channel].at(0).neg_peak_rms = neg_rms->GetRMS();
}

void Analyze::set_run_info(int subrun, int channel) {
    vector_struct[subrun][channel].push_back(wave);
}

//----------------------------------------------------------------------------
//----                          get definitions                           ----
//----------------------------------------------------------------------------
int Analyze::get_pos_npeaks(int subrun,int channel){ // number of positive peaks only
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.pos_npeaks;
}

int Analyze::get_neg_npeaks(int subrun,int channel){ // number of positive peaks only
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.neg_npeaks;
}

Int_t Analyze::get_pos_peak_x(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.pos_peak_x.at(element);
}

Int_t Analyze::get_neg_peak_x(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.neg_peak_x.at(element);
}

Double_t Analyze::get_pos_peak_y(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    //if(DBG) cout << "\tGetting positive peak: " << w.pos_peak_y.at(element) << "\t at bin: " << w.pos_peak_x.at(element) << endl;
    if (element < w.pos_peak_x.size()) return w.pos_peak_y.at(element);
    else {
        cout << "requested element is outside positive peak matrix size (Analyze::get_pos_peak_y)" << endl;
        return 0;
    }
}

Double_t Analyze::get_neg_peak_y(int subrun,int channel,int element){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    //if(DBG) cout << "\tGetting negative peak: " << w.neg_peak_y.at(element) << "\t at bin: " << w.neg_peak_x.at(element) << endl;
    if (element < w.neg_peak_x.size()) return w.neg_peak_y.at(element);
    else {
        cout << "requested element is outside positive peak matrix size (Analyze::get_neg_peak_y)" << endl;
        return 0;
    }
}

Double_t Analyze::get_baseline(int subrun, int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    return w.baseline;
}

Double_t Analyze::get_pos_peak_mean(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting positive peak mean: " << w.pos_peak_mean << endl;
    return w.pos_peak_mean;
}

Double_t Analyze::get_neg_peak_mean(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting negative peak mean: " << w.neg_peak_rms << endl;
    return w.neg_peak_mean;
}

Double_t Analyze::get_pos_peak_rms(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting positive peak rms: " << w.pos_peak_rms << endl;
    return w.pos_peak_rms;
}

Double_t Analyze::get_neg_peak_rms(int subrun,int channel){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting negative peak rms: " << w.neg_peak_rms << endl;
    return w.neg_peak_rms;
}

Int_t Analyze::get_risetime(int subrun, int channel, std::vector<unsigned short> wf[64][128]){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    float ratio;
    Double_t y1 = (w.pos_peak_y.at(1)-w.baseline), y2;
    Int_t risetime;
    
    for (int i = 0; i < 10; i++) {
        y2 = (wf[subrun][channel].at(w.pos_peak_x.at(1)-i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.01) break;
        risetime = i+1;
        
        if(DBG) cout << "ratio: " << ratio << endl;
        if(DBG) cout <<"x1: "<< w.pos_peak_x.at(1) <<" y1: "<<y1<< " x2: " << risetime << " y2: " << y2 << endl;
    }
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting risetime: " << risetime << endl;
    vector_struct[subrun][channel].at(0).rise_time = risetime;
    return risetime;
}

Int_t Analyze::get_width(int subrun, int channel, std::vector<unsigned short> wf[64][128]){
    waveform_struct w = vector_struct[subrun][channel].at(0);
    float ratio;
    Double_t y1 = (w.pos_peak_y.at(0)-w.baseline), y2;
    Int_t width, L, R;
    
    for (int i = 0; i < 10; i++) {
        y2 = (wf[subrun][channel].at(w.pos_peak_x.at(0)-i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.02) break;
        L = i+1;
        
        if(DBG) cout << "L ratio: " << ratio << endl;
        if(DBG) cout <<"x1: "<< w.pos_peak_x.at(0) <<" y1: "<<y1<< " x2: " << L << " y2: " << y2 << endl;
    }
    
    for (int i = 0; i < 10; i++) {
        y2 = (wf[subrun][channel].at(w.pos_peak_x.at(0)+i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.02) break;
        R = i+1;
        
        if(DBG) cout << "R ratio: " << ratio << endl;
        if(DBG) cout <<"x1: "<< w.pos_peak_x.at(0) <<" y1: "<<y1<< " x2: " << R << " y2: " << y2 << endl;
    }
    
    width = L + R + 1;
    
    //if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting risetime: " << risetime << endl;
    vector_struct[subrun][channel].at(0).width = width;
    return width;
}




