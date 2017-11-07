
const int DBG = 0;

class Analyze {
    
    public:
    
    struct waveform_struct{
        int subrun, channel, pos_npeaks, neg_npeaks, wf_elements;
        Double_t pos_peak_mean, pos_peak_rms, neg_peak_mean, neg_peak_rms, rise_time, width; //peak stats
        Double_t baseline, baseline_rms;
        std::vector<Int_t> pos_peak_x;
        std::vector<Int_t> pos_peak_y;
        std::vector<Int_t> neg_peak_x;
        std::vector<Int_t> neg_peak_y;
        
        std::vector<Double_t> fit_gain;
        std::vector<Double_t> fit_risetime;
        std::vector<Double_t> fit_quality;
    };
    waveform_struct wave;
    std::vector<waveform_struct> vector_of_struct[64][128][5]; //storage of structures
    
    //
    // set
    //
    void set_run_channel(int, int);
    void set_run_info(int,int,int); //initializes vector of structures for specific subruns and channels

    void set_baseline(int,int,int,std::vector<unsigned short> wf[64][128][10]);
    void set_baseline_rms(int,int,int,std::vector<unsigned short> wf[64][128][10]);
    
    // set peak mean & rms values
    void set_peaks(int,int,int,std::vector<unsigned short> wf[64][128][10],Int_t); // determine pos & neg peaks

    void set_pos_peak_mean(int,int,int);
    void set_neg_peak_mean(int,int,int);
    void set_pos_peak_rms(int,int,int);
    void set_neg_peak_rms(int,int,int);
    
    void set_fit_gn(int,int,int,int,int);
    void set_fit_rt(int,int,int,int,double);
    void set_fit_quality(int,int,int,int,double);
    
    //
    // get
    //
    int get_pos_npeaks(int,int,int); //full number of peaks (pos & neg)?
    int get_neg_npeaks(int,int,int); //full number of peaks (pos & neg)?

    // peak locations & amplitudes
    Int_t get_pos_peak_x(int,int,int,int);
    Int_t get_neg_peak_x(int,int,int,int);
    Double_t get_pos_peak_y(int,int,int,int);
    Double_t get_neg_peak_y(int,int,int,int);
    
    Double_t get_baseline(int,int,int);
    Double_t get_baseline_rms(int,int,int);

    // set peak mean & rms values
    Double_t get_pos_peak_mean(int,int,int);
    Double_t get_neg_peak_mean(int,int,int);
    Double_t get_pos_peak_rms(int,int,int);
    Double_t get_neg_peak_rms(int,int,int);
    
    // peak quality
    Int_t get_risetime(int,int,int,int,std::vector<unsigned short> wf[64][128][10]);
    Int_t get_width(int,int,int,int,std::vector<unsigned short> wf[64][128][10]);
    
    Int_t get_fit_gn(int,int,int,int);
    Double_t get_fit_rt(int,int,int,int);
    Double_t get_fit_quality(int,int,int,int);
};

//----------------------------------------------------------------------------
//----                          set definitions                           ----
//----------------------------------------------------------------------------
void Analyze::set_run_channel(int run, int channel){//1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}

void Analyze::set_baseline(int subrun, int channel, int file, std::vector<unsigned short> wf[64][128][10]){
    TString histName;
    histName.Form("Baseline MEAN (Subrun: %d, Channel: %d)", subrun, channel);
    TH1F *base_hist = new TH1F("","",5000,1,5000);
    for(int s = 0; s < 100; s++) base_hist->Fill(wf[subrun][channel][file].at(s));
    Double_t baseline = base_hist->GetMean();
    vector_of_struct[subrun][channel][file].at(0).baseline = baseline;
}

void Analyze::set_baseline_rms(int subrun, int channel, int file, std::vector<unsigned short> wf[64][128][10]){
    TString histName;
    histName.Form("Baseline RMS (Subrun: %d, Channel: %d)", subrun, channel);
    TH1F *base_hist = new TH1F("","",5000,1,5000);
    for(int s = 0; s < 100; s++) base_hist->Fill(wf[subrun][channel][file].at(s));
    Double_t baseline_rms = base_hist->GetRMS();
    vector_of_struct[subrun][channel][file].at(0).baseline_rms = baseline_rms;
}

void Analyze::set_peaks(int subrun, int channel, int file, std::vector<unsigned short> wf[64][128][10], Int_t interval){
    
    TH1F *wfH = new TH1F("","ADC Signal",wf[subrun][channel][file].size(),0,wf[subrun][channel][file].size());
    
    for(int s = 0; s < wf[subrun][channel][file].size(); s++) wfH->SetBinContent(s , wf[subrun][channel][file].at(s));

    int stop = (wf[subrun][channel][file].size()/interval);
    
    if(DBG) cout << "\n\tSetting Low & High Peaks (Subrun: " << subrun << " ,Channel: " << channel << ")\n" << endl;
    if(DBG) cout << "\tX_High\tHigh_Peak\tX_Low\tLow_Peak"<< endl;
    if(DBG) cout << "\t======\t=========\t=====\t========"<< endl;
    
    Int_t baseline = vector_of_struct[subrun][channel][file].at(0).baseline;
    for (int i = 0; i < stop; i++) {
        wfH->GetXaxis()->SetRange(interval*(i),interval*(i+1));
        
        Int_t x_min = wfH->GetMinimumBin()-1; // off by +1 for some reason...
        Int_t y_min = wfH->GetMinimum();
        
        Int_t x_max = wfH->GetMaximumBin()-1; // off by +1 for some reason...
        Int_t y_max = wfH->GetMaximum();
        
        if (y_min<(baseline-100)) {
            vector_of_struct[subrun][channel][file].at(0).neg_peak_x.push_back(x_min);
            vector_of_struct[subrun][channel][file].at(0).neg_peak_y.push_back(std::abs(y_min-baseline));
        }
        
        if (y_max>(baseline+100)) {
            vector_of_struct[subrun][channel][file].at(0).pos_peak_x.push_back(x_max);
            vector_of_struct[subrun][channel][file].at(0).pos_peak_y.push_back(std::abs(y_max-baseline));
        }
        
    }
    
    if(DBG) for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).pos_peak_x.size(); i++) cout << "pos. x peak:\t" << vector_of_struct[subrun][channel][file].at(0).pos_peak_x.at(i) << endl;
    if(DBG) for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).pos_peak_y.size(); i++) cout << "pos. y peak:\t" << vector_of_struct[subrun][channel][file].at(0).pos_peak_y.at(i) << endl;
    
    if(DBG) for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).neg_peak_x.size(); i++) cout << "neg. x peak:\t" << vector_of_struct[subrun][channel][file].at(0).neg_peak_x.at(i) << endl;
    if(DBG) for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).neg_peak_y.size(); i++) cout << "neg. y peak:\t" << vector_of_struct[subrun][channel][file].at(0).neg_peak_y.at(i) << endl;
    
    wfH->GetXaxis()->SetRange(1,wf[subrun][channel][file].size());
    
    vector_of_struct[subrun][channel][file].at(0).pos_npeaks = vector_of_struct[subrun][channel][file].at(0).pos_peak_y.size();
    vector_of_struct[subrun][channel][file].at(0).neg_npeaks = vector_of_struct[subrun][channel][file].at(0).neg_peak_y.size();
    
}

void Analyze::set_pos_peak_mean(int subrun, int channel, int file){
    int sum = 0;
    TH1F *hpos = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).pos_npeaks; i++) hpos->Fill(vector_of_struct[subrun][channel][file].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak mean: " << hpos->GetMean() << endl;
    vector_of_struct[subrun][channel][file].at(0).pos_peak_mean = hpos->GetMean();
}

void Analyze::set_neg_peak_mean(int subrun, int channel, int file){
    int sum = 0;
    TH1F *hneg = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).neg_npeaks; i++) hneg->Fill(vector_of_struct[subrun][channel][file].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak mean: " << hneg->GetMean() << endl;
    vector_of_struct[subrun][channel][file].at(0).neg_peak_mean = hneg->GetMean();
}

void Analyze::set_pos_peak_rms(int subrun, int channel, int file){
    int sum = 0;
    TH1F *pos_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).pos_npeaks; i++) pos_rms->Fill(vector_of_struct[subrun][channel][file].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak rms: " << pos_rms->GetRMS() << endl;
    vector_of_struct[subrun][channel][file].at(0).pos_peak_rms = pos_rms->GetRMS();
}

void Analyze::set_neg_peak_rms(int subrun, int channel, int file){
    int sum = 0;
    TH1F *neg_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_of_struct[subrun][channel][file].at(0).neg_npeaks; i++) neg_rms->Fill(vector_of_struct[subrun][channel][file].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak rms: " << neg_rms->GetRMS() << endl;
    vector_of_struct[subrun][channel][file].at(0).neg_peak_rms = neg_rms->GetRMS();
}

void Analyze::set_fit_gn(int subrun, int channel, int file, int peak, int gain){
    vector_of_struct[subrun][channel][file].at(0).fit_gain.resize(100); //100 elements
    vector_of_struct[subrun][channel][file].at(0).fit_gain.at(peak) = gain;
}

void Analyze::set_fit_rt(int subrun, int channel, int file, int peak, double risetime){
    vector_of_struct[subrun][channel][file].at(0).fit_risetime.resize(100);
    vector_of_struct[subrun][channel][file].at(0).fit_risetime.at(peak) = risetime;
}

void Analyze::set_fit_quality(int subrun, int channel, int file, int peak, double quality){
    vector_of_struct[subrun][channel][file].at(0).fit_quality.resize(100);
    vector_of_struct[subrun][channel][file].at(0).fit_quality.at(peak) = quality;
}

void Analyze::set_run_info(int subrun, int channel, int file) {
    vector_of_struct[subrun][channel][file].push_back(wave);
}

//----------------------------------------------------------------------------
//----                          get definitions                           ----
//----------------------------------------------------------------------------
int Analyze::get_pos_npeaks(int subrun, int channel, int file){ // number of positive peaks only
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    return w.pos_npeaks;
}

int Analyze::get_neg_npeaks(int subrun,int channel,int file){ // number of positive peaks only
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    return w.neg_npeaks;
}

Int_t Analyze::get_pos_peak_x(int subrun,int channel,int file,int element){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if(DBG) if (element>w.pos_npeaks) cout << "\nAnalyze::get_pos_peak_x: outside of peaks_x matrix (Error)" << endl;
    if(DBG) if (element>w.pos_npeaks) cout << "\tSubrun: " << subrun << " Channel: " << channel << "\n" << endl;
    if (element>w.pos_npeaks) return 0;
    
    return w.pos_peak_x.at(element);
}

Int_t Analyze::get_neg_peak_x(int subrun,int channel,int file,int element){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if(DBG) if (element>w.neg_npeaks) cout << "\nAnalyze::get_neg_peak_x: outside of peaks_x matrix (Error)" << endl;
    if(DBG) if (element>w.neg_npeaks) cout << "\tSubrun: " << subrun << " Channel: " << channel << "\n" << endl;
    if (element>w.neg_npeaks) return 0;
    
    return w.neg_peak_x.at(element);
}

Double_t Analyze::get_pos_peak_y(int subrun,int channel,int file,int element){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if (element < w.pos_peak_x.size()) return w.pos_peak_y.at(element);
    else {
        //cout << "requested element is outside positive peak matrix size (Analyze::get_pos_peak_y)" << endl;
        return 0;
    }
}

Double_t Analyze::get_neg_peak_y(int subrun,int channel,int file,int element){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if (element < w.neg_peak_x.size()) return w.neg_peak_y.at(element);
    else {
        cout << "requested element is outside positive peak matrix size (Analyze::get_neg_peak_y)" << endl;
        return 0;
    }
}

Double_t Analyze::get_baseline(int subrun, int channel, int file){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    return w.baseline;
}

Double_t Analyze::get_baseline_rms(int subrun, int channel, int file){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    return w.baseline_rms;
}

Double_t Analyze::get_pos_peak_mean(int subrun,int channel,int file){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting positive peak mean: " << w.pos_peak_mean << endl;
    return w.pos_peak_mean;
}

Double_t Analyze::get_neg_peak_mean(int subrun,int channel,int file){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting negative peak mean: " << w.neg_peak_rms << endl;
    return w.neg_peak_mean;
}

Double_t Analyze::get_pos_peak_rms(int subrun,int channel,int file){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting positive peak rms: " << w.pos_peak_rms << endl;
    return w.pos_peak_rms;
}

Double_t Analyze::get_neg_peak_rms(int subrun,int channel,int file){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting negative peak rms: " << w.neg_peak_rms << endl;
    return w.neg_peak_rms;
}

Int_t Analyze::get_risetime(int subrun, int channel, int file, int peak, std::vector<unsigned short> wf[64][128][10]){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    //if (element>w.pos_npeaks) cout << "\nAnalyze::get_risetime: outside of peaks matrix (Error)" << endl;
    //if (element>w.pos_npeaks) cout << "\tSubrun: " << subrun << " Channel: " << channel << "\n" << endl;
    if (peak >= w.pos_npeaks) return 0;
    if (w.pos_peak_x.at(peak) == 0) return 0;
    if (w.pos_peak_y.at(peak) == 0) return 0;
    
    float ratio;
    Double_t y1 = (w.pos_peak_y.at(peak));
    Double_t y2;
    Int_t risetime;
    
    for (int i = 0; i < 10; i++) {
        y2 = (wf[subrun][channel][file].at(w.pos_peak_x.at(peak)-i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.005) break;
        risetime = i;
        
        if(DBG) cout << "\nratio: " << ratio << endl;
        if(DBG) cout << "x1: "<< w.pos_peak_x.at(peak) <<" y1: "<<y1<< " x2: " << risetime << " y2: " << y2 << endl;
    }
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting risetime: " << risetime << endl;
    vector_of_struct[subrun][channel][file].at(0).rise_time = risetime;
    return risetime;
}

Int_t Analyze::get_width(int subrun, int channel, int file, int peak, std::vector<unsigned short> wf[64][128][10]){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    
    if (w.pos_npeaks <= peak) return 0;
    if (w.pos_peak_x.at(peak) == 0) return 0;
    if (w.pos_peak_y.at(peak) == 0) return 0;
    
    float ratio;
    Double_t y1 = (w.pos_peak_y.at(peak)), y2;
    Int_t width, L, R;
    
    for (int i = 0; i < 25; i++) {
        y2 = (wf[subrun][channel][file].at(w.pos_peak_x.at(peak)-i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.001) break;
        L = i;
        
        if(DBG) cout << "L ratio: " << ratio << endl;
        if(DBG) cout << "x1: " << w.pos_peak_x.at(peak) << " y1: " <<y1<< " x2: " << L << " y2: " << y2 << endl;
    }
    
    for (int i = 0; i < 25; i++) {
        y2 = (wf[subrun][channel][file].at(w.pos_peak_x.at(peak)+i)-w.baseline);
        ratio = y2/y1;
        if (ratio<0.001) break;
        R = i;
        
        if(DBG) cout << "R ratio: " << ratio << endl;
        if(DBG) cout <<"x1: "<< w.pos_peak_x.at(peak) <<" y1: "<<y1<< " x2: " << R << " y2: " << y2 << endl;
    }
    
    width = L + R;
    
    //if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tgetting risetime: " << risetime << endl;
    vector_of_struct[subrun][channel][file].at(0).width = width;
    
    return width;
}

Int_t Analyze::get_fit_gn(int subrun, int channel, int file, int peak){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if (w.fit_gain.size()==0) return 0;
    Int_t gain = w.fit_gain.at(peak);
    return gain;
}

Double_t Analyze::get_fit_rt(int subrun, int channel, int file, int peak){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if (w.fit_risetime.size()==0) return 0;
    Double_t risetime = w.fit_risetime.at(peak);
    return risetime;
}

Double_t Analyze::get_fit_quality(int subrun, int channel, int file, int peak){
    waveform_struct w = vector_of_struct[subrun][channel][file].at(0);
    if (w.fit_quality.size()==0) {
        cout << "nothing inside: get_fit_quality: ";
        return 0;
    }
    Double_t quality = w.fit_quality.at(peak);
    return quality;
}
