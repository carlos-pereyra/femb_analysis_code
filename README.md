
# FEMB ADC Analysis Script Documentation

### Peak.h

contains analysis header (analyze.h) and script (test.cpp)

#### Class Analyze

##### set_run_channel
```c++
void Analyze::set_run_channel(int run, int channel){//1. declare for which subrun and channel this information pertains to
    wave.subrun = run;
    wave.channel = channel;
}
```
##### set_baseline
```c++
void Analyze::set_baseline(int subrun, int channel, std::vector<unsigned short> wf[64][128]){
    TString histName;
    histName.Form("Baseline Data (Subrun: %d, Channel: %d)", subrun, channel); // rename!!! -cp
    
    TH1F *base_hist = new TH1F("","",5000,1,5000); //histogram subrun 1 baseline #1
    
    for(int s = 0; s < wf[subrun][channel].size(); s++) base_hist->Fill(wf[subrun][channel].at(s));
    
    Double_t baseline = base_hist->GetMean(); //Get Baseline Mean
    vector_struct[subrun][channel].at(0).baseline = baseline;
    
}
```
##### set_peaks
```c++
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
```
##### set_pos_peak_mean
```c++
void Analyze::set_pos_peak_mean(int subrun, int channel){
    int sum = 0;
    TH1F *hpos = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).pos_npeaks; i++) hpos->Fill(vector_struct[subrun][channel].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak mean: " << hpos->GetMean() << endl;
    vector_struct[subrun][channel].at(0).pos_peak_mean = hpos->GetMean();
}
```
##### set_neg_peak_mean
```c++
void Analyze::set_neg_peak_mean(int subrun, int channel){
    int sum = 0;
    TH1F *hneg = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).neg_npeaks; i++) hneg->Fill(vector_struct[subrun][channel].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak mean: " << hneg->GetMean() << endl;
    vector_struct[subrun][channel].at(0).neg_peak_mean = hneg->GetMean();
}
```
##### set_pos_peak_rms
```c++
void Analyze::set_pos_peak_rms(int subrun, int channel){
    int sum = 0;
    TH1F *pos_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).pos_npeaks; i++) pos_rms->Fill(vector_struct[subrun][channel].at(0).pos_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting positive peak rms: " << pos_rms->GetRMS() << endl;
    vector_struct[subrun][channel].at(0).pos_peak_rms = pos_rms->GetRMS();
}
```
##### set_neg_peak_rms
```c++
void Analyze::set_neg_peak_rms(int subrun, int channel){
    int sum = 0;
    TH1F *neg_rms = new TH1F("","",6000,1,6000); // expected max peak = 4096
    for (int i = 0; i < vector_struct[subrun][channel].at(0).neg_npeaks; i++) neg_rms->Fill(vector_struct[subrun][channel].at(0).neg_peak_y.at(i));
    if(DBG) cout << "subrun: " << subrun << " channel: " << channel << "\tsetting negative peak rms: " << neg_rms->GetRMS() << endl;
    vector_struct[subrun][channel].at(0).neg_peak_rms = neg_rms->GetRMS();
}
```
##### set_run_info
```c++
void Analyze::set_run_info(int subrun, int channel) {
    vector_struct[subrun][channel].push_back(wave);
}
```
