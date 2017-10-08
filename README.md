
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

```
##### set_pos_peak_mean
```c++

```
##### set_neg_peak_mean
```c++

```
##### set_pos_peak_rms
```c++

```
##### set_neg_peak_rms
```c++

```
##### set_run_info
