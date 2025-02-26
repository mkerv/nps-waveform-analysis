// This macro is used for the fit of the production runs not the elastic runs

 #include <iostream>
 #include <TObject.h>
 #include <TROOT.h>
 #include "TCanvas.h"
 #include "TGraphErrors.h"
 #include "TF1.h"
 #include "TH1F.h"
 #include "TH2F.h"
 #include "TLegend.h"
 #include "TMath.h"
 #include "TLorentzVector.h"
 #include "TChain.h"
 #include "TTree.h"
 #include "TStyle.h"
 #include "TSystem.h"
 #include "TRandom.h"
 #include "TFile.h"
 #include "TMatrixD.h"
 #include "TLatex.h"
 #include "TVirtualFFT.h"
 #include "TStopwatch.h"
 #include <Math/Interpolator.h>
 #include <fstream>
#include <TObjectTable.h>

 
const int ntime=110; // number of time samples for each fADC channel (100 for runs 55 and 56)
const int nfADC=7; //number of fADC modules (16 channels each) used in this run
const int nchannel=16; //number of channels in each fADC module
const int ncol=30; //number of calo colomns tested in this run (3 for runs 55 and 56)
const int nlin=36; //number of calo blocks in each colomn
const int nblocks=ncol*nlin; //number of tested calo blocks
const int maxpulses=4; //nb maximal de pulses dans la fentre ntime 
const int nsignals=nblocks*maxpulses;
const int maxwfpulses=12; //nb maximal de pulses que la wfa peut trouver
const int ntemp=56; //nb of temp sensors

Int_t wfnpulse[nblocks];
Double_t wfampl[nblocks][maxwfpulses]; 
Double_t wftime[nblocks][maxwfpulses];
Double_t timeref[nblocks];
Double_t timerefacc=-4.5;//car les bons pulses elastiques poussent vers 45.5 (4*ns) dans la wf alors que ceux des production runs sont vers 35.5 (4ns)
Double_t timerefacc2=0;// en (-4ns) car le pic de timewf pousse vers -22ns


ROOT::Math::Interpolator *interpolation[nblocks]; //function used for the interpolation
TH1F *hsig_i[nblocks];//histogram showing the waveform for block i and a given event
TF1 *finter[nblocks];


// FIT FUNCTION

Double_t func(Double_t *x, Double_t *par){ 
 Int_t j= (int)(par[0]);
 Double_t val=0;
 for(Int_t p=0;p<maxwfpulses;p++){
 if(x[0]-par[2+2*p]>1&&x[0]-par[2+2*p]<109)val+=par[3+2*p]*interpolation[j]->Eval(x[0]-par[2+2*p]);
 }
 return val+par[1];
}

///////////////////////////////////////////////////////////// FIT FUNCTION USED AND SCHEMATICS /////////////////////////////////////////////////////////////


void Fitwf(Int_t bn){
  
//Detect the pulses

wfnpulse[bn]=0;
Int_t good=0;
 
for(Int_t p=0;p<maxwfpulses;p++){ 

  //cout << "Setting parameters for pulse " << p << " with index " << 2+2*p << " and " << 3+2*p << endl;
 
wfampl[bn][p]=-1;
 }


for (Int_t it=0;it<ntime-6;it++){

 
if (!finter[bn]) {
   cout << " block=" << bn << " finter is nullptr" <<endl;
    return;
}

  
if (!hsig_i[bn]) {
    cout<< ", block=" << bn << " hsig_i is nullptr" <<endl;
    return;
}

  
// Condition over the number of samples in the pulse finding scheme
  
if(hsig_i[bn]->GetBinContent(it+1)<hsig_i[bn]->GetBinContent(it+2)&&hsig_i[bn]->GetBinContent(it+2)<hsig_i[bn]->GetBinContent(it+3)&&hsig_i[bn]->GetBinContent(it+3)<=hsig_i[bn]->GetBinContent(it+4)&&hsig_i[bn]->GetBinContent(it+4)>=hsig_i[bn]->GetBinContent(it+5)){

  
if(hsig_i[bn]->GetBinContent(it+4)>0){

//check if we exceeded the number of pulses
  
if (wfnpulse[bn] >= maxwfpulses) {
    cout << "Warning: wfnpulse[" << bn << "] exceeded maxwfpulses!" << endl;
    //  wfnpulse[bn] = maxwfpulses - 1; // Prevent overflow
}

 wfampl[bn][wfnpulse[bn]]=hsig_i[bn]->GetBinContent(it+4);// get the amplitude of the pulse found

 wftime[bn][wfnpulse[bn]]=hsig_i[bn]->GetBinCenter(it+4);// get the time of the pulse found

// flag for the good pulse
 
 if(TMath::Abs(wftime[bn][wfnpulse[bn]]-timeref[bn]-timerefacc)<4.1){
   good=1;
 }

 wfnpulse[bn]++;

// to prevent overflow
 if(wfnpulse[bn]==maxwfpulses){
   wfnpulse[bn]=maxwfpulses-1;
   it=ntime;
 }

 it+=4;

 }// end of the condition over (hsig_i[bn]->GetBinContent(it+4)>0){
 }// end of the condition over samples
 }// end of loop over it

 
 if(wfnpulse[bn]>maxwfpulses-2){cout<<"Warning : nb de pulses excessivement eleve dans la wf du bloc "<<bn<<endl;}

 //Adjust the parameters of the fit function



 //File<<"c_0"<<endl;
 
for(Int_t p=0;p<maxwfpulses;p++)
  {
    // File<<"c1"<<endl;
    finter[bn]->FixParameter(2+2*p,0.);
    //File<<"c2"<<endl;
    finter[bn]->FixParameter(3+2*p,0.);
  }

if(wfnpulse[bn]>0&&good==1){

for(Int_t p=0;p<TMath::Min(maxwfpulses,wfnpulse[bn]);p++){
  
finter[bn]->ReleaseParameter(2+2*p);
finter[bn]->ReleaseParameter(3+2*p);

finter[bn]->SetParameter(2+2*p,wftime[bn][p]-timeref[bn]);
finter[bn]->SetParameter(3+2*p,wfampl[bn][p]);

finter[bn]->SetParLimits(2+2*p,wftime[bn][p]-timeref[bn]-3,wftime[bn][p]-timeref[bn]+3);
finter[bn]->SetParLimits(3+2*p,wfampl[bn][p]*0.2,wfampl[bn][p]*3);
}}

//File<<"c3"<<endl;
 
if(wfnpulse[bn]>0&&good==0){
  
for(Int_t p=0;p<TMath::Min(maxwfpulses,wfnpulse[bn]);p++){

finter[bn]->ReleaseParameter(2+2*p);
finter[bn]->ReleaseParameter(3+2*p);

finter[bn]->SetParameter(2+2*p,wftime[bn][p]-timeref[bn]);

//Check to debug
if (wfampl[bn][p] > 0) { // Ensure it's positive before multiplying
    finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn][p] * 0.2, wfampl[bn][p] * 3);
} else {
    cout << "Warning: wfampl[" << bn << "][" << p << "] is non-positive!" << endl;
    // finter[bn]->SetParLimits(3 + 2 * p, 0.05, 10); // Set safe fallback limits
}

finter[bn]->SetParameter(3+2*p,wfampl[bn][p]);
 
finter[bn]->SetParLimits(2+2*p,wftime[bn][p]-timeref[bn]-3,wftime[bn][p]-timeref[bn]+3);
finter[bn]->SetParLimits(3+2*p,wfampl[bn][p]*0.2,wfampl[bn][p]*3);
}
 
//File<<"c4"<<endl;

//On recherche quand meme un eventuel pulse en temps 

finter[bn]->ReleaseParameter(2+2*wfnpulse[bn]);
finter[bn]->ReleaseParameter(3+2*wfnpulse[bn]);

finter[bn]->SetParameter(2+2*wfnpulse[bn],timerefacc);
finter[bn]->SetParameter(3+2*wfnpulse[bn],2);

finter[bn]->SetParLimits(2+2*wfnpulse[bn],timerefacc-4,timerefacc+4);
finter[bn]->SetParLimits(3+2*wfnpulse[bn],0.05,10);

wfnpulse[bn]++;

}
 
//File<<"c5"<<endl;

if(wfnpulse[bn]==0){

for(Int_t p=0;p<1;p++){

finter[bn]->ReleaseParameter(2+2*p);
finter[bn]->ReleaseParameter(3+2*p);
 
finter[bn]->SetParameter(2+2*p,timerefacc);
finter[bn]->SetParameter(3+2*p,2);
 
finter[bn]->SetParLimits(2+2*p,timerefacc-4,timerefacc+4);
finter[bn]->SetParLimits(3+2*p,0.05,10);

}
 
wfnpulse[bn]++;

}

//File<<"c6"<<endl;
 
finter[bn]->SetParameter(1,0.);
finter[bn]->SetParLimits(1,-100,100.);

//File<<"c7"<<endl;
//hsig_i[bn]->Fit(Form("finter_%d",bn),"Q","",TMath::Max(0.,wftime[bn][0]-20),ntime);
//File << "Checkpoint: Before fitting function for bn = " << bn << endl;

hsig_i[bn]->Fit(Form("finter_%d",bn),"Q","",0.,ntime);
 
//File<<"c8"<<endl;
}

/////////////////////////////////////////////////////////////MAIN FUNCTION ///////////////////////////////////////////////////////////////////////

void TEST_2(int run, int seg){

   int nthreads = 8;
 
   ROOT::EnableImplicitMT(nthreads);



std::cout << "Implicit MT enabled: " << ROOT::IsImplicitMTEnabled() << std::endl;
std::cout << "Number of threads: " << ROOT::GetThreadPoolSize() << std::endl;

  
   gROOT->Initialize();

    // Before processing events, check the list of objects in memory
    std::cout << "Objects in memory before event processing:" << std::endl;
    gObjectTable->Print();




  
  if(run==2016){timerefacc=-10;}//keep this condition and might add more runs if that happens again

  // timerefacc = timerefacc + 10;
  // timerefacc2 = timerefacc2 + 10 * 4;

  //int end_event, start_event;
  //const int evts_per_job = 250000;
  Double_t dt=4.; //time bin (sample) width (4 ns), the total time window is then ntime*dt
  const int nslots=1104; //nb maximal de slots dans tous les fADC
  const int Ndata=nslots*(ntime+1+1); //(1104 slots fADC au total mais pas tous utilises y compris 2 PM scintillateurs?) should correspond to Ndata.NPS.cal.fly.adcSampWaveform variable in the root file
  Double_t ADCtomV=1000./4096; //1000 mV correspond a 4096 canaux
  Double_t integtopC= (1./4096)*(dt*1.e3/50); // (1 V / 4096 adc channels) * (4000 ps time sample / 50 ohms input resistance) = 0.020 pc/channel 

TStopwatch t;
t.Start();

// where the cosmic pulse is expected to be
Int_t binmin;
Int_t binmax;
long nentries;

// Input root file

 TChain *tree=new TChain("T");

 //TString filename = Form("/volatile/hallc/nps/nps-ana/ROOTfiles/COIN/PROD/nps_hms_coin_%d_%d_1_-1.root", run, seg);

 // TString filename = Form("/cache/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", run, seg);

   TString filename = Form("/cache/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", run, seg);//new files with the HMS calibration
 
 //TString filename = Form("/volatile/hallc/nps/wassim/nps_hms_coin_%d_%d_1_-1.root", run, seg);//new files with the HMS calibration just for the first setting runs

 //TString filename = Form("/mss/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", run, seg);
 
 TFile *file = TFile::Open(filename);

 if (file && !file->IsZombie()) {
   tree->Add(filename);
   nentries = tree->GetEntries();
   cout << "There are:  " << nentries << " in this segment.\n";
   file->Close();
 } else {
   // File doesn't exist or is corrupted, move to the next rootfile
   cout << "File does not exist !! " << filename << endl;
 }

//Read reference waveforms

const int nbparameters=2+2*maxwfpulses;
cout<<"nbparameters = "<<nbparameters<<endl;
Double_t x[ntime],y[ntime];
ifstream filewf;
Double_t dum1,ymax;

for(Int_t i=0; i<nblocks; i++){

  //  filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/ref_wf_%d.txt",i));
  //filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/ref_wf_%d.txt",i));
      
    if (run>6183 && run<7500) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6171-6183/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>6168 && run<6171) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>5236 && run<6151) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>5208 && run<5217) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5183-5208/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>3898 && run<5183) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>2920 && run<3883) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2900-2920/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>2885 && run<2900) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2875-2885/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>2871 && run<2875) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2855-2871/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>1982 && run<2855) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/ref_wf_%d.txt",i));
    } 
    else if (run>1560 && run<1961) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1423-1511/RWF/ref_wf_%d.txt",i));
    }
      
   //comparison with the mean of the ref shapes
  //filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/fit_elastic_runs/results_elastics/refwf/final_refwf_block_%d.txt",i)); Done
  
timeref[i]=-1.e6;
if(filewf.is_open()){
  filewf>>timeref[i]>>dum1;//used for my ref shapes (be careful in switching on the mean method !!!!!!!!!!!!!!!!!)
  
for (Int_t it=0;it<ntime;it++){filewf>>x[it]>>y[it];}
ymax=0;
for (Int_t it=0;it<ntime;it++){if(y[it]>ymax){ymax=y[it];timeref[i]=x[it];}}
interpolation[i] = new ROOT::Math::Interpolator(ntime, ROOT::Math::Interpolation::kCSPLINE);
interpolation[i]->SetData(ntime, x, y);
finter[i]= new TF1(Form("finter_%d",i),func,-200,ntime+200,nbparameters);
finter[i]->FixParameter(0,i); //numero de bloc
//for(Int_t p=0;p<maxwfpulses;p++){cout<<"p ="<<p<<endl;finter[i]->FixParameter(2+2*p,0.);finter[i]->FixParameter(3+2*p,0.);}
finter[i]->SetLineColor(4);
finter[i]->SetNpx(1100);
}
filewf.close();
//cout<<i<<" "<<timeref[i]<<endl;
 }

Double_t corr_time_HMS;
//Read tdc_offset_param (needed to determine HMS corrections to the timing)//For now, it is just one file
Float_t tdcoffset[nblocks];
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6171-6183/fit_e_runs/fit_elastic_runs/results_elastics/refwf/tdc_offset_param.txt");//Done
ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/tdc_offset_param.txt");//Done
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/tdc_offset_param.txt");//Done
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5183-5208/fit_e_runs/RWF/tdc_offset_param.txt");//Done
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/RWF/tdc_offset_param.txt");//Done
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2900-2920/RWF/tdc_offset_param.txt");//Done
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2875-2885/RWF/tdc_offset_param.txt");//Done
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2855-2871/RWF/tdc_offset_param.txt");//Done 
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/tdc_offset_param.txt");//Done
//ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1423-1511/RWF/tdc_offset_param.txt");//Done

 for(Int_t i=0;i<nblocks;i++){filetdc>>tdcoffset[i];}

tree->SetBranchStatus("*", 0); // Deactivate all branches

// Activate only the needed branches
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampWaveform", 1);
tree->SetBranchStatus("NPS.cal.fly.adcSampWaveform", 1);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcCounter", 1);
tree->SetBranchStatus("NPS.cal.fly.adcCounter", 1);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseAmp", 1);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseAmp", 1);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseInt", 1);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseInt", 1);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPed", 1);
tree->SetBranchStatus("NPS.cal.fly.adcSampPed", 1);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTime", 1);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseTime", 1);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTimeRaw", 1);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseTimeRaw", 1);
tree->SetBranchStatus("T.hms.hT1_tdcTime", 1);
tree->SetBranchStatus("T.hms.hT2_tdcTime", 1);
tree->SetBranchStatus("T.hms.hT3_tdcTime", 1);
tree->SetBranchStatus("T.hms.hTRIG1_tdcTime", 1);
tree->SetBranchStatus("T.hms.hTRIG2_tdcTime", 1);
tree->SetBranchStatus("T.hms.hTRIG3_tdcTime", 1);
tree->SetBranchStatus("T.hms.hTRIG4_tdcTime", 1);
tree->SetBranchStatus("T.hms.hTRIG5_tdcTime", 1);
tree->SetBranchStatus("T.hms.hTRIG6_tdcTime", 1);
tree->SetBranchStatus("beta", 1);
tree->SetBranchStatus("cernpeSum", 1);
tree->SetBranchStatus("caletracknorm", 1);
tree->SetBranchStatus("caletottracknorm", 1);
tree->SetBranchStatus("caletotnorm", 1);
tree->SetBranchStatus("H.react.x", 1);
tree->SetBranchStatus("H.react.y", 1);
tree->SetBranchStatus("H.react.z", 1);
tree->SetBranchStatus("H.gtr.dp", 1);
tree->SetBranchStatus("H.gtr.th", 1);
tree->SetBranchStatus("H.gtr.ph", 1);
tree->SetBranchStatus("H.gtr.px", 1);
tree->SetBranchStatus("H.gtr.py", 1);
tree->SetBranchStatus("H.gtr.pz", 1);
tree->SetBranchStatus("H.cer.npeSum", 1);
tree->SetBranchStatus("H.cal.etottracknorm", 1);


// Input variables

Int_t NSampWaveForm; tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampWaveform",&NSampWaveForm);//should correspond to Ndata
Double_t SampWaveForm[Ndata]; tree->SetBranchAddress("NPS.cal.fly.adcSampWaveform",&SampWaveForm); // Contain the waveforms for the tested calo blocks

Int_t  NadcCounter;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcCounter",&NadcCounter) ;// total number of pulses above theshold (could have up to 4 pulses per block)
Double_t adcCounter[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcCounter",&adcCounter) ;

Int_t NadcSampPulseAmp;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseAmp",&NadcSampPulseAmp) ;
Double_t adcSampPulseAmp[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseAmp",&adcSampPulseAmp) ;

Int_t NadcSampPulseInt;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseInt",&NadcSampPulseInt) ;
Double_t adcSampPulseInt[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseInt",&adcSampPulseInt) ;

Int_t NadcSampPulsePed;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPed",&NadcSampPulsePed) ;
Double_t adcSampPulsePed[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPed",&adcSampPulsePed) ;

Int_t NadcSampPulseTime;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseTime",&NadcSampPulseTime) ;
Double_t adcSampPulseTime[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseTime",&adcSampPulseTime) ;

Int_t NadcSampPulseTimeRaw;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseTimeRaw",&NadcSampPulseTimeRaw) ;
Double_t adcSampPulseTimeRaw[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseTimeRaw",&adcSampPulseTimeRaw) ;

Double_t hT1_tdcTime;tree->SetBranchAddress("T.hms.hT1_tdcTime",&hT1_tdcTime) ;
Double_t hT2_tdcTime;tree->SetBranchAddress("T.hms.hT2_tdcTime",&hT2_tdcTime) ;
Double_t hT3_tdcTime;tree->SetBranchAddress("T.hms.hT3_tdcTime",&hT3_tdcTime) ;
Double_t hTRIG1_tdcTime;tree->SetBranchAddress("T.hms.hTRIG1_tdcTime",&hTRIG1_tdcTime) ;
Double_t hTRIG2_tdcTime;tree->SetBranchAddress("T.hms.hTRIG2_tdcTime",&hTRIG2_tdcTime) ;
Double_t hTRIG3_tdcTime;tree->SetBranchAddress("T.hms.hTRIG3_tdcTime",&hTRIG3_tdcTime) ;
Double_t hTRIG4_tdcTime;tree->SetBranchAddress("T.hms.hTRIG4_tdcTime",&hTRIG4_tdcTime) ;
Double_t hTRIG5_tdcTime;tree->SetBranchAddress("T.hms.hTRIG5_tdcTime",&hTRIG5_tdcTime) ;
Double_t hTRIG6_tdcTime;tree->SetBranchAddress("T.hms.hTRIG6_tdcTime",&hTRIG6_tdcTime) ;
Double_t beta;tree->SetBranchAddress("beta",&beta) ;
Double_t cernpeSum;tree->SetBranchAddress("cernpeSum",&cernpeSum) ;
Double_t caletracknorm;tree->SetBranchAddress("caletracknorm",&caletracknorm) ;
Double_t caletottracknorm;tree->SetBranchAddress("caletottracknorm",&caletottracknorm) ;
Double_t caletotnorm;tree->SetBranchAddress("caletotnorm",&caletotnorm) ;
Double_t vx;tree->SetBranchAddress("H.react.x",&vx) ;
Double_t vy;tree->SetBranchAddress("H.react.y",&vy) ;
Double_t vz;tree->SetBranchAddress("H.react.z",&vz) ;
Double_t dp;tree->SetBranchAddress("H.gtr.dp",&dp) ;
Double_t th;tree->SetBranchAddress("H.gtr.th",&th) ;
Double_t ph;tree->SetBranchAddress("H.gtr.ph",&ph) ;
Double_t px;tree->SetBranchAddress("H.gtr.px",&px) ;
Double_t py;tree->SetBranchAddress("H.gtr.py",&py) ;
Double_t pz;tree->SetBranchAddress("H.gtr.pz",&pz) ;
Double_t cernpe;tree->SetBranchAddress("H.cer.npeSum",&cernpe) ;
Double_t caltracknorm;tree->SetBranchAddress("H.cal.etottracknorm",&caltracknorm) ;
Double_t hcnps_intlk_cz_t_back[ntemp];
Double_t hcnps_intlk_cz_t_front[ntemp];



// Enable only the required branches
for (int i = 0; i < ntemp; ++i) {
    tree->SetBranchStatus(Form("hcnps_intlk_cz_t_back_%i", i+1), 1);
    tree->SetBranchStatus(Form("hcnps_intlk_cz_t_front_%i", i+1), 1);

    tree->SetBranchAddress(Form("hcnps_intlk_cz_t_back_%i", i+1), &hcnps_intlk_cz_t_back[i]);
    tree->SetBranchAddress(Form("hcnps_intlk_cz_t_front_%i", i+1), &hcnps_intlk_cz_t_front[i]);
}

 cout<<"Nb of entries in rootfiles = "<<nentries<<endl;

// Output rootfile

Int_t tracage=0; //1 to see the waveforms event by event

TFile *fout;

if (tracage == 0) {

  TString rootfilePath = Form("nps_production_%d_%d.root", run, seg);
  
    cout << "Rootfile location: " << rootfilePath << endl;
    fout = new TFile(rootfilePath, "recreate");
}

TTree *treeout=new TTree("T","Tree organized");
treeout->SetAutoFlush(-300000000);
 
// Output variables

Double_t signal[nblocks][ntime];//treeout->Branch("signal",&signal,Form("signal[%d][%d]/D",nblocks,ntime));//if we want to save the raw waveforms in the output rootfile
Double_t ampl2[nblocks];//treeout->Branch("ampl2",&ampl2,Form("ampl2[%d]/D",nblocks));//amplitude of the pulse relatively to background
Double_t ampl[nblocks];treeout->Branch("ampl",&ampl,Form("ampl[%d]/D",nblocks));//amplitude of the pulse (a comparer avec Sampampl)
Double_t amplwf[nblocks];treeout->Branch("amplwf",&amplwf,Form("amplwf[%d]/D",nblocks));//amplwfitude of the pulse (a comparer avec Sampamplwf)
//treeout->Branch("wfampl",&wfampl,Form("wfampl[%d][%d]/D",nblocks,maxwfpulses));
//treeout->Branch("wftime",&wftime,Form("wftime[%d][%d]/D",nblocks,maxwfpulses));
Double_t Npulse[nblocks]; 
treeout->Branch("wfnpulse",&wfnpulse,Form("wfnpulse[%d]/I",nblocks));//number of pulses in one block 
Double_t Sampampl[nblocks];treeout->Branch("Sampampl",&Sampampl,Form("Sampampl[%d]/D",nblocks));//amplitude of the 1st pulse 
Double_t Samptime[nblocks];treeout->Branch("Samptime",&Samptime,Form("Samptime[%d]/D",nblocks));//time of the 1st pulse 'ns)
Double_t Sampped[nblocks];
Double_t Sampener[nblocks];
Double_t ener[nblocks];
Double_t time[nblocks];
Double_t timewf[nblocks];treeout->Branch("timewf",&timewf,Form("timewf[%d]/D",nblocks));//timewf position of the pulse maximum
Double_t chi2[nblocks];treeout->Branch("chi2",&chi2,Form("chi2[%d]/D",nblocks));//time of the 1st pulse 'ns)
Double_t bkg[nblocks];
Double_t enertot;treeout->Branch("enertot",&enertot,"enertot/D");// sum of all ener[i] : total energy deposited in the calo (!! energy calibration is not done yet)
Double_t integ[nblocks];
Double_t integtot;treeout->Branch("integtot",&integtot,"integtot/D");// sum of all integ[i]
Double_t noise[nblocks];
Double_t larg50[nblocks];//treeout->Branch("larg50",&larg50,Form("larg50[%d]/D",nblocks));//RMS of waveforms relatively to zero
Double_t larg90[nblocks];//treeout->Branch("larg90",&larg90,Form("larg90[%d]/D",nblocks));//RMS of waveforms relatively to zero
Int_t pres[nblocks];treeout->Branch("pres",&pres,Form("pres[%d]/I",nblocks));//RMS of waveforms relatively to zero
Int_t nrun;treeout->Branch("nrun",&nrun,"nrun/I");// run number
treeout->Branch("seg",&seg,"seg/I"); 
//treeout->Branch("job_number",&job_number,"job_number/I"); 
Int_t event;treeout->Branch("event",&event,"event/I");// event number in input rootfiles
treeout->Branch("beta",&beta,"beta/D") ;
treeout->Branch("cernpeSum",&cernpeSum,"cernpeSum/D") ;
treeout->Branch("caletracknorm",&caletracknorm,"caletracknorm/D") ;
treeout->Branch("caletottracknorm",&caletottracknorm,"caletottracknorm/D") ;
treeout->Branch("caletotnorm",&caletotnorm,"caletotnorm/D") ;
treeout->Branch("corr_time_HMS",&corr_time_HMS,"corr_time_HMS/D");
treeout->Branch("dp",&dp,"dp/D");
treeout->Branch("th",&th,"th/D");
treeout->Branch("ph",&ph,"ph/D");
treeout->Branch("px",&px,"px/D");
treeout->Branch("py",&py,"py/D");
treeout->Branch("pz",&pz,"pz/D");
treeout->Branch("vx",&vx,"vx/D");
treeout->Branch("vy",&vy,"vy/D");
treeout->Branch("vz",&vz,"vz/D");
 
Double_t temp_back[ntemp];
Double_t temp_front[ntemp];

treeout->Branch("temp_back",&temp_back, Form("temp_back[%i]/D", ntemp));
treeout->Branch("temp_front",&temp_front, Form("temp_front[%i]/D", ntemp));

//Other variables
Double_t sigmax[nblocks];
Int_t bloc,nsamp,ns,ilin,icol,ilinc,icolc,inp;
Double_t max50,max90,min50,min90;
Int_t nsampwf=0;
Int_t ndataprob=0;
Int_t nb=0;

//Histograms
for(Int_t i=0;i<nblocks;i++){
hsig_i[i] = new TH1F(Form("hsig_i%d",i),Form("hsig_i%d",i),ntime,0,ntime); 
hsig_i[i]->SetLineColor(1);
hsig_i[i]->SetLineWidth(2);
hsig_i[i]->GetXaxis()->SetTitle("Time (4 ns)");
hsig_i[i]->GetYaxis()->SetTitle("(mV)");
hsig_i[i]->GetYaxis()->SetLabelSize(0.05);
}


TLatex* tex = new TLatex();
tex->SetTextSize(0.015);
TH1::AddDirectory(kFALSE);

//Read the mean time positions of cosmic pulses (this is determined by the macro analyse_wassim.C)
Float_t timemean[nblocks],timemean2[nblocks];
for(Int_t ii=0;ii<nblocks;ii++){timemean[ii]=50;timemean2[ii]=150;} //1st pass analysis, when the macro analyse_wassim.C is not executed yet (30 is the default value)
//2nd pass analysis, when the macro analyse_wassim.C is already executed. If not, comment the following lines
ifstream filetime("filetime.txt");
Float_t dum;

TH1F *h1time = new TH1F("h1time","pulse (>20mV) shift (4*ns units) relatively to elastic refwf (all found pulses included)",200,-50,50);
TH1F *h2time = new TH1F("h2time","pulse (>20mV) time (ns) (all found pulses included)",200,-100,100);
 
/////// ANALYZE //////////////////////////////////////////////////////////

//start_event = job_number * evts_per_job;
//end_event = (job_number+1)* evts_per_job;
//if(end_event > nentries) end_event = nentries;

	  
   for(Int_t evt=0; evt<1000; evt++){
     //if (evt % 100 == 0) {
    tree->GetEntry(evt);
   

    event=evt;
    nrun=run;

    //  bool skip_this_event = false;  // Flag to indicate skipping the event where nsamp different than ntime
    
    // for(Int_t i=nb_of_runs-1;i>-1;i--){if(evt<entries[i])nrun=runnumber[i];}

    if (evt%100==0) {
      cout << " Entry = " << evt <<"/"<<nentries<<"  run="<<nrun<<"  cpu time="<<t.RealTime()<<endl;
      t.Continue();
    }

if(TMath::Abs(th)<0.08&&TMath::Abs(ph)<0.04&&TMath::Abs(dp)<10){
      
    if(NSampWaveForm>Ndata){
      nsampwf++;
      cout<<"!!!! NSampWaveForm problem  "<<evt<<"  "<<NSampWaveForm<<" "<<Ndata<<endl;
    }

    if(NSampWaveForm<=Ndata){ //NSampWaveForm must be <= Ndata (otherwise correct Ndata value)

      for(Int_t i=0;i<nblocks;i++){

	pres[i]=0;

	for(Int_t j=0;j<ntime;j++){
	  
	  signal[i][j]=0;}
	
      }//liste de presence des blocs! pres=0 si bloc absent, pres=1 s'il est present
 
//Extract the data from the complex variable SampWaveForm[] (NPS.cal.fly.adcSampWaveform)
    
    ns=0; //ns represent for a given event the element number of the NPS.cal.fly.adcSampWaveform variable 
//    nb=0;
    while (ns < NSampWaveForm) {
      
    bloc=SampWaveForm[ns];ns++; //bloc number (actually the slot number)
    nsamp=SampWaveForm[ns];ns++; //time samples (should be equal to ntime (100))
    
if(bloc==2000)bloc=1080; //modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators
if(bloc==2001)bloc=1081; //modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators
 
if(bloc<0||bloc>nslots-0.5){
  cout<<"slot number problem "<<evt<<" "<<bloc<<endl;ns=NSampWaveForm+1;
 }//to exit the while()


if(bloc>-0.5&&bloc<nslots){// that's what we expect!
  
pres[bloc]=1;
 
for (Int_t it=0;it<nsamp;it++){
  
  if (nsamp==ntime){
    
if(bloc>-0.5&&bloc<nblocks){
  signal[bloc][it]=SampWaveForm[ns];hsig_i[bloc]->SetBinContent(it+1,signal[bloc][it]);
 }
 
ns++;

}}
 
}//fin if(bloc number is good)
 
}//fin while()



    

//Calculate the output variables of the rootfile


    
if(pres[454]==0)nb++;
 
//initialisation
 
enertot=0;
integtot=0;
 
for(Int_t i=0;i<nblocks;i++){
  timewf[i]=-100;
  amplwf[i]=-100;
  chi2[i]=-1;
  ener[i]=0;
  integ[i]=0;
  noise[i]=0;
  bkg[i]=0;
  sigmax[i]=-100;
  Sampampl[i]=-100;
  Sampped[i]=-100;
  Samptime[i]=-100;
  Sampener[i]=-100;
  Npulse[i]=0;
 }
 
for(Int_t i=0;i<nblocks;i++){
  wfnpulse[i]=-100;

  for(Int_t p=0;p<maxwfpulses;p++){
    wftime[i][p]=-100;
    wfampl[i][p]=-100;
  }}
 
for(Int_t i=0;i<nblocks;i++){
  hsig_i[i]->SetLineColor(1);
 }


 
//Read the hcana calculated variables

//if(!(NadcCounter==NadcSampPulseAmp&&NadcSampPulseInt==NadcSampPulsePed&&NadcSampPulseInt==NadcSampPulseTime)){cout<<"!!!!! Problem Ndata !!!!!! "<<evt<<endl;ndataprob++;}
 
corr_time_HMS=0;
 
for(Int_t iNdata = 0; iNdata < NadcCounter; iNdata++){
  
if(adcCounter[iNdata] ==2000) adcCounter[iNdata]=1080;// nouveau numero du scintillateur attribue par Malek
if(adcCounter[iNdata] ==2001) adcCounter[iNdata]=1081;// nouveau numero du scintillateur attribue par Malek


//determine HMS time correction

 
 if(iNdata==0){
   corr_time_HMS=adcSampPulseTime[iNdata]-(adcSampPulseTimeRaw[iNdata]/16.)-tdcoffset[(int)(adcCounter[iNdata])];
 }
 
if(iNdata!=0){
  
  if(TMath::Abs(corr_time_HMS-(adcSampPulseTime[iNdata]-adcSampPulseTimeRaw[iNdata]/16.-tdcoffset[(int)(adcCounter[iNdata])]))>0.001){
    cout<<"problem HMS time correction event "<<evt<<endl;
  }}

 if(!(adcCounter[iNdata] >= 0 && adcCounter[iNdata] <nblocks+2)/*&&adcCounter[iNdata]!=196*/){
   cout<<"****** Problem adcCounter ******* "<<evt<<" "<<iNdata<<" "<<adcCounter[iNdata]<<endl;
 }
 
if(adcCounter[iNdata] >= 0 && adcCounter[iNdata] <nblocks){
  
Npulse[(int)(adcCounter[iNdata])]+=1;
hsig_i[(int)(adcCounter[iNdata])]->SetLineColor(2);
 
if(Npulse[(int)(adcCounter[iNdata])]==1){
  
Sampampl[(int)(adcCounter[iNdata])]=adcSampPulseAmp[iNdata];
Samptime[(int)(adcCounter[iNdata])]=adcSampPulseTime[iNdata];
Sampener[(int)(adcCounter[iNdata])]=adcSampPulseInt[iNdata];
Sampped[(int)(adcCounter[iNdata])]=adcSampPulsePed[iNdata];
}
 
if(Npulse[(int)(adcCounter[iNdata])]>1){
  
if(TMath::Abs(Samptime[(int)(adcCounter[iNdata])]-timemean2[(int)(adcCounter[iNdata])])>TMath::Abs(adcSampPulseTime[iNdata]-timemean2[(int)(adcCounter[iNdata])])){//si le 2eme pulse est plus proche du temps de reference que le 1er pulse alors on prend le 2eme
  
Sampampl[(int)(adcCounter[iNdata])]=adcSampPulseAmp[iNdata];
Samptime[(int)(adcCounter[iNdata])]=adcSampPulseTime[iNdata];
Sampener[(int)(adcCounter[iNdata])]=adcSampPulseInt[iNdata];
Sampped[(int)(adcCounter[iNdata])]=adcSampPulsePed[iNdata];
 
}}
}}

// Fit de la wf

for(Int_t i=0;i<nblocks;i++){
  
icol=i%ncol;
ilin=(int)(i/ncol);
 
if(pres[i]==1){ //For this setting?

for(Int_t it=0;it<nsamp;it++){

     hsig_i[i]->SetBinError(it+1,TMath::Sqrt(TMath::Abs(hsig_i[i]->GetBinContent(it+1)*4.096))/4.096);

     if(hsig_i[i]->GetBinContent(it+1)<1.){
       
     hsig_i[i]->SetBinError(it+1,TMath::Sqrt(TMath::Abs(1.*4.096))/4.096);}
     
 }//end loop over it

 if (nsamp==ntime){
Fitwf(i);//We call the fit here
 }

  // **Add this check to skip the block if finter[i] is invalid**
  if (!finter[i]) {
       cout << "Skipping block " << i << " due to missing finter[" << i << "] at event " << evt << endl;
       continue;
        }

 
chi2[i]=finter[i]->GetChisquare()/finter[i]->GetNDF();

for(Int_t p=0;p<TMath::Max(wfnpulse[i],1);p++){

wftime[i][p]=finter[i]->GetParameter(2+2*p)*dt+corr_time_HMS;//temps du pulse en ns

wfampl[i][p]=finter[i]->GetParameter(3+2*p); //amplitude du pulse en ns

if(wfampl[i][p]>20){

h2time->Fill(wftime[i][p],1.);

h1time->Fill(finter[i]->GetParameter(2+2*p),1.);
 
}// Fill the time spectrums

if(p==0){

timewf[i]=wftime[i][p];
amplwf[i]=wfampl[i][p];
 
 }

 
//if(p>0){if(TMath::Abs(wftime[i][p]-timerefacc/dt)<TMath::Abs(timewf[i]-timerefacc/dt)){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}}//pour prendre le pulse dont le temps est le plus proche de timerefacc

//New modification based on the production runs
 
if(p>0){
  
if(TMath::Abs(wftime[i][p]-timerefacc2*dt)<TMath::Abs(timewf[i]-timerefacc2*dt)){
         timewf[i]=wftime[i][p];
	 amplwf[i]=wfampl[i][p];
          }
          }//pour prendre le pulse dont le temps est le plus proche de timerefacc


 }//end of loop over maxpulses
 }// end of condition over col,lin
 }//end of loop over les blocs

 
//Calculation

for(Int_t i=0;i<nblocks;i++){

binmin=(int)(timemean[i]-25);
binmax=(int)(timemean[i]+59);       

for(Int_t it=0;it<nsamp;it++){
  
integ[i]+=signal[i][it];
integtot+=signal[i][it];

 if(it>binmin&&it<binmax){//cosmic pulse window
   
ener[i]+=signal[i][it];
enertot+=signal[i][it];
}
 
if(!(it>binmin&&it<binmax)){
  bkg[i]+=signal[i][it];}//background window
 
if(signal[i][it]>sigmax[i]){
  time[i]=it;
  sigmax[i]=signal[i][it];
  ampl[i]=signal[i][it];
}//determine the pulse maximum
}

ener[i]-=bkg[i]*(binmax-binmin-1)/(nsamp-(binmax-binmin-1));//subtract the bkg contribution (we have to normalize since the window widths are not the same)

bkg[i]=bkg[i]/(nsamp-(binmax-binmin-1));//mean value of bkg

for(Int_t it=0;it<nsamp;it++){
  
if(!(it>binmin&&it<binmax)){
  noise[i]+=(signal[i][it]-bkg[i])*(signal[i][it]-bkg[i])/(nsamp-(binmax-binmin-1));
}// RMS of the bkg
}
noise[i]=TMath::Sqrt(noise[i]);
}


 
// Calculation of signal widths
for(Int_t i=0;i<nblocks;i++){

ampl2[i]=ampl[i]-bkg[i];//amplitude of the pulse relatively to bkg

max50=0;
max90=50;
min50=100;
min90=100;
  
for(Int_t it=time[i];it<nsamp;it++){
  //aller vers la droite du maximum
  if((signal[i][it]-bkg[i])>=ampl2[i]*0.5){max50=it;}
  if((signal[i][it]-bkg[i])>=ampl2[i]*0.1){max90=it;}
}
 
for(Int_t it=time[i];it>-0.5;it--){//aller vers la gauche du maximum

if((signal[i][it]-bkg[i])>=ampl2[i]*0.5){min50=it;}
if((signal[i][it]-bkg[i])>=ampl2[i]*0.1){min90=it;}
}
 
larg50[i]=max50-min50;
larg90[i]=max90-min90;

}	      

if(evt==2){cout<<"SampAmp="<<Sampampl[7]<<"  ampl="<<ampl[7]<<endl;}

    treeout->Fill();

		
    //std::cout << "we just filled the tree" << std::endl;
    //gObjectTable->Print();




    
}// end if(NSampWaveForm<=Ndata)
    
}//end of the cut over HMS

     // After processing events, check the list of objects in memory

 if(evt%1000==0) { std::cout << "Objects in memory After event processing:" << std::endl;
   gObjectTable->Print();}


} //end for(evt)
	
/////////////////////////////////////////////////////////////////////////////////////////
 
//Removed the chunk responsable for the event per event check for now

/////// Write the output files //////////////////////////////////////////////////////////

   for (int i = 0; i < nblocks; i++) {
    delete hsig_i[i];         // Delete histogram objects
    delete finter[i];         // Delete fit function objects
    delete interpolation[i];  // Delete interpolator objects
}


// Then safely write and delete at the end
fout->cd();
h1time->Write("h1time");
h2time->Write("h2time");
delete h1time;
delete h2time;
fout->Write();
fout->Close();

   
cout<<"fin de l'analyse"<<endl;
cout<<nsampwf<<" "<<ndataprob<<" "<<nb<<endl;

/*
  fout->cd();
  fout->Write();
  h1time->Write("h1time");
  h2time->Write("h2time");
  histo->Write("histo");
  delete histo;
  delete h1time;
  delete h2time;
  fout->Close();
*/


}
