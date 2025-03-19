 #include "TCanvas.h"
 #include "TROOT.h"
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
 #include <iostream>
 #include <fstream>
const int ntime=110; // number of time samples for each fADC channel (100 for runs 55 and 56)
const int nfADC=7; //number of fADC modules (16 channels each) used in this run
const int nchannel=16; //number of channels in each fADC module
const int ncol=30; //number of calo colomns tested in this run (3 for runs 55 and 56)
const int nlin=36; //number of calo blocks in each colomn
const int nblocks=ncol*nlin; //number of tested calo blocks
const int maxpulses=4; //nb maximal de pulses dans la fentre ntime 
const int nsignals=nblocks*maxpulses;
Double_t dt=4.; //time bin (sample) width (4 ns), the total time window is then ntime*dt
const int ntemp=56; //nb of temp sensors

ROOT::Math::Interpolator *interpolation[nblocks];
TH1F *hsig_i[nblocks];//histogram showing the waveform for block i and a given event
Int_t wfnpulse[nblocks],wfnpulsei[nblocks],wfnpulsesup[nblocks];
const int maxwfpulses=12; //nb maximal de pulses que la wfa peut trouver
Double_t wfampl[nblocks][maxwfpulses]; 
Double_t wftime[nblocks][maxwfpulses]; 
Int_t wfnpulse_sev[nblocks][maxwfpulses];
TF1 *finter[nblocks]; 
Double_t timeref[nblocks];
Double_t calo_dist=3;// calo distance to the target (m)
Double_t timerefacc=(calo_dist-9.5)/(3.e8*1.e-9*4);//time offset (4ns units) due to the change of calo_dist

Int_t cpulse;
Float_t corr_time_HMS;
Float_t cortime[nblocks];

// La fonction du fit de la wf
Double_t func(Double_t *x, Double_t *par){ 
 Int_t j= (int)(par[0]);
 Double_t val=0;
 for(Int_t p=0;p<maxwfpulses;p++){
 if(x[0]-par[2+2*p]>1&&x[0]-par[2+2*p]<109)val+=par[3+2*p]*interpolation[j]->Eval(x[0]-par[2+2*p]);
 }
 return val+par[1];
}


//////////////////////////////////////////////////////////////////////////////

void npulsewf(Int_t bn, Int_t minsamp, Int_t maxsamp){

//Detect the pulses

wfnpulse[bn]=0;
Int_t good=0;

for(Int_t p=0;p<maxwfpulses;p++){wfampl[bn][p]=-1;}
for (Int_t it=minsamp;it<maxsamp-6;it++){
if(hsig_i[bn]->GetBinContent(it+1)<hsig_i[bn]->GetBinContent(it+2)&&hsig_i[bn]->GetBinContent(it+2)<hsig_i[bn]->GetBinContent(it+3)&&hsig_i[bn]->GetBinContent(it+3)<=hsig_i[bn]->GetBinContent(it+4)&&hsig_i[bn]->GetBinContent(it+4)>=hsig_i[bn]->GetBinContent(it+5)){
if(hsig_i[bn]->GetBinContent(it+4)>0){
wfampl[bn][wfnpulse[bn]]=hsig_i[bn]->GetBinContent(it+4);
wftime[bn][wfnpulse[bn]]=hsig_i[bn]->GetBinCenter(it+4);
if(TMath::Abs(wftime[bn][wfnpulse[bn]]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt-timerefacc)<4.){good=1;cpulse=wfnpulse[bn];}
//cout<<"detection bloc "<<bn<<" pulse="<<wfnpulse[bn]<<" time="<<wftime[bn][wfnpulse[bn]]<<" (4ns) ampl="<<wfampl[bn][wfnpulse[bn]]<<" reftime= "<<timeref[bn]<<" (4ns) diff="<<TMath::Abs(wftime[bn][wfnpulse[bn]]-timeref[bn])<<" good="<<good<<endl;
wfnpulse[bn]++;
if(wfnpulse[bn]==maxwfpulses){wfnpulse[bn]=maxwfpulses-1;it=ntime;}
it+=4;
}}}
if(wfnpulse[bn]>maxwfpulses-2){cout<<"Warning : nb de pulses excessivement eleve dans la wf du bloc "<<bn<<endl;}
}


//////////////////////////////////////////////////////////////////////////////

void Fitwf(Int_t bn, Int_t minsamp, Int_t maxsamp){

//Detect the pulses

wfnpulse[bn]=0;
wfnpulsesup[bn]=0;
Int_t good=0;

for(Int_t p=0;p<maxwfpulses;p++){wfampl[bn][p]=-1;}
for (Int_t it=minsamp;it<maxsamp-6;it++){
if(hsig_i[bn]->GetBinContent(it+1)<hsig_i[bn]->GetBinContent(it+2)&&hsig_i[bn]->GetBinContent(it+2)<hsig_i[bn]->GetBinContent(it+3)&&hsig_i[bn]->GetBinContent(it+3)<=hsig_i[bn]->GetBinContent(it+4)&&hsig_i[bn]->GetBinContent(it+4)>=hsig_i[bn]->GetBinContent(it+5)){
if(hsig_i[bn]->GetBinContent(it+4)>0){
wfampl[bn][wfnpulse[bn]]=hsig_i[bn]->GetBinContent(it+4);
wftime[bn][wfnpulse[bn]]=hsig_i[bn]->GetBinCenter(it+4);

//check if the pulse can be found with a more severe condition
if(hsig_i[bn]->GetBinContent(it)<hsig_i[bn]->GetBinContent(it+1)&&hsig_i[bn]->GetBinContent(it+1)<hsig_i[bn]->GetBinContent(it+2)&&hsig_i[bn]->GetBinContent(it+2)<hsig_i[bn]->GetBinContent(it+3)&&hsig_i[bn]->GetBinContent(it+3)<=hsig_i[bn]->GetBinContent(it+4)&&hsig_i[bn]->GetBinContent(it+4)>=hsig_i[bn]->GetBinContent(it+5)&&hsig_i[bn]->GetBinContent(it+5)>=hsig_i[bn]->GetBinContent(it+6)){wfnpulse_sev[bn][wfnpulse[bn]]=1;}

if(TMath::Abs(wftime[bn][wfnpulse[bn]]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt-timerefacc)<4.){good=1;cpulse=wfnpulse[bn];}
wfnpulse[bn]++;
if(wfnpulse[bn]==maxwfpulses){wfnpulse[bn]=maxwfpulses-1;it=ntime;}
it+=4;
}}}
if(wfnpulse[bn]>maxwfpulses-2){cout<<"Warning : nb de pulses excessivement eleve dans la wf du bloc "<<bn<<endl;}

//Adjust the parameters of the fit function
for(Int_t p=0;p<maxwfpulses;p++){finter[bn]->FixParameter(2+2*p,0.);finter[bn]->FixParameter(3+2*p,0.);}
if(wfnpulse[bn]>0&&good==1){
for(Int_t p=0;p<TMath::Min(maxwfpulses,wfnpulse[bn]);p++){
finter[bn]->ReleaseParameter(2+2*p);finter[bn]->ReleaseParameter(3+2*p);
finter[bn]->SetParameter(2+2*p,wftime[bn][p]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt);finter[bn]->SetParameter(3+2*p,wfampl[bn][p]);
finter[bn]->SetParLimits(2+2*p,wftime[bn][p]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt-3,wftime[bn][p]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt+3);finter[bn]->SetParLimits(3+2*p,wfampl[bn][p]*0.2,wfampl[bn][p]*3);
}}
if(wfnpulse[bn]>0&&good==0){
for(Int_t p=0;p<TMath::Min(maxwfpulses,wfnpulse[bn]);p++){
finter[bn]->ReleaseParameter(2+2*p);finter[bn]->ReleaseParameter(3+2*p);
finter[bn]->SetParameter(2+2*p,wftime[bn][p]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt);finter[bn]->SetParameter(3+2*p,wfampl[bn][p]);
finter[bn]->SetParLimits(2+2*p,wftime[bn][p]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt-3,wftime[bn][p]-timeref[bn]+(corr_time_HMS-cortime[bn])/dt+3);finter[bn]->SetParLimits(3+2*p,wfampl[bn][p]*0.2,wfampl[bn][p]*3);
}

//Search for an eventual pulse
cpulse=wfnpulse[bn];
finter[bn]->ReleaseParameter(2+2*wfnpulse[bn]);finter[bn]->ReleaseParameter(3+2*wfnpulse[bn]);
finter[bn]->SetParameter(2+2*wfnpulse[bn],timerefacc+(corr_time_HMS-cortime[bn])/dt);finter[bn]->SetParameter(3+2*wfnpulse[bn],2);
finter[bn]->SetParLimits(2+2*wfnpulse[bn],timerefacc+(corr_time_HMS-cortime[bn])/dt-4,timerefacc+(corr_time_HMS-cortime[bn])/dt+4);finter[bn]->SetParLimits(3+2*wfnpulse[bn],0.05,10);
wfnpulse[bn]++;wfnpulsesup[bn]++;
}
if(wfnpulse[bn]==0){
cpulse=wfnpulse[bn];
for(Int_t p=0;p<1;p++){
finter[bn]->ReleaseParameter(2+2*p);finter[bn]->ReleaseParameter(3+2*p);
finter[bn]->SetParameter(2+2*p,timerefacc+(corr_time_HMS-cortime[bn])/dt);finter[bn]->SetParameter(3+2*p,2);
finter[bn]->SetParLimits(2+2*p,timerefacc+(corr_time_HMS-cortime[bn])/dt-4,timerefacc+(corr_time_HMS-cortime[bn])/dt+4);finter[bn]->SetParLimits(3+2*p,0.05,10);
}
wfnpulse[bn]++;wfnpulsesup[bn]++;
}
finter[bn]->SetParameter(1,0.);finter[bn]->SetParLimits(1,-100,100.);
hsig_i[bn]->Fit(Form("finter_%d",bn),"Q","",minsamp,maxsamp);
}


//////////////////////////////////////////////////////////////////////////////////////////


void TEST_2(Int_t runnumber, Double_t calodist, int seg){

// Constants (valid for cosmic run 901 taken in Sep 2023)
const int nslots=1104; //nb maximal de slots dans tous les fADC
const int Ndata=nslots*(ntime+1+1); //(1104 slots fADC au total mais pas tous utilises y compris 2 PM scintillateurs?) should correspond to Ndata.NPS.cal.fly.adcSampWaveform variable in the root file
Double_t ADCtomV=1000./4096; //1000 mV correspond a 4096 canaux
Double_t integtopC= (1./4096)*(dt*1.e3/50); // (1 V / 4096 adc channels) * (4000 ps time sample / 50 ohms input resistance) = 0.020 pc/channel 

TStopwatch t;
t.Start();
calo_dist=calodist;

// where the elastic pulse is expected to be
Int_t binmin=30;
Int_t binmax=109;

// Input root file
TChain *tree=new TChain("T");
//tree->Add(Form("rootfiles/nps_hms_coin_%d_0_1_-1.root",runnumber));
TString filename = Form("/cache/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", runnumber, seg);//new files with the HMS calibration
//TString filename = "../nps_hms_coin_1993_1000.root";// test rootfile to Malek
//Check if the rootfile is not found and check if we have a sensible number of events in the segment
 
 TFile *file = TFile::Open(filename);
 
 Int_t nentries=tree->GetEntries();

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
Int_t preswf[nblocks];
for(Int_t i=0;i<nblocks;i++){
  
    if (runnumber>6183 && runnumber<7500) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6171-6183/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>6168 && runnumber<6171) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>5236 && runnumber<6151) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>5208 && runnumber<5217) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5183-5208/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>3898 && runnumber<5183) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>2920 && runnumber<3883) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2900-2920/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>2885 && runnumber<2900) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2875-2885/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>2871 && runnumber<2875) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2855-2871/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>1982 && runnumber<2855) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>1560 && runnumber<1961) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1423-1511/RWF/ref_wf_%d.txt",i));
    }
    
timeref[i]=-1.e6;
preswf[i]=0;
if(filewf.is_open()){
filewf>>timeref[i]>>dum1;
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
preswf[i]=1;
}
filewf.close();
//cout<<i<<" "<<timeref[i]<<endl;
}

//Check if all ref wf are present
Int_t probrefwf=0;
for(Int_t i=0;i<nblocks;i++){if(preswf[i]==0){probrefwf++;cout<<"REF WF for block "<<i<<" is missing !"<<endl;}}

//Read tdc_offset_param (needed to determine HMS corrections to the timing)
Float_t tdcoffset[nblocks];
ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/tdc_offset_param.txt");//MUST be the SAME file used by HCANA 
for(Int_t i=0;i<nblocks;i++){filetdc>>tdcoffset[i];}

//Lecture des corrections timing 
ifstream filetime;
Float_t dum;
Float_t rmstime[nblocks];
filetime.open("filetime_step_i.txt");
for(Int_t i=0;i<nblocks;i++){filetime>>dum>>cortime[i]>>dum>>rmstime[i]>>dum;}


// Input variables
// T->Print(); to see the type (int, double, ...) of each variable in the rootfile
tree->SetBranchStatus("*", false);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampWaveform",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampWaveform",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcCounter",true);
tree->SetBranchStatus("NPS.cal.fly.adcCounter",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseAmp",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseAmp",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseInt",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseInt",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPed",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPed",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTime",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseTime",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTimeRaw",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseTimeRaw",true);
tree->SetBranchStatus("T.hms.hT1_tdcTime",true);
tree->SetBranchStatus("T.hms.hT2_tdcTime",true);
tree->SetBranchStatus("T.hms.hT3_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG1_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG2_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG3_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG4_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG5_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG6_tdcTime",true);
tree->SetBranchStatus("H.react.x",true);
tree->SetBranchStatus("H.react.y",true);
tree->SetBranchStatus("H.react.z",true);
tree->SetBranchStatus("H.gtr.dp",true);
tree->SetBranchStatus("H.gtr.th",true);
tree->SetBranchStatus("H.gtr.ph",true);
tree->SetBranchStatus("H.gtr.px",true);
tree->SetBranchStatus("H.gtr.py",true);
tree->SetBranchStatus("H.gtr.pz",true);
tree->SetBranchStatus("H.cer.npeSum",true);
tree->SetBranchStatus("H.cal.etottracknorm",true);


/////
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
if(tracage==0) fout=new TFile(Form("rootfiles/nps_elaswf_%i.root",runnumber),"recreate");

TTree *treeout=new TTree("T","Tree organized");

// Output variables
Double_t signal[nblocks][ntime];//treeout->Branch("signal",&signal,Form("signal[%d][%d]/D",nblocks,ntime));//if we want to save the raw waveforms in the output rootfile
Double_t ampl2[nblocks];//treeout->Branch("ampl2",&ampl2,Form("ampl2[%d]/D",nblocks));//amplitude of the pulse relatively to background
Double_t ampl[nblocks];treeout->Branch("ampl",&ampl,Form("ampl[%d]/D",nblocks));//amplitude of the pulse (a comparer avec Sampampl)
Double_t amplwf[nblocks];treeout->Branch("amplwf",&amplwf,Form("amplwf[%d]/D",nblocks));//amplwfitude of the pulse (a comparer avec Sampamplwf)
Double_t baselinewf[nblocks];treeout->Branch("baselinewf",&baselinewf,Form("baselinewf[%d]/D",nblocks));//baselinewfitude of the pulse (a comparer avec Sampbaselinewf)
treeout->Branch("wfampl",&wfampl,Form("wfampl[%d][%d]/D",nblocks,maxwfpulses));
treeout->Branch("wftime",&wftime,Form("wftime[%d][%d]/D",nblocks,maxwfpulses));
Double_t Npulse[nblocks];treeout->Branch("Npulse",&Npulse,Form("Npulse[%d]/D",nblocks));//number of pulses in one block 
treeout->Branch("wfnpulse",&wfnpulse,Form("wfnpulse[%d]/I",nblocks));//number of pulses in one block 
treeout->Branch("wfnpulse_sev",&wfnpulse_sev,Form("wfnpulse_sev[%d][%d]/I",nblocks,maxwfpulses)); // 1 if the pulse can be detected with a more severe condition
treeout->Branch("wfnpulsesup",&wfnpulsesup,Form("wfnpulsesup[%d]/I",nblocks));//1 if a supl fit was performed in the pulse region, 0 otherwise 
treeout->Branch("wfnpulsei",&wfnpulsei,Form("wfnpulsei[%d]/I",nblocks));//number of initial pulses in one block found by the algo 
Int_t wfnpulse10[nblocks];treeout->Branch("wfnpulse10",&wfnpulse10,Form("wfnpulse10[%d]/I",nblocks));//number of pulses (>10mV) in one block 
Int_t wfnpulse1[nblocks];treeout->Branch("wfnpulse1",&wfnpulse1,Form("wfnpulse1[%d]/I",nblocks));//number of pulses (>1mV) in one block 
Double_t Sampampl[nblocks];treeout->Branch("Sampampl",&Sampampl,Form("Sampampl[%d]/D",nblocks));//amplitude of the 1st pulse 
Double_t Samptime[nblocks];treeout->Branch("Samptime",&Samptime,Form("Samptime[%d]/D",nblocks));//time of the 1st pulse 'ns)
Double_t Sampped[nblocks];treeout->Branch("Sampped",&Sampped,Form("Sampped[%d]/D",nblocks));//ped of the 1st pulse 'ns)
Double_t Sampener[nblocks];//treeout->Branch("Sampener",&Sampener,Form("Sampener[%d]/D",nblocks));//integral of the 1st pulse 
Double_t ener[nblocks];treeout->Branch("ener",&ener,Form("ener[%d]/D",nblocks));//integral of waveforms in the pulse window - integral of waveform in the backgroung window
Double_t time[nblocks];//treeout->Branch("time",&time,Form("time[%d]/D",nblocks));//time position of the pulse maximum
Double_t timewf[nblocks];treeout->Branch("timewf",&timewf,Form("timewf[%d]/D",nblocks));//timewf position of the pulse 
Double_t chi2[nblocks];treeout->Branch("chi2",&chi2,Form("chi2[%d]/D",nblocks));//time of the 1st pulse 'ns)
Double_t bkg[nblocks];treeout->Branch("bkg",&bkg,Form("bkg[%d]/D",nblocks));//integral of waveforms in the backgroung window
Double_t enertot;treeout->Branch("enertot",&enertot,"enertot/D");// sum of all ener[i] : total energy deposited in the calo (!! energy calibration is not done yet)
Double_t integ[nblocks];//treeout->Branch("integ",&integ,Form("integ[%d]/D",nblocks));//integral of waveforms in [0,dt*ntime ns] window
Double_t integtot;//treeout->Branch("integtot",&integtot,"integtot/D");// sum of all integ[i]
Double_t noise[nblocks];treeout->Branch("noise",&noise,Form("noise[%d]/D",nblocks));//RMS of waveforms relatively to zero
Double_t larg50[nblocks];//treeout->Branch("larg50",&larg50,Form("larg50[%d]/D",nblocks));//RMS of waveforms relatively to zero
Double_t larg90[nblocks];//treeout->Branch("larg90",&larg90,Form("larg90[%d]/D",nblocks));//RMS of waveforms relatively to zero
Int_t pres[nblocks];treeout->Branch("pres",&pres,Form("pres[%d]/I",nblocks));//1 if the waveform data is available
treeout->Branch("preswf",&preswf,Form("preswf[%d]/I",nblocks));//1 if the ref shapes is available
Int_t nrun;treeout->Branch("nrun",&nrun,"nrun/I");// run number
treeout->Branch("seg",&seg,"seg/I"); 
Int_t event;treeout->Branch("event",&event,"event/I");// event number in input rootfiles
treeout->Branch("dp",&dp,"dp/D");
treeout->Branch("th",&th,"th/D");
treeout->Branch("ph",&ph,"ph/D");
treeout->Branch("px",&px,"px/D");
treeout->Branch("py",&py,"py/D");
treeout->Branch("pz",&pz,"pz/D");
treeout->Branch("vx",&vx,"vx/D");
treeout->Branch("vy",&vy,"vy/D");
treeout->Branch("vz",&vz,"vz/D");
treeout->Branch("corr_time_HMS",&corr_time_HMS,"corr_time_HMS/D");
Double_t cputime;treeout->Branch("cputime",&cputime,"cputime/D");

//Other variables
Double_t sigmax[nblocks];
Int_t bloc,nsamp,ns,ilin,icol,ilinc,icolc,inp;
Double_t max50,max90,min50,min90;

//Problem variables
Int_t nsampwf=0;
Int_t ndataprob=0;
Int_t probnsamp=0;
Int_t probslot=0;
Int_t probHMStime=0;

//Histograms
for(Int_t i=0;i<nblocks;i++){
hsig_i[i] = new TH1F(Form("hsig_i%d",i),Form("hsig_i%d",i),ntime,0,ntime); 
hsig_i[i]->SetLineColor(1);
hsig_i[i]->SetLineWidth(2);
hsig_i[i]->GetXaxis()->SetTitle("Time (4 ns)");
hsig_i[i]->GetYaxis()->SetTitle("(mV)");
hsig_i[i]->GetYaxis()->SetLabelSize(0.05);
}


// Variables du tracage
TCanvas *c1 = new TCanvas("c1","c1",1300,1000);c1->Divide(ncol,nlin,0,0);
TCanvas *c11 = new TCanvas("c11","c11",1000,800);

TLatex* tex = new TLatex();
tex->SetTextSize(0.015);
TH1::AddDirectory(kFALSE);
TH1F* histo= new TH1F("histo","histo",100,0,500);

Float_t timemean2[nblocks];
for(Int_t ii=0;ii<nblocks;ii++){timemean2[ii]=170+timerefacc*dt;} //where the coincidence pulse is expected (ns units)

TH1F *h1time = new TH1F("h1time","h1time",200,-50,50);
TH1F *h2time = new TH1F("h2time","h2time",200,-100,100);


/////// ANALYZE //////////////////////////////////////////////////////////

for(Int_t evt=0; evt<100000; evt++){
  
    tree->GetEntry(evt);
    event=evt;
    nrun=runnumber;
    // std::cout<<" timerefacc  =  "<<timerefacc<<endl;
    if (evt%1000==0) {cout << " Entry = " << evt <<"/"<<nentries<<"  run="<<nrun<<"  cpu time="<<t.RealTime()<<endl;t.Continue();}
    
    if(TMath::Abs(th)<0.08&&TMath::Abs(ph)<0.04&&TMath::Abs(dp)<10&&evt%100==0){//still have the 1% here
       
    if(NSampWaveForm>Ndata){nsampwf++; cout<<"!!!! NSampWaveForm problem  "<<evt<<"  "<<NSampWaveForm<<" "<<Ndata<<endl;}

    if(NSampWaveForm<=Ndata){ //NSampWaveForm must be <= Ndata (otherwise correct Ndata value)
    for(Int_t i=0;i<nblocks;i++){pres[i]=0;for(Int_t j=0;j<ntime;j++)signal[i][j]=0;}//liste de presence des blocs! pres=0 si bloc absent, pres=1 s'il est present
 
//Extract the data from the complex variable SampWaveForm[] (NPS.cal.fly.adcSampWaveform)
    ns=0; //ns represent for a given event the element number of the NPS.cal.fly.adcSampWaveform variable 
    while (ns < NSampWaveForm) {
    bloc=SampWaveForm[ns];ns++; //bloc number (actually the slot number)
    nsamp=SampWaveForm[ns];ns++; //time samples (should be equal to ntime (100))
    histo->Fill(nsamp,1);
if(bloc==2000)bloc=1080; //modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators
if(bloc==2001)bloc=1081; //modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators
if(bloc<0||bloc>nslots-0.5){cout<<"slot number problem "<<evt<<" "<<bloc<<endl;ns=NSampWaveForm+1;}//to exit the while() 
if(nsamp!=ntime){cout<<"nsamp problem "<<evt<<" "<<nsamp<<" "<<bloc<<endl;probnsamp++;} //just a check!
if(bloc>-0.5&&bloc<nslots){// that's what we expect!
//cout<<"lecture"<<endl;
pres[bloc]=1;
for (Int_t it=0;it<nsamp;it++){
if(bloc>-0.5&&bloc<nblocks){signal[bloc][it]=SampWaveForm[ns];hsig_i[bloc]->SetBinContent(it+1,signal[bloc][it]);}
ns++;
}
}//fin if(bloc number is good)
    }//fin while()


//initialisation
enertot=0;integtot=0;
for(Int_t i=0;i<nblocks;i++){timewf[i]=-100;amplwf[i]=-100;;baselinewf[i]=-300;chi2[i]=-1;ener[i]=0;integ[i]=0;noise[i]=0;bkg[i]=0;sigmax[i]=-100;Sampampl[i]=-100;Sampped[i]=-100;Samptime[i]=-100;Sampener[i]=-100;Npulse[i]=0;}
for(Int_t i=0;i<nblocks;i++){wfnpulse[i]=-100;wfnpulsesup[i]=-100;wfnpulsei[i]=-100;wfnpulse10[i]=-100;wfnpulse1[i]=-100;for(Int_t p=0;p<maxwfpulses;p++){wftime[i][p]=-300;wfampl[i][p]=-100;wfnpulse_sev[i][p]=0;}}
for(Int_t i=0;i<nblocks;i++){hsig_i[i]->SetLineColor(1);}

//Read the hcana calculated variables
if(!(NadcCounter==NadcSampPulseAmp&&NadcSampPulseInt==NadcSampPulsePed&&NadcSampPulseInt==NadcSampPulseTime)){cout<<"!!!!! Problem Ndata !!!!!! "<<evt<<endl;ndataprob++;}
corr_time_HMS=0;
for(Int_t iNdata = 0; iNdata < NadcCounter; iNdata++){
if(adcCounter[iNdata] ==2000) adcCounter[iNdata]=1080;// nouveau numero du scintillateur attribue par Malek
if(adcCounter[iNdata] ==2001) adcCounter[iNdata]=1081;// nouveau numero du scintillateur attribue par Malek


//determine HMS time correction
if(iNdata==0)corr_time_HMS=adcSampPulseTime[iNdata]-adcSampPulseTimeRaw[iNdata]/16.-tdcoffset[(int)(adcCounter[iNdata])];
if(iNdata!=0){if(TMath::Abs(corr_time_HMS-(adcSampPulseTime[iNdata]-adcSampPulseTimeRaw[iNdata]/16.-tdcoffset[(int)(adcCounter[iNdata])]))>0.001){cout<<"problem HMS time correction event "<<evt<<endl;probHMStime++;}}

if(!(adcCounter[iNdata] >= 0 && adcCounter[iNdata] <nblocks+2)/*&&adcCounter[iNdata]!=196*/)cout<<"****** Problem adcCounter ******* "<<evt<<" "<<iNdata<<" "<<adcCounter[iNdata]<<endl;
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

// Fit of the wf
for(Int_t i=0;i<nblocks;i++){
if(pres[i]==1&&preswf[i]==1&&nsamp==ntime){
for(Int_t it=0;it<nsamp;it++){
hsig_i[i]->SetBinError(it+1,TMath::Sqrt(TMath::Abs(hsig_i[i]->GetBinContent(it+1)*4.096/2.))/4.096);
if(hsig_i[i]->GetBinContent(it+1)<1.)hsig_i[i]->SetBinError(it+1,TMath::Sqrt(TMath::Abs(1.*4.096/2.))/4.096);
}
cpulse=-100;
npulsewf(i,1,110);
wfnpulsei[i]=wfnpulse[i];
Fitwf(i,1,110);
wfnpulse1[i]=0;
for(Int_t p=0;p<TMath::Max(wfnpulse[i],1);p++){if(wfampl[i][p]>1&&p!=cpulse)wfnpulse1[i]++;}
if(tracage==1)cout<<i<<"  cpulse ="<<cpulse<<" wfnpulse1= "<<wfnpulse1[i]<<endl;
chi2[i]=finter[i]->GetChisquare()/finter[i]->GetNDF();
wfnpulse10[i]=0;
baselinewf[i]=finter[i]->GetParameter(1);
for(Int_t p=0;p<TMath::Max(wfnpulse[i],1);p++){
wftime[i][p]=finter[i]->GetParameter(2+2*p)*dt+corr_time_HMS-cortime[i]-timerefacc*dt;//temps du pulse en ns
wfampl[i][p]=finter[i]->GetParameter(3+2*p); //amplitude du pulse en mV
if(wfampl[i][p]>10)wfnpulse10[i]++;
if(wfampl[i][p]>20){h2time->Fill(wftime[i][p],1.);h1time->Fill(finter[i]->GetParameter(2+2*p)-timerefacc+corr_time_HMS/dt,1.);}
if(p==0){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}
if(p>0){if(TMath::Abs(wftime[i][p])<TMath::Abs(timewf[i])){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}}//pour prendre le pulse dont le temps est le plus proche de timerefacc
}}}

//Calculation
for(Int_t i=0;i<nblocks;i++){
for(Int_t it=0;it<nsamp;it++){ 
integ[i]+=signal[i][it];
integtot+=signal[i][it];
if(it>binmin&&it<binmax){//cosmic pulse window
//if(it>47&&it<84){//cosmic pulse window
ener[i]+=signal[i][it];
enertot+=signal[i][it];
}
if(!(it>binmin&&it<binmax)){bkg[i]+=signal[i][it];}//background window
if(signal[i][it]>sigmax[i]){time[i]=it;sigmax[i]=signal[i][it];ampl[i]=signal[i][it];}//determine the pulse maximum
}

ener[i]-=bkg[i]*(binmax-binmin-1)/(nsamp-(binmax-binmin-1));//subtract the bkg contribution (we have to normalize since the window widths are not the same)
bkg[i]=bkg[i]/(nsamp-(binmax-binmin-1));//mean value of bkg
for(Int_t it=0;it<nsamp;it++){ 
if(!(it>binmin&&it<binmax)){noise[i]+=(signal[i][it]-bkg[i])*(signal[i][it]-bkg[i])/(nsamp-(binmax-binmin-1));}// RMS of the bkg
}
noise[i]=TMath::Sqrt(noise[i]);
}
// Calculation of signal widths
for(Int_t i=0;i<nblocks;i++){
ampl2[i]=ampl[i]-bkg[i];//amplitude of the pulse relatively to bkg
max50=0;max90=50;min50=100;min90=100;
for(Int_t it=time[i];it<nsamp;it++){//aller vers la droite du maximum
if((signal[i][it]-bkg[i])>=ampl2[i]*0.5)max50=it;
if((signal[i][it]-bkg[i])>=ampl2[i]*0.1)max90=it;
}
for(Int_t it=time[i];it>-0.5;it--){//aller vers la gauche du maximum
if((signal[i][it]-bkg[i])>=ampl2[i]*0.5)min50=it;
if((signal[i][it]-bkg[i])>=ampl2[i]*0.1)min90=it;
}
larg50[i]=max50-min50;
larg90[i]=max90-min90;
}
    cputime=t.RealTime();
    treeout->Fill();

    }// end if(NSampWaveForm<=Ndata)
}

} //end for(evt)

/////////////////////////////////////////////////////////////////////////////////////////
Double_t h1timemax=h1time->GetBinCenter(h1time->FindFirstBinAbove(h1time->GetMaximum()*0.99));
Double_t h2timemax=h2time->GetBinCenter(h2time->FindFirstBinAbove(h2time->GetMaximum()*0.99));



/////// Write the output files //////////////////////////////////////////////////////////

cout<<"fin de l'analyse"<<endl;
ofstream fileprob(Form("fileprob_%d.txt",runnumber));
cout<<probrefwf<<" "<<nsampwf<<" "<<ndataprob<<" "<<probnsamp<<" "<<probslot<<" "<<probHMStime<<" "<<h1timemax<<" "<<h2timemax<<endl;
cout<<cputime<<endl;
fileprob<<probrefwf<<" "<<nsampwf<<" "<<ndataprob<<" "<<probnsamp<<" "<<probslot<<" "<<probHMStime<<" "<<h1timemax<<" "<<h2timemax<<endl;
fileprob<<cputime<<endl;

  fout->cd();
  fout->Write();
  h1time->Write("h1time");
  h2time->Write("h2time");
  histo->Write("histo");
  fout->Close();



}
