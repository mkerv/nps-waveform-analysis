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
#include <ROOT/RDataFrame.hxx>
#include <vector>
#include <TH1.h> 
#include <ROOT/RDF/HistoModels.hxx> 
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "Fit/UnBinData.h"
#include "Fit/BinData.h"
#include "Fit/DataOptions.h"
#include "Math/WrappedMultiTF1.h"
#include <TKey.h>
#include "TThread.h"


template <typename T>
void checkType(const T &obj)
{
  static_assert(std::is_same<T, std::vector<Double_t>>::value, "Type is not correct!");
}

const int ntime = 110;           // number of time samples for each fADC channel (100 for runs 55 and 56)
const int nfADC = 7;             // number of fADC modules (16 channels each) used in this run
const int nchannel = 16;         // number of channels in each fADC module
const int ncol = 30;             // number of calo colomns tested in this run (3 for runs 55 and 56)
const int nlin = 36;             // number of calo blocks in each colomn
const int nblocks = ncol * nlin; // number of tested calo blocks
const int maxpulses = 4;         // nb maximal de pulses dans la fentre ntime
const int nsignals = nblocks * maxpulses;
const int maxwfpulses = 12; // nb maximal de pulses que la wfa peut trouver
const int ntemp = 56;       // nb of temp sensors
std::atomic<int> nFitFailures{0};

Double_t timeref[nblocks];
Float_t cortime[nblocks];
Int_t preswf[nblocks];
//Double_t timerefacc = -4.5; // car les bons pulses elastiques poussent vers 45.5 (4*ns) dans la wf alors que ceux des production runs sont vers 35.5 (4ns)
Double_t timerefacc = 0;   // en (-4ns) car le pic de timewf pousse vers -22ns

// ROOT::Math::Interpolator *interpolation[nblocks]; // function used for the interpolation
std::vector<std::vector<double>> interpX(nblocks),
                                interpY(nblocks);

///////////////////////////////////////////////////////////// FIT FUNCTION USED AND SCHEMATICS /////////////////////////////////////////////////////////////
bool FastCloneAndFilter(const TString &inName,
  const TString &outName) {
TFile *fin  = TFile::Open(inName,  "READ");
if (!fin || fin->IsZombie()) return false;
TFile *fout = TFile::Open(outName, "RECREATE");
if (!fout|| fout->IsZombie()) { fin->Close(); return false; }

// copy all non-T keys
for (auto k : *fin->GetListOfKeys()) {
TKey *key = static_cast<TKey*>(k);
TString name(key->GetName());
if (name == "T") continue;
//if (key->GetName() == "T") continue;
fout->cd();
key->ReadObj()->Write();
}
// clone T without the big branch
TTree *tin = (TTree*)fin->Get("T");
tin->SetBranchStatus("NPS.cal.fly.adcSampWaveform", 0);
TTree *tout = tin->CloneTree(-1, "fast");
fout->cd();
tout->Write("T");

fout->Close();
fin->Close();
return true;
}
/////////////////////////////////////////////////////////////MAIN FUNCTION ///////////////////////////////////////////////////////////////////////
void TEST_2(int run, int seg, int threads)
{
  TStopwatch t;
  t.Start();
  int nthreads = 6;// or any number
  nthreads = threads;
  

  // BUILD TCHAIN
  TChain chain("T");
  TString filename = Form("/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_%d_%d_1_-1.root", run, seg);
  //TString filename = Form("/mss/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", run, seg);
  //TString filename = Form("../nps_hms_coin_%d_%d_1_-1.root", run, seg);
  TFile *testOpen = TFile::Open(filename);
  if (!testOpen || testOpen->IsZombie())
  {
    std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
    return;
  }
  testOpen->Close();

TString outFile = Form("/volatile/hallc/nps/kerver/ROOTfiles/WF/nps_production_%d_%d_%d_interactive_infheap_test_sorted.root", run, seg,nthreads);

if (!FastCloneAndFilter(filename, outFile)) {
  std::cerr<<"Clone failed!\n";
  return;
}
cout <<"I/O Copy time = "<<t.RealTime()<<endl;
t.Continue();


  // ENABLE MT
  ROOT::EnableImplicitMT(nthreads);
  std::cout << "Implicit MT enabled: " << ROOT::IsImplicitMTEnabled() << "\n";
  std::cout << "Number of threads: " << ROOT::GetThreadPoolSize() << "\n";




  chain.Add(filename);
  chain.SetBranchStatus("*", 0);
// 2) …then enable _just_ the ones your analysis uses:
chain.SetBranchStatus("g.evnum",1);
chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampWaveform",1);
chain.SetBranchStatus("Ndata.NPS.cal.fly.adcCounter",1);
chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseAmp",1);
chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseInt",1);
chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPed",1);
chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTime",1);
chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTimeRaw",1);
chain.SetBranchStatus("NPS.cal.fly.adcSampWaveform",1);
chain.SetBranchStatus("NPS.cal.fly.adcCounter",1);
chain.SetBranchStatus("NPS.cal.fly.adcSampPulseAmp",1);
chain.SetBranchStatus("NPS.cal.fly.adcSampPulseInt",1);
chain.SetBranchStatus("NPS.cal.fly.adcSampPed",1);
chain.SetBranchStatus("NPS.cal.fly.adcSampPulseTime",1);
chain.SetBranchStatus("NPS.cal.fly.adcSampPulseTimeRaw",1);


  auto nEntries = chain.GetEntries(); 
  std::cout << "Chain has " << nEntries << " entries.\n";
  
  // Before processing events, check the list of objects in memory
  std::cout << "Objects in memory before event processing:" << std::endl;
  gObjectTable->Print();

  // Create an RDataFrame from the chain
  ROOT::RDataFrame df(chain);
  auto nEntriesDF = df.Count().GetValue();
  std::cout << "Dataframe has " << nEntriesDF << " entries.\n";

  // Define a new RDataFrame that processes only 1% of the events
  auto nEventsToProcess = nEntriesDF / 1000;
  //auto df1percent = df.Range(0, nEventsToProcess);
  //auto df1percent = df.Range(9999, 12000);

  Double_t dt = 4.;                                    // time bin (sample) width (4 ns), the total time window is then ntime*dt
  const int nslots = 1104;                             // nb maximal de slots dans tous les fADC
  const int Ndata = nslots * (ntime + 1 + 1);          //(1104 slots fADC au total mais pas tous utilises y compris 2 PM scintillateurs?) should correspond to Ndata.NPS.cal.fly.adcSampWaveform variable in the root file
  Double_t ADCtomV = 1000. / 4096;                     // 1000 mV correspond a 4096 canaux
  Double_t integtopC = (1. / 4096) * (dt * 1.e3 / 50); // (1 V / 4096 adc channels) * (4000 ps time sample / 50 ohms input resistance) = 0.020 pc/channel

  // Read reference waveforms
  const int nbparameters = 2 + 2 * maxwfpulses;
  cout << "nbparameters = " << nbparameters << endl;
  Double_t x[ntime], y[ntime];
  ifstream filewf;
  Double_t dum1, ymax;

  // Read tdc_offset_param (needed to determine HMS corrections to the timing)//For now, it is just one file
  Float_t tdcoffset[nblocks];

  ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/tdc_offset_param.txt"); // Done

  for (Int_t i = 0; i < nblocks; i++)
  {

    filetdc >> tdcoffset[i];

    if (run > 6183 && run < 7500)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6171-6183/fit_e_runs/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 6168 && run < 6171)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 5236 && run < 6151)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 5208 && run < 5217)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5183-5208/fit_e_runs/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 3898 && run < 5183)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 2920 && run < 3883)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2900-2920/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 2885 && run < 2900)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2875-2885/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 2871 && run < 2875)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2855-2871/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 1982 && run < 2855)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/ref_wf_%d.txt", i));
    }
    else if (run > 1560 && run < 1961)
    {
      filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1423-1511/RWF/ref_wf_%d.txt", i));
    }

    // comparison with the mean of the ref shapes
    // filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/fit_elastic_runs/results_elastics/refwf/final_refwf_block_%d.txt",i)); Done

    timeref[i] = -1.e6;
    preswf[i] = 0;
    //interpolation[i] = new ROOT::Math::Interpolator(ntime, ROOT::Math::Interpolation::kCSPLINE);

    ymax = 0;
    if (filewf.is_open()){
    filewf >> timeref[i] >> dum1;
    interpX[i].resize(ntime);
    interpY[i].resize(ntime);
    for (int it = 0; it < ntime; ++it) {
      filewf >> interpX[i][it] >> interpY[i][it];
              
                if (interpY[i][it] > ymax)
                {
                    ymax = interpY[i][it];
                    timeref[i] = interpX[i][it];
                }
            
    }
    preswf[i]=1;
   }
    filewf.close();

  }

     // Lecture des corrections timing
     ifstream filetime;
     Float_t dum;
     filetime.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/filetime_step_i.txt");
     for (Int_t i = 0; i < nblocks; i++)
     {
         filetime >> dum >> cortime[i] >> dum >> dum >> dum;
         if(cortime[i]==0){
          cortime[i]=-0.0000001;
         }
     }
     filetime.close();

  // Load only required branches into dataframe
  auto df2 = df.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")
  //  auto df2 = df1percent.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")
                 .Define("SampWaveForm", "NPS.cal.fly.adcSampWaveform")
                 .Define("NadcCounter", "Ndata.NPS.cal.fly.adcCounter")
                 .Define("adcCounter", "NPS.cal.fly.adcCounter")
                 .Define("NadcSampPulseAmp", "Ndata.NPS.cal.fly.adcSampPulseAmp")
                 .Define("adcSampPulseAmp", "NPS.cal.fly.adcSampPulseAmp")
                 .Define("NadcSampPulseInt", "Ndata.NPS.cal.fly.adcSampPulseInt")
                 .Define("adcSampPulseInt", "NPS.cal.fly.adcSampPulseInt")
                 .Define("NadcSampPulsePed", "Ndata.NPS.cal.fly.adcSampPed")
                 .Define("adcSampPulsePed", "NPS.cal.fly.adcSampPed")
                 .Define("NadcSampPulseTime", "Ndata.NPS.cal.fly.adcSampPulseTime")
                 .Define("adcSampPulseTime", "NPS.cal.fly.adcSampPulseTime")
                 .Define("NadcSampPulseTimeRaw", "Ndata.NPS.cal.fly.adcSampPulseTimeRaw")
                 .Define("adcSampPulseTimeRaw", "NPS.cal.fly.adcSampPulseTimeRaw")
                 .Define("evt", "g.evnum"); // define event number from row of rdataframe for troubleshoooting (threadsafe
  // Output rootfile
  // Other variables
  Int_t ilinc, icolc, inp;
  Int_t ndataprob = 0;

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.015);
  TH1::AddDirectory(kFALSE);



Double_t calodist = 9.5;

if (run > 1571 && run < 3667)
{
    calodist = 3.5;
}
else if (run > 3666 && run < 4632)
{
    calodist = 4.;
}
else if (run > 4635 && run < 4953)
{
    calodist = 6.;
}
else if (run > 4965 && run < 5344)
{
    calodist = 4.;
}
else if (run > 5354 && run < 5464)
{
    calodist = 3.;
}
else if (run> 5523 && run < 7013)
{
    calodist = 3.5;
}
timerefacc = (calodist - 9.5) / (3.e8 * 1.e-9 * 4);
  // Read the mean time positions of cosmic pulses (this is determined by the macro analyse_wassim.C)
  Float_t timemean2[nblocks];
  for (Int_t ii = 0; ii < nblocks; ii++)
  {
    timemean2[ii] = 170 + timerefacc * dt;
  } // 1st pass analysis, when the macro analyse_wassim.C is not executed yet (30 is the default value)
  // 2nd pass analysis, when the macro analyse_wassim.C is already executed. If not, comment the following lines

  ROOT::RDF::TH1DModel m_h1time("h1time", "pulse (>20mV) shift (4*ns units) relatively to elastic refwf (all found pulses included)", 200, -50, 50);
  ROOT::RDF::TH1DModel m_h2time("h2time", "pulse (>20mV) time (ns) (all found pulses included)", 200, -100, 100);
 
 //Fill with finter[blocknum]->GetParameter(2 + 2 * p), 1.) and wfttime[blocknumn][numpulses] 
  // auto h1time = df2.Histo1D(m_h1time ,"");
// auto h2time = df2.Histo1D(m_h2time ,"");

  /////// ANALYZE //////////////////////////////////////////////////////////

  // this is where sequential event loop was
  ////////Lambda funtion for the per-event wf analysis/////////
  auto analyze = [=](Int_t NSampWaveForm, const ROOT::VecOps::RVec<Double_t> &SampWaveForm, Double_t evt, Int_t NadcCounter, ROOT::VecOps::RVec<Double_t> &adcCounter, const ROOT::VecOps::RVec<Double_t> adcSampPulseTime, const ROOT::VecOps::RVec<Double_t> adcSampPulseTimeRaw, const ROOT::VecOps::RVec<Double_t> adcSampPulseAmp, const ROOT::VecOps::RVec<Double_t> adcSampPulseInt, const ROOT::VecOps::RVec<Double_t> adcSampPulsePed)
  {
    

  TStopwatch tlambda;
  tlambda.Start();

	  //TF1 *finter[nblocks];
    std::vector<std::unique_ptr<TF1>> finter(nblocks); //object is now thread-local
   // TH1F *hsig_i[nblocks];
   std::vector<TH1F*> hsig_i(nblocks, nullptr);
    //Int_t pres[nblocks] = {0};
    std::vector<Int_t> pres(nblocks, 0);
  //  Double_t signal[nblocks][ntime];
  std::vector<Double_t> signal(nblocks * ntime);

    Int_t bloc, nsamp;

    Double_t enertot = 0.;
    Double_t integtot = 0.;
    Double_t corr_time_HMS = 0.;

    std::vector<Double_t> timewf(nblocks, -100);
    std::vector<Double_t> amplwf(nblocks, -100);
    std::vector<Double_t> chi2(nblocks, -100.0);
  //  std::vector<Double_t> h1time(nblocks, -1.0);
  //  std::vector<Double_t> h2time(nblocks, -1.0);
    std::vector<Double_t> h1time;
    std::vector<Double_t> h2time;
    Double_t ener[nblocks];
    Double_t integ[nblocks];
    Double_t noise[nblocks];
    Double_t bkg[nblocks];
    Double_t sigmax[nblocks];
    std::vector<Double_t> Sampampl(nblocks, -100);
    Double_t Sampped[nblocks];
    std::vector<Double_t> Samptime(nblocks, -100);
    Double_t Sampener[nblocks];
    Double_t Npulse[nblocks];
    Int_t nsampwf = 0;

    std::vector<Int_t> wfnpulse(nblocks, -100);

    

// Now create an outer RVec that contains nblocks copies of tmp:
ROOT::RVec<Double_t> wfampl(nblocks * maxwfpulses, -100.0);
ROOT::RVec<Double_t> wftime(nblocks * maxwfpulses, -100.0);


    // not used but here they are anyway
    Double_t ampl2[nblocks];
    std::vector<Double_t> ampl(nblocks, -100);
    Double_t time[nblocks];
    Double_t larg50[nblocks];
    Double_t larg90[nblocks];
    Double_t max50, max90, min50, min90;
    // where the cosmic pulse is expected to be
    Int_t binmin;
    Int_t binmax;

    // The fit function need to be moved into the scope of analyze to capture all variables (wfnpulse)
    auto Fitwf = [=, &wfnpulse, &wfampl, &finter, &hsig_i, &wftime, &corr_time_HMS, &chi2](Double_t evt, Int_t bn)
    {
     // std::cout << ">>> Fitwf called for evt="<<evt<<" bn="<<bn<<std::endl;

  // Build a thread‐local Interpolator
  auto interpPtr = std::make_shared<ROOT::Math::Interpolator>(
    ntime,
    ROOT::Math::Interpolation::kCSPLINE
 );

   // copy your raw data in:
   interpPtr->SetData(ntime,
    interpX[bn].data(),
    interpY[bn].data());



      // Clone the interpolator for this block
     // ROOT::Math::Interpolator localInterp = *interpolation[bn];
      // (this copy brings along all the precomputed spline data)


      auto func = [interpPtr] (Double_t *x, Double_t *par) mutable {
        Int_t j = Int_t(par[0]); // block index
        Double_t val = 0;
        for (Int_t p = 0; p < maxwfpulses; ++p) {
         Double_t dt0 = x[0] - par[2 + 2*p];
        //  Double_t raw = x[0] - par[2 + 2*p];
        //  Double_t dt0 = std::min(std::max(raw, 0.0), double(ntime-1));
          if (dt0 > 1 && dt0 < ntime - 1) {
            val += par[3 + 2*p] * interpPtr->Eval(dt0);
          }
        }
        return val + par[1];
      };


      std::string fname = Form("finter_bn%d_ptr%p_evt%.0f", bn, (void*)finter[bn].get(), evt);
      finter[bn] = std::make_unique<TF1>(fname.c_str(), func, -200, ntime + 200, nbparameters);
      finter[bn]->FixParameter(0, bn);
      finter[bn]->SetNpx(1100);
      // Detect the pulses

      wfnpulse[bn] = 0;
      Int_t good = 0;


      if (!finter[bn])
      {
        cout << " block=" << bn << " finter is nullptr" << endl;
        return;
      }
      if (!hsig_i[bn])
      {
        cout << ", block=" << bn << " hsig_i is nullptr" << endl;
        return;
      }
      double maxamp=0.0;
      for (Int_t it = 1; it < ntime - 6; it++) // loop over number of samples(110) in an event for a given block
      {
        if(maxamp<hsig_i[bn]->GetBinContent(it)){
          maxamp=hsig_i[bn]->GetBinContent(it);
          }
      

        // Condition over the number of samples in the pulse finding scheme
        if (hsig_i[bn]->GetBinContent(it) < hsig_i[bn]->GetBinContent(it + 1) && hsig_i[bn]->GetBinContent(it + 1) < hsig_i[bn]->GetBinContent(it + 2) && hsig_i[bn]->GetBinContent(it + 2) < hsig_i[bn]->GetBinContent(it + 3) && hsig_i[bn]->GetBinContent(it + 3) <= hsig_i[bn]->GetBinContent(it + 4) && hsig_i[bn]->GetBinContent(it + 4) >= hsig_i[bn]->GetBinContent(it + 5) && hsig_i[bn]->GetBinContent(it + 5) >= hsig_i[bn]->GetBinContent(it + 6))
       //above is new condition of 5 increasing samples and 2 decresaing. below is old: 4 up 1 down.
        // if (hsig_i[bn]->GetBinContent(it + 1) < hsig_i[bn]->GetBinContent(it + 2) && hsig_i[bn]->GetBinContent(it + 2) < hsig_i[bn]->GetBinContent(it + 3) && hsig_i[bn]->GetBinContent(it + 3) <= hsig_i[bn]->GetBinContent(it + 4) && hsig_i[bn]->GetBinContent(it + 4) >= hsig_i[bn]->GetBinContent(it + 5))
        {
          if(hsig_i[bn]->GetBinContent(it + 4) > 2.0)
          //temp always true to match wassims's 5amp criteria
          //if (hsig_i[bn]->GetBinContent(it + 4) > 0)
          {
            // check if we exceeded the number of pulses
            if (wfnpulse[bn] >= maxwfpulses)
            {
              cout << "Warning: wfnpulse[" << bn << "] exceeded maxwfpulses!" << endl;
                wfnpulse[bn] = maxwfpulses - 1; // Prevent overflow
            }

            wfampl[bn * maxwfpulses + wfnpulse[bn]] = hsig_i[bn]->GetBinContent(it + 4); // get the amplitude of the pulse found

            wftime[bn * maxwfpulses + wfnpulse[bn]] = hsig_i[bn]->GetBinCenter(it + 4); // get the time of the pulse found

            // flag for the good pulse
           //if (TMath::Abs(wftime[bn][wfnpulse[bn]] - timeref[bn] +(corr_time_HMS - cortime[bn]) / dt - timerefacc) < 4.)
            if (TMath::Abs(wftime[bn * maxwfpulses + wfnpulse[bn]] - timeref[bn]) < 4.)
            {
              good = 1;
             //cpulse = wfnpulse[bn];
            }
             //   cout<<"detection bloc "<<bn<<" pulse="<<wfnpulse[bn]<<" time="<<wftime[bn][wfnpulse[bn]]<<" (4ns) ampl="<<wfampl[bn][wfnpulse[bn]]<<" reftime= "<<timeref[bn]<<" (4ns) diff="<<TMath::Abs(wftime[bn][wfnpulse[bn]]-timeref[bn])<<" good="<<good<<endl;

            wfnpulse[bn]++;
           // cout<<"detection bloc "<<bn<<" pulse="<<wfnpulse[bn]<<" time="<<wftime[bn][wfnpulse[bn]-1]<<" (4ns) ampl="<<wfampl[bn][wfnpulse[bn]-1]<<" reftime= "<<timeref[bn]<<" (4ns) diff="<<TMath::Abs(wftime[bn][wfnpulse[bn]-1]-timeref[bn])<<" good="<<good<<endl;

            // to prevent overflow
            if (wfnpulse[bn] == maxwfpulses)
            {
              wfnpulse[bn] = maxwfpulses - 1;
              it = ntime;
            }
            maxamp=hsig_i[bn]->GetBinContent(it+4);
            it += 4;

          } // end of the condition over (hsig_i[bn]->GetBinContent(it+4)>0){
        } // end of the condition over samples
      } // end of loop over it

      if (maxamp < 2.0)
      {
       // cout << "Warning : All samples < 2.0mV. Skipping " << bn << endl;
        return;
      }


      if (wfnpulse[bn] > maxwfpulses - 2)
      {
        cout << "Warning : excessively high number of pulses in the block wf " << bn << endl;
      }

      // Adjust the parameters of the fit function

      for (Int_t p = 0; p < maxwfpulses; p++)
      {
        finter[bn]->FixParameter(2 + 2 * p, 0.);
        finter[bn]->FixParameter(3 + 2 * p, 0.);
      }


      if (wfnpulse[bn] > 0 && good == 1)
      {
        for (Int_t p = 0; p < TMath::Min(maxwfpulses, wfnpulse[bn]); p++)
        {
          finter[bn]->ReleaseParameter(2 + 2 * p);
          finter[bn]->ReleaseParameter(3 + 2 * p);

    //finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn]);
          //finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt);   
          finter[bn]->SetParameter(2 + 2 * p, wftime[bn * maxwfpulses + p ]- timeref[bn]);   
          finter[bn]->SetParameter(3 + 2 * p, wfampl[bn * maxwfpulses + p ]);


          //finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] - 3, wftime[bn][p] - timeref[bn] + 3);
          //finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt - 3, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt + 3);
          finter[bn]->SetParLimits(2 + 2 * p, wftime[bn * maxwfpulses + p ]- timeref[bn] - 3, wftime[bn * maxwfpulses + p ]- timeref[bn] + 3);

          finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn * maxwfpulses + p ] * 0.2, wfampl[bn * maxwfpulses + p ] * 3);   

        }
      }

      if (wfnpulse[bn] > 0 && good == 0)
      {
        for (Int_t p = 0; p < TMath::Min(maxwfpulses, wfnpulse[bn]); p++)
        {
          finter[bn]->ReleaseParameter(2 + 2 * p);
          finter[bn]->ReleaseParameter(3 + 2 * p);
          finter[bn]->SetParameter(2 + 2 * p, wftime[bn * maxwfpulses + p ] - timeref[bn]);
         // finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt);
          finter[bn]->SetParameter(3 + 2 * p, wfampl[bn * maxwfpulses + p ]);
          finter[bn]->SetParLimits(2 + 2 * p, wftime[bn * maxwfpulses + p ] - timeref[bn] - 3, wftime[bn * maxwfpulses + p ] - timeref[bn] + 3);
          //finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt - 3, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt + 3);

          finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn * maxwfpulses + p ] * 0.2, wfampl[bn * maxwfpulses + p ] * 3);
        }

        // On recherche quand meme un eventuel pulse en temps
        //cpulse = wfnpulse[bn];

        finter[bn]->ReleaseParameter(2 + 2 * wfnpulse[bn]);
        finter[bn]->ReleaseParameter(3 + 2 * wfnpulse[bn]);
         finter[bn]->SetParameter(2 + 2 * wfnpulse[bn], 0.);
       // finter[bn]->SetParameter(2 + 2 * wfnpulse[bn], timerefacc);
        //finter[bn]->SetParameter(2 + 2 * wfnpulse[bn], timerefacc + (corr_time_HMS - cortime[bn]) / dt);
        finter[bn]->SetParameter(3 + 2 * wfnpulse[bn], 2);
         // finter[bn]->SetParLimits(2 + 2 * wfnpulse[bn],  - 1, 1);
      //  finter[bn]->SetParLimits(2 + 2 * wfnpulse[bn], timerefacc - 4, timerefacc + 4);
        //finter[bn]->SetParLimits(2 + 2 * wfnpulse[bn], timerefacc + (corr_time_HMS - cortime[bn]) / dt - 4, timerefacc + (corr_time_HMS - cortime[bn]) / dt + 4);

        finter[bn]->SetParLimits(3 + 2 * wfnpulse[bn], 0.05, 10);

        wfnpulse[bn]++;
       
      }

      if (wfnpulse[bn] == 0)
      {
       // cout<<"NO PULSE "<<bn<<" pulse="<<wfnpulse[bn]<<" time="<<wftime[bn][wfnpulse[bn]]<<" (4ns) ampl="<<wfampl[bn][wfnpulse[bn]]<<" reftime= "<<timeref[bn]<<" (4ns) diff="<<TMath::Abs(wftime[bn][wfnpulse[bn]]-timeref[bn])<<" good="<<good<<endl;

        for (Int_t p = 0; p < 1; p++)
        {
          finter[bn]->ReleaseParameter(2 + 2 * p);
          finter[bn]->ReleaseParameter(3 + 2 * p);

          finter[bn]->SetParameter(2 + 2 * p, 0);
          //finter[bn]->SetParameter(2 + 2 * p, timerefacc);
          //finter[bn]->SetParameter(2 + 2 * p, timerefacc + (corr_time_HMS - cortime[bn]) / dt);
          finter[bn]->SetParameter(3 + 2 * p, 2);

         // finter[bn]->SetParLimits(2 + 2 * p,  - 1, 1);
          //finter[bn]->SetParLimits(2 + 2 * p, timerefacc - 4, timerefacc + 4);
          //finter[bn]->SetParLimits(2 + 2 * p, timerefacc + (corr_time_HMS - cortime[bn]) / dt - 4, timerefacc + (corr_time_HMS - cortime[bn]) / dt + 4);

          finter[bn]->SetParLimits(3 + 2 * p, 0.05, 10);

       // finter[bn]->FixParameter(2 + 2 * p, timerefacc + (corr_time_HMS - cortime[bn]) / dt);
       finter[bn]->FixParameter(2 + 2 * p, 0);

          finter[bn]->FixParameter(3 + 2 * p, 2);
          finter[bn]->ReleaseParameter(1);
        }

        wfnpulse[bn]++;
      }

      finter[bn]->SetParameter(1, 0.);
      finter[bn]->SetParLimits(1, -100, 100.);
/*
      std::cout<<"  [DEBUG] before final fix: wfnpulse="<<wfnpulse[bn]
         <<"  maxwfpulses="<<maxwfpulses<<std::endl;
*/


      for (int p = wfnpulse[bn]; p < maxwfpulses; ++p) {
        finter[bn]->FixParameter(2 + 2*p, 0.);
        finter[bn]->FixParameter(3 + 2*p, 0.);
      }

//Prepare binned data from the histogram:
int nBins = hsig_i[bn]->GetNbinsX();
ROOT::Fit::BinData data(nBins, 1);
for (int ib = 1; ib <= nBins; ++ib) {
    double x[1] = { hsig_i[bn]->GetBinCenter(ib) };
    double y    = hsig_i[bn]->GetBinContent(ib);
    double err  = hsig_i[bn]->GetBinError(ib);
    data.Add(x, y, err);
}




//Wrap TF1 into a IModelFunction via WrappedMultiTF1:
ROOT::Math::WrappedMultiTF1 wfunc(*finter[bn], finter[bn]->GetNdim());
//Configure fitter
ROOT::Fit::Fitter fitter;
auto &cfg = fitter.Config();
//fitter.Config().SetMinimizer("Minuit2", "Migrad");
//cfg.SetMinimizer("Minuit2", "Migrad");
cfg.SetMinimizer("Minuit2", "MINIMIZE");
//auto &mopts = fitter.Config().MinimizerOptions();
auto &mopts = cfg.MinimizerOptions();
mopts.SetStrategy(0);
mopts.SetPrintLevel(0);
mopts.SetMaxIterations(1000);
fitter.SetFunction(wfunc, false);


cfg.CreateParamsSettings(wfunc);

int nFloat = 2 + 2*wfnpulse[bn];
for (int ip = nFloat; ip < nbparameters; ++ip)
cfg.ParSettings(ip).Fix();  
for (int ip = 1; ip < nFloat; ++ip)
cfg.ParSettings(ip).Release();

//if(good == 0 && wfnpulse[bn]>1){
  if(good == 0){
  //cout<<"SET LIMITS HERE "<<wfnpulse[bn]<< endl;  
  cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).SetLowerLimit(0.05);
  cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).SetUpperLimit(10);
  if(wfnpulse[bn]>1)cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).Fix();        // pin it at your seed = 0.05

   finter[bn]->SetParameter(3 + 2*(wfnpulse[bn]-1), 2.);
  }

cfg.ParSettings(0).Fix();  



bool ok = fitter.LeastSquareFit(data);
if (!ok) {
  std::cerr<<"Failed once for event "<<evt<<", block "<<bn<<"\n";
}

cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).Release();  
if(good == 0 && wfnpulse[bn]>1 && ok){
  //if a real pulse is found save the fit result first before releasing the ref pulse params
  auto &tempresult = fitter.Result();
  chi2[bn] = tempresult.Chi2()/tempresult.Ndf();
      
  for (Int_t p = 0; p < wfnpulse[bn]; ++p) {
    wfampl[bn * maxwfpulses + p ] = tempresult.Parameter(3 + 2*p);
    wftime[bn * maxwfpulses + p ]  = tempresult.Parameter(2 + 2*p)*dt                     // convert bins → ns
               + corr_time_HMS                // add HMS correction
               - cortime[bn]                  // subtract block‐by‐block cable delay
               - timerefacc*dt;               // subtract your reference‐time offset
}  

// cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).Release(); 
  ok = fitter.LeastSquareFit(data);
  }


if (!ok) {
  // retry just this one with a tougher configuration
 // cout<<"FAILED Again, retry"<<endl;
  cfg.MinimizerOptions().SetStrategy(2);
 // mopts.SetPrintLevel(1);
  cfg.MinimizerOptions().SetMaxIterations(5000);
  ok = fitter.LeastSquareFit(data);
}
if (!ok) {
//  std::cerr<<"STILL Fit failed for event "<<evt<<", block "<<bn<<"\n";
 //failcount++;
 nFitFailures.fetch_add(1, std::memory_order_relaxed);
 if(good == 0 && wfnpulse[bn]>1) return;
 for (int p = 0; p < wfnpulse[bn]; p++)
  {
  wftime[bn * maxwfpulses + p ] = -1000;
  wfampl[bn * maxwfpulses + p ] = -1000;
  }
  

  return;
}

//Extract parameters to arrays
auto &result = fitter.Result();

      for (Int_t p = 0; p < wfnpulse[bn]; ++p) {
        // 1) the fitted bin shift (same units your x[] was in)
        double binOff = result.Parameter(2 + 2*p);
    
    
        // 2) what bin *would* that be in the original histogram?
        double binIndex = binOff + timeref[bn];       // timeref[] is itself in bins
    
        // 3) uncorrected time in ns
        double uncorTime = binIndex * dt;             // dt = 4 ns/bin
    
        // 4) fully corrected time in ns
        double corrTime = uncorTime + (corr_time_HMS - cortime[bn]);
    
        // 5) amplitude
        wfampl[bn * maxwfpulses + p ] = result.Parameter(3 + 2*p);
        wftime[bn * maxwfpulses + p ]  = binOff*dt                     // convert bins → ns
                   + corr_time_HMS                // add HMS correction
                   - cortime[bn]                  // subtract block‐by‐block cable delay
                   - timerefacc*dt;               // subtract your reference‐time offset

    }
  unsigned npar = result.NPar();
  for(unsigned ip=0; ip<npar; ++ip) {
    finter[bn]->SetParameter(ip, result.Parameter(ip));
  }

  double tempchi2   = result.Chi2();   // total χ² from the fit :contentReference[oaicite:0]{index=0}
unsigned int ndf = result.Ndf(); // degrees of freedom :contentReference[oaicite:1]{index=1}

chi2[bn] = tempchi2/ndf;
      
      // Old Fit
     // hsig_i[bn]->Fit(finter[bn].get(), "Q", "", 1., ntime);

      for(int p=0; p<26;p++){
        //std::cout<<p<<" : "<<finter[bn]->GetParameter(p)<<endl;
     //   std::cout<<p<<" : "<<result.Parameter(3 + 2*p)<<endl;
    }

      for(int np=0;np<wfnpulse[bn];np++){
       //std::cout<<wfampl[bn][np]<<"  "<<wftime[bn][np]<<endl;
      }

    };

    if (NSampWaveForm > Ndata)
    {
      nsampwf++;
      cout << "!!!! NSampWaveForm problem  " << evt << "  " << NSampWaveForm << " " << Ndata << endl;
    }

    if (NSampWaveForm <= Ndata) // NSampWaveForm must be <= Ndata (otherwise correct Ndata value)
    {
      for (Int_t i = 0; i < nblocks; i++)
      {
        // Initialize per-block Arrays and Histos for WFA
        ener[i] = 0;
        integ[i] = 0;
        noise[i] = 0;
        bkg[i] = 0;
        sigmax[i] = -100;
        Sampped[i] = -100;
        Sampener[i] = -100;
        Npulse[i] = 0;


      // avoid threads creating hsig with same name at same time(before object deletion)
      TThread* me = TThread::Self();
      Long_t rootTid = me ? me->GetId() : -1;
      TString hname = Form("hsig_i_evt%.0f_thr%ld_blk%d", evt, rootTid, i);
        hsig_i[i] = new TH1F(hname, hname, ntime, 0, ntime);
        //hsig_i[i] = new TH1F(Form("hsig_i%d", i), Form("hsig_i%d", i), ntime, 0, ntime);
        hsig_i[i]->SetLineColor(1);
        hsig_i[i]->SetLineWidth(2);
        hsig_i[i]->GetXaxis()->SetTitle("Time (4 ns)");
        hsig_i[i]->GetYaxis()->SetTitle("(mV)");
        hsig_i[i]->GetYaxis()->SetLabelSize(0.05);



      } // liste de presence des blocs! pres=0 si bloc absent, pres=1 s'il est present
      std::fill(signal.begin(), signal.end(), 0.);
      // Extract the data from the variable SampWaveForm[] (NPS.cal.fly.adcSampWaveform)

      int ns = 0; // ns represent for a given event the element number of the NPS.cal.fly.adcSampWaveform variable
      while (ns < NSampWaveForm)
      {
        bloc = SampWaveForm[ns];
        ns++; // bloc number (actually the slot number)
        nsamp = SampWaveForm[ns];
        ns++; // time samples (should be equal to ntime (100))

        if (bloc == 2000)
          bloc = 1080; // modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators
        if (bloc == 2001)
          bloc = 1081; // modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators

        if (bloc < 0 || bloc > nslots - 0.5)
        {
          cout << "slot number problem " << evt << " " << bloc << endl;
          // ns = NSampWaveForm + 1;
          break;
        } // to exit the while()

        if (bloc > -0.5 && bloc < nslots)
        { // that's what we expect!

          pres[bloc] = 1;

          //wassims updated code doesnt have this requirement
          //if (nsamp == ntime)
          //{
            for (Int_t it = 0; it < nsamp; it++)
            {
              if (bloc > -0.5 && bloc < nblocks)
              {
                //signal[bloc][it] = SampWaveForm[ns];
                signal[ bloc * ntime + it ] = SampWaveForm[ns];
                //hsig_i[bloc]->SetBinContent(it + 1, signal[bloc][it]);
                hsig_i[bloc]->SetBinContent(it + 1, signal[bloc*ntime + it]);
              }
              ns++;
            }
         // }
        } // fin if(bloc number is good)
      } // fin while()

      // Read the hcana calculated variables

      // if(!(NadcCounter==NadcSampPulseAmp&&NadcSampPulseInt==NadcSampPulsePed&&NadcSampPulseInt==NadcSampPulseTime)){cout<<"!!!!! Problem Ndata !!!!!! "<<evt<<endl;ndataprob++;}

      for (Int_t iNdata = 0; iNdata < NadcCounter; iNdata++)
      {
        if (adcCounter[iNdata] == 2000)
          adcCounter[iNdata] = 1080; // nouveau numero du scintillateur attribue par Malek
        if (adcCounter[iNdata] == 2001)
          adcCounter[iNdata] = 1081; // nouveau numero du scintillateur attribue par Malek

        // determine HMS time correction
        if (iNdata == 0)
        {
          corr_time_HMS = adcSampPulseTime[iNdata] - (adcSampPulseTimeRaw[iNdata] / 16.) - tdcoffset[(int)(adcCounter[iNdata])];
        }
        if (iNdata != 0)
        {
          if (TMath::Abs(corr_time_HMS - (adcSampPulseTime[iNdata] - adcSampPulseTimeRaw[iNdata] / 16. - tdcoffset[(int)(adcCounter[iNdata])])) > 0.001)
          {
            //cout << "problem HMS time correction event " << evt <<"  "<< corr_time_HMS <<" "<<tdcoffset[(int)(adcCounter[iNdata])]<< endl;
          }
        }

        if (!(adcCounter[iNdata] >= 0 && adcCounter[iNdata] < nblocks + 2) /*&&adcCounter[iNdata]!=196*/)
        {
          cout << "****** Problem adcCounter ******* " << evt << " " << iNdata << " " << adcCounter[iNdata] << endl;
        }
        if (adcCounter[iNdata] >= 0 && adcCounter[iNdata] < nblocks)
        {
          Npulse[(int)(adcCounter[iNdata])] += 1;
          hsig_i[(int)(adcCounter[iNdata])]->SetLineColor(2); // why are we doing this here??

          if (Npulse[(int)(adcCounter[iNdata])] == 1)
          {
            Sampampl[(int)(adcCounter[iNdata])] = adcSampPulseAmp[iNdata];
            Samptime[(int)(adcCounter[iNdata])] = adcSampPulseTime[iNdata];
            Sampener[(int)(adcCounter[iNdata])] = adcSampPulseInt[iNdata];
            Sampped[(int)(adcCounter[iNdata])] = adcSampPulsePed[iNdata];
          }
          if (Npulse[(int)(adcCounter[iNdata])] > 1)
          {
            if (TMath::Abs(Samptime[(int)(adcCounter[iNdata])] - timemean2[(int)(adcCounter[iNdata])]) > TMath::Abs(adcSampPulseTime[iNdata] - timemean2[(int)(adcCounter[iNdata])]))
            { // if the 2nd pulse is closer to the reference time than the 1st pulse then we take the 2nd
              Sampampl[(int)(adcCounter[iNdata])] = adcSampPulseAmp[iNdata];
              Samptime[(int)(adcCounter[iNdata])] = adcSampPulseTime[iNdata];
              Sampener[(int)(adcCounter[iNdata])] = adcSampPulseInt[iNdata];
              Sampped[(int)(adcCounter[iNdata])] = adcSampPulsePed[iNdata];
            }
          }
        }
      }

      // Fit de la wf
          for (Int_t i = 0; i < nblocks; i++)
      {
        if (pres[i] == 1 && preswf[i]==1 ) // if this block is present during event
        {
          for (Int_t it = 0; it < nsamp; it++)
          {
            hsig_i[i]->SetBinError(it + 1, TMath::Sqrt(TMath::Abs(hsig_i[i]->GetBinContent(it + 1) * 4.096/2.)) / 4.096);

            if (hsig_i[i]->GetBinContent(it + 1) < 1.)
            {
              hsig_i[i]->SetBinError(it + 1, TMath::Sqrt(TMath::Abs(1. * 4.096/2.)) / 4.096);
            }
          }

          //if (!hsig_i[i]) continue; 
          //if (wfnpulse[i] == 0) continue; 
            Fitwf(evt,i); // We call the fit here
            finter[i]->SetLineColor(kBlue);
            finter[i]->SetLineWidth(2);
    

          // **Add this check to skip the block if finter[i] is invalid**
          if (!finter[i])
          {
            cout << "Skipping block " << i << " due to missing finter[" << i << "] at event " << evt << endl;
            continue;
          }
          //cout << i << "  cpulse =" << cpulse << " wfnpulse1= " << wfnpulse1[i] << endl;
          //chi2[i] = finter[i]->GetChisquare() / finter[i]->GetNDF();

          for (Int_t p = 0; p < TMath::Max(wfnpulse[i], 1); p++)
          {

           // wftime[i][p] = finter[i]->GetParameter(2 + 2 * p) * dt + corr_time_HMS; // temps du pulse en ns
           // wftime[i][p] = finter[i]->GetParameter(2 + 2 * p) * dt + corr_time_HMS - cortime[i] - timerefacc * dt;

           // wfampl[i][p] = finter[i]->GetParameter(3 + 2 * p); // amplitude du pulse en ns

            if (wfampl[i * maxwfpulses + p ] > 20)
            {

              /// fix later //// diagnostic histos
              
               // h2time->Fill(wftime[i][p], 1.);
               //h2time[i]=wftime[i][p];
               h2time.push_back(wftime[i * maxwfpulses + p ]);

               // h1time->Fill(finter[i]->GetParameter(2 + 2 * p), 1.);
               // h1time[i]=finter[i]->GetParameter(2 + 2 * p);

               h1time.push_back(finter[i]->GetParameter(2 + 2 * p) - timerefacc + corr_time_HMS / dt);
               //h1time.push_back(finter[i]->GsetParameter(2 + 2 * p));
              
            } // Fill the time spectrums

            if (p == 0)
            {
              timewf[i] = wftime[i * maxwfpulses + p ];
              amplwf[i] = wfampl[i * maxwfpulses + p ];
            }

            // if(p>0){if(TMath::Abs(wftime[i][p]-timerefacc/dt)<TMath::Abs(timewf[i]-timerefacc/dt)){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}}//pour prendre le pulse dont le temps est le plus proche de timerefacc

            // New modification based on the production runs

            if (p > 0)
            {
              if (TMath::Abs(wftime[i * maxwfpulses + p ]) < TMath::Abs(timewf[i]))
              {
                timewf[i] = wftime[i * maxwfpulses + p ];
                amplwf[i] = wfampl[i * maxwfpulses + p ];
              }
            } // pour prendre le pulse dont le temps est le plus proche de timerefacc

          } // end of loop over maxpulses
        } // end of condition over col,lin
      } // end of loop over les blocs

      // Calculation of some output branch variables
      for (Int_t i = 0; i < nblocks; i++)
      {

        binmin = 30;
        binmax = 109;

        for (Int_t it = 0; it < nsamp; it++)
        {

          //integ[i] += signal[i][it];
         //integtot += signal[i][it];
         integ[i] += signal[i*ntime + it];
         integtot += signal[i*ntime + it];

          if (it > binmin && it < binmax)
          { // cosmic pulse window

            //ener[i] += signal[i][it];
           // enertot += signal[i][it];
           ener[i] += signal[i*ntime + it];
           enertot += signal[i*ntime + it];
          }

          if (!(it > binmin && it < binmax))
          {
            //bkg[i] += signal[i][it];
            bkg[i] += signal[i*ntime + it];

          } // background window

         // if (signal[i][it] > sigmax[i])
         if (signal[i*ntime + it] > sigmax[i])
          {
            time[i] = it;
           // sigmax[i] = signal[i][it];
            //ampl[i] = signal[i][it];
            sigmax[i] = signal[i*ntime + it];
            ampl[i] = signal[i*ntime + it];

          } // determine the pulse maximum
        }

        // doesnt seem to be used??
        ener[i] -= bkg[i] * (binmax - binmin - 1) / (nsamp - (binmax - binmin - 1)); // subtract the bkg contribution (we have to normalize since the window widths are not the same)

        bkg[i] = bkg[i] / (nsamp - (binmax - binmin - 1)); // mean value of bkg

        for (Int_t it = 0; it < nsamp; it++)
        {

          if (!(it > binmin && it < binmax))
          {
            //noise[i] += (signal[i][it] - bkg[i]) * (signal[i][it] - bkg[i]) / (nsamp - (binmax - binmin - 1));
            noise[i] += (signal[i*ntime+it] - bkg[i]) * (signal[i*ntime+it] - bkg[i]) / (nsamp - (binmax - binmin - 1));
          } // RMS of the bkg
        }
        noise[i] = TMath::Sqrt(noise[i]);

        // Calculation of signal widths
        ampl2[i] = ampl[i] - bkg[i]; // amplitude of the pulse relatively to bkg

        max50 = 0;
        max90 = 50;
        min50 = 100;
        min90 = 100;

        for (Int_t it = time[i]; it < nsamp; it++)
        { // aller vers la droite du maximum
          //if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.5)
          if ((signal[i*ntime + it] - bkg[i]) >= ampl2[i]*0.5)
          {
            max50 = it;
          }
          //if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.1)
          if ((signal[i*ntime + it] - bkg[i]) >= ampl2[i]*0.1)
          {
            max90 = it;
          }
        }
        for (Int_t it = time[i]; it > -0.5; it--)
        { // aller vers la gauche du maximum
          //if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.5)
          if ((signal[i*ntime + it] - bkg[i]) >= ampl2[i]*0.5)
          {
            min50 = it;
          }
          //if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.1)
          if ((signal[i*ntime + it] - bkg[i]) >= ampl2[i]*0.1)
          {
            min90 = it;
          }
        }
        // not used anywhere??
        larg50[i] = max50 - min50;
        larg90[i] = max90 - min90;
      }

      if (evt == 2)
      {
        cout << "SampAmp=" << Sampampl[7] << "  ampl=" << ampl[7] << endl;
      }

      // Lambda returns a tuple of all the arrays to be written to dataframe. this replaces filling branches in ttree
      // treeout->Fill();
if((int)evt % 1000 == 0){
  cout <<" Entry = "<< evt <<"  cpu time="<<tlambda.RealTime()<<endl;
  tlambda.Continue();

       //gObjectTable->Print();

      }

    } // end if(NSampWaveForm<=Ndata)

    //Diagnostic Code for hsig_i
    /*
    int nonzero =0;
    if ((int)evt % 1000 == 0) {
      TFile *fout = new TFile(Form("/volatile/hallc/nps/kerver/ROOTfiles/debug_run_%d_evt_%lld.root",run, (long long)evt), "RECREATE");
      for (int i = 0; i < nblocks; ++i) {
        if (hsig_i[i] && hsig_i[i]->GetEntries() > 0) nonzero++;
        if (hsig_i[i]){
          std::cout << "evt=" << evt << ", block=" << i
          << ", integral=" << hsig_i[i]->Integral() << std::endl;
         hsig_i[i]->Write(Form("hsig_block_%d", i));
        }
      }
      fout->Close();
      delete fout;
    }
      
    //end of diagnostic snippet 
*/

    for (int i = 0; i < nblocks; ++i) {
      if (hsig_i[i]) {
      delete hsig_i[i];
      hsig_i[i] = nullptr;
      }
  }
  
//std::cout << "evt=" << evt << ": " << nonzero << " non-empty histograms\n";
    tlambda.Stop();
    return std::make_tuple(chi2, ampl, amplwf, wfnpulse, Sampampl, Samptime, timewf, enertot, integtot, pres, corr_time_HMS, h1time, h2time, wfampl, wftime);
  }; // End of lambda (event loop)

  /////////////////////////////////////////////////////////////////////////////////////////

  // Removed the chunk responsable for the event per event check for now

  /////// Write the output files //////////////////////////////////////////////////////////

  // Make output dataframe with columns as the output tuple and all the arrays within
  auto df_final = df2.Define("tuple", analyze, {"NSampWaveForm", "SampWaveForm", "evt", "NadcCounter", "adcCounter", "adcSampPulseTime", "adcSampPulseTimeRaw", "adcSampPulseAmp", "adcSampPulseInt", "adcSampPulsePed"});

  // Using auto casued compile error for lambda arguements, had to list them explicitly to work?
  df_final = df_final.Define("chi2", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
   std::vector<Double_t> output = std::get<0>(tuple);
   return output; },{"tuple"});
  df_final = df_final.Define("ampl", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
    std::vector<Double_t> output = std::get<1>(tuple);
    return output; },{"tuple"});
  df_final = df_final.Define("amplwf", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
    std::vector<Double_t> output = std::get<2>(tuple);
    return output; },{"tuple"});
    df_final = df_final.Define("wfnpulse", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
      std::vector<Int_t> output = std::get<3>(tuple);
      return output; },{"tuple"});  
      df_final = df_final.Define("Sampampl", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
        std::vector<Double_t> output = std::get<4>(tuple);
        return output; },{"tuple"});
        df_final = df_final.Define("Samptime", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
          std::vector<Double_t> output = std::get<5>(tuple);
          return output; },{"tuple"});  
          df_final = df_final.Define("timewf", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
            std::vector<Double_t> output = std::get<6>(tuple);
            return output; },{"tuple"});  
            df_final = df_final.Define("enertot", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
              Double_t output = std::get<7>(tuple);
              return output; },{"tuple"});  
              df_final = df_final.Define("integtot", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
                Double_t output = std::get<8>(tuple);
                return output; },{"tuple"});  
                df_final = df_final.Define("pres", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
                  std::vector<Int_t> output = std::get<9>(tuple);
                  return output; },{"tuple"});  
                  df_final = df_final.Define("corr_time_HMS", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
                    Double_t output = std::get<10>(tuple);
                    return output; },{"tuple"});
                    df_final = df_final.Define("h1time", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
                      std::vector<Double_t> output = std::get<11>(tuple);
                      return output; },{"tuple"}); 
                      df_final = df_final.Define("h2time", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
                        std::vector<Double_t> output = std::get<12>(tuple);
                        return output; },{"tuple"}); 
                        df_final = df_final.Define("wfampl", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>   > &tuple){ 
                          ROOT::RVec<Double_t> output = std::get<13>(tuple);
                          return output; },{"tuple"});  
                          df_final = df_final.Define("wftime", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,ROOT::RVec<Double_t> , ROOT::RVec<Double_t>    > &tuple){ 
                            ROOT::RVec<Double_t> output = std::get<14>(tuple);
                            return output; },{"tuple"});    
                                                  
                    
 auto h_h1time = df_final.Histo1D(m_h1time ,"h1time");
 auto h_h2time = df_final.Histo1D(m_h2time ,"h2time");

  // Save the dataframe to the output ROOT file
  //TString rootfilePath = Form("/volatile/hallc/nps/kerver/ROOTfiles/WF/nps_production_%d_%d_%d_interactive_allwf_0.05fix_2mvthresh.root", run, seg,nthreads);
 // TString rootfilePath = Form("../nps_production_wf_%d_%d_%d.root", run, seg,nthreads);
  auto columnNames = df_final.GetColumnNames();
  for (const auto &col : columnNames) {
   // std::cout << col << std::endl;
}

///////////
cout<<"About to Snapshot"<<endl;
TString tempout = Form("../nps_production_%d_%d_%d_temp.root", run, seg,nthreads);
ROOT::RDF::RSnapshotOptions opts;
opts.fMode = "RECREATE";  
//opts.fBasketSize = 512 * 1024;        // 512 KB per basket
df_final.Snapshot("WF", tempout, {"chi2","ampl","amplwf","wfnpulse","Sampampl","Samptime","timewf","enertot","integtot","pres","corr_time_HMS","h1time","h2time","evt","wfampl","wftime"},opts);
t.Stop();
std::cout
  <<"=== Snapshot done. Elapsed real time: "<<t.RealTime()
  <<" s, CPU time: "<<t.CpuTime()
  <<" s ===\n";
t.Continue();

// 1) Open and tune baskets
TFile *fin = TFile::Open(tempout, "READ");
if (!fin || fin->IsZombie()) {
  std::cerr<<"ERROR: could not open "<<tempout<<" for basket-tuning\n";
  return;
}
//TTree *tin = (TTree*)fin->Get("WF");
TTree* tin = static_cast<TTree*>(fin->Get("WF"));
if (!tin) {
  std::cerr<<"ERROR: WF not found in "<<tempout<<"\n";
  fin->Close();
  return;
}
tin->SetBasketSize("wfampl", 512*1024);   // 512 KB
tin->SetBasketSize("wftime", 512*1024);

tin->BuildIndex("evt");


  // 3) Clone into your main file (atomic overwrite)
  TFile *fout = TFile::Open(outFile, "UPDATE");
  if (!fout || fout->IsZombie()) {
    std::cerr<<"Cannot open output file "<<outFile<<"\n";
    fin->Close();
    return;
  }
  fout->cd();
  auto tout = tin->CloneTree(-1);
    // 4) (Optional) tune output baskets before write
tout->SetBasketSize("wfampl", 512*1024);
tout->SetBasketSize("wftime", 512*1024);
  tout->BuildIndex("evt");


std::cout
  <<"=== Sorting done. Elapsed real time: "<<t.RealTime()
  <<" s, CPU time: "<<t.CpuTime()
  <<" s ===\n";
t.Continue();


  tout->Write("WF", TObject::kOverwrite);
  fout->Close();
  fin->Close();

 
  cout << "fin de l'analyse" << endl;
  //cout << "Failed fits: "<<failcount<<endl;
  std::cout << "Total failed fits: " << nFitFailures.load() << "\n";
  t.Stop();
  t.Print();
}
