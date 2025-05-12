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
int failcount=0;

Double_t timeref[nblocks];
Float_t cortime[nblocks];
Int_t preswf[nblocks];
//Double_t timerefacc = -4.5; // car les bons pulses elastiques poussent vers 45.5 (4*ns) dans la wf alors que ceux des production runs sont vers 35.5 (4ns)
Double_t timerefacc = 0;   // en (-4ns) car le pic de timewf pousse vers -22ns

// ROOT::Math::Interpolator *interpolation[nblocks]; // function used for the interpolation
std::vector<std::vector<double>> interpX(nblocks),
                                interpY(nblocks);

///////////////////////////////////////////////////////////// FIT FUNCTION USED AND SCHEMATICS /////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////MAIN FUNCTION ///////////////////////////////////////////////////////////////////////

void wfMT1event(int run, int seg, int selectevent, int selectblock)
{

  // FIT FUNCTION redefined as a lambda function for MT use
  /*
  auto func = [=](Double_t *x, Double_t *par)
  {
    if (!par){
      cout << "NULL PARAM POINTER!" <<endl;
      return 0.0;  // not checking if par[0] is nullptr
    } 

    Int_t j = (int)(par[0]);
    Double_t val = 0;
    for (Int_t p = 0; p < maxwfpulses; p++)
    {
      if (x[0] - par[2 + 2 * p] > 1 && x[0] - par[2 + 2 * p] < 109)
        val += par[3 + 2 * p] * interpolation[j]->Eval(x[0] - par[2 + 2 * p]);
    }
    return val + par[1];
  };
*/
  TStopwatch t;
  t.Start();

  // ENABLE MT
  const int nthreads = 6;// or any number
  //ROOT::EnableImplicitMT(nthreads);
  std::cout << "Implicit MT enabled: " << ROOT::IsImplicitMTEnabled() << "\n";
  std::cout << "Number of threads: " << ROOT::GetThreadPoolSize() << "\n";

  // BUILD TCHAIN
  TChain chain("T");
  TString filename = Form("/cache/hallc/c-nps/analysis/pass2/replays/production/nps_hms_coin_%d_%d_1_-1.root", run, seg);
  //TString filename = Form("/mss/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", run, seg);
  //TString filename = Form("../nps_hms_coin_%d_%d_1_-1.root", run, seg);
  TFile *testOpen = TFile::Open(filename);
  if (!testOpen || testOpen->IsZombie())
  {
    std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
    return;
  }
  testOpen->Close();
  chain.Add(filename);
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
  auto df1percent = df.Range(selectevent-1, selectevent);

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

    // filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/ref_wf_%d.txt",i));
    // filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/ref_wf_%d.txt",i));

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

/*   
    // Fill array with reference waveform (x,y) for each block. set reftime @ ymax
    if (filewf.is_open())
    {
      filewf >> timeref[i] >> dum1; // used for my ref shapes (be careful in switching on the mean method !!!!!!!!!!!!!!!!!)

      ymax = 0;
      for (Int_t it = 0; it < ntime; it++)
      {
        filewf >> x[it] >> y[it];
        if (y[it] > ymax)
        {
          ymax = y[it];
          timeref[i] = x[it];
        }
      }
      interpolation[i]->SetData(ntime, x, y);
      preswf[i]=1;
    }
    filewf.close();
    */
    // cout<<i<<" "<<timeref[i]<<endl;
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
  // auto df2 = df.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")
  auto df2 = df1percent.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")
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
                 .Define("hT1_tdcTime", "T.hms.hT1_tdcTime")
                 .Define("hT2_tdcTime", "T.hms.hT2_tdcTime")
                 .Define("hT3_tdcTime", "T.hms.hT3_tdcTime")
                 .Define("hTRIG1_tdcTime", "T.hms.hTRIG1_tdcTime")
                 .Define("hTRIG2_tdcTime", "T.hms.hTRIG2_tdcTime")
                 .Define("hTRIG3_tdcTime", "T.hms.hTRIG3_tdcTime")
                 .Define("hTRIG4_tdcTime", "T.hms.hTRIG4_tdcTime")
                 .Define("hTRIG5_tdcTime", "T.hms.hTRIG5_tdcTime")
                 .Define("hTRIG6_tdcTime", "T.hms.hTRIG6_tdcTime")
                 //.Define("beta", "beta")
                // .Define("caletracknorm", "caletracknorm")
                 //.Define("caletottracknorm", "caletottracknorm")
                // .Define("caletotnorm", "caletotnorm")
                 .Define("vx", "H.react.x")
                 .Define("vy", "H.react.y")
                 .Define("vz", "H.react.z")
                 .Define("dp", "H.gtr.dp")
                 .Define("th", "H.gtr.th")
                 .Define("ph", "H.gtr.ph")
                 .Define("px", "H.gtr.px")
                 .Define("py", "H.gtr.py")
                 .Define("pz", "H.gtr.pz")
                 .Define("cernpe", "H.cer.npeSum")
                 .Define("caltracknorm", "H.cal.etottracknorm")
                 //.Define("evt", "rdfentry_") // define event number from row of rdataframe for troubleshoooting (threadsafe)
                 .Define("evt", "g.evnum"); // define event number from row of rdataframe for troubleshoooting (threadsafe)
                 //.Filter("TMath::Abs(H.gtr.th) < 0.08")
                 //.Filter("TMath::Abs(H.gtr.ph) < 0.04")
                 //.Filter("TMath::Abs(H.gtr.dp) < 10");

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

  /////TODO: convert all needed treeout histos to df.Histo1D() and add to dataframe
  //TH1F *h1time = new TH1F("h1time", "pulse (>20mV) shift (4*ns units) relatively to elastic refwf (all found pulses included)", 200, -50, 50);
  //TH1F *h2time = new TH1F("h2time", "pulse (>20mV) time (ns) (all found pulses included)", 200, -100, 100);


  ROOT::RDF::TH1DModel m_h1time("h1time", "pulse (>20mV) shift (4*ns units) relatively to elastic refwf (all found pulses included)", 200, -50, 50);
  ROOT::RDF::TH1DModel m_h2time("h2time", "pulse (>20mV) time (ns) (all found pulses included)", 200, -100, 100);
 
 //Fill with finter[blocknum]->GetParameter(2 + 2 * p), 1.) and wfttime[blocknumn][numpulses] 
  // auto h1time = df2.Histo1D(m_h1time ,"");
// auto h2time = df2.Histo1D(m_h2time ,"");

  /////// ANALYZE //////////////////////////////////////////////////////////

  // this is where sequential event loop was

  // HMS CUTS
/*
  df2 = df2.Filter("TMath::Abs(H.gtr.th) < 0.08")
            .Filter("TMath::Abs(H.gtr.ph) < 0.04")
            .Filter("TMath::Abs(H.gtr.dp) < 10");

*/
  ////////Lambda funtion for the per-event wf analysis/////////
  auto analyze = [=](Int_t NSampWaveForm, const ROOT::VecOps::RVec<Double_t> &SampWaveForm, Double_t evt, Int_t NadcCounter, ROOT::VecOps::RVec<Double_t> &adcCounter, const ROOT::VecOps::RVec<Double_t> adcSampPulseTime, const ROOT::VecOps::RVec<Double_t> adcSampPulseTimeRaw, const ROOT::VecOps::RVec<Double_t> adcSampPulseAmp, const ROOT::VecOps::RVec<Double_t> adcSampPulseInt, const ROOT::VecOps::RVec<Double_t> adcSampPulsePed)
  {


    TStopwatch tlambda;
    tlambda.Start();

    //TF1 *finter[nblocks];
    std::vector<std::unique_ptr<TF1>> finter(nblocks); //object is now thread-local
    TH1F *hsig_i[nblocks];
    //Int_t pres[nblocks] = {0};
    std::vector<Int_t> pres(nblocks, 0);
    Double_t signal[nblocks][ntime];
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
    //Double_t wfampl[nblocks][maxwfpulses];
    //Double_t wftime[nblocks][maxwfpulses];
std::vector<std::vector<Double_t>> wfampl(nblocks,
  std::vector<Double_t>(maxwfpulses, -100));
std::vector<std::vector<Double_t>> wftime(nblocks,
  std::vector<Double_t>(maxwfpulses, -100));

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
           TStopwatch tfit;
           tfit.Start();
  // 1) Build a thread‐local Interpolator
      auto interpPtr = std::make_shared<ROOT::Math::Interpolator>(
        ntime,
        ROOT::Math::Interpolation::kCSPLINE
     );
    
       // copy your raw data in:
       interpPtr->SetData(ntime,
        interpX[bn].data(),
        interpY[bn].data());
    
    
    
          // 2) Clone the interpolator for this block
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

      /*
      for (Int_t p = 0; p < maxwfpulses; p++)
      {
        // cout << "Setting parameters for pulse " << p << " with index " << 2+2*p << " and " << 3+2*p << endl;
        wfampl[bn][p] = -1;
      }
*/

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
              //  wfnpulse[bn] = maxwfpulses - 1; // Prevent overflow
            }

            wfampl[bn][wfnpulse[bn]] = hsig_i[bn]->GetBinContent(it + 4); // get the amplitude of the pulse found

            wftime[bn][wfnpulse[bn]] = hsig_i[bn]->GetBinCenter(it + 4); // get the time of the pulse found

            // flag for the good pulse
            if (TMath::Abs(wftime[bn][wfnpulse[bn]] - timeref[bn] +(corr_time_HMS - cortime[bn]) / dt - timerefacc) < 4.)
            {
              good = 1;
             //cpulse = wfnpulse[bn];
            }
            cout<<"detection bloc "<<bn<<" pulse="<<wfnpulse[bn]<<" time="<<wftime[bn][wfnpulse[bn]]<<" (4ns) ampl="<<wfampl[bn][wfnpulse[bn]]<<" reftime= "<<timeref[bn]<<" (4ns) diff="<<TMath::Abs(wftime[bn][wfnpulse[bn]]-timeref[bn])<<" good="<<good<<endl;

            wfnpulse[bn]++;

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
        cout << "Warning : All samples < 2.0mV. Max Amp = "<<maxamp<<". Skipping " << bn << endl;
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

// ── STEP 1: Floor any zero‐error bins so χ² is finite ──
 /*
int nBins = hsig_i[bn]->GetNbinsX();
for (int ib = 1; ib <= nBins; ++ib) {
   if (hsig_i[bn]->GetBinError(ib) == 0)
   hsig_i[bn]->SetBinError(ib, 1.0);
}
   */
      if (wfnpulse[bn] > 0 && good == 1)
      {
        for (Int_t p = 0; p < TMath::Min(maxwfpulses, wfnpulse[bn]); p++)
        {
          finter[bn]->ReleaseParameter(2 + 2 * p);
          finter[bn]->ReleaseParameter(3 + 2 * p);

          //finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn]);
          //finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt);   
          finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p]- timeref[bn]);   
          finter[bn]->SetParameter(3 + 2 * p, wfampl[bn][p]);


          //finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] - 3, wftime[bn][p] - timeref[bn] + 3);
          //finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt - 3, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt + 3);
          finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p]- timeref[bn] - 3, wftime[bn][p]- timeref[bn] + 3);

          finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn][p] * 0.2, wfampl[bn][p] * 3);

        }
      }

      if (wfnpulse[bn] > 0 && good == 0)
      {
        for (Int_t p = 0; p < TMath::Min(maxwfpulses, wfnpulse[bn]); p++)
        {
          finter[bn]->ReleaseParameter(2 + 2 * p);
          finter[bn]->ReleaseParameter(3 + 2 * p);
          finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn]);
          // finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt);
          finter[bn]->SetParameter(3 + 2 * p, wfampl[bn][p]);
          finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] - 3, wftime[bn][p] - timeref[bn] + 3);
          //finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt - 3, wftime[bn][p] - timeref[bn] + (corr_time_HMS - cortime[bn]) / dt + 3);

          finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn][p] * 0.2, wfampl[bn][p] * 3);

          std::cout<<"SET PARAM "<<2 + 2 * p<<" : "<<finter[bn]->GetParameter(2 + 2 * p)<<" FROM: "<<wftime[bn][p]<<" "<<timeref[bn]<<" "<<corr_time_HMS<<" "<<cortime[bn]<<" "<<dt<<endl;

        }

        // On recherche quand meme un eventuel pulse en temps
        //cpulse = wfnpulse[bn];

        finter[bn]->ReleaseParameter(2 + 2 * wfnpulse[bn]);
        finter[bn]->ReleaseParameter(3 + 2 * wfnpulse[bn]);
        finter[bn]->SetParameter(2 + 2 * wfnpulse[bn], timerefacc - timerefacc);
        //finter[bn]->SetParameter(2 + 2 * wfnpulse[bn], timerefacc + (corr_time_HMS - cortime[bn]) / dt - timeref[bn]);
        finter[bn]->SetParameter(3 + 2 * wfnpulse[bn], 2);
        //finter[bn]->SetParLimits(2 + 2 * wfnpulse[bn], timerefacc - timerefacc - 4, timerefacc - timerefacc + 4);
        //finter[bn]->SetParLimits(2 + 2 * wfnpulse[bn], timerefacc + (corr_time_HMS - cortime[bn]) / dt - 4, timerefacc + (corr_time_HMS - cortime[bn]) / dt + 4);

        //finter[bn]->SetParLimits(3 + 2 * wfnpulse[bn], 0.05, 10);

        wfnpulse[bn]++;

      }

      if (wfnpulse[bn] == 0)
      {

        for (Int_t p = 0; p < 1; p++)
        {
          finter[bn]->ReleaseParameter(2 + 2 * p);
          finter[bn]->ReleaseParameter(3 + 2 * p);

          finter[bn]->SetParameter(2 + 2 * p, 0);
          //finter[bn]->SetParameter(2 + 2 * p, timerefacc + (corr_time_HMS - cortime[bn]) / dt);
          finter[bn]->SetParameter(3 + 2 * p, 2);

       //   finter[bn]->SetParLimits(2 + 2 * p, timerefacc - 4, timerefacc + 4);
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
      for (int p = wfnpulse[bn]; p < maxwfpulses; ++p) {
        finter[bn]->FixParameter(2 + 2*p, 0.);
        finter[bn]->FixParameter(3 + 2*p, 0.);
      }

      // hsig_i[bn]->Fit(Form("finter_%d",bn),"Q","",TMath::Max(0.,wftime[bn][0]-20),ntime);
      // File << "Checkpoint: Before fitting function for bn = " << bn << endl;

      //hsig_i[bn]->Fit(Form("finter_%d", bn), "Q", "", 0., ntime);
             
// 3) Prepare binned data from the histogram:
int nBins = hsig_i[bn]->GetNbinsX();
ROOT::Fit::BinData data(nBins, 1);
for (int ib = 1; ib <= nBins; ++ib) {
    double x[1] = { hsig_i[bn]->GetBinCenter(ib) };
    double y    = hsig_i[bn]->GetBinContent(ib);
    double err  = hsig_i[bn]->GetBinError(ib);
   //if(bn==898) cout<<"BIN ERRER: "<<err<<endl;
    data.Add(x, y, err);
}




for(int p=0; p<26;p++){
  //std::cout<<p<<" : "<<finter[bn]->GetParameter(p)<<endl;
//  std::cout<<p<<" : "<<finter[bn]->GetParameter(2 + 2*p)<<endl;
 // std::cout<<p<<" : "<<finter[bn]->GetParameter(3 + 2*p)<<endl;
}




// 4) Wrap TF1 into a IModelFunction via WrappedMultiTF1:
ROOT::Math::WrappedMultiTF1 wfunc(*finter[bn], finter[bn]->GetNdim());
// 5) Configure fitter
ROOT::Fit::Fitter fitter;
auto &cfg = fitter.Config();
//fitter.Config().SetMinimizer("Minuit2", "Migrad");
cfg.SetMinimizer("Minuit2", "Migrad");
//auto &mopts = fitter.Config().MinimizerOptions();s
auto &mopts = cfg.MinimizerOptions();
mopts.SetStrategy(0);
mopts.SetPrintLevel(0);
mopts.SetMaxIterations(1000);
fitter.SetFunction(wfunc, /*useGradient=*/false);


// tell the Fitter about your TF1’s parameters
cfg.CreateParamsSettings(wfunc);

int nFloat = 2 + 2*wfnpulse[bn];
// now fix & release exactly the slots you want:
//cfg.FixParameter(0);


for (int ip = nFloat; ip < nbparameters; ++ip)
cfg.ParSettings(ip).Fix();  
for (int ip = 1; ip < nFloat; ++ip)
cfg.ParSettings(ip).Release();




if(good == 0 && wfnpulse[bn]>1){
  cout<<"SET LIMITS HERE "<<wfnpulse[bn]<< endl;  
  cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).SetLowerLimit(0.05);
  cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).SetUpperLimit(10);
  cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).Fix();        // pin it at your seed = 0.05
   finter[bn]->SetParameter(3 + 2*(wfnpulse[bn]-1), 0.05);
  }


cfg.ParSettings(0).Fix();  
// lock the fake‐pulse amplitude to never go below 0.05:



//cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).Fix();
//cfg.ParSettings(2 + 2*(wfnpulse[bn]-1)).Fix();
//cfg.UseDataErrors(true);
//cfg.UseFunctionErrors(false);


if(bn==selectblock){
for(int p=0; p<26;p++){
  std::cout<<"BEOFRE FIT "<<p<<" : "<<finter[bn]->GetParameter(p)<<endl;
}
}


int refParIdx = 3 + 2*(wfnpulse[bn]-1);
std::cout
  << " about to fit ref amp par["<<refParIdx<<"] = "
  << finter[bn]->GetParameter(refParIdx)
  << "  bounds = ["
  << cfg.ParSettings(refParIdx).LowerLimit() << ","
  << cfg.ParSettings(refParIdx).UpperLimit() << "]\n";


// 6) run the fit
bool ok = fitter.LeastSquareFit(data);
if (!ok) {
  cout<<"FAILED ONCE"<<endl;
}

if(good == 0 && wfnpulse[bn]>1 && ok){
  //if a real pulse is found save the fit result first before releasing the ref pulse params
  for (Int_t p = 0; p < wfnpulse[bn]; ++p) {
      wfampl[bn][p] = finter[bn]->GetParameter(3 + 2*p);
      wftime[bn][p]  = finter[bn]->GetParameter(2 + 2*p)*dt                     // convert bins → ns
               + corr_time_HMS                // add HMS correction
               - cortime[bn]                  // subtract block‐by‐block cable delay
               - timerefacc*dt;               // subtract your reference‐time offset
}

  cfg.ParSettings(3 + 2*(wfnpulse[bn]-1)).Release();  
  ok = fitter.LeastSquareFit(data);
  }

/*
for(int p=0; p<26;p++){
  //std::cout<<p<<" : "<<finter[bn]->GetParameter(p)<<endl;
  std::cout<<p<<" : "<<finter[bn]->GetParameter(2 + 2*p)<<endl;
  std::cout<<p<<" : "<<finter[bn]->GetParameter(3 + 2*p)<<endl;
}

*/

if (!ok) {
  //std::cerr<<"Fit failed for event "<<evt<<", block "<<bn<<"\n";
  for(int p=0; p<26;p++){
    //std::cout<<p<<" : "<<finter[bn]->GetParameter(p)<<endl;
  //  std::cout<<p<<" : "<<finter[bn]->GetParameter(2 + 2*p)<<endl;
  //  std::cout<<p<<" : "<<finter[bn]->GetParameter(3 + 2*p)<<endl;
  }
  //return;
}
if (!ok) {
  // retry just this one with a tougher configuration
  cout<<"FAILED Again, retry"<<endl;
  cfg.MinimizerOptions().SetStrategy(2);
  mopts.SetPrintLevel(1);
  cfg.MinimizerOptions().SetMaxIterations(5000);
  ok = fitter.LeastSquareFit(data);
}
if (!ok) {
  std::cerr<<"STILL Fit failed for event "<<evt<<", block "<<bn<<"\n";
 failcount++;
 if(good == 0 && wfnpulse[bn]>1) return;
  for (int p = 0; p < wfnpulse[bn]; p++)
  {
  wftime[bn][p] = -1000;
  wfampl[bn][p] = -1000;
  }
  
  for(int p=0; p<26;p++){
    //std::cout<<p<<" : "<<finter[bn]->GetParameter(p)<<endl;
   // std::cout<<p<<" : "<<finter[bn]->GetParameter(2 + 2*p)<<endl;
   // std::cout<<p<<" : "<<finter[bn]->GetParameter(3 + 2*p)<<endl;
   // std::cout<<wfampl[bn][p]<<"  "<<wftime[bn][p]<<endl;
  }
  
  return;
}

// 7) Extract parameters to your arrays:
auto &result = fitter.Result();

/*
    // Extract the fit parameters directly into your arrays: (replaces line 770/772)
    for (Int_t p = 0; p < wfnpulse[bn]; ++p) {
     // wftime[bn][p] = result.Parameter(2 + 2*p) * dt + corr_time_HMS - cortime[bn] - timerefacc * dt;
    wftime[bn][p] = result.Parameter(2 + 2*p) * dt + corr_time_HMS - cortime[bn] + timerefacc * dt;
     // wftime[bn][p] = result.Parameter(2 + 2*p);

      wfampl[bn][p] = result.Parameter(3 + 2*p);
  }
*/

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
    wfampl[bn][p] = result.Parameter(3 + 2*p);
    wftime[bn][p]  = binOff*dt                     // convert bins → ns
               + corr_time_HMS                // add HMS correction
               - cortime[bn]                  // subtract block‐by‐block cable delay
               - timerefacc*dt;               // subtract your reference‐time offset


    std::cout
      << "AFTER FIT: bloc " << bn
      << "  pulse="   << p
      << "  binOff="  << binOff
      << "  binIdx="  << binIndex
      << "  time_unc="<< uncorTime  << " ns"
      << "  time_cor="<< corrTime   << " ns"
      << "  ampl="    << wfampl[bn][p]       << " mV"
       << " wftime="    << wftime[bn][p]       << " mV"
      << std::endl;
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
     if(bn==selectblock){
      for(int p=0; p<26;p++){
        std::cout<<p<<" : "<<result.Parameter(p)<<endl;
    }
  }

      for(int np=0;np<wfnpulse[bn];np++){
        std::cout<<wfampl[bn][np]<<" "<<ntime<<endl;
    }
    cout <<" Entry = "<< evt <<" block: "<<bn<<"  fit time="<<tfit.RealTime()<<endl;
    tfit.Stop();
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

        //finter[i] = new TF1(Form("finter_%d", i), func, -200, ntime + 200, nbparameters);
        //finter[i] = std::make_unique<TF1>(Form("finter_%d", i),func, -200, ntime + 200,nbparameters); //now a thread-local object
        //finter[i]->FixParameter(0, i); // numero de bloc
        //finter[i]->SetLineColor(4);
        //finter[i]->SetNpx(1100);
/*
        for (Int_t p = 0; p < maxwfpulses; p++)
        {
          wftime[i][p] = -100;
          wfampl[i][p] = -100;
        }
*/
        hsig_i[i] = new TH1F(Form("hsig_i%d", i), Form("hsig_i%d", i), ntime, 0, ntime);
        hsig_i[i]->SetLineColor(1);
        hsig_i[i]->SetLineWidth(2);
        hsig_i[i]->GetXaxis()->SetTitle("Time (4 ns)");
        hsig_i[i]->GetYaxis()->SetTitle("(mV)");
        hsig_i[i]->GetYaxis()->SetLabelSize(0.05);

        for (Int_t j = 0; j < ntime; j++)
        {
          signal[i][j] = 0;
        }

      } // liste de presence des blocs! pres=0 si bloc absent, pres=1 s'il est present

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
                signal[bloc][it] = SampWaveForm[ns];
                hsig_i[bloc]->SetBinContent(it + 1, signal[bloc][it]);
              }
              ns++;
            }
         // }
        } // fin if(bloc number is good)
      } // fin while()

      // Read the hcana calculated variables

      // if(!(NadcCounter==NadcSampPulseAmp&&NadcSampPulseInt==NadcSampPulsePed&&NadcSampPulseInt==NadcSampPulseTime)){cout<<"!!!!! Problem Ndata !!!!!! "<<evt<<endl;ndataprob++;}

      cout << "current event = " << evt << "   NadcCounter=" << NadcCounter << endl;
      for (Int_t iNdata = 0; iNdata < NadcCounter; iNdata++)
      {
        if (adcCounter[iNdata] == 2000)
          adcCounter[iNdata] = 1080; // nouveau numero du scintillateur attribue par Malek
        if (adcCounter[iNdata] == 2001)
          adcCounter[iNdata] = 1081; // nouveau numero du scintillateur attribue par Malek
        cout<<"NDATA "<<iNdata<<endl;
        // determine HMS time correction
        if (iNdata == 0)
        {
          corr_time_HMS = adcSampPulseTime[iNdata] - (adcSampPulseTimeRaw[iNdata] / 16.) - tdcoffset[(int)(adcCounter[iNdata])];
          cout << "HMS time correction event " << evt <<"  "<< corr_time_HMS <<" "<<tdcoffset[(int)(adcCounter[iNdata])]<< endl;
        }
        if (iNdata != 0)
        {
          if (TMath::Abs(corr_time_HMS - (adcSampPulseTime[iNdata] - adcSampPulseTimeRaw[iNdata] / 16. - tdcoffset[(int)(adcCounter[iNdata])])) > 0.001)
          {
            cout << "problem HMS time correction event " << evt <<"  "<< corr_time_HMS <<" "<<tdcoffset[(int)(adcCounter[iNdata])]<< endl;
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

          cout << "BEFORE FIT HMS time correction event " << evt <<"  "<< corr_time_HMS << endl;
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
         // chi2[i] = finter[i]->GetChisquare() / finter[i]->GetNDF();
          cout<<"CHI2: "<<chi2[i]<<endl;

          for (Int_t p = 0; p < TMath::Max(wfnpulse[i], 1); p++)
          {

           // wftime[i][p] = finter[i]->GetParameter(2 + 2 * p) * dt + corr_time_HMS; // temps du pulse en ns
           // wftime[i][p] = finter[i]->GetParameter(2 + 2 * p) * dt + corr_time_HMS - cortime[i] - timerefacc * dt;

           // wfampl[i][p] = finter[i]->GetParameter(3 + 2 * p); // amplitude du pulse en ns

            if (wfampl[i][p] > 20)
            {

              /// fix later //// diagnostic histos
              
               // h2time->Fill(wftime[i][p], 1.);
               //h2time[i]=wftime[i][p];
               h2time.push_back(wftime[i][p]);

               // h1time->Fill(finter[i]->GetParameter(2 + 2 * p), 1.);
               // h1time[i]=finter[i]->GetParameter(2 + 2 * p);

               h1time.push_back(finter[i]->GetParameter(2 + 2 * p) - timerefacc + corr_time_HMS / dt);
               //h1time.push_back(finter[i]->GetParameter(2 + 2 * p));
              
            } // Fill the time spectrums

            if (p == 0)
            {
              timewf[i] = wftime[i][p];
              amplwf[i] = wfampl[i][p];
            }

            // if(p>0){if(TMath::Abs(wftime[i][p]-timerefacc/dt)<TMath::Abs(timewf[i]-timerefacc/dt)){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}}//pour prendre le pulse dont le temps est le plus proche de timerefacc

            // New modification based on the production runs

            if (p > 0)
            {
              if (TMath::Abs(wftime[i][p]) < TMath::Abs(timewf[i]))
              {
                timewf[i] = wftime[i][p];
                amplwf[i] = wfampl[i][p];
              }
            } // pour prendre le pulse dont le temps est le plus proche de timerefacc
            //cout<<"AFTER FIT: bloc "<<i<<" pulse="<<p<<" time="<<wftime[i][p]<<" (4ns) ampl="<<wfampl[i][p]<<" reftime= "<<timeref[i]<< " time uncorrected=" << (wftime[i][p] - corr_time_HMS + cortime[i] + timerefacc * dt) / dt + timeref[i] << "(4ns)"<<endl;

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

          integ[i] += signal[i][it];
          integtot += signal[i][it];

          if (it > binmin && it < binmax)
          { // cosmic pulse window

            ener[i] += signal[i][it];
            enertot += signal[i][it];
          }

          if (!(it > binmin && it < binmax))
          {
            bkg[i] += signal[i][it];
          } // background window

          if (signal[i][it] > sigmax[i])
          {
            time[i] = it;
            sigmax[i] = signal[i][it];
            ampl[i] = signal[i][it];
          } // determine the pulse maximum
        }

        // doesnt seem to be used??
        ener[i] -= bkg[i] * (binmax - binmin - 1) / (nsamp - (binmax - binmin - 1)); // subtract the bkg contribution (we have to normalize since the window widths are not the same)

        bkg[i] = bkg[i] / (nsamp - (binmax - binmin - 1)); // mean value of bkg

        for (Int_t it = 0; it < nsamp; it++)
        {

          if (!(it > binmin && it < binmax))
          {
            noise[i] += (signal[i][it] - bkg[i]) * (signal[i][it] - bkg[i]) / (nsamp - (binmax - binmin - 1));
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
          if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.5)
          {
            max50 = it;
          }
          if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.1)
          {
            max90 = it;
          }
        }
        for (Int_t it = time[i]; it > -0.5; it--)
        { // aller vers la gauche du maximum
          if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.5)
          {
            min50 = it;
          }
          if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.1)
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
    
    int nonzero =0;

    //  TFile *fout = new TFile(Form("/volatile/hallc/nps/kerver/ROOTfiles/WF/waveforms/debug_run_%d_evt_%lld.root",run, (long long)evt), "RECREATE");
      for (int i = 0; i < nblocks; ++i) {
        //if (hsig_i[i] && hsig_i[i]->GetEntries() > 0) nonzero++;
       // if (hsig_i[i]){
        //  std::cout << "evt=" << evt << ", block=" << i
       //   << ", integral=" << hsig_i[i]->Integral() << std::endl;
       //  hsig_i[i]->Write(Form("hsig_block_%d", i));


         if(i==selectblock){
         TCanvas *c2 = new TCanvas("c2", "c2", 1000, 800);
         hsig_i[i]->SetDirectory(0);  
         c2->cd();  
         c2->Update();
         hsig_i[i]->Draw();
         finter[i]->Draw("same");
         for(int p=0; p<26;p++){
          std::cout<<p<<" : "<<finter[i]->GetParameter(p)<<endl;
      }

         c2->Update();
         TString outName = Form("hist_fit_run%d_evt%lld.png", run, (long long)evt);
         c2->SaveAs(outName);
             for (Int_t p = 0; p < TMath::Max(wfnpulse[i], 1); p++)
             {
   
                 cout << "pulse n:" << p
                      << " time initial=" << (wftime[i][p] - corr_time_HMS + cortime[i] + timerefacc * dt) / dt + timeref[i] << "(4ns)"
                      << " time corrected=" << wftime[i][p] << "(ns)"
                      << " ampl=" << wfampl[i][p] << " (mV)" << endl;
             }
             //break;
   
         }
        //}


        }

       // fout->Close();
     // delete fout;      
      
      
    
      
    //end of diagnostic snippet 


    for (int i = 0; i < nblocks; ++i) {
   //   delete hsig_i[i];
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
  df_final = df_final.Define("chi2", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
   std::vector<Double_t> output = std::get<0>(tuple);
   return output; },{"tuple"});
  df_final = df_final.Define("ampl", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
    std::vector<Double_t> output = std::get<1>(tuple);
    return output; },{"tuple"});
  df_final = df_final.Define("amplwf", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
    std::vector<Double_t> output = std::get<2>(tuple);
    return output; },{"tuple"});
    df_final = df_final.Define("wfnpulse", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
      std::vector<Int_t> output = std::get<3>(tuple);
      return output; },{"tuple"});  
      df_final = df_final.Define("Sampampl", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
        std::vector<Double_t> output = std::get<4>(tuple);
        return output; },{"tuple"});
        df_final = df_final.Define("Samptime", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
          std::vector<Double_t> output = std::get<5>(tuple);
          return output; },{"tuple"});  
          df_final = df_final.Define("timewf", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
            std::vector<Double_t> output = std::get<6>(tuple);
            return output; },{"tuple"});  
            df_final = df_final.Define("enertot", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
              Double_t output = std::get<7>(tuple);
              return output; },{"tuple"});  
              df_final = df_final.Define("integtot", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
                Double_t output = std::get<8>(tuple);
                return output; },{"tuple"});  
                df_final = df_final.Define("pres", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
                  std::vector<Int_t> output = std::get<9>(tuple);
                  return output; },{"tuple"});  
                  df_final = df_final.Define("corr_time_HMS", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
                    Double_t output = std::get<10>(tuple);
                    return output; },{"tuple"});
                    df_final = df_final.Define("h1time", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
                      std::vector<Double_t> output = std::get<11>(tuple);
                      return output; },{"tuple"}); 
                      df_final = df_final.Define("h2time", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
                        std::vector<Double_t> output = std::get<12>(tuple);
                        return output; },{"tuple"}); 
                        df_final = df_final.Define("wfampl", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
                          std::vector<std::vector<Double_t>> output = std::get<13>(tuple);
                          return output; },{"tuple"});  
                          df_final = df_final.Define("wftime", [](const std::tuple<std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Int_t>,std::vector<Double_t>,std::vector<Double_t>,std::vector<Double_t>,Double_t,Double_t,std::vector<Int_t>,Double_t,std::vector<Double_t>,std::vector<Double_t>,std::vector<std::vector<Double_t>>,std::vector<std::vector<Double_t>> > &tuple){ 
                            std::vector<std::vector<Double_t>> output = std::get<14>(tuple);
                            return output; },{"tuple"});    
                                                  
                    
 auto h_h1time = df_final.Histo1D(m_h1time ,"h1time");
 auto h_h2time = df_final.Histo1D(m_h2time ,"h2time");

  // Save the dataframe to the output ROOT file
  TString rootfilePath = Form("/volatile/hallc/nps/kerver/ROOTfiles/WF/waveforms/nps_production_%d_%d_%d_interactive_allwf.root", run, seg,nthreads);
  auto columnNames = df_final.GetColumnNames();
  for (const auto &col : columnNames) {
   // std::cout << col << std::endl;
}

df_final.Snapshot("treeout", rootfilePath, {"chi2","ampl","amplwf","wfnpulse","Sampampl","Samptime","timewf","enertot","integtot","pres","corr_time_HMS","h1time","h2time","evt","wfampl","wftime"});
/*
    // Open output file and save histogram
    TFile outFile(rootfilePath, "UPDATE");
    h_h1time->Write();  // Write histogram to the output file
    h_h2time->Write();  // Write histogram to the output file
    outFile.Close();
*/
  // Then safely write and delete at the end
  /*
    fout->cd();
    h1time->Write("h1time");
    h2time->Write("h2time");
    delete h1time;
    delete h2time;
    fout->Write();
    fout->Close();
  */
  // TTree was never written to output file???

  cout << "fin de l'analyse" << endl;
  t.Stop();
  t.Print();
}
