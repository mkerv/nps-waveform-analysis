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
#include "TSpectrum.h"
#include "TLine.h"
#include "TError.h"
#include <thread>
#include <mutex>
template <typename T>
void checkType(const T &obj)
{
  static_assert(std::is_same<T, std::vector<Double_t>>::value, "Type is not correct!");
}
static std::mutex s_spectrum_mutex;

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
std::atomic<int> nFitSucceeds{0};

// Matched Filter constants
const int mfleft = 5;                     // start bin of the refwf used as a filter = max of the refwf - mfleft
const int mfright = 5;                    // end bin of the refwf used as a filter = max of the refwf + mfright
const int mfwidth = mfleft + mfright + 1; // width of the refwf used as a filter (4ns units)
const int mfstart = 10;                   // Search for peaks in the wf between bins mfstart and mfend (4ns units)
const int mfend = 100;                    // Search for peaks in the wf between bins mfstart and mfend (4ns units)
const double specthres = 0.02;            // peaks with amplitude less than specthres*highest_peak are discarded when TSpectrum::Search() is called
const double mfthres = 1.5;               // peaks with amplitude less than mfthres in the mf of the wf are discarded
const double trig_thres = 10;             // threshold (mV) on sum of 3x3 wf to allow the fit of the central wf
const int coinc_width = 20;               // (4ns units) the trig_thres will be applied on sum of 3x3 wf in the region {expected coinc time +/- coinc_width} instead of all wf range
Double_t mfyref[nblocks][mfwidth];        // part of the refwf that will be used in the MF
Double_t mfint[nblocks];

Double_t timeref[nblocks];
Float_t cortime[nblocks];
Int_t preswf[nblocks];
// Double_t timerefacc = -4.5; // car les bons pulses elastiques poussent vers 45.5 (4*ns) dans la wf alors que ceux des production runs sont vers 35.5 (4ns)
Double_t timerefacc = 0; // en (-4ns) car le pic de timewf pousse vers -22ns

// ROOT::Math::Interpolator *interpolation[nblocks]; // function used for the interpolation
std::vector<std::vector<double>> interpX(nblocks),
    interpY(nblocks);

///////////////////////////////////////////////////////////// FIT FUNCTION USED AND SCHEMATICS /////////////////////////////////////////////////////////////
bool FastCloneAndFilter(const TString &inName,
                        const TString &outName)
{
  TFile *fin = TFile::Open(inName, "READ");
  if (!fin || fin->IsZombie())
    return false;
  TFile *fout = TFile::Open(outName, "RECREATE");
  if (!fout || fout->IsZombie())
  {
    fin->Close();
    return false;
  }

  // copy all non-T keys
  for (auto k : *fin->GetListOfKeys())
  {
    TKey *key = static_cast<TKey *>(k);
    TString name(key->GetName());
    if (name == "T")
      continue;
    // if (key->GetName() == "T") continue;
    fout->cd();
    key->ReadObj()->Write();
  }
  // clone T without the big branch
  TTree *tin = (TTree *)fin->Get("T");
  tin->SetBranchStatus("NPS.cal.fly.adcSampWaveform", 0);
  TTree *tout = tin->CloneTree(-1, "fast");
  fout->cd();
  tout->Write("T");

  fout->Close();
  fin->Close();
  return true;
}
///////////////////////////////////////////////////////////// Peak Finding Function /////////////////////////////////////////////////////////////
void FindPulsesMF(int bn,
                  const Double_t fullSigArr[],    // length = nblocks * ntime
                  const std::vector<Int_t> &pres, // size = nblocks
                  Double_t minsignal,
                  const Double_t mfyref[][mfwidth], // [nblocks][mfwidth]
                  const Double_t mfint[],           // length = nblocks
                  Double_t timeref_bin,             // in sample‐units
                  Double_t timerefacc,              // in sample‐units
                  Int_t &wfnpulse_out,              // output: how many pulses found
                  std::vector<Double_t>& wftime_out, // size = nblocks*maxwfpulses
                  std::vector<Double_t>& wfampl_out,  // size = nblocks*maxwfpulses
                  Int_t offset
)
{

  if (pres[bn] == 0)
  {
    wfnpulse_out = 0;
    return;
  }

  // 2) Build the matched‐filter array mfVals[0..ntime−1]:
  std::array<Double_t, ntime> mfVals;
  mfVals.fill(0.0);
  Double_t mfmin = 1.0e6;

  // Convolution: for it in [mfleft .. ntime−mfright−1],
  //   mfVals[it] = sum_{jt=0..mfwidth−1} ( (sigArr[it + jt − mfright] − minsignal) * reversed_kernel )
  // where reversed_kernel = mfyref[bn][mfwidth−1−jt], normalized by mfint[bn].
  for (int it = mfleft; it < ntime - mfright; ++it)
  {
    Double_t acc = 0.0;
    for (int jt = 0; jt < mfwidth; ++jt)
    {
      Double_t raw = fullSigArr[bn * ntime + (it + jt - mfright)];
      Double_t delta = raw - minsignal;
      Double_t kern = mfyref[bn][mfwidth - 1 - jt];
      acc += (delta * kern) / mfint[bn];
    }
    mfVals[it] = acc;
    if (acc < mfmin)
      mfmin = acc;
  }
  // Subtract the minimum so mfVals ⩾ 0:
  for (int it = mfleft; it < ntime - mfright; ++it)
  {
    mfVals[it] -= mfmin;
  }

  TString hname = Form("hMF_blk%d_thr%ld", bn, (long)std::hash<std::thread::id>{}(std::this_thread::get_id()));
  std::unique_ptr<TH1F> hMF(new TH1F(hname, hname, ntime, 0, ntime));
  for (int i = 0; i < ntime; ++i)
  {
    // why did we shift by1
    hMF->SetBinContent(i + 1, mfVals[i]);
  }

  // — suppress ROOT warnings for just this call: all the tspectrum finding too many pulses on noise
  //auto old = gErrorIgnoreLevel;
 // gErrorIgnoreLevel = kError; // ignore kInfo and kWarning
  Int_t npeaks;
    // lock _only_ the Search() itself
    std::lock_guard<std::mutex> lock(s_spectrum_mutex);
    TSpectrum spec(maxwfpulses);
    npeaks = spec.Search(hMF.get(), 2, "nobackground,nodraw", specthres);

 // gErrorIgnoreLevel = old;
  // Loop over TSpectrum’s found peaks, record those within [mfstart..mfend], amplitude > mfthres:
  for (int ip = 0; ip < npeaks && wfnpulse_out < maxwfpulses; ++ip)
  {
    Double_t xpos = spec.GetPositionX()[ip] - 2.0; // shift by 2 bins, just like the single‐threaded code
    Double_t ypos = spec.GetPositionY()[ip];
    if (xpos > std::max(mfstart, 0) && xpos < std::min(mfend, ntime - 1) && ypos > mfthres)
    {
      int ti = static_cast<int>(std::round(xpos));
      // raw amplitude = | raw_waveform(ti) − minsignal |
      Double_t rawAmp = std::abs(fullSigArr[bn * ntime + ti] - minsignal);
      wftime_out[offset+wfnpulse_out]=xpos;
      wfampl_out[offset+wfnpulse_out]=rawAmp;
      //wftime_out.push_back(xpos);
      //wfampl_out.push_back(rawAmp);
      wfnpulse_out++;
    }
  }

  if (wfnpulse_out > maxwfpulses - 2)
  {
    std::cout << "Warning: high number of pulses in block " << bn
              << " (wfnpulse=" << wfnpulse_out << ")\n";
  }

  return;
}
/////////////////////////////////////////////////////////////Function to check cluster energy threshold///////////////////////////////////////////////////////////////////////
bool PassClusterThreshold(int bn,
                          const Double_t fullSigArr[],    // length = nblocks*ntime
                          const std::vector<Int_t> &pres, // size = nblocks
                          int ncol,
                          int nlin,
                          int nblocks,
                          int ntime,
                          Double_t timeref_bin, // in sample units
                          Double_t timerefacc,  // in sample units
                          Double_t trig_thres,  // mV threshold
                          int coinc_width)      // in sample units
{

  // Compute the “center time” around which we look for the 3×3 sum:
  Double_t center = timeref_bin + timerefacc;

  int row = bn / ncol;
  int col = bn % ncol;

  // For each time‐bin it ∈ [0..ntime−1], sum this block + all 8 neighbors (if present).
  Double_t globalMin = 1e6;
  Double_t maxInWindow = -1e6;

  for (int it = 0; it < ntime; ++it)
  {
    // Start with the raw sample from block `bn`:
    Double_t sum3x3 = fullSigArr[bn * ntime + it];

    // Look at all eight offsets (Δrow, Δcol):
    static const int dR[8] = {0, 0, +1, -1, +1, +1, -1, -1};
    static const int dC[8] = {+1, -1, 0, 0, +1, -1, +1, -1};

    for (int k = 0; k < 8; ++k)
    {
      int nr = row + dR[k];
      int nc = col + dC[k];
      if (nr < 0 || nr >= nlin || nc < 0 || nc >= ncol)
        continue;
      int nb = nr * ncol + nc;
      if (pres[nb] == 1)
      {
        sum3x3 += fullSigArr[nb * ntime + it];
      }
    }
    if (sum3x3 < globalMin)
    {
      globalMin = sum3x3;
    }
    // If this time‐bin is inside the “coincidence window,” track its max:
    if (std::abs(double(it) - center) < coinc_width)
    {
      if (sum3x3 > maxInWindow)
      {
        maxInWindow = sum3x3;
      }
    }
  }

  // 4) After scanning all it: if (maxInWindow - globalMin) > trig_thres → pass
  return ((maxInWindow - globalMin) > trig_thres);
}

/////////////////////////////////////////////////////////////MAIN FUNCTION ///////////////////////////////////////////////////////////////////////
void TEST_2(int run, int seg, int threads)
{
  TStopwatch t;
  t.Start();
  int nthreads = 6; // or any number
  nthreads = threads;

  // BUILD TCHAIN
  TChain chain("T");
  TString filename = Form("/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_%d_%d_1_-1.root", run, seg);
  // TString filename = Form("/mss/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", run, seg);
  // TString filename = Form("../nps_hms_coin_%d_%d_1_-1.root", run, seg);
  TFile *testOpen = TFile::Open(filename);
  if (!testOpen || testOpen->IsZombie())
  {
    std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
    return;
  }
  testOpen->Close();

  TString outFile = Form("/volatile/hallc/nps/kerver/ROOTfiles/WF/nps_production_%d_%d_%d_interactive_tspec.root", run, seg, nthreads);

  if (!FastCloneAndFilter(filename, outFile))
  {
    std::cerr << "Clone failed!\n";
    return;
  }
  cout << "I/O Copy time = " << t.RealTime() << endl;
  t.Continue();

  // ENABLE MT
  gErrorIgnoreLevel = kError;
  ROOT::EnableImplicitMT(nthreads);
  std::cout << "Implicit MT enabled: " << ROOT::IsImplicitMTEnabled() << "\n";
  std::cout << "Number of threads: " << ROOT::GetThreadPoolSize() << "\n";

  chain.Add(filename);
  chain.SetBranchStatus("*", 0);

  chain.SetBranchStatus("g.evnum", 1);
  chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampWaveform", 1);
  chain.SetBranchStatus("Ndata.NPS.cal.fly.adcCounter", 1);
  chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseAmp", 1);
  chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseInt", 1);
  chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPed", 1);
  chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTime", 1);
  chain.SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTimeRaw", 1);
  chain.SetBranchStatus("NPS.cal.fly.adcSampWaveform", 1);
  chain.SetBranchStatus("NPS.cal.fly.adcCounter", 1);
  chain.SetBranchStatus("NPS.cal.fly.adcSampPulseAmp", 1);
  chain.SetBranchStatus("NPS.cal.fly.adcSampPulseInt", 1);
  chain.SetBranchStatus("NPS.cal.fly.adcSampPed", 1);
  chain.SetBranchStatus("NPS.cal.fly.adcSampPulseTime", 1);
  chain.SetBranchStatus("NPS.cal.fly.adcSampPulseTimeRaw", 1);

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
  // auto df1percent = df.Range(0, nEventsToProcess);
  // auto df1percent = df.Range(10, 50);

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

    ymax = 0;
    if (filewf.is_open())
    {
      filewf >> timeref[i] >> dum1;
      interpX[i].resize(ntime);
      interpY[i].resize(ntime);
      for (int it = 0; it < ntime; ++it)
      {
        filewf >> interpX[i][it] >> interpY[i][it];

        if (interpY[i][it] > ymax)
        {
          ymax = interpY[i][it];
          timeref[i] = interpX[i][it];
        }
      }
      mfint[i] = 0;
      for (Int_t it = 0; it < ntime; it++)
      {
        if (TMath::Abs(timeref[i] - interpX[i][it]) < 0.001)
        {
          for (Int_t jt = 0; jt < mfwidth; jt++)
          {
            mfyref[i][jt] = interpY[i][it + jt - mfleft];
            mfint[i] += mfyref[i][jt];
          }
        }
      }
      preswf[i] = 1;
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
    if (cortime[i] == 0)
    {
      cortime[i] = -0.0000001;
    }
  }
  filetime.close();

  // Load only required branches into dataframe
  auto df2 = df.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")
  //              auto df2 = df1percent.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")
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
  else if (run > 5523 && run < 7013)
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

  /////// ANALYZE //////////////////////////////////////////////////////////

  // this is where sequential event loop was
  ////////Lambda funtion for the per-event wf analysis/////////
  auto analyze = [=](Int_t NSampWaveForm, const ROOT::VecOps::RVec<Double_t> &SampWaveForm, Double_t evt, Int_t NadcCounter, ROOT::VecOps::RVec<Double_t> &adcCounter, const ROOT::VecOps::RVec<Double_t> adcSampPulseTime, const ROOT::VecOps::RVec<Double_t> adcSampPulseTimeRaw, const ROOT::VecOps::RVec<Double_t> adcSampPulseAmp, const ROOT::VecOps::RVec<Double_t> adcSampPulseInt, const ROOT::VecOps::RVec<Double_t> adcSampPulsePed)
  {
    bool fitFailed = false;

    TStopwatch tlambda;
    tlambda.Start();

    std::vector<std::unique_ptr<TF1>> finter(nblocks); // object is now thread-local
    std::vector<Int_t> pres(nblocks, 0);
    std::vector<Double_t> signal(nblocks * ntime);
    std::vector<Double_t> minsignal(nblocks, 1e6);
    std::array<Double_t, ntime> Err{};

    Int_t bloc, nsamp;

    Double_t enertot = 0.;
    Double_t integtot = 0.;
    Double_t corr_time_HMS = 0.;

    std::vector<Double_t> timewf(nblocks, -100);
    std::vector<Double_t> amplwf(nblocks, -100);
    std::vector<Double_t> chi2(nblocks, -100.0);
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

    std::vector<Int_t> wfnpulse(nblocks, 0);

    //ROOT::RVec<Double_t> wfampl;
    //ROOT::RVec<Double_t> wftime;
    //wfampl.reserve(nblocks * maxwfpulses);
    //wftime.reserve(nblocks * maxwfpulses);
    std::vector<Double_t> wfampl, wftime;
wfampl.resize(nblocks*maxwfpulses, -999);
wftime.resize(nblocks*maxwfpulses, -999);
std::vector<Int_t> blockOffset(nblocks + 1, 0);
Int_t currentOffset = 0;


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
    auto Fitwf = [=, &blockOffset, &fitFailed, &signal, &wfnpulse, &wfampl, &finter, &wftime, &corr_time_HMS, &chi2](Double_t evt, Int_t bn, const Double_t Err_arr[])
    {

    // Dont fit pedestal with a constant function if wnpulse from tpsectrum==0  
    if (wfnpulse[bn] == 0) {
        chi2[bn] = -100.0;
        return;
    }


      // Build a thread‐local Interpolator
      auto interpPtr = std::make_shared<ROOT::Math::Interpolator>(
          ntime,
          ROOT::Math::Interpolation::kCSPLINE);

      // copy your raw data in:
      interpPtr->SetData(ntime,
                         interpX[bn].data(),
                         interpY[bn].data());

      auto func = [interpPtr, bn, &wfnpulse](Double_t *x, Double_t *par) mutable
      {
        Double_t val = par[0];
        for (Int_t p = 0; p < wfnpulse[bn]; ++p)
        {
          Double_t dt0 = x[0] - par[1 + 2 * p];
          //  Double_t raw = x[0] - par[2 + 2*p];
          //  Double_t dt0 = std::min(std::max(raw, 0.0), double(ntime-1));
          if (dt0 > 1 && dt0 < ntime - 1)
          {
            val += par[2 + 2 * p] * interpPtr->Eval(dt0);
          }
        }
        return val;
      };

      std::string fname = Form("finter_bn%d_ptr%p_evt%.0f", bn, (void *)finter[bn].get(), evt);
      finter[bn] = std::make_unique<TF1>(fname.c_str(), func, -200, ntime + 200, wfnpulse[bn] * 2 + 1);

      finter[bn]->SetNpx(1100);

      if (!finter[bn])
      {
        cout << " block=" << bn << " finter is nullptr" << endl;
        return;
      }

      /////////////////OLD WF PEAK FINDING WITH HSIG////////////

      if (wfnpulse[bn] > maxwfpulses - 2)
      {
        //cout << "Warning : excessively high number of pulses in the block wf " << bn << endl;
      }

      // Adjust the parameters of the fit function
      if (wfnpulse[bn] > 0)
      {
        for (Int_t p = 0; p < TMath::Min(maxwfpulses, wfnpulse[bn]); p++)
        {
          finter[bn]->ReleaseParameter(1 + 2 * p);
          finter[bn]->ReleaseParameter(2 + 2 * p);
          finter[bn]->SetParameter(1 + 2 * p, wftime[blockOffset[bn] + p] - timeref[bn]);
          finter[bn]->SetParameter(2 + 2 * p, wfampl[blockOffset[bn] + p]);
          finter[bn]->SetParLimits(1 + 2 * p, wftime[blockOffset[bn] + p] - timeref[bn] - 4, wftime[blockOffset[bn] + p] - timeref[bn] + 4);
          finter[bn]->SetParLimits(2 + 2 * p, wfampl[blockOffset[bn] + p] * 0.2, wfampl[blockOffset[bn] + p] * 5);
        }
      }

      finter[bn]->SetParameter(0, 0.);
      finter[bn]->SetParLimits(0, -100, 100.);
      double pedestal = 0;
      for (int i = 0; i < 20; ++i)
      {
        pedestal += signal[bn * ntime + i];
      }
      pedestal /= 20;
      finter[bn]->SetParameter(0, pedestal);

      // Prepare binned data from the histogram:
      ROOT::Fit::BinData data(ntime, /*nDim=*/1);
      for (int ib = mfstart; ib < mfend; ++ib)
      {
        double x[1] = {static_cast<double>(ib)};
        double y = signal[bn * ntime + ib];
        // double err  = std::sqrt(std::abs(y * 4.096 / 2.0)) / 4.096;
        double err = Err_arr[ib];
        data.Add(x, y, err);
      }

      // Wrap TF1 into a IModelFunction via WrappedMultiTF1:
      ROOT::Math::WrappedMultiTF1 wfunc(*finter[bn], finter[bn]->GetNdim());
      // Configure fitter
      ROOT::Fit::Fitter fitter;
      auto &cfg = fitter.Config();
      // fitter.Config().SetMinimizer("Minuit2", "Migrad");
      cfg.SetMinimizer("Minuit2", "Migrad");
      // cfg.SetMinimizer("Minuit2", "MINIMIZE");
      // cfg.SetMinimizer("Fumili","");
      // auto &mopts = fitter.Config().Mi  nimizerOptions();
      auto &mopts = cfg.MinimizerOptions();
      mopts.SetStrategy(1);
      mopts.SetPrintLevel(0);
      mopts.SetMaxIterations(1000);
      cfg.CreateParamsSettings(wfunc);

      /*
      for (int p = 0; p < wfnpulse[bn]; ++p) {
        int iA   = 2 + 2*p;
        double seed = wfampl[bn*maxwfpulses + p];
        double lo   = seed * 0.2;
        double hi   = seed * 5.0;
       // std::cout << "[debug] about to set param" << iA
       //           << " limits = [" << lo << "," << hi << "]\n";
        auto &ps = cfg.ParSettings(iA);
        //ps.SetLowerLimit(lo);
        //ps.SetUpperLimit(hi);
      ps.SetLimits(seed * 0.2,    // finite lower
                   seed * 5.0);   // finite upper
                   ps.Fix();


        // now immediately print what Minuit2 believes:
       // std::cout << "[debug] param" << iA <<" "<<ps.IsFixed()
        //          << " now has limits = ["
         //         << (ps.HasLowerLimit() ? std::to_string(ps.LowerLimit()) : "-inf")
          //        << ","
           //       << (ps.HasUpperLimit() ? std::to_string(ps.UpperLimit()) : "+inf")
            //      << "]\n";
      }
      */

      /*
      int npar2 = 1 + 2*wfnpulse[bn];

      for (int i = 0; i < npar2; ++i) {
        auto &ps = cfg.ParSettings(i);
       // std::cout<< "param"<<i
        //  <<"  value="<< ps.Value()
         // <<"  step=" << ps.StepSize()
         // <<"  limits=["<< (ps.HasLowerLimit()?std::to_string(ps.LowerLimit()):"-inf")
          //              <<","<< (ps.HasUpperLimit()?std::to_string(ps.UpperLimit()):"+inf") <<"]"
         // << std::endl;
      }
      */

      fitter.SetFunction(wfunc, false);

      /*
      cout<<"pre-fit params: ";
      for (int ip = 0; ip < finter[bn]->GetNpar(); ++ip) {
        cout<<finter[bn]->GetParameter(ip)<<" , ";
      }
      cout<<endl;
      */
      bool ok = fitter.LeastSquareFit(data);
      if (ok)
      {
        nFitSucceeds.fetch_add(1, std::memory_order_relaxed);
      }

      if (!ok)
      {
        // retry just this one with a tougher configuration
        // cout<<"FAILED Again, retry"<<endl;
        cfg.MinimizerOptions().SetStrategy(2);
        // mopts.SetPrintLevel(1);
        cfg.MinimizerOptions().SetMaxIterations(5000);
        ok = fitter.LeastSquareFit(data);
        if (ok)
        {
          nFitSucceeds.fetch_add(1, std::memory_order_relaxed);
        }
      }
      if (!ok)
      {
        // std::cerr<<std::fixed<<"STILL Fit failed for event "<<evt<<", block "<<bn<<"\n";
        fitFailed = true;
        nFitFailures.fetch_add(1, std::memory_order_relaxed);
        for (int p = 0; p < wfnpulse[bn]; p++)
        {
          // change time back to corrected in ns
          wftime[blockOffset[bn] + p] = (wftime[blockOffset[bn] + p] - timeref[bn]) * dt // convert bins → ns
                                        + corr_time_HMS                                  // add HMS correction
                                        - cortime[bn]                                    // subtract block‐by‐block cable delay
                                        - timerefacc * dt;

          // cout<< wfampl[blockOffset[bn] + p ]<<"  @ " <<wftime[blockOffset[bn] + p ] <<endl;
        }
        chi2[bn] = -100.;
        return;
      }

      // Extract parameters to arrays
      auto &result = fitter.Result();

      for (Int_t p = 0; p < wfnpulse[bn]; ++p)
      {
        // 1) the fitted bin shift (same units your x[] was in)
        double binOff = result.Parameter(1 + 2 * p);

        // 2) what bin *would* that be in the original histogram?
        double binIndex = binOff + timeref[bn]; // timeref[] is itself in bins

        // 3) uncorrected time in ns
        double uncorTime = binIndex * dt; // dt = 4 ns/bin

        // 4) fully corrected time in ns
        double corrTime = uncorTime + (corr_time_HMS - cortime[bn]);

        // 5) amplitude
        wfampl[blockOffset[bn] + p] = result.Parameter(2 + 2 * p);
        wftime[blockOffset[bn] + p] = binOff * dt        // convert bins → ns
                                      + corr_time_HMS    // add HMS correction
                                      - cortime[bn]      // subtract block‐by‐block cable delay
                                      - timerefacc * dt; // subtract your reference‐time offset
        // if(wfampl[blockOffset[bn]+p]<0.0)  cout<<"HERE!!!!!"<< wfampl[blockOffset[bn] + p ]<<"   +  "<<binOff<<endl;
      }
      unsigned npar = result.NPar();
      for (unsigned ip = 0; ip < npar; ++ip)
      {
        finter[bn]->SetParameter(ip, result.Parameter(ip));
      }

      double tempchi2 = result.Chi2(); // total χ² from the fit :contentReference[oaicite:0]{index=0}
      unsigned int ndf = result.Ndf(); // degrees of freedom :contentReference[oaicite:1]{index=1}

      chi2[bn] = tempchi2 / ndf;
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

      } // liste de presence des blocs! pres=0 si bloc absent, pres=1 s'il est present
      std::fill(signal.begin(), signal.end(), 0.);
      // Extract the data from the variable SampWaveForm[] (NPS.cal.fly.adcSampWaveform)

      int ns = 0; // ns represent for a given event the element number of the NPS.cal.fly.adcSampWaveform variable
      while (ns < NSampWaveForm)
      {
        bloc = SampWaveForm[ns];
        ns++; // bloc number (actually the slot number)
        nsamp = SampWaveForm[ns];
        ns++; // time samples (should be equal to ntime (110))

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

          for (Int_t it = 0; it < nsamp; it++)
          {
            if (bloc > -0.5 && bloc < nblocks)
            {
              signal[bloc * ntime + it] = SampWaveForm[ns];
              minsignal[bloc] = TMath::Min(minsignal[bloc], signal[bloc * ntime + it]);
            }
            ns++;
          }
        } // fin if(bloc number is good)
      } // fin while()

      // Read the hcana calculated variables

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
            // cout << "problem HMS time correction event " << evt <<"  "<< corr_time_HMS <<" "<<tdcoffset[(int)(adcCounter[iNdata])]<< endl;
          }
        }

        if (!(adcCounter[iNdata] >= 0 && adcCounter[iNdata] < nblocks + 2) /*&&adcCounter[iNdata]!=196*/)
        {
          cout << "****** Problem adcCounter ******* " << evt << " " << iNdata << " " << adcCounter[iNdata] << endl;
        }
        if (adcCounter[iNdata] >= 0 && adcCounter[iNdata] < nblocks)
        {
          Npulse[(int)(adcCounter[iNdata])] += 1;

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
        if (pres[i] == 1 && preswf[i] == 1) // if this block is present during event
        {
          for (Int_t it = 0; it < nsamp; it++)
          {
            Double_t y = signal[i * ntime + it];
            Double_t e = std::sqrt(std::abs(y * 4.096 / 2.)) / 4.096;

            if (e < 1.)
            {
              e = std::sqrt(std::abs(1.0 * 4.096 / 2.)) / 4.096;
            }
            Err[it] = e;
          }
          // Get the number of pulses from tspecrtum
         // blockOffset[i] = wftime.size();
          blockOffset[i] = currentOffset;
          FindPulsesMF(i, signal.data(), pres, minsignal[i], mfyref, mfint, timeref[i], timerefacc, wfnpulse[i], wftime, wfampl,blockOffset[i]);
          currentOffset += wfnpulse[i];
          bool okToFit = PassClusterThreshold(i, signal.data(), pres, ncol, nlin, nblocks, ntime, timeref[i], timerefacc, trig_thres, coinc_width);

          if (okToFit)
          {
            // cout<<evt<<" passed cluster threshold for block "<<i<<endl;
            Fitwf(evt, i, Err.data());
            if (finter[i])
            {
              finter[i]->SetLineColor(kBlue);
              finter[i]->SetLineWidth(2);
            }
            else
            {
              cout << "Failed to create TF1 in Fitwf for block " << i << ", event " << evt << endl;
            }
          }
          else
          {
            // We decided NOT to fit this block, so explicitly leave finter[i] == nullptr
            // and skip every attempt to use it.

            continue;
          }

          for (Int_t p = 0; p < wfnpulse[i]; p++)
          {

            if (wfampl[blockOffset[i] + p] > 20)
            {
              h2time.push_back(wftime[blockOffset[i] + p]);
              h1time.push_back(finter[i]->GetParameter(1 + 2 * p) - timerefacc + corr_time_HMS / dt);
              // h1time.push_back(finter[i]->GsetParameter(2 + 2 * p));

            } // Fill the time spectrums

            if (p == 0)
            {
              timewf[i] = wftime[blockOffset[i] + p];
              amplwf[i] = wfampl[blockOffset[i] + p];
            }

            // if(p>0){if(TMath::Abs(wftime[i][p]-timerefacc/dt)<TMath::Abs(timewf[i]-timerefacc/dt)){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}}//pour prendre le pulse dont le temps est le plus proche de timerefacc

            // New modification based on the production runs

            if (p > 0)
            {
              if (TMath::Abs(wftime[blockOffset[i] + p]) < TMath::Abs(timewf[i]))
              {
                timewf[i] = wftime[blockOffset[i] + p];
                amplwf[i] = wfampl[blockOffset[i] + p];
              }
            } // pour prendre le pulse dont le temps est le plus proche de timerefacc

          } // end of loop over maxpulses
        } // end of condition over col,lin
      } // end of loop over les blocs
      //blockOffset[nblocks] = wftime.size(); // after loop over blocks ensure the last offset is filled
      blockOffset[nblocks] = currentOffset;  //  total number of real pulses


      // Calculation of some output branch variables
      for (Int_t i = 0; i < nblocks; i++)
      {

        binmin = 30;
        binmax = 109;

        for (Int_t it = 0; it < nsamp; it++)
        {

          integ[i] += signal[i * ntime + it];
          integtot += signal[i * ntime + it];

          if (it > binmin && it < binmax)
          { // cosmic pulse window
            ener[i] += signal[i * ntime + it];
            enertot += signal[i * ntime + it];
          }

          if (!(it > binmin && it < binmax))
          {
            bkg[i] += signal[i * ntime + it];

          } // background window

          // if (signal[i][it] > sigmax[i])
          if (signal[i * ntime + it] > sigmax[i])
          {
            time[i] = it;
            sigmax[i] = signal[i * ntime + it];
            ampl[i] = signal[i * ntime + it];

          } // determine the pulse maximum
        }

        // doesnt seem to be used??
        ener[i] -= bkg[i] * (binmax - binmin - 1) / (nsamp - (binmax - binmin - 1)); // subtract the bkg contribution (we have to normalize since the window widths are not the same)

        bkg[i] = bkg[i] / (nsamp - (binmax - binmin - 1)); // mean value of bkg

        for (Int_t it = 0; it < nsamp; it++)
        {

          if (!(it > binmin && it < binmax))
          {
            noise[i] += (signal[i * ntime + it] - bkg[i]) * (signal[i * ntime + it] - bkg[i]) / (nsamp - (binmax - binmin - 1));
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
          // if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.5)
          if ((signal[i * ntime + it] - bkg[i]) >= ampl2[i] * 0.5)
          {
            max50 = it;
          }
          // if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.1)
          if ((signal[i * ntime + it] - bkg[i]) >= ampl2[i] * 0.1)
          {
            max90 = it;
          }
        }
        for (Int_t it = time[i]; it > -0.5; it--)
        { // aller vers la gauche du maximum
          // if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.5)
          if ((signal[i * ntime + it] - bkg[i]) >= ampl2[i] * 0.5)
          {
            min50 = it;
          }
          // if ((signal[i][it] - bkg[i]) >= ampl2[i] * 0.1)
          if ((signal[i * ntime + it] - bkg[i]) >= ampl2[i] * 0.1)
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

      if ((int)evt % 1000 == 0)
      {
        cout << " Entry = " << evt << "  cpu time=" << tlambda.RealTime() << endl;
        tlambda.Continue();

        // gObjectTable->Print();
      }

    } // end if(NSampWaveForm<=Ndata)

    ////////temp plotting on single thread
    // Only produce the PDF if this event is one you care about:
    // 1) collect exactly those blocks that had ≥1 peak:
    if (false)
    {
      std::vector<int> activeBlocks;
      for (int bn = 0; bn < nblocks; ++bn)
      {
        if (wfnpulse[bn] > 0 && finter[bn])
        {
          activeBlocks.push_back(bn);
        }
      }

      if (activeBlocks.empty())
      {
        std::cout << std::fixed << "[evt=" << evt << "] no fitted pulses found → skipping PDF.\n";
      }
      else
      {
        int Nactive = activeBlocks.size();
        int ncol_small = std::ceil(std::sqrt(double(Nactive)));
        int nrow_small = std::ceil(double(Nactive) / double(ncol_small));

        // Make a canvas sized so each pad is ~300×300 pixels
        TCanvas *c1 = new TCanvas("c1", "Fits for event", ncol_small * 300, nrow_small * 300);
        c1->Divide(ncol_small, nrow_small, 0.001, 0.001);
        gStyle->SetOptStat(0);

        // Keep pointers so we can delete them after printing:
        std::vector<TH1F *> keepHistos;
        keepHistos.reserve(Nactive);

        for (int idx = 0; idx < Nactive; ++idx)
        {
          int bn = activeBlocks[idx];
          int padIndex = idx + 1; // pads are 1–based
          c1->cd(padIndex);

          // Debug print:
          //  std::cout << "[evt=" << evt << "] drawing block " << bn
          //           << "  wfnpulse=" << wfnpulse[bn]
          //          << "  TF1_ptr=" << (void*)finter[bn].get() << "\n";

          // 2) build a histogram for the raw waveform of block `bn`
          TH1F *hraw = new TH1F(
              Form("hraw_evt%.0f_blk%03d", evt, bn),
              Form("Block %d (evt=%.0f)", bn, evt),
              ntime, 0.0, double(ntime));
          for (int it = 0; it < ntime; ++it)
          {
            hraw->SetBinContent(it + 1, signal[bn * ntime + it]);
          }
          hraw->SetLineColor(kBlack);
          hraw->Draw("hist");
          hraw->SetMinimum(hraw->GetMinimum() * 1.1 - 1.0);
          hraw->SetMaximum(hraw->GetMaximum() * 1.1 + 1.0);

          gPad->Update();
          keepHistos.push_back(hraw);

          // 3) overlay the TF1 fit, but force the visible x‐range to [0, ntime]:
          if (finter[bn])
          {
            // make sure it has at least one free parameter (otherwise it's “flat”)
            if (finter[bn]->GetNpar() > 1)
            {
              // If the seed amplitude is too small, bump it up to 2.0
              Double_t ampSeed = finter[bn]->GetParameter(3);
              // cout<<"fit params: ";
              for (int ip = 0; ip < finter[bn]->GetNpar(); ++ip)
              {
                // cout<<finter[bn]->GetParameter(ip)<<" , ";
              }
              // cout<<endl;

              finter[bn]->SetLineColor(kBlue);
              finter[bn]->SetLineWidth(2);
              // restrict the drawing range to [0..ntime]:
              finter[bn]->SetRange(0.0, double(ntime));
              finter[bn]->Draw("same");
            }
            else
            {
              std::cout << "  → TF1[" << bn << "] exists but has <=1 par ⇒ nothing to draw\n";
            }
          }
          else
          {
            std::cout << "  → finter[" << bn << "] is nullptr, skipping fit‐curve\n";
          }

          // 4) draw vertical red lines for each TSpectrum‐found peak
          for (int p = 0; p < wfnpulse[bn]; ++p)
          {
            Double_t t_ns = wftime[blockOffset[bn] + p];
            // convert that “corrected ns‐time” back to “bin index”:
            Double_t binOff = (t_ns + cortime[bn] + timerefacc * dt - corr_time_HMS) / dt + timeref[bn];

            // Debug: print out each candidate binOff
            //  std::cout << "    → peak["<<p<<"] t_ns="<<t_ns
            //            << " ⇒ binOff="<<binOff << "\n";

            // Now draw a TLine at x=binOff
            if (binOff >= 0.0 && binOff <= ntime)
            {
              Double_t ylo = hraw->GetMinimum();
              Double_t yhi = hraw->GetMaximum();
              if (yhi <= ylo)
              {
                // If flat, give yourself a small nonzero range so the vertical line is visible:
                ylo = 0;
                yhi = ylo + 1;
              }
              TLine ln(binOff, ylo, binOff, yhi);
              ln.SetLineColor(kRed);
              ln.SetLineStyle(2);
              ln.DrawClone("same");
              // TLine *ln = new TLine(binOff, ylo, binOff, yhi);

              // ln->SetLineColor(kRed);
              // ln->SetLineStyle(2);
              // ln->Draw("same");
            }
            else
            {
              std::cout << "       (binOff out of range [0,"
                        << ntime << "], skipping)\n";
            }
          } // end loop over peaks

          // 5) label the pad with block index
          TLatex tex;
          tex.SetNDC();
          tex.SetTextSize(0.04);
          tex.SetTextColor(kBlue);
          tex.DrawLatex(0.02, 0.90, Form("blk %d", bn));

          // leave hraw alive until after we print the canvas…
        } // end for each active `bn`

        // 6) finally update+print the canvas to a single‐page PDF
        c1->Update();
        TString pdfName = Form("figures/fits_run%d_evt%.0f.pdf", run, evt);
        c1->Print(pdfName);

        // 7) clean up
        for (auto h : keepHistos)
        {
          delete h;
        }
        delete c1;
      } // end if(activeBlocks non‐empty)

    } // if for plotting
    ///////// end of plotting

    tlambda.Stop();
return std::make_tuple(
  chi2, ampl, amplwf, wfnpulse,
  Sampampl, Samptime, timewf,
  enertot, integtot, pres,
  corr_time_HMS, h1time, h2time,
  ROOT::RVec<Double_t>(wfampl.begin(), wfampl.begin() + blockOffset[nblocks]),
  ROOT::RVec<Double_t>(wftime.begin(), wftime.begin() + blockOffset[nblocks])
);



  }; // End of lambda (event loop)

  /////// Write the output files //////////////////////////////////////////////////////////

  // Make output dataframe with columns as the output tuple and all the arrays within
  auto df_final = df2.Define("tuple", analyze, {"NSampWaveForm", "SampWaveForm", "evt", "NadcCounter", "adcCounter", "adcSampPulseTime", "adcSampPulseTimeRaw", "adcSampPulseAmp", "adcSampPulseInt", "adcSampPulsePed"});

  // Using auto casued compile error for lambda arguements, had to list them explicitly to work?
  df_final = df_final.Define("chi2", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
   std::vector<Double_t> output = std::get<0>(tuple);
   return output; }, {"tuple"});
  df_final = df_final.Define("ampl", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
    std::vector<Double_t> output = std::get<1>(tuple);
    return output; }, {"tuple"});
  df_final = df_final.Define("amplwf", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
    std::vector<Double_t> output = std::get<2>(tuple);
    return output; }, {"tuple"});
  df_final = df_final.Define("wfnpulse", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
      std::vector<Int_t> output = std::get<3>(tuple);
      return output; }, {"tuple"});
  df_final = df_final.Define("Sampampl", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
        std::vector<Double_t> output = std::get<4>(tuple);
        return output; }, {"tuple"});
  df_final = df_final.Define("Samptime", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
          std::vector<Double_t> output = std::get<5>(tuple);
          return output; }, {"tuple"});
  df_final = df_final.Define("timewf", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
            std::vector<Double_t> output = std::get<6>(tuple);
            return output; }, {"tuple"});
  df_final = df_final.Define("enertot", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
              Double_t output = std::get<7>(tuple);
              return output; }, {"tuple"});
  df_final = df_final.Define("integtot", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
                Double_t output = std::get<8>(tuple);
                return output; }, {"tuple"});
  df_final = df_final.Define("pres", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
                  std::vector<Int_t> output = std::get<9>(tuple);
                  return output; }, {"tuple"});
  df_final = df_final.Define("corr_time_HMS", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
                    Double_t output = std::get<10>(tuple);
                    return output; }, {"tuple"});
  df_final = df_final.Define("h1time", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
                      std::vector<Double_t> output = std::get<11>(tuple);
                      return output; }, {"tuple"});
  df_final = df_final.Define("h2time", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
                        std::vector<Double_t> output = std::get<12>(tuple);
                        return output; }, {"tuple"});
  df_final = df_final.Define("wfampl", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
                          ROOT::RVec<Double_t> output = std::get<13>(tuple);
                          return output; }, {"tuple"});
  df_final = df_final.Define("wftime", [](const std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Int_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, Double_t, Double_t, std::vector<Int_t>, Double_t, std::vector<Double_t>, std::vector<Double_t>, ROOT::RVec<Double_t>, ROOT::RVec<Double_t>> &tuple)
                             { 
                            ROOT::RVec<Double_t> output = std::get<14>(tuple);
                            return output; }, {"tuple"});

  auto h_h1time = df_final.Histo1D(m_h1time, "h1time");
  auto h_h2time = df_final.Histo1D(m_h2time, "h2time");

  // Save the dataframe to the output ROOT file
  // TString rootfilePath = Form("/volatile/hallc/nps/kerver/ROOTfiles/WF/nps_production_%d_%d_%d_interactive_allwf_0.05fix_2mvthresh.root", run, seg,nthreads);
  // TString rootfilePath = Form("../nps_production_wf_%d_%d_%d.root", run, seg,nthreads);
  auto columnNames = df_final.GetColumnNames();
  for (const auto &col : columnNames)
  {
    // std::cout << col << std::endl;
  }

  ///////////
  cout << "About to Snapshot" << endl;
  TString tempout = Form("../nps_production_%d_%d_%d_temp.root", run, seg, nthreads);
  ROOT::RDF::RSnapshotOptions opts;
  opts.fMode = "RECREATE";
  // opts.fBasketSize = 512 * 1024;        // 512 KB per basket
  df_final.Snapshot("WF", tempout, {"chi2", "ampl", "amplwf", "wfnpulse", "Sampampl", "Samptime", "timewf", "enertot", "integtot", "pres", "corr_time_HMS", "h1time", "h2time", "evt", "wfampl", "wftime"}, opts);
  t.Stop();
  std::cout
      << "=== Snapshot done. Elapsed real time: " << t.RealTime()
      << " s, CPU time: " << t.CpuTime()
      << " s ===\n";
  t.Continue();

  // 1) Open and tune baskets
  TFile *fin = TFile::Open(tempout, "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "ERROR: could not open " << tempout << " for basket-tuning\n";
    return;
  }
  // TTree *tin = (TTree*)fin->Get("WF");
  TTree *tin = static_cast<TTree *>(fin->Get("WF"));
  if (!tin)
  {
    std::cerr << "ERROR: WF not found in " << tempout << "\n";
    fin->Close();
    return;
  }
  tin->BuildIndex("evt");

  // 3) Clone into your main file (atomic overwrite)
  TFile *fout = TFile::Open(outFile, "UPDATE");
  if (!fout || fout->IsZombie())
  {
    std::cerr << "Cannot open output file " << outFile << "\n";
    fin->Close();
    return;
  }
  fout->cd();
  auto tout = tin->CloneTree(-1);
  tout->BuildIndex("evt");

  std::cout
      << "=== Sorting done. Elapsed real time: " << t.RealTime()
      << " s, CPU time: " << t.CpuTime()
      << " s ===\n";
  t.Continue();

  tout->Write("WF", TObject::kOverwrite);
  fout->Close();
  fin->Close();

  cout << "fin de l'analyse" << endl;
  // cout << "Failed fits: "<<failcount<<endl;
  std::cout << "Total failed fits: " << nFitFailures.load() << " total fits succeed:" << nFitSucceeds.load() << "\n";
  t.Stop();
  t.Print();
}
