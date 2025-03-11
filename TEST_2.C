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

Double_t timeref[nblocks];
Double_t timerefacc = -4.5; // car les bons pulses elastiques poussent vers 45.5 (4*ns) dans la wf alors que ceux des production runs sont vers 35.5 (4ns)
Double_t timerefacc2 = 0;   // en (-4ns) car le pic de timewf pousse vers -22ns

ROOT::Math::Interpolator *interpolation[nblocks]; // function used for the interpolation
// TH1F *hsig_i[nblocks];                            // histogram showing the waveform for block i and a given event
// TF1 *finter[nblocks];

// FIT FUNCTION - redfined in main function as a lmbda function for MT
/*
Double_t func(Double_t *x, Double_t *par)
{
  Int_t j = (int)(par[0]);
  Double_t val = 0;
  for (Int_t p = 0; p < maxwfpulses; p++)
  {
    if (x[0] - par[2 + 2 * p] > 1 && x[0] - par[2 + 2 * p] < 109)
      val += par[3 + 2 * p] * interpolation[j]->Eval(x[0] - par[2 + 2 * p]);
  }
  return val + par[1];
}
  */

///////////////////////////////////////////////////////////// FIT FUNCTION USED AND SCHEMATICS /////////////////////////////////////////////////////////////

auto Fitwf = [&](Int_t bn)
{
  // Detect the pulses

  wfnpulse[bn] = 0;
  Int_t good = 0;

  for (Int_t p = 0; p < maxwfpulses; p++)
  {

    // cout << "Setting parameters for pulse " << p << " with index " << 2+2*p << " and " << 3+2*p << endl;

    wfampl[bn][p] = -1;
  }

  for (Int_t it = 0; it < ntime - 6; it++) // loop over number of samples(110) in an event for a given block
  {

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

    // Condition over the number of samples in the pulse finding scheme

    if (hsig_i[bn]->GetBinContent(it + 1) < hsig_i[bn]->GetBinContent(it + 2) && hsig_i[bn]->GetBinContent(it + 2) < hsig_i[bn]->GetBinContent(it + 3) && hsig_i[bn]->GetBinContent(it + 3) <= hsig_i[bn]->GetBinContent(it + 4) && hsig_i[bn]->GetBinContent(it + 4) >= hsig_i[bn]->GetBinContent(it + 5))
    {

      if (hsig_i[bn]->GetBinContent(it + 4) > 0)
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

        if (TMath::Abs(wftime[bn][wfnpulse[bn]] - timeref[bn] - timerefacc) < 4.1)
        {
          good = 1;
        }

        wfnpulse[bn]++;

        // to prevent overflow
        if (wfnpulse[bn] == maxwfpulses)
        {
          wfnpulse[bn] = maxwfpulses - 1;
          it = ntime;
        }

        it += 4;

      } // end of the condition over (hsig_i[bn]->GetBinContent(it+4)>0){
    } // end of the condition over samples
  } // end of loop over it

  if (wfnpulse[bn] > maxwfpulses - 2)
  {
    cout << "Warning : excessively high number of pulses in the block wf " << bn << endl;
  }

  // Adjust the parameters of the fit function

  // File<<"c_0"<<endl;

  for (Int_t p = 0; p < maxwfpulses; p++)
  {
    // File<<"c1"<<endl;
    finter[bn]->FixParameter(2 + 2 * p, 0.);
    // File<<"c2"<<endl;
    finter[bn]->FixParameter(3 + 2 * p, 0.);
  }

  if (wfnpulse[bn] > 0 && good == 1)
  {

    for (Int_t p = 0; p < TMath::Min(maxwfpulses, wfnpulse[bn]); p++)
    {

      finter[bn]->ReleaseParameter(2 + 2 * p);
      finter[bn]->ReleaseParameter(3 + 2 * p);

      finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn]);
      finter[bn]->SetParameter(3 + 2 * p, wfampl[bn][p]);

      finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] - 3, wftime[bn][p] - timeref[bn] + 3);
      finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn][p] * 0.2, wfampl[bn][p] * 3);
    }
  }

  // File<<"c3"<<endl;

  if (wfnpulse[bn] > 0 && good == 0)
  {

    for (Int_t p = 0; p < TMath::Min(maxwfpulses, wfnpulse[bn]); p++)
    {

      finter[bn]->ReleaseParameter(2 + 2 * p);
      finter[bn]->ReleaseParameter(3 + 2 * p);

      finter[bn]->SetParameter(2 + 2 * p, wftime[bn][p] - timeref[bn]);

      // Check to debug
      if (wfampl[bn][p] > 0)
      { // Ensure it's positive before multiplying
        finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn][p] * 0.2, wfampl[bn][p] * 3);
      }
      else
      {
        cout << "Warning: wfampl[" << bn << "][" << p << "] is non-positive!" << endl;
        // finter[bn]->SetParLimits(3 + 2 * p, 0.05, 10); // Set safe fallback limits
      }

      finter[bn]->SetParameter(3 + 2 * p, wfampl[bn][p]);

      finter[bn]->SetParLimits(2 + 2 * p, wftime[bn][p] - timeref[bn] - 3, wftime[bn][p] - timeref[bn] + 3);
      finter[bn]->SetParLimits(3 + 2 * p, wfampl[bn][p] * 0.2, wfampl[bn][p] * 3);
    }

    // File<<"c4"<<endl;

    // On recherche quand meme un eventuel pulse en temps

    finter[bn]->ReleaseParameter(2 + 2 * wfnpulse[bn]);
    finter[bn]->ReleaseParameter(3 + 2 * wfnpulse[bn]);

    finter[bn]->SetParameter(2 + 2 * wfnpulse[bn], timerefacc);
    finter[bn]->SetParameter(3 + 2 * wfnpulse[bn], 2);

    finter[bn]->SetParLimits(2 + 2 * wfnpulse[bn], timerefacc - 4, timerefacc + 4);
    finter[bn]->SetParLimits(3 + 2 * wfnpulse[bn], 0.05, 10);

    wfnpulse[bn]++;
  }

  // File<<"c5"<<endl;

  if (wfnpulse[bn] == 0)
  {

    for (Int_t p = 0; p < 1; p++)
    {

      finter[bn]->ReleaseParameter(2 + 2 * p);
      finter[bn]->ReleaseParameter(3 + 2 * p);

      finter[bn]->SetParameter(2 + 2 * p, timerefacc);
      finter[bn]->SetParameter(3 + 2 * p, 2);

      finter[bn]->SetParLimits(2 + 2 * p, timerefacc - 4, timerefacc + 4);
      finter[bn]->SetParLimits(3 + 2 * p, 0.05, 10);
    }

    wfnpulse[bn]++;
  }

  // File<<"c6"<<endl;

  finter[bn]->SetParameter(1, 0.);
  finter[bn]->SetParLimits(1, -100, 100.);

  // File<<"c7"<<endl;
  // hsig_i[bn]->Fit(Form("finter_%d",bn),"Q","",TMath::Max(0.,wftime[bn][0]-20),ntime);
  // File << "Checkpoint: Before fitting function for bn = " << bn << endl;

  hsig_i[bn]->Fit(Form("finter_%d", bn), "Q", "", 0., ntime);

  // File<<"c8"<<endl;
}

/////////////////////////////////////////////////////////////MAIN FUNCTION ///////////////////////////////////////////////////////////////////////

void
TEST_2(int run, int seg)
{

  // FIT FUNCTION redefined as a lambda function for MT use
  auto func = [=](Double_t *x, Double_t *par)
  {
    Int_t j = (int)(par[0]);
    Double_t val = 0;
    for (Int_t p = 0; p < maxwfpulses; p++)
    {
      if (x[0] - par[2 + 2 * p] > 1 && x[0] - par[2 + 2 * p] < 109)
        val += par[3 + 2 * p] * interpolation[j]->Eval(x[0] - par[2 + 2 * p]);
    }
    return val + par[1];
  };

  TStopwatch t;
  t.Start();

  // ENABLE MT
  const int nthreads = 8; // or any number
  ROOT::EnableImplicitMT(nthreads);
  std::cout << "Implicit MT enabled: " << ROOT::IsImplicitMTEnabled() << "\n";
  std::cout << "Number of threads: " << ROOT::GetThreadPoolSize() << "\n";

  // BUILD TCHAIN
  TChain chain("T");
  TString fname = Form("/cache/hallc/c-nps/analysis/online/replays/nps_hms_coin_%d_%d_1_-1.root", run, seg);
  TFile *testOpen = TFile::Open(filename);
  if (!testOpen || testOpen->IsZombie())
  {
    std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
    return;
  }
  testOpen->Close();
  chain.Add(fname);
  auto nEntries = chain.GetEntries(); // replaced wassims nentries
  std::cout << "Chain has " << nEntries << " entries.\n";

  // Before processing events, check the list of objects in memory
  std::cout << "Objects in memory before event processing:" << std::endl;
  gObjectTable->Print();

  // Create an RDataFrame from the chain
  ROOT::RDataFrame df(chain);

  // Flatten the waveforms to handle multiple per entry
  // auto dfFlattened = df.Define("single_waveform", "Flatten(waveform)");

  if (run == 2016)
  {
    timerefacc = -10;
  } // keep this condition and might add more runs if that happens again

  // timerefacc = timerefacc + 10;
  // timerefacc2 = timerefacc2 + 10 * 4;

  // int end_event, start_event;
  // const int evts_per_job = 250000;
  Double_t dt = 4.;                                    // time bin (sample) width (4 ns), the total time window is then ntime*dt
  const int nslots = 1104;                             // nb maximal de slots dans tous les fADC
  const int Ndata = nslots * (ntime + 1 + 1);          //(1104 slots fADC au total mais pas tous utilises y compris 2 PM scintillateurs?) should correspond to Ndata.NPS.cal.fly.adcSampWaveform variable in the root file
  Double_t ADCtomV = 1000. / 4096;                     // 1000 mV correspond a 4096 canaux
  Double_t integtopC = (1. / 4096) * (dt * 1.e3 / 50); // (1 V / 4096 adc channels) * (4000 ps time sample / 50 ohms input resistance) = 0.020 pc/channel

  // where the cosmic pulse is expected to be
  Int_t binmin;
  Int_t binmax;

  // Read reference waveforms
  const int nbparameters = 2 + 2 * maxwfpulses;
  cout << "nbparameters = " << nbparameters << endl;
  Double_t x[ntime], y[ntime];
  ifstream filewf;
  Double_t dum1, ymax; // dont seem to be used anywhere?

  for (Int_t i = 0; i < nblocks; i++)
  {

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
    interpolation[i] = new ROOT::Math::Interpolator(ntime, ROOT::Math::Interpolation::kCSPLINE);
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
    }
    filewf.close();
    // Moved inside lamnda function. need an interpolation function array for each thread
    /*
    finter[i] = new TF1(Form("finter_%d", i), func, -200, ntime + 200, nbparameters);
    finter[i]->FixParameter(0, i); // numero de bloc
    finter[i]->SetLineColor(4);
    finter[i]->SetNpx(1100);
  */
    // cout<<i<<" "<<timeref[i]<<endl;
  }

  // Read tdc_offset_param (needed to determine HMS corrections to the timing)//For now, it is just one file
  Float_t tdcoffset[nblocks];
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6171-6183/fit_e_runs/fit_elastic_runs/results_elastics/refwf/tdc_offset_param.txt");//Done
  ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/tdc_offset_param.txt"); // Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/tdc_offset_param.txt");//Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5183-5208/fit_e_runs/RWF/tdc_offset_param.txt");//Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/RWF/tdc_offset_param.txt");//Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2900-2920/RWF/tdc_offset_param.txt");//Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2875-2885/RWF/tdc_offset_param.txt");//Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2855-2871/RWF/tdc_offset_param.txt");//Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/tdc_offset_param.txt");//Done
  // ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1423-1511/RWF/tdc_offset_param.txt");//Done

  for (Int_t i = 0; i < nblocks; i++)
  {
    filetdc >> tdcoffset[i];
  }

  // Load only required branches into dataframe
  auto df = df.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")
                .Define("SampWaveForm", "NPS.cal.fly.adcSampWaveform")
                .Define("NadcCounter", "Ndata.NPS.cal.fly.adcCounter")
                .Define("adcCounter", "NPS.cal.fly.adcCounter")
                .Define("NadcSampPulseAmp", "Ndata.NPS.cal.fly.adcSampPulseAmp")
                .Define("adcSampPulseAmp", "NPS.cal.fly.adcSampPulseAmp")
                .Define("NadcSampPulseInt", "Ndata.NPS.cal.fly.adcSampPulseInt")
                .Define("adcSampPulseInt", "NPS.cal.fly.adcSampPulseInt")
                .Define("NadcSampPulsePed", "Ndata.NPS.cal.fly.adcSampPulsePed")
                .Define("adcSampPulsePed", "NPS.cal.fly.adcSampPulsePed")
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
                .Define("beta", "beta")
                .Define("cernpeSum", "cernpeSum")
                .Define("caletracknorm", "caletracknorm")
                .Define("caletottracknorm", "caletottracknorm")
                .Define("caletotnorm", "caletotnorm")
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
                .Define("caltracknorm", "H.cal.etottracknorm");
  .Define("evt", "rdfentry_"); // define event number from row of rdataframe for troubleshoooting (threadsafe)

  // Output rootfile

  Int_t tracage = 0; // 1 to see the waveforms event by event

  TFile *fout;

  if (tracage == 0)
  {

    TString rootfilePath = Form("nps_production_%d_%d.root", run, seg);

    cout << "Rootfile location: " << rootfilePath << endl;
    fout = new TFile(rootfilePath, "recreate");
  }

  TTree *treeout = new TTree("T", "Tree organized");
  treeout->SetAutoFlush(-300000000);

  /////TODO: convert all needed treeout branches to df.Histo1D() then ->write() at end of code
  
  // Output variables


  // Double_t signal[nblocks][ntime]; // treeout->Branch("signal",&signal,Form("signal[%d][%d]/D",nblocks,ntime));//if we want to save the raw waveforms in the output rootfile
  // no global arays. needs to be created within lamda function for it to be thread safe.

  // Double_t ampl2[nblocks]; // treeout->Branch("ampl2",&ampl2,Form("ampl2[%d]/D",nblocks));//amplitude of the pulse relatively to background
  // Double_t ampl[nblocks];
  treeout->Branch("ampl", &ampl, Form("ampl[%d]/D", nblocks)); // amplitude of the pulse (a comparer avec Sampampl)
  // Double_t amplwf[nblocks];
  treeout->Branch("amplwf", &amplwf, Form("amplwf[%d]/D", nblocks)); // amplwfitude of the pulse (a comparer avec Sampamplwf)
  // treeout->Branch("wfampl",&wfampl,Form("wfampl[%d][%d]/D",nblocks,maxwfpulses));
  // treeout->Branch("wftime",&wftime,Form("wftime[%d][%d]/D",nblocks,maxwfpulses));
  treeout->Branch("wfnpulse", &wfnpulse, Form("wfnpulse[%d]/I", nblocks)); // number of pulses in one block
  treeout->Branch("Sampampl", &Sampampl, Form("Sampampl[%d]/D", nblocks)); // amplitude of the 1st pulse
  treeout->Branch("Samptime", &Samptime, Form("Samptime[%d]/D", nblocks)); // time of the 1st pulse 'ns)
  // Double_t ener[nblocks];
  //  Double_t time[nblocks];
  // Double_t timewf[nblocks];
  treeout->Branch("timewf", &timewf, Form("timewf[%d]/D", nblocks)); // timewf position of the pulse maximum
  treeout->Branch("chi2", &chi2, Form("chi2[%d]/D", nblocks));       // time of the 1st pulse 'ns)
  treeout->Branch("enertot", &enertot, "enertot/D");                 // sum of all ener[i] : total energy deposited in the calo (!! energy calibration is not done yet)
  treeout->Branch("integtot", &integtot, "integtot/D");              // sum of all integ[i]
  // Double_t larg50[nblocks];                                          // treeout->Branch("larg50",&larg50,Form("larg50[%d]/D",nblocks));//RMS of waveforms relatively to zero
  // Double_t larg90[nblocks];                                          // treeout->Branch("larg90",&larg90,Form("larg90[%d]/D",nblocks));//RMS of waveforms relatively to zero
  // Int_t pres[nblocks]; //no glbal arays. needs to be created within lamda function for it to be thread safe.
  // treeout->Branch("pres", &pres, Form("pres[%d]/I", nblocks)); // RMS of waveforms relatively to zero
  treeout->Branch("seg", &seg, "seg/I");
  // treeout->Branch("job_number",&job_number,"job_number/I");
  // Int_t event;
  treeout->Branch("event", &event, "event/I"); // event number in input rootfiles
  treeout->Branch("beta", &beta, "beta/D");
  treeout->Branch("cernpeSum", &cernpeSum, "cernpeSum/D");
  treeout->Branch("caletracknorm", &caletracknorm, "caletracknorm/D");
  treeout->Branch("caletottracknorm", &caletottracknorm, "caletottracknorm/D");
  treeout->Branch("caletotnorm", &caletotnorm, "caletotnorm/D");
  treeout->Branch("corr_time_HMS", &corr_time_HMS, "corr_time_HMS/D");
  treeout->Branch("dp", &dp, "dp/D");
  treeout->Branch("th", &th, "th/D");
  treeout->Branch("ph", &ph, "ph/D");
  treeout->Branch("px", &px, "px/D");
  treeout->Branch("py", &py, "py/D");
  treeout->Branch("pz", &pz, "pz/D");
  treeout->Branch("vx", &vx, "vx/D");
  treeout->Branch("vy", &vy, "vy/D");
  treeout->Branch("vz", &vz, "vz/D");

  // Other variables
  Int_t ilin, icol, ilinc, icolc, inp;
  Int_t nsampwf = 0;
  Int_t ndataprob = 0;
  Int_t nb = 0;

  // Histograms
  /* Moved inside threadsafe lambda function for now. not thew way to do it. should define df.histo1D 
  for (Int_t i = 0; i < nblocks; i++)
  {
    hsig_i[i] = new TH1F(Form("hsig_i%d", i), Form("hsig_i%d", i), ntime, 0, ntime);
    hsig_i[i]->SetLineColor(1);
    hsig_i[i]->SetLineWidth(2);
    hsig_i[i]->GetXaxis()->SetTitle("Time (4 ns)");
    hsig_i[i]->GetYaxis()->SetTitle("(mV)");
    hsig_i[i]->GetYaxis()->SetLabelSize(0.05);
  }
    */

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.015);
  TH1::AddDirectory(kFALSE);

  // Read the mean time positions of cosmic pulses (this is determined by the macro analyse_wassim.C)
  Float_t timemean[nblocks], timemean2[nblocks];
  for (Int_t ii = 0; ii < nblocks; ii++)
  {
    timemean[ii] = 50;
    timemean2[ii] = 150;
  } // 1st pass analysis, when the macro analyse_wassim.C is not executed yet (30 is the default value)
  // 2nd pass analysis, when the macro analyse_wassim.C is already executed. If not, comment the following lines
  ifstream filetime("filetime.txt");
  Float_t dum;

  TH1F *h1time = new TH1F("h1time", "pulse (>20mV) shift (4*ns units) relatively to elastic refwf (all found pulses included)", 200, -50, 50);
  TH1F *h2time = new TH1F("h2time", "pulse (>20mV) time (ns) (all found pulses included)", 200, -100, 100);

  /////// ANALYZE //////////////////////////////////////////////////////////

  // this is where sequential event loop was

  // HMS CUTS
  auto d2 = df.Filter(TMath::Abs(th) < 0.08)
                .Filter(Math::Abs(ph) < 0.04)
                .Filter(TMath::Abs(dp) < 10);

  ////////Lambda funtion for the per-event wf analysis/////////
  auto analyze = [=](Int_t NSampWaveForm, const std::vector<Int_t> &SampWaveForm, Int_t evt, Int_t NadcCounter, const std::vector<Int_t> &adcCounter, const std::vector<Int_t> adcSampPulseTime, const std::vector<Int_t> adcSampPulseTimeRaw)
  {
    TF1 *finter[nblocks];
    TH1F *hsig_i[nblocks];
    Int_t pres[nblocks];
    Double_t signal[nblocks][ntime];
    Int_t bloc, nsamp;
    if (NSampWaveForm > Ndata)
    {
      nsampwf++;
      cout << "!!!! NSampWaveForm problem  " << evt << "  " << NSampWaveForm << " " << Ndata << endl;
    }

    if (NSampWaveForm <= Ndata)
    { // NSampWaveForm must be <= Ndata (otherwise correct Ndata value)

      for (Int_t i = 0; i < nblocks; i++)
      {

        hsig_i[i] = new TH1F(Form("hsig_i%d", i), Form("hsig_i%d", i), ntime, 0, ntime);
        hsig_i[i]->SetLineColor(1);
        hsig_i[i]->SetLineWidth(2);
        hsig_i[i]->GetXaxis()->SetTitle("Time (4 ns)");
        hsig_i[i]->GetYaxis()->SetTitle("(mV)");
        hsig_i[i]->GetYaxis()->SetLabelSize(0.05);

        pres[i] = 0;

        for (Int_t j = 0; j < ntime; j++)
        {
          signal[i][j] = 0;
        }

      } // liste de presence des blocs! pres=0 si bloc absent, pres=1 s'il est present

      // Extract the data from the complex variable SampWaveForm[] (NPS.cal.fly.adcSampWaveform)

      int ns = 0; // ns represent for a given event the element number of the NPS.cal.fly.adcSampWaveform variable
      //    nb=0;
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

          if (nsamp == ntime)
          {
            for (Int_t it = 0; it < nsamp; it++)
            {

              if (bloc > -0.5 && bloc < nblocks)
              {
                signal[bloc][it] = SampWaveForm[ns];
                hsig_i[bloc]->SetBinContent(it + 1, signal[bloc][it]);
              }

              ns++;
            }
          }

        } // fin if(bloc number is good)

      } // fin while()

      // Calculate the output variables of the rootfile

      if (pres[454] == 0)
        nb++;

      // initialisation

      Double_t enertot = 0.;
      Double_t integtot = 0.;
      Double_t corr_time_HMS = 0.;

      Double_t timewf[nblocks];
      Double_t amplwf[nblocks];
      Double_t chi2[nblocks];
      Double_t ener[nblocks];
      Double_t integ[nblocks];
      Double_t noise[nblocks];
      Double_t bkg[nblocks];
      Double_t sigmax[nblocks];
      Double_t Sampampl[nblocks];
      Double_t Sampped[nblocks];
      Double_t Samptime[nblocks];
      Double_t Sampener[nblocks];
      Double_t Npulse[nblocks];

      Int_t wfnpulse[nblocks];
      Double_t wfampl[nblocks][maxwfpulses];
      Double_t wftime[nblocks][maxwfpulses];
      // not used but here they are anyway
      Double_t ampl2[nblocks];
      Double_t ampl[nblocks];
      Double_t time[nblocks];
      Double_t larg50[nblocks];
      Double_t larg90[nblocks];
      Double_t max50, max90, min50, min90;

      for (Int_t i = 0; i < nblocks; i++)
      {
        timewf[i] = -100;
        amplwf[i] = -100;
        chi2[i] = -1;
        ener[i] = 0;
        integ[i] = 0;
        noise[i] = 0;
        bkg[i] = 0;
        sigmax[i] = -100;
        Sampampl[i] = -100;
        Sampped[i] = -100;
        Samptime[i] = -100;
        Sampener[i] = -100;
        Npulse[i] = 0;

        finter[i] = new TF1(Form("finter_%d", i), func, -200, ntime + 200, nbparameters);
        finter[i]->FixParameter(0, i); // numero de bloc
        finter[i]->SetLineColor(4);
        finter[i]->SetNpx(1100);
      }

      for (Int_t i = 0; i < nblocks; i++)
      {
        wfnpulse[i] = -100;

        for (Int_t p = 0; p < maxwfpulses; p++)
        {
          wftime[i][p] = -100;
          wfampl[i][p] = -100;
        }
      }
      // End per-event intializations

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
            cout << "problem HMS time correction event " << evt << endl;
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

        // icol = i % ncol; //never used anywhere?
        // ilin = (int)(i / ncol);

        if (pres[i] == 1)
        { // For this setting?

          for (Int_t it = 0; it < nsamp; it++)
          {

            hsig_i[i]->SetBinError(it + 1, TMath::Sqrt(TMath::Abs(hsig_i[i]->GetBinContent(it + 1) * 4.096)) / 4.096);

            if (hsig_i[i]->GetBinContent(it + 1) < 1.)
            {

              hsig_i[i]->SetBinError(it + 1, TMath::Sqrt(TMath::Abs(1. * 4.096)) / 4.096);
            }

          } // end loop over it

          if (nsamp == ntime)
          {
            Fitwf(i); // We call the fit here
          }

          // **Add this check to skip the block if finter[i] is invalid**
          if (!finter[i])
          {
            cout << "Skipping block " << i << " due to missing finter[" << i << "] at event " << evt << endl;
            continue;
          }

          chi2[i] = finter[i]->GetChisquare() / finter[i]->GetNDF();

          for (Int_t p = 0; p < TMath::Max(wfnpulse[i], 1); p++)
          {

            wftime[i][p] = finter[i]->GetParameter(2 + 2 * p) * dt + corr_time_HMS; // temps du pulse en ns

            wfampl[i][p] = finter[i]->GetParameter(3 + 2 * p); // amplitude du pulse en ns

            if (wfampl[i][p] > 20)
            {

              /// fix later //// diagnostic histos
              /*
                h2time->Fill(wftime[i][p], 1.);

                h1time->Fill(finter[i]->GetParameter(2 + 2 * p), 1.);
              */
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

              if (TMath::Abs(wftime[i][p] - timerefacc2 * dt) < TMath::Abs(timewf[i] - timerefacc2 * dt))
              {
                timewf[i] = wftime[i][p];
                amplwf[i] = wfampl[i][p];
              }
            } // pour prendre le pulse dont le temps est le plus proche de timerefacc

          } // end of loop over maxpulses
        } // end of condition over col,lin
      } // end of loop over les blocs

      // can all these differenet loops over block# be combined?????

      // Calculation

      for (Int_t i = 0; i < nblocks; i++)
      {

        binmin = (int)(timemean[i] - 25);
        binmax = (int)(timemean[i] + 59);

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
      }

      // Calculation of signal widths
      for (Int_t i = 0; i < nblocks; i++)
      {

        ampl2[i] = ampl[i] - bkg[i]; // amplitude of the pulse relatively to bkg

        max50 = 0;
        max90 = 50;
        min50 = 100;
        min90 = 100;

        for (Int_t it = time[i]; it < nsamp; it++)
        {
          // aller vers la droite du maximum
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

      treeout->Fill();

      // std::cout << "we just filled the tree" << std::endl;
      // gObjectTable->Print();

    } // end if(NSampWaveForm<=Ndata)
  }

  // end of the cut over HMS

  // After processing events, check the list of objects in memory

  if (evt % 1000 == 0)
  {
    std::cout << "Objects in memory After event processing:" << std::endl;
    gObjectTable->Print();
  }

  // end for(evt)

  /////////////////////////////////////////////////////////////////////////////////////////

  // Removed the chunk responsable for the event per event check for now

  /////// Write the output files //////////////////////////////////////////////////////////

  for (int i = 0; i < nblocks; i++)
  {
    delete hsig_i[i];        // Delete histogram objects
    delete finter[i];        // Delete fit function objects
    delete interpolation[i]; // Delete interpolator objects
  }

  // Then safely write and delete at the end
  fout->cd();
  h1time->Write("h1time");
  h2time->Write("h2time");
  delete h1time;
  delete h2time;
  fout->Write();
  fout->Close();

  cout << "fin de l'analyse" << endl;
  cout << nsampwf << " " << ndataprob << " " << nb << endl;

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
