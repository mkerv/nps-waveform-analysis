
# NPS Waveform Fitting

## Overview
`npsWF.C` is a multi-threaded implementation of waveform fitting code orginally prepared by Wassim Hamdi(wassim@jlab.org) and Malek Mazouz (mazouz@jlab.org) to fit data from the NPS rg1a experiment. The code is written in ROOT/C++ and relies heavily on ROOT's RDataframes for the multi-threading. The code does the following:
1. Reads in raw ADC waveforms and per-block reference waveforms  
2. Applies a matched filter to the raw waveforms  
3. Finds candidate peaks via **TSpectrum**  
4. Fits each pulse using **Minuit2** (Fitter API)  
5. Writes per-pulse time/amplitude into a new `WF` TTree in your output root file  

## Prerequisites
- **ROOT ≥ 6.30.04** 
- load halla module **nps_replay/5.28.24**
  - pythia 6.4.28, pythia 8.311, evio 5.3, analyzer 1.7.12, hcana 1.1, NPSlib fddb8e9_1.1.0

- https://github.com/mkerv/nps-waveform-analysis


## Usage

```bash
# With run, segment, threads, diagnostics
root -l -b -q '.X npsWF.C+(run, seg, nThreads, makeDiagnostics)'
```

## Command-Line Arguments

| Position | Name            | Type | Default | Description                                      |
| -------: | --------------- | ---- | ------- | ------------------------------------------------ |
|        1 | run             | int  | 0       | Run number                                       |
|        2 | seg             | int  | 0       | Segment ID                                       |
|        3 | nThreads        | int  | 4       | Number of threads               |
|        4 | makeDiagnostics | bool | false   | If true, per-event diagnostic PDFs are generated |

## Workflow Chronology
1. **Copy original ROOT file**: `inputFile` → `outputFile` via `CloneTree("fast")` preserves all TTrees & objects from the original input Replay file **except** the `NPS.cal.fly.adcSampWaveform` branch, which contained the raw waveform data that is not longer needed. 
2. **Enable multithreading**: `ROOT::EnableImplicitMT(nThreads)`  
3. **Populate RDataFrame**: Read in input file as a TChain, and create dataframe from it ` ROOT::RDataFrame df(chain);` 
    * Define columns in the dataframe for **only** the branches needed for waveform fitting via `df.Define("NSampWaveForm", "Ndata.NPS.cal.fly.adcSampWaveform")`
4. **Read in timing offsets**: A per-block timing offset is read in and saved from file :`/work/hallc/nps/wassim/ANALYSIS/TEST/filetime_step_i.txt`
5. **Read in reference waveforms**: Per-block reference waveform data is read in from text files. The text files selected are run number dependent. The units of the data are amplitude in mV and time in 4ns bins.
6. **Call Analyze**: `analyze` is a lambda function containing the rest of the logic of the code. When the dataframe calls `Snapshot` all code in this lambda function is executed for one dataframe "row" (an event) on one thread. All varaibles and methods in this function will be local to one thread.
7. **Get raw ADC data**: `vector<Double_t> signal(nblocks * 110)` is filled with the raw ADC values for all time steps(0 to 110) for all blocks.
8. **Compute HMS time correction**:  
9. **Loop over blocks**: the rest of the following code logic is executed in a loop over all 1080 blocks
10. **Find pulses**: `FindPulsesMF` is called to find all pulses in a given block.
    1.  Build a `Matched Filter` from the block's raw adc data and it's reference wavefrom data. See Section Matched Filter for details.
    2.  Call `TSpectrum.Search()` to find all peaks in the matched filter above a 1.5mV threshold. The max number of peaks found is limited to 12. See https://root.cern.ch/doc/master/classTSpectrum.html#a5ba181a45b117934b48c4ef5f78d0b2b for information on the search alogorithm
    3.  `FindPulsesMF` returns `vector<Double_t> wfampl` and `vector<Double_t> wftime` which contain the amplitude and time bin(in 4ns units) for each pulse found.
11. **Cluster threshold check**: `PassClusterThreshold` is called, it does the following:
    1.  Sums the ADC values for the active block and it's surrounding 8 blocks at every time bin.
    2.  Finds the maximum ADC sum value within a `coinc_width = 20` bin window around the reference waveform time. If that ADC sum is < `trig_thres = 10` then the block skipped and is not fitted.

12. **Waveform fit**: Call `Fitwf` 
    1.  Create a `ROOT::Interpolator` with cubic splines over `ntime=110` points.filled with the reference waveform data points. This turns the descrete reference waveform data into a continuous differentiable function.
    2. Define a model function:  
    $f(x;\,\mathbf{p}) \;=\; p_0 \;+\;\sum_{n=0}^{N-1} A_{n}\;\bigl[\mathrm{ref}(x - t_{n})\bigr]$

       * Where $p_{0}$ is the constant pedestal. For each pulse $n$ :  $t_{n}$ is the time offset (in ADC-bin units) from the block’s reference pulse. $A_{n}$ is the amplitude of that pulse, and ref(t) returns the reference amplitude scale at time t normalized to 1.
    4. Wrap the model function in a `ROOT::TF1` and set the initial parameters of the model function as the `wfampl` and `wftime - timeref` for all the pulses found in `FindPulsesMF`. Note that the time is set as a distance away from the reference pulse time
    5. Creat a `ROOT::Fit::Fitter` and configure it with your desired minimizer. By default we use Minuit2's Migrad. See https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.html
    6. Provide the fitter with the raw binned data, the interpolated reference waveform, and the model. Call `fitter.LeastSquareFit(data)` to perform the fit. 
         $\chi^2 \;=\;\sum_i \frac{\bigl(y_i - f(x_i;p)\bigr)^2}{\sigma_i^2}$        $
    7. If the minimizer fails to converge, the fit will fail and reattempt with more iterations and a higher minimizer "stategy". If it fails once again the output results will sipmly be the time/amp of the raw data as found by `TSpectrum`
    8. The fit results are stored in `RVec<std::vector<Double_t>> wfampl` , `RVec<std::vector<Double_t>> wftime` , and `vector<std::vector<Double_t>> chi2`
        - Note: Before saving the fit result, the time is converted from 4ns bins to a time in ns with a refernce time offset and per-block time correction. See Section: Time Correction for more details
  
13. **Exit the Block Loop**: after the fit some diagnostic variables are computed and a debug option can be toggled to print PDFs of fitted waveforms. This will create plots for **all** fitted events, so only enable if you know what youre doing.  
 
14. **Exit analyze function**: The analyze lambda function now completes and returns a tuple of all the relevant output variables.  
15. **Snapshot**: when `df.Snapshot` is called all the logic in `analyze` listed above is performed for every event in the dataframe. Each on it's own thread. All the variables from the output tuple from `analyze` are added as thier own columns of the dataframe. Snapshot writes out the dataframe into `TTree("WF")` in a temporary rootfile.
16. **Build event index**: the event order as it appears in the doutput TTree have been shuffled  as a result from the asyncrounous computing across threads. `BuildIndex("evt")` re-indexes the TTree based on the global event number `evt`.   
17. **Merge**:  the `WF` tree is cloned into the `outputFile` and the temp file deleted

## Important Variables

| Name            | Type                  | Default | Description                                 |
| --------------- | --------------------- | ------- | ------------------------------------------- |

## Functions & Lambdas
- **`void npsWF(run, seg, threads)`**  
  Main function: reads input file, reads in TDC offsets, enables relevant branches, enables MT, orchestrates the workflow. Variables declared here are **global** and shared on all threads!  
- **`bool FastCloneAndFilter(const TString &inName,const TString &outName)`**: 
  Clones all the TTrees, branches, and histograms from the input file to and new output file and removes the `NPS.cal.fly.adcSampWaveform` branch. Returns bool true if clone succeeds.

- **`void analyze(...)` lambda**:
  This function is what is called **per-event** when the RDataframe calls Snapshot. It is ran on it's own thread for each event and every variable in it's scope is **thread-local**. A sort of *main* function for each thread. Loops over all blocks and calls all the following functions. Returns a tuple containing all the final output variables that will be written to the dataframe in `npsWF` 
- **`void FindPulsesMF(...)`**:
      Applied the matched filter, calls TSpectrum to find peaks, and updates `wfnpulses[]`, `wfampl[]`, and `wftime[]` for the found peaks.
- **`bool PassClusterThreshold(...)`**:
    Sums raw ADC values for 8 surrounding blocks to check that the max of the sum meets 10mV threshold. Returns bool: if true → okay to fit , if false → skip block

- **`void Fitwf(evt, blocknum, err[])` lambda**:  
Where the fit occurs: Creates interpolator of reference waveform data. Initializes TF1 fit function parameters from the `FindPulsesMF(...)` output. Creates the `Minuit2` minimizer and performs a least square fit on the waveform data. Returns void and updates `chi2[]`, `wfampl[]`, and `wftime[]` for that block's fit.

## Output
- **`../nps_production_%d_%d_%d_wf.root`**, where %d are run, seg, nthreads
  - Preserves all original input TTrees & objects  
  - Adds `WF` TTree with branches:  
    

| Branch        | Type                    | Description                                                                  |
| ------------- | ----------------------- | ---------------------------------------------------------------------------- |
| evt           | `Double_t`                 | Global event number                                |
| wftime        | `RVec<Double_t>` | Pulse time of each fitted pulse  per block (ns) See Section: Notes                                          |
| wfampl        | `RVec<Double_t>` | Waveform amplitude of each fitted pulse (mV) See Section: Notes                |
| chi2          | `std::vector<Double_t>` | χ² of the fit                                               |
| pres          | `std::vector<Int_t>`                 | Presence flag: 1 if this block has data/pulses, 0 if empty                   |
| corr_time_HMS | `Double_t` | HMS time correction (ns), computed from raw TDC times and channel offsets    |
| h1time        | `std::vector<Double_t>` |                                          |
| h2time        | `std::vector<Double_t>` |                                        |
| wfnpulse      | `std::vector<Int_t>`                 | Number of pulses found by TSpectrum in this block                            |
| timewf        | `std::vector<Double_t>` | Peak time of the first fitted pulse in a block (in 4ns ADC-bin units)                                  |
| amplwf        | `std::vector<Double_t>` | Peak amplitude of the first fitted pulse in a block (mV)|
| ampl          | `std::vector<Double_t>` | The maximum value of all ADC samples in each block       (mV)                        |
| Sampampl      | `std::vector<Double_t>` |                            |
| Samptime      | `std::vector<Double_t>` |                        |
| Sampener      | `std::vector<Double_t>` |                             |
| Sampped       | `std::vector<Double_t>` |                                 |
| enertot       | `Double_t` | Total integrated raw-waveform adc values for this event within a time window of (30,109) bins (sum of samples for all blocks in units of mV)         |
| integtot      | `Double_t` | Total integrated raw-waveform adc values for this event (sum of samples for all blocks in units of mV)         |

 ## Important Notes
- **Output event ordering**:  
  ROOT's Rdataframe gives each thread multiple events to analyze in "buckets". Some thread's buckets will finish before others. This results in the output TTree's events being in a different order than they came in (i.e. different from the ordering of events in the HMS tree T). To mitigate this the TTree has been "re-indexed" via `BuildIndex("g.evnum");` where `g.evnum` is the global event number from the original input file. For the end user, all this means is that looping over events in any analysis code will need a slightly different syntax than most comonly used. An example of this is given in the next Section.
- **Flattened Output Vectors**: The 2 output vectors `wfampl` and `wftime` are of size=nblocks*npulses. Because the number of pulses and number of active blocks per event change, these are dynamic sized vectors. If 1 of 1080 blocks is not active, then values for that block are not included in the array. To effectively parse these values the user should utilize the `wfnpulse` vector, which **always** has 1080 elements indicating the number of pulses per block. So for example, if a block number 4 in wfnpulse==0, then the user should know that the values in wfampl and wftime correspond to pulses for blocks 0,1,2,3, and then will jump to pulses for block 5.
- **Failed Fits**
  With the default theshold values and minimizer configuration, it seems about 1-2% of blocks fail to fit. This can happen for a variety of reasons that are still being looked into. Most common are very noisy blocks where TSpectrum found peaks above threshold, but the minimizer never converged. The value of the threholds and minimizer configuration can be tweaked to reduce the number of fails, but this drastically increases CPU time. If a fit fails, the `wfampl[]` and `wftime[]` will be set to the values of the raw data as found by TSpectrum and the `chi2` will be set to -100 to flag these events as failed.


## Output Usage Example
A very simple example code is included on the Github. This is how a basic analysis code would read events from the TTree T and ensure they correspond to the same events in the new TTree WF.

```cpp
void basicanalysiscode(const char *fileName)
{
 //Normal boiler plate is the same
    TFile *f1 = TFile::Open(fileName, "READ");
    TTree *t1 = (TTree *)f2->Get("WF");

 //Get the global event number
    t1->SetBranchAddress("evt", &evtVal);

 //This is where things are slightly different
    TVirtualIndex *vIdx = t1->GetTreeIndex();
    TTreeIndex *idx = dynamic_cast<TTreeIndex *>(vIdx);

    Long64_t *indexArray = idx->GetIndex();

 //Now loop over all entries
    Long64_t n = t1->GetEntries();
    for (Long64_t i = 0; i < n; ++i)
    {
        // Get entrys at the indexArray[i] , not at i !
        t2->GetEntry(indexArray[i]);
  
        std::cout << std::fixed << "sorted[" << i << "] → original Entry=" << indexArray[i] << ", evt=" << evtVal << "\n";
    }
}
```


## Extending & Customizing
- **Tweaking Hardcoded Values**: There are several hardcoded values that users are free to tweak and recompile:
  - `specthres = 0.02` peaks with amplitude less than specthres*highest_peak are discarded when TSpectrum::Search() is called.
  - `mfthres = 1.5`    peaks with amplitude less than mfthres in the matched filter are discarded.
  - `trig_thres = 10`  threshold (mV) of sum of 3x3 wf to allow the fit of the central wf.
  - `coinc_width = 20`  (4ns units) the trig_thres will be applied on sum of 3x3 wf in this region around the the coin_time.
  - `cfg.SetMinimizer("Minuit2", "Migrad")` this is the minimizer selected. Users can easily change minimizer here.
  - `maxwfpulses` the maximum number of waveform peaks TSpectrum is allowed to find.

- **Add a new output branch**: To add a new output branch to the TTree WF you would need to add/change some code in the following:  
  1. Declare and initialize the variable **within** the scope of `analyze()` function. Outside this scope will not be thread-safe.
     - per-block specific variables should be vectors, RVecs, or arrays of size 1080. 
  2. Append the variable to the return tuple in `analyze()`: `return std::make_tuple(myBranch,...);`
  3. Insert the variable into the dataframe `df.Define("myBranch", [](myBranchType)){"tuple"}` at the end of `npsWF()` before the snapshot is called.
  4. Append `"myBranch"` to the field list in `RDataFrame::Snapshot()`  
- **Running on Subset of Events**: Unfortunately the function to select a subset of events from RDataframes breaks multi-threading. So if you want to run the code only for a subset of events from a root file you are currently limited to 1-thread. The following changes to the code should be made.
    1. Comment out `ROOT::EnableImplicitMT(nthreads);`
    2. Uncomment out `//auto df = df.Range(0, 5000);` and change values to your event subset.
  
- **Enable diagnostics**: set `makeDiagnostics` to `true`.
  - **Warning**: This part of the code is a crude way of printing out plots of full fitted waveforms. For just one given event this mean creating ~25 hisograms and ~25 fits functions and writing it to a file. One pdf file is created *per event*. This mode should only be activated in the single-threaded mode with only a very limited range of events.


