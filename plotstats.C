#include <TFile.h>
#include <TTree.h>
#include <iostream>

void plotstats(const char *fileName)
{
    TFile *f2 = TFile::Open(fileName, "READ");
    TTree *t2 = (TTree *)f2->Get("WF");
    TTree *t1 = (TTree *)f2->Get("T");

    Double_t evtVal = 0;
    Double_t evtVal1 = 0;
    t2->SetBranchAddress("evt", &evtVal);
    t1->SetBranchAddress("g.evnum", &evtVal1);


    TVirtualIndex *vIdx = t2->GetTreeIndex();
    TTreeIndex *idx = dynamic_cast<TTreeIndex *>(vIdx);
    if (!idx)
    {
        std::cerr << "No TTreeIndex found in the output tree.\n";
        f2->Close();
        return;
    }

    // Get the index array (Long64_t*), not Int_t*
    Long64_t *indexArray = idx->GetIndex();
    Long64_t n = t2->GetEntries();
    Double_t lastevent = 0.;

    for (Long64_t i = 0; i < n; ++i)
    {
        // indexArray[i] is the original-entry number for the i-th smallest evt
        t2->GetEntry(indexArray[i]);
        t1->GetEntry(i);
        if (true)
        {
            //std::cout << std::fixed << "sorted[" << i << "] → original Entry=" << indexArray[i] << ", evt=" << evtVal << "\n";
                        std::cout << std::fixed << "sorted[" << i << "] → original Entry=" << indexArray[i] << ", evt=" << evtVal << " From T TREE "<< ", evt=" << evtVal1 << "\n";

        }

        if (i > 0 && evtVal != (lastevent + 1.0))
        {
            cout << "WRONG " << evtVal << " != " << (lastevent+1.0) << endl;
        }
        lastevent = evtVal;
    }
}