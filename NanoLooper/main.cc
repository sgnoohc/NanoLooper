#include "main.h"
#include "Nano.h" // Contains the definition of the "nt" object that reads the NanoAOD
#include "rooutil.h" // Contains the definitions of various help functions that start with "RooUtil"

int main()
{

    // Input file path (comma separated string).
    TString input_file_list_tstring = "/hadoop/cms/store/user/phchang/VBSHWWSignalGeneration/VBSWWH_C2V_4p5_RunIIAutumn18NanoAOD_VBSWWH_C2V_4p5_v3_ext1/merged/output.root";

    // Name of the TTree in the input root files.
    TString input_tree_name = "Events";

    // Create the TChain that holds the TTree's of the NanoAOD ntuples.
    TChain* events_tchain = RooUtil::FileUtil::createTChain(input_tree_name, input_file_list_tstring);

    // Max number of events to loop. (-1 = loop over all)
    int n_events = -1;

    // Create a looper that will handle event looping
    RooUtil::Looper<Nano> looper;

    // Initializer the looper
    looper.init(events_tchain, &nt, n_events);

    // Create your output root file
    // TFile* output_file = new TFile("output.root", "recreate");

    // Create your histograms
    // TH1F* ...
    // TH1F* ...
    // TH1F* ...
    // TH1F* ...

    while (looper.nextEvent())
    {
        // Do your analysis here
        // float mjj = blah blah;
    }

    // Done with the analysis
    // output_file->cd();
    // hist->Write();
    // ...
    // ...
    // ...

    return 0;
}
