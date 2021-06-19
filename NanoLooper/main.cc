#include "main.h"
#include "Nano.h" // Contains the definition of the "nt" object that reads the NanoAOD
#include "MCTools.h" // Contains the definition of the dumpGenParticleInfos();
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
    int n_events = 1; // Looping over first event only
    // int n_events = -1; // To loop over all events

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

        // Accessing the generator level information through branches with names that start with "GenPart_"

        // To access the TTree content we use a "nt" object which is already taking care of various "SetBranchAddress" behind the scenes.
        // By simply calling nt instance's methods with the same name as the TTree branch names, one can access the content of the TTree.
        for (int ipart = 0; ipart < nt.GenPart_pt().size(); ++ipart)
        {
            std::cout << "ipart = " << ipart << " Transverse momentum (pt)    = " << nt.GenPart_pt()[ipart] << std::endl;
            std::cout << "ipart = " << ipart << " Pseudorapditiy      (eta)   = " << nt.GenPart_eta()[ipart] << std::endl;
            std::cout << "ipart = " << ipart << " Azimuthal angle     (phi)   = " << nt.GenPart_phi()[ipart] << std::endl;
            std::cout << "ipart = " << ipart << " PDGID               (pdgid) = " << nt.GenPart_pdgId()[ipart] << std::endl;
        }

        // Alternatively, one can print the same (and more) information using tools developed by our research group's (and colleagues).
        dumpGenParticleInfos();

    }

    // Done with the analysis
    // output_file->cd();
    // hist->Write();
    // ...
    // ...
    // ...

    return 0;
}
