#include "main.h"
#include "Nano.h" // Contains the definition of the "nt" object that reads the NanoAOD
#include "MCTools.h" // Contains the definition of the dumpGenParticleInfos();
#include "rooutil.h" // Contains the definitions of various help functions that start with "RooUtil"
int main()
{
    // Create your output root file
    TFile* output_file = new TFile("/home/users/joytzphysics/NanoLooper/NanoLooper/wzh.root", "recreate");
    // pt
    TH1F* ptZ = new TH1F("ptZ", "Transverse momentum of Z bosons", 1080, 0, 1800);
    TH1F* ptW = new TH1F("ptW", "Transverse momentum of W bosons", 1080, 0, 1800);
    TH1F* ptH = new TH1F("ptH", "Transverse momentum of H bosons", 1080, 0, 1800);
    // VBF quarks
    TH1F* deltaEta = new TH1F("deltaEta", "Delta Eta of VBF quarks", 1080, 0, 10);
    TH1F* massVBF = new TH1F("massVBF", "Invariant mass of VBF quark system", 1080, 0, 3500);
    // Cross section
    Double_t cs = 0.018 * 137;

    // Input Path
    for (unsigned i = 0; i < 26; i++)
    {
        TString input = "/hadoop/cms/store/user/phchang/VBSHWWSignalGeneration/UL18_VBSWZH_incl_C2V_4_Azure_v1/merged/output_" + std::to_string(i) + ".root";

        // Name of the TTree in the input root files.
        TString input_tree_name = "Events";

        // Create the TChain that holds the TTree's of the NanoAOD ntuples.
        TChain* events_tchain = RooUtil::FileUtil::createTChain(input_tree_name, input);

        // Max number of events to loop. (-1 = loop over all)
        // int n_events = 1; // Looping over first event only
        int n_events = -1; // To loop over all events

        // Create a looper that will handle event looping
        RooUtil::Looper<Nano> looper;

        nt.SetYear(2018);
        // Initializer the looper
        looper.init(events_tchain, &nt, n_events);


        // Loop through events
        while (looper.nextEvent())
        {
            // Do your analysis here
            // float mjj = blah blah;

            // Select Leptons for the event

            bool isvbswzh = nt.GenPart_status()[2] == 23;

            if (isvbswzh)
            {
                // LV ordered [NObjects];
                const LV& ijet = nt.GenPart_p4()[2];
                const LV& jjet = nt.GenPart_p4()[3];
                const LV& jet0 = ijet.pt() > jjet.pt() ? ijet : jjet;
                const LV& jet1 = ijet.pt() > jjet.pt() ? jjet : ijet;
                const LV& W = nt.GenPart_p4()[4];
                const LV& Z = nt.GenPart_p4()[5];
                const LV& H = nt.GenPart_p4()[6];
                // Fill histograms
                massVBF->Fill((ijet + jjet).M());
                deltaEta->Fill(TMath::Abs(ijet.Eta()-jjet.Eta()));
                ptZ->Fill(Z.Pt());
                ptW->Fill(W.Pt());
                ptH->Fill(H.Pt());
            }
        }
        // Alternatively, one can print the same (and more) information using tools developed by our research group's (and colleagues).
        // dumpGenParticleInfos();

    }

    // Done with the analysis
    output_file->cd();
    massVBF->Scale(cs / (massVBF->Integral()));
    massVBF->Write();
    deltaEta->Scale(cs / (deltaEta->Integral()));
    deltaEta->Write();
    ptH->Scale(cs / (ptH->Integral()));
    ptH->Write();
    ptW->Scale(cs / (ptW->Integral()));
    ptW->Write();
    ptZ->Scale(cs / (ptZ->Integral()));
    ptZ->Write();

    output_file->Close();


    return 0;
}
