#include "main.h"
#include "Nano.h" // Contains the definition of the "nt" object that reads the NanoAOD
#include "MCTools.h" // Contains the definition of the dumpGenParticleInfos();
#include "rooutil.h" // Contains the definitions of various help functions that start with "RooUtil"

//=================================================================================================
// Key analysis variables
// NOTE: variables have "_" ending
//=================================================================================================
namespace Analysis
{

    //_______________________________________________________
    // Per event weight (normalized to 1/fb of data)
    Double_t scale1fb_;

    //_______________________________________________________
    // Integrated Luminosity
    Double_t lumi_;

    //_______________________________________________________
    // Year of the data
    Int_t year_;

    //_______________________________________________________
    // Flags used in the anlaysis
    bool is_vbs_wzh_;

    //_______________________________________________________
    // Signal's generator level kinematics
    LV gen_ijet_;
    LV gen_jjet_;
    LV gen_jet0_;
    LV gen_jet1_;
    LV gen_W_;
    LV gen_Z_;
    LV gen_H_;

    //_______________________________________________________
    // Reconstructed objects
    std::vector<LV> elecs_;
    std::vector<LV> muons_;
    std::vector<LV> jets_;
    std::vector<LV> bjets_;

    //_______________________________________________________
    // Set per event weight (normalized to 1/fb of data)
    void setScale1fb(Double_t scale1fb) { scale1fb_ = scale1fb; }

    //_______________________________________________________
    // Set integrated luminosity
    void setLumi(Double_t lumi) { lumi_ = lumi; }

    //_______________________________________________________
    // Set the year of the data
    void setYear(Int_t year) { year_ = year; }

    //_______________________________________________________
    // Clear all analysis variables
    void clearAnalysisVariables()
    {
        // Flags used in the anlaysis
        is_vbs_wzh_ = false;

        // Signal's generator level kinematics
        gen_ijet_ = LV();
        gen_jjet_ = LV();
        gen_jet0_ = LV();
        gen_jet1_ = LV();
        gen_W_ = LV();
        gen_Z_ = LV();
        gen_H_ = LV();

        // Reconstructed objects
        elecs_.clear();
        muons_.clear();
        jets_.clear();
        bjets_.clear();

    }

    //_______________________________________________________
    // Compute signal generator level kinematics
    void computeSignalGeneratorLevelKinematics()
    {
        is_vbs_wzh_ = nt.GenPart_status()[2] == 23;
        if (is_vbs_wzh_)
        {
            // LV ordered [NObjects];
            gen_ijet_ = nt.GenPart_p4()[2];
            gen_jjet_ = nt.GenPart_p4()[3];
            gen_jet0_ = gen_ijet_.pt() > gen_jjet_.pt() ? gen_ijet_ : gen_jjet_;
            gen_jet1_ = gen_ijet_.pt() > gen_jjet_.pt() ? gen_jjet_ : gen_ijet_;
            gen_W_ = nt.GenPart_p4()[4];
            gen_Z_ = nt.GenPart_p4()[5];
            gen_H_ = nt.GenPart_p4()[6];
        }
    }

}

//=================================================================================================
// Histograms
// Naming convention: "h_NAME_" (i.e. Starts with "h_" and ends with "_")
//=================================================================================================
namespace Hist
{
    // pt
    TH1F* h_gen_ptZ_;
    TH1F* h_gen_ptW_;
    TH1F* h_gen_ptH_;
    // VBF quarks
    TH1F* h_gen_deltaEta_;
    TH1F* h_gen_massVBF_;

    //_______________________________________________________
    // Book the histograms
    void bookHistograms()
    {
        // Pt's of the bosons
        h_gen_ptZ_ = new TH1F("ptZ", "Transverse momentum of Z bosons", 1080, 0, 1800);
        h_gen_ptW_ = new TH1F("ptW", "Transverse momentum of W bosons", 1080, 0, 1800);
        h_gen_ptH_ = new TH1F("ptH", "Transverse momentum of H bosons", 1080, 0, 1800);
        // VBF quarks kinematics
        h_gen_deltaEta_ = new TH1F("deltaEta", "Delta Eta of VBF quarks", 1080, 0, 10);
        h_gen_massVBF_ = new TH1F("massVBF", "Invariant mass of VBF quark system", 1080, 0, 3500);
    }

    //_______________________________________________________
    // Fill the histograms
    void fillHistograms()
    {
        // There are some cases where the event weight is negative
        Double_t sign_of_the_weight = ((nt.Generator_weight() > 0) - (nt.Generator_weight() < 0));
        Double_t per_event_weight = sign_of_the_weight * Analysis::scale1fb_ * Analysis::lumi_;

        // Fill histograms
        h_gen_massVBF_->Fill((Analysis::gen_ijet_ + Analysis::gen_jjet_).M(), per_event_weight);
        h_gen_deltaEta_->Fill(TMath::Abs(Analysis::gen_ijet_.Eta() - Analysis::gen_jjet_.Eta()), per_event_weight);
        h_gen_ptZ_->Fill(Analysis::gen_Z_.Pt(), per_event_weight);
        h_gen_ptW_->Fill(Analysis::gen_W_.Pt(), per_event_weight);
        h_gen_ptH_->Fill(Analysis::gen_H_.Pt(), per_event_weight);
    }

    //_______________________________________________________
    // Write the histograms
    void writeHistograms(TFile* ofile)
    {
        // Done with the analysis
        ofile->cd();
        h_gen_massVBF_->Write();
        h_gen_deltaEta_->Write();
        h_gen_ptH_->Write();
        h_gen_ptW_->Write();
        h_gen_ptZ_->Write();

    }

}


//=================================================================================================
// Main function
//=================================================================================================
int main()
{
    // Create your output root file
    TFile* output_file = new TFile("/home/users/joytzphysics/NanoLooper/NanoLooper/wzh.root", "recreate");

    // Create Histograms
    Hist::bookHistograms();

    // Set scale 1fb (the per event weight normalized for 1/fb)
    Analysis::setScale1fb(0.018);

    // Set the luminosity
    Analysis::setLumi(137); // TODO: Update properly in the future. For now it's a Placeholder!

    // Set the luminosity
    Analysis::setYear(2018); // TODO: Update properly in the future. For now it's a Placeholder!

    // Input Path (RooUtil::Looper can accept comma separated list)
    std::vector<TString> input_file_paths;
    for (unsigned i = 0; i < 26; i++)
    {
        TString input_temp = "/hadoop/cms/store/user/phchang/VBSHWWSignalGeneration/UL18_VBSWZH_incl_C2V_4_Azure_v1/merged/output_" + std::to_string(i) + ".root";
        input_file_paths.push_back(input_temp);
    }

    // Use RooUtil::StringUtil to concatenate the input file paths
    TString input = RooUtil::StringUtil::join(input_file_paths, ",");

    // Name of the TTree in the input root files.
    TString input_tree_name = "Events";

    // Create the TChain that holds the TTree's of the NanoAOD ntuples.
    TChain* events_tchain = RooUtil::FileUtil::createTChain(input_tree_name, input);

    // Max number of events to loop. (-1 = loop over all)
    // int n_events = 1; // Looping over first event only
    int n_events = -1; // To loop over all events

    // Create a looper that will handle event looping
    RooUtil::Looper<Nano> looper;

    nt.SetYear(Analysis::year_);
    // Initializer the looper
    looper.init(events_tchain, &nt, n_events);


    // Loop through events
    while (looper.nextEvent())
    {

        // Beginning of every event we should clear the analysis variables
        Analysis::clearAnalysisVariables();

        // Compute signal generator level kinematics
        Analysis::computeSignalGeneratorLevelKinematics();

        // Fill the histograms
        Hist::fillHistograms();


    }

    // Write out the histograms
    Hist::writeHistograms(output_file);

    // Close the file
    output_file->Close();


    return 0;
}
