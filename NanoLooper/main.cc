#include "main.h"
#include "Nano.h" // Contains the definition of the "nt" object that reads the NanoAOD
#include "MCTools.h" // Contains the definition of the dumpGenParticleInfos();
#include "ElectronSelections.h" // Contains the definitions of ttH::electronID
#include "MuonSelections.h" // Contains the definitions of ttH::muonID
#include "Config.h" // Contains the definitions of gconf
#include "rooutil.h" // Contains the definitions of various help functions that start with "RooUtil"

//=================================================================================================
// Objects
//=================================================================================================
namespace Obj
{
    //_______________________________________________________
    // Electron data structure
    struct Elec
    {
        LV p4;
        int isTight;
        int jetIdx; // Overlapping jet index
    };

    //_______________________________________________________
    // Muon data structure
    struct Muon
    {
        LV p4;
        int isTight;
        int jetIdx; // Overlapping jet index
    };

    //_______________________________________________________
    // Jet data structure
    struct Jet
    {
        LV p4;
        float isBtagScore;
        int isBtagTight;
        int isBtagLoose;
    };
}

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
    // Per event weight (normalized * lumi)
    Double_t wgt_;

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
    std::vector<Obj::Elec> elecs_;
    std::vector<Obj::Muon> muons_;
    std::vector<Obj::Jet> jets_;

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
    // Set up the NanoCORE's common configuration service tool
    void setConfig()
    {
        gconf.GetConfigs(nt.year());
    }

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

    }

    //_______________________________________________________
    // Compute the event weight for this event
    void computeEventWeight()
    {
        // There are some cases where the event weight is negative
        Double_t sign_of_the_weight = ((nt.Generator_weight() > 0) - (nt.Generator_weight() < 0));
        Double_t per_event_weight = sign_of_the_weight * Analysis::scale1fb_ * Analysis::lumi_;
        wgt_ = per_event_weight;
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

    //_______________________________________________________
    // Select electrons
    void selectElectrons()
    {
        // Select electrons
        for (unsigned int iel = 0; iel < nt.Electron_pt().size(); ++iel)
        {
            if (ttH::electronID(iel, ttH::IDveto, nt.year()))
            {
                Obj::Elec this_elec;
                this_elec.p4 = nt.Electron_p4()[iel];
                this_elec.isTight = ttH::electronID(iel, ttH::IDtight, nt.year());
                this_elec.jetIdx = nt.Electron_jetIdx()[iel];
                elecs_.push_back(this_elec);
            }
        }
    }

    //_______________________________________________________
    // Select muons
    void selectMuons()
    {
        // Select muons
        for (unsigned int imu = 0; imu < nt.Muon_pt().size(); ++imu)
        {
            if (ttH::muonID(imu, ttH::IDveto, nt.year()))
            {
                Obj::Muon this_muon;
                this_muon.p4 = nt.Muon_p4()[imu];
                this_muon.isTight = ttH::muonID(imu, ttH::IDtight, nt.year());
                this_muon.jetIdx = nt.Muon_jetIdx()[imu];
                muons_.push_back(this_muon);
            }
        }
    }

    //_______________________________________________________
    // Select jets
    void selectJets()
    {
        // Loop over the jets
        for (unsigned int ijet = 0; ijet < nt.Jet_pt().size(); ++ijet)
        {

            // Read jet p4
            const LV& jet_p4 = nt.Jet_p4()[ijet];

            if (nt.year() == 2016)
            {
                if (nt.Jet_jetId()[ijet] < 1) // For 2016 apparently it's >= 1
                    continue;
            }
            else
            {
                if (nt.Jet_jetId()[ijet] < 2) // "Tight" ID requirement while for others >= 2
                    continue;
            }

            // if (nt.year() == 2017)
            // {
            //     if (nt.Jet_puId()[ijet] < 7) // For 2017 "111" (pass tight)
            //         continue;
            // }

            // Overlap check against good electrons
            bool isOverlap = false;
            for (unsigned int ilep = 0; ilep < elecs_.size(); ++ilep)
            {
                if (elecs_[ilep].jetIdx == (int) ijet)
                {
                    isOverlap = true;
                    break;
                }
            }

            // Overlap check against good muons
            for (unsigned int ilep = 0; ilep < muons_.size(); ++ilep)
            {
                if (muons_[ilep].jetIdx == (int) ijet)
                {
                    isOverlap = true;
                    break;
                }
            }

            // Then skip
            if (isOverlap)
                continue;

            // We keep jets above 20 GeV only
            if (not (jet_p4.pt() > 20.))
                continue;

            Obj::Jet this_jet;
            this_jet.p4 = nt.Jet_p4()[ijet];
            this_jet.isBtagScore = nt.Jet_btagDeepFlavB()[ijet];
            this_jet.isBtagTight = nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_tight;
            this_jet.isBtagLoose = nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_loose;
            jets_.push_back(this_jet);

        }

    }

    //_______________________________________________________
    // Run the analysis for this event
    void runAlgorithms()
    {
        // Beginning of every event we should clear the analysis variables
        Analysis::clearAnalysisVariables();

        // Compute per event weight
        Analysis::computeEventWeight();

        // Compute signal generator level kinematics
        Analysis::computeSignalGeneratorLevelKinematics();

        // Select electrons and muons
        Analysis::selectElectrons();
        Analysis::selectMuons();
        Analysis::selectJets();

    }

}

//=================================================================================================
// Cutflows
// Naming convention: "h_NAME_" (i.e. Starts with "h_" and ends with "_")
//=================================================================================================
namespace Cutflow
{
    // cutflow
    TH1F* h_cutflow_;

    //_______________________________________________________
    // Cutflow
    enum Cuts
    {
        kNoSelection = 0,
        kTwoLightLeptons,
        kAtLeastFourPt20Jets,
        kNCuts,
    };

    //_______________________________________________________
    // Book the histograms
    void bookCutflow()
    {
        // Cutflow histogram
        h_cutflow_ = new TH1F("h_cutflow", "Cutflow of the selections", Cuts::kNCuts, 0, Cuts::kNCuts);
    }

    //_______________________________________________________
    // Fill cutflow
    void fillCutflow(Cuts cutstage)
    {
        h_cutflow_->Fill(cutstage, Analysis::wgt_);
    } 

    //_______________________________________________________
    // write cutflow
    void writeCutflow(TFile* ofile)
    {
        ofile->cd();
        h_cutflow_->Write();
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
        h_gen_ptZ_ = new TH1F("h_gen_ptZ", "Transverse momentum of Z bosons", 1080, 0, 1800);
        h_gen_ptW_ = new TH1F("h_gen_ptW", "Transverse momentum of W bosons", 1080, 0, 1800);
        h_gen_ptH_ = new TH1F("h_gen_ptH", "Transverse momentum of H bosons", 1080, 0, 1800);
        // VBF quarks kinematics
        h_gen_deltaEta_ = new TH1F("h_gen_deltaEta", "Delta Eta of VBF quarks", 1080, 0, 10);
        h_gen_massVBF_ = new TH1F("h_gen_massVBF", "Invariant mass of VBF quark system", 1080, 0, 3500);
    }

    //_______________________________________________________
    // Fill the histograms
    void fillHistograms()
    {
        // Fill histograms
        h_gen_massVBF_->Fill((Analysis::gen_ijet_ + Analysis::gen_jjet_).M(), Analysis::wgt_);
        h_gen_deltaEta_->Fill(TMath::Abs(Analysis::gen_ijet_.Eta() - Analysis::gen_jjet_.Eta()), Analysis::wgt_);
        h_gen_ptZ_->Fill(Analysis::gen_Z_.Pt(), Analysis::wgt_);
        h_gen_ptW_->Fill(Analysis::gen_W_.Pt(), Analysis::wgt_);
        h_gen_ptH_->Fill(Analysis::gen_H_.Pt(), Analysis::wgt_);
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

    // Create Cutflow Histogram
    Cutflow::bookCutflow();

    // Set scale 1fb (the per event weight normalized for 1/fb)
    Analysis::setScale1fb(0.018);

    // Set the luminosity
    Analysis::setLumi(137); // TODO: Update properly in the future. For now it's a Placeholder!

    // Set the luminosity
    Analysis::setYear(2018); // TODO: Update properly in the future. For now it's a Placeholder!

    // Set up year dependent configuration of the analysis that are POG specific from CMS
    Analysis::setConfig();

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

    // TODO: Update properly in the future. For now, this is an AdHoc solution.
    nt.SetYear(Analysis::year_);

    // Initializer the looper
    looper.init(events_tchain, &nt, n_events);

    // Loop through events
    while (looper.nextEvent())
    {

        // Run the analysis algorithms (selections, computing variables, etc.)
        Analysis::runAlgorithms();

        // Fill the "counter" histogram
        Cutflow::fillCutflow(Cutflow::Cuts::kNoSelection);

        // Cut #1: Require that there are exactly two leptons
        if (not (Analysis::elecs_.size() + Analysis::muons_.size() == 2)) { continue; }

        // Fill the "counter" histogram
        Cutflow::fillCutflow(Cutflow::Cuts::kTwoLightLeptons);

        // Cut #2: Require that there are at least 4 pt > 20 GeV jets
        if (not (Analysis::jets_.size() >= 4)) { continue; }

        // Fill the "counter" histogram
        Cutflow::fillCutflow(Cutflow::Cuts::kAtLeastFourPt20Jets);

        // All cuts have passed
        // Now fill the histograms
        Hist::fillHistograms();

    }

    // Write out the histograms
    Hist::writeHistograms(output_file);

    // Write out the cutflow histogram
    Cutflow::writeCutflow(output_file);

    // Close the file
    output_file->Close();


    return 0;
}
