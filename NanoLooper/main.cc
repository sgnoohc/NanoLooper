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
        int pdgid;
    };

    //_______________________________________________________
    // Muon data structure
    struct Muon
    {
        LV p4;
        int isTight;
        int jetIdx; // Overlapping jet index
        int pdgid;
    };


    //_______________________________________________________
    // Jet data structure
    struct Jet
    {
        LV p4;
        int JetIdx; // Overlapping jet index
        float hbbScore;
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
    // Flags used in the analysis
    bool is_vbs;
    bool is_hbb;
//    bool is_wjj;

    //_______________________________________________________
    // Signal's generator level kinematics
    LV gen_ijet_;
    LV gen_jjet_;
    LV gen_jet0_;
    LV gen_jet1_;
    std::vector<LV> gen_W_;
    std::vector<LV> gen_Z_;
    LV gen_H_;
    LV gen_b0_;
    LV gen_b1_;
    LV gen_Wb0_;
    LV gen_Wb1_;

    //_______________________________________________________
    // Reconstructed objects
    std::vector<Obj::Elec> elecs_;
    std::vector<Obj::Muon> muons_;
    std::vector<Obj::Jet> jets_;
   std::vector<LV> leptons_;
    std::vector<LV> VBFjets_;
    std::vector<LV> bjets_;
    std::vector<LV> Wjets_;

    //________________________________________________________
    // Fat Jets reconstruction
    std::vector<Obj::Jet> fatJets_;
    LV hbbFatJet_;
 

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
        is_vbs = false;
        is_hbb = false;
//        is_wjj = false;

        // Signal's generator level kinematics
        gen_ijet_ = LV();
        gen_jjet_ = LV();
        gen_jet0_ = LV();
        gen_jet1_ = LV();
        gen_H_ = LV();
        gen_b0_ = LV();
        gen_b1_ = LV();
        gen_Wb0_ = LV();
        gen_Wb1_ = LV();
        

        // Reconstructed objects
        gen_W_.clear();
	gen_Z_.clear();
        elecs_.clear();
        muons_.clear();
        jets_.clear();
        fatJets_.clear();
        leptons_.clear();
        VBFjets_.clear();
        bjets_.clear();
        Wjets_.clear();
        hbbFatJet_ = LV(); 

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
        is_vbs = nt.GenPart_status()[2] == 23;
        if (is_vbs)
        {
            // LV ordered [NObjects];
            gen_ijet_ = nt.GenPart_p4()[2];
            gen_jjet_ = nt.GenPart_p4()[3];
            gen_jet0_ = gen_ijet_.pt() > gen_jjet_.pt() ? gen_ijet_ : gen_jjet_;
            gen_jet1_ = gen_ijet_.pt() > gen_jjet_.pt() ? gen_jjet_ : gen_ijet_;
	    for (unsigned int i = 0; i < 7; ++i) {
		if (abs(nt.GenPart_pdgId()[i] == 24)) { gen_W_.push_back(nt.GenPart_p4()[i]);}
		if (abs(nt.GenPart_pdgId()[i] == 23)) { gen_Z_.push_back(nt.GenPart_p4()[i]);}
		if (abs(nt.GenPart_pdgId()[i] == 25)) { gen_H_ = nt.GenPart_p4()[i];} 
	    }
            // Select b quarks and hbb
            std::vector<LV> h_decay_p4s;
            std::vector<int> h_decay_pdgIds;
            // std::vector<int> h_decay_statuses;
            for (unsigned int igen = 0; igen < nt.GenPart_p4().size(); ++igen)
            {
                int imom = nt.GenPart_genPartIdxMother()[igen];
                if (abs(nt.GenPart_pdgId()[imom]) == 25 and nt.GenPart_status()[imom] == 62)
                {
                    // int pdgid = nt.GenPart_pdgId()[igen];
                    // int status = nt.GenPart_status()[igen];
                    h_decay_p4s.push_back(nt.GenPart_p4()[igen]);
                    h_decay_pdgIds.push_back(nt.GenPart_pdgId()[igen]);
                    // h_decay_statuses.push_back(nt.GenPart_status()[igen]);
                }
            }
            if (abs(h_decay_pdgIds[0]) == 5)
            {
                is_hbb = true;
            }
            if (is_hbb) 
            {
                gen_b0_ = h_decay_p4s[0].pt() > h_decay_p4s[1].pt() ? h_decay_p4s[0] : h_decay_p4s[1];
                gen_b1_ = h_decay_p4s[0].pt() > h_decay_p4s[1].pt() ? h_decay_p4s[1] : h_decay_p4s[0];
            }

            // Select jet quarks and wjj
            std::vector<LV> w_decay_p4s;
            std::vector<int> w_decay_pdgIds;
            // std::vector<int> w_decay_statuses;
            for (unsigned int igen = 0; igen < nt.GenPart_p4().size(); ++igen)
            {
                int imom = nt.GenPart_genPartIdxMother()[igen];
                if (abs(nt.GenPart_pdgId()[imom]) == 24 and nt.GenPart_status()[imom] == 62)
                {
                    // int pdgid = nt.GenPart_pdgId()[igen];
                    // int status = nt.GenPart_status()[igen];
                    w_decay_p4s.push_back(nt.GenPart_p4()[igen]);
                    w_decay_pdgIds.push_back(nt.GenPart_pdgId()[igen]);
                    // w_decay_statuses.push_back(nt.GenPart_status()[igen]);
                }
            }
//            if (abs(w_decay_pdgIds[0]) >= 1 and abs(w_decay_pdgIds[0]) <= 5 and (w_decay_pdgIds[0] * w_decay_pdgIds[1]) < 0){
//                is_wjj = true;
//            }
//            if (is_wjj) 
//            {
//                if (w_decay_p4s.size() == 2) {
//                    gen_Wb0_ = w_decay_p4s[0].pt() > w_decay_p4s[1].pt() ? w_decay_p4s[0] : w_decay_p4s[1];
//                    gen_Wb1_ = w_decay_p4s[0].pt() > w_decay_p4s[1].pt() ? w_decay_p4s[1] : w_decay_p4s[0];
//                }
//            }
//        }
    }

    

    //_______________________________________________________
    // Select electrons
    void selectElectrons()
    {
        for (unsigned int iel = 0; iel < nt.Electron_pt().size(); ++iel)
        {
            if (ttH::electronID(iel, ttH::IDveto, nt.year()))
            {
                Obj::Elec this_elec;
                this_elec.p4 = nt.Electron_p4()[iel];
                this_elec.isTight = ttH::electronID(iel, ttH::IDtight, nt.year());
                this_elec.jetIdx = nt.Electron_jetIdx()[iel];
                this_elec.pdgid = nt.Electron_pdgId()[iel];
                elecs_.push_back(this_elec);
                leptons_.push_back(this_elec.p4);
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
                this_muon.pdgid = nt.Muon_pdgId()[imu];
                muons_.push_back(this_muon);
                leptons_.push_back(this_muon.p4);
            }
        }
    }

    //_______________________________________________________
    // Select Fatjets
    void selectFatJets()
    {
        for (unsigned int ifatjet = 0; ifatjet < nt.FatJet_pt().size(); ++ifatjet)
        {
            Obj::Jet this_fatJet;
            // Overlap check against good electrons
            bool isOverlap = false;
            for (unsigned int ilep = 0; ilep < elecs_.size(); ++ilep)
            {
                if (RooUtil::Calc::DeltaR(elecs_[ilep].p4, nt.FatJet_p4()[ifatjet]) < 0.8)
                {
                    isOverlap = true;
                    break;
                }
            }

             // Overlap check against good muons
            for (unsigned int ilep = 0; ilep < muons_.size(); ++ilep)
            {
                if (RooUtil::Calc::DeltaR(muons_[ilep].p4, nt.FatJet_p4()[ifatjet]) < 0.8)
                {
                    isOverlap = true;
                    break;
                }
            }

            // Overlap 
            // Then skip
            if (isOverlap)
                continue;

            // Keep fat jets softdropmass above 40 GeV
            float boostedMass = nt.FatJet_msoftdrop()[ifatjet]; 
            if (not (boostedMass > 40.))
                continue;

            // p4
            TLorentzVector p4;
            p4.SetPtEtaPhiM(nt.FatJet_pt()[ifatjet], nt.FatJet_eta()[ifatjet], nt.FatJet_phi()[ifatjet], 
                            nt.FatJet_msoftdrop()[ifatjet]);
            this_fatJet.p4 = RooUtil::Calc::getLV(p4);
            this_fatJet.hbbScore = nt.FatJet_deepTagMD_HbbvsQCD()[ifatjet];
            this_fatJet.JetIdx = -999;
            fatJets_.push_back(this_fatJet);
        }
        
        // FatJet Selection
        float maxhbbscore = fatJets_[0].hbbScore;
        int maxHbbNo = 0;
        for (unsigned int ifatjet = 1; ifatjet < fatJets_.size(); ifatjet++) {
            if (fatJets_[ifatjet].hbbScore >= maxhbbscore) {
                maxHbbNo = ifatjet;
                maxhbbscore = fatJets_[ifatjet].hbbScore;
            }
        }
        hbbFatJet_ = fatJets_[maxHbbNo].p4;
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

            // Overlap check against good fatjets
            for (unsigned int ilep = 0; ilep < fatJets_.size(); ++ilep) {
                if (RooUtil::Calc::DeltaR(fatJets_[ilep].p4, jet_p4) <= 0.8) {
                    isOverlap = true;
                    break;
                }
            }

            // Overlap 
            // Then skip
            if (isOverlap)
                continue;

            // We keep jets above 30 GeV only
            if (not (jet_p4.pt() > 30.))
                continue;

            Obj::Jet this_jet;
            this_jet.p4 = nt.Jet_p4()[ijet];
            this_jet.isBtagScore = nt.Jet_btagDeepFlavB()[ijet];
            this_jet.isBtagTight = nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_tight;
            this_jet.isBtagLoose = nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_loose;
            jets_.push_back(this_jet);

        }
    }

        
    //___________________________________________________________
    // Select jet quarks
    void classifyJet()
    {
        double largestEta1 = (jets_[0].p4).Eta();
        double largestEta2 = (jets_[0].p4).Eta();
        int VBFIndex1 = 0;
        int VBFIndex2 = 0;
        for (unsigned int i = 0; i < jets_.size(); i++) {
            LV jeti = jets_[i].p4;
            if (jeti.Eta() > 0 && jeti.Eta() >= largestEta1) {
                VBFIndex1 = i;
                largestEta1 = jeti.Eta();
            }
            else if (jeti.Eta() < largestEta2 && jeti.Eta() < 0) {
                VBFIndex2 = i;
                largestEta2 = jeti.Eta();
            }
        } 
        VBFjets_.push_back(jets_[VBFIndex1].p4);
        VBFjets_.push_back(jets_[VBFIndex2].p4);
        double largestpT = 0;
        double largerpT = 0;
        int Wjet1 = 0;
        int Wjet2 = 0;
        for (unsigned int i = 0; i <jets_.size(); i++) {
            if (i != VBFIndex1 && i != VBFIndex2) {
                largestpT = (jets_[i].p4).Pt();
                largerpT = (jets_[i].p4).Pt();
                Wjet1 = i;
                Wjet2 = i;
                break;
            }
        }
        for (unsigned int j = 0; j <jets_.size(); j++) {
            if (j != VBFIndex1 && j != VBFIndex2) {
                double jetpT = (jets_[j].p4).Pt();
                if (jetpT >= largestpT) {
                    largerpT = largestpT;
                    Wjet2 = Wjet1;
                    largestpT = jetpT;
                    Wjet1 = j;
                }
                else if (jetpT < largestpT && jetpT > largerpT) {
                    Wjet2 = j;
                    largerpT = jetpT;
                }
            }
        }
        Wjets_.push_back(jets_[Wjet1].p4);
        Wjets_.push_back(jets_[Wjet2].p4);
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
        Analysis::selectFatJets();
        Analysis::selectJets();
        Analysis::classifyJet();

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
        kVbs,
        kHbb,
        kTwoLightLeptons,
        kOneHbbFatJet,
	kHbbScore,
        kAtLeastTwoPt30Jets,
	kMjj,
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

    // b quarks
    TH1F* h_gen_massbQsystem_;
    TH1F* h_gen_ptb0_;
    TH1F* h_gen_ptb1_;

    // leptons_
    TH1F* ptLep1;
    TH1F* ptLep2;
    TH1F* re_deltaEtaVBF;
    TH1F* etaLep1;
    TH1F* etaLep2;
    TH1F* phiLep1;
    TH1F* phiLep2;
    TH1F* massDiLep;
    TH1F* ptDiLep;
    TH1F* dRLep; 

    // jets
    TH1F* etaVBFjet1;
    TH1F* ptVBFjet1;
    TH1F* ptVBFjet2;
    TH1F* etaVBFjet2;
    TH1F* VBFjetMass;
    TH1F* ptB1;
    TH1F* ptB2;
    TH1F* BjetMass;
    TH1F* BjetMass_400cut;
    TH1F* ptWjet1;
    TH1F* ptWjet2;
    TH1F* WjetMass;
    
    // dR between reconstructed and gen-level
    TH1F* dRleadingVBF;
    TH1F* dRsubleadingVBF;
    TH1F* dRsubleadingbq;
    TH1F* dRleadingbq;
    TH1F* dRhiggs;

    // Recombination of VBF
    TH1F* recombdRleadingVBF;
    TH1F* recombdRleadingbq;

    // FatJets
    TH1F* ptHbb;
    TH1F* etaHbb;
    TH1F* massHbb;
    TH1F* nFatjets;

    // s-hat variable
    TH1F* massZH;
    TH1F* massZHzoom;

    //_______________________________________________________
    // Book the histograms
    void bookHistograms()
    {
        // Recombination of VBF
        recombdRleadingVBF = new TH1F("recombdRleadingVBF", "smallest delta R of leading VBF quarks between reco- and gen-level", 1080, 0, 5);
        recombdRleadingbq = new TH1F("recombdRleadingbq", "smallest delta R of leading bottom quarks between reco- and gen-level", 1080, 0, 5);
        // dR
        dRleadingVBF = new TH1F("dRleadingVBF", "delta R of leading VBF quarks", 1080, 0, 5);
        dRsubleadingbq = new TH1F("dRsubleadingbq", "delta R of subleading bottom quarks", 1080, 0, 5);
        dRsubleadingVBF = new TH1F("dRsubleadingVBF", "delta R of subleading VBF quarks", 1080, 0, 5);
        dRleadingbq = new TH1F("dRleadingbq", "delta R of leading bottom quarks", 1080, 0, 5);
        dRhiggs = new TH1F("dRhiggs", "delta R of gen higgs and tagged reco hbbfatjet", 1080, 0, 5);
        // Pt's of the bosons
        h_gen_ptZ_ = new TH1F("h_gen_ptZ", "Transverse momentum of Z bosons", 1080, 0, 1800);
        h_gen_ptW_ = new TH1F("h_gen_ptW", "Transverse momentum of W bosons", 1080, 0, 1800);
        h_gen_ptH_ = new TH1F("h_gen_ptH", "Transverse momentum of H bosons", 1080, 0, 1800);
        // VBF quarks kinematics
        h_gen_deltaEta_ = new TH1F("h_gen_deltaEta", "Delta Eta of VBF quarks", 1080, 0, 10);
        h_gen_massVBF_ = new TH1F("h_gen_massVBF", "Invariant mass of VBF quark system", 1080, 0, 3500);
        // bQ kinematics
        h_gen_massbQsystem_ = new TH1F("h_gen_massbQsystem", "Mass of bottom quark system", 1080, 0, 200);
        h_gen_ptb0_ = new TH1F("h_gen_ptb0", "pt of leading bottom quarks", 1080, 0, 400);
        h_gen_ptb1_ = new TH1F("h_gen_ptb1", "pt of subleading bottom quarks", 1080, 0, 400);
        // Leptons kinematics
        ptLep1 = new TH1F("ptLep1", "pt of leading leptons_", 1080, 0, 600);
        ptLep2 = new TH1F("ptLep2", "pt of subleading leptons_", 1080, 0, 600);
        etaLep1 = new TH1F("etaLep1", "eta of leading leptons_", 1080, -5, 5);
        etaLep2 = new TH1F("etaLep2", "eta of subleading leptons_", 1080, -5, 5);
        phiLep1 = new TH1F("phiLep1", "phi of leading leptons_", 1080, -5, 5);
        phiLep2 = new TH1F("phiLep2", "phi of subleading leptons_", 1080, -5, 5);
        massDiLep = new TH1F("massDiLep", "mass of dileptons", 1080, 0, 250);
        ptDiLep = new TH1F("ptDiLep", "pt of the dilepton system", 1080, 0, 1800);
        dRLep = new TH1F("dRLep", "delta R between two leptons_", 1080, 0, 5);

        // jets kinematics
        etaVBFjet1 = new TH1F("etaVBFjet1", "eta of leading VBF jet", 1080, -5, 5); 
        etaVBFjet2 = new TH1F("etaVBFjet2", "eta of subleading VBF jet", 1080, -5, 5); 
        re_deltaEtaVBF = new TH1F("re_deltaEtaVBF", "Delta Eta of VBF quarks (Reconstructed)", 1080, 0, 10);
        ptVBFjet1 = new TH1F("ptVBFjet1", "pt of leading VBF jets", 1080, 0, 500);
        ptVBFjet2 = new TH1F("ptVBFjet2", "pt of subleading VBF jets", 1080, 0, 500);
        VBFjetMass = new TH1F("VBFjetMass", "Invariant mass of VBF quark system", 1080, 0, 3500);
        ptB1 = new TH1F("ptB1", "pt of leading bottom quarks", 1080, 0, 400);
        ptB2 = new TH1F("ptB2", "pt of subleading bottom quarks", 1080, 0, 400);
        BjetMass = new TH1F("BjetMass", "Invariant mass of bottom quark system", 1080, 0, 300);
        BjetMass_400cut = new TH1F("BjetMass_400cut", "Invariant mass of b quark system after the cut", 1080, 0, 300);
        ptWjet1 = new TH1F("ptWjet1", "pt of leading jet quarks from W bosons", 1080, 0, 300);
        ptWjet2 = new TH1F("ptWjet2", "pt of subleading jet quarks from W bosons", 1080, 0, 300);
        WjetMass = new TH1F("WjetMass", "Invariant mass of jet quark system", 1080, 0, 100); 

        //______________________________________________________
        // fatJets Kinematics
        ptHbb = new TH1F("ptHbb", "pt of Hbb fatjet", 1080, 0, 1800);
        etaHbb = new TH1F("etaHbb", "eta of Hbb fatjet", 1080, -5, 5);
        massHbb = new TH1F("massHbb", "mass of Hbb fatjet", 1080, 0, 200);
        nFatjets = new TH1F("nFatjets", "Number of fatjets", 5, 0, 5);

        //_______________________________________________________
        // s-hat variable
        massZH = new TH1F("massZH", "Invariant mass of the ZH system", 1080, 0, 3500);
        massZHzoom = new TH1F("massZHzoom", "Invariant mass of the ZH system", 1080, 0, 500);
    }

    //_______________________________________________________
    // Fill fatjet Histograms
    void fillHbbFatJetKinematicHistograms()
    {
        ptHbb->Fill(Analysis::hbbFatJet_.Pt(), Analysis::wgt_);
        etaHbb->Fill(Analysis::hbbFatJet_.Eta(), Analysis::wgt_);
        massHbb->Fill(Analysis::hbbFatJet_.M(), Analysis::wgt_);
        nFatjets->Fill(Analysis::fatJets_.size(), Analysis::wgt_);
    }



    //_______________________________________________________
    // Fill Lepton histograms
    void fillLeptonsKinematicHistograms()
    {
        LV lep1 = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[0] : Analysis::leptons_[1];
        LV lep2 = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[1] : Analysis::leptons_[0];
        ptLep1->Fill(lep1.Pt(), Analysis::wgt_);
        ptLep2->Fill(lep2.Pt(), Analysis::wgt_);
        etaLep1->Fill(lep1.Eta(), Analysis::wgt_);
        etaLep2->Fill(lep2.Eta(), Analysis::wgt_);
        phiLep1->Fill(lep1.Phi(), Analysis::wgt_);
        phiLep2->Fill(lep2.Phi(), Analysis::wgt_);
        massDiLep->Fill((lep1+lep2).M(), Analysis::wgt_);
        ptDiLep->Fill((lep1+lep2).Pt(), Analysis::wgt_);
        dRLep->Fill(RooUtil::Calc::DeltaR(lep1, lep2), Analysis::wgt_);
    }



    void fillJetsKinematicHistograms()
    {
        LV wj1 = (Analysis::Wjets_[0]).Pt() > (Analysis::Wjets_[1]).Pt() ? Analysis::Wjets_[0] : Analysis::Wjets_[1];
        LV wj2 = (Analysis::Wjets_[0]).Pt() > (Analysis::Wjets_[1]).Pt() ? Analysis::Wjets_[1] : Analysis::Wjets_[0];
        LV VBF1 = (Analysis::VBFjets_[0]).Pt() > (Analysis::VBFjets_[1]).Pt() ? Analysis::VBFjets_[0] : Analysis::VBFjets_[1];
        LV VBF2 = (Analysis::VBFjets_[0]).Pt() > (Analysis::VBFjets_[1]).Pt() ? Analysis::VBFjets_[1] : Analysis::VBFjets_[0];
        ptWjet1->Fill(wj1.Pt(), Analysis::wgt_);
        ptWjet2->Fill(wj2.Pt(), Analysis::wgt_);
        WjetMass->Fill((wj1+wj2).M(), Analysis::wgt_);

        etaVBFjet1->Fill(VBF1.Eta(), Analysis::wgt_);
        etaVBFjet2->Fill(VBF2.Eta(), Analysis::wgt_);

        ptVBFjet1->Fill(VBF1.Pt(), Analysis::wgt_);
        ptVBFjet2->Fill(VBF2.Pt(), Analysis::wgt_);

        re_deltaEtaVBF->Fill(TMath::Abs(Analysis::VBFjets_[0].Eta()-Analysis::VBFjets_[1].Eta()), Analysis::wgt_);
        VBFjetMass->Fill((Analysis::VBFjets_[0] + Analysis::VBFjets_[1]).M());

        // Figuring out which combination of VBF jets is good
        float dR11 = RooUtil::Calc::DeltaR(VBF1, Analysis::gen_jet0_);
        float dR22 = RooUtil::Calc::DeltaR(VBF2, Analysis::gen_jet1_);
        float dR12 = RooUtil::Calc::DeltaR(VBF2, Analysis::gen_jet0_);
        float dR21 = RooUtil::Calc::DeltaR(VBF1, Analysis::gen_jet1_);

        LV VBF1_matched = dR11 + dR22 < dR12 + dR21 ? VBF1 : VBF2;
        LV VBF2_matched = dR11 + dR22 < dR12 + dR21 ? VBF2 : VBF1;

        // dR between reco- and gen- level
        dRleadingVBF->Fill(RooUtil::Calc::DeltaR(VBF1_matched, Analysis::gen_jet0_), Analysis::wgt_);
        dRsubleadingVBF->Fill(RooUtil::Calc::DeltaR(VBF2_matched, Analysis::gen_jet1_), Analysis::wgt_);
        dRhiggs->Fill(RooUtil::Calc::DeltaR(Analysis::hbbFatJet_, Analysis::gen_H_), Analysis::wgt_);

    }
    

    // Fill gen-level histograms
    void fillGenLevelHistograms()
    {
        h_gen_massVBF_->Fill((Analysis::gen_ijet_ + Analysis::gen_jjet_).M(), Analysis::wgt_);
        h_gen_deltaEta_->Fill(TMath::Abs(Analysis::gen_ijet_.Eta() - Analysis::gen_jjet_.Eta()), Analysis::wgt_);
        for (unsigned int j = 0; j < gen_Z_.size(); j++) {h_gen_ptZ_->Fill(Analysis::gen_Z_[j].Pt(), Analysis::wgt_);}
        for (unsigned int i = 0; i < gen_W_.size(); i++) { h_gen_ptW_->Fill(Analysis::gen_W_[i].Pt(), Analysis::wgt_);}
        h_gen_ptH_->Fill(Analysis::gen_H_.Pt(), Analysis::wgt_);
        h_gen_massbQsystem_->Fill((Analysis::gen_b0_+Analysis::gen_b1_).M(), Analysis::wgt_); 
        h_gen_ptb0_->Fill(Analysis::gen_b0_.Pt(), Analysis::wgt_);
        h_gen_ptb1_->Fill(Analysis::gen_b1_.Pt(), Analysis::wgt_);
    }

    // Fill s-hat kinematic variables
    void fillSHatHistograms()
    {
        massZH->Fill((Analysis::leptons_[0] + Analysis::leptons_[1] + Analysis::hbbFatJet_).M(), Analysis::wgt_);
        massZHzoom->Fill((Analysis::leptons_[0] + Analysis::leptons_[1] + Analysis::hbbFatJet_).M(), Analysis::wgt_);
    }
    

 
    //_______________________________________________________
    // Write the histograms
    void writeHistograms(TFile* ofile)
    {
        // Done with the analysis
        ofile->cd();
        // Gen-level
        h_gen_ptb0_->Write();
        h_gen_ptb1_->Write();
        h_gen_massbQsystem_->Write();
        h_gen_massVBF_->Write();
        h_gen_deltaEta_->Write();
        h_gen_ptH_->Write();
        h_gen_ptZ_->Write();
        h_gen_ptW_->Write();
        ptLep1->Write();
        ptLep2->Write();
        re_deltaEtaVBF->Write();
        etaLep1->Write();
        etaLep2->Write();
        phiLep1->Write();
        phiLep2->Write();
        massDiLep->Write();
        ptDiLep->Write();
        dRLep->Write();
        etaVBFjet1->Write();
        etaVBFjet2->Write();
        ptVBFjet2->Write();
        ptVBFjet1->Write();
        VBFjetMass->Write();
        ptB1->Write();
        ptB2->Write();
        BjetMass->Write();
        ptWjet1->Write();
        ptWjet2->Write();
        WjetMass->Write();
        dRleadingVBF->Write();
        dRsubleadingVBF->Write();
        dRsubleadingbq->Write();
        dRleadingbq->Write();
        recombdRleadingVBF->Write();
        recombdRleadingbq->Write();
        ptHbb->Write();
        etaHbb->Write();
        massHbb->Write();
        nFatjets->Write();
        massZH->Write();
        massZHzoom->Write();
        dRhiggs->Write();
    }

}




        
        
//====================================================================
// Main function
//====================================================================
int main() 
{
    // Create your output root file
    TFile* output_file = new TFile("VBSOSWWHoutput.root", "recreate");

    // Create Histograms
    Hist::bookHistograms();

    // Create Cutflow Histogram
    Cutflow::bookCutflow();

    // Set scale 1fb (the per event weight normalized for 1/fb)
    Analysis::setScale1fb(5.652/645745); // 11.23 fb / (# of events = 642498)

    // Set the luminosity
    Analysis::setLumi(137); // TODO: Update properly in the future. For now it's a Placeholder!

    // Set the year 
    Analysis::setYear(2018); // TODO: Update properly in the future. For now it's a Placeholder!

    // Set up year dependent configuration of the analysis that are POG specific from CMS
    Analysis::setConfig();

    // Input Path (RooUtil::Looper can accept comma separated list)
    std::vector<TString> input_file_paths;
    for (unsigned i = 0; i < 26; i++)
    {
        TString input_temp = "/home/users/joytzphysics/VBSOSWWH_C2V_3/output_" + std::to_string(i) + ".root";
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

    // Tracking cutflow number;
    int cut1_events = 0;
    int cut2_events = 0;
    int cut3_events = 0;
    int cut4_events = 0;
    int cut5_events = 0;
    int cut6_events = 0;

    // Loop through events
    while (looper.nextEvent())
    {
        // print out event details
//      dumpGenParticleInfos();

        // print out information
          
        // Run the analysis algorithms (selections, computing variables, etc.)
        Analysis::runAlgorithms();
        // Fill the "counter" histogram
        Cutflow::fillCutflow(Cutflow::Cuts::kNoSelection);

        // Cut#1: Check if vbs
        if (Analysis::is_vbs == false) { continue;}
        cut1_events ++;
        Cutflow::fillCutflow(Cutflow::Cuts::kVbs);
        // Cut#2: Check if hbb
        if (Analysis::is_hbb == false) { continue;}
        cut2_events ++;
        Cutflow::fillCutflow(Cutflow::Cuts::kHbb);
        // Cut#3: Check if wjj
//        if (Analysis::is_wjj == false) { continue;}
//        cut3_events ++;
//        Cutflow::fillCutflow(Cutflow::Cuts::kWjj);

        // Cut#3: Require that there are exactly two leptons
        if (not (Analysis::elecs_.size() + Analysis::muons_.size() == 2)) { continue; }
        // Cut#3: Require that the two leptons have OS (opposite-sign)
        int is_os = false;
        if (Analysis::elecs_.size() == 2)
        {
            is_os = Analysis::elecs_[0].pdgid * Analysis::elecs_[1].pdgid < 0;
        }
        else if (Analysis::muons_.size() == 2)
        {
            is_os = Analysis::muons_[0].pdgid * Analysis::muons_[1].pdgid < 0;
        }
        else if (Analysis::muons_.size() == 1 && Analysis::elecs_.size() == 1)
	{
	    is_os = Analysis::muons_[0].pdgid * Analysis::elecs_[0].pdgid < 0;
	}

        if (not (is_os)) { continue; }
	// Cut#3: Require that the leading leptons Pt >= 40, subleading >= 30
        LV leading = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[0] : Analysis::leptons_[1];
        LV subleading = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[1] : Analysis::leptons_[0];
        
	if (not (leading.Pt() >= 40 && subleading.Pt() >= 30)) {continue;}
        Cutflow::fillCutflow(Cutflow::Cuts::kTwoLightLeptons);
        cut3_events ++; 

        // Cut#4: Require at least one fatjet with softdropmass > 40 GeV
        if (not (Analysis::fatJets_.size() >= 1 ) ) { continue;}
        cut4_events ++;
        Cutflow::fillCutflow(Cutflow::Cuts::kOneHbbFatJet);
	// Cut#5: Require the Hbb score > 0.5
	int no_Hbb = 0;
	for (unsigned int ifatjet = 1; ifatjet < fatJets_.size(); ifatjet++) {
            if (Analysis::fatJets_[ifatjet].hbbScore >= 0.5) { no_Hbb += 1; }
	}
	if (no_Hbb = 0) { continue; }
	cut5_events ++;
	Cutflow::fillCutflow(Cutflow::Cuts::kHbbScore);
        // Cut#6: Require that there are at least 2 pt > 30 GeV jets
        
        if (not (Analysis::jets_.size() >= 2)) { continue; }
        cut6_events ++;
        Cutflow::fillCutflow(Cutflow::Cuts::kAtLeastTwoPt30Jets);

	// Cut#7: Require mjj > 500, deltaEtajj > 3
	if (not (Analysis::VBFjets_[0] + Analysis::VBFjets_[1]).M() > 500) { continue;}
	if (not (TMath::Abs(Analysis::VBFjets_[0].Eta()-Analysis::VBFjets_[1].Eta()) > 3)) { continue;}
	Cutflow::fillCutflow(Cutflow::Cuts::kMjj);

        // Fill Lepton Histogram;
        Hist::fillLeptonsKinematicHistograms();

        // All cuts have passed
        // Now fill the histograms
        Hist::fillGenLevelHistograms();
        Hist::fillSHatHistograms();
        Hist::fillJetsKinematicHistograms();
        Hist::fillHbbFatJetKinematicHistograms();
 
    }

    // Write out the histograms
    Hist::writeHistograms(output_file);

    // Write out the cutflow histogram
    Cutflow::writeCutflow(output_file);

    // Produce cutflow events number
    cout << "The first cut (vbswzh) produces " << cut1_events << " events." << endl;
    cout << "The second cut (hbb) produces " << cut2_events << " events." << endl;
    cout << "The third cut (wjj) produces " << cut3_events << " events." << endl;
    cout << "The fourth cut (two leptons) produces " << cut4_events << "events." << endl;
    cout << "The fifth cut (1 fatjet) produces " << cut5_events << " events." << endl;
    cout << "The sixth cut (6 jets) produces " << cut6_events << " events." << endl;

    // Close the file
    output_file->Close();


    return 0;
}
