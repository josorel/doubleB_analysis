#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"

void test_TMVA()
{
    // Load data
    TString label_model = "bjet";
    // TString label_data = "herwig_dijet_deepJet";
    // TString label_data = "herwig_dijet_pnet";
    TString label_data = "bjet";
    // TString label = "ttbar_highPU";
    // TString fname = "./saved_models/TMVA.root";
    // TString fname = "./data_root_ttbar_highPU_30_pt_700/data.root";
    TString fname = "/data_CMS/cms/sorel/SV_bjet_TMVA.root";
    std::cout << "Reading data from: " << fname << std::endl;
    TFile *fin = new TFile(fname);
    TTree *testTree = (TTree *) fin->Get(label_model + "_dataloader/TestTree");

    // Load the TMVA model
    TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

    Float_t SVdist_bdt;
    Float_t SVmass1_bdt;
    Float_t SVpt1_bdt;
    Float_t SVntrks1_bdt;
    Float_t SVmass2_bdt;
    Float_t SVpt2_bdt;
    Float_t SVntrks2_bdt;

    Int_t classID_bdt;

    Float_t SVdist_tree;
    Float_t SVmass1_tree;
    Float_t SVpt1_tree;
    Float_t SVntrks1_tree;
    Float_t SVmass2_tree;
    Float_t SVpt2_tree;
    Float_t SVntrks2_tree;

    Int_t classID_tree;
    Float_t weight_tree;

    std::vector<TString> variable_names_ = {"SVdist", "SVmass1", "SVpt1", "SVntrks1",
                                            "SVmass2", "SVpt2", "SVntrks2",
                                            };

    std::map<TString, float *> variables_bdt;
    variables_bdt["SVdist"] = &SVdist_bdt;
    variables_bdt["SVmass1"] = &SVmass1_bdt;
    variables_bdt["SVpt1"] = &SVpt1_bdt;
    variables_bdt["SVntrks1"] = &SVntrks1_bdt;
    variables_bdt["SVmass2"] = &SVmass2_bdt;
    variables_bdt["SVpt2"] = &SVpt2_bdt;
    variables_bdt["SVntrks2"] = &SVntrks2_bdt;

    std::map<TString, float *> variables_tree;
    variables_tree["SVdist"] = &SVdist_tree;
    variables_tree["SVmass1"] = &SVmass1_tree;
    variables_tree["SVpt1"] = &SVpt1_tree;
    variables_tree["SVntrks1"] = &SVntrks1_tree;
    variables_tree["SVmass2"] = &SVmass2_tree;
    variables_tree["SVpt2"] = &SVpt2_tree;
    variables_tree["SVntrks2"] = &SVntrks2_tree;

    for (TString var : variable_names_) {
        // std::cout << "var: " << var << std::endl;
        reader->AddVariable(var, variables_bdt[var]);
        testTree->SetBranchAddress(var, variables_tree[var]);
    }
    testTree->SetBranchAddress("classID", &classID_tree);
    testTree->SetBranchAddress("weight", &weight_tree);

    TString weightfile = label_model+"_dataloader/weights/TMVAClassification_BDTG.weights.xml";
    // TString weightfile = "/home/llr/cms/kalipoliti/C10630p1_miniAOD/src/RecoHI/HiJetAlgos/data/TMVAClassification_BDTG.weights.xml";
    // std::cout << "Reading weightfile from: " << weightfile << std::endl;
    reader->BookMVA("BDTG", weightfile);

    // Declare WP
    float wp = -0.6; // -0.3
    // float wp = 0.;

    // Get WP effS and effB
    float totalS = 0;
    float totalB = 0;
    float passS = 0;
    float passB = 0;

    // Fill histo for roc curve
    TH1F *h_discr_s = new TH1F("h_discr_s", "x=discr value, y=events with that value", 40, -1., 1.);
    TH1F *h_discr_b = new TH1F("h_discr_b", "x=discr value, y=events with that value", 40, -1., 1.);

    for (Long64_t ient = 0; ient < testTree->GetEntries(); ient++) {
        testTree->GetEntry(ient);

        weight_tree = 1.;

        // [DEBUG]
        // if (ient > 0) break;

        SVdist_bdt=SVdist_tree;
        SVmass1_bdt=SVmass1_tree;
        SVpt1_bdt=SVpt1_tree;
        SVntrks1_bdt=SVntrks1_tree;
        SVmass2_bdt=SVmass2_tree;
        SVpt2_bdt=SVpt2_tree;
        SVntrks2_bdt=SVntrks2_tree;

        classID_bdt = classID_tree;
        
        bool pass = false;
        float proba = reader->EvaluateMVA("BDTG");
        // std::cout << proba << std::endl;

        if (proba > wp) {
            pass = true;
        }

        // totalS += weight_tree;
        // if (pass) passS += weight_tree;
        // h_discr_s->Fill(proba, weight_tree);
        // totalB += weight_tree;
        // if (pass) passB += weight_tree;
        // h_discr_b->Fill(proba, weight_tree);

        // cout << "ID: " << classID_bdt << endl;

        if (classID_bdt == 0) {
            totalS += weight_tree;
            if (pass) passS += weight_tree;
            h_discr_s->Fill(proba, weight_tree);
        } else {
            totalB += weight_tree;
            if (pass) passB += weight_tree;
            h_discr_b->Fill(proba, weight_tree);
        }
    } // end signal tree loop

    // Save the histos
    TString foutName = label_data + "_roc.root";
    TFile *fout = new TFile(foutName, "recreate");
    
    for (auto h : {h_discr_s, h_discr_b}) {
        h->Write();
    }
    fout->Close();
    delete fout;

    std::cout << "for wp : " << wp 
              << "\n\tsignal efficiency = " << (float) passS / totalS
              << "\n\tbkg rejection = " << 1. - ((float) passB / totalB)
              << std::endl;
}