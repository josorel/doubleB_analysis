#include "tTree.h"
// #include "binning.h"

bool skipMC(Float_t jtpt, Float_t refpt, Float_t pthat) {
    if (!(refpt>0)) return true;    
    if (pthat<0.35*jtpt) return true;
    return false;
}

// Int_t findMostCommon(const std::vector<Int_t>& nums) {
//     std::unordered_map<Int_t, Int_t> frequencyMap;

//     // Count the frequency of each number
//     for (Int_t num : nums) {
//         frequencyMap[num]++;
//     }

//     // Variables to track the most common number and its count
//     Int_t mostCommon = nums[0];
//     Int_t maxCount = 0;
//     Int_t tieCount = 0;

//     // Find the number with the highest frequency
//     for (const auto& entry : frequencyMap) {
//         if (entry.second > maxCount) {
//             mostCommon = entry.first;
//             maxCount = entry.second;
//             tieCount = 0;  // Reset tie counter when a new highest frequency is found
//         } else if (entry.second == maxCount) {
//             tieCount++;  // Increment tie counter if there is a tie
//         }
//     }

//     // If there is more than one number with the highest frequency, return -999
//     if (tieCount > 0) {
//         return -999;
//     }

//     return mostCommon;
// }

// DISPLAY METHODS --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void displayMerges(const std::map<Float_t, std::vector<ROOT::Math::PtEtaPhiMVector>>& Merges) {
    std::cout << "Displaying Merges:" << std::endl;
    for (const auto& entry : Merges) {
        // Print the Pt() value (the key)
        std::cout << "Key (Pt): " << entry.first << std::endl;
        
        // Print the contents of the associated vector (the value)
        std::cout << "Tracks: " << std::endl;
        for (const auto& track : entry.second) {
            std::cout << "  Track Pt: " << track.Pt() << ", Eta: " << track.Eta() 
                      << ", Phi: " << track.Phi() << ", M: " << track.M() << std::endl;
        }
    }
}

void displayTrackVectors(const std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors) {
    // std::cout << "Displaying Track Vectors:" << std::endl;

    Float_t totalPt = 0;
    
    // Iterate over each track in the vector and display its properties
    for (size_t i = 0; i < trackVectors.size(); ++i) {
        const ROOT::Math::PtEtaPhiMVector& track = trackVectors[i];
        std::cout << "Track " << i + 1 << ": "
                  << "Pt: " << track.Pt() 
                  << ", Eta: " << track.Eta() 
                  << ", Phi: " << track.Phi() 
                  << ", M: " << track.M() 
                  << std::endl;
                  totalPt += track.Pt();
    }
}

// AGGREGATION AND DECLUSTERING METHODS

// // This method groups closest secondary vertices until only two remain
// void groupVertexes(std::vector<ROOT::Math::PtEtaPhiMVector>& vecFinalSecVtxs){

//     while (vecFinalSecVtxs.size() > 2){
//         double min_distance=std::numeric_limits<double>::infinity();
//         Int_t index1 = 0, index2 = 0;

//         //loop through given tracks to find those with smallest distance
//         for (Int_t i = 0; i < vecFinalSecVtxs.size(); ++i) {
//             for (Int_t j = i + 1; j < vecFinalSecVtxs.size(); ++j) {
//                 Float_t dist = std::pow((vecFinalSecVtxs[i].Eta()-vecFinalSecVtxs[j].Eta()),2)+std::pow(acos(cos(vecFinalSecVtxs[i].Phi()-vecFinalSecVtxs[j].Phi())),2);
//                if (dist < min_distance) {
//                     min_distance = dist;
//                     index1 = i;
//                     index2 = j;
//                 }
//             }
//         }

//         ROOT::Math::PtEtaPhiMVector newVector = vecFinalSecVtxs[index1]+vecFinalSecVtxs[index2];

//         // Add combined track
//         vecFinalSecVtxs.push_back(newVector);
//         // remove merged tracks
//         if (index1 > index2) std::swap(index1, index2);
//         vecFinalSecVtxs.erase(vecFinalSecVtxs.begin() + index2);
//         vecFinalSecVtxs.erase(vecFinalSecVtxs.begin() + index1);
//     }
// }

// // This method recreates B hadrons using secondary vertices
// void ReconstructBs(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, map<Int_t, ROOT::Math::PtEtaPhiMVector> finalSecVtxs){
//     std::vector<Int_t> indices;

//     // You can loop over the tracks of a specific jet using trkJetId
//     for (Int_t itrk = 0; itrk < trackVectors.size(); itrk++) {

//         //if it is from a b hadron, create vector v1 for track and set properties
//         if (t.trkMatchSta[itrk]>1){
//             ROOT::Math::PtEtaPhiMVector v1;

//             v1.SetEta(t.trkEta[itrk]);
//             v1.SetPt(t.trkPt[itrk]);
//             v1.SetPhi(t.trkPhi[itrk]);

//             if(std::abs(t.trkPdgId[itrk])==211){
//                 v1.SetM(0.139570);
//             }
//             if(std::abs(t.trkPdgId[itrk])==13){
//                 v1.SetM(0.105658);
//             }
//             if(std::abs(t.trkPdgId[itrk])==11){
//                 v1.SetM(0.000510);
//             }
//             if(std::abs(t.trkPdgId[itrk])==2212){
//                 v1.SetM(0.938272);
//             }
//             if(std::abs(t.trkPdgId[itrk])==321){
//                 v1.SetM(0.493677);
//             }

//             // erase from list of tracks
//             for (Int_t i = 0; i < trackVectors.size(); i++){
//                 if (trackVectors[i]==v1){
//                     indices.push_back(i);
//                 }
//             }

//             // add to map of reconstructed b hadrons
//             for (auto i : finalSecVtxs){
//                 if (t.trkEta[itrk]==i.first){
//                     i.second+=v1;
//                 }
//             }
//         }
//     }

//     //add reconstructed b hadrons to list of tracks
//     std::sort(indices.begin(),indices.end(), greater<int>());
//     for (Int_t i = 0; i < indices.size(); i++){
//         trackVectors.erase(trackVectors.begin() + indices[i]);
//     }
//     for (auto i : finalSecVtxs) trackVectors.push_back(i.second);
// }

// SUBJET CLUSTERING METHOD - groups tracks until only two remain (which are then our subjets)
void CambridgeAachen(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, std::map<Float_t, std::vector<ROOT::Math::PtEtaPhiMVector>>& Merges){

    while (trackVectors.size() > 2){
        double min_distance=std::numeric_limits<double>::infinity();
        Int_t index1 = 0, index2 = 0;

        //loop through given tracks to find those with smallest distance
        for (Int_t i = 0; i < trackVectors.size(); ++i) {
            for (Int_t j = i + 1; j < trackVectors.size(); ++j) {
                Float_t dist = std::pow((trackVectors[i].Rapidity()-trackVectors[j].Rapidity()),2)+std::pow(acos(cos(trackVectors[i].Phi()-trackVectors[j].Phi())),2);
                if (dist < min_distance) {
                    min_distance = dist;
                    index1 = i;
                    index2 = j;
                }
            }
        }

        //keep track of merges - needed for Soft Drop
        ROOT::Math::PtEtaPhiMVector newVector = trackVectors[index1]+trackVectors[index2];
        Merges.insert({newVector.Pt(), {trackVectors[index1],trackVectors[index2]}});

        for (Int_t k = 0;k<(Merges[newVector.Pt()]).size();k++){
            for (auto it = Merges.begin(); it != Merges.end(); ++it){
                if (abs((Merges[newVector.Pt()])[k].Pt()-(it->first))<pow(10,-2)){
                    std::vector<ROOT::Math::PtEtaPhiMVector> tempTracks = it->second;
                    Merges[newVector.Pt()].erase(Merges[newVector.Pt()].begin() + k);
                    Merges[newVector.Pt()].insert(Merges[newVector.Pt()].end(),tempTracks.begin(), tempTracks.end());
                }
            }
        }


        // Add combined track
        trackVectors.push_back(newVector);
        // remove merged tracks
        if (index1 > index2) std::swap(index1, index2);
        trackVectors.erase(trackVectors.begin() + index2);
        trackVectors.erase(trackVectors.begin() + index1);

    }
}

// SOFTDROP METHOD
Int_t SoftDrop(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors){
    if (trackVectors.size() < 2) return 2;

    // check if one of the subjets is only 10% of total Pt
    // SoftDrop is 1 if it fails and 0 if it passes
    Int_t softDrop = 1;
    if(trackVectors[0].Pt()<trackVectors[1].Pt()){
        if (trackVectors[0].Pt()/(trackVectors[0].Pt()+trackVectors[1].Pt())>0.1){
            softDrop = 0;
        }
    }
    else{
        if (trackVectors[1].Pt()/(trackVectors[0].Pt()+trackVectors[1].Pt())>0.1){
            softDrop = 0;
        }
    }
    return softDrop;
}

// MAIN METHOD --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void run_data_inclusive()
{   
    // Load tree with class 

    // TFile *_file0 = TFile::Open("/data_CMS/cms/kalipoliti/qcdMC/dijet/aggrTMVA_fixedGenBug/merged_HiForestMiniAOD.root"); 
    // TFile *_file0 = TFile::Open("/data_CMS/cms/kalipoliti/qcdMC/bjet/noAggr_withPNET/merged_HiForestMiniAOD.root"); 
    TFile *_file0 = TFile::Open("/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_withPFJetTrigger/merged_HiForestMiniAOD.root"); 
    TString fin_name = "/data_CMS/cms/kalipoliti/bJet2017G/LowEGJet/aggrTMVA_withPFJetTrigger/merged_HiForestMiniAOD.root";
    TString fout_name = "histogramsIdealSvtxDataLowNewBins.root";
    tTree t(fin_name);

    Double_t weight=1;

    // Declare histograms

    const Int_t jtpt_binsVectorSize = 4;
    Int_t jtpt_bins = jtpt_binsVectorSize - 1;
    Double_t jtpt_binsVector[jtpt_binsVectorSize] = {
        60., 
        80., 
        110.,
        150.
    };
    const Int_t logrg_binsVectorSize = 7;
    Int_t logrg_bins = logrg_binsVectorSize - 1;
    Double_t logrg_binsVector[logrg_binsVectorSize] = {
        // -1.2, 
        // -1.0,
        0., 
        0.25,
        0.5,
        0.75,
        1.0,
        1.3,
        3.0
    };
    const Int_t zg_binsVectorSize = 5;
    Int_t zg_bins = zg_binsVectorSize - 1;
    Double_t zg_binsVector[zg_binsVectorSize] = {
        // -0.1, 
        // -0.2,
        0.1, 
        // 0.15,
        // 0.15,
        0.2,
        // 0.25,
        // 0.25,
        0.3,
        // 0.35,
        0.4,
        // 0.45,
        0.5,
    };
    const Int_t msv_binsVectorSize = 8;
    Int_t msv_bins = msv_binsVectorSize - 1;
    Double_t msv_binsVector[msv_binsVectorSize] = {
        // -1.2, 
        0., 
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
        6.0,
        7.0,
    };
    const Int_t kt_binsVectorSize = 5;
    Int_t kt_bins = kt_binsVectorSize - 1;
    Double_t kt_binsVector[kt_binsVectorSize] = {
        // -1.0,
        0., 
        0.4,
        0.8,
        1.2,
        3,
    };

    TH3D *h_Zg_inclusive = new TH3D("h_InvMass_vs_Zg_Data", "svtx invariant mass vs zg Data", msv_bins, msv_binsVector ,zg_bins, zg_binsVector,jtpt_bins, jtpt_binsVector);
    TH3D *h_Kt_inclusive = new TH3D("h_InvMass_vs_Kt_Data", "svtx invariant mass vs kt Data", msv_bins, msv_binsVector ,kt_bins, kt_binsVector,jtpt_bins, jtpt_binsVector);
    TH3D *h_Rg_inclusive = new TH3D("h_InvMass_vs_Dr_Data", "svtx invariant mass vs dr Data", msv_bins, msv_binsVector ,logrg_bins, logrg_binsVector,jtpt_bins, jtpt_binsVector);

    // Turn off all branches and turn on only the interesting ones
    // Attention! If a branch is off, it will return bs without crashing 
    t.SetBranchStatus("*", 0);
    std::vector<TString> active_branches = {
        "weight",
        "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1",
        "nref", "jtpt", "jteta", "lumi", "jtptCh", "discr_particleNet_BvsAll", "pthat", // reco jet
        "refpt", "refeta", "jtHadFlav", "jtNbHad", //"sjt1Eta", "sjt2Eta", "sjt1Pt", "sjt2Pt",  // gen jet
        "ntrk", "trkPt", "trkMatchSta", "trkJetId", "trkSvtxId", "trkPdgId", "trkEta", "trkPhi", // tracks
        "jtNsvtx", "svtxpt", "svtxJetId", //"svtxm","svtxmcorr","svtxdls", "svtxdls2d", "svtxNtrk", "svtxnormchi2", "svtxpt", // SVs 
    };
    t.SetBranchStatus(active_branches, 1);
    
    std::cout << "Looping over events" << std::endl; 

    cout << "Events: " << t.GetEntries() << endl;
    for (Long64_t ient = 0; ient < t.GetEntries(); ient++) {

        // Print progress
        if (ient % 100000 == 0) {
            std::cout << "entry nb = " << ient << std::endl;
        }
 
        t.GetEntry(ient); 

        // Loop over jets  ===============================================================================================================
        for (Int_t ijet = 0; ijet < t.nref; ijet++) {

            //REQUIREMENTS**********************************************************************
            // Pt and eta cuts
            if (std::abs(t.jteta[ijet]) > 2.) continue;
            if ((t.jtpt[ijet]) < 60 || (t.jtpt[ijet]) > 150) continue; 

            // DATA REQS - use triggers 80+100 for high and 40+60 for low
            // if (t.HLT_HIAK4PFJet80_v1==0&&t.HLT_HIAK4PFJet100_v1==0) continue;
            if (t.HLT_HIAK4PFJet40_v1==0&&t.HLT_HIAK4PFJet60_v1==0) continue;
            if (t.HLT_HIAK4PFJet60_v1==1&&(t.HLT_HIAK4PFJet80_v1==1||t.HLT_HIAK4PFJet100_v1==1)) continue;
            if (t.HLT_HIAK4PFJet40_v1==1&&(t.HLT_HIAK4PFJet60_v1==1||t.HLT_HIAK4PFJet80_v1==1||t.HLT_HIAK4PFJet100_v1==1)) continue;

            if (t.HLT_HIAK4PFJet60_v1==1) weight=1;
            if (t.HLT_HIAK4PFJet40_v1==1) weight=33.917210;

            //***********************************************************************************

            //vector of TLorentz vectors for each track inside jet, used for subjet clustering
            std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;

            // You can loop over the tracks of a specific jet using trkJetId --------------------------------------------------------------
            for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
                if (t.trkJetId[itrk] != ijet) continue; 

                //create vectors for each track
                ROOT::Math::PtEtaPhiMVector v1;
                v1.SetEta(t.trkEta[itrk]);
                v1.SetPt(t.trkPt[itrk]);
                v1.SetPhi(t.trkPhi[itrk]);

                if(std::abs(t.trkPdgId[itrk])==211){
                    v1.SetM(0.139570);
                }
                if(std::abs(t.trkPdgId[itrk])==13){
                    v1.SetM(0.105658);
                }
                if(std::abs(t.trkPdgId[itrk])==11){
                    v1.SetM(0.000510);
                }
                if(std::abs(t.trkPdgId[itrk])==2212){
                    v1.SetM(0.938272);
                }
                if(std::abs(t.trkPdgId[itrk])==321){
                    v1.SetM(0.493677);
                }

                trackVectors.push_back(v1);
            }

            std::vector<ROOT::Math::PtEtaPhiMVector> initialVectors = trackVectors;

            //do not cluster into subjets if there are less than two tracks
            if (trackVectors.size()<2) continue;

            // DECLUSTERING: create two subjets from tracks in jet ------------------------------------------------------------------------------
            std::map<Float_t, std::vector<ROOT::Math::PtEtaPhiMVector>> Merges;
            CambridgeAachen(trackVectors, Merges);

            //SOFTDROP: make sure new tracks pass softdrop, if not go back to previous merges, redo cambridgeAachen, and test again ----------------------------
            Int_t failed = 0;
            while (SoftDrop(trackVectors)==1&&failed!=1){
                if(trackVectors[0].Pt()>trackVectors[1].Pt()){
                    for (auto it = Merges.begin(); it != Merges.end(); ++it) {
                        if (abs((it->first)-trackVectors[0].Pt())<pow(10,-1)){
                            trackVectors = it->second;
                        }
                    }
                }
                else{
                    for (auto it = Merges.begin(); it != Merges.end(); ++it) {
                        if (abs((it->first)-trackVectors[1].Pt())<pow(10,-1)){
                            trackVectors = it->second;
                        }
                    }
                    
                }
                if (trackVectors.size() == 2 && SoftDrop(trackVectors)==1) failed=1;
                if (SoftDrop(trackVectors)==2) failed=1;
                // If Soft Drop fails and only two tracks remain, the jet is discarded
                if (failed == 1) continue;
                Merges.clear();
                initialVectors = trackVectors;
                CambridgeAachen(trackVectors, Merges);
            }

            //do not add jet to histogram if softdrop failed in all cases ----------------------------------------
            if (failed == 1) continue;

            Double_t minPt = std::numeric_limits<double>::infinity();
            for (Int_t i=0; i < 2; i++){
                if (trackVectors[i].Pt() < minPt){
                    minPt = trackVectors[i].Pt();
                }
            }
            // KT CUT
            if (std::log(minPt*(ROOT::Math::VectorUtil::DeltaR(trackVectors[0], trackVectors[1])))<0) continue;

            // Record values for substructure variables
            Float_t Kt_value=0;
            Float_t Zg_value=0;
            Float_t Rg_value=0;
            Rg_value = std::log(0.4/(ROOT::Math::VectorUtil::DeltaR(trackVectors[0], trackVectors[1])));
            if(trackVectors[0].Pt()<trackVectors[1].Pt()){
                Kt_value=std::log(trackVectors[0].Pt()*(ROOT::Math::VectorUtil::DeltaR(trackVectors[0], trackVectors[1])));
                Zg_value=trackVectors[0].Pt()/(trackVectors[0].Pt()+trackVectors[1].Pt());
            }
            else{
                Kt_value=std::log(trackVectors[1].Pt()*(ROOT::Math::VectorUtil::DeltaR(trackVectors[0], trackVectors[1])));
                Zg_value=trackVectors[1].Pt()/(trackVectors[0].Pt()+trackVectors[1].Pt());
            }

            h_Zg_inclusive->Fill(Zg_value, t.jtpt[ijet], weight);
            h_Kt_inclusive->Fill(Kt_value,t.jtpt[ijet], weight);
            h_Rg_inclusive->Fill(Rg_value,t.jtpt[ijet], weight);

        }
    } // entry loop ---------------------------------------------------------------------------------------------------------------------------

    // Save histograms in a new file 
    std::cout << "\n(Re)creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    for (auto h : {h_Zg_inclusive,h_Kt_inclusive,h_Rg_inclusive,}){
        h->Write();
    }

    fout->Close();
}
