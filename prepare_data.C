#include "TFile.h"
#include "tTree.h"
#include "TString.h"

// #include "../../myMacros/C106X/tTree.h"

#include <iostream>

void groupVertexes(std::vector<ROOT::Math::PtEtaPhiMVector>& vecFinalSecVtxs, tTree& t, map<Int_t, ROOT::Math::PtEtaPhiMVector>& finalSecVtxs, map<Float_t, Int_t>& ntrkVector){
    // if (trackVectors.size() < 2) return;

    // in the end, we want 2 vectors corresponding to 2 subjets
    // displayTrackVectors(trackVectors);
    while (vecFinalSecVtxs.size() > 2){
        double min_distance=std::numeric_limits<double>::infinity();
        Int_t index1 = 0, index2 = 0;

        //loop through given tracks to find those with smallest distance
        for (Int_t i = 0; i < vecFinalSecVtxs.size(); ++i) {
            // cout << "Rapidity: " << trackVectors[i].Rapidity() << " eta: " << trackVectors[i].Eta() << endl;
            for (Int_t j = i + 1; j < vecFinalSecVtxs.size(); ++j) {
                Float_t dist = std::pow((vecFinalSecVtxs[i].Eta()-vecFinalSecVtxs[j].Eta()),2)+std::pow(acos(cos(vecFinalSecVtxs[i].Phi()-vecFinalSecVtxs[j].Phi())),2);
                // cout << dist << endl;
                if (dist < min_distance) {
                    min_distance = dist;
                    index1 = i;
                    index2 = j;
                }
            }
        }

        ROOT::Math::PtEtaPhiMVector newVector = vecFinalSecVtxs[index1]+vecFinalSecVtxs[index2];

        // Add combined track
        vecFinalSecVtxs.push_back(newVector);
        ntrkVector.insert({newVector.Pt(),ntrkVector[vecFinalSecVtxs[index1].Pt()]+ntrkVector[vecFinalSecVtxs[index2].Pt()]});
        if (vecFinalSecVtxs.size()!=ntrkVector.size()) break;
        std::map<Float_t,Int_t>::iterator it;
        it=ntrkVector.find(vecFinalSecVtxs[index1].Pt());
        ntrkVector.erase(it);
        it=ntrkVector.find(vecFinalSecVtxs[index2].Pt());
        ntrkVector.erase(it);
        // remove merged tracks
        if (index1 > index2) std::swap(index1, index2);
        vecFinalSecVtxs.erase(vecFinalSecVtxs.begin() + index2);
        vecFinalSecVtxs.erase(vecFinalSecVtxs.begin() + index1);

        // cout << "GROUPED! " << endl;
    }
}

vector<ROOT::Math::PtEtaPhiMVector> makeSvtxs(tTree& t, Int_t& ijet, Float_t& SVntrks1, Float_t& SVntrks2){
    std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;

    map<Int_t, std::vector<ROOT::Math::PtEtaPhiMVector>> secVtxs;
    std::vector<ROOT::Math::PtEtaPhiMVector> emptySVs;
    map<Float_t, Int_t> ntrkVector;

    // map<Int_t, std::vector<Int_t>> secVtxIDs;

    // You can loop over the tracks of a specific jet using trkJetId --------------------------------------------------------------
    for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
        if (t.trkJetId[itrk] != ijet) continue; 

        // cout << "STA: " << t.trkMatchSta[itrk] << " ID: " << t.trkSvtxId[itrk] << " Pt: " << t.trkPt[itrk] << endl;

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

        //remake secondary vertices
        if (t.trkSvtxId[itrk]<0) continue;

        // if (secVtxIDs.count(t.trkSvtxId[itrk])) secVtxIDs[t.trkSvtxId[itrk]].push_back(t.trkMatchSta[itrk]);
        // if (!(secVtxIDs.count(t.trkSvtxId[itrk]))) secVtxIDs.insert({t.trkSvtxId[itrk], {t.trkMatchSta[itrk]}});


        if (secVtxs.count(t.trkSvtxId[itrk])) secVtxs[t.trkSvtxId[itrk]].push_back(v1);
        if (!(secVtxs.count(t.trkSvtxId[itrk]))) secVtxs.insert({t.trkSvtxId[itrk], {v1}});

        // cout << "SEC VTX with ID " << t.trkSvtxId[itrk] << ": " << endl;
        // displayTrackVectors(secVtxs[t.trkSvtxId[itrk]]);
        // cout << "TRACKS not in svtx: " << endl;
        // displayTrackVectors(trackVectors);

        // displayTrackVectors(secVtxs[t.trkSvtxId[itrk]]);
    }
    if (secVtxs.size()<2) return emptySVs;

    // displayTrackVectors(trackVectors);

    // map<Int_t, Int_t> svtxStas;

    // for (auto i : secVtxIDs){
    //     // for (auto j : i.second) cout << j << endl;
    //     for (auto j : i.second){
    //         if (i.second.size()>1 && count(i.second.begin(), i.second.end(), j)<2){
    //             unpure++;
    //             // emergency=1;
    //             // cout << "here" << endl;
    //         }
    //     }
    //     Int_t commonSta = findMostCommon(i.second);
    //     if (commonSta==-999){
    //         Double_t maxPt = 0;
    //         for (Int_t j=0; j < secVtxs[i.first].size(); j++){
    //             if (secVtxs[i.first][j].Pt() > maxPt){
    //                 maxPt = secVtxs[i.first][j].Pt();
    //                 commonSta = secVtxIDs[i.first][j];
    //             }
    //         }
    //     }
    //     // cout << "Common sta: " << commonSta << endl;
    //     // cout << "status: " << commonSta << endl;
    //     svtxStas.insert({i.first, commonSta});
    // }

    // //remove secondary vertices tracks from list of tracks
    // vector<Int_t> indices;
    // for (auto i : secVtxs){
    //     for (Int_t j=0;j<trackVectors.size();j++){
    //         if (count(i.second.begin(), i.second.end(), trackVectors[j])) indices.push_back(j);
    //     }
    // }
    // std::sort(indices.begin(),indices.end(), greater<int>());
    // for (Int_t i = 0; i < indices.size(); i++){
    //     trackVectors.erase(trackVectors.begin() + indices[i]);
    // }

    map<Int_t, ROOT::Math::PtEtaPhiMVector> finalSecVtxs;
    for (auto i : secVtxs){
        ROOT::Math::PtEtaPhiMVector vtxSum;
        for (Int_t j=0; j<i.second.size(); j++) vtxSum+=i.second[j];
        finalSecVtxs.insert({i.first, vtxSum});
        // if (svtxStas[i.first]>1){
            // if (finalSecVtxs.count(svtxStas[i.first])){
            //     h_b2AngleBkgd->Fill(std::pow((finalSecVtxs[svtxStas[i.first]].Rapidity()-vtxSum.Rapidity()),2)+std::pow(acos(cos(finalSecVtxs[svtxStas[i.first]].Phi()-vtxSum.Phi())),2), t.weight);
            //     // finalSecVtxs[svtxStas[i.first]]+=vtxSum;
            //     filled=1;
            //     split++;
            //     finalSecVtxs.insert({9999, vtxSum});
            // }
        // }
    }


    // if (finalSecVtxs.size()<2) continue;
    
    vector<ROOT::Math::PtEtaPhiMVector> vecFinalSecVtxs;
    for (auto i : finalSecVtxs){
        vecFinalSecVtxs.push_back(i.second);
    }

    // displayTrackVectors(vecFinalSecVtxs);
    for (auto i : finalSecVtxs) ntrkVector.insert({(i.second).Pt(),t.svtxNtrk[i.first]});
    if (vecFinalSecVtxs.size()>2) groupVertexes(vecFinalSecVtxs, t, finalSecVtxs, ntrkVector);

    Float_t biggestPt=0;
    Float_t smallestPt=0;
    for (auto i : ntrkVector) if (i.first>biggestPt) biggestPt = i.first;
    for (auto i : ntrkVector) if (i.first<biggestPt) smallestPt = i.first;
    SVntrks1=ntrkVector[biggestPt];
    SVntrks2=ntrkVector[smallestPt];

    // cout << "GROUPED: " << endl;
    // displayTrackVectors(vecFinalSecVtxs);

    return vecFinalSecVtxs;
}

void prepare_data()
{
    // Setup 
    float ptMin = 60.;
    float ptMax = 150.;

    const float missing_value = -1000.;

    // TString indir = "/home/llr/cms/kalipoliti/C10630p1_miniAOD/src/HeavyIonsAnalysis/Configuration/test/rootf/";
    // TString label = "HiForestMiniAOD_HighPU_50000events_conmatch_newVars_truth";

    // TString indir = "/data_CMS/cms/kalipoliti/ttbarMC/highPU/aggrGenNoReco/";
    // TString label = "ttbar_highPU";
    TString indir = "/data_CMS/cms/kalipoliti/herwigMC/dijet/herwig_dijet_official_aggrTMVA/";
    TString label = "merged_HiForestMiniAOD";
    // TString fin = indir + "merged_HiForestMiniAOD.root";
    TString fin = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";

    TString odir = "/data_CMS/cms/sorel/";
    TString fout = odir + label + "_" + ptMin + "_pt_" + ptMax + ".root";   

    // Load tree
    tTree t = tTree(fin);

    // Activate branches of interest
    std::vector<TString> activeBranches = {"nref", "jteta", "jtpt", "jtHadFlav", "jtNbHad","jtNsvtx",
                                           "discr_particleNet_BvsAll",
                                           "ntrk", "trkJetId", "trkSvtxId", "trkPt", "trkEta", "trkPhi",
                                           "trkMatchSta", "trkPdgId",
                                           "nsvtx", "svtxJetId", "svtxNtrk", "svtxm", "svtxmcorr", "svtxpt",
                                           "svtxdl", "svtxdls", "svtxdl2d", "svtxdls2d", "svtxnormchi2",
                                           "weight"
                                           };
    t.SetBranchStatus("*", 0);
    t.SetBranchStatus(activeBranches, 1);

    // Create new file to store signal and background trees
    TFile *foutPtr = new TFile(fout, "recreate");
    std::cout << "(Re)creating " << fout << " file." << std::endl;

    TTree *TreeS = new TTree("TreeS", "Signal tree for MVA");
    TTree *TreeB = new TTree("TreeB", "Background tree for MVA");
    TTree *TreePU = new TTree("TreePU", "Pile-up tree for MVA");

    Float_t SVdist;
    
    Float_t SVmass1;
    Float_t SVpt1;
    Float_t SVntrks1;
    Float_t SVmass2;
    Float_t SVpt2;
    Float_t SVntrks2;

    Float_t jtpt;
    Float_t weight;

    TreeS->Branch("SVdist", &SVdist, "SVdist/F");
    TreeS->Branch("SVmass1", &SVmass1, "SVmass1/F");
    TreeS->Branch("SVpt1", &SVpt1, "SVpt1/F");
    TreeS->Branch("SVntrks1", &SVntrks1, "SVntrks1/F");
    TreeS->Branch("SVmass2", &SVmass2, "SVmass2/F");
    TreeS->Branch("SVpt2", &SVpt2, "SVpt2/F");
    TreeS->Branch("SVntrks2", &SVntrks2, "SVntrks2/F");
    TreeS->Branch("jtpt", &jtpt, "jtpt/F");
    TreeS->Branch("weight", &weight, "weight/F");

    TreeB->Branch("SVdist", &SVdist, "SVdist/F");
    TreeB->Branch("SVmass1", &SVmass1, "SVmass1/F");
    TreeB->Branch("SVpt1", &SVpt1, "SVpt1/F");
    TreeB->Branch("SVntrks1", &SVntrks1, "SVntrks1/F");
    TreeB->Branch("SVmass2", &SVmass2, "SVmass2/F");
    TreeB->Branch("SVpt2", &SVpt2, "SVpt2/F");
    TreeB->Branch("SVntrks2", &SVntrks2, "SVntrks2/F");
    TreeB->Branch("jtpt", &jtpt, "jtpt/F");
    TreeB->Branch("weight", &weight, "weight/F");

    TreePU->Branch("SVdist", &SVdist, "SVdist/F");
    TreePU->Branch("SVmass1", &SVmass1, "SVmass1/F");
    TreePU->Branch("SVpt1", &SVpt1, "SVpt1/F");
    TreePU->Branch("SVntrks1", &SVntrks1, "SVntrks1/F");
    TreePU->Branch("SVmass2", &SVmass2, "SVmass2/F");
    TreePU->Branch("SVpt2", &SVpt2, "SVpt2/F");
    TreePU->Branch("SVntrks2", &SVntrks2, "SVntrks2/F");
    TreePU->Branch("jtpt", &jtpt, "jtpt/F");
    TreePU->Branch("weight", &weight, "weight/F");


    
    int sig = 0;
    int bkg = 0;
    int pu = 0;

    int sigcut=0;
    int bkgdcut=0;

    for (Long64_t ient = 0; ient < t.GetEntries(); ient++) {
        // for debugging purposes 
        //if (ient < 217316) continue;
        // if (ient > 1) break;
            
        // Show progress
        if (ient % 1000000 == 0) {
            std::cout << "i = " << ient << std::endl;
        }

        t.GetEntry(ient);

        for (Int_t ijet = 0; ijet < t.nref; ijet++) {
            weight = t.weight;

            bool skipJet = false;
            skipJet |= (std::abs(t.jteta[ijet]) > 2);
            skipJet |= (t.jtpt[ijet] < ptMin || t.jtpt[ijet] > ptMax);
            bool isBjet = (t.jtHadFlav[ijet] == 5);
            skipJet |= (!isBjet);
            bool passBtag = (t.discr_particleNet_BvsAll[ijet] > 0.9);
            skipJet |= (t.jtNsvtx[ijet]<2);
            skipJet |= (!passBtag);

            if (skipJet) continue;

            vector <ROOT::Math::PtEtaPhiMVector> secondaryVs=makeSvtxs(t,ijet, SVntrks1, SVntrks2);
            if (secondaryVs.size()!=2) continue;

            ROOT::Math::PtEtaPhiMVector leadingSV;
            ROOT::Math::PtEtaPhiMVector subleadingSV;
            if (secondaryVs[0].Pt()>secondaryVs[1].Pt()){
                leadingSV = secondaryVs[0];
                subleadingSV = secondaryVs[1];
            }
            else{
                leadingSV = secondaryVs[1];
                subleadingSV = secondaryVs[0];
            }
            
            SVdist = ROOT::Math::VectorUtil::DeltaR(secondaryVs[0], secondaryVs[1]);
            SVmass1 = leadingSV.M();
            SVpt1 = leadingSV.Pt();
            SVmass2 = subleadingSV.M();
            SVpt2 = subleadingSV.Pt();

            // SVdls;
            // SVnormchi2;
            // SVntrks;

            Int_t nbh = t.jtNbHad[ijet];
            if (nbh == 2) {
                TreeS->Fill();
                sig++;
                if (SVdist<0.1) sigcut++;
            } else {
                if (nbh == 1) {
                    TreeB->Fill();
                    bkg++;
                    if (SVdist<0.1) bkgdcut++;
                } else {
                    // TreeB->Fill();
                    // TreeS->Fill();
                    TreePU->Fill();
                    pu++;
                }
            }

            // for (Int_t itrk = 0; itrk < t.ntrk; itrk++) {
            //     int ijet = t.trkJetId[itrk];
            //     bool skipJet = false;
            //     skipJet |= (std::abs(t.jteta[ijet]) > 2);
            //     skipJet |= (t.jtpt[ijet] < ptMin || t.jtpt[ijet] > ptMax);
            //     bool isBjet = (t.jtHadFlav[ijet] == 5);
            //     skipJet |= (!isBjet);
            //     // bool passBtag = (t.discr_deepFlavour_b[ijet] + t.discr_deepFlavour_bb[ijet] + t.discr_deepFlavour_lepb[ijet]) > 0.9;
            //     // passBtag &= (t.discr_deepFlavour_bb[ijet] < 0.2);
            //     bool passBtag = (t.discr_particleNet_BvsAll[ijet] > 0.9);
            //     skipJet |= (!passBtag);

            //     if (skipJet) continue;

            //     jtpt = t.jtpt[ijet];

            //     Float_t trkPt = t.trkPt[itrk];
            //     trkIp3dSig = t.trkIp3dSig[itrk];
            //     if (trkIp3dSig != trkIp3dSig) trkIp3dSig = missing_value; // nan
            //     // trkIp2dSig = t.trkIp2dSig[itrk];
            //     if (trkIp2dSig != trkIp2dSig) trkIp2dSig = missing_value; // nan
            //     trkDistToAxis = t.trkDistToAxis[itrk];
            //     if (trkDistToAxis != trkDistToAxis) trkDistToAxis = missing_value; // nan
            //     // trkDz = t.trkDz[itrk];
            //     if (trkDz != trkDz) trkDz = missing_value; // nan

            //     // if (std::abs(trkDz) > 0.6) continue;

            //     trkPtOverJet = trkPt / jtpt;
            //     trkPdgId = t.trkPdgId[itrk]; 
            //     if (std::abs(trkPdgId) == 11 || std::abs(trkPdgId) == 13) {
            //         trkIsLepton = 1;
            //     } else {
            //         trkIsLepton = 0;
            //     }
            
            //     int trkSvtxId = t.trkSvtxId[itrk];
            //     trkInSV = int(trkSvtxId >= 0);
            //     if (trkSvtxId >= 0) {
            //         svtxdls = t.svtxdls[trkSvtxId];
            //         svtxdls2d = t.svtxdls2d[trkSvtxId];
            //         svtxm = t.svtxm[trkSvtxId];
            //         svtxmcorr = t.svtxmcorr[trkSvtxId];
            //         svtxNtrk = t.svtxNtrk[trkSvtxId];
            //         svtxnormchi2 = t.svtxnormchi2[trkSvtxId];
            //         svtxTrkPtOverSv = trkPt / t.svtxpt[trkSvtxId];
            //     } else {
            //         svtxdls = missing_value;
            //         svtxdls2d = missing_value;
            //         svtxm = missing_value;
            //         svtxmcorr = missing_value;
            //         svtxNtrk = missing_value;
            //         svtxnormchi2 = missing_value;
            //         svtxTrkPtOverSv = missing_value;
            //     }

            //     Int_t sta = t.trkMatchSta[itrk];
            //     if (sta >= 100) {
            //         TreeS->Fill();
            //         sig++;
            //     } else {
            //         if (sta == 1) {
            //             TreeB->Fill();
            //             bkg++;
            //         } else {
            //             // TreeB->Fill();
            //             // TreeS->Fill();
            //             TreePU->Fill();
            //             pu++;
            //         }
            //     }

            // } // end track loop
        }
    } // end entry loop

    int total = sig + bkg + pu;

    std::cout << "sig / total : " << (float) sig / total << std::endl;
    std::cout << "bkg / total : " << (float) bkg / total << std::endl;
    std::cout << "pu / total : " << (float) pu / total << std::endl;

    std::cout << "bkg / sig : " << (float) bkg / sig << std::endl; 
    std::cout << "pu / sig : " << (float) pu / sig << std::endl; 

    std::cout << "sigcut / sig : " << (float) sigcut / sig << std::endl; 
    std::cout << "bkgdcut / bkg : " << (float) bkgdcut / bkg << std::endl; 

    TreeS->Write("", TObject::kOverwrite);
    TreeB->Write("", TObject::kOverwrite);
    TreePU->Write("", TObject::kOverwrite);

    foutPtr->Close();
    delete foutPtr;         
}


