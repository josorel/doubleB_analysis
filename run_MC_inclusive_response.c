#include "tTree.h"
#include "binning.h"
#include <random>

bool skipMC(Float_t jtpt, Float_t refpt, Float_t pthat) {
    if (!(refpt>0)) return true;    
    if (pthat<0.35*jtpt) return true;
    return false;
}

// std::pair<double,double> computePull(const ROOT::Math::PtEtaPhiMVector& jet, const std::vector<ROOT::Math::PtEtaPhiMVector>& constituents) {
//     double pullX = 0, pullY = 0;
//     double jetPt = jet.Pt();

//     for (const auto& c : constituents) {
//         double dEta = c.Rapidity() - jet.Rapidity();
//         double dPhi = TVector2::Phi_mpi_pi(c.Phi() - jet.Phi());
//         double r    = std::hypot(dEta, dPhi);      // = |r_i|
//         if (r == 0) continue;                      // skip exact center
//         pullX += c.Pt() * r * dEta;
//         pullY += c.Pt() * r * dPhi;
//     }

//     pullX /= jetPt;
//     pullY /= jetPt;
//     return {pullX, pullY};
// }

void fill_jk_resampling(std::vector<TH2D *> histos, double num, double x, double y, double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 histograms
        histos[i]->Fill(x,y,w);
    }
}

void fill_jk_resampling_response(std::vector<RooUnfoldResponse *> responses, double num, 
                                 double x_reco, double y_reco, 
                                 double x_gen, double y_gen, 
                                 double w) {
    /*
    histos: vector of size 10
    num:  double between 0 and 1
    */
    
    int numRescaled = (int) (num*10); // from 0 to 9
    for (int i=0; i<10; i++) {
        if (i==numRescaled) continue; // fill all but 1 responses
        responses[i]->Fill(x_reco,y_reco,x_gen,y_gen,w);
    }
}

void create_response_fast(Float_t& jtpt,Float_t& jtptCh,auto& logrg,auto& logkt,auto& zg,auto& logrg_gen,auto& logkt_gen,auto& zg_gen,Float_t& weight,Float_t& pthat,
                          Float_t& jtptCh_gen, Float_t& bpt, Float_t& bpt_gen, Float_t& jtpt_gen,
                          TString jer_opt="nom", TString jec_opt="nom")
{

    // random number generator for jackknife resampling
    const double range_from = 0;
    const double range_to = 1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distr(range_from, range_to);

    // Random number for jack-knife resampling
    double num = distr(generator);

    // if (jer_opt=="nom") jtpt = jtpt * jer_sf_nom;
    // else if (jer_opt=="up") jtpt = jtpt * jer_sf_up;
    // else if (jer_opt=="down") jtpt = jtpt * jer_sf_down;

    // double jec_fact = 0; // nominal
    // if (jec_opt=="up") jec_fact = 1;
    // if (jec_opt=="down") jec_fact = -1;
    // double newCorrection = 1 + (jec_fact * jec_unc);
    // jtpt = jtpt * newCorrection;

    // ---- Fix observable limits ----

    // Checks for logrg 
    if (logrg < 0. && logrg > -900.) return; // rg>0.4 -> skip for now
    if (logrg < -900.) logrg = -0.2; // SD-untagged bin range
    if (logkt < 0. && logkt > -900.) logrg = -0.2; // kt<1 -> SD-untagged bin range
    if (logrg >= 2.5) logrg = 2.499; // overflow bin range
    if (logrg_gen < 0. && logrg_gen>-900) logrg_gen = -900; // rg_gen>0.4 -> out of range but not SD-untagged
    if (logrg_gen < -900.) logrg_gen = -0.2; // SD-untagged bin range
    if (logkt_gen < 0. && logkt_gen > -900.) logrg_gen = -0.2; // kt_gen<1 -> SD-untagged bin range
    if (logrg_gen >= 2.5) logrg_gen = 2.499; // overflow bin range

    // Checks for zg
    if (zg < -900.) zg = 0.; // SD-untagged bin range
    if (logkt < 0. && logkt > -900.) zg = 0.; // kt<1 -> SD-untagged bin range
    if (zg >= 0.5) zg = 0.499; // zg=0.5 included in last bin
    if (zg_gen < -900.) zg_gen = 0.; // SD-untagged bin range
    if (logkt_gen < 0. && logkt_gen > -900.) zg_gen = 0.; // kt<1 -> SD-untagged bin range
    if (zg_gen >= 0.5) zg_gen = 0.499; // zg=0.5 included in last bin

    // Checks for kt
    if (std::abs(jtptCh) < 1e-4) return;
    // --------------------------

    // Check if pass cuts
    bool has_gen_match = (jtpt_gen > 0);

    bool reco_pass_cuts_rg=false;
    bool true_pass_cuts_rg=false;

    bool reco_pass_cuts_zg=false;
    bool true_pass_cuts_zg=false;

    bool reco_pass_cuts_kt=false;
    bool true_pass_cuts_kt=false;

    if (zg>0.0){
        reco_pass_cuts_rg=true;
        reco_pass_cuts_zg=true;
        reco_pass_cuts_kt=true;
    }
    if (zg_gen>0.0){
        true_pass_cuts_rg=true;
        true_pass_cuts_zg=true;
        true_pass_cuts_kt=true;
    }

    // bool reco_pass_cuts_rg = (jtpt < jtpt_max && jtpt >= jtpt_min && logrg < logrg_max && logrg >= logrg_min);
    // bool true_pass_cuts_rg = (jtpt_gen < jtpt_max && jtpt_gen >= jtpt_min && logrg_gen < logrg_max && logrg_gen >= logrg_min);

    // bool reco_pass_cuts_zg = (jtpt < jtpt_max && jtpt >= jtpt_min && zg < zg_max && zg >= zg_min);
    // bool true_pass_cuts_zg = (jtpt_gen < jtpt_max && jtpt_gen >= jtpt_min && zg_gen < zg_max && zg_gen >= zg_min);

    // bool reco_pass_cuts_kt = (jtpt < jtpt_max && jtpt >= jtpt_min && kt < kt_max && kt >= kt_min);
    // bool true_pass_cuts_kt = (jtpt_gen < jtpt_max && jtpt_gen >= jtpt_min && kt_gen < kt_max && kt_gen >= kt_min);

    // Fill histograms
    // if (!has_gen_match) {   
    //     // fill fakes
    //     return; 
    // } 
    
    // The rest of the histograms don;t include any fakes

    // fill rg histograms
    if (true_pass_cuts_rg) {
        fill_jk_resampling(histos_efficiency_denominator_rgpt, num, logrg_gen, jtpt_gen, weight);
        if (num<0.5) h_half0_efficiency_denominator_rgpt->Fill(logrg_gen, jtpt_gen, weight);
        else h_half1_efficiency_denominator_rgpt->Fill(logrg_gen, jtpt_gen, weight);
    }
    if (reco_pass_cuts_rg) {
        fill_jk_resampling(histos_purity_denominator_rgpt, num, logrg, jtpt, weight);
        if (num<0.5) h_half0_purity_denominator_rgpt->Fill(logrg, jtpt, weight);
        else h_half1_purity_denominator_rgpt->Fill(logrg, jtpt, weight);
    }
    if (true_pass_cuts_rg && reco_pass_cuts_rg) {
        fill_jk_resampling(histos_efficiency_numerator_rgpt, num, logrg_gen, jtpt_gen, weight);
        fill_jk_resampling(histos_purity_numerator_rgpt, num, logrg, jtpt, weight);
        fill_jk_resampling_response(responses_rgpt, num, logrg, jtpt, logrg_gen, jtpt_gen, weight);

        if (num<0.5) {
            h_half0_efficiency_numerator_rgpt->Fill(logrg_gen, jtpt_gen, weight);
            h_half0_purity_numerator_rgpt->Fill(logrg, jtpt, weight);
            response_half0_rgpt->Fill(logrg, jtpt, logrg_gen, jtpt_gen, weight);
        } else {
            h_half1_efficiency_numerator_rgpt->Fill(logrg_gen, jtpt_gen, weight);
            h_half1_purity_numerator_rgpt->Fill(logrg, jtpt, weight);
            response_half1_rgpt->Fill(logrg, jtpt, logrg_gen, jtpt_gen, weight);
        }
        response_full_rgpt->Fill(logrg, jtpt, logrg_gen, jtpt_gen, weight);
    }
    
    // fill zg histograms
    if (true_pass_cuts_zg) {
        fill_jk_resampling(histos_efficiency_denominator_zgpt, num, zg_gen, jtpt_gen, weight);
        if (num<0.5) h_half0_efficiency_denominator_zgpt->Fill(zg_gen, jtpt_gen, weight);
        else h_half1_efficiency_denominator_zgpt->Fill(zg_gen, jtpt_gen, weight);
    }
    if (reco_pass_cuts_zg) {
        fill_jk_resampling(histos_purity_denominator_zgpt, num, zg, jtpt, weight);
        if (num<0.5) h_half0_purity_denominator_zgpt->Fill(zg, jtpt, weight);
        else h_half1_purity_denominator_zgpt->Fill(zg, jtpt, weight);
    }
    if (true_pass_cuts_zg && reco_pass_cuts_zg) {
        fill_jk_resampling(histos_efficiency_numerator_zgpt, num, zg_gen, jtpt_gen, weight);
        fill_jk_resampling(histos_purity_numerator_zgpt, num, zg, jtpt, weight);
        fill_jk_resampling_response(responses_zgpt, num, zg, jtpt, zg_gen, jtpt_gen, weight);

        if (num<0.5) {
            h_half0_efficiency_numerator_zgpt->Fill(zg_gen, jtpt_gen, weight);
            h_half0_purity_numerator_zgpt->Fill(zg, jtpt, weight);
            response_half0_zgpt->Fill(zg, jtpt, zg_gen, jtpt_gen, weight);
        } else {
            h_half1_efficiency_numerator_zgpt->Fill(zg_gen, jtpt_gen, weight);
            h_half1_purity_numerator_zgpt->Fill(zg, jtpt, weight);
            response_half1_zgpt->Fill(zg, jtpt, zg_gen, jtpt_gen, weight);
        }
        response_full_zgpt->Fill(zg, jtpt, zg_gen, jtpt_gen, weight);
    }
    // fill kt histograms
    if (true_pass_cuts_kt) {
        fill_jk_resampling(histos_efficiency_denominator_ktpt, num, logkt_gen, jtpt_gen, weight);
        if (num<0.5) h_half0_efficiency_denominator_ktpt->Fill(logkt_gen, jtpt_gen, weight);
        else h_half1_efficiency_denominator_ktpt->Fill(logkt_gen, jtpt_gen, weight);
    }
    if (reco_pass_cuts_kt) {
        fill_jk_resampling(histos_purity_denominator_ktpt, num, logkt, jtpt, weight);
        if (num<0.5) h_half0_purity_denominator_ktpt->Fill(logkt, jtpt, weight);
        else h_half1_purity_denominator_ktpt->Fill(logkt, jtpt, weight);
    }
    if (true_pass_cuts_kt && reco_pass_cuts_kt) {
        fill_jk_resampling(histos_efficiency_numerator_ktpt, num, logkt_gen, jtpt_gen, weight);
        fill_jk_resampling(histos_purity_numerator_ktpt, num, logkt, jtpt, weight);
        fill_jk_resampling_response(responses_ktpt, num, logkt, jtpt, logkt_gen, jtpt_gen, weight);

        if (num<0.5) {
            h_half0_efficiency_numerator_ktpt->Fill(logkt_gen, jtpt_gen, weight);
            h_half0_purity_numerator_ktpt->Fill(logkt, jtpt, weight);
            response_half0_ktpt->Fill(logkt, jtpt, logkt_gen, jtpt_gen, weight);
        } else {
            h_half1_efficiency_numerator_ktpt->Fill(logkt_gen, jtpt_gen, weight);
            h_half1_purity_numerator_ktpt->Fill(logkt, jtpt, weight);
            response_half1_ktpt->Fill(logkt, jtpt, logkt_gen, jtpt_gen, weight);
        }
        response_full_ktpt->Fill(logkt, jtpt, logkt_gen, jtpt_gen, weight);
    }    
}

void purityEfficiencyHists(){
    // Create purity and efficiency histograms
    TH2D *h0_purity_rgpt, 
        *h1_purity_rgpt,
        *h2_purity_rgpt,
        *h3_purity_rgpt,
        *h4_purity_rgpt,
        *h5_purity_rgpt,
        *h6_purity_rgpt,
        *h7_purity_rgpt,
        *h8_purity_rgpt,
        *h9_purity_rgpt;

    std::vector<TH2D *> histos_purity_rgpt = {
        h0_purity_rgpt,
        h1_purity_rgpt,
        h2_purity_rgpt,
        h3_purity_rgpt,
        h4_purity_rgpt,
        h5_purity_rgpt,
        h6_purity_rgpt,
        h7_purity_rgpt,
        h8_purity_rgpt,
        h9_purity_rgpt,
    };

    TH2D *h0_efficiency_rgpt, 
        *h1_efficiency_rgpt,
        *h2_efficiency_rgpt,
        *h3_efficiency_rgpt,
        *h4_efficiency_rgpt,
        *h5_efficiency_rgpt,
        *h6_efficiency_rgpt,
        *h7_efficiency_rgpt,
        *h8_efficiency_rgpt,
        *h9_efficiency_rgpt;

    std::vector<TH2D *> histos_efficiency_rgpt = {
        h0_efficiency_rgpt,
        h1_efficiency_rgpt,
        h2_efficiency_rgpt,
        h3_efficiency_rgpt,
        h4_efficiency_rgpt,
        h5_efficiency_rgpt,
        h6_efficiency_rgpt,
        h7_efficiency_rgpt,
        h8_efficiency_rgpt,
        h9_efficiency_rgpt,
    };

    TH2D *h0_purity_zgpt, 
        *h1_purity_zgpt,
        *h2_purity_zgpt,
        *h3_purity_zgpt,
        *h4_purity_zgpt,
        *h5_purity_zgpt,
        *h6_purity_zgpt,
        *h7_purity_zgpt,
        *h8_purity_zgpt,
        *h9_purity_zgpt;

    std::vector<TH2D *> histos_purity_zgpt = {
        h0_purity_zgpt,
        h1_purity_zgpt,
        h2_purity_zgpt,
        h3_purity_zgpt,
        h4_purity_zgpt,
        h5_purity_zgpt,
        h6_purity_zgpt,
        h7_purity_zgpt,
        h8_purity_zgpt,
        h9_purity_zgpt,
    };

    TH2D *h0_efficiency_zgpt, 
        *h1_efficiency_zgpt,
        *h2_efficiency_zgpt,
        *h3_efficiency_zgpt,
        *h4_efficiency_zgpt,
        *h5_efficiency_zgpt,
        *h6_efficiency_zgpt,
        *h7_efficiency_zgpt,
        *h8_efficiency_zgpt,
        *h9_efficiency_zgpt;

    std::vector<TH2D *> histos_efficiency_zgpt = {
        h0_efficiency_zgpt,
        h1_efficiency_zgpt,
        h2_efficiency_zgpt,
        h3_efficiency_zgpt,
        h4_efficiency_zgpt,
        h5_efficiency_zgpt,
        h6_efficiency_zgpt,
        h7_efficiency_zgpt,
        h8_efficiency_zgpt,
        h9_efficiency_zgpt,
    };

    TH2D *h0_purity_ktpt, 
        *h1_purity_ktpt,
        *h2_purity_ktpt,
        *h3_purity_ktpt,
        *h4_purity_ktpt,
        *h5_purity_ktpt,
        *h6_purity_ktpt,
        *h7_purity_ktpt,
        *h8_purity_ktpt,
        *h9_purity_ktpt;

    std::vector<TH2D *> histos_purity_ktpt = {
        h0_purity_ktpt,
        h1_purity_ktpt,
        h2_purity_ktpt,
        h3_purity_ktpt,
        h4_purity_ktpt,
        h5_purity_ktpt,
        h6_purity_ktpt,
        h7_purity_ktpt,
        h8_purity_ktpt,
        h9_purity_ktpt,
    };

    TH2D *h0_efficiency_ktpt, 
        *h1_efficiency_ktpt,
        *h2_efficiency_ktpt,
        *h3_efficiency_ktpt,
        *h4_efficiency_ktpt,
        *h5_efficiency_ktpt,
        *h6_efficiency_ktpt,
        *h7_efficiency_ktpt,
        *h8_efficiency_ktpt,
        *h9_efficiency_ktpt;

    std::vector<TH2D *> histos_efficiency_ktpt = {
        h0_efficiency_ktpt,
        h1_efficiency_ktpt,
        h2_efficiency_ktpt,
        h3_efficiency_ktpt,
        h4_efficiency_ktpt,
        h5_efficiency_ktpt,
        h6_efficiency_ktpt,
        h7_efficiency_ktpt,
        h8_efficiency_ktpt,
        h9_efficiency_ktpt,
    };
    
    // initialize histograms
    for (int i=0; i<10; i++) {
        histos_purity_rgpt[i] = (TH2D *) histos_purity_numerator_rgpt[i]->Clone("h"+TString(Form("%d",i))+"_purity_rgpt");
        histos_purity_rgpt[i]->Divide(histos_purity_numerator_rgpt[i], histos_purity_denominator_rgpt[i], 1., 1., "b");
        histos_efficiency_rgpt[i] = (TH2D *) histos_efficiency_numerator_rgpt[i]->Clone("h"+TString(Form("%d",i))+"_efficiency_rgpt");
        histos_efficiency_rgpt[i]->Divide(histos_efficiency_numerator_rgpt[i], histos_efficiency_denominator_rgpt[i], 1., 1., "b");

        histos_purity_zgpt[i] = (TH2D *) histos_purity_numerator_zgpt[i]->Clone("h"+TString(Form("%d",i))+"_purity_zgpt");
        histos_purity_zgpt[i]->Divide(histos_purity_numerator_zgpt[i], histos_purity_denominator_zgpt[i], 1., 1., "b");
        histos_efficiency_zgpt[i] = (TH2D *) histos_efficiency_numerator_zgpt[i]->Clone("h"+TString(Form("%d",i))+"_efficiency_zgpt");
        histos_efficiency_zgpt[i]->Divide(histos_efficiency_numerator_zgpt[i], histos_efficiency_denominator_zgpt[i], 1., 1., "b");

        histos_purity_ktpt[i] = (TH2D *) histos_purity_numerator_ktpt[i]->Clone("h"+TString(Form("%d",i))+"_purity_ktpt");
        histos_purity_ktpt[i]->Divide(histos_purity_numerator_ktpt[i], histos_purity_denominator_ktpt[i], 1., 1., "b");
        histos_efficiency_ktpt[i] = (TH2D *) histos_efficiency_numerator_ktpt[i]->Clone("h"+TString(Form("%d",i))+"_efficiency_ktpt");
        histos_efficiency_ktpt[i]->Divide(histos_efficiency_numerator_ktpt[i], histos_efficiency_denominator_ktpt[i], 1., 1., "b");
    }

    // declare the per half purity + efficiency histograms
    TH2D *h_half0_purity_rgpt = (TH2D *) h_half0_purity_numerator_rgpt->Clone("h_half0_purity_rgpt");
    h_half0_purity_rgpt->Divide(h_half0_purity_numerator_rgpt, h_half0_purity_denominator_rgpt, 1., 1., "b");
    TH2D *h_half0_efficiency_rgpt = (TH2D *) h_half0_efficiency_numerator_rgpt->Clone("h_half0_efficiency_rgpt");
    h_half0_efficiency_rgpt->Divide(h_half0_efficiency_numerator_rgpt, h_half0_efficiency_denominator_rgpt, 1., 1., "b");

    TH2D *h_half1_purity_rgpt = (TH2D *) h_half1_purity_numerator_rgpt->Clone("h_half1_purity_rgpt");
    h_half1_purity_rgpt->Divide(h_half1_purity_numerator_rgpt, h_half1_purity_denominator_rgpt, 1., 1., "b");
    TH2D *h_half1_efficiency_rgpt = (TH2D *) h_half1_efficiency_numerator_rgpt->Clone("h_half1_efficiency_rgpt");
    h_half1_efficiency_rgpt->Divide(h_half1_efficiency_numerator_rgpt, h_half1_efficiency_denominator_rgpt, 1., 1., "b");

    TH2D *h_half0_purity_zgpt = (TH2D *) h_half0_purity_numerator_zgpt->Clone("h_half0_purity_zgpt");
    h_half0_purity_zgpt->Divide(h_half0_purity_numerator_zgpt, h_half0_purity_denominator_zgpt, 1., 1., "b");
    TH2D *h_half0_efficiency_zgpt = (TH2D *) h_half0_efficiency_numerator_zgpt->Clone("h_half0_efficiency_zgpt");
    h_half0_efficiency_zgpt->Divide(h_half0_efficiency_numerator_zgpt, h_half0_efficiency_denominator_zgpt, 1., 1., "b");

    TH2D *h_half1_purity_zgpt = (TH2D *) h_half1_purity_numerator_zgpt->Clone("h_half1_purity_zgpt");
    h_half1_purity_zgpt->Divide(h_half1_purity_numerator_zgpt, h_half1_purity_denominator_zgpt, 1., 1., "b");
    TH2D *h_half1_efficiency_zgpt = (TH2D *) h_half1_efficiency_numerator_zgpt->Clone("h_half1_efficiency_zgpt");
    h_half1_efficiency_zgpt->Divide(h_half1_efficiency_numerator_zgpt, h_half1_efficiency_denominator_zgpt, 1., 1., "b");

    TH2D *h_half0_purity_ktpt = (TH2D *) h_half0_purity_numerator_ktpt->Clone("h_half0_purity_ktpt");
    h_half0_purity_ktpt->Divide(h_half0_purity_numerator_ktpt, h_half0_purity_denominator_ktpt, 1., 1., "b");
    TH2D *h_half0_efficiency_ktpt = (TH2D *) h_half0_efficiency_numerator_ktpt->Clone("h_half0_efficiency_ktpt");
    h_half0_efficiency_ktpt->Divide(h_half0_efficiency_numerator_ktpt, h_half0_efficiency_denominator_ktpt, 1., 1., "b");

    TH2D *h_half1_purity_ktpt = (TH2D *) h_half1_purity_numerator_ktpt->Clone("h_half1_purity_ktpt");
    h_half1_purity_ktpt->Divide(h_half1_purity_numerator_ktpt, h_half1_purity_denominator_ktpt, 1., 1., "b");
    TH2D *h_half1_efficiency_ktpt = (TH2D *) h_half1_efficiency_numerator_ktpt->Clone("h_half1_efficiency_ktpt");
    h_half1_efficiency_ktpt->Divide(h_half1_efficiency_numerator_ktpt, h_half1_efficiency_denominator_ktpt, 1., 1., "b");

    // declare the full purity + efficiency histograms
    TH2D *h_full_purity_numerator_rgpt = (TH2D *) h_half0_purity_numerator_rgpt->Clone("h_full_purity_numerator_rgpt");
    h_full_purity_numerator_rgpt->Add(h_half1_purity_numerator_rgpt);
    TH2D *h_full_purity_denominator_rgpt = (TH2D *) h_half0_purity_denominator_rgpt->Clone("h_full_purity_denominator_rgpt");
    h_full_purity_denominator_rgpt->Add(h_half1_purity_denominator_rgpt);
    TH2D *h_full_purity_rgpt = (TH2D *) h_full_purity_numerator_rgpt->Clone("h_full_purity_rgpt");
    h_full_purity_rgpt->Divide(h_full_purity_numerator_rgpt, h_full_purity_denominator_rgpt, 1., 1., "b");

    TH2D *h_full_efficiency_numerator_rgpt = (TH2D *) h_half0_efficiency_numerator_rgpt->Clone("h_full_efficiency_numerator_rgpt");
    h_full_efficiency_numerator_rgpt->Add(h_half1_efficiency_numerator_rgpt);
    TH2D *h_full_efficiency_denominator_rgpt = (TH2D *) h_half0_efficiency_denominator_rgpt->Clone("h_full_efficiency_denominator_rgpt");
    h_full_efficiency_denominator_rgpt->Add(h_half1_efficiency_denominator_rgpt);
    TH2D *h_full_efficiency_rgpt = (TH2D *) h_full_efficiency_numerator_rgpt->Clone("h_full_efficiency_rgpt");
    h_full_efficiency_rgpt->Divide(h_full_efficiency_numerator_rgpt, h_full_efficiency_denominator_rgpt, 1., 1., "b");

    TH2D *h_full_purity_numerator_zgpt = (TH2D *) h_half0_purity_numerator_zgpt->Clone("h_full_purity_numerator_zgpt");
    h_full_purity_numerator_zgpt->Add(h_half1_purity_numerator_zgpt);
    TH2D *h_full_purity_denominator_zgpt = (TH2D *) h_half0_purity_denominator_zgpt->Clone("h_full_purity_denominator_zgpt");
    h_full_purity_denominator_zgpt->Add(h_half1_purity_denominator_zgpt);
    TH2D *h_full_purity_zgpt = (TH2D *) h_full_purity_numerator_zgpt->Clone("h_full_purity_zgpt");
    h_full_purity_zgpt->Divide(h_full_purity_numerator_zgpt, h_full_purity_denominator_zgpt, 1., 1., "b");

    TH2D *h_full_efficiency_numerator_zgpt = (TH2D *) h_half0_efficiency_numerator_zgpt->Clone("h_full_efficiency_numerator_zgpt");
    h_full_efficiency_numerator_zgpt->Add(h_half1_efficiency_numerator_zgpt);
    TH2D *h_full_efficiency_denominator_zgpt = (TH2D *) h_half0_efficiency_denominator_zgpt->Clone("h_full_efficiency_denominator_zgpt");
    h_full_efficiency_denominator_zgpt->Add(h_half1_efficiency_denominator_zgpt);
    TH2D *h_full_efficiency_zgpt = (TH2D *) h_full_efficiency_numerator_zgpt->Clone("h_full_efficiency_zgpt");
    h_full_efficiency_zgpt->Divide(h_full_efficiency_numerator_zgpt, h_full_efficiency_denominator_zgpt, 1., 1., "b");

    TH2D *h_full_purity_numerator_ktpt = (TH2D *) h_half0_purity_numerator_ktpt->Clone("h_full_purity_numerator_ktpt");
    h_full_purity_numerator_ktpt->Add(h_half1_purity_numerator_ktpt);
    TH2D *h_full_purity_denominator_ktpt = (TH2D *) h_half0_purity_denominator_ktpt->Clone("h_full_purity_denominator_ktpt");
    h_full_purity_denominator_ktpt->Add(h_half1_purity_denominator_ktpt);
    TH2D *h_full_purity_ktpt = (TH2D *) h_full_purity_numerator_ktpt->Clone("h_full_purity_ktpt");
    h_full_purity_ktpt->Divide(h_full_purity_numerator_ktpt, h_full_purity_denominator_ktpt, 1., 1., "b");

    TH2D *h_full_efficiency_numerator_ktpt = (TH2D *) h_half0_efficiency_numerator_ktpt->Clone("h_full_efficiency_numerator_ktpt");
    h_full_efficiency_numerator_ktpt->Add(h_half1_efficiency_numerator_ktpt);
    TH2D *h_full_efficiency_denominator_ktpt = (TH2D *) h_half0_efficiency_denominator_ktpt->Clone("h_full_efficiency_denominator_ktpt");
    h_full_efficiency_denominator_ktpt->Add(h_half1_efficiency_denominator_ktpt);
    TH2D *h_full_efficiency_ktpt = (TH2D *) h_full_efficiency_numerator_ktpt->Clone("h_full_efficiency_ktpt");
    h_full_efficiency_ktpt->Divide(h_full_efficiency_numerator_ktpt, h_full_efficiency_denominator_ktpt, 1., 1., "b");

    // Create output file
    TString fout_name = "responseMatrixBjetInclusiveMC_JERdown.root";
    std::cout << "Creating file: " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    // Write jk resampling histograms + responses
    for (int i=0; i<10; i++) {
        histos_purity_numerator_rgpt[i]->Write();
        histos_purity_denominator_rgpt[i]->Write();
        histos_purity_rgpt[i]->Write();

        histos_efficiency_numerator_rgpt[i]->Write();
        histos_efficiency_denominator_rgpt[i]->Write();
        histos_efficiency_rgpt[i]->Write();

        responses_rgpt[i]->Write();

        histos_purity_numerator_zgpt[i]->Write();
        histos_purity_denominator_zgpt[i]->Write();
        histos_purity_zgpt[i]->Write();

        histos_efficiency_numerator_zgpt[i]->Write();
        histos_efficiency_denominator_zgpt[i]->Write();
        histos_efficiency_zgpt[i]->Write();

        responses_zgpt[i]->Write();

        histos_purity_numerator_ktpt[i]->Write();
        histos_purity_denominator_ktpt[i]->Write();
        histos_purity_ktpt[i]->Write();

        histos_efficiency_numerator_ktpt[i]->Write();
        histos_efficiency_denominator_ktpt[i]->Write();
        histos_efficiency_ktpt[i]->Write();

        responses_ktpt[i]->Write();
    }

    // Write per half histograms 
    h_half0_purity_numerator_rgpt->Write();
    h_half0_purity_denominator_rgpt->Write();
    h_half0_purity_rgpt->Write();

    h_half0_efficiency_numerator_rgpt->Write();
    h_half0_efficiency_denominator_rgpt->Write();
    h_half0_efficiency_rgpt->Write();

    response_half0_rgpt->Write();

    h_half1_purity_numerator_rgpt->Write();
    h_half1_purity_denominator_rgpt->Write();
    h_half1_purity_rgpt->Write();

    h_half1_efficiency_numerator_rgpt->Write();
    h_half1_efficiency_denominator_rgpt->Write();
    h_half1_efficiency_rgpt->Write();

    response_half1_rgpt->Write();

    h_half0_purity_numerator_zgpt->Write();
    h_half0_purity_denominator_zgpt->Write();
    h_half0_purity_zgpt->Write();

    h_half0_efficiency_numerator_zgpt->Write();
    h_half0_efficiency_denominator_zgpt->Write();
    h_half0_efficiency_zgpt->Write();

    response_half0_zgpt->Write();

    h_half1_purity_numerator_zgpt->Write();
    h_half1_purity_denominator_zgpt->Write();
    h_half1_purity_zgpt->Write();

    h_half1_efficiency_numerator_zgpt->Write();
    h_half1_efficiency_denominator_zgpt->Write();
    h_half1_efficiency_zgpt->Write();

    response_half1_zgpt->Write();

    h_half0_purity_numerator_ktpt->Write();
    h_half0_purity_denominator_ktpt->Write();
    h_half0_purity_ktpt->Write();

    h_half0_efficiency_numerator_ktpt->Write();
    h_half0_efficiency_denominator_ktpt->Write();
    h_half0_efficiency_ktpt->Write();

    response_half0_ktpt->Write();

    h_half1_purity_numerator_ktpt->Write();
    h_half1_purity_denominator_ktpt->Write();
    h_half1_purity_ktpt->Write();

    h_half1_efficiency_numerator_ktpt->Write();
    h_half1_efficiency_denominator_ktpt->Write();
    h_half1_efficiency_ktpt->Write();

    response_half1_ktpt->Write();
    
    // Write full histograms 
    h_full_purity_numerator_rgpt->Write();
    h_full_purity_denominator_rgpt->Write();
    h_full_purity_rgpt->Write();

    h_full_efficiency_numerator_rgpt->Write();
    h_full_efficiency_denominator_rgpt->Write();
    h_full_efficiency_rgpt->Write();

    response_full_rgpt->Write();

    h_full_purity_numerator_zgpt->Write();
    h_full_purity_denominator_zgpt->Write();
    h_full_purity_zgpt->Write();

    h_full_efficiency_numerator_zgpt->Write();
    h_full_efficiency_denominator_zgpt->Write();
    h_full_efficiency_zgpt->Write();

    response_full_zgpt->Write();

    h_full_purity_numerator_ktpt->Write();
    h_full_purity_denominator_ktpt->Write();
    h_full_purity_ktpt->Write();

    h_full_efficiency_numerator_ktpt->Write();
    h_full_efficiency_denominator_ktpt->Write();
    h_full_efficiency_ktpt->Write();

    response_full_ktpt->Write();

    fout->Close();
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
    std::cout << "Displaying Track Vectors:" << std::endl;

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
    // std::cout << "TOTAL PT = " << totalPt << endl;
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

// // This method groups closest reco B hadrons (using GEN info) until only two remain
// void groupBs(map<int, ROOT::Math::PtEtaPhiMVector>& mp){

//     std::vector<ROOT::Math::PtEtaPhiMVector> groupedBs;
//     for (auto i : mp) groupedBs.push_back(i.second);

//     while (groupedBs.size() > 2){
//         double min_distance=std::numeric_limits<double>::infinity();
//         Int_t index1 = 0, index2 = 0;
//         //loop through given tracks to find those with smallest distance
//         for (Int_t i = 0; i < groupedBs.size(); ++i) {
//             for (Int_t j = i + 1; j < groupedBs.size(); ++j) {
//                 Float_t dist = std::pow((groupedBs[i].Eta()-groupedBs[j].Eta()),2)+std::pow(acos(cos(groupedBs[i].Phi()-groupedBs[j].Phi())),2);
//                 if (dist < min_distance) {
//                     min_distance = dist;
//                     index1 = i;
//                     index2 = j;
//                 }
//             }
//         }

//         ROOT::Math::PtEtaPhiMVector newVector = groupedBs[index1]+groupedBs[index2];

//         // Add combined track
//         groupedBs.push_back(newVector);
//         // remove merged tracks
//         if (index1 > index2) std::swap(index1, index2);
//         groupedBs.erase(groupedBs.begin() + index2);
//         groupedBs.erase(groupedBs.begin() + index1);

//     }

//     mp.clear();
//     Int_t extra=0;
//     for (auto i : groupedBs){
//         mp.insert({extra,i});
//         extra++;
//     }
// }

// // This method reconstructs B hadrons using GEN track info (track statuses)
// void ReconstructBsGen(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, tTree& t, map<int, ROOT::Math::PtEtaPhiMVector>& mp, Int_t& ijet){
//     std::vector<Int_t> indices;

//     // You can loop over the tracks of a specific jet using trkJetId
//     for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
//         if (t.refTrkJetId[itrk] != ijet) continue; 

//         //if it is from a b hadron, create vector v1 for track and set properties
//         if (t.refTrkSta[itrk]>1){

//             ROOT::Math::PtEtaPhiMVector v1;

//             v1.SetEta(t.refTrkEta[itrk]);
//             v1.SetPt(t.refTrkPt[itrk]);
//             v1.SetPhi(t.refTrkPhi[itrk]);

//             if(std::abs(t.refTrkPdgId[itrk])==211){
//                 v1.SetM(0.139570);
//             }
//             if(std::abs(t.refTrkPdgId[itrk])==13){
//                 v1.SetM(0.105658);
//             }
//             if(std::abs(t.refTrkPdgId[itrk])==11){
//                 v1.SetM(0.000510);
//             }
//             if(std::abs(t.refTrkPdgId[itrk])==2212){
//                 v1.SetM(0.938272);
//             }
//             if(std::abs(t.refTrkPdgId[itrk])==321){
//                 v1.SetM(0.493677);
//             }

//             // erase from list of tracks
//             for (Int_t i = 0; i < trackVectors.size(); i++){
//                 if (trackVectors[i]==v1){
//                     indices.push_back(i);
//                 }
//             }

//             // add to map of reconstructed b hadrons
//             if (mp.count(t.refTrkSta[itrk])){
//                 mp[t.refTrkSta[itrk]]+=v1;
//             }
//             else{
//                 mp.insert({t.refTrkSta[itrk], v1});
//             }
//         }
//     }

//     //add reconstructed b hadrons to list of tracks
//     std::sort(indices.begin(),indices.end(), greater<int>());
//     for (Int_t i = 0; i < indices.size(); i++){
//         trackVectors.erase(trackVectors.begin() + indices[i]);
//     }
    
//     for (auto i : mp) trackVectors.push_back(i.second);
// }

// SUBJET CLUSTERING METHOD - groups tracks until only two remain (which are then our subjets)
void CambridgeAachen(std::vector<ROOT::Math::PtEtaPhiMVector>& trackVectors, std::map<Float_t, std::vector<ROOT::Math::PtEtaPhiMVector>>& Merges){

    // While there are more than two tracks
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

// Method called to get substructure variable values using GEN information
vector<Float_t> getGenVals(Long64_t& ient, Int_t& ijet, tTree& t){

    t.GetEntry(ient); 

    // Pass eta and pt cuts using gen info
    if (std::abs(t.refeta[ijet]) > 2.) return {0.0,0.0,0.0};
    if ((t.refpt[ijet]) < 60 || (t.refpt[ijet]) > 150) return {0.0,0.0,0.0};

    //vector of TLorentz vectors for each track inside jet, used for subjet clustering
    std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;

    // Loop over tracks using track jet IDs --------------------------------------------------------------
    for (Int_t itrk = 0; itrk < t.nrefTrk; itrk++) {
        if (t.refTrkJetId[itrk] != ijet) continue; 

        //create vectors for each track
        ROOT::Math::PtEtaPhiMVector v1;
        v1.SetEta(t.refTrkEta[itrk]);
        v1.SetPt(t.refTrkPt[itrk]);
        v1.SetPhi(t.refTrkPhi[itrk]);

        if(std::abs(t.refTrkPdgId[itrk])==211){
            v1.SetM(0.139570);
        }
        if(std::abs(t.refTrkPdgId[itrk])==13){
            v1.SetM(0.105658);
        }
        if(std::abs(t.refTrkPdgId[itrk])==11){
            v1.SetM(0.000510);
        }
        if(std::abs(t.refTrkPdgId[itrk])==2212){
            v1.SetM(0.938272);
        }
        if(std::abs(t.refTrkPdgId[itrk])==321){
            v1.SetM(0.493677);
        }

        trackVectors.push_back(v1);
    }

    std::vector<ROOT::Math::PtEtaPhiMVector> initialVectors = trackVectors;

    //do not cluster into subjets if there are less than two tracks
    if (trackVectors.size()<2) return {0.0,0.0,0.0};

    // DECLUSTERING: create two subjets from tracks in jet ------------------------------------------------------------------------------
    std::map<Float_t, std::vector<ROOT::Math::PtEtaPhiMVector>> Merges;
    CambridgeAachen(trackVectors, Merges);

    //SOFTDROP: make sure new tracks pass softdrop, if not go back to previous merges, redo cambridgeAachen, and test again ----------------------------
    Int_t failed = 0;
    while (SoftDrop(trackVectors)==1){
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
        if (failed == 1) return {0.0,0.0,0.0};
        Merges.clear();
        initialVectors = trackVectors;
        CambridgeAachen(trackVectors, Merges);
    }

    Double_t minPt = std::numeric_limits<double>::infinity();
    for (Int_t i=0; i < 2; i++){
        if (trackVectors[i].Pt() < minPt){
            minPt = trackVectors[i].Pt();
        }
    }
    // KT CUT
    if (std::log(minPt*(ROOT::Math::VectorUtil::DeltaR(trackVectors[0], trackVectors[1])))<0) return {0.0,0.0,0.0};

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
    return {Zg_value,Kt_value,Rg_value};
}

vector<Float_t> getRecoVals(Long64_t& ient, Int_t& ijet, tTree& t){

    t.GetEntry(ient); 

    // Jet energy resolution 
    Float_t jtpt=t.jtpt[ijet];
    // if (t.jer_opt[ijet]=="nom") jtpt = jtpt * t.jer_sf_nom[ijet];
    jtpt = jtpt * t.jer_sf_down[ijet];
    // else if (t.jer_opt[ijet]=="down") jtpt = jtpt * t.jer_sf_down[ijet];

    // Reco eta and pt cuts
    if (std::abs(t.jteta[ijet]) > 2.) return {0.0,0.0,0.0};
    if ((jtpt) < 60 || (jtpt) > 150) return {0.0,0.0,0.0};

    //vector of TLorentz vectors for each track inside jet, used for subjet clustering
    std::vector<ROOT::Math::PtEtaPhiMVector> trackVectors;

    // Loop over the tracks of a specific jet using trkJetId --------------------------------------------------------------
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
    if (trackVectors.size()<2) return {0.0,0.0,0.0};

    // DECLUSTERING: create two subjets from tracks in jet ------------------------------------------------------------------------------
    std::map<Float_t, std::vector<ROOT::Math::PtEtaPhiMVector>> Merges;
    CambridgeAachen(trackVectors, Merges);

    //SOFTDROP: make sure new tracks pass softdrop, if not go back to previous merges, redo cambridgeAachen, and test again ----------------------------
    Int_t failed = 0;
    while (SoftDrop(trackVectors)==1){
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
        if (failed == 1) return {0.0,0.0,0.0};
        Merges.clear();
        initialVectors = trackVectors;
        CambridgeAachen(trackVectors, Merges);
    }

    Double_t minPt = std::numeric_limits<double>::infinity();
    for (Int_t i=0; i < 2; i++){
        if (trackVectors[i].Pt() < minPt){
            minPt = trackVectors[i].Pt();
        }
    }
    // KT CUT
    if (std::log(minPt*(ROOT::Math::VectorUtil::DeltaR(trackVectors[0], trackVectors[1])))<0) return {0.0,0.0,0.0};

    // // Final tracks lists each initial track inside of each subjet
    // map<Float_t, vector<ROOT::Math::PtEtaPhiMVector>> finalTracks;
    // for (Int_t k=0; k < 2; k++){
    //     for (auto i : Merges){
    //         if (abs(i.first-trackVectors[k].Pt())<0.01){
    //             vector<ROOT::Math::PtEtaPhiMVector> tempFinalTracks;
    //             for (auto j : i.second) tempFinalTracks.push_back(j);
    //             finalTracks.insert({i.first, tempFinalTracks});
    //         }
    //     }
    //     for (auto i : initialVectors) if (i==trackVectors[k]) finalTracks.insert({trackVectors[k].Pt(), {trackVectors[k]}});
    // }

    // Float_t maxPt = 0;
    // for (Int_t i=0; i < 2; i++){
    //     if (trackVectors[i].Pt() > maxPt){
    //         maxPt = trackVectors[i].Pt();
    //     }
    // }
    
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
    return {Zg_value,Kt_value,Rg_value};
}

// MAIN METHOD --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void run_MC_inclusive_response(){

    TFile *_file0 = TFile::Open("/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root"); 
    TString fin_name = "/data_CMS/cms/kalipoliti/qcdMC/bjet/aggrTMVA_fixedMassBug/merged_HiForestMiniAOD.root";
    TString fout_name = "histogramsSvtxGenBjetUnfoldingInclusiveMC_JERdown.root";
    tTree t(fin_name);

    // Declare histograms
    TH2D *h_Zg_inclusive = new TH2D("h_Zg_inclusive", "x=Zg inclusive", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_Kt_inclusive = new TH2D("h_Kt_inclusive", "x=Kt inclusive", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
    TH2D *h_Rg_inclusive = new TH2D("h_Rg_inclusive", "x=Rg inclusive", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);

    // Turn off all branches and turn on only the interesting ones
    // Attention! If a branch is off, it will return bs without crashing 
    t.SetBranchStatus("*", 0);
    std::vector<TString> active_branches = {
        "weight", "evt",
        "nref", "jtpt", "jteta", "lumi", "jtptCh", "discr_particleNet_BvsAll", "pthat", "ngen", "jtphi", "jtm", // reco jet
        // "sjt1Eta","sjt2Eta","sjt1Pt","sjt2Pt", "rsjt1Eta","rsjt2Eta","rsjt1Pt","rsjt2Pt", "rsjt1Y","rsjt2Y","rsjt1Phi","rsjt2Phi",
        "refpt", "refeta", "jtHadFlav", "jtNbHad", //"sjt1Eta", "sjt2Eta", "sjt1Pt", "sjt2Pt",  // gen jet
        "ntrk", "trkPt", "trkMatchSta", "trkJetId", "trkSvtxId", "trkPdgId", "trkEta", "trkPhi", // tracks
        "jtNsvtx", "svtxpt", "svtxJetId", //"svtxm","svtxmcorr","svtxdls", "svtxdls2d", "svtxNtrk", "svtxnormchi2", "svtxpt", // SVs 
        "nrefTrk", "refTrkJetId", "refTrkSta", "refTrkEta", "refTrkPt", "refTrkPhi", "refTrkPdgId",
        "jtBpt","refBpt", "refptCh", 
        "jer_opt", "jer_sf_up", "jer_sf_down",
        // "HLT_HIAK4PFJet80_v1", "HLT_HIAK4PFJet100_v1", "HLT_HIAK4PFJet60_v1", "HLT_HIAK4PFJet40_v1",
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

            // Jet energy resolution
            Float_t jtpt=t.jtpt[ijet];
            // if (t.jer_opt[ijet]=="nom") jtpt = jtpt * t.jer_sf_nom[ijet];
            // jtpt = jtpt * t.jer_sf_down[ijet];
            // else if (t.jer_opt[ijet]=="down") jtpt = jtpt * t.jer_sf_down[ijet];

            // Get values of substructure variables using gen and reco info
            vector<Float_t> genVals = getGenVals(ient, ijet,t);
            vector<Float_t> recoVals = getRecoVals(ient, ijet,t);

            create_response_fast(jtpt,t.jtptCh[ijet],recoVals[2],recoVals[1],recoVals[0],genVals[2],genVals[1],genVals[0],t.weight,t.pthat, t.refptCh[ijet], t.jtBpt[ijet], t.refBpt[ijet], t.refpt[ijet]);

            vector<Float_t> emptyCheck = {0.0,0.0,0.0};
            if (recoVals!=emptyCheck){
                h_Zg_inclusive->Fill(recoVals[0], jtpt, t.weight);
                h_Kt_inclusive->Fill(recoVals[1],jtpt, t.weight);
                h_Rg_inclusive->Fill(recoVals[2],jtpt, t.weight);
            }

        }
    } // entry loop ---------------------------------------------------------------------------------------------------------------------------

    // Save histograms in a new file 
    std::cout << "\n(Re)creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    for (auto h : {
        h_Zg_inclusive,h_Rg_inclusive,h_Kt_inclusive,
    }) {
        // h->Scale(1.0/h->Integral());
        h->Write();
    }

    fout->Close();
    
    // Create purity and efficiency histograms for corrections
    purityEfficiencyHists();
}
