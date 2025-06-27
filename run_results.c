#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldInvert.h"

void run_results(TString observable="zg", TString jer_opt="nom", TString jec_opt="nom")
{    
    // ---------------- Options for b    -------------------
    bool withSF = true;
    bool sfDown = false;
    bool sfUp = false;

    TString sample_sf = "pythia_PF40"; 
    TString label_sf = "aggrTMVA_XXT";

    // ---------------- Options for both -------------------
    TString suffix_in = "_jer_" + jer_opt + "_jec_" + jec_opt;

    TString sample_data = "data_PF40to100";
    TString label_data = "aggrTMVA_XXT";

    bool unfoldBayes = false;
    bool purityPythia = true;
    bool responsePythia = true;
    bool efficiencyPythia = true;
    bool btagEfficiencyPythia = false;

    TString pythia_unfolding = "pythia_PF40";
    TString herwig_unfolding = "herwig_official_PF40";

    TString sample_truth = "herwig_official_PF40"; // sample used for comparing to truth
    TString label_truth = "aggrTMVA_XXT";
    TString sample_unfolding = "herwig_official_PF40"; // sample used for unfolding
    TString label_unfolding = "aggrTMVA_XXT";

    bool inclusive = sample_unfolding.Contains("dijet");
    if (inclusive) withSF = false;

    TString suffix_out = suffix_in;
    if (withSF) suffix_out += "_withSF";
    if (sfUp) suffix_out += "Up";
    if (sfDown) suffix_out += "Down";

    int ibin_pt = 2;

    std::cout << "Options:"
              << "\n\tunfoldBayes:" << unfoldBayes
              << "\n\tibin_pt:" << ibin_pt
              << "\n\tinclusive:" << inclusive
              << "\n\t\twithSF:" << withSF
              << "\n\t\tsfUp:" << sfUp
              << "\n\t\tsfDown:" << sfDown
              << std::endl;

    TString fout_name = "unfoldedHistograms.root";
    // --------------------------------------------


    // ---------------- Plotting setup ------------
    gSystem->Load("libRooUnfold.so");
    gStyle->SetErrorX(0.5);

    TString xlabel;
    if (observable=="rg") xlabel = "ln(R/R_{g})";
    else if (observable=="zg") xlabel = "z_{g}";
    else if (observable=="zpt") xlabel = "ln(k_{T})";
    TString ylabel = "1/N dN/d" + xlabel;

    // Float_t text_size = 200.;
    // gStyle->SetTextSize(text_size);
    // gStyle->SetLegendTextSize(text_size);
    // gStyle->SetLabelSize(text_size, "XYZ");
    // gStyle->SetTitleSize(text_size, "XYZ");
    // --------------------------------------------


    // ----------- Grab data -----------

    TString fname_data = "histogramsDataFittedNewBinsWithLow"+observable+"_new.root";
    // TString fname_data = "histogramsDataFittedNewBins"+observable+".root";
    std::cout << "Getting data from " << fname_data << std::endl;
    TFile *fin_data = new TFile(fname_data);
    TH2D *h_data = (TH2D *) fin_data->Get("h_data_2d");
    TH2D *h_sig_fraction = (TH2D *) fin_data->Get("signalFractions");

    TString fname_data_inclusive = "histogramsSvtxGenDataHighInclusive_new.root";
    TFile *fin_data_inclusive = new TFile(fname_data_inclusive);
    TH2D *h_data_inclusive = (TH2D *) fin_data_inclusive->Get("h_" + observable + "_inclusive");
    TString fname_data_inclusive2 = "histogramsSvtxGenDataLowInclusive_new.root";
    TFile *fin_data_inclusive2 = new TFile(fname_data_inclusive2);
    TH2D *h_data_inclusive2 = (TH2D *) fin_data_inclusive2->Get("h_" + observable + "_inclusive");
    h_data_inclusive->Add(h_data_inclusive2);
    // h_data->GetYaxis()->SetRange(0, h_data->GetNbinsY() + 1); // all mb bins

    // TH2D *h_sig_fraction;
    // TString fname_data = "responseMatrixBjetNewBins.root";
    // TFile *fin_data = new TFile(fname_data);
    // TH2D *h_data = (TH2D *) fin_data->Get("h_half0_purity_numerator_"+observable+"pt"); // true MC to compare w/ data after BOTH efficiency corrections

    Int_t isData=0;
    if (fname_data.Contains("Data")) isData=1;

    // TFile *file4 = TFile::Open("histogramsSvtxGenBjetUnfoldingBinned3.root");
    // TH3D *hist2b_Y_B_3D = (TH3D*)file4->Get("h_InvMass_vs_Zg_2b");
    // TH2D *h_data = (TH2D*) hist2b_Y_B_3D->Project3D("zy");

    // TFile *file4 = TFile::Open("histogramsSvtxGenBjetUnfoldingBinned.root");
    // TH3D *hist2b_Y_B_3D = (TH3D*)file4->Get("h_InvMass_vs_Zg_2b");
    // TH2D *h_data = (TH2D *) hist2b_Y_B_3D->Project3D("zy");

    int nbins_x = h_data->GetNbinsX();
    int nbins_pt = h_data->GetNbinsY();
    int dim = nbins_pt*nbins_x;

    int ibin_x_min = 1;
    int ibin_x_max = nbins_x;

    double pt_min = h_data->GetYaxis()->GetBinLowEdge(ibin_pt);
    double pt_max = h_data->GetYaxis()->GetBinUpEdge(ibin_pt);

    double x_min = h_data->GetXaxis()->GetBinLowEdge(1);
    double x_max = h_data->GetXaxis()->GetBinUpEdge(nbins_x);
    
    // TH2D *h_data_reco = (TH2D *) h_data->Project3D("zx");
    // h_data_reco->GetXaxis()->SetRange(1, h_data_reco->GetNbinsX());
    // h_data_reco->GetYaxis()->SetRange(1, h_data_reco->GetNbinsY());

    // Multiply histograms by signal fraction
    TH2D *h_data_after_fit = (TH2D *) h_data->Clone("h_data_after_fit");
    if (isData==1) h_data_after_fit->Multiply(h_sig_fraction);
    
    // if (!inclusive) {
    //     std::cout << "\t---->Multiplying by signal fraction" << std::endl;
    //     // Grab signal fraction from template fit
    //     TString fit_option = "_glued";
    //     TString fname_fit = "../template_fit/histos/fitted_parameters_RooFit_data_"+label_data+"_" + observable + fit_option + "_jer_nom_jec_nom.root";
    //     std::cout << "Getting signal fraction from " << fname_fit << std::endl;
    //     TFile *fin_fit = new TFile(fname_fit);
    //     TH2D *h_sig_fraction = (TH2D *) fin_fit->Get("h_sig_fraction");    
    //     h_data_after_fit->Multiply(h_sig_fraction);

    //     std::cout << "\t---->Dividing by SF" << std::endl;
    //     TString fin_sf_name = "../btag/histos/aggrTMVA_inclusive_"+observable+"_sfs.root";
    //     std::cout << "Getting SFs from " << fin_sf_name << std::endl;
    //     TFile *fin_sf = new TFile(fin_sf_name);
    //     TH2D *h_sf = (TH2D*) fin_sf->Get("h_eff_sf_3bins");

    //     if (sfUp||sfDown) {
    //         TFile *fin_sf_unc = new TFile("../btag/histos/aggrTMVA_inclusive_"+observable+"_sf_unc.root");
    //         if (sfUp) {
    //             TH2D *h_sf_unc = (TH2D *) fin_sf_unc->Get("h_sf_unc_up_3bins");
    //             h_sf->Add(h_sf_unc);
    //         } else if (sfDown) {
    //             TH2D *h_sf_unc = (TH2D *) fin_sf_unc->Get("h_sf_unc_down_3bins");
    //             h_sf->Add(h_sf_unc, -1.);
    //         }
    //     }  
    //     h_data_after_fit->Divide(h_sf);
    // }
    
    // Note: Result = unfold(raw * purity) * 1 / (efficiency)
    //       fakes are negligible

    // ---- Grab response matrix + corrections
    TString fname_pythia_unfolding = "responseMatrixBjet_new.root";
    // TString fname_pythia_unfolding2 = "responseMatrixDijetNewBins.root";
    std::cout << "Getting pythia response + corrections from : " << fname_pythia_unfolding << std::endl;
    TFile *fin_pythia_unfolding = new TFile(fname_pythia_unfolding);
    // TFile *fin_pythia_unfolding2 = new TFile(fname_pythia_unfolding2);
    TString fname_pythia_unfolding_inclusive = "responseMatrixDijetInclusiveMC_new.root";
    TFile *fin_pythia_unfolding_inclusive = new TFile(fname_pythia_unfolding_inclusive);

    // TString fname_herwig_unfolding = "./histos/"+herwig_unfolding+"_"+label_unfolding+"_response" + suffix_in + ".root";
    // std::cout << "Getting herwig response + corrections from : " << fname_herwig_unfolding << std::endl;
    // TFile *fin_herwig_unfolding = new TFile(fname_herwig_unfolding);

    TH2D *h_full_purity;
    TH2D *h_full_efficiency;
    TH2D *h_mc_reco;
    RooUnfoldResponse *response;
    TH2D *h_mc_true_no_eff;

    TH2D *h_full_purity_inclusive;
    TH2D *h_full_efficiency_inclusive;
    RooUnfoldResponse *response_inclusive;

    // get purity 
    if (purityPythia)
        h_full_purity = (TH2D *) fin_pythia_unfolding->Get("h_full_purity_"+observable+"pt"); // reconstruction purity correction
        h_full_purity_inclusive = (TH2D *) fin_pythia_unfolding_inclusive->Get("h_full_purity_"+observable+"pt");
    // else 
    //     h_full_purity = (TH2D *) fin_herwig_unfolding->Get("h_full_purity_"+observable+"pt"); // reconstruction purity correction

    // get efficiency
    if (efficiencyPythia)
        h_full_efficiency = (TH2D *) fin_pythia_unfolding->Get("h_full_efficiency_"+observable+"pt"); // reconsturction efficiency correction
        h_full_efficiency_inclusive = (TH2D *) fin_pythia_unfolding_inclusive->Get("h_full_efficiency_"+observable+"pt");
    // else 
    //     h_full_efficiency = (TH2D *) fin_herwig_unfolding->Get("h_full_efficiency_"+observable+"pt"); // reconsturction efficiency correction

    // get transfer matrix
    if (responsePythia) {
        h_mc_reco = (TH2D *) fin_pythia_unfolding->Get("h_full_purity_numerator_"+observable+"pt"); // reco MC to compare w/ data after purity correction
        response = (RooUnfoldResponse *) fin_pythia_unfolding->Get("response_full_"+observable+"pt"); // response 
        response_inclusive = (RooUnfoldResponse *) fin_pythia_unfolding_inclusive->Get("response_full_"+observable+"pt");
        // response2 = (RooUnfoldResponse *) fin_pythia_unfolding2->Get("response_full_"+observable+"pt");
        // response->Add(response2);
        h_mc_true_no_eff = (TH2D *) fin_pythia_unfolding->Get("h_full_efficiency_numerator_"+observable+"pt"); // true MC to compare w/ data BEFORE efficiency corrections

    } 
    // else {
    //     h_mc_reco = (TH2D *) fin_herwig_unfolding->Get("h_full_purity_numerator_"+observable+"pt"); // reco MC to compare w/ data after purity correction
    //     response = (RooUnfoldResponse *) fin_herwig_unfolding->Get("response_full_"+observable+"pt"); // response 
    //     h_mc_true_no_eff = (TH2D *) fin_herwig_unfolding->Get("h_full_efficiency_numerator_"+observable+"pt"); // true MC to compare w/ data BEFORE efficiency corrections
    // }
    // TEST FSR DOWN
    // TFile *fin_pythia_var = new TFile("./histos/dijet_PF40_aggrTMVA_inclusive_response_jer_nom_jec_nom_pythia_variations.root");
    // h_full_purity = (TH2D *) fin_pythia_var->Get("h_FSRdown_purity_"+observable+"pt");

    // ---- Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    std::cout << "\t---->Condition number nominal = " << cond_number
              << std::endl;

    // ---- Grab the truth level MC ---- 
    TString fname_response_truth = "responseMatrixBjet_new.root";
    TString fname_response_truth2 = "responseMatrixDijet_new.root";
    std::cout << "Getting truth from : " << fname_response_truth << std::endl;
    TFile *fin_response_truth = new TFile(fname_response_truth);
    TFile *fin_response_truth2 = new TFile(fname_response_truth2);
    TH2D *h_mc_true = (TH2D *) fin_response_truth->Get("h_full_efficiency_denominator_"+observable+"pt"); // true MC to compare w/ data after BOTH efficiency corrections
    // TH2D *h_mc_true2 = (TH2D *) fin_response_truth2->Get("h_full_efficiency_denominator_"+observable+"pt"); // true MC to compare w/ data after BOTH efficiency corrections
    // h_mc_true->Add(h_mc_true2);

    //------- Apply purity correction
    std::cout << "\t---->Multiplying data by purity" << std::endl;
    TH2D *h_data_purity_corrected = (TH2D *) h_data_after_fit->Clone("h_data_purity_corrected");
    if (isData==1) h_data_purity_corrected->Multiply(h_full_purity);
    TH2D *h_data_purity_corrected_inclusive= (TH2D *) h_data_inclusive->Clone("h_data_purity_corrected_inclusive");

    // ---- Unfold

    std::cout << "\t---->Unfolding" << std::endl;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    TH2D *h_data_unfolded;
    TMatrixD covariance_matrix_before_unfolding(dim,dim);
    TMatrixD covariance_matrix_after_unfolding(dim,dim);
    if (unfoldBayes) {
        Int_t niter =  100;
        RooUnfoldBayes unfold(response, h_data_purity_corrected, niter);
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment);
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    } else {
        RooUnfoldInvert unfold(response, h_data_purity_corrected);
        h_data_unfolded = (TH2D *) unfold.Hreco(errorTreatment);
        covariance_matrix_before_unfolding = unfold.GetMeasuredCov();
        covariance_matrix_after_unfolding = unfold.Ereco();
    }    
    RooUnfoldInvert unfold(response_inclusive, h_data_purity_corrected_inclusive);
    TH2D *h_data_unfolded_inclusive = (TH2D *) unfold.Hreco(errorTreatment);

    h_data_unfolded->SetName("h_data_unfolded");

    // TH2D *h_data_unfolded = (TH2D *) h_data_purity_corrected->Clone("h_data_unfolded");

    // ---- Fold back
    std::cout << "\t---->Refolding" << std::endl;
    TH2D *h_data_refolded = (TH2D *) response->ApplyToTruth(h_data_unfolded, "h_data_refolded");

    // ---- Apply efficiency correction 
    std::cout << "\t---->Dividing by recostruction efficiency" << std::endl;
    TH2D *h_data_efficiency_corrected = (TH2D *) h_data_unfolded->Clone("h_data_efficiency_corrected");
    h_data_efficiency_corrected->Divide(h_full_efficiency);
    TH2D *h_data_efficiency_corrected_inclusive = (TH2D *) h_data_unfolded_inclusive->Clone("h_data_efficiency_corrected_inclusive");
    h_data_efficiency_corrected_inclusive->Divide(h_full_efficiency_inclusive);

    // ---- Final corrections
    TH2D *h_data_fully_corrected = (TH2D *) h_data_efficiency_corrected->Clone("h_data_fully_corrected");
    TH2D *h_data_fully_corrected_inclusive = (TH2D *) h_data_efficiency_corrected_inclusive->Clone("h_data_fully_corrected_inclusive");
    TH2D *h_eff;
    TH2D *h_eff_withSF;
    if (!inclusive&&false) {
        // ---- Grab b tagging efficiency correction 
        TString fname_b_tag_eff;
        if (btagEfficiencyPythia) 
            fname_b_tag_eff = "../btag/histos/"+pythia_unfolding+"_"+label_sf+"_" +observable + "_efficiency" + suffix_out + ".root";
        else 
            fname_b_tag_eff = "../btag/histos/"+herwig_unfolding+"_"+label_sf+"_" +observable + "_efficiency" + suffix_out + ".root";
        std::cout << "Getting b tagging efficiency from: " << fname_b_tag_eff << std::endl;
        TFile *fin_b_tag_eff = new TFile(fname_b_tag_eff);
        h_eff = (TH2D *) fin_b_tag_eff->Get("h_eff");
        h_eff_withSF = (TH2D *) fin_b_tag_eff->Get("h_eff_withSF");
        std::cout << "\t---->Dividing by b tagging efficiency (with SF)" << std::endl;
        // if (withSF) h_data_fully_corrected->Divide(h_eff_withSF);
        // else h_data_fully_corrected->Divide(h_eff);
        h_data_fully_corrected->Divide(h_eff);

        // // ----- SF refolding test 
        // if (withSF&&false) {
        //     std::cout << "\t---->SF refolding test" << std::endl;
        //     // undo the efficiency correction but leave the SF
        //     TH2D *h_data_unfolded_with_sf = (TH2D *) h_data_fully_corrected->Clone("h_data_unfolded_with_sf");
        //     h_data_unfolded_with_sf->Multiply(h_eff);
        //     h_data_unfolded_with_sf->Multiply(h_full_efficiency);
        //     // and refold
        //     TH2D *h_data_refolded_with_sf = (TH2D *) response->ApplyToTruth(h_data_unfolded_with_sf, "h_data_refolded_with_sf");
        //     // fix the unc 
        //     for (int ibin_x=1; ibin_x<=nbins_x; ibin_x++) {
        //         for (int ibin_pt=1; ibin_pt<=nbins_pt; ibin_pt++) {
        //             double stat_unc = h_data_purity_corrected->GetBinError(ibin_x, ibin_pt);
        //             h_data_refolded->SetBinError(ibin_x, ibin_pt, stat_unc);
        //             h_data_refolded_with_sf->SetBinError(ibin_x, ibin_pt, stat_unc);
        //         }
        //     }  
        //     // compare to other refolded => ratio should be the reco SF
        //     TH2D *h_refolded_ratio = (TH2D *) h_data_refolded->Clone("h_refolded_ratio");
        //     h_refolded_ratio->Divide(h_data_refolded_with_sf);

        //     TH1D *h_refolded_ratio_1d = (TH1D *) h_refolded_ratio->ProjectionX("h_refolded_ratio_1d", 2, 2);
        //     h_refolded_ratio_1d->SetTitle("SF from refolding");

        //     TString fname_sfs = "../btag/histos/variations/aggrTMVA_inclusive_"+observable+"_sfs_nom.root";
        //     std::cout << "Getting reco sfs from: " << fname_sfs << std::endl;
        //     TFile *fin_sfs = new TFile(fname_sfs);
        //     TH2D *h_sf = (TH2D *) fin_sfs->Get("h_eff_sf")->Clone("h_sf");
        //     TH1D *h_sf_1d = (TH1D *) h_sf->ProjectionX("h_sf_1d", 1, 1);
        //     h_sf_1d->SetMarkerColor(kRed);
        //     h_sf_1d->SetLineColor(kRed);
        //     h_sf_1d->SetTitle("SF reco");
        //     h_sf_1d->GetXaxis()->SetTitle(xlabel);
        //     h_sf_1d->GetYaxis()->SetTitle("b tagging efficiency SF");
        //     h_sf_1d->SetMinimum(0.95);
        //     h_sf_1d->SetMaximum(1.3);

        //     TCanvas *c_sf_closure = new TCanvas("c_sf_closure", "", 800, 600);
        //     TPad *top_pad = new TPad("top_pad", "", 0., 0.33, 1., 1.);
        //     TPad *bot_pad = new TPad("bot_pad", "", 0., 0., 1., 0.33);
        //     bot_pad->SetTopMargin(0.01);
        //     bot_pad->SetBottomMargin(0.3);
        //     top_pad->SetBottomMargin(0.01);

        //     top_pad->cd();
        //     h_sf_1d->Draw("pe1 same");
        //     h_refolded_ratio_1d->Draw("pe1 same");
        //     auto leg_sf = top_pad->BuildLegend();
        //     leg_sf->SetHeader("100 < p_{T}^{jet} < 120 (GeV)");
        //     drawHeader();

        //     bot_pad->cd();
        //     TH1D *h_sf_ratio_1d = (TH1D *) h_refolded_ratio_1d->Clone("h_sf_ratio_1d");
        //     h_sf_ratio_1d->Divide(h_sf_1d);
        //     h_sf_ratio_1d->GetYaxis()->SetTitle("SF refolded / SF reco");
        //     h_sf_ratio_1d->GetYaxis()->SetNdivisions(6);
        //     h_sf_ratio_1d->GetXaxis()->SetTitle(xlabel);
        //     h_sf_ratio_1d->GetXaxis()->SetTitleOffset(4.);
        //     h_sf_ratio_1d->Draw("pe1");

        //     TLine *line = new TLine(x_min, 1., x_max, 1.);
        //     line->SetLineWidth(2.); 
        //     line->SetLineStyle(kDashed);
        //     line->SetLineColor(kGray);
        //     line->Draw();

        //     c_sf_closure->cd();
        //     top_pad->Draw();
        //     bot_pad->Draw();

        //     c_sf_closure->Draw();
        //     // c_sf_closure->Print("plots_an/"+sample_unfolding+"_SF_closure_test_"+observable+".png");
        // }
    } // end if !inclusive

    // ---- Quantitative bottomline test
    if (true) {
        // for (int col=0; col<dim; col++) {
        //     int row1=2;
        //     int row2=3;
        //     covariance_matrix_before_unfolding[row1][col]=covariance_matrix_before_unfolding[row2][col];
        // }
        TMatrixD covariance_matrix_before_unfolding_inverted = covariance_matrix_before_unfolding;
        covariance_matrix_before_unfolding_inverted.Invert();
        TMatrixD covariance_matrix_after_unfolding_inverted = covariance_matrix_after_unfolding;
        covariance_matrix_after_unfolding_inverted.Invert();
        // ---- sanity checks ----
        // TMatrixD cc1(dim, dim);
        // cc1.Mult(covariance_matrix_before_unfolding, covariance_matrix_before_unfolding_inverted);
        // cc1.Mult(covariance_matrix_before_unfolding_inverted, covariance_matrix_before_unfolding);
        // cc1.Mult(covariance_matrix_after_unfolding, covariance_matrix_after_unfolding_inverted);
        // cc1.Mult(covariance_matrix_after_unfolding_inverted, covariance_matrix_after_unfolding);
        // cc1.Draw("colz");
        // ------------------------
        TVectorD ydiff_before_unfolding(dim); // data-MC after purity correction
        TVectorD ydiff_after_unfolding(dim); // data-MC after unfolding before efficiency correction
        TVectorD ystat_before_unfolding(dim); // data stat unc
        for (int i=1; i<=nbins_x; i++) {
            for (int j=1; j<=nbins_pt; j++) {
                // x:i=1,...,n ; y:j=1,...,m ; from matrix(i,j) to vector(k) => k=(i-1)+[(j-1)*n]
                ydiff_before_unfolding[i-1+nbins_x*(j-1)]=(h_data_purity_corrected->GetBinContent(i,j)-h_mc_reco->GetBinContent(i,j)); 
                ydiff_after_unfolding[i-1+nbins_x*(j-1)]=(h_data_unfolded->GetBinContent(i,j)-h_mc_true_no_eff->GetBinContent(i,j)); 
                ystat_before_unfolding[i-1+nbins_x*(j-1)]=h_data_purity_corrected->GetBinError(i,j); 
            }
        }

        TVectorD cov_times_ydiff_before_unfolding(dim); // [yn/sn^2 ... y1/s1^2] 
        TVectorD cov_times_ydiff_after_unfolding(dim); // [yn/sn^2 ... y1/s1^2]
        for ( int i=dim-1; i>=0; i--) {
            TMatrixDRow row_before_unfolding(covariance_matrix_before_unfolding_inverted,i);
            TMatrixDRow row_after_unfolding(covariance_matrix_after_unfolding_inverted,i);
            double sum_before_unfolding=0;
            double sum_after_unfolding=0;
            for ( int j=0; j<dim; j++) {
                sum_before_unfolding += row_before_unfolding[j]*ydiff_before_unfolding[j];
                sum_after_unfolding += row_after_unfolding[j]*ydiff_after_unfolding[j];
            }
            cov_times_ydiff_before_unfolding[dim-i-1] = sum_before_unfolding;
            cov_times_ydiff_after_unfolding[dim-i-1] = sum_after_unfolding;
        }

        double chi2_before_unfolding = 0.;
        double naive_chi2_before_unfolding = 0.;
        double chi2_after_unfolding = 0.;
        double test_chi2 = 0.;
        for ( int i=0; i<dim; i++) {
            chi2_before_unfolding += ydiff_before_unfolding[i]*cov_times_ydiff_before_unfolding[dim-i-1];
            chi2_after_unfolding += ydiff_after_unfolding[i]*cov_times_ydiff_after_unfolding[dim-i-1];
            naive_chi2_before_unfolding += std::pow(ydiff_before_unfolding[i], 2) / std::pow(ystat_before_unfolding[i], 2);

            // std::cout << ydiff_before_unfolding[i]*cov_times_ydiff_before_unfolding[dim-i-1] 
            //           << ", " 
            //           << ydiff_after_unfolding[i]*cov_times_ydiff_after_unfolding[dim-i-1] << std::endl;
        }
        // std::cout << "ndof = " << 
        std::cout << Form("naive_chi2_before_unfolding = %f", naive_chi2_before_unfolding) << std::endl;
        std::cout << Form("chi2_before_unfolding = %f", chi2_before_unfolding) << std::endl;
        std::cout << Form("chi2_after_unfolding = %f", chi2_after_unfolding) << std::endl;
        std::cout << Form("chi2_after_unfolding / chi2_before_unfolding = %f", chi2_after_unfolding / chi2_before_unfolding) << std::endl;
    }

    // ---- Graphical bottomline test
    if (true) {
        std::cout << "Performing graphical bottomline test" << std::endl;
        TH1D *h_mc_reco_1d = (TH1D *) h_mc_reco->ProjectionX("h_mc_reco_1d", ibin_pt, ibin_pt);
        TH1D *h_mc_true_1d = (TH1D *) h_mc_true->ProjectionX("h_mc_true_1d", ibin_pt, ibin_pt);
        TH1D *h_data_purity_corrected_1d = (TH1D *) h_data_purity_corrected->ProjectionX("h_data_purity_corrected_1d", ibin_pt, ibin_pt);
        TH1D *h_data_fully_corrected_1d = (TH1D *) h_data_fully_corrected->ProjectionX("h_data_fully_corrected_1d", ibin_pt, ibin_pt);
        TH1D *h_data_refolded_1d = (TH1D *) h_data_refolded->ProjectionX("h_data_refolded_1d", ibin_pt, ibin_pt);

        TFile *file10 = TFile::Open("histogramsSvtxGenDijetUnfoldingInclusiveMC_new.root");
        TH2D *h_X_inclusive2d = (TH2D*)file10->Get("h_" + observable + "_inclusive");
        TH1D *h_X_inclusive = (TH1D *) h_X_inclusive2d->ProjectionX("h_X_inclusive", ibin_pt, ibin_pt);
        TFile *file11 = TFile::Open("histogramsSvtxGenBjetInclusiveDoubleB.root");
        TH1D *h_X_2B = (TH1D*)file11->Get("h_" + observable + "_2b");
        TH1D *h_data_fully_corrected_inclusive_1d = (TH1D *) h_data_fully_corrected_inclusive->ProjectionX("h_data_fully_corrected_inclusive_1d", ibin_pt, ibin_pt);

        TFile *fileA = TFile::Open("/grid_mnt/data__data.polcms/cms/kalipoliti/forJoseph/herwig_official_PF40_aggrTMVA_XXT_unfolded_histograms_zg_jer_nom_jec_nom_forJoseph.root");
        TH2D *data_singleB = (TH2D*)fileA->Get("h_data_unfolded");
        // data_singleB->GetXaxis()->SetRangeUser(0,2.5);
        TH1D *data_singleB_1d = (TH1D *) data_singleB->ProjectionX("data_singleB_1d", 2, 2);

        double ymax = 0.;
        for (auto h : {
                    h_mc_reco_1d, 
                    h_mc_true_1d,
                    h_data_purity_corrected_1d,
                    h_data_fully_corrected_1d,
                    h_data_refolded_1d,
                    h_X_inclusive,
                    h_X_2B,
                    h_data_fully_corrected_inclusive_1d,
                    data_singleB_1d,
                    }) {
                        h->GetXaxis()->SetRange(ibin_x_min, ibin_x_max);
                        h->Scale(1/h->Integral(), "width");
                        if (h->GetMaximum()>ymax) ymax =  h->GetMaximum();
                    }

        gStyle->SetOptTitle(0);

        TLegend *leg = new TLegend(0.5,0.6,0.85,0.8);
        if (observable=="zpt") leg = new TLegend(0.5,0.1,0.85,0.3);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetMargin(0.15);
        leg->SetHeader(Form("%.0f < p_{T}^{jet} < %.0f (GeV)", pt_min, pt_max));

        h_data_purity_corrected_1d->SetMarkerColor(kBlack);
        h_data_purity_corrected_1d->SetMarkerStyle(kFullCircle);
        h_data_purity_corrected_1d->SetMarkerSize(1);
        // h_data_fully_corrected_1d->GetYaxis()->SetRangeUser(0., ymax*1.1);
        // h_data_fully_corrected_1d->GetYaxis()->SetRangeUser(0., 0.8);
        h_data_purity_corrected_1d->GetYaxis()->SetTitle(ylabel);
        // leg->AddEntry(h_data_purity_corrected_1d, "Detector level data (tagged double b jets)", "pe1");

        h_mc_reco_1d->SetMarkerColor(kRed);
        h_mc_reco_1d->SetLineColor(kRed);
        h_mc_reco_1d->SetMarkerStyle(kFullTriangleUp);
        h_mc_reco_1d->SetMarkerSize(1);
        // leg->AddEntry(h_mc_reco_1d, "Detector level MC (tagged double b jets)", "pe1");

        h_data_fully_corrected_1d->SetMarkerColor(kBlue);
        h_data_fully_corrected_1d->SetLineColor(kBlue);
        h_data_fully_corrected_1d->SetMarkerStyle(kOpenCross);
        h_data_fully_corrected_1d->SetMarkerSize(1);
        leg->AddEntry(h_data_fully_corrected_1d, "Unfolded data (double b jets)", "pe1");

        h_mc_true_1d->SetMarkerColor(kRed);
        h_mc_true_1d->SetLineColor(kRed);
        h_mc_true_1d->SetMarkerStyle(kFullTriangleUp);
        h_mc_true_1d->SetMarkerSize(1);
        leg->AddEntry(h_mc_true_1d, "Particle level MC (double b jets)", "pe1");

        h_data_refolded_1d->SetMarkerColor(kBlue);
        h_data_refolded_1d->SetLineColor(kBlue);
        h_data_refolded_1d->SetMarkerStyle(kFullCross);
        h_data_refolded_1d->SetMarkerSize(1);
        // leg->AddEntry(h_data_refolded_1d, "Refolded data", "pe1");

        TCanvas *c_unfold = new TCanvas("c_unfold", "", 800, 600);
        TPad *pad1 = new TPad("pad1", "", 0., 0., 1., 0.3);
        TPad *pad2 = new TPad("pad2", "", 0., 0.3, 1., 1.);
        pad1->SetTopMargin(0.01);
        pad1->SetBottomMargin(0.3);
        pad2->SetBottomMargin(0.01);

        pad2->cd();

        TString drawOption="pe1 same";
        THStack *hs=new THStack();
        h_data_fully_corrected_1d->SetFillColor(kGray);
        h_data_fully_corrected_1d->SetMarkerStyle(kFullCircle);
        h_data_fully_corrected_1d->SetMarkerColor(kBlack);

        h_data_fully_corrected_inclusive_1d->SetMarkerStyle(kOpenCircle);
        h_data_fully_corrected_inclusive_1d->SetMarkerColor(kBlue);
        h_data_fully_corrected_inclusive_1d->SetFillColor(kGray);
        // h_data_fully_corrected_inclusive_1d->SetFillStyle(3001);

        // data_singleB_1d->SetFillColor(kRed);
        data_singleB_1d->SetMarkerStyle(kOpenTriangleUp);
        data_singleB_1d->SetMarkerColor(kRed);
        data_singleB_1d->SetLineColor(kRed);
        // data_singleB_1d->SetFillStyle(3001);

        leg->AddEntry(h_data_fully_corrected_inclusive_1d, "Unfolded data (inclusive jets)", "pe1");
        // leg->AddEntry(data_singleB_1d, "Unfolded data (single B jets)", "pe1");
        // h_data_purity_corrected_1d->SetLineColor(kBlack);
        // h_data_purity_corrected_1d->SetFillColor(kBlack);
        // h_data_purity_corrected_1d->SetFillStyle(3001);
        h_data_fully_corrected_1d->GetYaxis()->SetTitle(ylabel);
        hs->Add(h_data_fully_corrected_1d,"e");
        hs->Add(h_data_fully_corrected_1d,"e2");
        hs->Add(h_data_fully_corrected_inclusive_1d,"e");
        hs->Add(h_data_fully_corrected_inclusive_1d,"e2");
        data_singleB_1d->SetFillColor(0);
        // hs->Add(data_singleB_1d,"hist e");
        // hs->Add(data_singleB_1d,"e2");
        // hs->Add(h_data_purity_corrected_1d, "e");
        // hs->Add(h_data_purity_corrected_1d, "e2");
        // hs->Add(h_data_purity_corrected_1d, "hist");
        // hs->Add(h_mc_reco_1d, drawOption);
        // hs->Add(h_mc_reco_1d, "hist");
        // h_mc_true_1d->SetLineStyle(kDashed);
        hs->Add(h_mc_true_1d, drawOption);
        hs->Add(h_mc_true_1d, "hist");
        // // hs->Add(h_data_refolded_1d, drawOption);
        // // hs->Add(h_data_refolded_1d, "hist");

        h_X_inclusive->SetMarkerStyle(kOpenTriangleUp);
        h_X_inclusive->SetMarkerColor(kRed);
        h_X_inclusive->SetLineColor(kRed);
        leg->AddEntry(h_X_inclusive, "Particle level MC (inclusive jets)", "pe1");

        h_X_2B->SetMarkerStyle(kFullTriangleUp);
        h_X_2B->SetMarkerColor(kBlue);
        h_X_2B->SetLineColor(kBlue);
        TH1D *h_X_2B2 = (TH1D *) h_X_2B->Clone("h_X_2B2");
        h_X_2B2->SetFillColor(kBlue);
        h_X_2B2->SetFillStyle(3001);

        h_X_inclusive->SetLineStyle(kDashed);
        hs->Add(h_X_inclusive,drawOption);
        hs->Add(h_X_inclusive,"hist");
        // hs->Add(h_X_2B,drawOption);
        // hs->Add(h_X_2B,"hist");
        // hs->Add(h_X_2B2,"e2");
        hs->Draw("nostack");

        leg->Draw();

        TLatex *test_info_text = new TLatex;
        test_info_text->SetNDC();
        // test_info_text->SetTextSize(text_size);
        if (isData==1) test_info_text->DrawLatex(0.65, 0.93, "BOTTOMLINE TEST");
        else test_info_text->DrawLatex(0.65, 0.93, "TRIVIAL TEST");
        // test_info_text->Draw();
        // drawHeader();    

        TLine *line = new TLine(x_min, 1., x_max, 1.);
        line->SetLineWidth(2.); 
        line->SetLineStyle(kDashed);
        line->SetLineColor(kGray);

        TLegend *leg_ratio = new TLegend(0.56, 0.75, 0.8, 0.925);
        leg_ratio->SetBorderSize(0);
        leg_ratio->SetFillStyle(0);

        TH1D *h_data_mc_reco_ratio = (TH1D *) h_data_purity_corrected_1d->Clone("h_data_mc_reco_ratio");
        h_data_mc_reco_ratio->Divide(h_mc_reco_1d);
        h_data_mc_reco_ratio->SetMarkerStyle(kFullCircle);
        h_data_mc_reco_ratio->SetMarkerColor(kBlack);
        h_data_mc_reco_ratio->SetLineColor(kBlack);
        h_data_mc_reco_ratio->SetMarkerSize(1);
        h_data_mc_reco_ratio->GetYaxis()->SetRangeUser(0.,2.);
        h_data_mc_reco_ratio->GetYaxis()->SetTitle("ratio");
        h_data_mc_reco_ratio->GetXaxis()->SetTitle(xlabel);
        h_data_mc_reco_ratio->GetXaxis()->SetTitleOffset(3.5);
        h_data_mc_reco_ratio->GetYaxis()->SetNdivisions(8);
        // leg_ratio->AddEntry(h_data_mc_reco_ratio, "reco data / reco mc", "pe1");

        TH1D *h_data_mc_true_ratio = (TH1D *) h_data_fully_corrected_1d->Clone("h_data_mc_true_ratio");
        h_data_mc_true_ratio->Divide(h_mc_true_1d);
        h_data_mc_true_ratio->SetMarkerStyle(kFullCross);
        h_data_mc_true_ratio->SetMarkerColor(kBlack);
        h_data_mc_true_ratio->SetLineColor(kBlack);
        h_data_mc_true_ratio->SetMarkerSize(1);
        leg_ratio->AddEntry(h_data_mc_true_ratio, "unfolded data / true mc (double b jets)", "pe1");

        TH1D *h_data_mc_true_ratio_inclusive = (TH1D *) h_data_fully_corrected_inclusive_1d->Clone("h_data_mc_true_ratio_inclusive");
        h_data_mc_true_ratio_inclusive->Divide(h_X_inclusive);
        h_data_mc_true_ratio_inclusive->SetMarkerStyle(kOpenCross);
        h_data_mc_true_ratio_inclusive->SetMarkerColor(kBlue);
        h_data_mc_true_ratio_inclusive->SetLineColor(kBlue);
        h_data_mc_true_ratio_inclusive->SetMarkerSize(1);
        h_data_mc_true_ratio_inclusive->GetYaxis()->SetRangeUser(0.,2.);
        leg_ratio->AddEntry(h_data_mc_true_ratio_inclusive, "unfolded data / true mc (inclusive jets)", "pe1");

        pad1->cd();
        h_data_mc_true_ratio->SetStats(0);
        // h_data_mc_reco_ratio->Draw(drawOption);
        h_data_mc_true_ratio->Draw(drawOption);
        h_data_mc_true_ratio_inclusive->Draw(drawOption);
        leg_ratio->Draw();
        line->Draw();

        c_unfold->cd();
        pad1->Draw();
        pad2->Draw();
    }

    // // ---- Covariance after unfolding
    // TMatrixD covariance_matrix = unfold.Ereco();
    // TH2D *covariance_histogram = new TH2D(covariance_matrix);
    // covariance_histogram->GetXaxis()->SetTitle("Particle level " + xlabel + " * p_{T}^{jet} bins");
    // covariance_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + " * p_{T}^{jet} bins");
    // TCanvas *c_cov = new TCanvas("c_cov", "covariance", 800, 600);
    // covariance_histogram->Draw("colz");

    // // ---- Correlation after unfolding
    // TH2D *correlation_histogram = (TH2D *) covariance_histogram->Clone("correlation_histogram");
    // for (int ibinx=0; ibinx <= correlation_histogram->GetNbinsX()+1; ibinx++) {
    //     for (int ibiny=0; ibiny <= correlation_histogram->GetNbinsY()+1; ibiny++) {
    //         double covxy = covariance_histogram->GetBinContent(ibinx, ibiny);
    //         double covxx = covariance_histogram->GetBinContent(ibinx, ibinx);
    //         double covyy = covariance_histogram->GetBinContent(ibiny, ibiny);
    //         correlation_histogram->SetBinContent(ibinx, ibiny, covxy/std::sqrt(covxx*covyy));
    //     }
    // }
    // correlation_histogram->GetXaxis()->SetTitle("Particle level " + xlabel + " * p_{T}^{jet} bins");
    // correlation_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + " * p_{T}^{jet} bins");
    // correlation_histogram->GetZaxis()->SetRangeUser(-1.,1.);
    // TCanvas *c_cor = new TCanvas("c_cor", "correlation", 800, 600);
    // correlation_histogram->Draw("colz");

    // for (int ibin_pt=1; ibin_pt<=nbins_pt; ibin_pt++) {
    //     double y = ibin_pt*nbins_x;
    //     // std::cout << nbins_pt*nbins_x << std::endl;
    //     TLine *line_cory = new TLine(0, y, nbins_pt*nbins_x, y);
    //     TLine *line_corx = new TLine(y, 0, y, nbins_pt*nbins_x);

    //     c_cov->cd();
    //     line_cory->Draw();
    //     line_corx->Draw();

    //     c_cor->cd();
    //     line_cory->Draw();
    //     line_corx->Draw();
    // }

    // if (sample_unfolding.Contains("herwig")) drawHeaderHerwig();
    // else drawHeaderSimulation();
    // c_cor->Draw();
    // // c_cor->Print("plots_an/"+sample_unfolding+"_"+label_unfolding+"_correlation_"+observable+".png");

    std::cout << "Creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    h_data_fully_corrected->SetName("h_data_unfolded");
    h_data_fully_corrected->Write();

    h_data_after_fit->SetName("h_data_singleb");
    h_data_after_fit->Write();

    h_mc_reco->SetName("h_mc_reco_singleb");
    h_mc_reco->Write();

    // c_unfold->BuildLegend();

    // // JACK-KNIFE RESAMPLING FOR MC STAT UNC
    // for (int i=0; i<10; i++) {
    //     std::cout << "JK resampling i=" << i << std::endl;
    //     TH2D *h_purity_jk = (TH2D *) fin_response_unfolding->Get("h"+TString(Form("%d",i))+"_purity_"+observable+"pt");
    //     TH2D *h_efficiency_jk = (TH2D *) fin_response_unfolding->Get("h"+TString(Form("%d",i))+"_efficiency_"+observable+"pt");
    //     RooUnfoldResponse *response_jk = (RooUnfoldResponse *) fin_response_unfolding->Get("response"+TString(Form("%d",i))+"_"+observable+"pt");
 
    //     TH2D *h_purity_corrected_jk = (TH2D *) h_data_after_fit->Clone(Form("h_purity_corrected_jk_%d",i));
    //     h_purity_corrected_jk->Multiply(h_data_after_fit, h_purity_jk);

    //     RooUnfoldInvert unfold_jk(response_jk, h_purity_corrected_jk);
    //     TH2D *h_unfolded_jk = (TH2D *) unfold_jk.Hreco(errorTreatment);

    //     TH2D *h_efficiency_corrected_jk = (TH2D *) h_unfolded_jk->Clone(Form("h_efficiency_corrected_jk_%d",i));
    //     h_efficiency_corrected_jk->Divide(h_unfolded_jk, h_efficiency_jk);

    //     TH2D *h_fully_corrected = (TH2D *) h_efficiency_corrected_jk->Clone(Form("h_fully_corrected_jk_%d",i));
    //     if (!inclusive) h_fully_corrected->Divide(h_efficiency_corrected_jk, h_eff_withSF);

    //     h_fully_corrected->Write();
    // }

    fout->Close();
    delete fout;

    // gApplication -> Terminate(0);
}
