// #include "draw_utils.h"
// #include "RooUnfold.h"
#include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"
// #include "RooUnfoldBinByBin.h"
// #include "RooUnfoldInvert.h"
// #include "RooUnfoldSvd.h"
#include "TH2.h"

void draw_response(TString observable="zg")
{
    // Float_t font_size = 26.;
    // gStyle->SetPalette(57);
    // gStyle->SetPaintTextFormat(".2f"); 

    TString xlabel;
    if (observable=="rg") xlabel = "ln(0.4/R_{g})";
    else if (observable=="zg") xlabel = "z_{g}";
    else if (observable=="zpt") xlabel = "k_{T}";

    // TString fin_name = "responseMatrixBjetNewBins.root";   
    // TString fin_name2 = "responseMatrixDijetNewBins.root";  
    TString fin_name = "responseMatrixBjetMC_run.root";

    std::cout << "File in: " << fin_name << std::endl;
    TFile *fin = new TFile(fin_name);
    // TFile *fin2 = new TFile(fin_name2);
    TH2D *h_purity;
    TH2D *h_efficiency;
    TH2D *h_purity2;
    TH2D *h_efficiency2;

    h_purity = (TH2D *) fin->Get("h_full_purity_" + observable + "pt");
    // h_purity2 = (TH2D *) fin2->Get("h_full_purity_" + observable + "pt");
    // h_purity->Add(h_purity2);
    h_efficiency = (TH2D *) fin->Get("h_full_efficiency_" + observable + "pt");
    // h_efficiency2 = (TH2D *) fin2->Get("h_full_efficiency_" + observable + "pt");
    // h_efficiency->Add(h_efficiency2);
    RooUnfoldResponse *response = (RooUnfoldResponse *) fin->Get("response_full_" + observable + "pt");
    // RooUnfoldResponse *response2 = (RooUnfoldResponse *) fin2->Get("response_full_" + observable + "pt");
    // response->Add(response2);

    // Draw purity + efficiency 
    // for (auto h : {h_purity, h_efficiency}) {
    //     h->GetXaxis()->SetTitleSize(font_size);    
    //     h->GetXaxis()->SetLabelSize(font_size);
    //     h->GetXaxis()->SetTitle(xlabel);

    //     h->GetYaxis()->SetTitleSize(font_size);    
    //     h->GetYaxis()->SetLabelSize(font_size);
    //     h->GetYaxis()->SetTitle("p_{T}^{jet}");

    //     h->SetMarkerSize(800);
    // }

    TCanvas *c_purity = new TCanvas("c_purity", "purity", 800, 600);
    h_purity->GetZaxis()->SetTitle("Reconstruction purity");
    h_purity->Draw("colz text");
    h_purity->SetStats(0);
    // if (sample.Contains("herwig")) drawHeaderHerwig();
    // else drawHeaderSimulation();
    c_purity->Draw();
    // // c_purity->Print("plots_an/"+sample+"_"+label+"_purity_"+observable+".png");

    TCanvas *c_efficiency = new TCanvas("c_efficiency", "efficiency", 800, 600);
    h_efficiency->GetZaxis()->SetTitle("Reconstruction efficiency");
    h_efficiency->Draw("colz text");
    h_efficiency->SetStats(0);
    // if (sample.Contains("herwig")) drawHeaderHerwig();
    // else drawHeaderSimulation();
    c_efficiency->Draw();
    // // c_efficiency->Print("plots_an/"+sample+"_"+label+"_efficiency_"+observable+".png");

    // Draw the response matrix
    Int_t nbins_x = h_purity->GetNbinsX(); 
    Int_t nbins_pt = h_purity->GetNbinsY(); 

    TMatrixD response_matrix = response->Mresponse();
    TH2D *response_histogram = new TH2D(response_matrix);

    if (true) {
        TCanvas *c_response = new TCanvas("c_response", "response", 800, 600);
        // c_response->SetLogz();
        response_histogram->GetXaxis()->SetTitle("Detector level " + xlabel + " * p_{T}^{jet} bins");
        response_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + " * p_{T}^{jet} bins");
        response_histogram->GetZaxis()->SetTitle("Migration probability");
        response_histogram->SetStats(0);
        response_histogram->Draw("colz");

        for (int i = 1; i < nbins_pt; i++) {
            double coord = i * nbins_x;
            TLine *vline = new TLine(coord, 0, coord, nbins_x*nbins_pt);
            vline->SetLineColor(kBlack);
            vline->SetLineWidth(2);
            vline->Draw();

            TLine *hline = new TLine(0, coord, nbins_x*nbins_pt, coord);
            hline->SetLineColor(kBlack);
            hline->SetLineWidth(2);
            hline->Draw();
        }
        // if (sample.Contains("herwig")) drawHeaderHerwig();
        // else drawHeaderSimulation();
        // c_response->Print("plots_an/"+sample+"_"+label+"_response_"+observable+".png");
    }

    // // DRAW SPECIFIC PT BIN
    // if (false) {
    //     int ibin_pt = 2;
    //     TCanvas *c_response_i = new TCanvas(Form("c_response_%d",ibin_pt), "response", 800, 600);
    //     response_histogram->GetXaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
    //     response_histogram->GetYaxis()->SetRange(((ibin_pt-1)*nbins_x)+1,ibin_pt*nbins_x);
    //     response_histogram->GetXaxis()->SetTitle("Detector level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
    //     response_histogram->GetYaxis()->SetTitle("Particle level " + xlabel + Form(", p_{T}^{jet} bin = %d", ibin_pt));
    //     for (int i=1;i<=nbins_x+1;i++) {
    //         TString label = Form("%.2f",h_purity->GetXaxis()->GetBinLowEdge(i));
    //         response_histogram->GetYaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
    //         response_histogram->GetXaxis()->ChangeLabel(i, -1, -1, -1, -1, -1, label);
    //     }
    //     response_histogram->GetZaxis()->SetRangeUser(0, 0.7);

    //     response_histogram->Draw("colz");
    //     for (int i = 1; i < nbins_pt; i++) {
    //         double coord = i * nbins_x;
    //         TLine *vline = new TLine(coord, 0, coord, nbins_x*nbins_pt);
    //         vline->SetLineColor(kBlack);
    //         vline->SetLineWidth(2);
    //         vline->Draw();

    //         TLine *hline = new TLine(0, coord, nbins_x*nbins_pt, coord);
    //         hline->SetLineColor(kBlack);
    //         hline->SetLineWidth(2);
    //         hline->Draw();
    //     }
    //     // if (sample.Contains("herwig")) drawHeaderHerwig();
    //     // else drawHeaderSimulation();
    // }

    // Print condition number
    TDecompSVD *svd= new TDecompSVD(response->Mresponse());  // response is a RooUnfold response object, svd is the singular value decomposition (SVD) matrix. the response->Mresponse() returns the normalized migration matrix
    auto singular_values = svd->GetSig(); //this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
    double cond_number = singular_values.Max() / singular_values.Min();
    for (int i=0; i < singular_values.GetNrows(); i++) {
        double val = singular_values[i];
        // std::cout << val << std::endl;
    }
    std::cout << "Largest value = " << singular_values.Max() 
              << "\nSmallest value = " << singular_values.Min()
              << "\nCondition number = " << cond_number
              << std::endl;

}