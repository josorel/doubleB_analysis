// #include "root_utils.h"
// #include "../../myMacros/C106X/draw_utils.h"
// #include "../../myMacros/C106X/cms_palette.h"

const Float_t text_size = 0.02;

void draw_info()
{
    // TLatex *prelim = new TLatex;
    // prelim->SetNDC();
    // prelim->SetTextSize(text_size);
    // prelim->SetTextAlign(12);
    // prelim->DrawLatex(0.15, 0.95, "#bf{CMS} #it{Simulation Supplementary}");

    TLatex *lumi = new TLatex;
    lumi->SetNDC();
    lumi->SetTextSize(text_size);
    lumi->SetTextAlign(32);
    lumi->DrawLatex(0.95, 0.95, "PYTHIA8 CP5 (pp 5.02 TeV)");
}

void draw_ROC_TMVA()
{
    // ---- Setup 
    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetTextSize(text_size);
    gStyle->SetLegendTextSize(text_size);
    gStyle->SetLabelSize(text_size, "XYZ");
    gStyle->SetTitleSize(text_size, "XYZ");
    
    // TString label = "qcd_bjet";
    TString label = "bjet";
    TString finName = label + "_roc.root";
    TFile *fin = new TFile (finName);

    TH1F *h_discr_s = (TH1F *) fin->Get("h_discr_s");
    TH1F *h_discr_b = (TH1F *) fin->Get("h_discr_b");

    THStack *h_discr = new THStack("h_discr", "");
    h_discr->SetTitle("; discriminator value; n_{tracks}");
    TLegend *leg_discr = new TLegend(0.4, 0.7, 0.6, 0.8);

    // Count roc graphs 
    Int_t nbins = h_discr_s->GetNbinsX();

    TGraph *gr_roc = new TGraph(nbins);
    gr_roc->SetLineColor(kBlue);
    gr_roc->SetLineWidth(3);
    gr_roc->GetXaxis()->SetTitle("Signal efficiency");
    gr_roc->GetYaxis()->SetTitle("Background rejection");

    Float_t total_sig = h_discr_s->Integral(1, nbins);
    Float_t total_bkg = h_discr_b->Integral(1, nbins);
    // std::cout << nbins << std::endl;
    // gr_roc->SetPoint(nbins, 1., 0.);
    for (Int_t ibin = 1; ibin <= nbins; ibin++) {
        float sig_pass = h_discr_s->Integral(ibin, nbins);
        float sig_eff = sig_pass / total_sig;

        float bkg_pass = h_discr_b->Integral(ibin, nbins);
        float bkg_rej = 1 - (bkg_pass / total_bkg);

        // std::cout << "edge low " << h_discr_s->GetXaxis()->GetBinLowEdge(ibin) 
        //           << ", edge high " << h_discr_s->GetXaxis()->GetBinUpEdge(ibin) 
        //           << std::endl;
        gr_roc->SetPoint(ibin-1, sig_eff, bkg_rej);

        float purity = sig_pass / (sig_pass + bkg_pass);

        // if (sig_eff > 0.92) 
            std::cout << "eff: " << sig_eff << ", bkg rej: " << bkg_rej << ", pur: " << purity << std::endl;
    }
    gr_roc->SetPoint(nbins, 0., 1.); // for discr > 1

    // format histos
    h_discr_s->SetLineWidth(2);
    h_discr_s->SetLineColor(kBlue);
    h_discr->Add(h_discr_s);
    leg_discr->AddEntry(h_discr_s, "Match to signal", "l");

    h_discr_b->SetLineWidth(2);
    h_discr_b->SetLineColor(kRed);
    h_discr->Add(h_discr_b);
    leg_discr->AddEntry(h_discr_b, "Match to background", "l");

    // TCanvas *c_discr = new TCanvas("c_discr", "", 1200, 1000);
    // h_discr->Draw("histo nostack");
    // leg_discr->Draw();
    // c_discr->Draw();
    // SetRealAspectRatio(c_discr);
    // c_discr->SetFixedAspectRatio();

    // TPaveText *info_bdt = new TPaveText(0., 1., 0.35, 1.08, "nb ndc");
    // info_bdt->SetTextSize(text_size);
    // info_bdt->SetFillStyle(0);
    // info_bdt->SetLineWidth(0);
    // info_bdt->SetBorderSize(0);
    // info_bdt->AddText("#it{Gradient BDT}");

    TPaveText *info_tracks = new TPaveText(0.5, 0.5, 0.9, 0.9, "nb ndc");
    info_tracks->SetTextSize(text_size);
    info_tracks->SetFillStyle(0);
    info_tracks->SetBorderSize(0);
    info_tracks->SetLineWidth(0);
    info_tracks->SetTextAlign(12);
    // info_tracks->AddText("#bf{Charged particle classification}");
    info_tracks->AddText("#bf{CMS}");
    info_tracks->AddText("#it{Simulation Supplementary}");
    info_tracks->AddText("");
    info_tracks->AddText("");
    info_tracks->AddText("#it{Gradient BDT}");
    info_tracks->AddText("");
    info_tracks->AddText("Charged particles,^{} p_{T}^{track} > 1 GeV");
    info_tracks->AddText("from b-tagged b jets");
    info_tracks->AddText("60 GeV < p_{T}^{jet} < 150 GeV, |#eta^{jet}| < 2");
    info_tracks->AddText("");
    info_tracks->AddText("Signal:");
    info_tracks->AddText("Jets with two SVs from two Bs");
    info_tracks->AddText("");
    info_tracks->AddText("Background:");
    info_tracks->AddText("Jets with two SVs from single Bs");
    // info_tracks->AddText("");
    // info_tracks->AddText("Gradient BDT");

    // std::cout << gr_roc->GetXaxis()->GetLabelSize() << std::endl;

    gr_roc->GetXaxis()->SetLabelSize(text_size);
    gr_roc->GetYaxis()->SetLabelSize(text_size);
    gr_roc->GetXaxis()->SetTitleSize(text_size);
    gr_roc->GetYaxis()->SetTitleSize(text_size);

    gr_roc->GetYaxis()->SetNdivisions(511);
    gr_roc->GetXaxis()->SetNdivisions(511);
    gr_roc->GetXaxis()->SetTitleOffset(1.2);

    TCanvas *c_roc = new TCanvas("c_roc", "", 900, 700);
    c_roc->SetRightMargin(0.05);
    c_roc->SetLeftMargin(0.15);
    c_roc->SetTopMargin(0.07);
    c_roc->SetBottomMargin(0.12);
    gr_roc->Draw("al");
    // info_top_left->Draw();
    // info_top_right->Draw();
    info_tracks->Draw();
    // info_bdt->Draw();

    TLine *linev = new TLine(0., 1., 1.1, 1.);
    linev->SetLineStyle(kDashed);
    linev->SetLineColorAlpha(kBlack, 0.9);
    linev->Draw("same");

    TLine *line = new TLine(1., 0., 1., 1.1);
    line->SetLineStyle(kDashed);
    line->SetLineColorAlpha(kBlack, 0.9);
    // line->Draw("same");

    gr_roc->Draw("l same");

    c_roc->SetRealAspectRatio();
    // drawHeaderSimulation();
    draw_info();
    c_roc->Draw("Same");
    // c_roc->SetGrid();
    TString c_roc_name = "./plots/" + label + "_TMVA_roc.png";
    c_roc->Print(c_roc_name);  
    TString c_roc_name_pdf = "./plots/" + label + "_TMVA_roc.pdf";
    c_roc->Print(c_roc_name_pdf);    

    // Double_t x1 = c_roc->GetX1();
    // Double_t y1 = c_roc->GetY1();
    // Double_t x2 = c_roc->GetX2();
    // Double_t y2 = c_roc->GetY2();

    // std::cout << x1 << " " << x2 << " " << y1 << " " << y2 << std::endl;

}