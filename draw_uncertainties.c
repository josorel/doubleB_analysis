#include <cmath>

double calculateStandardDeviation(const vector<Double_t>& arr)
{
    double sum = 0.0, mean, standardDeviation = 0.0;

    int size = arr.size();

    // Calculate the sum of elements in the vector
    for (int i = 0; i < size; ++i) {
        sum += arr[i];
    }

    // Calculate the mean
    mean = sum / size;

    // Calculate the sum of squared differences from the
    // mean
    for (int i = 0; i < size; ++i) {
        standardDeviation += pow(arr[i] - mean, 2);
    }

    // Calculate the square root of the variance to get the
    // standard deviation
    return sqrt(standardDeviation / size);
}


void draw_uncertainties(){

    TCanvas *c1 = new TCanvas("c1","c1",600,600);

    TString observable="zg";

    TFile *file0 = TFile::Open("unfoldedHistograms_Central"+observable+".root");
    
    TFile *file1 = TFile::Open("unfoldedHistograms_var1"+observable+".root");
    TFile *file2 = TFile::Open("unfoldedHistograms_var2"+observable+".root");
    TFile *file3 = TFile::Open("unfoldedHistograms_var3"+observable+".root");
    TFile *file4 = TFile::Open("unfoldedHistograms_var4"+observable+".root");
    TFile *file5 = TFile::Open("unfoldedHistograms_var5"+observable+".root");
    TFile *file6 = TFile::Open("unfoldedHistograms_var6"+observable+".root");
    TFile *file7 = TFile::Open("unfoldedHistograms_var7"+observable+".root");
    TFile *file8 = TFile::Open("unfoldedHistograms_var8"+observable+".root");
    TFile *file9 = TFile::Open("unfoldedHistograms_var9"+observable+".root");
    TFile *file10 = TFile::Open("unfoldedHistograms_var10"+observable+".root");

    TFile *file11 = TFile::Open("unfoldedHistograms_NoGlue"+observable+".root");

    TFile *file12 = TFile::Open("unfoldedHistograms_GlueX2"+observable+".root");

    TFile *file13 = TFile::Open("unfoldedHistograms_TrkEff"+observable+".root");

    TFile *file14 = TFile::Open("unfoldedHistograms_JER"+observable+".root");

    TFile *file15 = TFile::Open("unfoldedHistograms_Herwig"+observable+".root");

    TFile *file16 = TFile::Open("unfoldedHistograms_HerwigTemplate"+observable+".root");

    TH1D *hist0 = (TH1D*)file0->Get("h_data_unfolded_var");
    TH1D *hist1 = (TH1D*)file1->Get("h_data_unfolded_var");
    TH1D *hist2 = (TH1D*)file2->Get("h_data_unfolded_var");
    TH1D *hist3 = (TH1D*)file3->Get("h_data_unfolded_var");
    TH1D *hist4 = (TH1D*)file4->Get("h_data_unfolded_var");
    TH1D *hist5 = (TH1D*)file5->Get("h_data_unfolded_var");
    TH1D *hist6 = (TH1D*)file6->Get("h_data_unfolded_var");
    TH1D *hist7 = (TH1D*)file7->Get("h_data_unfolded_var");
    TH1D *hist8 = (TH1D*)file8->Get("h_data_unfolded_var");
    TH1D *hist9 = (TH1D*)file9->Get("h_data_unfolded_var");
    TH1D *hist10 = (TH1D*)file10->Get("h_data_unfolded_var");

    TH1D *hist11 = (TH1D*)file11->Get("h_data_unfolded_var");
    TH1D *hist12 = (TH1D*)file12->Get("h_data_unfolded_var");
    TH1D *hist13 = (TH1D*)file13->Get("h_data_unfolded_var");
    TH1D *hist14 = (TH1D*)file14->Get("h_data_unfolded_var");
    TH1D *hist15 = (TH1D*)file15->Get("h_data_unfolded_var");
    TH1D *hist16 = (TH1D*)file16->Get("h_data_unfolded_var");

    TH1D *hist20 = (TH1D *) hist0->Clone("hist20");
    for (Int_t i=1; i<hist0->GetNbinsX()+1;i++){
        hist20->SetBinContent(i,hist0->GetBinError(i));
    }

    hist1->Add(hist0,-1);
    hist2->Add(hist0,-1);
    hist3->Add(hist0,-1);
    hist4->Add(hist0,-1);
    hist5->Add(hist0,-1);
    hist6->Add(hist0,-1);
    hist7->Add(hist0,-1);
    hist8->Add(hist0,-1);
    hist9->Add(hist0,-1);
    hist10->Add(hist0,-1);

    hist1->Divide(hist0);
    hist2->Divide(hist0);
    hist3->Divide(hist0);
    hist4->Divide(hist0);
    hist5->Divide(hist0);
    hist6->Divide(hist0);
    hist7->Divide(hist0);
    hist8->Divide(hist0);
    hist9->Divide(hist0);
    hist10->Divide(hist0);

    for (auto h: {hist1,hist2,hist3,hist4,hist5,hist6,hist7,hist8,hist9,hist10}
    ){
        for (Int_t i=1; i<h->GetNbinsX()+1;i++) h->SetBinContent(i,abs(h->GetBinContent(i)));
    }

    map <Int_t, vector<Double_t>> statunc;
    vector<Double_t> final_statunc;
    for (Int_t i=0;i<hist0->GetNbinsX();i++){
        for (auto h: {hist1,hist2,hist3,hist4,hist5,hist6,hist7,hist8,hist9,hist10}
        ){
            if (statunc.count(i)) statunc[i].push_back(h->GetBinContent(i+1));
            else statunc.insert({i,{h->GetBinContent(i+1)}});
        }
    }
    for (auto i : statunc) final_statunc.push_back(calculateStandardDeviation(i.second)*10/9);

    TH1D *hist30 = (TH1D *) hist0->Clone("hist30");
    for (Int_t i=1; i<hist0->GetNbinsX()+1;i++){
        hist30->SetBinContent(i,final_statunc[i-1]);
    }

    hist11->Add(hist0,-1);
    hist12->Add(hist0,-1);
    hist13->Add(hist0,-1);
    hist14->Add(hist0,-1);
    hist15->Add(hist0,-1);
    hist16->Add(hist0,-1);

    hist11->Divide(hist0);
    hist12->Divide(hist0);
    hist13->Divide(hist0);
    hist14->Divide(hist0);
    hist15->Divide(hist0);
    hist16->Divide(hist0);

    hist11->Add(hist12);

    hist0->SetLineColor(kBlack);
    hist1->SetLineColor(kYellow);
    hist2->SetLineColor(kOrange);
    hist3->SetLineColor(kMagenta);
    hist4->SetLineColor(kRed-9);
    hist5->SetLineColor(kViolet-7);
    hist6->SetLineColor(kBlue);
    hist7->SetLineColor(kCyan);
    hist8->SetLineColor(kGreen);
    hist9->SetLineColor(kGreen-8);
    hist10->SetLineColor(kCyan-5);

    hist11->SetLineColor(kRed);
    hist12->SetLineColor(kYellow);
    hist13->SetLineColor(kCyan);
    hist14->SetLineColor(kMagenta);
    hist15->SetLineColor(kRed-9);
    hist16->SetLineColor(kViolet-7);
    hist30->SetLineColor(kGreen);

    hist11->SetMarkerStyle(kOpenCircle);
    hist13->SetMarkerStyle(kFullCircle);
    hist14->SetMarkerStyle(kOpenTriangleUp);
    hist15->SetMarkerStyle(kFullTriangleUp);
    hist16->SetMarkerStyle(kFullCross);
    hist20->SetMarkerStyle(0);
    hist30->SetMarkerStyle(kOpenCross);

    hist11->SetMarkerColor(kRed);
    hist13->SetMarkerColor(kCyan);
    hist14->SetMarkerColor(kMagenta);
    hist15->SetMarkerColor(kRed-9);
    hist16->SetMarkerColor(kViolet-7);
    hist30->SetMarkerColor(kGreen);

    // hist0->SetTitle("All Events");
    // hist1->SetTitle("Variation 1");
    // hist2->SetTitle("Variation 2");
    // hist3->SetTitle("Variation 3");
    // hist4->SetTitle("Variation 4");
    // hist5->SetTitle("Variation 5");
    // hist6->SetTitle("Variation 6");
    // hist7->SetTitle("Variation 7");
    // hist8->SetTitle("Variation 8");
    // hist9->SetTitle("Variation 9");
    // hist10->SetTitle("Variation 10");

    hist11->SetTitle("Light and charm misd. rate");
    // hist12->SetTitle("");
    hist13->SetTitle("Track efficiency");
    hist14->SetTitle("Jet energy resolution");
    hist15->SetTitle("Response matrix model dependence");
    hist16->SetTitle("B-jet fraction model dependence");
    hist30->SetTitle("Response matrix stat.");

    hist20->SetTitle("Statistical uncertainty");

    hist0->SetFillColor(0);
    hist1->SetFillColor(0);
    hist2->SetFillColor(0);
    hist3->SetFillColor(0);
    hist4->SetFillColor(0);
    hist5->SetFillColor(0);
    hist6->SetFillColor(0);
    hist7->SetFillColor(0);
    hist8->SetFillColor(0);
    hist9->SetFillColor(0);
    hist10->SetFillColor(0);

    hist20->SetStats(0);
    hist20->SetLineColor(kGray);
    hist20->SetFillColor(kGray);
    // hist20->SetFillStyle(kDashed);
    // for (Int_t i=1; i<hist20->GetNbinsX()+1;i++) hist20->SetBinError(i,0);
    hist20->Draw("hist");
    // hist20->Draw("pe1 Same");

    for (auto h: {hist1,hist2,hist3,hist4,hist5,hist6,hist7,hist8,hist9,hist10,hist13,hist15,hist30, hist11, hist12,hist14,hist16}//11,12,14,16
    ){
        h->Multiply(hist0);
        h->SetFillColor(0);
        for (Int_t i=1; i<h->GetNbinsX()+1;i++) h->SetBinContent(i,abs(h->GetBinContent(i)));
        for (Int_t i=1; i<h->GetNbinsX()+1;i++) h->SetBinError(i,0);
    }

    // hist13->SetStats(0);
    hist0->SetMarkerStyle(kFullCircle);
    hist0->SetMarkerColor(kBlack);
    // hist0->Draw("p1");
    // hist1->Draw("hist Same");
    // hist2->Draw("hist Same");
    // hist3->Draw("hist Same");
    // hist4->Draw("hist Same");
    // hist5->Draw("hist Same");
    // hist6->Draw("hist Same");
    // hist7->Draw("hist Same");
    // hist8->Draw("hist Same");
    // hist9->Draw("hist Same");
    // hist10->Draw("hist Same");
    hist11->Draw("p L Same");
    // hist12->Draw("hist Same");
    hist13->Draw("p L Same");
    hist14->Draw("p L Same");
    hist15->Draw("p L Same");
    hist16->Draw("p L Same");
    hist30->Draw("p L Same");

    c1->BuildLegend();
}


