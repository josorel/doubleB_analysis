
void run_templateFit(TString observable="Zg"){

    TString fout_name = "histogramsDataFittedNewBinsWithLow" + observable + "_test.root";
    if (observable=="Kt") fout_name = "histogramsDataFittedNewBinsWithLowkt_test.root";
    if (observable=="Zg") fout_name = "histogramsDataFittedNewBinsWithLowzg_test.root";
    if (observable=="Rg") fout_name = "histogramsDataFittedNewBinsWithLowrg_test.root";

    TCanvas *c1 = new TCanvas("c1","c1",600,600);

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
    const Int_t kt_binsVectorSize = 5;
    Int_t kt_bins = kt_binsVectorSize - 1;
    Double_t kt_binsVector[kt_binsVectorSize] = {
        // -1.0,
        0., 
        0.4,
        0.8,
        1.2,
        3.0
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

    TH2D *signalFractions;
    if (observable=="Kt") signalFractions = new TH2D("signalFractions","x=signalFractions",kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
    if (observable=="Zg") signalFractions = new TH2D("signalFractions","x=signalFractions",zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
    if (observable=="Rg") signalFractions = new TH2D("signalFractions","x=signalFractions",logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);

    TFile *file = TFile::Open("histogramsIdealSvtxDataHighNewBins.root");

    TH3D *hist2b_Y_A_3D = (TH3D*)file->Get("h_InvMass_vs_" + observable + "_Data");

    TFile *file3 = TFile::Open("histogramsIdealSvtxDataLowNewBins.root");

    TH3D *hist2b_Y_A2_3D = (TH3D*)file3->Get("h_InvMass_vs_" + observable + "_Data");

    hist2b_Y_A_3D->Add(hist2b_Y_A2_3D);

    TFile *file2 = TFile::Open("histogramsSvtxGenBjetUnfoldingMC_var0.root");

    TH3D *hist2b_Y_B2_3D = (TH3D*)file2->Get("h_InvMass_vs_" + observable + "_2b");
    TH3D *hist1b_Y_B2_3D = (TH3D*)file2->Get("h_InvMass_vs_" + observable + "_1b");
    TH3D *hist0b_Y_B2_3D = (TH3D*)file2->Get("h_InvMass_vs_" + observable + "_0b");
    TH3D *histMoreb_Y_B2_3D = (TH3D*)file2->Get("h_InvMass_vs_" + observable + "_Moreb");

    TFile *file4 = TFile::Open("histogramsSvtxGenDijetUnfoldingMC_var0.root");

    TH3D *hist2b_Y_B_3D = (TH3D*)file4->Get("h_InvMass_vs_" + observable + "_2b");
    TH3D *hist1b_Y_B_3D = (TH3D*)file4->Get("h_InvMass_vs_" + observable + "_1b");
    TH3D *hist0b_Y_B_3D = (TH3D*)file4->Get("h_InvMass_vs_" + observable + "_0b");
    TH3D *histMoreb_Y_B_3D = (TH3D*)file4->Get("h_InvMass_vs_" + observable + "_Moreb");

    hist2b_Y_B_3D->Add(hist2b_Y_B2_3D);
    hist1b_Y_B_3D->Add(hist1b_Y_B2_3D);
    hist0b_Y_B_3D->Add(hist0b_Y_B2_3D);
    histMoreb_Y_B_3D->Add(histMoreb_Y_B2_3D);

    hist2b_Y_B_3D->Scale(hist2b_Y_B2_3D->Integral()/hist2b_Y_B_3D->Integral());
    hist1b_Y_B_3D->Scale(hist1b_Y_B2_3D->Integral()/hist1b_Y_B_3D->Integral());
    hist0b_Y_B_3D->Scale(hist0b_Y_B2_3D->Integral()/hist0b_Y_B_3D->Integral());
    histMoreb_Y_B_3D->Scale(histMoreb_Y_B2_3D->Integral()/histMoreb_Y_B_3D->Integral());

    for (Int_t jptbin = 1; jptbin <= 3; jptbin++){
        hist2b_Y_A_3D->GetZaxis()->SetRange(jptbin,jptbin);
        TH2D *hist2b_Y_A = (TH2D*) hist2b_Y_A_3D->Project3D("xy");

        hist2b_Y_B_3D->GetZaxis()->SetRange(jptbin,jptbin);
        TH2D *hist2b_Y_B = (TH2D*) hist2b_Y_B_3D->Project3D("xy");
        hist1b_Y_B_3D->GetZaxis()->SetRange(jptbin,jptbin);
        TH2D *hist1b_Y_B = (TH2D*) hist1b_Y_B_3D->Project3D("xy");
        hist0b_Y_B_3D->GetZaxis()->SetRange(jptbin,jptbin);
        TH2D *hist0b_Y_B = (TH2D*) hist0b_Y_B_3D->Project3D("xy");
        histMoreb_Y_B_3D->GetZaxis()->SetRange(jptbin,jptbin);
        TH2D *histMoreb_Y_B = (TH2D*) histMoreb_Y_B_3D->Project3D("xy");

        for (Int_t ibin_y = 1; ibin_y <= hist2b_Y_A->GetNbinsY()+1; ibin_y++) {

            std::cout << "\n\n\n\n -------- Beginning the fit of (" << ibin_y << ")\n\n\n\n" << std::endl; 
    
            TH1D *h_data_mSv = (TH1D *) hist2b_Y_A->ProjectionY(Form("h_data_mSv_%d", ibin_y), ibin_y, ibin_y);
            TH1D *h_2B_mSv = (TH1D *) hist2b_Y_B->ProjectionY(Form("h_2B_mSv_%d", ibin_y), ibin_y, ibin_y);
            TH1D *h_1B_mSv = (TH1D *) hist1b_Y_B->ProjectionY(Form("h_1B_mSv_%d", ibin_y), ibin_y, ibin_y);
            TH1D *h_0B_mSv = (TH1D *) hist0b_Y_B->ProjectionY(Form("h_0B_mSv_%d", ibin_y), ibin_y, ibin_y);
            TH1D *h_MoreB_mSv = (TH1D *) histMoreb_Y_B->ProjectionY(Form("h_MoreB_mSv_%d", ibin_y), ibin_y, ibin_y);
    
            double int0 = h_2B_mSv->Integral(1, 7);
            double int1 = h_1B_mSv->Integral(1, 7);
            double int2 = h_0B_mSv->Integral(1, 7);
            double int3 = h_MoreB_mSv->Integral(1, 7);
    
            double bkg_fraction = (int1+int3+int2) / (int0 + int1 + int3+int2);
            cout << "bkg_fraction: " << bkg_fraction << endl;
    
            // Create the observable
            RooRealVar mSv(Form("mSv_%d", ibin_y), "mSv", 0, 7);

            h_1B_mSv->Add(h_0B_mSv);
            h_1B_mSv->Add(h_MoreB_mSv);
            h_1B_mSv->Add(h_0B_mSv);
            h_1B_mSv->Add(h_MoreB_mSv);
    
            // Create the RooDataHist object for the observed data + templates
            RooDataHist *dh_data_mSv = new RooDataHist(Form("dh_data_mSv_%d", ibin_y), "dh_data_mSv", mSv, RooFit::Import(*h_data_mSv));
            RooDataHist *dh_sig_mSv = new RooDataHist(Form("dh_sig_mSv_%d", ibin_y), "dh_sig_mSv", mSv, RooFit::Import(*h_2B_mSv));
            RooDataHist *dh_bkg_mSv = new RooDataHist(Form("dh_bkg_mSv_%d", ibin_y), "dh_bkg_mSv", mSv, RooFit::Import(*h_1B_mSv));
    
            // Create the RooHistPdf objects for the template PDFs
            RooHistPdf sig_template(Form("sig_template_%d", ibin_y), "sig_template", mSv, *dh_sig_mSv);
            RooHistPdf bkg_template(Form("bkg_template_%d", ibin_y), "bkg_template", mSv, *dh_bkg_mSv);
            RooArgList template_list(sig_template, bkg_template, "template_list");
    
            // Create the RooRealVar for the fit parameter (e.g., fraction of template A)
            RooRealVar sig_fraction_val(Form("sig_fraction_val_%d", ibin_y),"sig_fraction_val",1-bkg_fraction,0.,1.);
    
            // Create the composite PDF using a linear combination of the template PDFs
            RooAddPdf model0(Form("mode10_%d", ibin_y), "model0", template_list, sig_fraction_val, true);
            RooFitResult* result = model0.fitTo(*dh_data_mSv, RooFit::SumW2Error(true), RooFit::Save(), RooFit::CloneData(true), RooFit::PrintLevel(2), RooFit::Strategy(1), RooFit::Minos(false)); // result is already given a unique name
            Int_t status = result->status();
            result->Print();
            // mode10.plotOn(frame, DrawOption("L"));
    
            // std::cout << "covariance matrix:" << std::endl;
            // (result->covarianceMatrix().Print());
    
            // Get the fitted parameter values
            double a = sig_fraction_val.getValV();
            double da = sig_fraction_val.getError();
    
            signalFractions->SetBinContent(ibin_y,jptbin,sig_fraction_val.getValV());
            signalFractions->SetBinError(ibin_y,jptbin,sig_fraction_val.getError());
            
        }
    }

    hist2b_Y_A_3D->GetZaxis()->SetRange(1,3);
    TH2D *h_data_2d_proj = (TH2D *) hist2b_Y_A_3D->Project3D("zy");
    TH2D *h_data_2d = (TH2D *) h_data_2d_proj->Clone("h_data_2d");

    std::cout << "\n(Re)creating file " << fout_name << std::endl;
    TFile *fout = new TFile(fout_name, "recreate");

    TH1D *h_data_1d = (TH1D *) h_data_2d->ProjectionX("h_data_1d",2,2);
    TH1D *sf1d = (TH1D *) signalFractions->ProjectionX("sf1d",2,2);

    sf1d->Write();
    h_data_1d->Write();
    h_data_2d->Write();
    signalFractions->Write();

    signalFractions->Draw("colz text");

    fout->Close();
}