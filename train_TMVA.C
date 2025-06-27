
void train_TMVA() 
{
    TMVA::Tools::Instance();

    // Choose methods
    std::map<std::string,int> Use;
    Use["BDTG"] = 1; // uses Gradient Boost

    std::cout << std::endl;
    std::cout << "==> Start TMVAClassification" << std::endl;

    TFile *input(0);
    // TString label = "ttbar_highPU";
    TString label = "bjet";
    // TString fname = "./data_root_" + label + "_80_pt_140/data.root";
    // TString fname = "./ntuples/" + label + "_30_pt_700.root";
    TString fname = "/data_CMS/cms/sorel/merged_HiForestMiniAOD_60_pt_150.root";
    if (!gSystem->AccessPathName( fname )) {
        input = TFile::Open( fname ); // check if file in local directory exists
    }
    if (!input) {
        std::cout << "ERROR: could not open data file" << std::endl;
        exit(1);
    }
    std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;


    TTree *signalTree     = (TTree*)input->Get("TreeS");
    TTree *background     = (TTree*)input->Get("TreeB");

    // TTree *signalTree     = (TTree*)input->Get("TreeS");
    // TTree *background     = (TTree*)input->Get("TreeB");

    TString outfileName("/data_CMS/cms/sorel/SV_" + label + "_TMVA.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                            "!V:!Silent:Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    TMVA::DataLoader *dataloader=new TMVA::DataLoader(label + "_dataloader");

    dataloader->AddVariable("SVdist",'F');

    dataloader->AddVariable("SVmass1",'F');
    dataloader->AddVariable("SVpt1",'F');
    dataloader->AddVariable("SVntrks1",'F');
    dataloader->AddVariable("SVmass2",'F');
    dataloader->AddVariable("SVpt2",'F');
    dataloader->AddVariable("SVntrks2",'F');
    

    // dataloader->AddSpectator( "trkIsLepton", 'I' );
    // dataloader->AddSpectator( "trkInSV", 'I' );
    // dataloader->AddSpectator( "jtpt", 'F' );

    Double_t signalWeight = 1.0;
    Double_t backgroundWeight = 1.0;

    dataloader->AddSignalTree    ( signalTree,     signalWeight );
    dataloader->AddBackgroundTree( background, backgroundWeight );

    dataloader->SetSignalWeightExpression("weight");
    dataloader->SetBackgroundWeightExpression("weight");


    // TCut mycuts = "trkInSV==1"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    // TCut mycutb = "trkInSV==1"; // for example: TCut mycutb = "abs(var1)<0.5";

    TCut mycuts = "";
    TCut mycutb = mycuts;

    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
        "nTrain_Signal=150000:nTrain_Background=150000:nTest_Signal=150000:nTest_Background=150000:SplitMode=Random:NormMode=NumEvents:!V" );

    // Book methods
    if (Use["BDTG"]) // Gradient Boost
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
            "!H:!V:NTrees=1000:MinNodeSize=5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=10:MaxDepth=5" );

    // Train
    factory->TrainAllMethods();

    // Test
    factory->TestAllMethods();

    // Evaluate
    factory->EvaluateAllMethods();

    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete dataloader;

}