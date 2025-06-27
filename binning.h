#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldInvert.h"
#include "TH2.h"

// jet pt bins : no underflow or overflow
const Int_t jtpt_binsVectorSize = 4;
Int_t jtpt_bins = jtpt_binsVectorSize - 1;
Double_t jtpt_binsVector[jtpt_binsVectorSize] = {
    60., 
    80., 
    110.,
    150.
};
Double_t jtpt_min = jtpt_binsVector[0];
Double_t jtpt_max = jtpt_binsVector[jtpt_bins];

// jet pt bins w/ underflow and overflow
const Int_t jtpt_wUO_binsVectorSize = 11;
Int_t jtpt_wUO_bins = jtpt_wUO_binsVectorSize - 1;
Double_t jtpt_wUO_binsVector[jtpt_wUO_binsVectorSize] = {
    5,
    50,
    60,
    70,
    80., 
    // 90,
    100., 
    // 110,
    120.,
    // 130,
    140.,
    150,
    200,
    5000,
};
Double_t jtpt_wUO_min = jtpt_wUO_binsVector[0];
Double_t jtpt_wUO_max = jtpt_wUO_binsVector[jtpt_wUO_bins];

// ln(0.4/rg) bins : 1st bin untagged (unphysical)
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
Double_t logrg_min = logrg_binsVector[0];
Double_t logrg_max = logrg_binsVector[logrg_bins];

// ln(0.4/rg) bins w/ underflow (rg>0.4 or logrg<0) : 1st bin underflow, 2nd bin untagged (unphysical)
const Int_t logrg_wU_binsVectorSize = 10;
Int_t logrg_wU_bins = logrg_wU_binsVectorSize - 1;
Double_t logrg_wU_binsVector[logrg_wU_binsVectorSize] = {
    -2., // SD-untagged or kT<1 => use -1.1 as logrg_underflow value
    -1., // rg>0.4
    0., 
    // 0.15,
    0.3,
    // 0.45,
    0.6,
    // 0.75,
    0.9,
    // 1.05,
    1.2,
    // 1.4,
    1.6,
    // 1.8,
    2.1,
    2.5
};
Double_t logrg_wU_min = logrg_wU_binsVector[0];
Double_t logrg_wU_max = logrg_wU_binsVector[logrg_wU_bins];


// zg bins : 1st bin untagged (unphysical)
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
Double_t zg_min = zg_binsVector[0];
Double_t zg_max = zg_binsVector[zg_bins];

// kt bins (bpt/jtptCh)
const Int_t kt_binsVectorSize = 5;
Int_t kt_bins = kt_binsVectorSize - 1;
Double_t kt_binsVector[kt_binsVectorSize] = {
    // -1.0,
    0., 
    0.4,
    0.8,
    1.2,
    3.,
};
Double_t kt_min = kt_binsVector[0];
Double_t kt_max = kt_binsVector[kt_bins];
// const Int_t kt_binsVectorSize = 7;
// Int_t kt_bins = kt_binsVectorSize - 1;
// Double_t kt_binsVector[kt_binsVectorSize] = {
//     0., 
//     0.35,
//     0.55,
//     0.7,
//     0.8,
//     0.9,
//     1.
// };
// Double_t kt_min = kt_binsVector[0];
// Double_t kt_max = kt_binsVector[kt_bins];

// kt bins (bpt/jtptCh) w/ underflow (no charged particles in the jet)
const Int_t kt_wU_binsVectorSize = 8;

Int_t kt_wU_bins = kt_wU_binsVectorSize - 1;
Double_t kt_wU_binsVector[kt_wU_binsVectorSize] = {
    -0.1, // jtptCh==0
    0., 
    0.35,
    0.55,
    0.7,
    0.8,
    0.9,
    1.
};
Double_t kt_wU_min = kt_wU_binsVector[0];
Double_t kt_wU_max = kt_wU_binsVector[kt_wU_bins];

// Double_t jtpt;
// Double_t jtptCh;
// Double_t logrg;
// Double_t logkt;
// Double_t zg;
// Double_t mb;
// Double_t bpt;
// Double_t jtpt_gen;
// Double_t jtptCh_gen;
// Double_t logrg_gen;
// Double_t logkt_gen;
// Double_t zg_gen;
// Double_t mb_gen;
// Double_t bpt_gen;
// Double_t weight;
// Double_t pthat;
// Double_t jer_sf_nom;
// Double_t jer_sf_up;
// Double_t jer_sf_down;
// Double_t jec_unc;  

// Declare histograms
TH2D *h_half0_purity_numerator_rgpt = new TH2D("h_half0_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_purity_denominator_rgpt = new TH2D("h_half0_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_efficiency_numerator_rgpt = new TH2D("h_half0_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_efficiency_denominator_rgpt = new TH2D("h_half0_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_purity_numerator_rgpt = new TH2D("h_half1_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_purity_denominator_rgpt = new TH2D("h_half1_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_efficiency_numerator_rgpt = new TH2D("h_half1_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_efficiency_denominator_rgpt = new TH2D("h_half1_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);

RooUnfoldResponse *response_half0_rgpt = new RooUnfoldResponse(h_half0_purity_denominator_rgpt, h_half0_efficiency_denominator_rgpt, "response_rgpt_half0", "response for 2d: rg and jet pt"); 
RooUnfoldResponse *response_half1_rgpt = new RooUnfoldResponse(h_half0_purity_denominator_rgpt, h_half0_efficiency_denominator_rgpt, "response_rgpt_half1", "response for 2d: rg and jet pt"); 
RooUnfoldResponse *response_full_rgpt = new RooUnfoldResponse(h_half0_purity_denominator_rgpt, h_half0_efficiency_denominator_rgpt, "response_full_rgpt", "response for 2d: rg and jet pt"); 

TH2D *h_half0_purity_numerator_zgpt = new TH2D("h_half0_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_purity_denominator_zgpt = new TH2D("h_half0_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_efficiency_numerator_zgpt = new TH2D("h_half0_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_efficiency_denominator_zgpt = new TH2D("h_half0_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_purity_numerator_zgpt = new TH2D("h_half1_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_purity_denominator_zgpt = new TH2D("h_half1_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_efficiency_numerator_zgpt = new TH2D("h_half1_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_efficiency_denominator_zgpt = new TH2D("h_half1_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);

RooUnfoldResponse *response_half0_zgpt = new RooUnfoldResponse(h_half0_purity_denominator_zgpt, h_half0_efficiency_denominator_zgpt, "response_zgpt_half0", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response_half1_zgpt = new RooUnfoldResponse(h_half0_purity_denominator_zgpt, h_half0_efficiency_denominator_zgpt, "response_zgpt_half1", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response_full_zgpt = new RooUnfoldResponse(h_half0_purity_denominator_zgpt, h_half0_efficiency_denominator_zgpt, "response_full_zgpt", "response for 2d: zg and jet pt"); 

TH2D *h_half0_purity_numerator_ktpt = new TH2D("h_half0_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_purity_denominator_ktpt = new TH2D("h_half0_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_efficiency_numerator_ktpt = new TH2D("h_half0_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half0_efficiency_denominator_ktpt = new TH2D("h_half0_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_purity_numerator_ktpt = new TH2D("h_half1_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_purity_denominator_ktpt = new TH2D("h_half1_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_efficiency_numerator_ktpt = new TH2D("h_half1_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h_half1_efficiency_denominator_ktpt = new TH2D("h_half1_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);

RooUnfoldResponse *response_half0_ktpt = new RooUnfoldResponse(h_half0_purity_denominator_ktpt, h_half0_efficiency_denominator_ktpt, "response_ktpt_half0", "response for 2d: kt and jet pt"); 
RooUnfoldResponse *response_half1_ktpt = new RooUnfoldResponse(h_half0_purity_denominator_ktpt, h_half0_efficiency_denominator_ktpt, "response_ktpt_half1", "response for 2d: kt and jet pt"); 
RooUnfoldResponse *response_full_ktpt = new RooUnfoldResponse(h_half0_purity_denominator_ktpt, h_half0_efficiency_denominator_ktpt, "response_full_ktpt", "response for 2d: kt and jet pt"); 

// rg jk resampling histogram definiton
TH2D *h0_purity_numerator_rgpt = new TH2D("h"+TString(Form("%d",0))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",1))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",2))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",3))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h4_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",4))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h5_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",5))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h6_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",6))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h7_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",7))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h8_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",8))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h9_purity_numerator_rgpt= new TH2D("h"+TString(Form("%d",9))+"_purity_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 

std::vector<TH2D *> histos_purity_numerator_rgpt = {
    h0_purity_numerator_rgpt, 
    h1_purity_numerator_rgpt, 
    h2_purity_numerator_rgpt, 
    h3_purity_numerator_rgpt, 
    h4_purity_numerator_rgpt, 
    h5_purity_numerator_rgpt, 
    h6_purity_numerator_rgpt, 
    h7_purity_numerator_rgpt, 
    h8_purity_numerator_rgpt, 
    h9_purity_numerator_rgpt,
};

TH2D *h0_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",0))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",1))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",2))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",3))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h4_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",4))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h5_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",5))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h6_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",6))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h7_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",7))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h8_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",8))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h9_purity_denominator_rgpt=new TH2D("h"+TString(Form("%d",9))+"_purity_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_purity_denominator_rgpt = {
    h0_purity_denominator_rgpt, 
    h1_purity_denominator_rgpt, 
    h2_purity_denominator_rgpt, 
    h3_purity_denominator_rgpt, 
    h4_purity_denominator_rgpt, 
    h5_purity_denominator_rgpt, 
    h6_purity_denominator_rgpt, 
    h7_purity_denominator_rgpt, 
    h8_purity_denominator_rgpt, 
    h9_purity_denominator_rgpt,
};

TH2D *h0_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",0))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",1))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h2_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",2))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);
TH2D *h3_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",3))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h4_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",4))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h5_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",5))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h6_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",6))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h7_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",7))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h8_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",8))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h9_efficiency_numerator_rgpt=new TH2D("h"+TString(Form("%d",9))+"_efficiency_numerator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_efficiency_numerator_rgpt = {
    h0_efficiency_numerator_rgpt, 
    h1_efficiency_numerator_rgpt, 
    h2_efficiency_numerator_rgpt, 
    h3_efficiency_numerator_rgpt, 
    h4_efficiency_numerator_rgpt, 
    h5_efficiency_numerator_rgpt, 
    h6_efficiency_numerator_rgpt, 
    h7_efficiency_numerator_rgpt, 
    h8_efficiency_numerator_rgpt, 
    h9_efficiency_numerator_rgpt,
};

TH2D *h0_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",0))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",1))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",2))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",3))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h4_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",4))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h5_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",5))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h6_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",6))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h7_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",7))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h8_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",8))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h9_efficiency_denominator_rgpt= new TH2D("h"+TString(Form("%d",9))+"_efficiency_denominator_rgpt", "x=logrg, y=jtpt", logrg_bins, logrg_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_efficiency_denominator_rgpt = {
    h0_efficiency_denominator_rgpt, 
    h1_efficiency_denominator_rgpt, 
    h2_efficiency_denominator_rgpt, 
    h3_efficiency_denominator_rgpt, 
    h4_efficiency_denominator_rgpt, 
    h5_efficiency_denominator_rgpt, 
    h6_efficiency_denominator_rgpt, 
    h7_efficiency_denominator_rgpt, 
    h8_efficiency_denominator_rgpt, 
    h9_efficiency_denominator_rgpt,
};

RooUnfoldResponse *response0_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[0],histos_efficiency_denominator_rgpt[0], "response"+TString(Form("%d",0))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response1_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[1],histos_efficiency_denominator_rgpt[1], "response"+TString(Form("%d",1))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response2_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[2],histos_efficiency_denominator_rgpt[2], "response"+TString(Form("%d",2))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response3_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[3],histos_efficiency_denominator_rgpt[3], "response"+TString(Form("%d",3))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response4_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[4],histos_efficiency_denominator_rgpt[4], "response"+TString(Form("%d",4))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response5_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[5],histos_efficiency_denominator_rgpt[5], "response"+TString(Form("%d",5))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response6_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[6],histos_efficiency_denominator_rgpt[6], "response"+TString(Form("%d",6))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response7_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[7],histos_efficiency_denominator_rgpt[7], "response"+TString(Form("%d",7))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response8_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[8],histos_efficiency_denominator_rgpt[8], "response"+TString(Form("%d",8))+"_rgpt", "response for 2d: rg and jet pt");
RooUnfoldResponse *response9_rgpt= new RooUnfoldResponse(histos_purity_denominator_rgpt[9],histos_efficiency_denominator_rgpt[9], "response"+TString(Form("%d",9))+"_rgpt", "response for 2d: rg and jet pt");

std::vector<RooUnfoldResponse *> responses_rgpt = {
    response0_rgpt, 
    response1_rgpt,
    response2_rgpt,
    response3_rgpt,
    response4_rgpt,
    response5_rgpt,
    response6_rgpt,
    response7_rgpt,
    response8_rgpt,
    response9_rgpt,
};

// zg jk resampling histogram definiton
TH2D *h0_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",0))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",1))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h2_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",2))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h3_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",3))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h4_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",4))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h5_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",5))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h6_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",6))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h7_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",7))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h8_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",8))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h9_purity_numerator_zgpt= new TH2D("h"+TString(Form("%d",9))+"_purity_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_purity_numerator_zgpt = {
    h0_purity_numerator_zgpt, 
    h1_purity_numerator_zgpt, 
    h2_purity_numerator_zgpt, 
    h3_purity_numerator_zgpt, 
    h4_purity_numerator_zgpt, 
    h5_purity_numerator_zgpt, 
    h6_purity_numerator_zgpt, 
    h7_purity_numerator_zgpt, 
    h8_purity_numerator_zgpt, 
    h9_purity_numerator_zgpt,
};

TH2D *h0_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",0))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",1))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",2))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",3))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h4_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",4))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h5_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",5))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h6_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",6))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h7_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",7))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h8_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",8))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h9_purity_denominator_zgpt= new TH2D("h"+TString(Form("%d",9))+"_purity_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_purity_denominator_zgpt = {
    h0_purity_denominator_zgpt, 
    h1_purity_denominator_zgpt, 
    h2_purity_denominator_zgpt, 
    h3_purity_denominator_zgpt, 
    h4_purity_denominator_zgpt, 
    h5_purity_denominator_zgpt, 
    h6_purity_denominator_zgpt, 
    h7_purity_denominator_zgpt, 
    h8_purity_denominator_zgpt, 
    h9_purity_denominator_zgpt,
};

TH2D *h0_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",0))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",1))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",2))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",3))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h4_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",4))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h5_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",5))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h6_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",6))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h7_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",7))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h8_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",8))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h9_efficiency_numerator_zgpt= new TH2D("h"+TString(Form("%d",9))+"_efficiency_numerator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_efficiency_numerator_zgpt = {
    h0_efficiency_numerator_zgpt, 
    h1_efficiency_numerator_zgpt, 
    h2_efficiency_numerator_zgpt, 
    h3_efficiency_numerator_zgpt, 
    h4_efficiency_numerator_zgpt, 
    h5_efficiency_numerator_zgpt, 
    h6_efficiency_numerator_zgpt, 
    h7_efficiency_numerator_zgpt, 
    h8_efficiency_numerator_zgpt, 
    h9_efficiency_numerator_zgpt,
};

TH2D *h0_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",0))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",1))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",2))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",3))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h4_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",4))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h5_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",5))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h6_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",6))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h7_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",7))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h8_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",8))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h9_efficiency_denominator_zgpt= new TH2D("h"+TString(Form("%d",9))+"_efficiency_denominator_zgpt", "x=zg, y=jtpt", zg_bins, zg_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_efficiency_denominator_zgpt = {
    h0_efficiency_denominator_zgpt, 
    h1_efficiency_denominator_zgpt, 
    h2_efficiency_denominator_zgpt, 
    h3_efficiency_denominator_zgpt, 
    h4_efficiency_denominator_zgpt, 
    h5_efficiency_denominator_zgpt, 
    h6_efficiency_denominator_zgpt, 
    h7_efficiency_denominator_zgpt, 
    h8_efficiency_denominator_zgpt, 
    h9_efficiency_denominator_zgpt,
};

RooUnfoldResponse *response0_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[0], histos_efficiency_denominator_zgpt[0], "response"+TString(Form("%d",0))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response1_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[1], histos_efficiency_denominator_zgpt[1], "response"+TString(Form("%d",1))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response2_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[2], histos_efficiency_denominator_zgpt[2], "response"+TString(Form("%d",2))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response3_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[3], histos_efficiency_denominator_zgpt[3], "response"+TString(Form("%d",3))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response4_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[4], histos_efficiency_denominator_zgpt[4], "response"+TString(Form("%d",4))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response5_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[5], histos_efficiency_denominator_zgpt[5], "response"+TString(Form("%d",5))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response6_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[6], histos_efficiency_denominator_zgpt[6], "response"+TString(Form("%d",6))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response7_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[7], histos_efficiency_denominator_zgpt[7], "response"+TString(Form("%d",7))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response8_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[8], histos_efficiency_denominator_zgpt[8], "response"+TString(Form("%d",8))+"_zgpt", "response for 2d: zg and jet pt"); 
RooUnfoldResponse *response9_zgpt=new RooUnfoldResponse(histos_purity_denominator_zgpt[9], histos_efficiency_denominator_zgpt[9], "response"+TString(Form("%d",9))+"_zgpt", "response for 2d: zg and jet pt");

std::vector<RooUnfoldResponse *> responses_zgpt = {
    response0_zgpt, 
    response1_zgpt,
    response2_zgpt,
    response3_zgpt,
    response4_zgpt,
    response5_zgpt,
    response6_zgpt,
    response7_zgpt,
    response8_zgpt,
    response9_zgpt,
};

// kt jk resampling histogram definiton
TH2D *h0_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",0))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",1))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",2))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",3))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h4_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",4))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h5_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",5))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h6_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",6))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h7_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",7))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h8_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",8))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h9_purity_numerator_ktpt= new TH2D("h"+TString(Form("%d",9))+"_purity_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_purity_numerator_ktpt = {
    h0_purity_numerator_ktpt, 
    h1_purity_numerator_ktpt, 
    h2_purity_numerator_ktpt, 
    h3_purity_numerator_ktpt, 
    h4_purity_numerator_ktpt, 
    h5_purity_numerator_ktpt, 
    h6_purity_numerator_ktpt, 
    h7_purity_numerator_ktpt, 
    h8_purity_numerator_ktpt, 
    h9_purity_numerator_ktpt,
};

TH2D *h0_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",0))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",1))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h2_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",2))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h3_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",3))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h4_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",4))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h5_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",5))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h6_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",6))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h7_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",7))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h8_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",8))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);  
TH2D *h9_purity_denominator_ktpt=new TH2D("h"+TString(Form("%d",9))+"_purity_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_purity_denominator_ktpt = {
    h0_purity_denominator_ktpt, 
    h1_purity_denominator_ktpt, 
    h2_purity_denominator_ktpt, 
    h3_purity_denominator_ktpt, 
    h4_purity_denominator_ktpt, 
    h5_purity_denominator_ktpt, 
    h6_purity_denominator_ktpt, 
    h7_purity_denominator_ktpt, 
    h8_purity_denominator_ktpt, 
    h9_purity_denominator_ktpt,
};

TH2D *h0_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",0))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",1))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h2_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",2))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h3_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",3))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h4_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",4))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h5_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",5))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h6_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",6))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h7_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",7))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h8_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",8))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h9_efficiency_numerator_ktpt= new TH2D("h"+TString(Form("%d",9))+"_efficiency_numerator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_efficiency_numerator_ktpt = {
    h0_efficiency_numerator_ktpt, 
    h1_efficiency_numerator_ktpt, 
    h2_efficiency_numerator_ktpt, 
    h3_efficiency_numerator_ktpt, 
    h4_efficiency_numerator_ktpt, 
    h5_efficiency_numerator_ktpt, 
    h6_efficiency_numerator_ktpt, 
    h7_efficiency_numerator_ktpt, 
    h8_efficiency_numerator_ktpt, 
    h9_efficiency_numerator_ktpt,
};

TH2D *h0_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",0))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h1_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",1))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h2_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",2))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h3_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",3))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h4_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",4))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h5_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",5))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h6_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",6))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h7_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",7))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h8_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",8))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector); 
TH2D *h9_efficiency_denominator_ktpt=new TH2D("h"+TString(Form("%d",9))+"_efficiency_denominator_ktpt", "x=kt, y=jtpt", kt_bins, kt_binsVector, jtpt_bins, jtpt_binsVector);

std::vector<TH2D *> histos_efficiency_denominator_ktpt = {
    h0_efficiency_denominator_ktpt, 
    h1_efficiency_denominator_ktpt, 
    h2_efficiency_denominator_ktpt, 
    h3_efficiency_denominator_ktpt, 
    h4_efficiency_denominator_ktpt, 
    h5_efficiency_denominator_ktpt, 
    h6_efficiency_denominator_ktpt, 
    h7_efficiency_denominator_ktpt, 
    h8_efficiency_denominator_ktpt, 
    h9_efficiency_denominator_ktpt,
};

RooUnfoldResponse *response0_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[0], histos_efficiency_denominator_ktpt[0], "response"+TString(Form("%d",0))+"_ktpt", "response for 2d: kt and jet pt"); 
RooUnfoldResponse *response1_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[1], histos_efficiency_denominator_ktpt[1], "response"+TString(Form("%d",1))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response2_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[2], histos_efficiency_denominator_ktpt[2], "response"+TString(Form("%d",2))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response3_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[3], histos_efficiency_denominator_ktpt[3], "response"+TString(Form("%d",3))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response4_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[4], histos_efficiency_denominator_ktpt[4], "response"+TString(Form("%d",4))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response5_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[5], histos_efficiency_denominator_ktpt[5], "response"+TString(Form("%d",5))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response6_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[6], histos_efficiency_denominator_ktpt[6], "response"+TString(Form("%d",6))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response7_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[7], histos_efficiency_denominator_ktpt[7], "response"+TString(Form("%d",7))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response8_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[8], histos_efficiency_denominator_ktpt[8], "response"+TString(Form("%d",8))+"_ktpt", "response for 2d: kt and jet pt");
RooUnfoldResponse *response9_ktpt= new RooUnfoldResponse(histos_purity_denominator_ktpt[9], histos_efficiency_denominator_ktpt[9], "response"+TString(Form("%d",9))+"_ktpt", "response for 2d: kt and jet pt");

std::vector<RooUnfoldResponse *> responses_ktpt = {
    response0_ktpt, 
    response1_ktpt,
    response2_ktpt,
    response3_ktpt,
    response4_ktpt,
    response5_ktpt,
    response6_ktpt,
    response7_ktpt,
    response8_ktpt,
    response9_ktpt,
};

TH1D *h_Zg_idealAg = new TH1D("h_Zg_idealAg", "x=Zg 2b", zg_bins, zg_binsVector);
TH1D *h_Kt_idealAg = new TH1D("h_Kt_idealAg", "x=Kt 2b", kt_bins, kt_binsVector);
TH1D *h_Rg_idealAg = new TH1D("h_Rg_idealAg", "x=Rg 2b", logrg_bins, logrg_binsVector);

TH1D *h_Zg_2BidealAg = new TH1D("h_Zg_2BidealAg", "x=Zg 2b", zg_bins, zg_binsVector);
TH1D *h_Kt_2BidealAg = new TH1D("h_Kt_2BidealAg", "x=Kt 2b", kt_bins, kt_binsVector);
TH1D *h_Rg_2BidealAg = new TH1D("h_Rg_2BidealAg", "x=Rg 2b", logrg_bins, logrg_binsVector);

TH1D *h_Zg_2BrecoAg = new TH1D("h_Zg_2BrecoAg", "x=Zg 2b", zg_bins, zg_binsVector);
TH1D *h_Kt_2BrecoAg = new TH1D("h_Kt_2BrecoAg", "x=Kt 2b", kt_bins, kt_binsVector);
TH1D *h_Rg_2BrecoAg = new TH1D("h_Rg_2BrecoAg", "x=Rg 2b", logrg_bins, logrg_binsVector);

TH1D *h_Zg_recoAg = new TH1D("h_Zg_recoAg", "x=Zg 2b", zg_bins, zg_binsVector);
TH1D *h_Kt_recoAg = new TH1D("h_Kt_recoAg", "x=Kt 2b", kt_bins, kt_binsVector);
TH1D *h_Rg_recoAg = new TH1D("h_Rg_recoAg", "x=Rg 2b", logrg_bins, logrg_binsVector);