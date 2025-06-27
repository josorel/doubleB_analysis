Interested in the **Jet Substructure of Unresolved Gluon to b-b Jets** ?
-

This code uses ROOT C++. Packages to install may include RooUnfold, RooFit. 
In this repository, you will find important files to conduct the analysis, such as:

- `run_MC_inclusive_response.c` - creates MC distributions and response matrix for inclusive jets
- `run_MC_response.c` - creates MC distributions and response matrix for double B jets
- `run_data_inclusive.c` - creates data distributions for inclusive jets
- `run_data.c` - creates data distributions for double B jets
- `run_templateFit.c` - performs template fit on data using MC distributions, generates signal fractions
- `run_results.c` - applies signal fractions from template fit and unfolding/corrections, then can display various results
<br />

- `draw_response.C` - Draws response matrices as well as efficiency and purity hists
- `draw_uncertainties.c` - Can be used to draw systematic uncertainties
<br />

- `binning.h` - has bins for substructure variables as well as histograms used in running other files
- `tTree.h` - has template for tree generated in files
<br />

The following files are used for BDT purposes:
- `prepare_data.C`
- `test_TMVA.C`
- `train_TMVA.C`
- `draw_ROC_TMVA.C`

