/**
 *  @file  ubcc1pi_standalone/Macros/Macros.h
 *
 *  @brief The header file for the macros
 */

#ifndef UBCC1PI_STANDALONE_MACROS_MACROS
#define UBCC1PI_STANDALONE_MACROS_MACROS

#include "ubcc1pi_standalone/Objects/Config.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

    // /**
    //  *  @brief  Macro to print the configuration
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PrintConfig(const Config &config = Config());

    /**
     *  @brief  Process nTuples for unfolded analysis
     *
     *  @param  config the input configuration
     */
    void Analyzer(const Config &config = Config());

    /**
     *  @brief  Process nTuples for unfolded analysis
     *
     *  @param  config the input configuration
     */
    void AnalyzerTest(const Config &config = Config());

    /**
     *  @brief  Process nTuples for unfolded analysis
     *
     *  @param  config the input configuration
     */
    void AnalyzerTest2(const Config &config = Config());

    /**
     *  @brief  Check that the numucc filter and the CC inclusive cuts included in the analysis have the same effect
     *
     *  @param  config the input configuration
     */
    void TestCCInclusiveSelection(const Config &config = Config());

    /**
     *  @brief  Check that the selection is the same when PeLEE ntuples are used
     *
     *  @param  config the input configuration
     */
    void TestSelectionsPeLEEVsUbcc1pi(const Config &config = Config());

    /**
     *  @brief  Check that the sideband selection is the same when PeLEE ntuples are used
     *
     *  @param  config the input configuration
     */
    void TestSidebandSelectionsPeLEEVsUbcc1pi(const Config &config = Config());

    /**
     *  @brief  Count the POT for the overlays and dirt samples and print the results
     *
     *  @param  config the input configuration
     */
    void CountPOT(const Config &config = Config());

    /**
     *  @brief  Get the run-subrun list that's needed to run Zarco's POT / trigger / spill counting tool
     *
     *  @param  config the input configuration
     */
    void GetRunSubrunList(const Config &config = Config());

    // /**
    //  *  @brief  Look at the truth information related to signal events
    //  *
    //  *  @param  config the input configuration
    //  */
    // void TruthStudy(const Config &config = Config());

    /**
     *  @brief  Look at the truth information relating to secondary interactions
     *
     *  @param  config the input configuration
     */
    void SecondaryInteractionsStudy(const Config &config = Config());

    // /**
    //  *  @brief  Get the range->momentum curve fit parameters
    //  *
    //  *  @param  config the input configuration
    //  */
    // void FitRangeCurves(const Config &config = Config());

    // /**
    //  *  @brief  Determine the momentum thresholds to apply
    //  *
    //  *  @param  config the input configuration
    //  */
    // void MomentumThresholdsStudy(const Config &config = Config());

    /**
     *  @brief  Look at the accuracy of the muon PID from the CC inclusive as a function of the muon kinematics
     *
     *  @param  config the input configuration
     */
    void CCInclusiveMuonPIDStudy(const Config &config = Config());

    /**
     *  @brief  Demonstrate the impact of track direction on the PID, and how it can be improved
     *
     *  @param  config the input configuration
     */
    void MultiPlanePIDDemo(const Config &config = Config());

    /**
     *  @brief  Plot the variables that form the input to the particle ID BDTs after the CC inclusive selection
     *
     *  @param  config the input configuration
     */
    void PlotInputVariables(const Config &config = Config());

    // /**
    //  *  @brief  Plot the variables that form the input to the particle ID BDTs after the CC inclusive selection broken down for pions
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PlotPionInputVariables(const Config &config = Config());

    /**
     *  @brief  Get the correlation plots between the input variables to the BDT
     *
     *  @param  config the input configuration
     */
    void GetCorrelationPlots(const Config &config = Config());

    /**
     *  @brief  Train the BDTs
     *
     *  @param  config the input configuration
     */
    void TrainBDTs(const Config &config = Config());

    /**
     *  @brief  Run the N-1 BDT study
     *
     *  @param  config the input configuration
     */
    void NMinusOneBDTStudy(const Config &config = Config());


    /**
     *  @brief  Print the normalisations of the input files
     *
     *  @param  config the input configuration
     */
    void PrintFileNormalisations(const Config &config = Config());

    // /**
    //  *  @brief  Make data/MC plots for one iteration of the N-1 BDT study
    //  *
    //  *  @param  config the input configuration
    //  */
    // void NMinusOneBDTDataMCStudy(const Config &config = Config());

    // /**
    // *  @brief  Run multiple iterations of the N-1 BDT study following a hard-coded ordering
    // *
    // *  @param  config the input configuration
    // */
    // void NMinusOneBDTStudyFull(const Config &config = Config());

    // /**
    //  *  @brief  Study the impact of the muon PID
    //  *
    //  *  @param  config the input configuration
    //  */
    // void MuonPIDStudy(const Config &config = Config());

    /**
     *  @brief  Make the event selection table
     *
     *  @param  config the input configuration
     */
    void MakeEventSelectionTable(const Config &config = Config());

    /**
     *  @brief  Make the event selection table
     *
     *  @param  config the input configuration
     */
    void MakeEventSelectionTableRecycle(const Config &config = Config());

    /**
     *  @brief  Execute the unfolding steps in the xsec_analyser
     *
     *  @param  config the input configuration
     */
    void Unfold(const Config &config = Config());

    /**
     *  @brief  Make the sideband event selection table
     *
     *  @param  config the input configuration
     */
    void MakeSidebandEventSelectionTable(const Config &config = Config());

    /**
     *  @brief  Get the properties of pion candidates from CC0pi events in the CC1pi selection
     *
     *  @param  config the input configuration
     */
    void MakeCC0piProtonSelectionTable(const Config &config = Config());

    /**
     *  @brief  Make the PID table
     *
     *  @param  config the input configuration
     */
    void MakeSelectedPIDTable(const Config &config = Config());

    /**
     *  @brief  Make the PID table
     *
     *  @param  config the input configuration
     */
    void MakeSidebandSelectedPIDTable(const Config &config = Config());

    /**
     *  @brief  Make the selection efficiency plots
     *
     *  @param  config the input configuration
     */
    void MakeEventSelectionEfficiencyPlots(const Config &config = Config());

    /**
     *  @brief  Plot the variables that are used by the selection at the point of usage
     *
     *  @param  config the input configuration
     */
    void PlotEventSelectionCuts(const Config &config = Config());

    /**
     *  @brief  Plot the selected events after each cut
     *
     *  @param  config the input configuration
     */
    void PlotEventSelectionCutsAlternate(const Config &config = Config());

        /**
     *  @brief  Plot the variables that are used by the sideband selection at the point of usage
     *
     *  @param  config the input configuration
     */
    void PlotSidebandEventSelectionCuts(const Config &config = Config());

    // /**
    //  *  @brief  Plot the reconstructed variables for the muon candidate
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PlotMuonRecoVariables(const Config &config = Config());

    /**
     *  @brief  Plot the reconstructed variables for the highest-energy proton candidate
     *
     *  @param  config the input configuration
     */
    void PlotProtonVariables(const Config &config = Config());

    /**
     *  @brief  Querry .bin objects created during fit
     *
     *  @param  config the input configuration
     */
    void QuerryFitObject(const Config &config = Config());

    /**
     *  @brief  Querry weights of input files
     *
     *  @param  config the input configuration
     */
    void QuerryWeights(const Config &config = Config());

    // /**
    //  *  @brief  Plot stacked by interaction the true momentum for the highest-energy MC protons
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PlotProtonMomentumByInteraction(const Config &config = Config());

    // /**
    //  *  @brief  Plot wiggliness vs recoonstructed momentum
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PlotEBRequests(const Config &config = Config());

    // /**
    //  *  @brief  Plot the recontructed variables
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PlotReconstructedVariables(const Config &config = Config());

    // /**
    // *  @brief  Dump the truth and reconstructed information about selected events to a ROOT file
    // *
    // *  @param  config the input configuration
    // */
    // void DumpSelectedEventInfo(const Config &config = Config());

    // /**
    //  *  @brief  Make the plots to motivate the binning choice
    //  *
    //  *  @param  config the input configuration
    //  */
    // void MakeBinningPlots(const Config &config = Config());

    // /**
    //  *  @brief  Plot the flux distribution
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PlotFlux(const Config &config = Config());

    // /**
    //  *  @brief  Plot the variations of the flux for each systematic parameter
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PlotFluxVariations(const Config &config = Config());

    // /**
    //  *  @brief  Extract the cross-sections
    //  *
    //  *  @param  config the input configuration
    //  */
    // void ExtractXSecs(const Config &config = Config());

    /**
     *  @brief  Print universe weights for background CC0pi
     *
     *  @param  config the input configuration
     */
    void PrintUniverseWeights(const Config &config = Config());


    // /**
    //  *  @brief  Extract the Nuwro fake data cross-sections
    //  *
    //  *  @param  config the input configuration
    //  */
    // void ExtractNuWroXSecs(const Config &config = Config());

    /**
     *  @brief  Extract the Nuwro fake data cross-sections
     *
     *  @param  config the input configuration
     */
    void ExtractNuWroXSecs2(const Config &config = Config());

    // /**
    //  *  @brief  Extract the cross-sections
    //  *
    //  *  @param  config the input configuration
    //  */
    // void ExtractXSecsOld(const Config &config = Config());

    // /**
    //  *  @brief  Extract the sideband fit
    //  *
    //  *  @param  config the input configuration
    //  */
    // void ExtractSidebandFit(const Config &config = Config());

    /**
     *  @brief  Extract the sideband fit of Nuwro to MC (GENIE)
     *
     *  @param  config the input configuration
     */
    void ExtractNuWroSidebandFit(const Config &config = Config());

    /**
     *  @brief  Extract the sideband fit of Nuwro to MC (GENIE)
     *
     *  @param  config the input configuration
     */
    void ExtractNuWroSidebandFit2(const Config &config = Config());

    /**
     *  @brief  Extract the sideband fit of Nuwro to MC (GENIE)
     *
     *  @param  config the input configuration
     */
    void ExtractNuWroSidebandFit3(const Config &config = Config());

    /**
     *  @brief  Extract the sideband fit of Nuwro to MC (GENIE)
     *
     *  @param  config the input configuration
     */
    void MakeNuWroXSecTruthDistributionPlots(const Config &config = Config());

    /**
     *  @brief  Extract the sideband fit of Nuwro to MC (GENIE)  with all (BNB) systematic uncertainties
     *
     *  @param  config the input configuration
     */
    void ExtractNuWroAllSystSidebandFit(const Config &config = Config());

    /**
     *  @brief  Extract the Nuwro fake data cross-sections with all (BNB) systematic uncertainties
     *
     *  @param  config the input configuration
     */
    void ExtractNuWroAllSystXSecs(const Config &config = Config());

    /**
     *  @brief  Extract the sideband without fitting
     *
     *  @param  config the input configuration
     */
    //void ExtractSideband(const Config &config = Config());

    /**
     *  @brief  Print out detector variation xsec-values
     *
     *  @param  config the input configuration
     */
    void PrintDetVar(const Config &config = Config());

    /**
     *  @brief  Save the sideband fit parameters as a numpy readable txt file
     *
     *  @param  config the input configuration
     */
    void ParamsToTxt(const Config &config = Config());

    /**
     *  @brief  Compare xsec and sideband selection
     *
     *  @param  config the input configuration
     */
    void SelectionComparison(const Config &config = Config());

    // /**
    //  *  @brief  Generate the CC0Pi selection
    //  *
    //  *  @param  config the input configuration
    //  */
    // void ExtractCC0PiNormalisation(const Config &config = Config());

    // /**
    //  *  @brief  Fit sideband MC to data
    //  *
    //  *  @param  config the input configuration
    //  */
    // void MakeSidebandTemplateFit(const Config &config = Config());

    // /**
    //  *  @brief  Print a summary of the uncertainties for the total cross-section
    //  *
    //  *  @param  config the input configuration
    //  */
    // void PrintUncertaintiesSummary(const Config &config = Config());

    /**
    *  @brief  Produce the plots for the previously extracted cross-section data
    *
    *  @param  config the input configuration
    */
    void MakeXSecPlots(const Config &config = Config());

    /**
    *  @brief  Produce the plots for the previously extracted cross-section data
    *
    *  @param  config the input configuration
    */
    void MakeNuWroXSecPlots(const Config &config = Config());

    /**
    *  @brief  Produce the plots for the previously extracted sideband fit
    *
    *  @param  config the input configuration
    */
    void MakeSidebandFitPlots(const Config &config = Config());

    /*
     *  @brief  Make plots for the CC0pi1p sideband sample
     *
     *  @param  config the input configuration
     */
    void MakeSidebandSamplePlots(const Config &config = Config());

    /*
     *  @brief  Make plots for the CC0pi1p sideband sample
     *
     *  @param  config the input configuration
     */
    void MakeSidebandSamplePlotsV2(const Config &config = Config());

    /*
     *  @brief  Optimize the sideband cuts
     *
     *  @param  config the input configuration
     */
    void OptimizeSidebandCuts(const Config &config = Config());

    /*
     *  @brief  Confusion plots for muon and pion candidates
     *
     *  @param  config the input configuration
     */
    void PlotMuonConfusion(const Config &config = Config());

    // /*
    //  *  @brief  Optimize the proton momentum threshold
    //  *
    //  *  @param  config the input configuration
    //  */
    // void OptimizeProtonMomentumCut(const Config &config = Config());

    /*
     *  @brief  Make plots for the sideband parameters
     *
     *  @param  config the input configuration
     */
    void MakeSidebandParameterPlots(const Config &config = Config());

    /**
     *  @brief  Optimize the sideband truth bin edges
     *
     *  @param  config the input config
     */
    void OptimizeSidebandBinEdges(const Config &config = Config());

    /**
     *  @brief  Check if the pion selected proton in CC0pi selection is also the true leading proton
     *
     *  @param  config the input config
     */
    void CheckCC0piProtonMomentum(const Config &config = Config());

    /**
     *  @brief  Top-level macro that runs the full analysis chain
     *
     *  @param  config the input config
     */
    void RunFullAnalysis(const Config &config = Config());

} // namespace ubcc1pi_macros

#endif
