#include <cstdlib>
#include <stdexcept>
#include <string>

#include <vector>
#include <map>

{
    // Add boost library support
    auto boost_inc = gSystem->Getenv("BOOST_INC");
    auto boost_lib = gSystem->Getenv("BOOST_LIB");
    //gSystem->AddIncludePath(("-I " + std::string(boost_inc) + "/").c_str());
    gInterpreter->AddIncludePath((std::string(boost_inc) + "/boost/archive/").c_str());
    gInterpreter->AddIncludePath((std::string(boost_inc) + "/boost/serialization/").c_str());
    //gSystem->AddIncludePath(("  -I " + std::string(boost_inc) + "/boost/archive/").c_str());
    //gSystem->AddIncludePath(("  -I " + std::string(boost_inc) + "/boost/serialization/").c_str());

    gSystem->Load((std::string(boost_lib) + "/libboost_filesystem.so").c_str());
    gSystem->Load((std::string(boost_lib) + "/libboost_serialization.so").c_str());
    gSystem->Load((std::string(boost_lib) + "/libboost_wserialization.so").c_str());

    // Generate any required dictionaries
    gInterpreter->GenerateDictionary("std::map<std::string, std::vector<float> >", "map;string;vector");

    // Add ubcc1pi_standalone to the include path
    const auto ubcc1piStandaloneDir = std::getenv("UBCC1PI_STANDALONE_DIR");
    if (!ubcc1piStandaloneDir)
        throw std::logic_error("The environment variable UBCC1PI_STANDALONE_DIR hasn't been set!");

    std::cout << "Adding " << ubcc1piStandaloneDir << " to include path" << std::endl;
    gSystem->AddIncludePath(("  -I" + std::string(ubcc1piStandaloneDir) + "/../").c_str());

    std::cout << "Adding " << ubcc1piStandaloneDir << "/ubsmear/inc/ to include path" << std::endl;
    gSystem->AddIncludePath(("  -I" + std::string(ubcc1piStandaloneDir) + "/ubsmear/inc/").c_str());


    // Compile the code using root
    for (const auto &file : std::vector<std::string>({

        // Interface
        "Interface/Event.cxx",
        "Interface/EventPeLEE.cxx",
        "Interface/EventXSec.cxx",
        "Interface/Subrun.cxx",
        "Interface/SubrunPeLEE.cxx",
        "Interface/SubrunXSec.cxx",

        // Objects
        "Objects/FileReader.cxx",
        "Objects/TreeWriter.cxx",

        // Helpers
        "Helpers/AnalysisHelper.cxx",
        "Helpers/BDTHelper.cxx",
        "Helpers/PlottingHelper.cxx",
        "Helpers/SelectionHelper.cxx",
        "Helpers/NormalisationHelper.cxx",
        "Helpers/ExtractionHelper.cxx",
        "Helpers/FittingHelper.cxx",
        "Helpers/CrossSectionHelper.cxx",

        // Macros
        "Macros/Analyzer.cxx",
        // "Macros/AnalyzerTest.cxx",
        "Macros/AnalyzerTest2.cxx",
        // "Macros/TestCCInclusiveSelection.cxx",
        // "Macros/TestSelectionsPeLEEVsUbcc1pi.cxx",
        // "Macros/TestSidebandSelectionsPeLEEVsUbcc1pi.cxx",
        // "Macros/PrintConfig.cxx",
        "Macros/CountPOT.cxx",
        "Macros/GetRunSubrunList.cxx",
        "Macros/PrintFileNormalisations.cxx",
        // "Macros/TruthStudy.cxx",
        // "Macros/FitRangeCurves.cxx",
        // "Macros/SecondaryInteractionsStudy.cxx",
        // "Macros/MomentumThresholdsStudy.cxx",
        // "Macros/CCInclusiveMuonPIDStudy.cxx",
        // "Macros/MultiPlanePIDDemo.cxx",
        // "Macros/PlotPionInputVariables.cxx",
        "Macros/PlotInputVariables.cxx",
        // "Macros/GetCorrelationPlots.cxx",
        "Macros/TrainBDTs.cxx",
        // // "Macros/NMinusOneBDTStudy.cxx",
        // // "Macros/NMinusOneBDTDataMCStudy.cxx",
        // // "Macros/NMinusOneBDTStudyFull.cxx",
        // // "Macros/MuonPIDStudy.cxx",
        "Macros/MakeEventSelectionTable.cxx",
        "Macros/MakeEventSelectionTableRecycle.cxx",
        "Macros/MakeSidebandEventSelectionTable.cxx",
        "Macros/MakeCC0piProtonSelectionTable.cxx",
        // "Macros/MakeEventSelectionEfficiencyPlots.cxx",
        // "Macros/PlotEventSelectionCutsAlternate.cxx",
        "Macros/PlotEventSelectionCuts.cxx",
        "Macros/Unfold.cxx",
        "Macros/PlotSidebandEventSelectionCuts.cxx",
        // // "Macros/MakeSelectedPIDTable.cxx",
        // // "Macros/PlotMuonRecoVariables.cxx",
        // // "Macros/PlotReconstructedVariables.cxx",
        // // "Macros/PlotProtonVariables.cxx",
        "Macros/OptimizeSidebandCuts.cxx",
        "Macros/PlotMuonConfusion.cxx",
        // // "Macros/OptimizeProtonMomentumCut.cxx",
        // // "Macros/PlotProtonMomentumByInteraction.cxx",
        // // "Macros/DumpSelectedEventInfo.cxx",
        // // "Macros/MakeBinningPlots.cxx",
        // // "Macros/PlotFlux.cxx",
        // // "Macros/PlotFluxVariations.cxx",
        // // "Macros/ExtractXSecs.cxx",
        // // "Macros/ExtractNuWroXSecs.cxx",
        // "Macros/ExtractNuWroXSecs2.cxx",
        // // "Macros/ExtractXSecsOld.cxx",
        // // "Macros/ExtractSidebandFit.cxx",
        // // "Macros/ExtractNuWroSidebandFit.cxx",
        // // "Macros/ExtractNuWroSidebandFit2.cxx",
        // "Macros/ExtractNuWroSidebandFit3.cxx",
        // "Macros/ExtractNuWroAllSystSidebandFit.cxx",
        // "Macros/ExtractNuWroAllSystXSecs.cxx",
        // "Macros/MakeNuWroXSecTruthDistributionPlots.cxx",
        // //"Macros/ExtractSideband.cxx",
        // // "Macros/PrintUniverseWeights.cxx",
        // // "Macros/PrintUncertaintiesSummary.cxx",
        // // "Macros/MakeXSecPlots.cxx",
        // "Macros/MakeNuWroXSecPlots.cxx",
        // "Macros/MakeSidebandFitPlots.cxx",
        "Macros/MakeSidebandSamplePlots.cxx",
        "Macros/MakeSidebandSamplePlotsV2.cxx",
        // "Macros/MakeSidebandParameterPlots.cxx",
        // // "Macros/PrintDetVar.cxx",
        // // "Macros/ParamsToTxt.cxx",
        // // "Macros/SelectionComparison.cxx",
        // // "Macros/MakeSidebandSelectedPIDTable.cxx",
        // // "Macros/QuerryFitObject.cxx",
        // // "Macros/QuerryWeights.cxx",
        // // "Macros/CheckCC0piProtonMomentum.cxx",
        // "Macros/PlotEBRequests.cxx",
        // "Macros/RunFullAnalysis.cxx",
        // "Macros/OptimizeSidebandBinEdges.cxx"
    }))
    {
        std::cout << "Compiling " << file << std::endl;
        gROOT->ProcessLine((".L " + std::string(ubcc1piStandaloneDir) + "/" + file + "+").c_str());
    }

    // Setting the global plotting environment
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetTitleFont(133, "XYZ");
    gStyle->SetTitleSize(18, "XYZ");
    gStyle->SetTitleFontSize(18);

    gStyle->SetLabelFont(133, "XYZ");
    gStyle->SetLabelSize(28, "XYZ");

    TGaxis::SetMaxDigits(3);

    gROOT->ForceStyle();
}