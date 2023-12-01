/**
 *  @file  ubcc1pi_standalone/Macros/MakeEventSelectionTable.cxx
 *
 *  @brief The implementation file of the MakeEventSelectionTable macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeEventSelectionTable(const Config &config)
{
    // Setup the input files
    // std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we define a vector of tuples with 4 entries
    //   - First, the sample type (e.g. overlay)
    //   - Second, the path to the input file
    //   - Third, the normalisation factor to apply to all events in that file

    // std::cout<<"##########################################\nUSING NUWRO AS DATA & Only CC0pi!\n##########################################"<<std::endl;
    // for (const auto fileRun: config.global.runs)
    // {
    //     if(fileRun == 1)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysfilePath, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 1));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun1.nuWrofilePath, 1.f);
    //     }
    //     else if(fileRun == 2)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysfilePath, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 2));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun2.nuWrofilePath, 1.f);
    //     }
    //     else if(fileRun == 3)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysfilePath, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 3));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun3.nuWrofilePath, 1.f);
    //     }
    //     else throw std::logic_error("PlotEventSelectionCuts - Invalid run number");
    // }

    // for (const auto fileRun: config.global.runs)
    // {
    //     if(fileRun == 1)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysfilePath, NormalisationHelper::GetOverlaysNormalisation(config, 1));
    //         inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun1.dirtfilePath, NormalisationHelper::GetDirtNormalisation(config, 1));
    //         inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun1.dataEXTfilePath, NormalisationHelper::GetDataEXTNormalisation(config, 1));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun1.dataBNBfilePath, 1.f);
    //     }
    //     else if(fileRun == 2)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysfilePath, NormalisationHelper::GetOverlaysNormalisation(config, 2));
    //         inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun2.dirtfilePath, NormalisationHelper::GetDirtNormalisation(config, 2));
    //         inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun2.dataEXTfilePath, NormalisationHelper::GetDataEXTNormalisation(config, 2));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun2.dataBNBfilePath, 1.f);
    //     }
    //     else if(fileRun == 3)
    //     {
    //         inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysfilePath, NormalisationHelper::GetOverlaysNormalisation(config, 3));
    //         inputData.emplace_back(AnalysisHelper::Dirt,    config.filesRun3.dirtfilePath, NormalisationHelper::GetDirtNormalisation(config, 3));
    //         inputData.emplace_back(AnalysisHelper::DataEXT, config.filesRun3.dataEXTfilePath, NormalisationHelper::GetDataEXTNormalisation(config, 3));
    //         inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun3.dataBNBfilePath, 1.f);
    //     }
    //     else throw std::logic_error("ExtractSidebandFit - Invalid run number");
    // }


    std::cout<<"DEBUG MakeEventSelectionTable -2"<<std::endl;
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get list of good runs
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::vector<int> goodRuns;
    std::ifstream file(config.global.goodRunListFile);
    int run;
    while (file >> run) {
        goodRuns.push_back(run);
    }

    std::cout<<"DEBUG MakeEventSelectionTable -1"<<std::endl;

    // Set up the selection
    auto selection = SelectionHelper::GetDefaultSelection2(true);

    // Setup the event counter
    AnalysisHelper::EventCounter counter;

    // Loop over the input files
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        std::cout<<"DEBUG MakeEventSelectionTable 0"<<std::endl;
        if(sampleType == AnalysisHelper::NuWro) continue;
        if(sampleType == AnalysisHelper::DetectorVariation) continue;
        if(!useThisFile) continue;
        std::cout << "Processing file - " << filePath << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);
        const auto isMC = (sampleType != AnalysisHelper::DataEXT && sampleType != AnalysisHelper::DataBNB);

        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        // if (isMC) readerPeLEE.EnableSystematicBranches(); // Todo: Is this correct/optimal?
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

        
    

        // Loop over the events
        std::cout << "### Only processing 10\% of events ###" << std::endl;
        for (unsigned int i = 0; i < nEvents/10; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);

            readerPeLEE.LoadEvent(i);

            const auto run = pEventPeLEE->metadata.run();
            const auto isGoodRun = (isDataBNB || isDataEXT) ? std::find(goodRuns.begin(), goodRuns.end(), run) != goodRuns.end() : true; // Apply good runs cuts to data
            // if(!isGoodRun) continue;
            if(!isGoodRun)
            {
                // std::cout << "\tDEBUG - fileName: " << fileName << std::endl;
                continue;
            }

            Event event(*pEventPeLEE, false);
            // Event event(*pEventPeLEE, isMC); // todo remove this and use next line
            // Event event(*pEventPeLEE);
            const auto pEvent = std::make_shared<Event>(event);

            // const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold); // todo remove this
            // if(!isTrueCC0Pi) continue;

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            // std::cout << weight << " " << AnalysisHelper::GetNominalEventWeight(pEvent) << " " << normalisation << std::endl;

            // For overlay & dirt samples only, count all events before any selection cuts
            // ATTN The EXT and BNB data has been filtered through the CC inclusive selection, so we can't get the count before any cuts were applied
            if (sampleType == AnalysisHelper::Overlay || sampleType == AnalysisHelper::Dirt)
            {
                counter.CountEvent("all", sampleType, pEvent, weight);
            }

            // Run the selection
            const auto &[passedAllCuts, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Check which of the cuts the event passed
            for (const auto &cutName : selection.GetCuts())
            {
                // If the event didn't pass this cut, then move on
                if (!SelectionHelper::IsCutPassed(cutsPassed, cutName))
                    break;

                // Count the event, and use the cut name as the "tag"
                counter.CountEvent(cutName, sampleType, pEvent, weight);
            }
        }
    }

    // Print the cuts that were used
    FormattingHelper::PrintLine();
    std::cout << "Cuts" << std::endl;
    FormattingHelper::PrintLine();

    FormattingHelper::Table table({"Cut", "", "Value"});
    for (const auto &cutName : selection.GetCuts())
    {
        table.AddEmptyRow();
        table.SetEntry("Cut", cutName);

        if (selection.CutHasValue(cutName))
            table.SetEntry("Value", selection.GetCutValue(cutName));
    }

    table.WriteToFile("eventSelection_cuts.md");

    // Print the breakdown of the event counts
    FormattingHelper::PrintLine();
    std::cout << "Summary" << std::endl;
    FormattingHelper::PrintLine();
    counter.PrintBreakdownSummary("eventSelection_summary.md");

    FormattingHelper::PrintLine();
    std::cout << "Details" << std::endl;
    FormattingHelper::PrintLine();

    const auto nEntriesToPrint = 10u;
    counter.PrintBreakdownDetails("eventSelection_details.md", nEntriesToPrint);
}

} // namespace ubcc1pi_macros
