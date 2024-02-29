/**
 *  @file  ubcc1pi_standalone/Macros/MakeSidebandEventSelectionTable.cxx
 *
 *  @brief The implementation file of the MakeSidebandEventSelectionTable macro
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

void MakeSidebandEventSelectionTable(const Config &config)
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get list of good runs
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::vector<int> goodRuns;
    std::ifstream file(config.global.goodRunListFile);
    int run;
    while (file >> run) {
        goodRuns.push_back(run);
    }

    // Set up the selection
    // auto selection = SelectionHelper::GetCC0piSelectionModifiedPeLEE(-0.48f, 0.12f, true);
    auto selection = SelectionHelper::GetCC0piSelectionModifiedPeLEE(-0.45f, 0.12f, true);
    // auto selection = SelectionHelper::GetCC0piSelectionModifiedAgainPeLEE(-0.45f, 0.12f, true);
    // auto selection = SelectionHelper::GetCC0piSelectionOriginalPeLEE(-0.5f, 0.2f, false);

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
        std::cout << "### Only processing 5\% of events ###" << std::endl;
        for (unsigned int i = 0; i < nEvents/20; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            const auto pEvent = std::make_shared<Event>(event);

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;


            // For overlay & dirt samples only, count all events before any selection cuts
            // ATTN The EXT and BNB data has been filtered through the CC inclusive selection, so we can't get the count before any cuts were applied
            if (sampleType == AnalysisHelper::Overlay || sampleType == AnalysisHelper::Dirt)
            {
                counter.CountEvent("all", sampleType, pEvent, weight, true, config.global.protonMomentumThreshold);
            }

            // Run the selection
            const auto &[passedAllCuts, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            if(assignedPdgCodes.size() != pEvent->reco.particles.size())
                throw std::logic_error("Assigned pdg codes vector has different size to the number of reconstructed particles");


            // Check which of the cuts the event passed
            for (const auto &cutName : selection.GetCuts())
            {
                // If the event didn't pass this cut, then move on
                if (!SelectionHelper::IsCutPassed(cutsPassed, cutName))
                    break;

                // Count the event, and use the cut name as the "tag"
                counter.CountEvent(cutName, sampleType, pEvent, weight, true, config.global.protonMomentumThreshold);
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

    table.WriteToFile("eventSelection_cuts_sidband.md");

    // Print the breakdown of the event counts
    FormattingHelper::PrintLine();
    std::cout << "Summary" << std::endl;
    FormattingHelper::PrintLine();
    counter.PrintBreakdownSummary("eventSelection_summary_sidband.md");

    FormattingHelper::PrintLine();
    std::cout << "Details" << std::endl;
    FormattingHelper::PrintLine();

    const auto nEntriesToPrint = 10u;
    counter.PrintBreakdownDetails("eventSelection_details_sidband.md", nEntriesToPrint);
}

} // namespace ubcc1pi_macros