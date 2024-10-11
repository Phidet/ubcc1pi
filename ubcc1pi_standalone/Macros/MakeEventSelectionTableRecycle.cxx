/**
 *  @file  ubcc1pi_standalone/Macros/MakeEventSelectionTableRecycle.cxx
 *
 *  @brief The implementation file of the MakeEventSelectionTableRecycle macro
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

void MakeEventSelectionTableRecycle(const Config &config)
{

    std::cout<<"DEBUG MakeEventSelectionTableRecycle -2"<<std::endl;
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get list of good runs
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::vector<int> goodRuns;
    std::ifstream file(config.global.goodRunListFile);
    int run;
    while (file >> run) {
        goodRuns.push_back(run);
    }

    std::cout<<"DEBUG MakeEventSelectionTableRecycle -1"<<std::endl;

    // Set up the selection
    // auto selection = SelectionHelper::GetDefaultSelection2(true);

    // Setup the event counter
    AnalysisHelper::EventCounter counter;

    // Check if the sizes of input and output files are identical
    if (config.input.files.size() != config.output.files.size())
        throw std::runtime_error("Input and output file sizes are not identical");

    // Loop over the input files
    for (unsigned int f = 0; f < config.input.files.size(); ++f)
    {
        const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] = config.input.files.at(f);
        const auto &[fileRunOut, normalisationOut, sampleTypeOut, useThisFileOut, filePathOut] = config.output.files.at(f);

        // Check that the input and output sampleTypes and runs are the same
        if (fileRun != fileRunOut)
            throw std::runtime_error("Input and output file runs are not identical");
        if (sampleType != sampleTypeOut)
            throw std::runtime_error("Input and output file sample types are not identical");
        if(useThisFile!=useThisFileOut)
            throw std::runtime_error("Input and output file 'useThisFile' flag are not identical");

        std::cout<<"DEBUG MakeEventSelectionTableRecycle 0"<<std::endl;
        if(sampleType == AnalysisHelper::NuWro) continue;
        if(sampleType == AnalysisHelper::DetectorVariation) continue;
        if(!useThisFile) continue;
        std::cout << "Processing file - " << filePath << " and " << filePathOut << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);
        const auto isMC = (!isDataEXT && !isDataBNB);


        std::cout<<"DEBUG MakeEventSelectionTableRecycle 1: "<<filePath<<" isMC: "<<isMC<<std::endl;
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        if(sampleType == isOverlay) readerPeLEE.EnableSystematicBranches();
        std::cout<<"DEBUG MakeEventSelectionTableRecycle 1.1"<<std::endl;
        std::cout<<"DEBUG MakeEventSelectionTableRecycle 1.2"<<std::endl;
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        std::cout<<"DEBUG MakeEventSelectionTableRecycle 2"<<std::endl;

        FileReader<EventXSec, SubrunXSec> readerXSec(filePathOut, isMC);
        // if (isMC) readerXSec.EnableSystematicBranches(); // Todo: Is this correct/optimal?
        const auto pEventXSec = readerXSec.GetBoundEventAddress();
        std::cout<<"DEBUG MakeEventSelectionTableRecycle 3"<<std::endl;


        // Loop over the events
        unsigned int indexPeLEE = 0;
        unsigned int indexXSec = 0;
        // for (unsigned int i = 0; i < nEvents/100; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // std::cout << "### Only processing 1\% of events ###" << std::endl;
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto nEventsXSec = readerXSec.GetNumberOfEvents();
        while(indexPeLEE<nEvents)
        {
            // AnalysisHelper::PrintLoadingBar(indexPeLEE, nEvents);
            if(indexPeLEE%1000==0) std::cout<<"\rProgress: "<<indexPeLEE<<"/"<<nEvents<<" "<<indexXSec<<"/"<<nEventsXSec<<std::flush;

            readerPeLEE.LoadEvent(indexPeLEE);

            auto isGoodRun = true;
            auto run = -1;
            if(isDataBNB || isDataEXT)
            {
                run = pEventPeLEE->metadata.run();
                isGoodRun = std::find(goodRuns.begin(), goodRuns.end(), run) != goodRuns.end(); // Apply good runs cuts to data
            }

            if(!isGoodRun)
            {
                // std::cout << "Skipping bad run: " << run << std::endl;
                indexPeLEE++;
                continue;
            }

            readerXSec.LoadEvent(indexXSec);
            Event event(*pEventPeLEE, false);
            // Event event(*pEventPeLEE, isMC); // todo remove this and use next line
            // Event event(*pEventPeLEE);
            const auto pEvent = std::make_shared<Event>(event);
            
            indexPeLEE++;
            indexXSec++;

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


            if(pEventXSec->reco.passed_particleTrackScore()) counter.CountEvent("particleTrackScore", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_particleVertexDistance()) counter.CountEvent("particleVertexDistance", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_particleGeneration()) counter.CountEvent("particleGeneration", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_particleTrackLength()) counter.CountEvent("particleTrackLength", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_particleProtonChi2()) counter.CountEvent("particleProtonChi2", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_particleMuonChi2()) counter.CountEvent("particleMuonChi2", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_particleProtonChi2OverMuonChi2()) counter.CountEvent("particleProtonChi2OverMuonChi2", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_pandoraNuPDGIsNumu()) counter.CountEvent("pandoraNuPDGIsNumu", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_daughterVerticesContained()) counter.CountEvent("daughterVerticesContained", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_nuVertexFiducial()) counter.CountEvent("nuVertexFiducial", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_topologicalOrFlashMatch()) counter.CountEvent("topologicalOrFlashMatch", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_topologicalScoreCC()) counter.CountEvent("topologicalScoreCC", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_min2Tracks()) counter.CountEvent("min2Tracks", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_max1Uncontained()) counter.CountEvent("max1Uncontained", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_2NonProtons()) counter.CountEvent("2NonProtons", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_pionHasValiddEdx()) counter.CountEvent("pionHasValiddEdx", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_pionNotInGap()) counter.CountEvent("pionNotInGap", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_muonNotInGap()) counter.CountEvent("muonNotInGap", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_openingAngle()) counter.CountEvent("openingAngle", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_topologicalScore()) counter.CountEvent("topologicalScore", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_startNearVertex()) counter.CountEvent("startNearVertex", sampleType, pEvent, weight); else continue;
            if(pEventXSec->reco.passed_likelyGoldenPion()) counter.CountEvent("likelyGoldenPion", sampleType, pEvent, weight); else continue;

            // // Run the selection
            // const auto &[passedAllCuts, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // // Check which of the cuts the event passed
            // for (const auto &cutName : selection.GetCuts())
            // {
            //     // If the event didn't pass this cut, then move on
            //     if (!SelectionHelper::IsCutPassed(cutsPassed, cutName))
            //         break;

            //     // Count the event, and use the cut name as the "tag"
            //     counter.CountEvent(cutName, sampleType, pEvent, weight);
            // }
        }
        // if(indexPeLEE!=nEvents || indexXSec!=nEventsXSec)
        //     throw std::runtime_error("Issue with skipping bad runs; indexPeLEE: "+std::to_string(indexPeLEE)+" nEvents: "+std::to_string(nEvents)+" indexXSec: "+std::to_string(indexXSec)+" nEventsXSec: "+std::to_string(nEventsXSec)); 
    }

    std::vector<std::string> cutNames = {
        "particleTrackScore",
        "particleVertexDistance",
        "particleGeneration",
        "particleTrackLength",
        "particleProtonChi2",
        "particleMuonChi2",
        "particleProtonChi2OverMuonChi2",
        "pandoraNuPDGIsNumu",
        "daughterVerticesContained",
        "nuVertexFiducial",
        "topologicalOrFlashMatch",
        "topologicalScoreCC",
        "min2Tracks",
        "max1Uncontained",
        "2NonProtons",
        "pionHasValiddEdx",
        "pionNotInGap",
        "muonNotInGap",
        "openingAngle",
        "topologicalScore",
        "startNearVertex",
        "likelyGoldenPion"
    };

    // Print the cuts that were used
    FormattingHelper::PrintLine();
    std::cout << "Cuts" << std::endl;
    FormattingHelper::PrintLine();

    FormattingHelper::Table table({"Cut", "", "Value"});
    for (const auto &cutName : cutNames)
    {
        table.AddEmptyRow();
        table.SetEntry("Cut", cutName);

        // if (selection.CutHasValue(cutName))
        //     table.SetEntry("Value", selection.GetCutValue(cutName));
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
