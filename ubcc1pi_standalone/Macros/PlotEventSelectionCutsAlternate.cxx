/**
 *  @file  ubcc1pi_standalone/Macros/PlotEventSelectionCutsAlternate.cxx
 *
 *  @brief The implementation file of the PlotEventSelectionCutsAlternate macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotEventSelectionCutsAlternate(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    for (const auto run: config.global.runs)
    {
        if(run == 1)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 1));
        }
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 2));
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config, 3));
        }
        else throw std::logic_error("PlotEventSelectionCutsAlternate - Invalid run number");
    }

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetSelection(config.global.selection);
    const auto allCuts = selection.GetCuts();

    const std::string yLabelParticles = "Number of particles";
    const std::string yLabelEvents = "Number of events";

    PlottingHelper::MultiPlot passesCCInclusivePlot("passesCCInclusive", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot min2TracksPlot("min2Tracks", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot max1UncontainedPlot("max1Uncontained", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot twoNonProtonsPlot("twoNonProtons", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot pionHasValiddEdxPlot("pionHasValiddEdx", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot pionNotInGapPlot("pionNotInGap", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot muonNotInGapPlot("muonNotInGap", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot openingAnglePlot("openingAngle", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot topologicalScorePlot("topologicalScore", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    PlottingHelper::MultiPlot startNearVertexPlot("startNearVertex", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));
    // PlottingHelper::MultiPlot likelyGoldenPionPlot("likelyGoldenPion", yLabelEvents, PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f));


    // Define a function scoped to this macro to determine if a cut was passed
    auto passedCut = [&](const std::string &cut, const std::vector<std::string> &cutsPassed) -> bool {
        return (std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end());
    };

    // Define a function scoped to this macro to determine if this event formed the input to a given cut
    auto wasInputToCut = [&](const std::string &cut, const std::vector<std::string> &cutsPassed) -> bool {
        const auto &cutIter = std::find(allCuts.begin(), allCuts.end(), cut);
        if (cutIter == allCuts.end())
            throw std::invalid_argument("PlotEventSelecitonCuts - Unknown cut: " + cut);

        const auto cutIndex = std::distance(allCuts.begin(), cutIter);
        if (cutIndex == 0)
            throw std::invalid_argument("PlotEventSelecitonCuts - Cut: " + cut + " has no preceeding cut");

        // If we passed the preceeding cut, then this event was and input to the supplied cut
        const auto &preceedingCut = allCuts.at(cutIndex - 1);

        return passedCut(preceedingCut, cutsPassed);
    };

    // Loop over the events
    for (const auto [sampleType, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        FileReader reader(fileName);
        auto pEvent = reader.GetBoundEventAddress();

        const auto nEvents = reader.GetNumberOfEvents();

        std::cout<<"############################\nUsing "<<20<<"\% of events!\n############################"<<std::endl;
        for (unsigned int i = 0; i < nEvents/5; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // For brevity assign the particles a variable
            const auto &truthParticles = pEvent->truth.particles;
            const auto &recoParticles = pEvent->reco.particles;

            // Run the event selection and store which cuts are passed
            const auto &[passesGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Get the plot style of the event
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyleEvent = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            const auto isTrueCC1Pi = AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
            if(!isTrueCC1Pi) continue;

            const auto truthData = (
                isTrueCC1Pi
                    ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            if (wasInputToCut("min2Tracks", cutsPassed))
            {
                passesCCInclusivePlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("max1Uncontained", cutsPassed))
            {
                min2TracksPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("2NonProtons", cutsPassed))
            {
                max1UncontainedPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("pionHasValiddEdx", cutsPassed))
            {
                twoNonProtonsPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if(wasInputToCut("pionNotInGap", cutsPassed))
            {
                pionHasValiddEdxPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("muonNotInGap", cutsPassed))
            {
                pionNotInGapPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("openingAngle", cutsPassed)){
                muonNotInGapPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("topologicalScore", cutsPassed))
            {
                openingAnglePlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("startNearVertex", cutsPassed))
            {
                topologicalScorePlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            if (wasInputToCut("likelyGoldenPion", cutsPassed))
            {
                startNearVertexPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            }

            // if (wasInputToCut("", cutsPassed))
            // {
            //     likelyGoldenPionPlot.Fill(truthData.muonPhi, plotStyleEvent, weight);
            // }
        }
    }

    passesCCInclusivePlot.SaveAsStacked("1-plotEventSelectionCutsAlt_passesCCInclusive");
    min2TracksPlot.SaveAsStacked("2-plotEventSelectionCutsAlt_min2Tracks");
    max1UncontainedPlot.SaveAsStacked("3-plotEventSelectionCutsAlt_max1Uncontained");
    twoNonProtonsPlot.SaveAsStacked("4-plotEventSelectionCutsAlt_twoNonProtons");
    pionHasValiddEdxPlot.SaveAsStacked("5-plotEventSelectionCutsAlt_pionHasValiddEdx");
    pionNotInGapPlot.SaveAsStacked("6-plotEventSelectionCutsAlt_pionNotInGap");
    muonNotInGapPlot.SaveAsStacked("7-plotEventSelectionCutsAlt_muonNotInGap");
    openingAnglePlot.SaveAsStacked("8-plotEventSelectionCutsAlt_openingAngle");
    topologicalScorePlot.SaveAsStacked("9-plotEventSelectionCutsAlt_topologicalScore");
    startNearVertexPlot.SaveAsStacked("10-plotEventSelectionCutsAlt_startNearVertex");
    // likelyGoldenPionPlot.SaveAsStacked("11-plotEventSelectionCutsAlt_likelyGoldenPion");
}

} // namespace ubcc1pi_macros
