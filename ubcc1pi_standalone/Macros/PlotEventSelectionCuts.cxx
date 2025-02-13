/**
 *  @file  ubcc1pi_standalone/Macros/PlotEventSelectionCuts.cxx
 *
 *  @brief The implementation file of the PlotEventSelectionCuts macro
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

void PlotEventSelectionCuts(const Config &config)
{
    //
    // Setup the input files
    //
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData;

    // inputData.emplace_back(AnalysisHelper::Overlay, config.files.overlaysFileName, NormalisationHelper::GetOverlaysNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::Dirt,    config.files.dirtFileName,     NormalisationHelper::GetDirtNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::DataEXT, config.files.dataEXTFileName,  NormalisationHelper::GetDataEXTNormalisation(config));
    // inputData.emplace_back(AnalysisHelper::DataBNB, config.files.dataBNBFileName,  1.f);


std::cout<<"##########################################\nUSING NUWRO AS DATA & Only CC0pi!\n##########################################"<<std::endl;
 for (const auto run: config.global.runs)
    {
        if(run == 1)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun1.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 1));
            inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun1.nuWroFileName, 1.f);
        }
        else if(run == 2)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun2.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 2));
            inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun2.nuWroFileName, 1.f);
        }
        else if(run == 3)
        {
            inputData.emplace_back(AnalysisHelper::Overlay, config.filesRun3.overlaysFileName, NormalisationHelper::GetOverlaysNormalisationToNuWro(config, 3));
            inputData.emplace_back(AnalysisHelper::DataBNB, config.filesRun3.nuWroFileName, 1.f);
        }
        else throw std::logic_error("PlotEventSelectionCuts - Invalid run number");
    }

    //
    // Get the selection
    //
    auto selection = SelectionHelper::GetSelection(config.global.selection);
    const auto allCuts = selection.GetCuts();

    const auto protonBDTCut = selection.GetCutValue("2NonProtons");

    //
    // Setup the plots
    //

    const std::string yLabelParticles = "Number of particles";
    const std::string yLabelEvents = "Number of events";

    PlottingHelper::MultiPlot nTracksPlot("Number of tracks", yLabelEvents, 6u, 1, 7);
    PlottingHelper::MultiPlot isContainedPlot("Is contained?", yLabelEvents, 2u, 0, 2);
    PlottingHelper::MultiPlot uncontainedVertexDistPlot("Distance to vertex / cm", yLabelParticles, PlottingHelper::GenerateLogBinEdges(40u, 0.15f, 1000.f));
    PlottingHelper::MultiPlot nUncontainedPlot("Number of uncontained particles", yLabelEvents, 4u, 0, 4);
    PlottingHelper::MultiPlot protonBDTResponsePlot("Proton BDT response", yLabelParticles, 40u, -0.60f, 0.60f);
    PlottingHelper::MultiPlot recoProtonMomentumPlot("Reco Proton momentum", yLabelParticles, 40u, -0.f, 2.5f);
    PlottingHelper::MultiPlot trueProtonMomentumPlot("True Proton momentum", yLabelParticles, 40u, -0.f, 2.5f);
    PlottingHelper::MultiPlot nNonProtonsPlot("Number of non-protons", yLabelEvents, 6u, 0, 6);
    PlottingHelper::MultiPlot truncatedMeandEdxBeforePlot("Pion cos(theta) / rad", yLabelEvents, 40u, -1.f, 1.f);
    PlottingHelper::MultiPlot truncatedMeandEdxAfterPlot("Pion cos(theta) / rad", yLabelEvents, 40u, -1.f, 1.f);
    PlottingHelper::MultiPlot pionNotInGapBeforePlot("Pion phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot pionNotInGapAfterPlot("Pion phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot openingAnglePlot("Muon-Pion opening angle / rad", yLabelEvents, 40u, 0.f, 3.142f);
    PlottingHelper::MultiPlot topologicalScorePlot("TopologicalScore", yLabelEvents, 40u, 0.06f, 1.f);
    PlottingHelper::MultiPlot startNearVertexParticlePlot("Distance to vertex / cm", yLabelParticles, PlottingHelper::GenerateLogBinEdges(40u, 0.03f, 1000.f));
    PlottingHelper::MultiPlot startNearVertexEventPlot("Min distance to vertex / cm", yLabelEvents, PlottingHelper::GenerateLogBinEdges(40u, 0.15f, 1000.f));
    PlottingHelper::MultiPlot likelyGoldenPionParticlePlot("Golden pion BDT response", yLabelParticles, 40u, -0.55f, 0.4f);
    PlottingHelper::MultiPlot likelyGoldenPionEventPlot("Golden pion BDT response", yLabelEvents, 40u, -0.55f, 0.4f);

    // Set the bin labels where appropriate
    nTracksPlot.SetIntegerBinLabels();
    nUncontainedPlot.SetIntegerBinLabels();
    nNonProtonsPlot.SetIntegerBinLabels();
    isContainedPlot.SetBinLabels({"No", "Yes"});

    // Add the cut values
    nTracksPlot.AddCutLine(2);
    nUncontainedPlot.AddCutLine(2);
    protonBDTResponsePlot.AddCutLine(protonBDTCut);
    nNonProtonsPlot.AddCutLine(2);
    nNonProtonsPlot.AddCutLine(3);
    openingAnglePlot.AddCutLine(selection.GetCutValue("openingAngle"));
    topologicalScorePlot.AddCutLine(selection.GetCutValue("topologicalScore"));
    startNearVertexParticlePlot.AddCutLine(selection.GetCutValue("startNearVertex"));
    startNearVertexEventPlot.AddCutLine(selection.GetCutValue("startNearVertex"));
    likelyGoldenPionParticlePlot.AddCutLine(selection.GetCutValue("likelyGoldenPion"));
    likelyGoldenPionEventPlot.AddCutLine(selection.GetCutValue("likelyGoldenPion"));

    //
    // Setup the BDTs
    //
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;

    BDTHelper::BDT goldenPionBDT("goldenPion", goldenPionFeatureNames);
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

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
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            reader.LoadEvent(i);

            // For brevity assign the particles a variable
            const auto &truthParticles = pEvent->truth.particles;
            const auto &recoParticles = pEvent->reco.particles;

            // Run the event selection and store which cuts are passed
            const auto &[passesGoldenPionSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Get the plot style of the event
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyleEvent = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // if(PlottingHelper::PlotStyle::NumuCC0Pi != plotStyleEvent) continue;
            const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold); // todo remove this
            if(!isTrueCC0Pi) continue;

            // std::cout<<"truthParticles:"<<std::endl;
            // for(const auto &particle: truthParticles)
            //     std::cout<<particle.pdgCode()<<" "<<std::endl;
            // std::cout<<"recoParticles:"<<std::endl;
            // for(const auto &particle: recoParticles)
            //     std::cout<<particle.pdgCode()<<" "<<std::endl;

            // Fill the plots
            // ...

            if (wasInputToCut("min2Tracks", cutsPassed))
            {
                unsigned int nTracks = 0u;
                for (const auto &recoParticle : recoParticles)
                {
                    if (AnalysisHelper::HasTrackFit(recoParticle))
                        nTracks++;
                }

                nTracksPlot.Fill(static_cast<float>(nTracks), plotStyleEvent, weight);
            }

            if (wasInputToCut("max1Uncontained", cutsPassed))
            {
                unsigned int nUncontained = 0u;
                for (const auto &recoParticle : recoParticles)
                {
                    if (!AnalysisHelper::HasTrackFit(recoParticle))
                        continue;

                    const auto isContained = AnalysisHelper::IsContained(recoParticle);
                    nUncontained += isContained ? 0u : 1u;

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(recoParticle, sampleType, truthParticles, false, config.global.useAbsPdg);
                    isContainedPlot.Fill(isContained ? 1.f : 0.f, plotStyleParticle, weight);

                    if (!isContained)
                    {
                        const TVector3 start(recoParticle.startX(), recoParticle.startY(), recoParticle.startZ());
                        const auto &recoVertex = pEvent->reco.nuVertex();
                        const float vertexDist2 = (start - recoVertex).Mag2();
                        const auto vertexDist = std::pow(vertexDist2, 0.5f);

                        uncontainedVertexDistPlot.Fill(vertexDist, plotStyleParticle, weight);
                    }
                }

                nUncontainedPlot.Fill(static_cast<float>(nUncontained), plotStyleEvent, weight);
            }

            if (wasInputToCut("2NonProtons", cutsPassed))
            {
                const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
                unsigned int nNonProtons = 0u;

                for (unsigned int index = 0; index < recoParticles.size(); ++index)
                {
                    if (index == muonIndex)
                    {
                        nNonProtons++;
                        continue;
                    }

                    const auto &particle = recoParticles.at(index);
                    if (!AnalysisHelper::HasTrackFit(particle))
                        continue;

                    std::vector<float> features;
                    const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, protonFeatureNames, features);

                    if (!hasFeatures)
                        continue;

                    const auto protonBDTResponse = protonBDT.GetResponse(features);

                    if (protonBDTResponse < protonBDTCut)
                        nNonProtons++;

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg);
                    protonBDTResponsePlot.Fill(protonBDTResponse, plotStyleParticle, weight);
                    recoProtonMomentumPlot.Fill(AnalysisHelper::GetProtonMomentumFromRange(particle.range()), plotStyleParticle, weight);
                }

                // const auto recoData = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes);
                const auto truthData = AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
                // recoProtonMomentumPlot.Fill(recoData.protonMomentum, plotStyleEvent, weight);
                trueProtonMomentumPlot.Fill(truthData.protonMomentum, plotStyleEvent, weight);

                nNonProtonsPlot.Fill(static_cast<float>(nNonProtons), plotStyleEvent, weight);
            }

            if(wasInputToCut("pionHasValiddEdx", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passesGoldenPionSelection);
                truncatedMeandEdxBeforePlot.Fill(recoData.pionCosTheta, plotStyleEvent, weight);

                if (passedCut("pionHasValiddEdx", cutsPassed))
                {
                    truncatedMeandEdxAfterPlot.Fill(recoData.pionCosTheta, plotStyleEvent, weight);
                }
            }

            if (wasInputToCut("pionNotInGap", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passesGoldenPionSelection);
                pionNotInGapBeforePlot.Fill(recoData.pionPhi, plotStyleEvent, weight);

                if (passedCut("pionNotInGap", cutsPassed))
                {
                    pionNotInGapAfterPlot.Fill(recoData.pionPhi, plotStyleEvent, weight);
                }
            }

            if (wasInputToCut("openingAngle", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passesGoldenPionSelection);
                openingAnglePlot.Fill(recoData.muonPionAngle, plotStyleEvent, weight);
            }

            if (wasInputToCut("topologicalScore", cutsPassed))
            {
                topologicalScorePlot.Fill(pEvent->reco.selectedTopologicalScore(), plotStyleEvent, weight);
            }

            if (wasInputToCut("startNearVertex", cutsPassed))
            {
                float maxVertexDist = -std::numeric_limits<float>::max();

                for (const auto &particle : recoParticles)
                {
                    if (!AnalysisHelper::HasTrackFit(particle))
                        continue;

                    const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
                    const auto &recoVertex = pEvent->reco.nuVertex();
                    const float vertexDist2 = (start - recoVertex).Mag2();
                    const auto vertexDist = std::pow(vertexDist2, 0.5f);

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg);
                    startNearVertexParticlePlot.Fill(vertexDist, plotStyleParticle, weight);
                    maxVertexDist = std::max(maxVertexDist, vertexDist);
                }

                startNearVertexEventPlot.Fill(maxVertexDist, plotStyleEvent, weight);
            }

            if (wasInputToCut("likelyGoldenPion", cutsPassed))
            {
                // Get the pion candidate
                const auto pionIter = std::find_if(assignedPdgCodes.begin(), assignedPdgCodes.end(), [](const auto &pdgCode){ return pdgCode == 211; });
                if (pionIter == assignedPdgCodes.end())
                    throw std::logic_error("PlotEventSelectionCuts - No pion candidate found");

                const auto pion = recoParticles.at(std::distance(assignedPdgCodes.begin(), pionIter));

                std::vector<float> features;
                const auto hasFeatures = BDTHelper::GetBDTFeatures(pion, goldenPionFeatureNames, features);

                if (!hasFeatures)
                    throw std::logic_error("PlotEventSelectionCuts - Pion candidate doesn't have BDT features");

                const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(features);

                const auto plotStyleParticle = PlottingHelper::GetPlotStyle(pion, sampleType, truthParticles, false, config.global.useAbsPdg);
                likelyGoldenPionParticlePlot.Fill(goldenPionBDTResponse, plotStyleParticle, weight);
                likelyGoldenPionEventPlot.Fill(goldenPionBDTResponse, plotStyleEvent, weight);
            }
        }
    }
    const std::string prefix = "CC0pi";
    nTracksPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_min2Tracks_nTracks");
    isContainedPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_max1Uncontained_isContained");
    uncontainedVertexDistPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_max1Uncontained_vertexDist_uncontainedParticles", true);
    nUncontainedPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_max1Uncontained_nUncontained");
    protonBDTResponsePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_2NonProtons_protonBDTResponse");
    recoProtonMomentumPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_recoProtonMomentumPlot");
    trueProtonMomentumPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_trueProtonMomentumPlot");
    nNonProtonsPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_2NonProtons_nNonProtons");
    truncatedMeandEdxBeforePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionTruncatedMeandEdx-before");
    truncatedMeandEdxAfterPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionTruncatedMeandEdx-after");
    openingAnglePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_openingAngle_openingAngle");
    pionNotInGapBeforePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionNotInGap_pionPhi-before");
    pionNotInGapAfterPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionNotInGap_pionPhi-after");
    topologicalScorePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_topologicalScore_topologicalScore", false, false, true);
    startNearVertexParticlePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_startNearVertex_vertexDist_allParticles", true);
    startNearVertexEventPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_startNearVertex_vertexDist_furthestParticle", true);
    likelyGoldenPionParticlePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_likelyGoldenPion_goldenPionBDTResponse_particles");
    likelyGoldenPionEventPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_likelyGoldenPion_goldenPionBDTResponse_events");
}

} // namespace ubcc1pi_macros
