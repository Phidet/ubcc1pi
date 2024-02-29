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
    auto selection = SelectionHelper::GetDefaultSelection2(true); // todo decide on final selection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const auto allCuts = selection.GetCuts();

    const auto protonBDTCut = selection.GetCutValue("2NonProtons");
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get list of good runs
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::vector<int> goodRuns;
    std::ifstream file(config.global.goodRunListFile);
    int run;
    while (file >> run) {
        goodRuns.push_back(run);
    }

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
    PlottingHelper::MultiPlot nNonProtonsPlot("Number of non-protons", yLabelEvents, 5u, 1, 6);
    PlottingHelper::MultiPlot truncatedMeandEdxPlot("Pion truncated dE/dx / MeV/cm", yLabelEvents, 40u, 0.f, 4.5f);
    // PlottingHelper::MultiPlot truncatedMeandEdxBeforePlot("Pion cos(theta) / rad", yLabelEvents, 40u, -1.f, 1.f);
    // PlottingHelper::MultiPlot truncatedMeandEdxAfterPlot("Pion cos(theta) / rad", yLabelEvents, 40u, -1.f, 1.f);
    PlottingHelper::MultiPlot pionNotInGapBeforePlot("Pion phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot pionNotInGapAfterPlot("Pion phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot muonNotInGapBeforePlot("Muon phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot muonNotInGapAfterPlot("Muon phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
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
    truncatedMeandEdxPlot.AddCutLine(selection.GetCutValue("pionHasValiddEdx"));

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
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        if(!useThisFile) continue;        
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);
        const auto isMC = (sampleType != AnalysisHelper::DataBNB) && (sampleType != AnalysisHelper::DataEXT);

        if(sampleType != AnalysisHelper::Overlay && sampleType != AnalysisHelper::Dirt && sampleType != AnalysisHelper::DataBNB && sampleType != AnalysisHelper::DataEXT) continue;

        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        if (isMC) readerPeLEE.EnableSystematicBranches(); // Todo: Is this correct/optimal?
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

        // Loop over the events
        // std::cout << "### Only processing 5\% of events ###" << std::endl;
        for (unsigned int i = 0; i < nEvents; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);

            readerPeLEE.LoadEvent(i);

            const auto run = pEventPeLEE->metadata.run();
            const auto isGoodRun = (isDataBNB || isDataEXT) ? std::find(goodRuns.begin(), goodRuns.end(), run) != goodRuns.end() : true; // Apply good runs cuts to data
            if(!isGoodRun)
            {
                // std::cout << "DEBUG - bad run: "<<run<<std::endl;
                continue;
            }

            Event event(*pEventPeLEE, true); // here we cut out the generation > 2 particles
            const auto pEvent = std::make_shared<Event>(event);

            // For brevity assign the particles a variable
            const auto &truthParticles = pEvent->truth.particles;
            const auto &recoParticles = pEvent->reco.particles;

            // Run the event selection and store which cuts are passed
            // std::cout<<"DEBUG Point 0"<<std::endl;
            const auto &[passesGoldenPionSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            // std::cout<<"DEBUG Point 1"<<std::endl;

            // Get the plot style of the event
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyleEvent = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg, true);
            // std::cout<<"DEBUG Point 2"<<std::endl;

            // if(PlottingHelper::PlotStyle::NumuCC0Pi != plotStyleEvent) continue;

            // const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold); // todo remove this
            // if(!isTrueCC0Pi) continue;

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

            // std::cout<<"DEBUG Point 3"<<std::endl;

            if (wasInputToCut("max1Uncontained", cutsPassed))
            {
                unsigned int nUncontained = 0u;
                for (const auto &recoParticle : recoParticles)
                {
                    // std::cout<<"DEBUG Point 3.1"<<std::endl;
                    if (!AnalysisHelper::HasTrackFit(recoParticle))
                        continue;

                    const auto isContained = AnalysisHelper::IsContained(recoParticle);
                    nUncontained += isContained ? 0u : 1u;

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(recoParticle, sampleType, truthParticles, false, config.global.useAbsPdg, true);
                    isContainedPlot.Fill(isContained ? 1.f : 0.f, plotStyleParticle, weight);

                    // std::cout<<"DEBUG Point 3.2"<<std::endl;
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

            // std::cout<<"DEBUG Point 3.3"<<std::endl;

            if (wasInputToCut("2NonProtons", cutsPassed))
            {
                // const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
                const auto muonIndex = AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 13);

                unsigned int nNonProtons = 0u;

                for (unsigned int index = 0; index < recoParticles.size(); ++index)
                {
                    // std::cout<<"DEBUG Point 3.4"<<std::endl;
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
                    
                    // std::cout<<"DEBUG Point 3.5"<<std::endl;
                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg, true);
                    protonBDTResponsePlot.Fill(protonBDTResponse, plotStyleParticle, weight);
                    recoProtonMomentumPlot.Fill(AnalysisHelper::GetProtonMomentumFromRange(particle.range()), plotStyleParticle, weight);
                    // std::cout<<"DEBUG Point 3.6"<<std::endl;
                }

                // // std::cout<<"DEBUG Point 3.7"<<std::endl;
                // // const auto recoData = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes);
                // const auto truthData = AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
                // // recoProtonMomentumPlot.Fill(recoData.protonMomentum, plotStyleEvent, weight);
                // trueProtonMomentumPlot.Fill(truthData.protonMomentum, plotStyleEvent, weight);
                // // std::cout<<"DEBUG Point 3.8"<<std::endl;

                nNonProtonsPlot.Fill(static_cast<float>(nNonProtons), plotStyleEvent, weight);
                // // std::cout<<"DEBUG Point 3.9"<<std::endl;
            }

            // std::cout<<"DEBUG Point 4"<<std::endl;

            if(wasInputToCut("pionHasValiddEdx", cutsPassed))
            {
                // Get the pion candidate
                const auto pionIter = std::find_if(assignedPdgCodes.begin(), assignedPdgCodes.end(), [](const auto &pdgCode){ return pdgCode == 211; });
                if (pionIter == assignedPdgCodes.end())
                    throw std::logic_error("PlotEventSelectionCuts - No pion candidate found");

                const auto pion = recoParticles.at(std::distance(assignedPdgCodes.begin(), pionIter));
                // const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passesGoldenPionSelection);
                truncatedMeandEdxPlot.Fill(pion.truncatedMeandEdx(), plotStyleEvent, weight);
                // truncatedMeandEdxBeforePlot.Fill(pion.truncatedMeandEdx(), plotStyleEvent, weight);

                // if (passedCut("pionHasValiddEdx", cutsPassed))
                // {
                //     truncatedMeandEdxAfterPlot.Fill(pion.truncatedMeandEdx(), plotStyleEvent, weight);
                // }
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

            if (wasInputToCut("muonNotInGap", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passesGoldenPionSelection);
                muonNotInGapBeforePlot.Fill(recoData.muonPhi, plotStyleEvent, weight);

                if (passedCut("muonNotInGap", cutsPassed))
                {
                    muonNotInGapAfterPlot.Fill(recoData.muonPhi, plotStyleEvent, weight);
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

            // std::cout<<"DEBUG Point 5"<<std::endl;

            if (wasInputToCut("startNearVertex", cutsPassed))
            {
                float maxVertexDist = -std::numeric_limits<float>::max();

                for (const auto &particle : recoParticles)
                {
                    // if(particle.reco.generation>2) continue;
                    if (!AnalysisHelper::HasTrackFit(particle))
                        continue;

                    const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
                    const auto &recoVertex = pEvent->reco.nuVertex();
                    const float vertexDist2 = (start - recoVertex).Mag2();
                    const auto vertexDist = std::pow(vertexDist2, 0.5f);

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg, true);
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

                const auto plotStyleParticle = PlottingHelper::GetPlotStyle(pion, sampleType, truthParticles, false, config.global.useAbsPdg, true);
                likelyGoldenPionParticlePlot.Fill(goldenPionBDTResponse, plotStyleParticle, weight);
                likelyGoldenPionEventPlot.Fill(goldenPionBDTResponse, plotStyleEvent, weight);
            }

            // std::cout<<"DEBUG Point 6"<<std::endl;
        }
    }

    const std::string prefix = "CC1pi";//"CC0pi";
    nTracksPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_min2Tracks_nTracks");
    isContainedPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_max1Uncontained_isContained");
    uncontainedVertexDistPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_max1Uncontained_vertexDist_uncontainedParticles", true);
    nUncontainedPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_max1Uncontained_nUncontained");
    protonBDTResponsePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_2NonProtons_protonBDTResponse");
    recoProtonMomentumPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_recoProtonMomentumPlot");
    trueProtonMomentumPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_trueProtonMomentumPlot");
    nNonProtonsPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_2NonProtons_nNonProtons");
    truncatedMeandEdxPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionTruncatedMeandEdx");
    // truncatedMeandEdxBeforePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionTruncatedMeandEdx-before");
    // truncatedMeandEdxAfterPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionTruncatedMeandEdx-after");
    openingAnglePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_openingAngle_openingAngle");
    pionNotInGapBeforePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionNotInGap_pionPhi-before");
    pionNotInGapAfterPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_pionNotInGap_pionPhi-after");
    muonNotInGapBeforePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_muonNotInGap_muonPhi-before");
    muonNotInGapAfterPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_muonNotInGap_muonPhi-after");
    topologicalScorePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_topologicalScore_topologicalScore", false, false, true);
    startNearVertexParticlePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_startNearVertex_vertexDist_allParticles", true);
    startNearVertexEventPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_startNearVertex_vertexDist_furthestParticle", true);
    likelyGoldenPionParticlePlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_likelyGoldenPion_goldenPionBDTResponse_particles");
    likelyGoldenPionEventPlot.SaveAsStacked("plotEventSelectionCuts_" + prefix + "_likelyGoldenPion_goldenPionBDTResponse_events");
}

} // namespace ubcc1pi_macros
