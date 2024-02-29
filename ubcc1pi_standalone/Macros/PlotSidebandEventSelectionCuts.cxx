/**
 *  @file  ubcc1pi_standalone/Macros/PlotSidebandEventSelectionCuts.cxx
 *
 *  @brief The implementation file of the PlotSidebandEventSelectionCuts macro
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

unsigned int GetLeadingProtonIndex(const std::vector<Event::Reco::Particle> &recoParticles, const std::vector<int> assignedPdgCodes)
{
    auto leadingProtonIndex = -std::numeric_limits<int>::max();
    auto highestMom = -std::numeric_limits<float>::max();
    auto highestProtonMom = -std::numeric_limits<float>::max();
    auto highestRange = -std::numeric_limits<float>::max();
    // std::cout<<"DEBUG GetLeadingProtonIndex - recoParticles.size(): "<<recoParticles.size()<<std::endl;
    for (unsigned int index = 0; index < recoParticles.size(); ++index)
    {
        const auto recoParticle = recoParticles.at(index);
        const auto pdgCode = assignedPdgCodes.at(index);
        // std::cout<<"DEBUG GetLeadingProtonIndex - pdgCode: "<<pdgCode<<std::endl;
        // std::cout<<"DEBUG GetLeadingProtonIndex - range: "<<range<<std::endl;
        if(pdgCode!=2212) continue;
        if (!AnalysisHelper::HasTrackFit(recoParticle))
            continue;
        std::vector<float> features;
        const auto hasFeatures = BDTHelper::GetBDTFeatures(recoParticle, BDTHelper::ProtonBDTFeatureNames, features);
        if (!hasFeatures)
            continue;
        auto range = recoParticle.range();
        const auto protonmom = AnalysisHelper::GetProtonMomentumFromRange(range);
        const auto momentum = AnalysisHelper::GetPionMomentumFromRange(range); //ATTN: Use pion momentum here for reco values
        // std::cout<<"DEBUG GetLeadingProtonIndex - momentum: "<<momentum<<std::endl;
        // std::cout<<"DEBUG GetLeadingProtonIndex - truncatedMeandEdx() V2: "<<recoParticle.truncatedMeandEdx()<<std::endl;
        if (momentum>highestMom)
        {
            if(range<=highestRange)
                throw std::logic_error("PlotEventSelectionCuts - Proton momentum is not monotonically increasing with range");
            if(protonmom<=highestProtonMom)
                throw std::logic_error("PlotEventSelectionCuts - Proton momentum is not monotonically increasing with range");
            leadingProtonIndex = index;
            highestMom=momentum;
            highestRange=range;
            highestProtonMom=protonmom;
        }
    }
    if(leadingProtonIndex<0)
        throw std::logic_error("PlotEventSelectionCuts - No leading proton candidate found");

    // std::cout<<"DEBUG GetLeadingProtonIndex - leadingProtonIndex: "<<leadingProtonIndex<<std::endl;
    return leadingProtonIndex;
}

void PlotSidebandEventSelectionCuts(const Config &config)
{
    // auto selection = SelectionHelper::GetCC0piSelectionModifiedPeLEE(-0.48f, 0.12f, true);
    auto selection = SelectionHelper::GetCC0piSelectionModifiedAgainPeLEE(-0.48f, 0.12f, true);
    const auto allCuts = selection.GetCuts();

    const auto protonBDTCut = selection.GetCutValue("1NonProton");
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

    // min2Tracks
    PlottingHelper::MultiPlot nTracksPlot("Number of tracks", yLabelEvents, 6u, 1, 7);

    PlottingHelper::MultiPlot isContainedPlot("Is contained?", yLabelEvents, 2u, 0, 2);
    // max1Uncontained
    PlottingHelper::MultiPlot uncontainedVertexDistPlot("Distance to vertex / cm", yLabelParticles, PlottingHelper::GenerateLogBinEdges(40u, 0.15f, 1000.f));
    PlottingHelper::MultiPlot nUncontainedPlot("Number of uncontained particles", yLabelEvents, 4u, 0, 4);
    // 1NonProton && atLeast1Proton
    PlottingHelper::MultiPlot protonBDTResponsePlot("Proton BDT response", yLabelParticles, 40u, -0.60f, 0.60f);
    PlottingHelper::MultiPlot nNonProtonsPlot("Number of non-protons", yLabelEvents, 5u, 1, 6);
    PlottingHelper::MultiPlot nProtonsPlot("Number of protons", yLabelEvents, 6u, 0, 6);
    // muonLikeProton
    // PlottingHelper::MultiPlot muonLikeProtonPlot("Muon-like proton", yLabelEvents, 2u, 0, 2);
    PlottingHelper::MultiPlot muonLikeProtonBeforePlot("Muon-like proton", yLabelParticles, 40u, 0.f, 0.6f);
    PlottingHelper::MultiPlot muonLikeProtonAfterPlot("Muon-like proton", yLabelParticles, 40u, 0.f, 0.6f);
    // barelyResemblingProton
    // PlottingHelper::MultiPlot barelyResemblingProtonPlot("Barely resembling proton", yLabelEvents, 2u, 0, 2);
    PlottingHelper::MultiPlot barelyResemblingProtonBeforePlot("Barely resembling proton", yLabelParticles, 40u, 0.f, 0.6f);
    PlottingHelper::MultiPlot barelyResemblingProtonAfterPlot("Barely resembling proton", yLabelParticles, 40u, 0.f, 0.6f);

    // protonHasValiddEdx
    PlottingHelper::MultiPlot truncatedMeandEdxPlot("Leading proton as pion truncated dE/dx / MeV / cm", yLabelEvents, 40u, 0.f, 4.5f);
    // PlottingHelper::MultiPlot truncatedMeandEdxBeforePlot("Pion cos(theta) / rad", yLabelEvents, 40u, -1.f, 1.f);
    // PlottingHelper::MultiPlot truncatedMeandEdxAfterPlot("Pion cos(theta) / rad", yLabelEvents, 40u, -1.f, 1.f);

    // muonNotInGap
    PlottingHelper::MultiPlot muonNotInGapBeforePlot("Muon phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot muonNotInGapAfterPlot("Muon phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    // protonNotInGap
    PlottingHelper::MultiPlot protonNotInGapBeforePlot("Leading proton phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    PlottingHelper::MultiPlot protonNotInGapAfterPlot("Leading proton phi / rad", yLabelEvents, 40u, -3.142f, 3.142f);
    // openingAngle
    PlottingHelper::MultiPlot openingAnglePlot("Muon-Proton opening angle / rad", yLabelEvents, 40u, 0.f, 3.142f);
    // topologicalScore
    PlottingHelper::MultiPlot topologicalScorePlot("TopologicalScore", yLabelEvents, 40u, 0.06f, 1.f);
    // startNearVertex
    PlottingHelper::MultiPlot startNearVertexParticlePlot("Distance to vertex / cm", yLabelParticles, PlottingHelper::GenerateLogBinEdges(40u, 0.03f, 1000.f));
    PlottingHelper::MultiPlot startNearVertexEventPlot("Min distance to vertex / cm", yLabelEvents, PlottingHelper::GenerateLogBinEdges(40u, 0.15f, 1000.f));
    // likelyGoldenProton
    PlottingHelper::MultiPlot likelyGoldenProtonParticlePlot("Golden pion BDT response to leading proton", yLabelParticles, 40u, -0.55f, 0.4f);
    PlottingHelper::MultiPlot likelyGoldenProtonEventPlot("Golden pion BDT response to leading proton", yLabelEvents, 40u, -0.55f, 0.4f);

    PlottingHelper::MultiPlot resultGeneric("Final generic total selection", yLabelEvents, 1u, 0, 1);
    PlottingHelper::MultiPlot resultGolden("Final golden total selection", yLabelEvents, 1u, 0, 1);

    // min2Tracks
    // max1Uncontained
    // 1NonProton
    // atLeast1Proton
    // muonLikeProton
    // barelyResemblingProton
    // protonHasValiddEdx
    // muonNotInGap
    // protonNotInGap
    // openingAngle
    // topologicalScore
    // startNearVertex
    // likelyGoldenProton


    // Set the bin labels where appropriate
    nTracksPlot.SetIntegerBinLabels();
    nUncontainedPlot.SetIntegerBinLabels();
    nNonProtonsPlot.SetIntegerBinLabels();
    nProtonsPlot.SetIntegerBinLabels();
    isContainedPlot.SetBinLabels({"No", "Yes"});

    // muonLikeProtonPlot.SetBinLabels({"No", "Yes"});
    // barelyResemblingProtonPlot.SetBinLabels({"No", "Yes"});

    // Add the cut values
    nTracksPlot.AddCutLine(2);
    nUncontainedPlot.AddCutLine(2);
    protonBDTResponsePlot.AddCutLine(protonBDTCut);
    nNonProtonsPlot.AddCutLine(1);
    nNonProtonsPlot.AddCutLine(2);
    nProtonsPlot.AddCutLine(1);
    openingAnglePlot.AddCutLine(selection.GetCutValue("openingAngle"));
    topologicalScorePlot.AddCutLine(selection.GetCutValue("topologicalScore"));
    startNearVertexParticlePlot.AddCutLine(selection.GetCutValue("startNearVertex"));
    startNearVertexEventPlot.AddCutLine(selection.GetCutValue("startNearVertex"));
    likelyGoldenProtonParticlePlot.AddCutLine(selection.GetCutValue("likelyGoldenProton"));
    likelyGoldenProtonEventPlot.AddCutLine(selection.GetCutValue("likelyGoldenProton"));
    truncatedMeandEdxPlot.AddCutLine(selection.GetCutValue("protonHasValiddEdx"));

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
        std::cout << "### Only processing 5\% of events ###" << std::endl;
        for (unsigned int i = 0; i < nEvents/20; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            const auto &[passesGoldenProtonSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Get the plot style of the event
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            const auto plotStyleEvent = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            // if(PlottingHelper::PlotStyle::NumuCC0Pi != plotStyleEvent) continue;
            // const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
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

            if (wasInputToCut("max1Uncontained", cutsPassed))
            {
                unsigned int nUncontained = 0u;
                for (const auto &recoParticle : recoParticles)
                {
                    if (!AnalysisHelper::HasTrackFit(recoParticle))
                        continue;

                    const auto isContained = AnalysisHelper::IsContained(recoParticle);
                    nUncontained += isContained ? 0u : 1u;

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(recoParticle, sampleType, truthParticles, false, config.global.useAbsPdg, true);
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

            if (wasInputToCut("1NonProton", cutsPassed))
            {
                // const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
                const auto muonIndex = AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 13);
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
                        nNonProtons++;;

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg, true);
                    protonBDTResponsePlot.Fill(protonBDTResponse, plotStyleParticle, weight);
                }

                nNonProtonsPlot.Fill(static_cast<float>(nNonProtons), plotStyleEvent, weight);
            }

            if (wasInputToCut("atLeast1Proton", cutsPassed))
            {
                // const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
                const auto muonIndex = AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 13);
                unsigned int nProtons = 0u;

                for (unsigned int index = 0; index < recoParticles.size(); ++index)
                {
                    if (index == muonIndex)
                        continue;

                    const auto &particle = recoParticles.at(index);
                    if (!AnalysisHelper::HasTrackFit(particle))
                        continue;

                    std::vector<float> features;
                    const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, protonFeatureNames, features);

                    if (!hasFeatures)
                        continue;

                    const auto protonBDTResponse = protonBDT.GetResponse(features);

                    if (protonBDTResponse > protonBDTCut)
                        nProtons++;
                }

                nProtonsPlot.Fill(static_cast<float>(nProtons), plotStyleEvent, weight);
            }

            // if(wasInputToCut("muonLikeProton", cutsPassed))
            // {
            //     const auto passesCut = wasInputToCut("barelyResemblingProtonPlot", cutsPassed);
            //     muonLikeProtonPlot.Fill(passesCut, plotStyleEvent, weight);
            // }

            // if(wasInputToCut("BarelyResemblingProton", cutsPassed))
            // {
            //     const auto passesCut = wasInputToCut("protonHasValiddEdx", cutsPassed);
            //     barelyResemblingProtonPlot.Fill(passesCut, plotStyleEvent, weight);
            // }


            if (wasInputToCut("muonLikeProton", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, passesGoldenProtonSelection);
                muonLikeProtonBeforePlot.Fill(recoData.protonMomentum, plotStyleEvent, weight);

                if (passedCut("muonLikeProton", cutsPassed))
                {
                    muonLikeProtonAfterPlot.Fill(recoData.protonMomentum, plotStyleEvent, weight);
                }
            }

            if (wasInputToCut("barelyResemblingProton", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, passesGoldenProtonSelection);
                barelyResemblingProtonBeforePlot.Fill(recoData.protonMomentum, plotStyleEvent, weight);

                if (passedCut("barelyResemblingProton", cutsPassed))
                {
                    barelyResemblingProtonAfterPlot.Fill(recoData.protonMomentum, plotStyleEvent, weight);
                }
            }



            if(wasInputToCut("protonHasValiddEdx", cutsPassed))
            {
                const auto leadingProtonIndex = GetLeadingProtonIndex(recoParticles, assignedPdgCodes);
                const auto &proton = recoParticles.at(leadingProtonIndex);
                truncatedMeandEdxPlot.Fill(proton.truncatedMeandEdx(), plotStyleEvent, weight);
            }

            if (wasInputToCut("muonNotInGap", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, passesGoldenProtonSelection);
                muonNotInGapBeforePlot.Fill(recoData.muonPhi, plotStyleEvent, weight);

                if (passedCut("muonNotInGap", cutsPassed))
                {
                    muonNotInGapAfterPlot.Fill(recoData.muonPhi, plotStyleEvent, weight);
                }
            }

            if (wasInputToCut("protonNotInGap", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, passesGoldenProtonSelection);
                protonNotInGapBeforePlot.Fill(recoData.protonPhi, plotStyleEvent, weight);

                if (passedCut("protonNotInGap", cutsPassed))
                {
                    protonNotInGapAfterPlot.Fill(recoData.protonPhi, plotStyleEvent, weight);
                }
            }

            if (wasInputToCut("openingAngle", cutsPassed))
            {
                const auto recoData = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, passesGoldenProtonSelection);
                openingAnglePlot.Fill(recoData.muonProtonAngle, plotStyleEvent, weight);
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

                    const auto plotStyleParticle = PlottingHelper::GetPlotStyle(particle, sampleType, truthParticles, false, config.global.useAbsPdg, true);
                    startNearVertexParticlePlot.Fill(vertexDist, plotStyleParticle, weight);
                    maxVertexDist = std::max(maxVertexDist, vertexDist);
                }

                startNearVertexEventPlot.Fill(maxVertexDist, plotStyleEvent, weight);
            }

            if (wasInputToCut("likelyGoldenProton", cutsPassed))
            {
                // // Get the pion candidate
                // const auto pionIter = std::find_if(assignedPdgCodes.begin(), assignedPdgCodes.end(), [](const auto &pdgCode){ return pdgCode == 211; });
                // if (pionIter == assignedPdgCodes.end())
                //     throw std::logic_error("PlotSidebandEventSelectionCuts - No pion candidate found");

                // const auto pion = recoParticles.at(std::distance(assignedPdgCodes.begin(), pionIter));
                const auto leadingProtonIndex = GetLeadingProtonIndex(recoParticles, assignedPdgCodes);
                const auto &proton = recoParticles.at(leadingProtonIndex);
                std::vector<float> features;
                const auto hasFeatures = BDTHelper::GetBDTFeatures(proton, goldenPionFeatureNames, features);

                if (!hasFeatures)
                    throw std::logic_error("PlotSidebandEventSelectionCuts - Pion candidate doesn't have BDT features");

                const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(features);

                const auto plotStyleParticle = PlottingHelper::GetPlotStyle(proton, sampleType, truthParticles, false, config.global.useAbsPdg, true);
                likelyGoldenProtonParticlePlot.Fill(goldenPionBDTResponse, plotStyleParticle, weight);
                likelyGoldenProtonEventPlot.Fill(goldenPionBDTResponse, plotStyleEvent, weight);
            }

            if(wasInputToCut("likelyGoldenProton", cutsPassed))
            {
                resultGeneric.Fill(0.5f, plotStyleEvent, weight);
            }

            if(passesGoldenProtonSelection)
            {
                resultGolden.Fill(0.5f, plotStyleEvent, weight);
            }
        }
    }
    const std::string prefix = "CC0pi";
    nTracksPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_min2Tracks_nTracks");
    isContainedPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_max1Uncontained_isContained");
    uncontainedVertexDistPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_max1Uncontained_vertexDist_uncontainedParticles", true);
    nUncontainedPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_max1Uncontained_nUncontained");
    protonBDTResponsePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_1NonProton_protonBDTResponse");
    nNonProtonsPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_1NonProton_nNonProtons");
    nProtonsPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_AtLeast1Proton_nProtons");
    // muonLikeProtonPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_MuonLikeProton_muonLikeProton");
    muonLikeProtonBeforePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_MuonLikeProton_muonLikeProton_before");
    muonLikeProtonAfterPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_MuonLikeProton_muonLikeProton_after");
    // barelyResemblingProtonPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_BarelyResemblingProton_barelyResemblingProton");
    barelyResemblingProtonBeforePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_BarelyResemblingProton_barelyResemblingProton_before");
    barelyResemblingProtonAfterPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_BarelyResemblingProton_barelyResemblingProton_after");
    truncatedMeandEdxPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_leadingProtonTruncatedMeandEdx");
    // truncatedMeandEdxBeforePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_protonTruncatedMeandEdx-before");
    // truncatedMeandEdxAfterPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_protonTruncatedMeandEdx-after");
    openingAnglePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_openingAngle_openingAngle");
    muonNotInGapBeforePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_muonNotInGap_muonPhi-before");
    muonNotInGapAfterPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_muonNotInGap_muonPhi-after");
    protonNotInGapBeforePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_protonNotInGap_protonPhi-before");
    protonNotInGapAfterPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_protonNotInGap_protonPhi-after");
    topologicalScorePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_topologicalScore_topologicalScore", false, false, true);
    startNearVertexParticlePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_startNearVertex_vertexDist_allParticles", true);
    startNearVertexEventPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_startNearVertex_vertexDist_furthestParticle", true);
    likelyGoldenProtonParticlePlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_likelyGoldenProton_goldenPionBDTResponse_particles");
    likelyGoldenProtonEventPlot.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_likelyGoldenProton_goldenPionBDTResponse_events");
    resultGeneric.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_result_generic");
    resultGolden.SaveAsStacked("plotSidebandEventSelectionCuts_" + prefix + "_result_golden");
}

} // namespace ubcc1pi_macros
