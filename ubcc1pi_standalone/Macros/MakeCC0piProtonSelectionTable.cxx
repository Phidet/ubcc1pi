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

void MakeCC0piProtonSelectionTable(const Config &config)
{
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get list of good runs
    // -------------------------------------------------------------------------------------------------------------------------------------
    // std::vector<int> goodRuns;
    // std::ifstream file(config.global.goodRunListFile);
    // int run;
    // while (file >> run) {
    //     goodRuns.push_back(run);
    // }

    // Set up the selection
    auto selection = SelectionHelper::GetDefaultSelection2(true);

    // Setup the event counter
    AnalysisHelper::EventCounter counter;

    float numSelectedCC0piEvents = 0.f;
    float numSelectedCC0piEventsWithProton = 0.f; // Whether the pion candidate is a proton (rather than the muon)
    float numSelectedCC0piEventsWithLeadingMomentumProton = 0.f; // The number of events for which the pion candidate is the leading momentum proton
    float numSelectedCC0piEventsWithLowestBDTProton = 0.f; // The number of events for which the pion candidate is the lowest BDT proton

    const auto pProtonBDT = std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames);

    // Loop over the input files
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        if(!useThisFile) continue;
        std::cout << "Processing file - " << filePath << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);

        if(!isOverlay) continue; // Only process overlay files

        const auto isMC = true;

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

            // const auto run = pEventPeLEE->metadata.run();
            // const auto isGoodRun = (isDataBNB || isDataEXT) ? std::find(goodRuns.begin(), goodRuns.end(), run) != goodRuns.end() : true; // Apply good runs cuts to data
            // // if(!isGoodRun) continue;
            // if(!isGoodRun)
            // {
            //     // std::cout << "\tDEBUG - fileName: " << fileName << std::endl;
            //     continue;
            // }

            Event event(*pEventPeLEE, false);
            const auto pEvent = std::make_shared<Event>(event);

            const auto isTrueCC0pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            if(!isTrueCC0pi) continue; // Skip non-CC0Pi events

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            if(assignedPdgCodes.size() != pEvent->reco.particles.size())
                throw std::logic_error("Assigned pdg codes vector has different size to the number of reconstructed particles");
            const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            if(!passedGenericSelection) continue; // Skip events that don't pass the generic selection

            const auto &trueParticles = pEvent->truth.particles;
            double maxProtonMomentum = 0.0;
            auto leadingProtonByMomentumIndex = std::numeric_limits<unsigned int>::max();
            for (unsigned int index = 0; index < trueParticles.size(); ++index)
            {
                const auto &trueParticle = trueParticles.at(index);
                if (std::abs(trueParticle.pdgCode()) == 2212 && trueParticle.momentum() > maxProtonMomentum)
                {
                    maxProtonMomentum = trueParticle.momentum();
                    leadingProtonByMomentumIndex = index;
                }
            }

            if(leadingProtonByMomentumIndex==std::numeric_limits<unsigned int>::max())
                throw std::logic_error("Leading proton index not set for leadingProtonByMomentumIndex.");

            const auto &recoParticles = pEvent->reco.particles;
            const auto falsePionCandidateIndex = AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211);
            unsigned int falsePionTruthMatchIndex = std::numeric_limits<unsigned int>::max();
            try
            {
                falsePionTruthMatchIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticles.at(falsePionCandidateIndex), trueParticles, true);
            }
            catch(const std::exception& e)
            {
                std::cout<<"try failed for falsePionTruthMatchIndex"<<std::endl;
                continue;
            }

            if(falsePionTruthMatchIndex==std::numeric_limits<unsigned int>::max())
                throw std::logic_error("falsePionTruthMatchIndex not set");

            auto lowestProtonBDT = std::numeric_limits<float>::max();
            auto lowestProtonBDTMatchedIndex = std::numeric_limits<unsigned int>::max();
            for(unsigned int index = 0; index < recoParticles.size(); ++index)
            {
                // std::cout<<"index: "<<index<<std::endl;
                const auto &recoParticle = recoParticles.at(index);
                if(recoParticle.generation() != 2) continue; // This is needed because the PeLEE ntuples contain all particles
                unsigned int matchedTruthParticleIndex = std::numeric_limits<unsigned int>::max();
                try
                {                
                    matchedTruthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, trueParticles, true);
                }
                catch(const std::exception& e)
                {
                    std::cout<<"try failed for matchedTruthParticleIndex"<<std::endl;
                    continue;
                }

                const auto matchedTruthParticle = trueParticles.at(matchedTruthParticleIndex);
                if(matchedTruthParticle.pdgCode() != 2212) continue;

                // Get run the proton BDT
                std::vector<float> features;
                const auto hasFeatures = BDTHelper::GetBDTFeatures(recoParticle, BDTHelper::ProtonBDTFeatureNames, features);
                if (!hasFeatures) continue;
                const auto bdtResponse = pProtonBDT->GetResponse(features);

                if(bdtResponse < lowestProtonBDT)
                {
                    lowestProtonBDT = bdtResponse;
                    lowestProtonBDTMatchedIndex = matchedTruthParticleIndex;
                }
            }

            // if(lowestProtonBDTMatchedIndex==std::numeric_limits<unsigned int>::max())
            //     continue;
            //     throw std::logic_error("lowestProtonBDTMatchedIndex not set");

            const auto falsePionTruthMatch = trueParticles.at(falsePionTruthMatchIndex);
            if(falsePionTruthMatch.pdgCode() == 2212)
            {
                numSelectedCC0piEventsWithProton += weight;
            }
            // else
            // {
            //     std::cout<<"falsePionTruthMatch.pdgCode(): "<<falsePionTruthMatch.pdgCode()<<std::endl;
            // }

            if(falsePionTruthMatchIndex == leadingProtonByMomentumIndex) numSelectedCC0piEventsWithLeadingMomentumProton += weight;
            if(falsePionCandidateIndex == lowestProtonBDTMatchedIndex) numSelectedCC0piEventsWithLowestBDTProton += weight;
            numSelectedCC0piEvents += weight; 

        } // end of event loop
    } // end of file loop

    std::cout<<"Results: "<<std::endl;
    std::cout<<"numSelectedCC0piEvents: "<<numSelectedCC0piEvents<<std::endl;
    const auto fractionProtonIsPionCandidate = numSelectedCC0piEventsWithProton/numSelectedCC0piEvents;
    std::cout<<"\tFraction of which the pion candidate is a proton: "<<fractionProtonIsPionCandidate<<std::endl;
    const auto fractionLargestProtonMomentumIsPionCandidate = numSelectedCC0piEventsWithLeadingMomentumProton/numSelectedCC0piEvents;
    std::cout<<"\tFraction of which the largest momentum proton is the pion candidate: "<<fractionLargestProtonMomentumIsPionCandidate<<std::endl;
    const auto fractionLowestBDTProtonIsPionCandidate = numSelectedCC0piEventsWithLowestBDTProton/numSelectedCC0piEvents;
    std::cout<<"\t fraction of which the lowest BDT proton is the pion candidate: "<<fractionLowestBDTProtonIsPionCandidate<<std::endl;

} // end of macro

} // namespace ubcc1pi_macros