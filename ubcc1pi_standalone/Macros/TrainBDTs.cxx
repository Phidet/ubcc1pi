/**
 *  @file  ubcc1pi_standalone/Macros/TrainBDTs.cxx
 *
 *  @brief The implementation file of the TrainBDTs macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
// #include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void TrainBDTs(const Config &config)
{
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Extract the CC1Pi events tha pass the pre-selection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::cout << "Finding CC1Pi events" << std::endl;
    auto ccInclusiveSelection = SelectionHelper::GetCCInclusiveSelection2(true); // todo decide on final selection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    std::cout << "DEBUG Point 0" << std::endl;
    std::vector<std::pair<unsigned int, unsigned int>> cc1PiEventIndices;
    // for (const auto &[sampleType, fileRun, fileName, normalisation] : inputData)
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {   
        std::cout << "DEBUG Point 1" << std::endl;
        if(!useThisFile || sampleType != AnalysisHelper::Overlay) continue;
        // Read the input file
        std::cout << "DEBUG Point 2" << std::endl;
        std::cout<<"Using file: "<<filePath<<std::endl;
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, true);
        std::cout<<"WARNING: Only using 5\% of events!"<<std::endl; // TODO: Remove this line
        const auto nEvents = readerPeLEE.GetNumberOfEvents()/20;
        std::cout << "DEBUG Point 3" << std::endl;
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        std::cout << "DEBUG Point 4" << std::endl;
        for (unsigned int eventIndex = 0; eventIndex < nEvents; ++eventIndex)
        {
            AnalysisHelper::PrintLoadingBar(eventIndex, nEvents);
            readerPeLEE.LoadEvent(eventIndex);
            Event event(*pEventPeLEE, true);// true or false decides wether to cut generation!=2 particles
            const auto pEvent = std::make_shared<Event>(event);

            // Event must be true CC1Pi
            if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
                continue;

            // Event must pass the CCInclusive selection
            const auto &[passesCCInclusive, cutsPassed, assignedPdgCodes] = ccInclusiveSelection.Execute(pEvent);
            if (!passesCCInclusive)
                continue;

            cc1PiEventIndices.emplace_back(fileRun, eventIndex);
        }
        std::cout << "DEBUG Point 5" << std::endl;
    }
    std::cout << "DEBUG Point 6" << std::endl;
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Randomly choose the training events
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const auto nCC1PiEvents = cc1PiEventIndices.size();
    const auto nTrainingEvents = static_cast<unsigned int>(std::floor(static_cast<float>(nCC1PiEvents) * config.trainBDTs.trainingFraction));
    BDTHelper::EventShuffler shuffler(nCC1PiEvents, nTrainingEvents);
    std::cout << "Found " << nCC1PiEvents << " CC1Pi events passing CC inclusive seleciton. Using " << nTrainingEvents << " for training." << std::endl;
    if(nCC1PiEvents==0) throw std::runtime_error("No CC1Pi events found. Check input files and selection.");

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Setup the BDTs to train
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDTFactory goldenPionBDTFactory("goldenPion", goldenPionFeatureNames);
    BDTHelper::BDTFactory protonBDTFactory("proton", protonFeatureNames);
    BDTHelper::BDTFactory muonBDTFactory("muon", muonFeatureNames);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Fill the BDT training and testing entries
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {   
        if(!useThisFile || sampleType != AnalysisHelper::Overlay) continue;
        // Read the input file
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, true);
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

        std::cout << "Filling the "<<nCC1PiEvents<<" BDT entries for file "<<filePath<<std::endl;
        for (unsigned int i = 0; i < nCC1PiEvents; ++i)
        {
            const auto run = cc1PiEventIndices.at(i).first;
            if(fileRun != run) continue;

            AnalysisHelper::PrintLoadingBar(i, nCC1PiEvents);

            const auto eventIndex = cc1PiEventIndices.at(i).second;
            const auto isTrainingEvent = shuffler.IsTrainingEvent(i);
            readerPeLEE.LoadEvent(i);
            Event event(*pEventPeLEE, true); // true or false decides wether to cut generation!=2 particles
            const auto pEvent = std::make_shared<Event>(event);

            const auto truthParticles = pEvent->truth.particles;
            const auto recoParticles = pEvent->reco.particles;
            const float eventWeight = AnalysisHelper::GetNominalEventWeight(pEvent)*normalisation;

            for (const auto &recoParticle : recoParticles)
            {
                if(recoParticle.generation() > 2) throw std::runtime_error("Found reco particle with generation > 2");

                // Only use contained particles for training
                if (!AnalysisHelper::HasTrackFit(recoParticle) || !AnalysisHelper::IsContained(recoParticle))
                    continue;

                // Determine the true origin of the reco particle
                bool isExternal = true;
                int truePdgCode = -std::numeric_limits<int>::max();
                bool trueIsGolden = false;
                float trueMomentum = -std::numeric_limits<float>::max();
                float completeness = -std::numeric_limits<float>::max();

                if(recoParticle.pdgBacktracked.IsSet()) std::cout<<"DEBUG recoParticle.pdgBacktracked(): "<<recoParticle.pdgBacktracked()<<std::endl;
                else std::cout<<"DEBUG recoParticle.pdgBacktracked() not set"<<std::endl;

                try
                {
                    const auto truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles);
                    std::cout<<"DEBUG truthParticleIndex: "<<truthParticleIndex<<std::endl;
                    const auto truthParticle = truthParticles.at(truthParticleIndex);

                    isExternal = false;
                    truePdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
                    trueMomentum = truthParticle.momentum();
                    trueIsGolden = AnalysisHelper::IsGolden(truthParticle);

                    std::cout<<"DEBUG truePdgCode: "<<truePdgCode<<" trueMomentum: "<<trueMomentum<<" trueIsGolden: "<<trueIsGolden<<std::endl;

                    // completeness = recoParticle.truthMatchCompletenesses().at(truthParticleIndex);
                    completeness = recoParticle.completenessBacktracked();
                }
                catch (const std::exception &e) 
                {
                    std::cout << "Exception caught in TrainBDTs.cxx: " << e.what() << std::endl;
                }

                // Only use good matches for training
                if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                {
                    std::cout << "DEBUG: Skipping due to values - onlyGoodTruthMatches: " << config.trainBDTs.onlyGoodTruthMatches 
                              << ", isExternal: " << isExternal 
                              << ", completeness: " << completeness << std::endl;
                    continue;
                }

                std::cout<<"DEBUG: Passed onlyGoodTruthMatches"<<std::endl;

                // Extract the features
                std::vector<float> goldenPionFeatures;
                const auto areAllFeaturesAvailableGoldenPion = BDTHelper::GetBDTFeatures(recoParticle, goldenPionFeatureNames, goldenPionFeatures, true); // true for debugging

                std::vector<float> protonFeatures;
                const auto areAllFeaturesAvailableProton = BDTHelper::GetBDTFeatures(recoParticle, protonFeatureNames, protonFeatures, true); // true for debugging

                std::vector<float> muonFeatures;
                const auto areAllFeaturesAvailableMuon = BDTHelper::GetBDTFeatures(recoParticle, muonFeatureNames, muonFeatures, true); // true for debugging

                // Define the weight
                const auto completenessWeight = (config.trainBDTs.weightByCompleteness ? (isExternal ? 1.f : completeness) : 1.f);
                std::cout<<"DEBUG completenessWeight: "<<completenessWeight<<" isExternal: "<<isExternal<<std::endl;
                const auto weight = eventWeight * completenessWeight;

                if (areAllFeaturesAvailableGoldenPion)
                {
                    const bool isGoldenPion = !isExternal && truePdgCode == 211 && trueIsGolden;
                    std::cout<<"DEBUG Added golden pion entry - isGoldenPion: "<<isGoldenPion<<" !isExternal: "<<!isExternal<<" truePdgCode: "<<truePdgCode<<" trueIsGolden: "<<trueIsGolden<<std::endl;
                    goldenPionBDTFactory.AddEntry(goldenPionFeatures, isGoldenPion, isTrainingEvent, weight);
                }
                else
                {
                    std::cout<<"DEBUG Failed to add golden pion entry with true pdg code: "<<truePdgCode<<std::endl;
                }

                if (areAllFeaturesAvailableProton)
                {
                    const bool isProton = !isExternal && truePdgCode == 2212;
                    std::cout<<"DEBUG Added proton entry - isProton: "<<isProton<<" !isExternal: "<<!isExternal<<" truePdgCode: "<<truePdgCode<<std::endl;
                    protonBDTFactory.AddEntry(protonFeatures, isProton, isTrainingEvent, weight);
                }
                else
                {
                    std::cout<<"DEBUG Failed to add proton entry with true pdg code: "<<truePdgCode<<std::endl;
                }

                if (areAllFeaturesAvailableMuon)
                {
                    const bool isMuon = !isExternal && truePdgCode == 13;
                    std::cout<<"DEBUG Added muon entry - isMuon: "<<isMuon<<" !isExternal: "<<!isExternal<<" truePdgCode: "<<truePdgCode<<std::endl;
                    muonBDTFactory.AddEntry(muonFeatures, isMuon, isTrainingEvent, weight);
                }
                else
                {
                    std::cout<<"DEBUG Failed to add muon entry with true pdg code: "<<truePdgCode<<std::endl;
                }
            }
        }
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Optimize the BDTs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (config.trainBDTs.shouldOptimize)
    {
        std::cout << "Optimizing golden pion BDT" << std::endl;
        goldenPionBDTFactory.OptimizeParameters();

        std::cout << "Optimizing proton BDT" << std::endl;
        protonBDTFactory.OptimizeParameters();

        std::cout << "Optimizing muon BDT" << std::endl;
        muonBDTFactory.OptimizeParameters();
    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Train and test the BDTs
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::cout << "Training and testing golden pion BDT" << std::endl;
    goldenPionBDTFactory.TrainAndTest();

    std::cout << "Training and testing proton BDT" << std::endl;
    protonBDTFactory.TrainAndTest();

    std::cout << "Training and testing muon BDT" << std::endl;
    muonBDTFactory.TrainAndTest();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Make plots
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!config.trainBDTs.shouldMakePlots)
        return;

    const std::string yLabel = "Fraction of reco particles";
    PlottingHelper::MultiPlot goldenPionBDTPlot("Golden pion BDT response", yLabel, 50, -0.9f, 0.45f);
    PlottingHelper::MultiPlot protonBDTPlot("Proton BDT response", yLabel, 50, -0.8f, 0.7f);
    PlottingHelper::MultiPlot muonBDTPlot("Muon BDT response", yLabel, 50, -0.9f, 0.6f);

    // Using the newly trained BDT weight files, setup up a BDT for evaluation
    BDTHelper::BDT goldenPionBDT("goldenPion", goldenPionFeatureNames);
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

    std::cout << "Making BDT training plots" << std::endl;
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {   
        if(!useThisFile || sampleType != AnalysisHelper::Overlay) continue;
        // Read the input file
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, true);
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

        std::cout << "Filling the BDT entries" << std::endl;

        for (unsigned int i = 0; i < nCC1PiEvents; ++i)
        {
            const auto run = cc1PiEventIndices.at(i).first;
            if(fileRun != run) continue;

            AnalysisHelper::PrintLoadingBar(i, nCC1PiEvents);
            const auto eventIndex = cc1PiEventIndices.at(i).second;
            const auto isTrainingEvent = shuffler.IsTrainingEvent(i);
            readerPeLEE.LoadEvent(i);
            Event event(*pEventPeLEE, true);// true or false decides wether to cut generation!=2 particles
            const auto pEvent = std::make_shared<Event>(event);

            const auto truthParticles = pEvent->truth.particles;
            const auto recoParticles = pEvent->reco.particles;
            const float eventWeight = AnalysisHelper::GetNominalEventWeight(pEvent)*normalisation;

            for (const auto &recoParticle : recoParticles)
            {
                // Only use contained particles for testing
                if (!AnalysisHelper::HasTrackFit(recoParticle) || !AnalysisHelper::IsContained(recoParticle))
                    continue;

                // Determine the true origin of the reco particle
                bool isExternal = true;
                int truePdgCode = -std::numeric_limits<int>::max();
                bool trueIsGolden = false;
                float completeness = -std::numeric_limits<float>::max();

                try
                {
                    const auto truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles);
                    const auto truthParticle = truthParticles.at(truthParticleIndex);

                    isExternal = false;
                    truePdgCode = config.global.useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
                    trueIsGolden = AnalysisHelper::IsGolden(truthParticle);
                    // completeness = recoParticle.truthMatchCompletenesses().at(truthParticleIndex);
                    completeness = recoParticle.completenessBacktracked();
                }
                catch (const std::exception &e)
                {
                    std::cout << "Exception caught in TrainBDTs.cxx (Point 2): " << e.what() << std::endl;
                }
                // Only use good matches for testing
                if (config.trainBDTs.onlyGoodTruthMatches && (isExternal || completeness < 0.5f))
                    continue;

                std::cout << "DEBUG: Passed good matches test" << std::endl;

                // Fill to the plots
                const auto style = PlottingHelper::GetPlotStyle(recoParticle, AnalysisHelper::Overlay, truthParticles, isTrainingEvent, config.global.useAbsPdg, true);
                std::cout << "DEBUG style: " << style << std::endl;

                // For these plots skip neutrons
                if (style == PlottingHelper::Other || style == PlottingHelper::OtherPoints)
                    continue;

                std::cout << "DEBUG: Passed neutron skip test" << std::endl;

                if (style == PlottingHelper::Electron || style == PlottingHelper::ElectronPoints ||
                    style == PlottingHelper::Photon || style == PlottingHelper::PhotonPoints) continue; // Skip shower particles

                std::cout << "DEBUG: Passed shower particles skip test" << std::endl;

                // Extract the features
                std::vector<float> goldenPionFeatures;
                const auto areAllFeaturesAvailableGoldenPion = BDTHelper::GetBDTFeatures(recoParticle, goldenPionFeatureNames, goldenPionFeatures);

                std::cout << "DEBUG: Extracted goldenPionFeatures" << std::endl;

                std::vector<float> protonFeatures;
                const auto areAllFeaturesAvailableProton = BDTHelper::GetBDTFeatures(recoParticle, protonFeatureNames, protonFeatures);

                std::cout << "DEBUG: Extracted protonFeatures" << std::endl;

                std::vector<float> muonFeatures;
                const auto areAllFeaturesAvailableMuon = BDTHelper::GetBDTFeatures(recoParticle, muonFeatureNames, muonFeatures);

                std::cout << "DEBUG: Extracted muonFeatures" << std::endl;

                if (areAllFeaturesAvailableGoldenPion)
                {
                    const auto goldenPionBDTResponse = goldenPionBDT.GetResponse(goldenPionFeatures);
                    goldenPionBDTPlot.Fill(goldenPionBDTResponse, style, eventWeight);
                }

                std::cout << "DEBUG: Filled goldenPionBDTPlot" << std::endl;

                if (areAllFeaturesAvailableProton)
                {
                    const auto protonBDTResponse = protonBDT.GetResponse(protonFeatures);
                    protonBDTPlot.Fill(protonBDTResponse, style, eventWeight);
                }

                std::cout << "DEBUG: Filled protonBDTPlot" << std::endl;

                if (areAllFeaturesAvailableMuon)
                {
                    const auto muonBDTResponse = muonBDT.GetResponse(muonFeatures);
                    muonBDTPlot.Fill(muonBDTResponse, style, eventWeight);
                }

                std::cout << "DEBUG: Filled muonBDTPlot" << std::endl;
            }
        }
    }

    // Save the plots
    goldenPionBDTPlot.SaveAs("goldenPionBDTResponse");
    protonBDTPlot.SaveAs("protonBDTResponse");
    muonBDTPlot.SaveAs("muonBDTResponse");
}

} // namespace ubcc1pi_plots
