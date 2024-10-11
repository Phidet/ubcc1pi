/**
 *  @file  ubcc1pi_standalone/Macros/PlotReconstructedVariablesRun4Validation.cxx
 *
 *  @brief The implementation file of the PlotReconstructedVariablesRun4Validation macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotReconstructedVariablesRun4Validation(const Config &config)
{
    //
    // Setup the input files
    //

    // std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData123;
    std::vector< std::tuple<AnalysisHelper::SampleType, std::string, float> > inputData1,inputData2,inputData3,inputData4;

    const auto normalisationFactor_123vs4 = (config.normsRun1.dataBNBTor875WCut + config.normsRun2.dataBNBTor875WCut + config.normsRun3.dataBNBTor875WCut)/ config.normsRun4.dataBNBTor875WCut;
    std::cout<<"normalisation factor - Runs 123 vs 4: "<<normalisationFactor_123vs4<<std::endl;
    std::cout<<"config.normsRun1.dataBNBTor875WCut: "<<config.normsRun1.dataBNBTor875WCut<<std::endl;
    std::cout<<"config.normsRun2.dataBNBTor875WCut: "<<config.normsRun2.dataBNBTor875WCut<<std::endl;
    std::cout<<"config.normsRun3.dataBNBTor875WCut: "<<config.normsRun3.dataBNBTor875WCut<<std::endl;
    std::cout<<"config.normsRun4.dataBNBTor875WCut: "<<config.normsRun4.dataBNBTor875WCut<<std::endl;

    inputData1.emplace_back(AnalysisHelper::DataBNB, config.filesRun1.dataBNBFileName, 1.f/normalisationFactor_123vs4);//config.normsRun4.dataBNBTor875WCut/config.normsRun1.dataBNBTor875WCut);
    inputData2.emplace_back(AnalysisHelper::DataBNB, config.filesRun2.dataBNBFileName, 1.f/normalisationFactor_123vs4);//config.normsRun4.dataBNBTor875WCut/config.normsRun2.dataBNBTor875WCut);
    inputData3.emplace_back(AnalysisHelper::DataBNB, config.filesRun3.dataBNBFileName, 1.f/normalisationFactor_123vs4);//config.normsRun4.dataBNBTor875WCut/config.normsRun3.dataBNBTor875WCut);
    inputData4.emplace_back(AnalysisHelper::DataBNB, config.filesRun4.dataBNBFileName, 1.f);

    //
    // Setup the plots
    //

    const std::string yLabel = "Number of particles";


    PlottingHelper::MultiPlot muonMomentumPlot("Muon momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonCosThetaPlot("Muon cos(theta)", yLabel, 50u, config.global.muonCosTheta.min, config.global.muonCosTheta.max, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonPhiPlot("Muon phi / rad", yLabel, 50u, config.global.muonPhi.min, config.global.muonPhi.max, true, config.global.axisTitles);

    PlottingHelper::MultiPlot muonLogBragg_pToMIPPlot("Muon LogBragg_pToMIP", yLabel, 50u, -9.f, 7.f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonLogBragg_piToMIPPlot("Muon LogBragg_piToMIP", yLabel, 50u, -4.f, 6.f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonTruncMeandEdxPlot("Muon TruncMeandEdx", yLabel, 50u, 0.f, 10.f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonProtonForwardPlot("Muon ProtonForward", yLabel, 50u, 0.35f, 0.65f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonMuonForwardPlot("Muon MuonForward", yLabel, 50u, 0.35f, 0.65f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonNDescendentsPlot("Muon NDescendents", yLabel, 4u, 0, 4, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonNSpacePointsNearEndPlot("Muon NSpacePointsNearEnd", yLabel, 50u, 0, 100, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonWigglinessPlot("Muon Wiggliness", yLabel, 50u, 0.f, 0.004f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot muonTrackScorePlot("Muon TrackScore", yLabel, 50u, 0.f, 1.f, true, config.global.axisTitles);

    // PlottingHelper::MultiPlot muonMomentumParticlePlot("Muon momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot muonCosThetaParticlePlot("Muon cos(theta)", yLabel, 50u, config.global.muonCosTheta.min, config.global.muonCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot muonPhiParticlePlot("Muon phi / rad", yLabel, 50u, config.global.muonPhi.min, config.global.muonPhi.max, true, config.global.axisTitles);

    PlottingHelper::MultiPlot pionMomentumPlot("Pion momentum / GeV", yLabel, 50u, 0.f, 0.5f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionCosThetaPlot("Pion cos(theta)", yLabel, 50u, config.global.pionCosTheta.min, config.global.pionCosTheta.max, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionPhiPlot("Pion phi / rad", yLabel, 50u, config.global.pionPhi.min, config.global.pionPhi.max, true, config.global.axisTitles);

    PlottingHelper::MultiPlot pionLogBragg_pToMIPPlot("Pion LogBragg_pToMIP", yLabel, 50u, -9.f, 7.f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionLogBragg_piToMIPPlot("Pion LogBragg_piToMIP", yLabel, 50u, -4.f, 6.f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionTruncMeandEdxPlot("Pion TruncMeandEdx", yLabel, 50u, 0.f, 10.f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionProtonForwardPlot("Pion ProtonForward", yLabel, 50u, 0.35f, 0.65f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionMuonForwardPlot("Pion MuonForward", yLabel, 50u, 0.35f, 0.65f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionNDescendentsPlot("Pion NDescendents", yLabel, 4u, 0, 4, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionNSpacePointsNearEndPlot("Pion NSpacePointsNearEnd", yLabel, 50u, 0, 100, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionWigglinessPlot("Pion Wiggliness", yLabel, 50u, 0.f, 0.004f, true, config.global.axisTitles);
    PlottingHelper::MultiPlot pionTrackScorePlot("Pion TrackScore", yLabel, 50u, 0.f, 1.f, true, config.global.axisTitles);


    // PlottingHelper::MultiPlot pionMomentumParticlePlot("Pion momentum / GeV", yLabel, 50u, 0.f, 0.5f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot pionCosThetaParticlePlot("Pion cos(theta)", yLabel, 50u, config.global.pionCosTheta.min, config.global.pionCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot pionPhiParticlePlot("Pion phi / rad", yLabel, 50u, config.global.pionPhi.min, config.global.pionPhi.max, true, config.global.axisTitles);

    PlottingHelper::MultiPlot muonPionAnglePlot("Muon-pion opening angle / rad", yLabel, 50u, config.global.muonPionAngle.min, config.global.muonPionAngle.max, true, config.global.axisTitles);
    PlottingHelper::MultiPlot nProtonsPlot("Proton multiplicity", yLabel, 5u, 0, 5, true, config.global.axisTitles);

    // PlottingHelper::MultiPlot protonMomentumPlot("Proton momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot protonCosThetaPlot("Proton cos(theta)", yLabel, 50u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot protonPhiPlot("Proton phi / rad", yLabel, 50u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);

    // PlottingHelper::MultiPlot protonMomentumParticlePlot("Proton momentum / GeV", yLabel, 50u, 0.f, 2.f, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot protonCosThetaParticlePlot("Proton cos(theta)", yLabel, 50u, config.global.protonCosTheta.min, config.global.protonCosTheta.max, true, config.global.axisTitles);
    // PlottingHelper::MultiPlot protonPhiParticlePlot("Proton phi / rad", yLabel, 50u, config.global.protonPhi.min, config.global.protonPhi.max, true, config.global.axisTitles);


    // logBragg_pToMIP, logBragg_piToMIP, truncMeandEdx, protonForward, muonForward, nDescendents, nSpacePointsNearEnd, wiggliness, trackScore
    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    std::cout << "DEBUG - featureNames.size() = " << featureNames.size() << std::endl;

    //
    // Get the selection
    //
    // std::cout<<"DEBUG - Point 0"<<std::endl;
    auto selection = SelectionHelper::GetDefaultSelection();
    // std::cout<<"DEBUG - Point 1"<<std::endl;
    const auto run4Data = std::make_pair(PlottingHelper::BNBData, inputData4); // Photon to make multiplot ratio work without changes
    const auto run3Data = std::make_pair(PlottingHelper::Quaternary, inputData3);//Secondary
    const auto run2Data = std::make_pair(PlottingHelper::Quaternary, inputData2);//Tertiary
    const auto run1Data = std::make_pair(PlottingHelper::Quaternary, inputData1);//Quaternary

    // Loop over the events
    for(const auto &[plotStyle, inputData] : {run1Data, run2Data, run3Data, run4Data})
    {
        // std::cout<<"DEBUG - Point 2"<<std::endl;
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

                // Run the event selection and store which cuts are passed
                const auto &[isSelectedGolden, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
                const auto isSelectedGeneric = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());
                // std::cout<<"DEBUG - Point 4"<<std::endl;

                // Only use events that at least pass the generic selection
                if (!isSelectedGeneric)
                    continue;

                // Get the truth and reco analysis data
                // const auto normalisation = AnalysisHelper::GetNominalEventnormalisation(pEvent) * normalisation;
                const auto recoData = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, isSelectedGolden);

                // const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

                // Get the true origin of the selected muon and pion candidates
                const auto &recoParticles = pEvent->reco.particles;
                // const auto &truthParticles = pEvent->truth.particles;

                // std::cout<<"DEBUG - Point 4"<<std::endl;
                // const auto proton = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 2212));
                const auto muon = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 13));
                const auto pion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211));

                std::vector<float> featuresMuon, featuresPion;
                BDTHelper::GetBDTFeatures(muon, featureNames, featuresMuon);
                BDTHelper::GetBDTFeatures(pion, featureNames, featuresPion);

                if(!featuresMuon.empty())
                {
                    muonLogBragg_pToMIPPlot.Fill(featuresMuon.at(0), plotStyle, normalisation);
                    muonLogBragg_piToMIPPlot.Fill(featuresMuon.at(1), plotStyle, normalisation);
                    muonTruncMeandEdxPlot.Fill(featuresMuon.at(2), plotStyle, normalisation);
                    muonProtonForwardPlot.Fill(featuresMuon.at(3), plotStyle, normalisation);
                    muonMuonForwardPlot.Fill(featuresMuon.at(4), plotStyle, normalisation);
                    muonNDescendentsPlot.Fill(featuresMuon.at(5), plotStyle, normalisation);
                    muonNSpacePointsNearEndPlot.Fill(featuresMuon.at(6), plotStyle, normalisation);
                    muonWigglinessPlot.Fill(featuresMuon.at(7), plotStyle, normalisation);
                    muonTrackScorePlot.Fill(featuresMuon.at(8), plotStyle, normalisation);
                }

                if(!featuresPion.empty())
                {
                    pionLogBragg_pToMIPPlot.Fill(featuresPion.at(0), plotStyle, normalisation);
                    pionLogBragg_piToMIPPlot.Fill(featuresPion.at(1), plotStyle, normalisation);
                    pionTruncMeandEdxPlot.Fill(featuresPion.at(2), plotStyle, normalisation);
                    pionProtonForwardPlot.Fill(featuresPion.at(3), plotStyle, normalisation);
                    pionMuonForwardPlot.Fill(featuresPion.at(4), plotStyle, normalisation);
                    pionNDescendentsPlot.Fill(featuresPion.at(5), plotStyle, normalisation);
                    pionNSpacePointsNearEndPlot.Fill(featuresPion.at(6), plotStyle, normalisation);
                    pionWigglinessPlot.Fill(featuresPion.at(7), plotStyle, normalisation);
                    pionTrackScorePlot.Fill(featuresPion.at(8), plotStyle, normalisation);
                }


                // const auto protonPlotStyle = PlottingHelper::GetPlotStyle(proton, sampleType, truthParticles, false, config.global.useAbsPdg);
                // const auto muonPlotStyle = PlottingHelper::GetPlotStyle(muon, sampleType, truthParticles, false, config.global.useAbsPdg);
                // const auto pionPlotStyle = PlottingHelper::GetPlotStyle(pion, sampleType, truthParticles, false, config.global.useAbsPdg);

                // Fill the plots
                muonMomentumPlot.Fill(recoData.muonMomentum, plotStyle, normalisation);
                muonCosThetaPlot.Fill(recoData.muonCosTheta, plotStyle, normalisation);
                muonPhiPlot.Fill(recoData.muonPhi, plotStyle, normalisation);

                // muonMomentumParticlePlot.Fill(recoData.muonMomentum, muonPlotStyle, normalisation);
                // muonCosThetaParticlePlot.Fill(recoData.muonCosTheta, muonPlotStyle, normalisation);
                // muonPhiParticlePlot.Fill(recoData.muonPhi, muonPlotStyle, normalisation);

                // protonMomentumPlot.Fill(recoData.protonMomentum, plotStyle, normalisation);
                // protonCosThetaPlot.Fill(recoData.protonCosTheta, plotStyle, normalisation);
                // protonPhiPlot.Fill(recoData.protonPhi, plotStyle, normalisation);

                // protonMomentumParticlePlot.Fill(recoData.protonMomentum, muonPlotStyle, normalisation);
                // protonCosThetaParticlePlot.Fill(recoData.protonCosTheta, muonPlotStyle, normalisation);
                // protonPhiParticlePlot.Fill(recoData.protonPhi, muonPlotStyle, normalisation);

                if (recoData.hasGoldenPion)
                {
                    pionMomentumPlot.Fill(recoData.pionMomentum, plotStyle, normalisation);
                    // pionMomentumParticlePlot.Fill(recoData.pionMomentum, pionPlotStyle, normalisation);
                }

                pionCosThetaPlot.Fill(recoData.pionCosTheta, plotStyle, normalisation);
                pionPhiPlot.Fill(recoData.pionPhi, plotStyle, normalisation);

                // pionCosThetaParticlePlot.Fill(recoData.pionCosTheta, pionPlotStyle, normalisation);
                // pionPhiParticlePlot.Fill(recoData.pionPhi, pionPlotStyle, normalisation);

                muonPionAnglePlot.Fill(recoData.muonPionAngle, plotStyle, normalisation);
                nProtonsPlot.Fill(recoData.nProtons, plotStyle, normalisation);
            }
        }
    }

    muonMomentumPlot.SaveAsStacked("reco_muonMomentum",false,true,false,config.global.axisTitles);

    muonCosThetaPlot.SaveAsStacked("reco_muonCosTheta",false,true,false,config.global.axisTitles);
    muonPhiPlot.SaveAsStacked("reco_muonPhi",false,true,false,config.global.axisTitles);
    muonWigglinessPlot.SaveAsStacked("reco_muonWiggliness",false,true,false,config.global.axisTitles);

    muonLogBragg_pToMIPPlot.SaveAsStacked("reco_muonLogBragg_pToMIP",false,true,false,config.global.axisTitles);
    muonLogBragg_piToMIPPlot.SaveAsStacked("reco_muonLogBragg_piToMIP",false,true,false,config.global.axisTitles);
    muonTruncMeandEdxPlot.SaveAsStacked("reco_muonTruncMeandEdx",false,true,false,config.global.axisTitles);
    muonProtonForwardPlot.SaveAsStacked("reco_muonProtonForward",false,true,false,config.global.axisTitles);
    muonMuonForwardPlot.SaveAsStacked("reco_muonMuonForward",false,true,false,config.global.axisTitles);
    muonNDescendentsPlot.SaveAsStacked("reco_muonNDescendents",false,true,false,config.global.axisTitles);
    muonNSpacePointsNearEndPlot.SaveAsStacked("reco_muonNSpacePointsNearEnd",false,true,false,config.global.axisTitles);
    muonWigglinessPlot.SaveAsStacked("reco_muonWiggliness",false,true,false,config.global.axisTitles);
    muonTrackScorePlot.SaveAsStacked("reco_muonTrackScore",false,true,false,config.global.axisTitles);

    // muonMomentumParticlePlot.SaveAsStacked("reco_muonMomentum_particle",false,false,config.global.axisTitles);
    // muonCosThetaParticlePlot.SaveAsStacked("reco_muonCosTheta_particle",false,false,config.global.axisTitles);
    // muonPhiParticlePlot.SaveAsStacked("reco_muonPhi_particle",false,false,config.global.axisTitles);

    // protonMomentumPlot.SaveAsStacked("reco_protonMomentum",false,false,config.global.axisTitles);
    // protonCosThetaPlot.SaveAsStacked("reco_protonCosTheta",false,false,config.global.axisTitles);
    // protonPhiPlot.SaveAsStacked("reco_protonPhi",false,false,config.global.axisTitles);

    // protonMomentumParticlePlot.SaveAsStacked("reco_protonMomentum_particle",false,false,config.global.axisTitles);
    // protonCosThetaParticlePlot.SaveAsStacked("reco_protonCosTheta_particle",false,false,config.global.axisTitles);
    // protonPhiParticlePlot.SaveAsStacked("reco_protonPhi_particle",false,false,config.global.axisTitles);

    pionMomentumPlot.SaveAsStacked("reco_pionMomentum",false,true,false,config.global.axisTitles);
    pionCosThetaPlot.SaveAsStacked("reco_pionCosTheta",false,true,false,config.global.axisTitles);
    pionPhiPlot.SaveAsStacked("reco_pionPhi",false,true,false,config.global.axisTitles);
    pionWigglinessPlot.SaveAsStacked("reco_pionWiggliness",false,true,false,config.global.axisTitles);

    pionLogBragg_pToMIPPlot.SaveAsStacked("reco_pionLogBragg_pToMIP",false,true,false,config.global.axisTitles);
    pionLogBragg_piToMIPPlot.SaveAsStacked("reco_pionLogBragg_piToMIP",false,true,false,config.global.axisTitles);
    pionTruncMeandEdxPlot.SaveAsStacked("reco_pionTruncMeandEdx",false,true,false,config.global.axisTitles);
    pionProtonForwardPlot.SaveAsStacked("reco_pionProtonForward",false,true,false,config.global.axisTitles);
    pionMuonForwardPlot.SaveAsStacked("reco_pionMuonForward",false,true,false,config.global.axisTitles);
    pionNDescendentsPlot.SaveAsStacked("reco_pionNDescendents",false,true,false,config.global.axisTitles);
    pionNSpacePointsNearEndPlot.SaveAsStacked("reco_pionNSpacePointsNearEnd",false,true,false,config.global.axisTitles);
    pionWigglinessPlot.SaveAsStacked("reco_pionWiggliness",false,true,false,config.global.axisTitles);
    pionTrackScorePlot.SaveAsStacked("reco_pionTrackScore",false,true,false,config.global.axisTitles);

    // pionMomentumParticlePlot.SaveAsStacked("reco_pionMomentum_particle",false,false,config.global.axisTitles);
    // pionCosThetaParticlePlot.SaveAsStacked("reco_pionCosTheta_particle",false,false,config.global.axisTitles);
    // pionPhiParticlePlot.SaveAsStacked("reco_pionPhi_particle",false,false,config.global.axisTitles);

    muonPionAnglePlot.SaveAsStacked("reco_muonPionAngle",false,true,false,config.global.axisTitles);
    nProtonsPlot.SaveAsStacked("reco_nProtons",false,true,false,config.global.axisTitles);
}

} // namespace ubcc1pi_macros