/**
 *  @file  ubcc1pi_standalone/Macros/AnalyzerTest.cxx
 *
 *  @brief The implementation file of the AnalyzerTest macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// #include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
// #include "ubcc1pi_standalone/ubsmear/inc/ubsmear/Helpers/UBSmearingHelper.h"

// Boost libraries
// #include "binary_iarchive.hpp"
#include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"


#define COMPARETWOSINGLEVALUES(A, X, Y, Z) \
    if(X.Z.IsSet() != Y.Z.IsSet()) std::cout<<"! "<<#A<<" "<<#Z<<" - one value not set: "<<X.Z.IsSet()<<(X.Z.IsSet() ? std::string(" (")+std::to_string(X.Z())+std::string(")"):std::string(" "))<<" vs "<<Y.Z.IsSet()<<(Y.Z.IsSet() ? std::string(" (")+std::to_string(Y.Z())+std::string(")"):std::string(" "))<<std::endl; \
    else if(!X.Z.IsSet() && !Y.Z.IsSet()) std::cout<<"✓ "<<#A<<" "<<#Z<<" Both values not set."<<std::endl; \
    else {if(X.Z() - Y.Z() > 1e-3f ) std::cout<<"✕ "<<#A<<" "<<#Z<<" values not identical in PeLEE and ubcc1pi events: "<<X.Z()<<" vs "<<Y.Z()<<std::endl; \
    else std::cout<<"✓ "<<#A<<" "<<#Z<<" values identical in PeLEE and ubcc1pi events: "<<X.Z()<<" vs "<<Y.Z()<<std::endl;}

#define COMPAREANYTWOVALUES(A, X, Y, Z) \
    if(X.Z.IsSet() != Y.Z.IsSet()) std::cout<<"! "<<#A<<" "<<#Z<<" - one value not set: "<<X.Z.IsSet()<<" vs "<<Y.Z.IsSet()<<std::endl; \
    else if(!X.Z.IsSet() && !Y.Z.IsSet()) std::cout<<"✓ "<<#A<<" "<<#Z<<" Both values not set."<<std::endl; \
    else {if(X.Z() != Y.Z()) std::cout<<"✕ "<<#A<<" "<<#Z<<" values not identical in PeLEE and ubcc1pi events."<<std::endl; \
    else std::cout<<"✓ "<<#A<<" "<<#Z<<" values identical in PeLEE and ubcc1pi events."<<std::endl;}

#define CHECKANYVALUE(A, X, Z) \
    std::cout<<#A<<" "<<#Z<<" - value "; \
    if(X.Z.IsSet()) std::cout<<"set"<<std::endl; \
    else std::cout<<"not set"<<std::endl; \

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void AnalyzerTest(const Config &config)
{

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;
    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

        // The names of the cross-section kinematic parameters, and their binning information.
        // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
        { "muonCosTheta",  config.global.muonCosTheta,  true  },
        { "muonPhi",       config.global.muonPhi,       true  },
        { "muonMomentum",  config.global.muonMomentum,  true  },

        { "pionCosTheta",  config.global.pionCosTheta,  true  },
        { "pionPhi",       config.global.pionPhi,       true  },
        { "pionMomentum",  config.global.pionMomentum,  true  },

        { "muonPionAngle", config.global.muonPionAngle, true  },
        { "nProtons",      config.global.nProtons,      false }

    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));
    }

    // Set output precision
    std::cout << std::fixed;
    std::cout << std::setprecision(4);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the plots
    // -------------------------------------------------------------------------------------------------------------------------------------
    const std::string yLabelParticles = "Number of particles";
    const std::string yLabelEvents = "Number of events";
    // PlottingHelper::MultiPlot truth_particle_mc_pdg_plot("mc pdg", yLabelParticles,0u, 3000u, 10);
    // PlottingHelper::MultiPlot truth_particle_mc_endx_plot("mc endx", yLabelParticles,-5.f, 5.f, 10);
    // PlottingHelper::MultiPlot truth_particle_mc_endy_plot("mc endy", yLabelParticles,-5.f, 5.f, 10);
    // PlottingHelper::MultiPlot truth_particle_mc_endz_plot("mc endz", yLabelParticles,-5.f, 5.f, 10);
    // PlottingHelper::MultiPlot truth_particle_mc_px_plot("mc px", yLabelParticles,-5.f, 5.f, 10);
    // PlottingHelper::MultiPlot truth_particle_mc_py_plot("mc py", yLabelParticles,-5.f, 5.f, 10);
    // PlottingHelper::MultiPlot truth_particle_mc_pz_plot("mc pz", yLabelParticles,-5.f, 5.f, 10);
    // PlottingHelper::MultiPlot truth_particle_mc_E_plot("mc E", yLabelParticles,-5.f, 5.f, 10);

    // // Set the bin labels where appropriate
    // truth_particle_mc_pdg_plot.SetIntegerBinLabels();

    const auto ll = -100.f; // std::numeric_limits<float>::lowest();
    const auto ul = 100.f; // std::numeric_limits<float>::max();

    // PlottingHelper::MultiPlot truth_particle_pdgCode_plot("truth particle pdgCode", yLabelParticles, 10u, 0u, 3000u);
    // PlottingHelper::MultiPlot truth_particle_endX_plot("truth particle endX", yLabelParticles, 10u, ll, ul);
    // PlottingHelper::MultiPlot truth_particle_endY_plot("truth particle endY", yLabelParticles, 10u, ll, ul);
    // PlottingHelper::MultiPlot truth_particle_endZ_plot("truth particle endZ", yLabelParticles, 10u, ll, ul);
    // PlottingHelper::MultiPlot truth_particle_momentumX_plot("truth particle momentumX", yLabelParticles, 10u, ll, ul);
    // PlottingHelper::MultiPlot truth_particle_momentumY_plot("truth particle momentumY", yLabelParticles, 10u, ll, ul);
    // PlottingHelper::MultiPlot truth_particle_momentumZ_plot("truth particle momentumZ", yLabelParticles, 10u, ll, ul);
    // // PlottingHelper::MultiPlot truth_particle_momentum_plot("truth particle momentum", yLabelParticles, 10u, -5.f, 5.f);
    // PlottingHelper::MultiPlot truth_particle_energy_plot("truth particle energy", yLabelParticles, 10u, ll, ul);

    // // Set the bin labels where appropriate
    // truth_particle_pdgCode_plot.SetIntegerBinLabels();


    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section and for the sideband
    // -------------------------------------------------------------------------------------------------------------------------------------
    // ExtractionHelper::AnalysisValueMap getValue;
    // ExtractionHelper::AnalysisValueMap getSidebandValue;
    // ExtractionHelper::PopulateAnalysisValueMap(getValue, false);
    // ExtractionHelper::PopulateAnalysisValueMap(getSidebandValue, true); // true: creates sideband getters

    // std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-0.48f, barelyResemblingProtonValue=0.12f\n.........................................."<<std::endl;
    // auto sidebandSelection = SelectionHelper::GetCC0piSelectionModified(-0.48f, 0.12f);
    auto ubcc1piSelection = SelectionHelper::GetCC0piSelection();
    auto peleeSelection = SelectionHelper::GetCC0piSelectionOriginalPeLEE();

    // auto ubcc1piSelection = SelectionHelper::GetOriginalSelection(); // Relies on the pre-applied CC-inclusive filter
    // auto peleeSelection = SelectionHelper::GetDefaultSelection(); // Recreates cuts from the CC-inclusive filter in the selection

    ExtractionHelper::InputFileList inputData;
    // float totalExposurePOT;
    // ExtractionHelper::PopulateInputFileList(config, inputData, totalExposurePOT);

    // inputData.push_back(std::make_tuple(AnalysisHelper::Overlay,"",std::string("/pnfs/uboone/scratch/users/jdetje/v08_00_00_60/PeLee_Run4b_BNB_Overlay_for_ubcc1pi/0/0/0/57790885_0/neutrinoselection_filt_49039757-ae93-4b79-8111-a141719aaef3.root"),1.f)); // outdated
    // inputData.push_back(std::make_tuple(AnalysisHelper::Overlay,"",std::string("/pnfs/uboone/scratch/users/jdetje/v08_00_00_60/PeLee_Run4b_BNB_Overlay_for_ubcc1pi_2/0/0/0/60169221_0/neutrinoselection_filt_f34088d8-f55d-46e6-b783-8bf197e1e999.root"),1.f));

    // inputData.push_back(std::make_tuple(AnalysisHelper::Overlay,"",std::string("/pnfs/uboone/scratch/users/jdetje/v08_00_00_60/PeLee_Run4b_BNB_Overlay_for_ubcc1pi_4/0/0/8/66842097_8/neutrinoselection_filt_ce69f77c-3147-4279-a07d-799714e646d6.root"),1.f));

    // inputData.push_back(std::make_tuple(AnalysisHelper::Overlay,"",std::string("/uboone/app/users/jdetje/searchingfornues/files/steps2/PhysicsRun-2019_1_17_7_54_14-eventweight_20230414T155652_eventweight_extragenie1_20230414T164323_eventweight_extragenie2_20230414T165358_eventweight_extragenie3_20230414T170405_eventweight_extragenie4.root"),1.f));
    inputData.push_back(std::make_tuple(AnalysisHelper::Overlay, "", std::string("/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root"), 1.f));
    // inputData.push_back(std::make_tuple(AnalysisHelper::Overlay,"",std::string("/uboone/app/users/jdetje/searchingfornues/files/steps3/PhysicsRun_eventweight_20230428T142324_eventweight_extragenie1_20230428T151854_eventweight_extragenie2_20230428T152700_eventweight_extragenie3_20230428T153504_eventweight_extragenie4.root"),1.f));
    inputData.push_back(std::make_tuple(AnalysisHelper::DataBNB, "", std::string("/uboone/data/users/jdetje/ubcc1piVSpelee/ubcc1pi/ubcc1piAnalysis_0_4k.root"), 1.f));

    // Loop over the files
    FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(std::get<2>(inputData.at(0)), true);
    FileReader<Event, Subrun> reader(std::get<2>(inputData.at(1)), true);

    if(readerPeLEE.GetNumberOfEvents()!=reader.GetNumberOfEvents())
    {
        std::cout<<"ERROR: Number of events in PeLEE and ubcc1pi files are different - PeLEE: "<<readerPeLEE.GetNumberOfEvents()<<" vs ubcc1pi: "<<reader.GetNumberOfEvents()<<std::endl;
        throw std::logic_error("Unequal number of events in PeLEE and ubcc1pi files!");
    }

    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;
        // const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        // const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        // const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        // const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        // const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        // const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);

        std::cout << "DEBUG Point Y0" << std::endl;
        const auto isPeLEE    = (sampleType == AnalysisHelper::Overlay);
        std::cout << "DEBUG Point Y1" << std::endl;
        const auto isUbcc1pi  = (sampleType == AnalysisHelper::DataBNB);
        std::cout << "DEBUG Point Y2" << std::endl;

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        // FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(fileName);
        // FileReader<Event, Subrun> reader(fileName);
        // FileReader<Event, Subrun> reader(fileName);

        // if (isOverlay || isNuWro) // Todo: Is isNuWro needed ?
        //     reader.EnableSystematicBranches();

        std::cout << "DEBUG Point Y3" << std::endl;
        auto pEvent = reader.GetBoundEventAddress();
        std::cout << "DEBUG Point Y4" << std::endl;
        auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        std::cout << "DEBUG Point Y5" << std::endl;

        // Loop over the events in the file
        // const auto nEvents = isPeLEE ? readerPeLEE.GetNumberOfEvents() : reader.GetNumberOfEvents();
        const auto nEvents = std::min(readerPeLEE.GetNumberOfEvents(), reader.GetNumberOfEvents());

        const auto style = isPeLEE ? PlottingHelper::ExternalPoints : PlottingHelper::External;//BNBData;

        std::cout << "DEBUG Point Y6" << std::endl;
        // ******************************************************************************************
        // Check the input ntuples
        // ******************************************************************************************
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            std::cout << "DEBUG Point Y7" << std::endl;
            readerPeLEE.LoadEvent(i);
            std::cout << "DEBUG Point Y8" << std::endl;
            reader.LoadEvent(i);
            std::cout << "DEBUG Point Y9" << std::endl;

            // const auto event = isPeLEE ? static_cast<Event>(*pEventPeLEE) : *pEvent;
            const auto hasTruth = true;
            Event event(*pEventPeLEE);//, hasTruth);
            // event.Print();
            std::cout << "DEBUG Point Y10" << std::endl;

            // std::cout<<"("<<i<<")\nPeLEE: "<<"metadata.sub: "<<pEventPeLEE->metadata.sub()<<" "<<pEventPeLEE->metadata.run()<<std::endl;
            // std::cout<<"ubcc1pi: "<<"metadata.subRun: "<<event.metadata.subRun()<<" "<<event.metadata.run()<<std::endl;
            // std::cout<<"PeLEE: "<<"truth.ccnc: "<<pEventPeLEE->truth.ccnc()<<std::endl;
            // std::cout<<"ubcc1pi: "<<"truth.isCC: "<<event.truth.isCC()<<std::endl;

            // for(unsigned int j = 0; j < pEventPeLEE->truth.particles.size(); j++)
            // {
            //     truth_particle_mc_pdg_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_pdg()), style, 1.f); // weight = 1.f for now
            //     truth_particle_mc_endx_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_endx()), style, 1.f);
            //     truth_particle_mc_endy_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_endy()), style, 1.f);
            //     truth_particle_mc_endz_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_endz()), style, 1.f);
            //     truth_particle_mc_px_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_px()), style, 1.f);
            //     truth_particle_mc_py_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_py()), style, 1.f);
            //     truth_particle_mc_pz_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_pz()), style, 1.f);
            //     truth_particle_mc_E_plot.Fill(static_cast<float>(pEventPeLEE->truth.particles.at(j).mc_E()), style, 1.f);
            // }

            // ******************************************************************************************
            std::cout<<"\n~~~CHECKING METADATA~~~"<<std::endl;
            // ******************************************************************************************
            COMPARETWOSINGLEVALUES(Metadata, event.metadata, pEvent->metadata, run)
            COMPARETWOSINGLEVALUES(Metadata, event.metadata, pEvent->metadata, subRun)
            COMPARETWOSINGLEVALUES(Metadata, event.metadata, pEvent->metadata, event)
            // hasTruthInfo

            // ******************************************************************************************
            std::cout<<"\n~~~CHECKING TRUTH EVENT~~~"<<std::endl;
            // ******************************************************************************************
            const auto pTruth = &event.truth;
            const auto pTruthUbcc1pi = &pEvent->truth;

            COMPARETWOSINGLEVALUES(Truth, event.truth, pEvent->truth, splineEventWeight)
            COMPARETWOSINGLEVALUES(Truth, event.truth, pEvent->truth, genieTuneEventWeight)
            COMPAREANYTWOVALUES(Truth, event.truth, pEvent->truth, systParamNames)

            if(!pTruth->systParamFirstValueIndex.IsSet() || !pTruthUbcc1pi->systParamFirstValueIndex.IsSet()) std::cout<<"Truth systParamFirstValueIndex - at least one value not set: "<<pTruth->systParamFirstValueIndex.IsSet()<<" vs "<<pTruthUbcc1pi->systParamFirstValueIndex.IsSet()<<std::endl;
            else if(pTruth->systParamFirstValueIndex() != pTruthUbcc1pi->systParamFirstValueIndex())
            {
                std::cout<<"Truth systParamFirstValueIndex values not identical in PeLEE and ubcc1pi events."; // Might need to use weightsFlux, weightsGenie, weightsReint at some point to sort out the indices
                for(const auto &v : pTruth->systParamFirstValueIndex()) std::cout<<v<<std::endl;
                std::cout<<"\n"<<std::endl;
                for(const auto &v : pTruthUbcc1pi->systParamFirstValueIndex()) std::cout<<v<<std::endl;
            }

            COMPAREANYTWOVALUES(Truth, event.truth, pEvent->truth, systParamValues)
            std::cout<<"PeLEE: ";
            if(event.truth.systParamValues.IsSet()) for (const auto & value : event.truth.systParamValues()) std::cout<<value<<" ";
            std::cout<<std::endl;
            std::cout<<"ubcc1pi: ";
            if(pEvent->truth.systParamValues.IsSet()) for (const auto & value : pEvent->truth.systParamValues()) std::cout<<value<<" ";
            std::cout<<std::endl;

            COMPARETWOSINGLEVALUES(Truth, event.truth, pEvent->truth, isCC)
            COMPARETWOSINGLEVALUES(Truth, event.truth, pEvent->truth, interactionMode)
            // interactionString = DebugHelper::GetInteractionString(eventPeLEE.truth.interaction(), true); // To do: Needs to be implemented in DebugHelper
            COMPARETWOSINGLEVALUES(Truth, event.truth, pEvent->truth, nuPdgCode)
            COMPARETWOSINGLEVALUES(Truth, event.truth, pEvent->truth, nuEnergy)


            if(!pTruth->nuVertex.IsSet() || !pTruthUbcc1pi->nuVertex.IsSet()) std::cout<<"Truth nuVertex - at least one value not set: "<<pTruth->nuVertex.IsSet()<<" vs "<<pTruthUbcc1pi->nuVertex.IsSet()<<std::endl;
            else
            {
                const auto sameVertex = (pTruth->nuVertex() - pTruthUbcc1pi->nuVertex()).Mag() < 1e-3f;
                if(!sameVertex)
                std::cout<<"Truth nuVertex values not identical in PeLEE and ubcc1pi events: "<<pTruth->nuVertex().X()<<" "<<pTruth->nuVertex().Y()<<" "<<pTruth->nuVertex().Z() \
                                                                              <<" vs "<<pTruthUbcc1pi->nuVertex().X()<<" "<<pTruthUbcc1pi->nuVertex().Y()<<" "<<pTruthUbcc1pi->nuVertex().Z()<<std::endl;
            }

            // nFinalStates
            // COMPAREANYTWOVALUES(Truth, event.truth, pEvent->truth, slicePurities)
            // COMPAREANYTWOVALUES(Truth, event.truth, pEvent->truth, sliceCompletenesses)

            // ******************************************************************************************
            std::cout<<"\n~~~CHECKING TRUTH PARTICLES~~~"<<std::endl;
            // ******************************************************************************************
            if(event.truth.particles.size()!=pEvent->truth.particles.size())
            {
                std::cout<<"ERROR: Number of truth particles in PeLEE and ubcc1pi events are different: "<<event.truth.particles.size()<<" vs "<<pEvent->truth.particles.size()<<std::endl;
                // continue; // todo: go back to throwing an exception
                // throw std::logic_error("Unequal number of truth particles in PeLEE and ubcc1pi events!");
            }

            for(unsigned long j = 0; j < std::max(event.truth.particles.size(), pEvent->truth.particles.size()); j++)
            {
                std::cout<<"-----"<<j<<"-----"<<std::endl;
                const auto j1 = std::min(j, event.truth.particles.size()-1);
                const auto j2 = std::min(j, pEvent->truth.particles.size()-1);

                const auto particle = event.truth.particles.at(j1);
                const auto particleUbcc1pi = pEvent->truth.particles.at(j2);
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, pdgCode)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, startX)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, startY)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, startZ)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, endX)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, endY)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, endZ)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, momentumX)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, momentumY)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, momentumZ)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, momentum)

                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, energy)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, endMomentum)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, isStopping)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, nElasticScatters)
                COMPARETWOSINGLEVALUES(Truth Particles, particle, particleUbcc1pi, nInelasticScatters)

                // truth_particle_pdgCode_plot.Fill(static_cast<float>(pParticle->pdgCode()), style, 1.f); // weight = 1.f for now
                // truth_particle_endX_plot.Fill(pParticle->endX(), style, 1.f);
                // truth_particle_endY_plot.Fill(pParticle->endY(), style, 1.f);
                // truth_particle_endZ_plot.Fill(pParticle->endZ(), style, 1.f);
                // truth_particle_momentumX_plot.Fill(pParticle->momentumX(), style, 1.f);
                // truth_particle_momentumY_plot.Fill(pParticle->momentumY(), style, 1.f);
                // truth_particle_momentumZ_plot.Fill(pParticle->momentumZ(), style, 1.f);
                // // truth_particle_momentum_plot
                // truth_particle_energy_plot.Fill(pParticle->energy(), style, 1.f);
            }

            // ******************************************************************************************
            std::cout<<"\n~~~CHECKING RECO EVENT~~~"<<std::endl;
            // ******************************************************************************************
            const auto pReco = &event.reco;
            const auto pRecoUbcc1pi = &pEvent->reco;
            COMPARETWOSINGLEVALUES(Reco, event.reco, pEvent->reco, passesCCInclusive)
            // COMPARETWOSINGLEVALUES(Reco, event.reco, pEvent->reco, nSlices)
            // hasSelectedSlice
            COMPARETWOSINGLEVALUES(Reco, event.reco, pEvent->reco, selectedTopologicalScore)
            // sliceTopologicalScores
            // sliceIsSelectedAsNu
            // hasNeutrino
            COMPARETWOSINGLEVALUES(Reco, event.reco, pEvent->reco, nuPdgCode)
            if(!pReco->nuVertex.IsSet() || !pRecoUbcc1pi->nuVertex.IsSet()) std::cout<<"Reco nuVertex - at least one value not set: "<<pReco->nuVertex.IsSet()<<" vs "<<pRecoUbcc1pi->nuVertex.IsSet()<<std::endl;
            else if((pReco->nuVertex() - pRecoUbcc1pi->nuVertex()).Mag() > 1e-3f) std::cout<<"Reco nuVertex values not identical - "<<pReco->nuVertex().X()<<" "<< pReco->nuVertex().Y()<<" "<< pReco->nuVertex().Z() \
                                                                     << " vs " << pRecoUbcc1pi->nuVertex().X()<<" "<<pRecoUbcc1pi->nuVertex().Y()<<" "<<pRecoUbcc1pi->nuVertex().Z()<<std::endl;
            if(!pReco->nuVertexNoSCC.IsSet() || !pRecoUbcc1pi->nuVertexNoSCC.IsSet()) std::cout<<"Reco nuVertexNoSCC - at least one value not set: "<<pReco->nuVertexNoSCC.IsSet()<<" vs "<<pRecoUbcc1pi->nuVertexNoSCC.IsSet()<<std::endl;
            else if((pReco->nuVertexNoSCC() - pRecoUbcc1pi->nuVertexNoSCC()).Mag() > 1e-3f) std::cout<<"Reco nuVertexNoSCC values not identical - "<<pReco->nuVertexNoSCC().X()<<" "<< pReco->nuVertexNoSCC().Y()<<" "<< pReco->nuVertexNoSCC().Z() \
                                                                     << " vs " << pRecoUbcc1pi->nuVertexNoSCC().X()<<" "<<pRecoUbcc1pi->nuVertexNoSCC().Y()<<" "<<pRecoUbcc1pi->nuVertexNoSCC().Z()<<std::endl;

            COMPARETWOSINGLEVALUES(Reco, event.reco, pEvent->reco, flashChi2)
            // ******************************************************************************************
            std::cout<<"\n~~~CHECKING RECO PARTICLES~~~"<<std::endl;
            // ******************************************************************************************
            if(event.reco.particles.size()!=pEvent->reco.particles.size())
            {
                std::cout<<"ERROR: Number of reco particles in PeLEE and ubcc1pi events are different: "<<event.reco.particles.size()<<" vs "<<pEvent->reco.particles.size()<<std::endl;
                // const auto particle = event.reco.particles.back();
                // CHECKANYVALUE(Reco Particles Test, particle, pdgCode)
                // std::cout<<"pdgCode: "<<particle.pdgCode()<<std::endl;
                // CHECKANYVALUE(Reco Particles Test, particle, nHitsU)
                // std::cout<<"nHitsU: "<<particle.nHitsU()<<std::endl;
                // CHECKANYVALUE(Reco Particles Test, particle, nHitsV)
                // std::cout<<"nHitsV: "<<particle.nHitsV()<<std::endl;
                // CHECKANYVALUE(Reco Particles Test, particle, nHitsW)
                // std::cout<<"nHitsW: "<<particle.nHitsW()<<std::endl;
                // CHECKANYVALUE(Reco Particles Test, particle, nDaughters)
                // std::cout<<"nDaughters: "<<particle.nDaughters()<<std::endl;
                // CHECKANYVALUE(Reco Particles Test, particle, nDescendents)
                // std::cout<<"nDescendents: "<<particle.nDescendents()<<std::endl;
                // CHECKANYVALUE(Reco Particles Test, particle, trackScore)
                // std::cout<<"trackScore: "<<particle.trackScore()<<std::endl;
                // CHECKANYVALUE(Reco Particles Test, particle, startX)
                // CHECKANYVALUE(Reco Particles Test, particle, startY)
                // CHECKANYVALUE(Reco Particles Test, particle, startZ)
                // CHECKANYVALUE(Reco Particles Test, particle, endX)
                // CHECKANYVALUE(Reco Particles Test, particle, endY)
                // CHECKANYVALUE(Reco Particles Test, particle, endZ)
                // CHECKANYVALUE(Reco Particles Test, particle, directionX)
                // CHECKANYVALUE(Reco Particles Test, particle, directionY)
                // CHECKANYVALUE(Reco Particles Test, particle, directionZ)
                // CHECKANYVALUE(Reco Particles Test, particle, yzAngle)
                // CHECKANYVALUE(Reco Particles Test, particle, xyAngle)
                // CHECKANYVALUE(Reco Particles Test, particle, xzAngle)
                // CHECKANYVALUE(Reco Particles Test, particle, range)
                // CHECKANYVALUE(Reco Particles Test, particle, wiggliness)
                // CHECKANYVALUE(Reco Particles Test, particle, nSpacePointsNearEnd)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardMuon)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardMuon)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardPion)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardPion)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardProton)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardProton)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodMIP)
                // CHECKANYVALUE(Reco Particles Test, particle, truncatedMeandEdxU)
                // CHECKANYVALUE(Reco Particles Test, particle, truncatedMeandEdxV)
                // CHECKANYVALUE(Reco Particles Test, particle, truncatedMeandEdxW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardProtonW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardMuonW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardPionW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardMuonW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardProtonU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardMuonU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardPionU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardMuonU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardProtonV)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardMuonV)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardPionV)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardMuonV)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardProtonW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardMuonW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardPionW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardMuonW)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardProtonU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardMuonU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardPionU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardMuonU)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardProtonV)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardMuonV)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodBackwardPionV)
                // CHECKANYVALUE(Reco Particles Test, particle, likelihoodForwardMuonV)

                // throw std::logic_error("Unequal number of reco particles in PeLEE and ubcc1pi events!");
            }

            for(unsigned long j = 0; j < std::max(event.reco.particles.size(), pEvent->reco.particles.size()); j++)
            {
                const auto j1 = std::min(event.reco.particles.size()-1, j);
                const auto j2 = std::min(pEvent->reco.particles.size()-1, j);
                std::cout<<"-----"<<j<<"("<<j1<<" vs "<<j2<<")"<<"-----"<<std::endl;
                const auto particle = event.reco.particles.at(j1);
                const auto particleUbcc1pi = pEvent->reco.particles.at(j2);

                // isCCInclusiveMuonCandidate
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, pdgCode)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, nHitsU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, nHitsV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, nHitsW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, nDaughters)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, nDescendents)
                // nDescendentHitsU
                // nDescendentHitsV
                // nDescendentHitsW
                // nHitsInLargestDescendent
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, trackScore)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, startX)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, startY)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, startZ)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, endX)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, endY)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, endZ)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, directionX)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, directionY)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, directionZ)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, yzAngle)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, xyAngle)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, xzAngle)
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, length)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, range)
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, transverseVertexDist)
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, longitudinalVertexDist)
                // mcsMomentumForwardMuon
                // mcsMomentumUncertaintyForwardMuon
                // mcsLogLikelihoodForwardMuon
                // mcsMomentumBackwardMuon
                // mcsMomentumUncertaintyBackwardMuon
                // mcsLogLikelihoodBackwardMuon
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, wiggliness)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, nSpacePointsNearEnd)

                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuon)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuon)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardPionU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardPionV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardPionW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardPion)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardPionU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardPionV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardPionW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardPion)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardProtonU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardProtonV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardProtonW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardProton)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardProtonU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardProtonV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardProtonW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardProton)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodMIPU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodMIPV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodMIPW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodMIP)

                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, truncatedMeandEdxU)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, truncatedMeandEdxV)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, truncatedMeandEdxW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, truncatedMeandEdx)

                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardProtonW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardPionW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardProtonU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardPionU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardProtonV) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonV) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardPionV) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonV) //todo remove

                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardProtonW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardPionW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonW) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardProtonU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardPionU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonU) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardProtonV) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardMuonV) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodBackwardPionV) //todo remove
                // COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, likelihoodForwardMuonV) //todo remove

                // truthMatchPurities
                // truthMatchCompletenesses
                // hasMatchedMCParticle
                // bestMatchedMCParticleIndex

                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, generation)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, chi2ForwardMuonW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, chi2ForwardProtonW)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, distance)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, vertexX)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, vertexY)
                COMPARETWOSINGLEVALUES(Reco Particles, particle, particleUbcc1pi, vertexZ)
            }

            // ******************************************************************************************
            // Check the selections
            // ******************************************************************************************
            const auto &[passedGoldenSelectionDefault, cutsPassedDefault, assignedPdgCodesDefault] = peleeSelection.Execute(pEvent);
            const auto passedCCInclusiveSelectionDefault = SelectionHelper::IsCutPassed(cutsPassedDefault, "topologicalScoreCC");
            const auto &[passedGoldenSelectionOriginal, cutsPassedOriginal, assignedPdgCodesOriginal] = ubcc1piSelection.Execute(pEvent);
            const auto passedCCInclusiveSelectionOriginal = SelectionHelper::IsCutPassed(cutsPassedOriginal, "passesCCInclusive");

            std::cout<<"Ubcc1pi event - passes CC-inclusive: ubcc1pi: "<<passedCCInclusiveSelectionOriginal<<" vs pelee: "<<passedCCInclusiveSelectionDefault<<std::endl;
            std::cout<<"Ubcc1pi event - passes golden selection: ubcc1pi: "<<passedGoldenSelectionOriginal<<" vs pelee: "<<passedGoldenSelectionDefault<<std::endl;
            // std::cout<<"    passed default cuts:";
            // for(const auto &cut : cutsPassedDefault) std::cout<<" "<<cut;
            // std::cout<<std::endl;

        }
        break;

        // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // // auto pFileIn = new TFile(fileName, "READ");
        // // auto pFileIn = std::make_shared<TFile>(fileName.c_str(), "READ"); // std::unique_ptr<TFile>(myFile(TFile::Open(fileName, "READ")));
        // // auto pTreeEvents = (TTree*)pFileIn->Get("events");
        // // auto pTreeSubruns = (TTree*)pFileIn->Get("subruns");
        // // // pFileIn->Close();
        // // auto pFileOut = std::make_shared<TFile>("/uboone/data/users/jdetje/ubcc1pi/ubcc1pi_stage2/test.root", "RECREATE");
        // // pFileOut->WriteObject(pTreeEvents, "events");
        // // pFileOut->WriteObject(pTreeSubruns, "subruns");
        // // auto pTree = (TTree*)pFileOut->Get("events");
        // // pFileIn->Close();

        // // std::filesystem::copy(fileName, "/uboone/data/users/jdetje/ubcc1pi/ubcc1pi_stage2/test.root", std::filesystem::copy_options::overwrite_existing); // copy file

        // const std::string path = "/uboone/data/users/jdetje/ubcc1pi/ubcc1pi_stage2/test.root";
        //          // If there is no accessible file, copy first; avoids copying when using this macro consecutively
        //  // ATTN: This does not check whether file is corrupted or not; if in doubt, clear output directory before running
        // // auto pFileOut = std::make_shared<TFile>(path.c_str(), "UPDATE");
        //         // if(!pFileOut || pFileOut->IsZombie())//
        // if(access( path.c_str(), F_OK ) == -1)
        // {
        //     // pFileOut->Close();
        //     auto pFileIn = std::make_shared<TFile>(fileName.c_str(), "READ"); // std::unique_ptr<TFile>(myFile(TFile::Open(fileName, "READ")));
        //     pFileIn->Cp(path.c_str());
        //     pFileIn->Close();
        //     // pFileOut = std::make_shared<TFile>(path.c_str(), "UPDATE");
        // }


        // auto pFileOut = std::make_shared<TFile>(path.c_str(), "UPDATE");
        // auto pTreeFriend = (TTree*)pFileOut->Get("treeFriend");
        // auto pTree = (TTree*)pFileOut->Get("events");
        // if(!pTreeFriend)
        // {
        //     pTree->RemoveFriend(pTreeFriend);
        //     pTreeFriend->Delete();
        // }

        // // bool o_isSelectedGenericCC0Pi, o_isSelectedGenericCC1Pi, o_isSelectedGoldenCC0Pi, o_isSelectedGoldenCC1Pi;
        // // pTree->Branch("isSelectedGenericCC0Pi", &o_isSelectedGenericCC0Pi);
        // // pTree->Branch("isSelectedGenericCC1Pi", &o_isSelectedGenericCC1Pi);
        // // pTree->Branch("isSelectedGoldenCC0Pi", &o_isSelectedGoldenCC0Pi);
        // // pTree->Branch("isSelectedGoldenCC1Pi", &o_isSelectedGoldenCC1Pi);

        // // bool o_isTrueCC0Pi, o_isTrueCC1Pi, o_isSignalCC0Pi, o_isSignalCC1Pi;
        // // if(isOverlay || isDetVar || isNuWro) // Only MC
        // // {
        // //     pTree->Branch("isTrueCC0Pi", &o_isTrueCC0Pi);
        // //     pTree->Branch("isTrueCC1Pi", &o_isTrueCC1Pi);
        // //     pTree->Branch("isSignalCC0Pi", &o_isSignalCC0Pi);
        // //     pTree->Branch("isSignalCC1Pi", &o_isSignalCC1Pi);
        // // }

        // std::vector<Bool_t> isSelectedGenericCC0PiVect(nEvents, false); //todo remove default value -  only for debugging with >100% of events!!!!!!!!!!!!!!!!!!
        // std::vector<Bool_t> isSelectedGenericCC1PiVect(nEvents, false);
        // std::vector<Bool_t> isSelectedGoldenCC0PiVect(nEvents, false);
        // std::vector<Bool_t> isSelectedGoldenCC1PiVect(nEvents, false);
        // std::vector<Bool_t> isTrueCC0PiVect(nEvents, false);
        // std::vector<Bool_t> isTrueCC1PiVect(nEvents, false);
        // std::vector<Bool_t> isSignalCC0PiVect(nEvents, false);
        // std::vector<Bool_t> isSignalCC1PiVect(nEvents, false);

        // for(unsigned h = 2*nEvents/3; h < nEvents; ++h)
        // {
        //     isSelectedGenericCC0PiVect.at(h) = true;
        //     isSignalCC1PiVect.at(h) = true;
        // }

        // auto b1 = pTree->Branch("isSignalCC1Pi", &isSignalCC1PiVect, "isSignalCC1Pi/O");


        // for (unsigned int i = 0; i < 5000/*nEvents*/; ++i) // todo change back to all events!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // {
        //     AnalysisHelper::PrintLoadingBar(i, int(nEvents));
        //     reader.LoadEvent(i);

        //     // ############################################################
        //     // ############### Run the sideband selection #################
        //     // ############################################################
        //     const auto &[passedGoldenSidebandSelection, cutsPassed, assignedPdgCodes] = sidebandSelection.Execute(pEvent);

        //     const auto passedGenericSidebandSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

        //     // o_isSelectedGenericCC0Pi = passedGenericSidebandSelection;
        //     // o_isSelectedGoldenCC0Pi =  passedGoldenSidebandSelection;
        //     isSelectedGenericCC0PiVect.at(i) = passedGenericSidebandSelection;
        //     isSelectedGoldenCC0PiVect.at(i) = passedGoldenSidebandSelection;

        //     // Get the reco analysis data (if available, otherwise set to dummy values)
        //     const auto recoSidebandData = (
        //         passedGenericSidebandSelection
        //             ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, passedGoldenSidebandSelection)
        //             : AnalysisHelper::GetDummyAnalysisData()
        //     );

        //     // Here we apply reco-level phase-space restrictions
        //     // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
        //     // min/max values supplied in the binning. If so, then reject the event.
        //     auto passesSidebandPhaseSpaceReco = false;
        //     if (passedGenericSidebandSelection)
        //     {
        //         // Start by assuming the event passes the phase-space cuts
        //         passesSidebandPhaseSpaceReco = true;

        //         // Check the value of the kinematic quantities are within the phase-space limits
        //         for (const auto &[name, minMax] : phaseSpaceMap)
        //         {
        //             const auto &[min, max] = minMax;
        //             const auto value = getSidebandValue.at(name)(recoSidebandData);

        //             if (value < min || value > max)
        //             {
        //                 passesSidebandPhaseSpaceReco = false;
        //                 break;
        //             }
        //         }
        //     }
        //     const auto isSelectedSidebandGeneric = passedGenericSidebandSelection && passesSidebandPhaseSpaceReco;
        //     const auto isSelectedSidebandGolden = passedGoldenSidebandSelection && passesSidebandPhaseSpaceReco;


        //     // ############################################################
        //     // #################### Run the main selection #####################
        //     // ############################################################
        //     const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
        //     const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

        //     // o_isSelectedGenericCC1Pi = passedGenericSelection;
        //     // o_isSelectedGoldenCC1Pi =  passedGoldenSelection;
        //     isSelectedGenericCC1PiVect.at(i) = passedGenericSelection;
        //     isSelectedGoldenCC1PiVect.at(i) = passedGoldenSelection;

        //     // Get the reco analysis data (if available, otherwise set to dummy values)
        //     const auto recoData = (
        //         passedGenericSelection
        //             ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
        //             : AnalysisHelper::GetDummyAnalysisData()
        //     );

        //     // Here we apply reco-level phase-space restrictions
        //     // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
        //     // min/max values supplied in the binning. If so, then reject the event.
        //     auto passesPhaseSpaceReco = false;
        //     if (passedGenericSelection)
        //     {
        //         // Start by assuming the event passes the phase-space cuts
        //         passesPhaseSpaceReco = true;

        //         // Check the value of the kinematic quantities are within the phase-space limits
        //         for (const auto &[name, minMax] : phaseSpaceMap)
        //         {
        //             const auto &[min, max] = minMax;
        //             const auto value = getValue.at(name)(recoData);

        //             if (value < min || value > max)
        //             {
        //                 passesPhaseSpaceReco = false;
        //                 break;
        //             }
        //         }
        //     }
        //     const auto isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
        //     const auto isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;


        //     if (isDataBNB) continue; // For BNB data that's all we need to do!

        //     // Determine if this is truly a CC0Pi event
        //     const auto isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
        //     // o_isTrueCC0Pi = isTrueCC0Pi;
        //     isTrueCC0PiVect.at(i) = isTrueCC0Pi;

        //     // Get the truth analysis data (if available, otherwise set to dummy values)
        //     const auto truthSidebandData = (
        //         (isTrueCC0Pi)
        //             ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
        //             : AnalysisHelper::GetDummyAnalysisData()
        //     );

        //     // Here we apply truth-level phase-space restrictions
        //     // For all true CC0Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
        //     // event is not classed as "signal"
        //     bool passesSidebandPhaseSpaceTruth = false;
        //     if (isTrueCC0Pi)
        //     {
        //         // Start by assuming the event passes the phase-space cuts
        //         passesSidebandPhaseSpaceTruth = true;

        //         // Check the value of the kinematic quantities are within the phase-space limits
        //         for (const auto &[name, minMax] : phaseSpaceMap)
        //         {
        //             // if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection // todo check this
        //             const auto &[min, max] = minMax;
        //             const auto value = getSidebandValue.at(name)(truthSidebandData);

        //             if (value < min || value > max)
        //             {
        //                 passesSidebandPhaseSpaceTruth = false;
        //                 break;
        //             }
        //         }
        //     }
        //     const auto isCC0PiSignal = isTrueCC0Pi && passesSidebandPhaseSpaceTruth;
        //     // o_isSignalCC0Pi = isCC0PiSignal;
        //     isSignalCC0PiVect.at(i) = isCC0PiSignal;


        //     const auto isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
        //     // o_isTrueCC1Pi = isTrueCC1Pi;
        //     isTrueCC1PiVect.at(i) = isTrueCC1Pi;

        //     // Get the truth analysis data (if available, otherwise set to dummy values)
        //     const auto truthData = (
        //         isTrueCC1Pi
        //             ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
        //             : AnalysisHelper::GetDummyAnalysisData()
        //     );

        //     // Here we apply truth-level phase-space restrictions
        //     // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
        //     // event is not classed as "signal"
        //     bool passesPhaseSpaceTruth = false;
        //     if (isTrueCC1Pi)
        //     {
        //         // Start by assuming the event passes the phase-space cuts
        //         passesPhaseSpaceTruth = true;

        //         // Check the value of the kinematic quantities are within the phase-space limits
        //         for (const auto &[name, minMax] : phaseSpaceMap)
        //         {
        //             const auto &[min, max] = minMax;
        //             const auto value = getValue.at(name)(truthData);

        //             if (value < min || value > max)
        //             {
        //                 passesPhaseSpaceTruth = false;
        //                 break;
        //             }
        //         }
        //     }

        //     const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruth;
        //     // o_isSignalCC1Pi = isCC1PiSignal;
        //     // isSignalCC1PiVect.at(i) = isCC1PiSignal;
        //     isSignalCC1PiVect.at(i) = true; //todo remove debugging only

        //     // pTree->Fill();
        // } // End of event-level iteration

        // Float_t isSelectedGenericCC0Pi = 44.14f;

        // // auto pB0 = pTree->GetBranch("isSelectedGenericCC0Pi");
        // // if(pB0)
        // // {
        // //     auto pLeaf = pTree->GetLeaf("isSelectedGenericCC0Pi");
        // //     pTree->GetListOfBranches()->Remove(pB0);
        // //     pTree->GetListOfLeaves()->Remove(pLeaf);
        // //     pTree->Write();
        // // }

        // // const auto pB1 = pTree->Branch("isSelectedGenericCC0Pi", &isSelectedGenericCC0Pi, "isSelectedGenericCC0Pi/F");
        // TTree treeFriend("treeFriend", "ubcc1pi");
        // const auto pB1 = treeFriend.Branch("isSelectedGenericCC0Pi", &isSelectedGenericCC0Pi, "isSelectedGenericCC0Pi/F");

        // // // Disable everything...
        // // pTree->SetBranchStatus("*", true);
        // // // ...but the branch we need
        // // pTree->SetBranchStatus("isSelectedGenericCC0Pi", true);
        // // pTree->SetBranchAddress("isSelectedGenericCC0Pi", &isSelectedGenericCC0Pi);
        // Long64_t nentries = pTree->GetEntries();
        // for (Long64_t k = 0; k < nentries; k++) {
        //     if(k%100==0) std::cout<<(1.0*k)/nentries<<std::endl;
        //     // pTree->GetEntry(k);
        //     if(k<nentries/3) isSelectedGenericCC0Pi = 41.3f; //isSelectedGenericCC0PiVect.at(k);
        //     else isSelectedGenericCC0Pi = 41.5f;
        //     pB1->Fill(); // Use Branch->BackFill() with newer ROOT versions
        // }

        // treeFriend.Write("", TObject::kOverwrite);
        // pTree->AddFriend(&treeFriend);
        // pTree->Write("", TObject::kOverwrite);
        // std::cout<<std::endl;

        // // // pTree->Branch("isSelectedGenericCC0Pi", &isSelectedGenericCC0PiVect, "isSelectedGenericCC0Pi/O")->Fill();
        // // pTree->Branch("isSelectedGenericCC1Pi", &isSelectedGenericCC1PiVect, "isSelectedGenericCC1Pi/O")->Fill();
        // // pTree->Branch("isSelectedGoldenCC0Pi", &isSelectedGoldenCC0PiVect, "isSelectedGoldenCC0Pi/O")->Fill();
        // // pTree->Branch("isSelectedGoldenCC1Pi", &isSelectedGoldenCC1PiVect, "isSelectedGoldenCC1Pi/O")->Fill();


        // // for(int l = 0; l < 50; l++)
        // // {
        // //     std::cout<<" "<<isTrueCC0PiVect.at(l);
        // // }
        // // std::cout<<"\n--------------------"<<std::endl;
        // // for(int l = 0; l < 50; l++)
        // // {
        // //     std::cout<<" "<<isTrueCC0PiVect.at(nEvents-1-l);
        // // }
        // // std::cout<<std::endl;

        // // if(isOverlay || isDetVar || isNuWro) // Only MC
        // // {
        // //     pTree->Branch("isTrueCC0Pi", &isTrueCC0PiVect, "isTrueCC0Pi/O")->Fill();
        // //     pTree->Branch("isTrueCC1Pi", &isTrueCC1PiVect, "isTrueCC1Pi/O")->Fill();
        // //     pTree->Branch("isSignalCC0Pi", &isSignalCC0PiVect, "isSignalCC0Pi/O")->Fill();
        // //     pTree->Branch("isSignalCC1Pi", &isSignalCC1PiVect, "isSignalCC1Pi/O")->Fill();
        // //     // b1->Fill();
        // // }

        // // pTree->Fill();

        // //////////////////////////////////////////////////////////////
        // // Add values to input ntuples and save them as new files
        // //////////////////////////////////////////////////////////////
        // pFileOut->Write("", TObject::kOverwrite);
        // pFileOut->Close();
        // break; //TODO Remove after first tests !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    } // End of file-level iterration


    // truth_particle_pdgCode_plot.SaveAs("truth_particle_pdgCode");
    // truth_particle_endX_plot.SaveAs("truth_particle_endX");
    // truth_particle_endY_plot.SaveAs("truth_particle_endY");
    // truth_particle_endZ_plot.SaveAs("truth_particle_endZ");
    // truth_particle_momentumX_plot.SaveAs("truth_particle_momentumX");
    // truth_particle_momentumY_plot.SaveAs("truth_particle_momentumY");
    // truth_particle_momentumZ_plot.SaveAs("truth_particle_momentumZ");
    // // truth_particle_momentum_plot
    // truth_particle_energy_plot.SaveAs("truth_particle_energy");

    std::cout<<"-----------------Done-----------------"<<std::endl;


}

} // namespace ubcc1pi_macros
