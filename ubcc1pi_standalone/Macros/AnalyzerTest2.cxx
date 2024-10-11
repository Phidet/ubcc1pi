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

void AnalyzerTest2(const Config &config)
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
    // auto ubcc1piSelection = SelectionHelper::GetCC0piSelection();
    // auto peleeSelection = SelectionHelper::GetCC0piSelectionOriginalPeLEE();

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
            Event event(*pEventPeLEE, true);
            // event.Print();
            std::cout << "DEBUG Point Y10" << std::endl;

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

            for(unsigned long j = 0; j < std::max(event.reco.particles.size(), pEvent->reco.particles.size()); j++)
            {
                const auto j1 = std::min(event.reco.particles.size()-1, j);
                const auto j2 = std::min(pEvent->reco.particles.size()-1, j);
                std::cout<<"-----"<<j<<"("<<j1<<" vs "<<j2<<")"<<"-----"<<std::endl;
                const auto particle = event.reco.particles.at(j1);
                const auto particleUbcc1pi = pEvent->reco.particles.at(j2);

                std::cout << (particle.mcsMuonMomentum.IsSet() ? std::to_string(particle.mcsMuonMomentum()) : "Not set") << " ";
                std::cout << (particleUbcc1pi.mcsMomentumForwardMuon.IsSet() ? std::to_string(particleUbcc1pi.mcsMomentumForwardMuon()) : "Not set") << " ";
                std::cout << (particleUbcc1pi.mcsMomentumBackwardMuon.IsSet() ? std::to_string(particleUbcc1pi.mcsMomentumBackwardMuon()) : "Not set") << " " << std::endl;            }

        }
        break;

    } // End of file-level iterration

    std::cout<<"-----------------Done-----------------"<<std::endl;


}

} // namespace ubcc1pi_macros
