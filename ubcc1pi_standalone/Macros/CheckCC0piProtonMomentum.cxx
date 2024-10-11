/**
 *  @file  ubcc1pi_standalone/Macros/CheckCC0piProtonMomentum.cxx
 *
 *  @brief The implementation file of the CheckCC0piProtonMomentum macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubsmear.h"


#include <fstream> //todo remove

// Boost libraries
#include "binary_iarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

const auto testingFractionXSecCheckProtonMom = 0.2f;


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void CheckCC0piProtonMomentum(const Config &config)
{

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the details of the systematic parameters to apply
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we read in the "dimensions" of the systematic parameters to apply. For multisim parameters (flux & xsec), this is a map from
    // each parameter name to the expected number of universes. For unisim parameters (detector variations), this is a map from the
    // identidier (i.e. name) of the detector variation sample, to the identifier of the corresponding central-value sample. Additionally
    // we set the number of bootstrap universes (for the MC stat uncertainty), the corresponding weights are generated for each event.
    CrossSectionHelper::CrossSection::SystParams systParams;
    systParams.fluxDimensions = config.extractXSecs.fluxDimensions;
    systParams.xsecDimensions = config.extractXSecs.xsecDimensions;
    systParams.reintDimensions = config.extractXSecs.reintDimensions;
    systParams.detVarDimensions = config.extractXSecs.detVarDimensions;
    systParams.potFracUncertainty = config.extractXSecs.potFracUncertainty;
    systParams.nBootstrapUniverses = config.extractXSecs.nBootstrapUniverses;
    systParams.nSidebandFitUniverses = config.extractXSecs.nSidebandFitUniverses;

    CrossSectionHelper::SystDimensionsMap emptySystDimMap = {};
    auto systParamsNuWro = systParams; // shallow copy sufficient since struct contains no pointers
    systParamsNuWro.potFracUncertainty = 0.0;
    systParamsNuWro.fluxDimensions = emptySystDimMap;
    systParamsNuWro.reintDimensions = emptySystDimMap;

    auto systParamsGenie = systParamsNuWro; // shallow copy sufficient since struct contains no pointers
    systParamsGenie.xsecDimensions = emptySystDimMap;
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup an object that holds the information about how we should scale an event rate to obtain a cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we specify the:
    // - Flux             [10^-10 cm^-2 POT^-1]
    // - Exposure POT     [10^20 POT]               (stored in config as [POT])
    // - Target density   [10^31 nucleons/cm^3]     (stored in config as [10^23 nucleons/cm^3])
    // - Fiducial volume  [cm^3]
    //
    // Hence the units of the eventual cross-seciton are:
    // - Cross-section    [Flux * Exposure * Target density * Fiducial volume]^-1 = [10^-41 cm^2 / nucleon]
    //
    // Here we use a FluxReweightor to specify the flux. For each universe of each flux systematic paramter, this is uses the ratio of the
    // total neutrino event rate (as a function of the true neutrino energy) to the same rate in the nominal simulation to reweight the
    // input neutrino flux distribution. The integrated flux is each universe is used to scale the selected event rate when calculating the
    // cross-section in that universe. For all non-flux parameters, the nominal integrated flux is used.

    // const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
    // const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);
    // std::map<std::string,std::map<std::map<std::string, CrossSectionHelper::CrossSection::ScalingData>>> scalingDataScaledMap;
    // for (const auto &dataTypeName : {"NuWro", "BNB"})
    // {
    //     for (const auto &selectionName : {"generic", "golden"})
    //     {
    //         for (const auto &name : {...})
    //         {
    //             scalingDataScaledMap[dataTypeName].emplace(selectionName, CrossSectionHelper::CrossSection::ScalingData);
    //             scalingDataScaledMap.at(dataTypeName).at(selectionName).pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
    //         }
    //     }
    // }

    // std::cout << "- Flux:            " << scalingData.pFluxReweightor->GetIntegratedNominalFlux() << " * 10^-10 cm^-2 POT^-1" << std::endl;
    // std::cout << "- Exposure:        " << scalingData.exposurePOT << " * 10^20 POT" << std::endl;
    // std::cout << "- Target density:  " << config.global.targetDensity << " * 10^31 nucleons/cm^3" << std::endl;
    // std::cout << "- Fiducial volume: " << AnalysisHelper::GetFiducialVolume() << " * cm^3" << std::endl;
    // std::cout << "- nTargets:        " << scalingData.nTargets << " * 10^31 nucleons" << std::endl;

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the event selection
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we are using the default CC1pi selection. The final cut of this selection (named "likelyGoldenPion") is used to identify events
    // in which the reconstructed & selected pion candidate is "golden". Here a golden pion is one that is: contained within the TPC,
    // doesn't undergo any scatters on the argon, and comes to rest (before any secondary interaction). If an event passes all the cuts of
    // this selection, it's said to pass the "golden" selection. If an event passes the penultimate cut, it is said to have passed the
    // "generic" selection. Events passing the generic selection are likely to be CC1pi, but may or may not have a golden pion. All events
    // that pass the golden selection also pass the generic selection - but not all events that pass the generic selection also pass the
    // golden selection.
    auto selection = SelectionHelper::GetDefaultSelection();
    // auto sidebandSelection = SelectionHelper::GetSelection("CC0pi");

    // We additionally make a map from each cross-section to the limits of the phase-space that we should consider. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;

    // ATTN the configuration allows the user to enable or disable each cross-section. If a cross-section has been disabled, then it won't
    // be added the xsecMapNuWro. However, the phaseSpaceMap always includes all kinematic paramters!

    std::map<std::string, std::map<std::string, std::map<std::string, CrossSectionHelper::CrossSection::ScalingData>>> scalingDataScaledMap;
    const auto fluxHistNames = CrossSectionHelper::GetNominalFluxHistNames(config.flux.nuPdgsSignal, config.flux.nuPdgToHistName, config.flux.nomHistPattern);
    const auto &[fluxBinEdges, fluxValues] = CrossSectionHelper::ReadNominalFlux(config.flux.fileName, fluxHistNames, config.flux.pot);
    CrossSectionHelper::CrossSection::ScalingData scalingDataUnscaled;
    // CrossSectionHelper::CrossSection::ScalingData scalingDataNuWroTruth;
    scalingDataUnscaled.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
    // scalingDataNuWroTruth.pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);

    // The names of the cross-section kinematic parameters, and their sideband binning information.
    std::map<std::string, std::vector<float>> sidebandBinningMap {
        { "muonCosTheta",  config.global.muonCosTheta.binEdges           },
        { "muonPhi",       config.global.muonPhi.binEdges                },
        { "muonMomentum",  config.global.sidebandMuonMomentum.binEdges   },
        { "pionCosTheta",  config.global.pionCosTheta.binEdges           },
        { "pionPhi",       config.global.pionPhi.binEdges                },
        { "pionMomentum",  config.global.sidebandPionMomentum.binEdges   },
        { "muonPionAngle", config.global.muonPionAngle.binEdges          },
        { "nProtons",      config.global.nProtons.binEdges               }};

    // Add the differential cross-sections
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
    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
    // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section and for the sideband
    // -------------------------------------------------------------------------------------------------------------------------------------
    ExtractionHelper::AnalysisValueMap getValue;
    ExtractionHelper::AnalysisValueMap getSidebandValue;
    ExtractionHelper::PopulateAnalysisValueMap(getValue, false);
    ExtractionHelper::PopulateAnalysisValueMap(getSidebandValue, true); // true: creates sideband getters

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the input files
    // -------------------------------------------------------------------------------------------------------------------------------------
    ExtractionHelper::InputFileList inputData;
    float totalExposurePOT;
    ExtractionHelper::PopulateInputFileList(config, inputData, totalExposurePOT);


    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events for CC1pi
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);

        if(!isOverlay) continue;

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        if (isOverlay || isNuWro) //todo check if nuwro is needed here
            reader.EnableSystematicBranches();

        auto pEvent = reader.GetBoundEventAddress();
        const auto nEvents = reader.GetNumberOfEvents();

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Loop over the events in the file
        std::cout<<"############################\nUsing "<<testingFractionXSecCheckProtonMom*100<<"\% of events!\n############################"<<std::endl;
        for (unsigned int i = 0; i < nEvents*testingFractionXSecCheckProtonMom; ++i) // todo change back to every event!!!!!!!!!!!!!!!!!!!!!!
        {
            AnalysisHelper::PrintLoadingBar(i, int(nEvents*testingFractionXSecCheckProtonMom));
            reader.LoadEvent(i);

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------
            // Run the selection
            const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedGenericSelection
                    ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            bool passesPhaseSpaceReco = false;
            if (passedGenericSelection)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceReco = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValue.at(name)(recoData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceReco = false;
                        break;
                    }
                }
            }
            const auto isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
            const auto isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;
            std::map<std::string, bool> isSelectedMap = {{"generic",isSelectedGeneric},{"golden",isSelectedGolden}};

            if(!isSelectedGolden) continue;

            // Determine if this is truly a CC1Pi event
            const auto isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto truthData = (
                isTrueCC1Pi
                    ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // Here we apply truth-level phase-space restrictions
            // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
            // event is not classed as "signal"
            bool passesPhaseSpaceTruth = false;
            if (isTrueCC1Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceTruth = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValue.at(name)(truthData);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruth;
            const auto isSignal = isCC1PiSignal;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event is signal, and apply any phase-space restrictions based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Determine if this is truly a CC0Pi event
            const auto isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);


            // Get the truth analysis data (if available, otherwise set to dummy values)
            const auto sidebandTruthData = (
                isTrueCC0Pi
                    ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            bool passesSidebandPhaseSpaceTruth = false;
            if (isTrueCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesSidebandPhaseSpaceTruth = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getSidebandValue.at(name)(sidebandTruthData);

                    if (value < min || value > max)
                    {
                        if (name == "pionMomentum") continue; // todo: Not really compatible with the pion momentum in the CC1pi selection
                        passesSidebandPhaseSpaceTruth = false;
                        break;
                    }
                }
            }

            const auto isCC0PiSignal = isTrueCC0Pi && passesSidebandPhaseSpaceTruth;
            if(!isCC0PiSignal) continue;

            const auto name = "pionMomentum";
            const auto recoValue = getValue.at(name)(recoData);
            const auto trueValue = getValue.at(name)(truthData);

            // Get the index of the proton in the truth particles that matches the pion in the reco particles
            if (assignedPdgCodes.size() != pEvent->reco.particles.size())
                throw std::logic_error("AnalysisHelper::GetRecoAnalysisDataCC0pi - The assigned PDG codes is the wrong size");
            const auto recoPionCandidate = pEvent->reco.particles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211));
            try
            {
                const auto protonAsPionTruthIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoPionCandidate, pEvent->truth.particles);
                const auto leadingProtonTruthIndex = AnalysisHelper::GetTrueLeadingProtonIndex(pEvent->truth, true, config.global.protonMomentumThreshold);
                if(protonAsPionTruthIndex == leadingProtonTruthIndex)
                {
                    std::cout<<"Yes!  ";
                }
                else
                {
                    std::cout<<"No!  ";
                }
                std::cout<<protonAsPionTruthIndex<<"("<<pEvent->truth.particles.at(protonAsPionTruthIndex).pdgCode()<<"/"<<pEvent->truth.particles.at(protonAsPionTruthIndex).momentum()<<")"<<" vs "<<leadingProtonTruthIndex<<"("<<pEvent->truth.particles.at(leadingProtonTruthIndex).momentum()<<")"<<std::endl;
            }
            catch (const std::exception &)
            {
                std::cout<<"No match found for proton as pion"<<std::endl;
                continue;
            }
        }
    }

    std::cout<<"------------- All Done -------------"<<std::endl;
}

} // namespace ubcc1pi_macros

