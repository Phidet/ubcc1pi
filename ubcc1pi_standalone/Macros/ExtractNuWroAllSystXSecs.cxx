/**
 *  @file  ubcc1pi_standalone/Macros/ExtractNuWroAllSystXSecs.cxx
 *
 *  @brief The implementation file of the ExtractNuWroAllSystXSecs macro
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

const auto testingFractionXSecSystAll = 1.f;


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void ExtractNuWroAllSystXSecs(const Config &config)
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
    // systParamsNuWro.potFracUncertainty = 0.0;
    // systParamsNuWro.fluxDimensions = emptySystDimMap;
    // systParamsNuWro.reintDimensions = emptySystDimMap;

    auto systParamsNuWroTruth = systParams; // shallow copy sufficient since struct contains no pointers
    systParamsNuWroTruth.potFracUncertainty = 0.0;
    systParamsNuWroTruth.fluxDimensions = emptySystDimMap;
    systParamsNuWroTruth.reintDimensions = emptySystDimMap;

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

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the cross-section objects
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Here we make a map from a name of the cross-section to the cross-section object itself. In this way, we can iterate through the
    // cross-section objects and reduce code-bloat. The first index is an identifier for the selection that's applied (generic or goldlen),
    // the second index is an identifier for the kinematic quantity that's relevant for the cross-section (e.g. muonMomentum), and the
    // mapped type is the cross-section object.
    typedef std::map<std::string, CrossSectionHelper::CrossSection> CrossSectionMap;
    typedef std::map<std::string, CrossSectionMap > NestedCrossSectionMap;

    NestedCrossSectionMap xsecMapBNBUnscaled;
    NestedCrossSectionMap xsecMapBNBScaled;
    NestedCrossSectionMap xsecMapNuWroScaled;
    NestedCrossSectionMap xsecMapNuWroUnscaled;
    NestedCrossSectionMap xsecMapGenieScaled;
    NestedCrossSectionMap xsecMapGenieUnscaled;
    NestedCrossSectionMap xsecMapNuWroTrue; // True NuWro cross-sections

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

    // std::cout<<"DEBUG Y-6"<<std::endl;
    // The names of the cross-section kinematic parameters, and their sideband binning information.
    std::map<std::string, std::vector<float>> sidebandBinningMap {
        { "total",         {-1.f, 1.f}                                   },
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

        // Here we calculate every cross-section using both the generic and golden selection. In the end we only use the golden selection for
        // the pion momentum, but we additionally apply it to other cross-sections as a cross-check.
        // Add the cross-section object to the map using the binning from the input configuration
        const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        for (const auto &selectionName : {"generic", "golden"})
        {
            // Don't setup a cross-section object if it's been disabled in the configuration
            if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at(name))
            {
                for (const auto &dataTypeName : {"NuWro", "BNB", "Genie"})
                {
                    const auto pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
                    scalingDataScaledMap[dataTypeName][selectionName].emplace(name, CrossSectionHelper::CrossSection::ScalingData(pFluxReweightor));
                    // scalingDataScaledMap.at(dataTypeName).at(selectionName).at(name).pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
                }
                xsecMapBNBUnscaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapBNBScaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParams, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapNuWroScaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsNuWro, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapNuWroUnscaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsNuWro, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapGenieScaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsGenie, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapGenieUnscaled[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsGenie, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
                xsecMapNuWroTrue[selectionName].emplace(name, CrossSectionHelper::CrossSection(systParamsNuWroTruth, extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
            }
        }
    }
    // std::cout<<"DEBUG Y-5"<<std::endl;
    // ATTN here we use the machinary for a differential cross-section, and treat the total cross-section as a single-bin measurement.
    // The "kinematic quantity" in this case is just a dummy parameter. Here we define a single bin with edges arbitrarily chosen to be
    // (-1 -> +1), and we request that the cross-section object does not apply bin-width scaling. When we fill this object, we will use
    // the dummy kinematic quantity with a value of 0 for all events. This is arbitrary, as long as it's within the bin edges we chose.
    // In this way the single bin contains all events. It's just a trick to avoid implementing extra logic for the total cross-section.
    // std::cout<<"DEBUG Y-4.6"<<std::endl;
    for (const auto &selectionName : {"generic", "golden"})
    {
        // std::cout<<"DEBUG Y-4.5"<<std::endl;
        // Don't setup a cross-section object if it's been disabled in the configuration
        if (config.extractXSecs.crossSectionIsEnabled.at(selectionName).at("total"))
        {
            // std::cout<<"DEBUG Y-4.4"<<std::endl;
            for (const auto &dataTypeName : {"NuWro", "BNB", "Genie"})
            {
                // std::cout<<"DEBUG Y-4.3"<<std::endl;
                const auto pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
                scalingDataScaledMap[dataTypeName][selectionName].emplace("total", CrossSectionHelper::CrossSection::ScalingData(pFluxReweightor));
                // scalingDataScaledMap.at(dataTypeName).at(selectionName).at("total").pFluxReweightor = std::make_shared<CrossSectionHelper::FluxReweightor>(fluxBinEdges, fluxValues, systParams.fluxDimensions);
            }
            // std::cout<<"DEBUG Y-4.2"<<std::endl;
            xsecMapBNBUnscaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapBNBScaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParams, {-1.f, 1.f}, false, false, false));
            xsecMapNuWroScaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsNuWro, {-1.f, 1.f}, false, false, false));
            xsecMapNuWroUnscaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsNuWro, {-1.f, 1.f}, false, false, false));
            xsecMapGenieScaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsGenie, {-1.f, 1.f}, false, false, false));
            xsecMapGenieUnscaled[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsGenie, {-1.f, 1.f}, false, false, false));
            xsecMapNuWroTrue[selectionName].emplace("total", CrossSectionHelper::CrossSection(systParamsNuWroTruth, {-1.f, 1.f}, false, false, false));
            // std::cout<<"DEBUG Y-4.1"<<std::endl;
        }
    }
    // std::cout<<"DEBUG Y-4"<<std::endl;

    // Check to see if we have any cross-sections enabled
    if (xsecMapNuWroScaled.empty())
    {
        std::cout << "All cross-sections have been disabled in the configuration! Nothing more to do" << std::endl;
        return;
    }

    // Print the names of the cross-sections we are going to extract
    std::cout << "The following cross-sections are enabled:" << std::endl;
    for (const auto &[selectionName, xsecs] : xsecMapNuWroScaled)
    {
        std::cout << " - Selection: " << selectionName << std::endl;
        for (const auto &entry : xsecs)
        {
            const auto &name = entry.first;
            std::cout << " - " << name << std::endl;
        }
    }

    // std::cout<<"DEBUG Y-3"<<std::endl;

    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Load the sideband weights from file
    // // -------------------------------------------------------------------------------------------------------------------------------------
    // // Loop over all cross-section objects

    typedef std::pair<std::vector<Double_t>,std::vector<Double_t>> paramAndErrorPair; // Todo: Improve code!
    //Parameters: dataTypeName selectionName name
    typedef std::map<std::string, std::map<std::string, std::map<std::string, paramAndErrorPair>>> nominalFitMap;
    //Parameters: dataTypeName selectionName name paramName (i.e. golden muonMomentum hadronProduction)
    typedef std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<paramAndErrorPair>>>>> universeFitMap;

    //Parameters: selectionName, name;
    nominalFitMap cc0piNominalConstraintMap;
    universeFitMap cc0piUniverseConstraintMap;

    std::ifstream ifs1("cc0piNominalConstraintMapAllSyst.bin", std::ios::binary);
    std::ifstream ifs2("cc0piUniverseConstraintMapAllSyst.bin", std::ios::binary);
    boost::archive::binary_iarchive iarch1(ifs1);
    boost::archive::binary_iarchive iarch2(ifs2);
    iarch1 >> cc0piNominalConstraintMap;
    iarch2 >> cc0piUniverseConstraintMap;
    ifs1.close();
    ifs2.close();

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

    // std::cout<<"DEBUG Y-2"<<std::endl;

    for (auto &[dataTypeName, selectionNameMap] : scalingDataScaledMap)
    {
        for (auto &[selectionName, nameMap] : selectionNameMap)
        {
            for (auto &[name, scalingData]: nameMap)
            {
                scalingData.exposurePOT = totalExposurePOT*testingFractionXSecSystAll; //todo check that this value is also correct for NuWro
                scalingData.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
            }
        }
    }
    scalingDataUnscaled.exposurePOT = totalExposurePOT*testingFractionXSecSystAll;
    scalingDataUnscaled.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();
    // scalingDataNuWroTruth.exposurePOT = totalExposurePOT*testingFractionXSecSystAll; //todo check this is correct and shouldn't use nuwro values
    // scalingDataNuWroTruth.nTargets = config.global.targetDensity * (1e-8) * AnalysisHelper::GetFiducialVolume();

    // std::cout<<"DEBUG Y-1"<<std::endl;
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Count the events for CC1pi
    // -------------------------------------------------------------------------------------------------------------------------------------
    // Loop over the files
    for (const auto &[sampleType, sampleName, fileName, normalisation] : inputData)
    {
        std::cout << "Reading input file: " << fileName << std::endl;

        // std::cout<<"DEBUG Y0"<<std::endl;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);

        // std::cout<<"DEBUG Y1"<<std::endl;
        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader reader(fileName);

        // std::cout<<"DEBUG Y1.1"<<std::endl;
        if (isOverlay || isNuWro) //todo check if nuwro is needed here
            reader.EnableSystematicBranches();

        // std::cout<<"DEBUG Y2"<<std::endl;
        auto pEvent = reader.GetBoundEventAddress();
        const auto nEvents = reader.GetNumberOfEvents();

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Loop over the events in the file
        std::cout<<"############################\nUsing "<<testingFractionXSecSystAll*100<<"\% of events!\n############################"<<std::endl;
        for (unsigned int i = 0; i < nEvents*testingFractionXSecSystAll; ++i) // todo change back to every event!!!!!!!!!!!!!!!!!!!!!!
        {
            // std::cout<<"DEBUG TYT - P 0"<<std::endl;
            AnalysisHelper::PrintLoadingBar(i, int(nEvents*testingFractionXSecSystAll));
            reader.LoadEvent(i);
            // std::cout<<"DEBUG TYT - P 1"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Work out if this event passed the selection and apply any additional phase-space cuts based on the input binning
            // -----------------------------------------------------------------------------------------------------------------------------
            // Run the selection
            const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
            const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

            // std::cout<<"DEBUG TYT - P 2"<<std::endl;
            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoData = (
                passedGenericSelection
                    ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            // std::cout<<"DEBUG TYT - P 3"<<std::endl;
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

            // std::cout<<"DEBUG TYT - P 4"<<std::endl;
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

            // std::cout<<"DEBUG TYT - P 5"<<std::endl;
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

            // std::cout<<"DEBUG TYT - P 6"<<std::endl;
            const auto isCC0PiSignal = isTrueCC0Pi && passesSidebandPhaseSpaceTruth;

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle 'data' (BNB & fake-data)
            // -----------------------------------------------------------------------------------------------------------------------------
            for (auto &[selectionName, xsecs] : xsecMapBNBUnscaled)
            {
                // Determine if we passed the relevant selection
                const auto isSelected = isSelectedMap.at(selectionName);
                // Only count events passing the selection
                if (!isSelected) continue;

                for (auto &[name, xsec] : xsecs)
                {
                    if (isDataBNB)
                    {
                        xsec.AddSelectedBNBDataEvent(getValue.at(name)(recoData)); //todo was getSidebandValue, check getValue is correct
                        xsecMapBNBScaled.at(selectionName).at(name).AddSelectedBNBDataEvent(getValue.at(name)(recoData)); //todo was getSidebandValue, check getValue is correct
                    }
                    else if (isNuWro)
                    {
                        xsecMapNuWroUnscaled.at(selectionName).at(name).AddWeightedSelectedDataEvent(getValue.at(name)(recoData), weight); //todo was getSidebandValue, check getValue is correct
                        xsecMapNuWroScaled.at(selectionName).at(name).AddWeightedSelectedDataEvent(getValue.at(name)(recoData), weight); //todo was getSidebandValue, check getValue is correct
                    }
                    else if (config.global.useGenieAsData && isOverlay)
                    {
                        xsecMapGenieUnscaled.at(selectionName).at(name).AddWeightedSelectedDataEvent(getValue.at(name)(recoData), weight);
                        xsecMapGenieScaled.at(selectionName).at(name).AddWeightedSelectedDataEvent(getValue.at(name)(recoData), weight);
                    }
                }
            }
            if (isDataBNB) continue; // For BNB data that's all we need to do!

            // std::cout<<"DEBUG TYT - P 7"<<std::endl;
            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle the detector variation samples (unisims)
            // -----------------------------------------------------------------------------------------------------------------------------
            if (isDetVar)
            {
                for (auto &[selectionName, xsecs] : xsecMapBNBUnscaled)
                {
                    // Determine if we passed the relevant selection
                    const auto isSelected = isSelectedMap.at(selectionName);

                    for (auto &[name, xsec] : xsecs)
                    {
                        const auto recoValue = getValue.at(name)(recoData);

                        if (isSignal)
                        {
                            // Handle signal events
                            const auto trueValue = getValue.at(name)(truthData);
                            xsec.AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                            xsecMapBNBScaled.at(selectionName).at(name).AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                            xsecMapNuWroScaled.at(selectionName).at(name).AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                            xsecMapNuWroUnscaled.at(selectionName).at(name).AddSignalEventDetVar(recoValue, trueValue, isSelected, weight, sampleName);
                        }
                        else if (isSelected)
                        {
                            // Handle selected background events
                            xsec.AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                            xsecMapBNBScaled.at(selectionName).at(name).AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                            xsecMapNuWroScaled.at(selectionName).at(name).AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                            xsecMapNuWroUnscaled.at(selectionName).at(name).AddSelectedBackgroundEventDetVar(recoValue, weight, sampleName);
                        }
                    }
                }

                // For detector variation samples, that's all we need to do!
                continue;
            }

            // -----------------------------------------------------------------------------------------------------------------------------
            // Handle all other events (i.e those from the nominal simulation): Overlays, dirt, EXT data
            // -----------------------------------------------------------------------------------------------------------------------------
            // Get the flux weights
            auto fluxWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.fluxDimensions, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions)
            );

            // std::cout<<"DEBUG TYT - P 8"<<std::endl;
            // Get the cross-section weights
            // ATTN here we optionally scale the cross-section weights down by the genieTuneEventWeight - this is done so we don't
            // double count this weight (once in the nominal event weight, and once in the xsec systematic event weights)
            // const auto xsecWeightsScaleFactor = 1.f;
            auto xsecWeightsScaleFactor = (isOverlay && config.extractXSecs.scaleXSecWeights) ? pEvent->truth.genieTuneEventWeight() : 1.f;
            // if (xsecWeightsScaleFactor<0.1f && isCC0PiSignal) std::cout<<"xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<< " - event: " << i <<" - true protonMultiplicity: "<< getSidebandValue.at("nProtons")(truthData) << std::endl;
            // xsecWeightsScaleFactor = std::max(xsecWeightsScaleFactor, 0.0001f);

            if (xsecWeightsScaleFactor<=std::numeric_limits<float>::epsilon())// || std::isinf(xsecWeightsScaleFactor))
            {
                if (weight>std::numeric_limits<float>::epsilon())
                {
                    std::cout<<"xsecWeightsScaleFactor is close to zero but nominal weight is not!"<<std::endl;
                    std::cout<<"Skipped event: "<<i<<" - weight: "<<weight<<" - xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<std::endl;
                    throw std::runtime_error("xsecWeightsScaleFactor is close to zero but nominal weight is not!");
                }
                continue;
            }

            auto xsecWeights = (
                isOverlay
                    ? CrossSectionHelper::ScaleWeightsMap(CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.xsecDimensions, config.extractXSecs.mutuallyExclusiveDimensions), xsecWeightsScaleFactor, config.extractXSecs.xsecUBGenieScalingMap)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions)
            );

            // for (auto &[paramName,weightVector] : xsecWeights)
            // {
            //     int i = 0;
            //     // std::cout<<"DEBUG xsecWeights Unscaled SuperPre isCC0PiSignal: "<<isCC0PiSignal<<" - paramName: "<<paramName<<" - weight: "<<weight<<" - ";
            //     for (auto &w : weightVector)
            //     {
            //         std::cout<<" "<<w;
            //         i++;
            //         if(i>10) break;
            //     }
            //     std::cout<<std::endl;
            //     break;
            // }

            // Get the reinteraction weights
            auto reintWeights = (
                isOverlay
                    ? CrossSectionHelper::GetWeightsMap(pEvent->truth, systParams.reintDimensions)//, config.extractXSecs.mutuallyExclusiveDimensions)
                    : CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions)
            );

            // std::cout<<"DEBUG TYT - P 9"<<std::endl;
            // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
            const auto useForFluxReweighting = isOverlay && AnalysisHelper::IsInActiveVolume(pEvent->truth.nuVertex()) &&
                std::find(config.flux.nuPdgsSignal.begin(), config.flux.nuPdgsSignal.end(), pEvent->truth.nuPdgCode()) != config.flux.nuPdgsSignal.end();
            if (useForFluxReweighting) scalingDataUnscaled.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);

            // Unit weights to be used for fake-data studies (NuWro/Genie)
            // const auto fluxUnitWeights = CrossSectionHelper::GetUnitWeightsMap(systParams.fluxDimensions);
            // const auto xsecUnitWeights = CrossSectionHelper::GetUnitWeightsMap(systParams.xsecDimensions);
            // const auto reintUnitWeights = CrossSectionHelper::GetUnitWeightsMap(systParams.reintDimensions);
            const auto sidebandUnitWeights = std::vector<float>(systParams.nSidebandFitUniverses, 1.f);
            CrossSectionHelper::SystFloatMap emptySystFloatMap = {};

            // std::cout<<"DEBUG K-1"<<std::endl;
            // std::cout<<"DEBUG TYT - P 9.1"<<std::endl;
            // Handle selected background events
            if (!isSignal && !isNuWro) // !isNuWro = isOverlay || isDataEXT || isDirt
            {
                for (auto &[selectionName, xsecs] : xsecMapNuWroScaled)
                {
                    // std::cout<<"DEBUG K-1"<<std::endl;
                    // Only use selected background events
                    if (!isSelectedMap.at(selectionName)) continue;
                    // std::cout<<"DEBUG TYT - P 9.2"<<std::endl;

                    for (auto &[name, xsec] : xsecs)
                    {
                        // std::cout<<"DEBUG TYT - P 9.3"<<std::endl;
                        // std::cout<<"DEBUG K0"<<std::endl;
                        const auto recoValue = getValue.at(name)(recoData);
                        const auto trueValue = getValue.at(name)(truthData);
                        // std::cout<<"DEBUG TYT - P 9.3.1 - name: "<<name<<std::endl;
                        const auto trueSidebandValue = getSidebandValue.at(name)(sidebandTruthData);
                        // std::cout<<"DEBUG TYT - P 9.3.1.1"<<std::endl;
                        const auto sidebandBinEdges = sidebandBinningMap.at(name);
                        // std::cout<<"DEBUG TYT - P 9.3.2"<<std::endl;

                        // Add CC0pi constraint
                        for (const auto &dataTypeName: {std::string("NuWro"), std::string("BNB"), std::string("Genie")})
                        {
                            // std::cout<<"DEBUG TYT - P 9.4"<<std::endl;
                            // std::cout<<"DEBUG K1"<<std::endl;
                            // std::cout<<"DEBUG K0"<<std::endl;
                            if (dataTypeName=="NuWro" && (!config.global.useNuWroAsData || !isOverlay)) continue; // use only overlay for nuwro
                            else if (dataTypeName=="Genie" && (!config.global.useGenieAsData || !isOverlay)) continue; // use only overlay for genie
                            else if (dataTypeName=="BNB" && (!config.global.useBNBAsData)) continue; // use overlay, dirt, extbnb for BNB
                            // std::cout<<"DEBUG K1.1"<<std::endl;
                            // std::cout<<"DEBUG K0.1 - dataTypeName: "<<dataTypeName<<std::endl;
                            // std::cout<<"DEBUG TYT - P 9.5"<<std::endl;
                            const auto cc0piNominalConstraintParam = cc0piNominalConstraintMap.at(dataTypeName).at(selectionName).at(name).first;
                            const auto cc0piNominalConstraintParamError = cc0piNominalConstraintMap.at(dataTypeName).at(selectionName).at(name).second;
                            // std::cout<<"DEBUG TYT - P 9.6"<<std::endl;

                            // std::cout<<"DEBUG K0.2"<<std::endl;
                            const auto nominalWeightFactor = isCC0PiSignal ? xsec.GetSidebandScaling(trueSidebandValue, cc0piNominalConstraintParam, sidebandBinEdges) : 1.f;

                            // std::cout<<"DEBUG K1"<<std::endl;
                            // Todo: arbitrary threshold in function. Fixed??? // Also change from vector<float> to SystFloatMap for consistency
                            const auto sidebandWeights = isCC0PiSignal ? xsec.GetSidebandParameterWeights(trueSidebandValue, cc0piNominalConstraintParam, cc0piNominalConstraintParamError, sidebandBinEdges) : sidebandUnitWeights;

                            if (dataTypeName=="Genie")
                            {
                                // std::cout<<"DEBUG TYT - P 9.7"<<std::endl;
                                // std::cout<<"DEBUG K1.0.1"<<" "<<selectionName<<" "<<name<<" - recoValue: "<<recoValue<<std::endl;
                                xsecMapGenieUnscaled.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, emptySystFloatMap, emptySystFloatMap, emptySystFloatMap, sidebandUnitWeights);
                                if (config.global.useCC0piConstraint) xsecMapGenieScaled.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, emptySystFloatMap, emptySystFloatMap, emptySystFloatMap, sidebandWeights, nominalWeightFactor);
                                // if (useForFluxReweighting) scalingDataScaledMap.at(dataTypeName).at(selectionName).at(name).pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights, nominalWeightFactor); //todo check if fluxWeights is correct here
                            }
                            else
                            {
                                // std::cout<<"DEBUG TYT - P 9.8"<<std::endl;
                                // std::cout<<"DEBUG K1.1"<<std::endl;
                                const auto cc0piUniverseConstraintVector = cc0piUniverseConstraintMap.at(dataTypeName).at(selectionName).at(name);
                                const auto scaledXSecWeights = isCC0PiSignal ? xsec.GetSidebandUniverseScaling(xsecWeights, trueSidebandValue, cc0piUniverseConstraintVector, sidebandBinEdges) : xsecWeights;


                                // for (auto &[paramName,weightVector] : xsecWeights)
                                // {
                                //     // std::cout<<"DEBUG xsecWeights extrema: "<<*std::min_element(weightVector.begin(), weightVector.end())<<" "<<*std::max_element(weightVector.begin(), weightVector.end())<<std::endl;
                                //     int i = 0;
                                //     // std::cout<<"DEBUG xsecWeights isCC0PiSignal: "<<isCC0PiSignal<<" - recoValue: "<<recoValue<<" - trueValue: "<<trueValue<<"- dataTypeName: "<<dataTypeName<<" - paramName: "<<paramName<<" - weight: "<<weight<<" - xsecWeightsScaleFactor: "<<xsecWeightsScaleFactor<<" - nominalWeightFactor: "<<nominalWeightFactor<<" | ";
                                //     for (auto &w : weightVector)
                                //     {
                                //         std::cout<<" "<<w;
                                //         i++;
                                //         if(i>10) break;
                                //     }
                                //     std::cout<<std::endl;
                                //     break;
                                // }

                                // for (auto &[paramName,weightVector] : scaledXSecWeights)
                                // {
                                //     // std::cout<<"DEBUG scaledXSecWeights extrema: "<<*std::min_element(weightVector.begin(), weightVector.end())<<" "<<*std::max_element(weightVector.begin(), weightVector.end())<<std::endl;
                                //     int i = 0;
                                //     // std::cout<<"DEBUG scaledXSecWeights - isCC0PiSignal: "<<isCC0PiSignal<<" | ";
                                //     for (auto &w : weightVector)
                                //     {
                                //         std::cout<<" "<<w;
                                //         i++;
                                //         if(i>10) break;
                                //     }
                                //     std::cout<<std::endl;
                                //     break;
                                // }
                                // std::cout<<"DEBUG K1.2 - trueSidebandValue: "<<trueSidebandValue<<" - dataTypeName: "<<dataTypeName<<" - name: "<<name<<" - isCC0PiSignal: "<<isCC0PiSignal<<std::endl;
                                const auto scaledFluxWeights = isCC0PiSignal ? xsec.GetSidebandUniverseScaling(fluxWeights, trueSidebandValue, cc0piUniverseConstraintVector, sidebandBinEdges) : fluxWeights;
                                // std::cout<<"DEBUG K1.2.1"<<std::endl;
                                const auto scaledReintWeights = isCC0PiSignal ? xsec.GetSidebandUniverseScaling(reintWeights, trueSidebandValue, cc0piUniverseConstraintVector, sidebandBinEdges) : reintWeights;
                                // std::cout<<"DEBUG K1.2.2"<<std::endl;

                                if (dataTypeName=="NuWro")
                                {
                                    // std::cout<<"DEBUG K1.3"<<std::endl;
                                    xsecMapNuWroUnscaled.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, sidebandUnitWeights);
                                    if (config.global.useCC0piConstraint) xsec.AddSelectedBackgroundEvent(recoValue, isDirt, weight, scaledFluxWeights, scaledXSecWeights, scaledReintWeights, sidebandWeights, nominalWeightFactor);
                                    // std::cout<<"DEBUG K1.4"<<std::endl;
                                    // if (useForFluxReweighting) scalingDataScaledMap.at(dataTypeName).at(selectionName).at(name).pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights, nominalWeightFactor); //todo check if fluxWeights is correct here
                                } else if (dataTypeName=="BNB")
                                {
                                    // std::cout<<"DEBUG K1.6"<<std::endl;
                                    xsecMapBNBUnscaled.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, fluxWeights, xsecWeights, reintWeights, sidebandUnitWeights);
                                    // std::cout<<"DEBUG K1.7"<<std::endl;
                                    if (config.global.useCC0piConstraint) xsecMapBNBScaled.at(selectionName).at(name).AddSelectedBackgroundEvent(recoValue, isDirt, weight, scaledFluxWeights, scaledXSecWeights, scaledReintWeights, sidebandWeights, nominalWeightFactor);
                                    // std::cout<<"DEBUG K1.8"<<std::endl;
                                    if (useForFluxReweighting) scalingDataScaledMap.at(dataTypeName).at(selectionName).at(name).pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, scaledFluxWeights, nominalWeightFactor);
                                    // std::cout<<"DEBUG K1.9"<<std::endl;
                                }
                                // std::cout<<"DEBUG TYT - P 9.9"<<std::endl;
                                // std::cout<<"DEBUG K1.10"<<std::endl;
                            }
                            // std::cout<<"DEBUG K2"<<std::endl;
                            // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
                            // if (useForFluxReweighting) scalingDataScaledMap.at(dataTypeName).at(selectionName).at(name).pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, scaledFluxWeights, nominalWeightFactor);
                        }
                    }
                }
            }
            // Handle other events
            else
            {
                // std::cout<<"DEBUG K3"<<std::endl;
                // Fill the flux-reweightor with all overlay events from a neutrinos of the desired flavour in the active volume
                if (useForFluxReweighting)
                {
                    for (auto &[dataTypeName, selectionNameMap] : scalingDataScaledMap)
                    {
                        for (auto &[selectionName, nameMap] : selectionNameMap)
                        {
                            for (auto &[name, scalingData]: nameMap)
                            {
                                scalingData.pFluxReweightor->AddEvent(pEvent->truth.nuEnergy(), weight, fluxWeights);
                            }
                        }
                    }
                }

                // std::cout<<"DEBUG K4"<<std::endl;
                // Handle signal events
                if (isSignal)
                {
                    for (auto &[selectionName, xsecs] : xsecMapNuWroScaled)
                    {
                        // Determine if we passed the relevant selection
                        const auto isSelected = isSelectedMap.at(selectionName);

                        for (auto &[name, xsec] : xsecs)
                        {
                            // std::cout<<"DEBUG K4.1"<<std::endl;
                            const auto recoValue = getValue.at(name)(recoData);
                            const auto trueValue = getValue.at(name)(truthData);
                            if (isNuWro)
                            {
                                xsecMapNuWroTrue.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, emptySystFloatMap, xsecWeights, emptySystFloatMap);
                            } else
                            {
                                // std::cout<<"DEBUG K4.2"<<std::endl;
                                if (isOverlay) // Genie MC for NuWro (no EXT or dirt events in NuWro)
                                {
                                    if (config.global.useNuWroAsData)
                                    {
                                        if (config.global.useCC0piConstraint) xsec.AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights);
                                        xsecMapNuWroUnscaled.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights);
                                    }
                                    if (config.global.useGenieAsData)
                                    {
                                        // std::cout<<"DEBUG K4.3"<<" "<<selectionName<<" "<<name<<" - recoValue: "<<recoValue<<" - trueValue: "<<trueValue<<" - isSelected: "<<isSelected<<std::endl;
                                        xsecMapGenieUnscaled.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, emptySystFloatMap, emptySystFloatMap, emptySystFloatMap);
                                        // std::cout<<"DEBUG K4.3.1"<<std::endl;
                                        if (config.global.useCC0piConstraint) xsecMapGenieScaled.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, emptySystFloatMap, emptySystFloatMap, emptySystFloatMap);
                                        // std::cout<<"DEBUG K4.4"<<std::endl;
                                    }
                                }
                                // std::cout<<"DEBUG K4.5"<<std::endl;
                                if (config.global.useBNBAsData) // Genie MC for Data (isOverlay || isDataEXT || isDirt)
                                {
                                    xsecMapBNBUnscaled.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights);
                                    if (config.global.useCC0piConstraint) xsecMapBNBScaled.at(selectionName).at(name).AddSignalEvent(recoValue, trueValue, isSelected, weight, fluxWeights, xsecWeights, reintWeights);
                                }
                                // std::cout<<"DEBUG K4.6"<<std::endl;
                            }
                        }
                    }
                }
                // std::cout<<"DEBUG K5"<<std::endl;
            }
        }
        // std::cout<<"DEBUG TYT - P 10"<<std::endl;
    }
    // Loop over all cross-section objects
    std::vector<std::tuple<std::string, std::string, NestedCrossSectionMap>> xsecMapVector;
    if (config.global.useNuWroAsData)
    {
        xsecMapVector.push_back(std::make_tuple(std::string("NuWro"), std::string("Unscaled"), xsecMapNuWroUnscaled));
        if (config.global.useCC0piConstraint) xsecMapVector.push_back(std::make_tuple(std::string("NuWro"), std::string("Scaled"), xsecMapNuWroScaled));
    }
    if (config.global.useGenieAsData)
    {
        xsecMapVector.push_back(std::make_tuple(std::string("Genie"), std::string("Unscaled"), xsecMapGenieUnscaled));
        if (config.global.useCC0piConstraint) xsecMapVector.push_back(std::make_tuple(std::string("Genie"), std::string("Scaled"), xsecMapGenieScaled));
    }
    if (config.global.useBNBAsData)
    {
        xsecMapVector.push_back(std::make_tuple(std::string("BNB"), std::string("Unscaled"), xsecMapBNBUnscaled));
        if (config.global.useCC0piConstraint) xsecMapVector.push_back(std::make_tuple(std::string("BNB"), std::string("Scaled"), xsecMapBNBScaled));
    }

    for (const auto&[dataTypeName, scaling, xsecMap]: xsecMapVector)
    {
        const auto postfix = dataTypeName + scaling;
        for (const auto &[selectionName, xsecs] : xsecMap)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                std::cout << "Processing cross-section: " << postfix << " - " << selectionName << " - " << name << std::endl;

                const auto scalingData = (scaling == "Unscaled") ? scalingDataUnscaled : scalingDataScaledMap.at(dataTypeName).at(selectionName).at(name);
                ExtractionHelper::SaveCrossSectionMatricies(xsec, scalingData, selectionName, name, postfix, false);
            }
        }
    }

    if (config.global.useNuWroAsData)
    {
        // Loop over all cross-section objects
        for (const auto &[selectionName, xsecs] : xsecMapNuWroTrue)
        {
            for (const auto &[name, xsec] : xsecs)
            {
                std::cout << "True: Processing cross-section: "<<selectionName<< " - " << name << std::endl;

                // Instead of using some extra object 'scalingDataNuWroTruth', scalingDataUnscaled is used.
                // It contains the right exposurePOT, nTargets and nominal flux for NuWroTruth.
                // The reweighted flux might be different, however, this is not needed here.
                const auto prediction = xsec.GetPredictedCrossSection(scalingDataUnscaled);

                std::cout << "Predicted NuWro cross-section (truth-space)" << std::endl;
                FormattingHelper::SaveMatrix(prediction, "xsec_" + selectionName + "_" + name + "_prediction_NuWroTruth.txt");

                std::cout << "NuWro smearing Matrix (reco-space rows, truth-space columns)" << std::endl;
                const auto smearingMatrix = xsec.GetSmearingMatrix();
                FormattingHelper::SaveMatrix(smearingMatrix, "xsec_" + selectionName + "_" + name + "_smearingMatrix_NuWroTruth.txt");

                std::cout << "NuWro smearing Matrix all selected (reco-space rows, truth-space columns)" << std::endl;
                const auto smearingMatrixAllSelected = xsec.GetSmearingMatrixAllSelected();
                FormattingHelper::SaveMatrix(smearingMatrixAllSelected, "xsec_" + selectionName + "_" + name + "_smearingMatrixAllSelected_NuWroTruth.txt");
            }
        }
    }

    std::cout<<"------------- All Done -------------"<<std::endl;
    // return;
}

} // namespace ubcc1pi_macros
