/**
 *  @file  ubcc1pi_standalone/Helpers/CrossSectionHelper.cxx
 *
 *  @brief The implementation of the cross section helper class
 */

#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/GeometryHelper.h"

#include <TFile.h>

#include <stdexcept>
#include <regex>

namespace ubcc1pi
{

CrossSectionHelper::FluxReweightor::FluxReweightor(const std::vector<float> &binEdges, const std::vector<float> &binValuesNominal, const SystDimensionsMap &fluxWeightsDimensions) :
    /// @cond Doxygen can't handle this initilizer list
    m_dimensions(fluxWeightsDimensions),
    m_pFluxNominal(CrossSectionHelper::GetTH1F(binEdges)),
    m_pSpectrumNominal(CrossSectionHelper::GetTH1F(binEdges)),
    m_spectrumVariations(CrossSectionHelper::GetSystTH1FMap(binEdges, fluxWeightsDimensions)),

    // Define the cache-function that gets the integrated flux variations
    // ATTN it might look a bit odd to some to define a function in the constructor! To be clear, here we are initalizing the member
    // variable named m_getIntegratedFluxVariationFunction. This is a SystCacheFunction that takes a lambda function as it's first paramter
    // and stores the logic for later. When the user calls GetIntegratedFluxVariation() we ask the SystCacheFunction to check if it has the
    // value we are looking for in it's cache. If the answer is yes, then we just return the cached value. If the answer is no, then we call
    // the logic (defined in the lambda function below), store the result in the cache for later and return it.
    m_getIntegratedFluxVariationFunction([&](const std::string &paramName, const unsigned int universeIndex)
    {
        // Setup a new histogram to hold the reweighted flux
        auto pReweightedFlux = CrossSectionHelper::GetTH1F(binEdges);

        // Get the event rate spectrum in this universe
        const auto pSpectrumVariation = m_spectrumVariations.at(paramName).at(universeIndex);

        // Loop over each bin
        const auto nBins = binEdges.size() - 1;
        for (unsigned iBin = 1; iBin <= nBins; ++iBin)
        {
            // Get bin value for the nominal spectrum and the specturm in this universe
            const auto nEventsNominal = m_pSpectrumNominal->GetBinContent(iBin);
            const auto nEventsVariation = pSpectrumVariation->GetBinContent(iBin);

            // Get the ratio of the variation to the nominal (if not possible use a unit weight)
            const auto weight = (nEventsNominal <= std::numeric_limits<float>::epsilon() ? 1.f : (nEventsVariation / nEventsNominal) );

            // Apply this weight to the nominal flux
            const auto nominalFlux = m_pFluxNominal->GetBinContent(iBin);
            const auto reweightedFlux = nominalFlux * weight;

            // Store the result
            pReweightedFlux->SetBinContent(iBin, reweightedFlux);
        }

        // Return the integrated flux
        return this->GetIntegratedFlux(pReweightedFlux);

    }, m_dimensions)
    /// @endcond
{
    // There should be one more bin edge than the number of bins
    if (binEdges.size() != binValuesNominal.size() + 1)
        throw std::invalid_argument("FluxReweightor::FluxReweightor - Number of bin edges doesn't match number of bins!");

    // Fill the nominal flux histogram
    for (unsigned int iBin = 1; iBin <= binValuesNominal.size(); ++iBin)
    {
        // ATTN here we shift the index by 1 as root enumerates bins from 1 and the vector enumerates from 0
        m_pFluxNominal->SetBinContent(iBin, binValuesNominal.at(iBin - 1));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::FluxReweightor::AddEvent(const float trueNuEnergy, const float unscaledNominalWeight, const SystFloatMap &fluxWeights, const float nominalWeightFactor)
{
    // Make sure we have a valid input map
    CrossSectionHelper::ValidateSystMap(fluxWeights, m_dimensions);

    // Here we clear any cached values for the integrated flux variations because we are modifying the event rate spectra
    m_getIntegratedFluxVariationFunction.ClearCache();

    // Fill the nominal spectrum
    m_pSpectrumNominal->Fill(trueNuEnergy, unscaledNominalWeight * nominalWeightFactor);

    // Fill the universes
    for (const auto &[paramName, weights] : fluxWeights)
    {
        auto &spectrumUniverses = m_spectrumVariations.at(paramName);

        for (unsigned int iUni = 0; iUni < weights.size(); ++iUni)
        {
            const auto universeWeight = unscaledNominalWeight * weights.at(iUni);
            spectrumUniverses.at(iUni)->Fill(trueNuEnergy, universeWeight);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelper::FluxReweightor::GetNominalFlux() const
{
    return m_pFluxNominal;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::FluxReweightor::GetIntegratedFlux(const std::shared_ptr<TH1F> &pFlux) const
{
    float total = 0.f;
    for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pFlux->GetNbinsX()); ++iBin)
    {
        total += pFlux->GetBinContent(iBin);
    }

    return total;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::FluxReweightor::GetIntegratedNominalFlux() const
{
    return this->GetIntegratedFlux(m_pFluxNominal);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float CrossSectionHelper::FluxReweightor::GetIntegratedFluxVariation(const std::string &paramName, const unsigned int universeIndex)
{
    // Call the cache-function to either calculate the return value or retrieve it from the cache
    return m_getIntegratedFluxVariationFunction(paramName, universeIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::CrossSection::CrossSection(const SystParams &systParams, const std::vector<float> &binEdges, const bool hasUnderflow, const bool hasOverflow, const bool scaleByBinWidth, const std::vector<float> &trueBinEdges) :
    m_systParams(systParams),
    m_binEdges(binEdges),
    m_metadata(binEdges, hasUnderflow, hasOverflow, scaleByBinWidth),
    m_scaleByBinWidth(scaleByBinWidth),
    m_pSignal_true_nom(CrossSectionHelper::GetTH1F(trueBinEdges.empty() ? binEdges : trueBinEdges)),
    m_pSignal_selected_recoTrue_nom(CrossSectionHelper::GetTH2F(binEdges, trueBinEdges.empty() ? binEdges : trueBinEdges)),
    m_pBackground_selected_reco_nom(CrossSectionHelper::GetTH1F(binEdges)),
    m_pBNBData_selected_reco(CrossSectionHelper::GetTH1F(binEdges))
    {
    // Make sure the bin edges are valid
    const int nAnalysisBins = (binEdges.size() - 1) - (hasUnderflow ? 1 : 0) - (hasOverflow ? 1 : 0);
    if (nAnalysisBins <= 0)
        throw std::invalid_argument("CrossSection::CrossSection - Insufficient bin edges provided!");


    std::cout<<"binEdges:";
    for(const auto &binEdge: binEdges)
    {
        std::cout<<" "<<binEdge;
    }
    std::cout<<std::endl;
    std::cout<<"trueBinEdges:";
    for(const auto &trueBinEdge: trueBinEdges)
    {
        std::cout<<" "<<trueBinEdge;
    }
    std::cout<<std::endl;

    // Setup a dimensions map for the miscellaneous ("misc") multisim parameters.
    // These are added in ad-hoc (i.e. they don't come from the weights stored in the truth information)
    const SystDimensionsMap miscDimensions = {

        // Here we apply a weight to each event using the Poisson bootstrap method to estimate the MC stats uncertainty
        { "bootstrap", systParams.nBootstrapUniverses },

        // Apply the normal weights
        { "sidebandWeights", systParams.nBootstrapUniverses },

        // Here we apply a 2-universe multisim in which we weight the dirt up and down by 100%
        { "dirt", 2 }
    };

    const auto binEdges2 = trueBinEdges.empty() ? binEdges : trueBinEdges;
    // Setup the maps to hold the histograms for the multisim parameters
    m_signal_true_multisims.emplace("flux", CrossSectionHelper::GetSystTH1FMap(binEdges2, systParams.fluxDimensions));
    m_signal_true_multisims.emplace("xsec", CrossSectionHelper::GetSystTH1FMap(binEdges2, systParams.xsecDimensions));
    m_signal_true_multisims.emplace("reint", CrossSectionHelper::GetSystTH1FMap(binEdges2, systParams.reintDimensions));
    m_signal_true_multisims.emplace("misc", CrossSectionHelper::GetSystTH1FMap(binEdges2, miscDimensions));

    m_signal_selected_recoTrue_multisims.emplace("flux", CrossSectionHelper::GetSystTH2FMap(binEdges, systParams.fluxDimensions, binEdges2));
    m_signal_selected_recoTrue_multisims.emplace("xsec", CrossSectionHelper::GetSystTH2FMap(binEdges, systParams.xsecDimensions, binEdges2));
    m_signal_selected_recoTrue_multisims.emplace("reint", CrossSectionHelper::GetSystTH2FMap(binEdges, systParams.reintDimensions, binEdges2));
    m_signal_selected_recoTrue_multisims.emplace("misc", CrossSectionHelper::GetSystTH2FMap(binEdges, miscDimensions, binEdges2));

    m_background_selected_reco_multisims.emplace("flux", CrossSectionHelper::GetSystTH1FMap(binEdges, systParams.fluxDimensions));
    m_background_selected_reco_multisims.emplace("xsec", CrossSectionHelper::GetSystTH1FMap(binEdges, systParams.xsecDimensions));
    m_background_selected_reco_multisims.emplace("reint", CrossSectionHelper::GetSystTH1FMap(binEdges, systParams.reintDimensions));
    m_background_selected_reco_multisims.emplace("misc", CrossSectionHelper::GetSystTH1FMap(binEdges, miscDimensions));

    // Setup the maps to hold the histograms for the unisim parameters
    m_signal_true_unisims.emplace("detector", CrossSectionHelper::GetSystUnisimTH1FMap(binEdges2, systParams.detVarDimensions));
    m_signal_selected_recoTrue_unisims.emplace("detector", CrossSectionHelper::GetSystUnisimTH2FMap(binEdges, systParams.detVarDimensions, binEdges2));
    m_background_selected_reco_unisims.emplace("detector", CrossSectionHelper::GetSystUnisimTH1FMap(binEdges, systParams.detVarDimensions));
    // ATTN we could in theory include other unisim parameters here using a different name. At the moment we only use "detector".
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSignalEvent(const float recoValue, const float trueValue, const bool isSelected, 
    const float nominalWeight, const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights, const SystFloatMap &reintWeights)
{
    const auto bootstrapWeights = CrossSectionHelper::GenerateBootstrapWeights(m_systParams.nBootstrapUniverses);

    const auto sidebandWeights = std::vector<float>(m_systParams.nBootstrapUniverses, 1.0);

    // Get the miscellaneous weights
    const SystFloatMap miscWeights = {

        // Generate the bootstrap weights
        { "bootstrap", bootstrapWeights },

        // Apply the normal weights
        { "sidebandWeights", sidebandWeights },

        // A signal event is never dirt - so just use a unit weight
        { "dirt", {1.f, 1.f} }
    };

    m_pSignal_true_nom->Fill(trueValue, nominalWeight);
    CrossSectionHelper::FillSystTH1FMap(trueValue, nominalWeight, fluxWeights, m_signal_true_multisims.at("flux"));
    CrossSectionHelper::FillSystTH1FMap(trueValue, nominalWeight, xsecWeights, m_signal_true_multisims.at("xsec"));
    CrossSectionHelper::FillSystTH1FMap(trueValue, nominalWeight, reintWeights, m_signal_true_multisims.at("reint"));
    CrossSectionHelper::FillSystTH1FMap(trueValue, nominalWeight, miscWeights, m_signal_true_multisims.at("misc"));

    if (!isSelected)
        return;

    // ATTN the reco value is only used if the event is selected
    m_pSignal_selected_recoTrue_nom->Fill(recoValue, trueValue, nominalWeight);
    CrossSectionHelper::FillSystTH2FMap(recoValue, trueValue, nominalWeight, fluxWeights, m_signal_selected_recoTrue_multisims.at("flux"));
    CrossSectionHelper::FillSystTH2FMap(recoValue, trueValue, nominalWeight, xsecWeights, m_signal_selected_recoTrue_multisims.at("xsec"));
    CrossSectionHelper::FillSystTH2FMap(recoValue, trueValue, nominalWeight, reintWeights, m_signal_selected_recoTrue_multisims.at("reint"));
    CrossSectionHelper::FillSystTH2FMap(recoValue, trueValue, nominalWeight, miscWeights, m_signal_selected_recoTrue_multisims.at("misc"));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSignalEventDetVar(const float recoValue, const float trueValue, const bool isSelected, const float nominalWeight, const std::string &paramName)
{
    if (!static_cast<bool>(m_signal_true_unisims.at("detector").count(paramName)))
    {
        std::cout<<"CrossSection::AddSignalEventDetVar - Unknown parameter or CV sample name \"" << paramName << "\""<<std::endl;
        throw std::invalid_argument("CrossSection::AddSignalEventDetVar - Unknown parameter or CV sample name \"" + paramName + "\"");
    }
    m_signal_true_unisims.at("detector").at(paramName)->Fill(trueValue, nominalWeight);

    if (!isSelected)
        return;

    m_signal_selected_recoTrue_unisims.at("detector").at(paramName)->Fill(recoValue, trueValue, nominalWeight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBackgroundEvent(const float recoValue, const bool isDirt, const float unscaledNominalWeight,
    const SystFloatMap &fluxWeights, const SystFloatMap &xsecWeights, const SystFloatMap &reintWeights, std::vector<float> sidebandWeights, const float nominalWeightFactor)
{
    // The weight to apply as a result of dirt
    //   - For non-dirt events, just use zero (i.e. no change from the nominal simulation)
    //   - For dirt event we here use +-100%
    const auto dirtWeightDelta = isDirt ? 1.f : 0.f;
    auto bootstrapWeights = CrossSectionHelper::GenerateBootstrapWeights(m_systParams.nBootstrapUniverses);
    

    if(sidebandWeights.empty())
    {    
        sidebandWeights = std::vector<float>(m_systParams.nBootstrapUniverses, 1.0);
    } else if(sidebandWeights.size() != m_systParams.nBootstrapUniverses)
    {
        std::cout<<"CrossSectionHelper::CrossSection::AddSelectedBackgroundEvent - Number of sideband universes does not match vector size."<<std::endl;
        throw std::logic_error("CrossSectionHelper::CrossSection::AddSelectedBackgroundEvent - Number of sideband universes does not match vector size.");
    }

    // std::cout<<"AddSelectedBackgroundEvent - sidebandWeights:";
    // for (const auto &s : sidebandWeights)
    //     std::cout<<" "<<s;
    // std::cout<<std::endl;

    // std::cout<<"AddSelectedBackgroundEvent - xsec_All_UBGenie:";
    // for (const auto &s : xsecWeights.at("All_UBGenie"))
    //     std::cout<<" "<<s;
    // std::cout<<std::endl;


    // Get the miscellaneous weights
    for (auto &b: bootstrapWeights) b *= nominalWeightFactor;
    const SystFloatMap miscWeights = {

        // Apply the bootstrap weights
        { "bootstrap", bootstrapWeights },

        // Apply the normal weights
        { "sidebandWeights", sidebandWeights }, // they already have weightFactors integrated

        // Apply the dirt weight
        { "dirt", {(1.f - dirtWeightDelta)*nominalWeightFactor, (1.f + dirtWeightDelta)*nominalWeightFactor} }
    };


    m_pBackground_selected_reco_nom->Fill(recoValue, unscaledNominalWeight*nominalWeightFactor);
    CrossSectionHelper::FillSystTH1FMap(recoValue, unscaledNominalWeight, fluxWeights, m_background_selected_reco_multisims.at("flux")); // universe weights are included in fluxWeights
    CrossSectionHelper::FillSystTH1FMap(recoValue, unscaledNominalWeight, xsecWeights, m_background_selected_reco_multisims.at("xsec")); // universe weights are included in xsecWeights
    CrossSectionHelper::FillSystTH1FMap(recoValue, unscaledNominalWeight, reintWeights, m_background_selected_reco_multisims.at("reint")); // universe weights are included in reintWeights
    CrossSectionHelper::FillSystTH1FMap(recoValue, unscaledNominalWeight, miscWeights, m_background_selected_reco_multisims.at("misc")); // nominalWeightFactor/universe weights are included in miscWeights
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBackgroundEventDetVar(const float recoValue, const float nominalWeight, const std::string &paramName)
{
    if (!static_cast<bool>(m_background_selected_reco_unisims.at("detector").count(paramName)))
        throw std::invalid_argument("CrossSection::AddSelectedBackgroundEventDetVar - Unknown parameter or CV sample name \"" + paramName + "\"");

    m_background_selected_reco_unisims.at("detector").at(paramName)->Fill(recoValue, nominalWeight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddSelectedBNBDataEvent(const float recoValue)
{
    m_pBNBData_selected_reco->Fill(recoValue);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::CrossSection::AddWeightedSelectedDataEvent(const float recoValue, const float weight = 1.f)
{
    m_pBNBData_selected_reco->Fill(recoValue, weight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::CrossSection::HasOverflow() const
{
    return m_metadata.HasOverflow();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool CrossSectionHelper::CrossSection::HasUnderflow() const
{
    return m_metadata.HasUnderflow();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBXSecMeta CrossSectionHelper::CrossSection::GetMetadata() const
{
    return m_metadata;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::CrossSection::GetBinEdges() const
{
    return m_binEdges;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBinWidths() const
{
    return m_metadata.GetBinWidths();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetCrossSection(const ubsmear::UBMatrix &selected, const ubsmear::UBMatrix &backgrounds, const float integratedFlux, const float exposurePOT, const float nTargets) const
{
    // Check the input event rate has the right dimensions
    if (selected.GetColumns() != 1)
    {
        std::cout << "CrossSection::GetCrossSection  - Input event rate isn't a column vector" << std::endl;
        throw std::invalid_argument("CrossSection::GetCrossSection - Input event rate isn't a column vector");
    }

    if (selected.GetRows() != m_metadata.GetNBins())
    {
        std::cout << "CrossSection::GetCrossSection - Input event rate has the wrong number of bins" << std::endl;
        throw std::invalid_argument("CrossSection::GetCrossSection - Input event rate has the wrong number of bins");
    }

    // Check the input backgrounds have the right dimensions
    if (backgrounds.GetColumns() != 1)
    {
        std::cout << "CrossSection::GetCrossSection - Input background rate isn't a column vector" << std::endl;
        throw std::invalid_argument("CrossSection::GetCrossSection - Input background rate isn't a column vector");
    }

    if (backgrounds.GetRows() != m_metadata.GetNBins())
    {
        std::cout << "CrossSection::GetCrossSection - Input background rate has the wrong number of bins" << std::endl;
        throw std::invalid_argument("CrossSection::GetCrossSection - Input background rate has the wrong number of bins");
    }

    // Get the normalisation factor
    const auto norm = integratedFlux * exposurePOT * nTargets;
    if (norm <= std::numeric_limits<float>::epsilon())
    {
        std::cout << "CrossSection::GetCrossSection - Product of flux, exposure and targets non-positive" << std::endl;
        throw std::invalid_argument("CrossSection::GetCrossSection - Product of flux, exposure and targets non-positive");
    }

    // Get the column vector of bin widths
    const auto binWidths = this->GetBinWidths();

    // Now get the flux-integrated, forward-folded cross-section by
    //   - Subtracting the backgrounds
    //   - Scaling by the normalisation factor
    //   - Scaling (element wise) by the bin-widths
    return ubsmear::ElementWiseOperation(
        (selected - backgrounds),
        binWidths * norm,
        [](const float numerator, const float denominator) {

            if (denominator <= std::numeric_limits<float>::epsilon())
            {
                std::cout << "CrossSection::GetCrossSection - Warning: denominator is zero" << std::endl;
                throw std::logic_error("CrossSection::GetCrossSection - Found a bin in which the denominator (flux * POT * nTargets * binWidth) is <= 0");
            }

            return numerator / denominator;
        }
    );
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrix(const std::shared_ptr<TH1F> &pSignal_true, const std::shared_ptr<TH2F> &pSignal_selected_recoTrue) const
{
    // Get the input histograms as a matrix
    const auto allTrue = CrossSectionHelper::GetMatrixFromHist(pSignal_true);
    const auto selectedRecoTrue = CrossSectionHelper::GetMatrixFromHist(pSignal_selected_recoTrue);

    // Check theh matrices have the expected sizes
    if (!ubsmear::UBMatrixHelper::IsSquare(selectedRecoTrue))
    {
        std::cout<<"CrossSection::GetSmearingMatrix - input reco-true matrix isn't square"<<std::endl;
        throw std::invalid_argument("CrossSection::GetSmearingMatrix - input reco-true matrix isn't square");
    }

    const auto nBins = allTrue.GetRows();
    if (selectedRecoTrue.GetRows() != nBins)
    {
        std::cout<<"CrossSection::GetSmearingMatrix - the input objects differ in number of bins"<<std::endl;
        throw std::invalid_argument("CrossSection::GetSmearingMatrix - the input objects differ in number of bins");
    }

    // Get the elements of the smearing matrix
    std::vector<float> elements;
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
        {
            // Get the total number of signal events in this true bin
            // const auto denominator = allTrue.At(iTrue, 0);
            auto denominator = allTrue.At(iTrue, 0);

            // ATTN. If we find a bin in which there are no signal events, then there's no sensible value for the smearing matrix element
            // Here we return null to indicate that the smearing matrix can't be found
            if (denominator <= std::numeric_limits<float>::epsilon())
            {
                if(iReco==0 && iTrue==0)
                { 
                    std::cout<<"DEBUG CrossSectionHelper GetSmearingMatrix (true bins) - denominator: " << denominator << " ";
                    for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
                    {
                        // Get the total number of signal events in this true bin
                        std::cout<<allTrue.At(iTrue, 0)<<" ";
                    }
                    std::cout<<std::endl;
                }
                denominator = 1.0;
                // return std::shared_ptr<ubsmear::UBMatrix>(nullptr);
            }

            elements.push_back(selectedRecoTrue.At(iReco, iTrue) / denominator);
        }
    }

    return std::make_shared<ubsmear::UBMatrix>(elements, nBins, nBins);
}


// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixAllSelected(const std::shared_ptr<TH2F> &pSignal_selected_recoTrue) const
{
    // Get the input histograms as a matrix
    // const auto allTrue = CrossSectionHelper::GetMatrixFromHist(pSignal_true);
    const auto allTrue = GetSignalSelectedTrue(pSignal_selected_recoTrue);
    const auto selectedRecoTrue = CrossSectionHelper::GetMatrixFromHist(pSignal_selected_recoTrue);

    // Check theh matrices have the expected sizes
    if (!ubsmear::UBMatrixHelper::IsSquare(selectedRecoTrue))
    {
        std::cout<<"CrossSection::GetSmearingMatrix - input reco-true matrix isn't square"<<std::endl;
        throw std::invalid_argument("CrossSection::GetSmearingMatrix - input reco-true matrix isn't square");
    }
    const auto nBins = allTrue.GetRows();
    if (selectedRecoTrue.GetRows() != nBins)
    {
        std::cout<<"CrossSection::GetSmearingMatrix - the input objects differ in number of bins"<<std::endl;
        throw std::invalid_argument("CrossSection::GetSmearingMatrix - the input objects differ in number of bins");
    }
    // Get the elements of the smearing matrix
    std::vector<float> elements;
    for (unsigned int iReco = 0; iReco < nBins; ++iReco)
    {
        for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
        {
            // Get the total number of signal events in this true bin
            // const auto denominator = allTrue.At(iTrue, 0);
            auto denominator = allTrue.At(iTrue, 0);

            // ATTN. If we find a bin in which there are no signal events, then there's no sensible value for the smearing matrix element
            // Here we return null to indicate that the smearing matrix can't be found
            if (denominator <= std::numeric_limits<float>::epsilon())
            {
                if(iReco==0 && iTrue==0)
                { 
                    std::cout<<"DEBUG CrossSectionHelper GetSmearingMatrixAllSelected (true bins) - denominator: " << denominator << " ";
                    for (unsigned int iTrue = 0; iTrue < nBins; ++iTrue)
                    {
                        // Get the total number of signal events in this true bin
                        std::cout<<allTrue.At(iTrue, 0)<<" ";
                    }
                    std::cout<<std::endl;
                }
                denominator = 1.0;
                // return std::shared_ptr<ubsmear::UBMatrix>(nullptr);
            }

            elements.push_back(selectedRecoTrue.At(iReco, iTrue) / denominator);
        }
    }

    return std::make_shared<ubsmear::UBMatrix>(elements, nBins, nBins);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSignalSelectedTrue(const std::shared_ptr<TH2F> &pSignal_selected_recoTrue) const
{
    // Get the input histogram as a matrix
    const auto recoTrueMatrix = CrossSectionHelper::GetMatrixFromHist(pSignal_selected_recoTrue);

    // Integrate over the reco bins
    std::vector<float> elements;
    for (unsigned int iTrue = 0; iTrue < recoTrueMatrix.GetColumns(); ++iTrue)
    {
        float sum = 0.f;
        for (unsigned int iReco = 0; iReco < recoTrueMatrix.GetRows(); ++iReco)
        {
            sum += recoTrueMatrix.At(iReco, iTrue);
        }

        elements.push_back(sum);
    }

    // Return the result as a column vector
    return ubsmear::UBMatrix(elements, recoTrueMatrix.GetColumns(), 1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSelectedBNBDataEvents() const
{
    return CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSelectedBackgroundEvents() const
{
    return CrossSectionHelper::GetMatrixFromHist(m_pBackground_selected_reco_nom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSelectedSignalEvents() const
{
    return this->GetSignalSelectedTrue(m_pSignal_selected_recoTrue_nom);
}

// // -----------------------------------------------------------------------------------------------------------------------------------------

// ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSignalSelectedTrueMatrix() const
// {
//     return CrossSectionHelper::GetMatrixFromHist(m_pSignal_selected_recoTrue_nom);
// }

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSignalEvents() const
{
    return CrossSectionHelper::GetMatrixFromHist(m_pSignal_true_nom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSection(const ScalingData &scalingData) const
{
    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get the number of predicted backgrounds in the nominal universe
    const auto backgrounds = CrossSectionHelper::GetMatrixFromHist(m_pBackground_selected_reco_nom);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Get the cross-section
    return this->GetCrossSection(selected, backgrounds, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetNuWroTrueCrossSection(const ScalingData &scalingData) const
{
    // Get the number of selected events in BNB data
    const auto signal = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(signal.GetRows(), 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Get the cross-section
    return this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}


// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetPredictedCrossSection(const ScalingData &scalingData) const
{
    // Get the number of signal events in the nominal simulation in true bins
    const auto signal = CrossSectionHelper::GetMatrixFromHist(m_pSignal_true_nom);

    // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(signal.GetRows(), 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Get the cross-section
    return this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetPredictedCrossSectionInUniverse(const std::string &group, const std::string &paramName, const unsigned int universeIndex, const ScalingData &scalingData) const
{
    // Get the number of signal events in the nominal simulation in true bins
    const auto signal = CrossSectionHelper::GetMatrixFromHist(m_signal_true_multisims.at(group).at(paramName).at(universeIndex));

    // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(signal.GetRows(), 1);

    // Get the integrated flux in the supplied universe (if it's not a flux parameter, then use the nominal universe)
    const auto integratedFlux = (
        group == "flux"
            ? scalingData.pFluxReweightor->GetIntegratedFluxVariation(paramName, universeIndex)
            : scalingData.pFluxReweightor->GetIntegratedNominalFlux()
    );

    // Get the cross-section
    return this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}


// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixNominal() const
{
    // Get the smearing matrix in the nominal universe
    return this->GetSmearingMatrix(m_pSignal_true_nom, m_pSignal_selected_recoTrue_nom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSmearingMatrix() const
{
    // Get the smearing matrix in the nominal universe
    const auto pSmearingMatrix = this->GetSmearingMatrixNominal();

    // Check we were able to find the smearing matrix
    if (!pSmearingMatrix)
        throw std::logic_error("CrossSection::GetSmearingMatrix - There's a bin with no signal events! Can't get the smearing matrix");

    return *pSmearingMatrix;
}

// // // -----------------------------------------------------------------------------------------------------------------------------------------

// ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSmearingMatrixInUniverseAllSelected(const std::string &group, const std::string &paramName, const unsigned int universeIndex) const
// {
//     // Get the smearing matrix in the nominal universe
//     const auto pSmearingMatrix = this->GetSmearingMatrixInUniverse(const std::string &group, const std::string &paramName, const unsigned int universeIndex);

//     // Check we were able to find the smearing matrix
//     if (!pSmearingMatrix)
//         throw std::logic_error("CrossSection::GetSmearingMatrixInUniverseWithCheck - There's a bin with no signal events! Can't get the smearing matrix");

//     return *pSmearingMatrix;
// }

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixInUniverseAllSelected(const std::string &group, const std::string &paramName, const unsigned int universeIndex) const
{
    return this->GetSmearingMatrixAllSelected(m_signal_selected_recoTrue_multisims.at(group).at(paramName).at(universeIndex));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetSmearingMatrixAllSelected() const
{
    // Get the smearing matrix in the nominal universe
    const auto pSmearingMatrix = this->GetSmearingMatrixAllSelected(m_pSignal_selected_recoTrue_nom);

    // Check we were able to find the smearing matrix
    if (!pSmearingMatrix)
    {
        std::cout<<"CrossSection::GetSmearingMatrix - There's a bin with no signal events! Can't get the smearing matrix"<<std::endl;
        throw std::logic_error("CrossSection::GetSmearingMatrix - There's a bin with no signal events! Can't get the smearing matrix");
    }
    return *pSmearingMatrix;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSectionStatUncertainty(const ScalingData &scalingData) const
{
    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get a vector containing the count uncertainties (Poisson) on each bin
    std::vector<float> elements;
    const auto nBins = selected.GetRows();
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
    {
        elements.push_back(AnalysisHelper::GetCountUncertainty(selected.At(iBin, 0)));
    }

    const ubsmear::UBMatrix selectedUncertainty(elements, nBins, 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Scale the uncertainty on the number of selected events by same factor as the cross-section.
    // ATTN to apply the scaling we reuse the GetCrossSection function but intead pass it the count uncertainties as "selected" a zero
    // vector as the "backgrounds"
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
    return this->GetCrossSection(selectedUncertainty, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovarianceMap CrossSectionHelper::CrossSection::GetBNBDataCrossSectionSystUncertainties(const ScalingData &scalingData) const
{
    // Setup the output map
    SystBiasCovarianceMap outputMap;

    // Add the multisim parameters
    // ATTN the use of m_signal_true_multisims is arbitrary, we could use any multisims map
    for (const auto &[group, map] : m_signal_true_multisims)
    {
        for (const auto &entry : map)
        {
            const auto &paramName = entry.first;
            outputMap[group].emplace(paramName, this->GetBNBDataCrossSectionDistributionParams(group, paramName, scalingData));
        }
    }

    // Add the unisim parameters
    // ATTN again the use of m_signal_true_unisims is arbitrary, we could use any unisims map
    for (const auto &[group, map] : m_signal_true_unisims)
    {
        // ATTN for now we only have detector parameters that are unisims, if others are added later then we would need to add a clause here
        // to get the relevent dimensions.
        if (group != "detector")
            throw std::logic_error("CrossSection::GetBNBDataCrossSectionSystUncertainties - Don't know dimensions of group: \"" + group + "\"");

        const auto dimensions = m_systParams.detVarDimensions;

        // Extract the central-value names from the dimensions object
        std::vector<std::string> cvNames;
        for (const auto &[paramName, cvName] : dimensions)
            cvNames.push_back(cvName);

        for (const auto &entry : map)
        {
            const auto &paramName = entry.first;

            // ATTN The paramName could either be a detector variation sample of the name of a CV sample.
            // Here we check if the paramName is for a CV sample and if so, we can skip it.
            const auto isCVName = (std::find(cvNames.begin(), cvNames.end(), paramName) != cvNames.end());
            if (isCVName)
                continue;

            // Get the central-value name that corresponds to this parameter
            const auto &cvName = dimensions.at(paramName);

            outputMap[group].emplace(paramName, this->GetBNBDataCrossSectionDistributionParamsUnisim(group, paramName, cvName, scalingData));
        }
    }

    // Get the special POT normalisation uncertainty
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSection(scalingData) );
    outputMap["misc"].emplace("POT", this->GetDistributionParamsNormalisation(pXSecNom, m_systParams.potFracUncertainty));

    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovarianceMap CrossSectionHelper::CrossSection::GetSmearingMatrixSystUncertainties() const
{
    // Setup the output map
    SystBiasCovarianceMap outputMap;

    std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 0"<<std::endl;

    // Add the multisim parameters
    // ATTN the use of m_signal_true_multisims is arbitrary, we could use any multisims map
    for (const auto &[group, map] : m_signal_true_multisims)
    {
        std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 1"<<std::endl;
        for (const auto &entry : map)
        {
            std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 2"<<std::endl;
            const auto &paramName = entry.first;
            outputMap[group].emplace(paramName, this->GetSmearingMatrixDistributionParams(group, paramName));
            std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 3 - paramName: "<<paramName<<std::endl;
        }
    }
    std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 4"<<std::endl;
    // Add the unisim parameters
    // ATTN again the use of m_signal_true_unisims is arbitrary, we could use any unisims map
    for (const auto &[group, map] : m_signal_true_unisims)
    {
        std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 5"<<std::endl;
        // ATTN for now we only have detector parameters that are unisims, if others are added later then we would need to add a clause here
        // to get the relevent dimensions.
        if (group != "detector")
        {
            std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 5.1"<<std::endl;
            throw std::logic_error("CrossSection::GetSmearingMatrixSystUncertainties - Don't know dimensions of group: \"" + group + "\"");
        }

        const auto dimensions = m_systParams.detVarDimensions;

        // Extract the central-value names from the dimensions object
        std::vector<std::string> cvNames;
        for (const auto &[paramName, cvName] : dimensions)
            cvNames.push_back(cvName);
        
        std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 6"<<std::endl;

        for (const auto &entry : map)
        {
            const auto &paramName = entry.first;

            // ATTN The paramName could either be a detector variation sample or the name of a CV sample.
            // Here we check if the paramName is for a CV sample and if so, we can skip it.
            const auto isCVName = (std::find(cvNames.begin(), cvNames.end(), paramName) != cvNames.end());
            if (isCVName)
                continue;

            // Get the central-value name that corresponds to this parameter
            const auto &cvName = dimensions.at(paramName);

            outputMap[group].emplace(paramName, this->GetSmearingMatrixDistributionParamsUnisim(group, paramName, cvName));
        }
        std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 7"<<std::endl;
    }

    std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 8"<<std::endl;
    // Get the special POT normalisation uncertainty
    // ATTN no uncertainty is assigned to the smearing matrix due to POT normalisation (use a factor of 0.f)
    const auto pSmearingNom = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixNominal());
    outputMap["misc"].emplace("POT", this->GetDistributionParamsNormalisation(pSmearingNom, 0.f));

    std::cout<<"### DEBUG - GetSmearingMatrixSystUncertainties - Point 9"<<std::endl;
    return outputMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

// CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetSidebandWeights(const ScalingData &scalingData) const
// {
//     // Get the number of universes for this systematic parameter
//     // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
//     const auto nUniverses = m_systParams.nBootstrapUniverses;

//     // Get the nominal cross-section
//     const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>(GetPredictedCrossSection(scalingData));

//     // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
//     const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(pXSecNom->GetRows(), 1);

//     // Get the integrated flux in the nominal universe
//     const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

//     // The bootstrap universes
//     const auto universes = m_signal_true_multisims.at("misc").at("bootstrap");

//     // Get the distribution parameters
//     return this->GetDistributionParams([&](const unsigned int &iUni)
//     {
//         // This is the function we call to get the value in each universe

//         // Get the number of signal events in the nominal simulation in true bins
//         const auto signal = CrossSectionHelper::GetMatrixFromHist(universes.at(iUni));

//         // // Get the cross-section
//         // return std::make_shared<ubsmear::UBMatrix> ( this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets) );

//     }, nUniverses, pXSecNom);
// }

// -----------------------------------------------------------------------------------------------------------------------------------------

// CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetSidebandWeights(std::vector<float> &x, std::vector<float> &y, std::vector<float> &errorY, std::vector<float> &S;)
// {
//     const auto nUniverses = m_systParams.nBootstrapUniverses;
//     if (nUniverses == 0)
//         throw std::logic_error("CrossSection::GetDisributionParams - No universes supplied");


//     const auto selectedEventsData = this->GetSelectedBNBDataEvents();
//     const auto smearingMatrixNominal = this->GetSmearingMatrix();
//     const auto selectedEventsBackgroundNominalReco = this->GetSelectedBackgroundEvents();
//     const auto selectedEventsSignalNominalTruth = this->GetSelectedSignalEvents();
//     const auto signalData = selectedEventsData - selectedEventsBackgroundNominalReco;

//     const auto selectedSignalUniverses = m_signal_selected_recoTrue_multisims.at("misc").at("bootstrap");
//     const auto selectedBackgroundUniverses = m_background_selected_reco_multisims.at("misc").at("bootstrap");
//     for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
//     {
//         const auto selectedSignalTruth = this->GetSignalSelectedTrue(selectedSignalUniverses.at(iUni));
//         const auto selectedBackgoundReco = CrossSectionHelper::GetMatrixFromHist(selectedBackgroundUniverses.at(iUni));
//     }



//     // Insist the nominal input is a column vector
//     if (pNominal->GetColumns() != 1)
//         throw std::logic_error("CrossSection::GetDisributionParams - Input nominal is not a column vector");

//     // Get the number of bins
//     const auto nBins = pNominal->GetRows();

//     // Setup some empty matrices to hold the mean vector and error matrices
//     auto meanSum = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
//     auto errorMatrixSum = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

//     // Loop over the universes
//     unsigned int nValidUniverses = 0u;
//     for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
//     {
//         //// BEGIN DEBUG
//         AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
//         //// END DEBUG

//         // Get the value in this universe
//         const auto pUniverse = func(iUni);

//         // Insist that we get a valid quanitiy in this universes
//         if (!pUniverse)
//             continue;

//         // Count the number of valid universes
//         nValidUniverses++;

//         // Loop over the bins
//         // ATTN for speed we don't expcitly check here if the universe has the same number of bins as the nominal. Instead we rely on the
//         // internal range checking of the ubsmear::UBMatrix class
//         for (unsigned int iBin = 0; iBin < nBins; ++iBin)
//         {
//             // Add up the universes
//             meanSum.SetElement(iBin, 0, meanSum.At(iBin, 0) + pUniverse->At(iBin, 0));
//             // Loop over the bins again
//             // ATTN the error matrix is symmetric so here we only set one half of the off-diagonals within the universe loop and then copy
//             // them over to the other half afterwards. This is done for performance reasons as it halves the number of inserts required
//             const auto diffI = pUniverse->At(iBin, 0) - pNominal->At(iBin, 0);
//             for (unsigned int jBin = 0; jBin <= iBin; ++jBin)
//             {
//                 const auto diffJ = pUniverse->At(jBin, 0) - pNominal->At(jBin, 0);

//                 // Add up the error matrix elements
//                 errorMatrixSum.SetElement(iBin, jBin, errorMatrixSum.At(iBin, jBin) + diffI*diffJ);
//             }
//         }
//     }

//     std::cout << "DEBUG - Valid universes: " << nValidUniverses << " / " << nUniverses << std::endl;

//     // Scale the sums by the number of universes
//     if (nValidUniverses == 0)
//         throw std::logic_error("CrossSection::GetDisributionParams - Desired quantity was invalid in all universes");

//     const auto scaleFactor = 1.f / static_cast<float>(nValidUniverses);
//     const auto mean = meanSum * scaleFactor;
//     const auto errorMatrix = errorMatrixSum * scaleFactor;

//     // Get the bias vector
//     const auto pBias = std::make_shared<ubsmear::UBMatrix>(mean - *pNominal);

//     // Get the covariance matrix
//     auto pCovarianceMatrix = std::make_shared<ubsmear::UBMatrix>(ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins));
//     for (unsigned int iBin = 0; iBin < nBins; ++iBin)
//     {
//         for (unsigned int jBin = 0; jBin <= iBin; ++jBin)
//         {
//             // ATTN here we use the fact that V_ij = E_ij - b_i*b_j
//             //   - E_ij = error matrix element
//             //   - V_ij = covariance matrix element
//             //   - b_i = bias vector element
//             //
//             // This is done instead of calculating the covariance directly as this only requires one pass through the universe loop.
//             const auto covariance = errorMatrix.At(iBin, jBin) - pBias->At(iBin, 0)*pBias->At(jBin, 0);

//             pCovarianceMatrix->SetElement(iBin, jBin, covariance);
//             pCovarianceMatrix->SetElement(jBin, iBin, covariance);
//         }
//     }

//     return { pBias, pCovarianceMatrix };
// }


// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetPredictedCrossSectionStatUncertainty(const ScalingData &scalingData) const
{
    // Get the number of universes for this systematic parameter
    // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
    const auto nUniverses = m_systParams.nBootstrapUniverses;

    // Get the nominal cross-section
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>(this->GetPredictedCrossSection(scalingData));

    // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(pXSecNom->GetRows(), 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // The bootstrap universes
    const auto universes = m_signal_true_multisims.at("misc").at("bootstrap");

    // Get the distribution parameters
    return this->GetDistributionParams([&](const unsigned int &iUni)
    {
        // This is the function we call to get the value in each universe

        // Get the number of signal events in the nominal simulation in true bins
        const auto signal = CrossSectionHelper::GetMatrixFromHist(universes.at(iUni));

        // Get the cross-section
        return std::make_shared<ubsmear::UBMatrix> ( this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets) );

    }, nUniverses, pXSecNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetPredictedSidebandCrossSectionStatUncertainty(const ScalingData &scalingData) const
{
    // Get the number of universes for this systematic parameter
    // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
    const auto nUniverses = m_systParams.nBootstrapUniverses;

    // Get the nominal cross-section
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>(this->GetPredictedCrossSection(scalingData));

    // // Get the number of backgrounds is zero as we are effectively applying a "perfect" selection
    const auto zeroVector = ubsmear::UBMatrixHelper::GetZeroMatrix(pXSecNom->GetRows(), 1);

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // The bootstrap universes
    const auto universesSignal = m_signal_true_multisims.at("misc").at("sidebandWeights");

    // // The bootstrap universes
    // const auto universesBackground = m_signal_true_multisims.at("misc").at("sidebandWeights");

    // Get the distribution parameters
    return this->GetDistributionParams([&](const unsigned int &iUni)
    {
        // This is the function we call to get the value in each universe

        // Get the number of signal events in the nominal simulation in true bins
        const auto signal = CrossSectionHelper::GetMatrixFromHist(universesSignal.at(iUni));
        // const auto background = CrossSectionHelper::GetMatrixFromHist(universesBackground.at(iUni));

        // Get the cross-section
        return std::make_shared<ubsmear::UBMatrix> ( this->GetCrossSection(signal, zeroVector, integratedFlux, scalingData.exposurePOT, scalingData.nTargets) );

    }, nUniverses, pXSecNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::map<std::string, ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetTotalErrorMatrixMap(const ScalingData &scalingData) const
{
    std::map<std::string, ubsmear::UBMatrix> totalErrorMatrixMap;
    const auto dataStatUncertainties = this->GetBNBDataCrossSectionStatUncertainty(scalingData);
    const auto dataSystBiasCovariances = this->GetBNBDataCrossSectionSystUncertainties(scalingData);
    const auto smearingMatrixSystBiasCovariances = this->GetSmearingMatrixSystUncertainties();

    std::map<std::string, std::map<std::string, ubsmear::UBMatrix> > errorMatrixMap;

    // -----------------------------------------------------------------------------------------------------------------------------
    // Read the matrices from disk
    // -----------------------------------------------------------------------------------------------------------------------------
    // Loop over the two quantities that have systematic uncertainties
    const auto data = this->GetSelectedBNBDataEvents();
    const auto smearingMatrix = this->GetSmearingMatrix();
    for (const std::string &quantity : {"data", "smearingMatrix"})
    {
        std::cout<<"data.GetRows(): "<<data.GetRows()<<" - smearingMatrix.GetRows(): "<<smearingMatrix.GetRows()<<std::endl;
        // Get the number of bins (for the smearing matrix there are N^2 bins when flattened)
        const auto nBins = (quantity == "data" ? data.GetRows() : std::pow(smearingMatrix.GetRows(), 2)); //todo find a better way to do this (just use getnbins() ???)

        // Setup an empty error matrix for this quantity
        auto errorMatrixTotalSum =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

        // For the data cross-section, we also have a stat uncertainty
        // For the smearing matrix, just use a zero vector
        const auto statUncertainties = (quantity == "data") ? dataStatUncertainties : ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);

        // To convert these stat uncertainties into an error matrix, we produce a diagonal matrix whose diagonal entries contain the
        // variances (i.e. square of the stat uncerainties)
        const auto statVariances = ubsmear::ElementWiseOperation(statUncertainties, statUncertainties, [](const auto &l, const auto &r) { return l * r; });
        const auto statErrorMatrix = ubsmear::UBMatrixHelper::GetDiagonalMatrix(statVariances);

        // Store this in the map
        errorMatrixMap[quantity].emplace("stat", statErrorMatrix);

        // Add this to the grand total error matrix
        errorMatrixTotalSum = errorMatrixTotalSum + statErrorMatrix;

        // Handle the multisim parameters
        for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystDimensionsMap>(
            {
                {"flux", m_systParams.fluxDimensions},
                {"xsec", m_systParams.xsecDimensions},
                {"reint", m_systParams.reintDimensions},
                {"misc", {
                    {"bootstrap", m_systParams.nBootstrapUniverses},
                    {"dirt", 2},
                    {"POT", 0}
                }}
            }))
        {
            // Setup an empty error matrix for this group
            auto errorMatrixTotal =  ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

            // Loop over the parameters in this group
            for (const auto &[paramName, nUniverses] : dimensions)
            {
                // Get the bias vector and covariance matrix
                const auto &[pBias, pCovariance] =  (quantity == "data") ? dataSystBiasCovariances.at(group).at(paramName) : smearingMatrixSystBiasCovariances.at(group).at(paramName);
                // getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_NuWro");
                const auto biasVector = *pBias;
                // getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_NuWro");
                const auto covarianceMatrix = *pCovariance;

                // Get the total error matrix from the bias and covariance
                const auto errorMatrix = CrossSectionHelper::GetErrorMatrix(biasVector, covarianceMatrix);

                // Add this error matrix to the total
                errorMatrixTotal = errorMatrixTotal + errorMatrix;
            }

            // Store this total in the map
            errorMatrixMap[quantity].emplace(group, errorMatrixTotal);

            // Add this total to the grand total error matrix
            errorMatrixTotalSum = errorMatrixTotalSum + errorMatrixTotal;
        }

        // Handle the unisim parameters
        for (const auto &[group, dimensions] : std::map<std::string, CrossSectionHelper::SystUnisimDimensionsMap>(
            {
                {"detector", m_systParams.detVarDimensions}
            }))
        {
            // Setup an empty error matrix for this group
            auto errorMatrixTotal = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

            // Loop over the parameters in this group
            for (const auto &[paramName, cvName] : dimensions)
            {
                // Get the bias vector and (dummy) covariance matrix
                const auto &[pBias, pCovariance] =  (quantity == "data") ? dataSystBiasCovariances.at(group).at(paramName) : smearingMatrixSystBiasCovariances.at(group).at(paramName);
                // getMatrixFunction(quantity + "_" + group + "_" + paramName + "_bias_NuWro");
                const auto biasVector = *pBias;
                // getMatrixFunction(quantity + "_" + group + "_" + paramName + "_covariance_NuWro");
                const auto covarianceMatrix = *pCovariance;

                // Get the total error matrix from the bias and covariance
                const auto errorMatrix = CrossSectionHelper::GetErrorMatrix(biasVector, covarianceMatrix);

                // Add this error matrix to the total
                errorMatrixTotal = errorMatrixTotal + errorMatrix;
            }

            // Store this total in the map
            errorMatrixMap[quantity].emplace(group, errorMatrixTotal);

            // Add this total to the grand total error matrix
            errorMatrixTotalSum = errorMatrixTotalSum + errorMatrixTotal;
        }

        // Add the grand total to the map
        totalErrorMatrixMap.emplace(quantity, errorMatrixTotalSum);
    }
    return totalErrorMatrixMap;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSectionInUniverse(const std::string &group, const std::string &paramName, const unsigned int universeIndex, const ScalingData &scalingData) const
{
    // ATTN for speed, this function (and other similar functions) doesn't explicitly check if the input group, paramName, universeIndex is
    // valid. Instead we rely on the internal range checking of the .at() method for STL containers. This is for peformance reasons as this
    // function gets called once per universe!

    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get the number of predicted backgrounds in the supplied universe
    const auto backgrounds = CrossSectionHelper::GetMatrixFromHist(m_background_selected_reco_multisims.at(group).at(paramName).at(universeIndex));

    // Get the integrated flux in the supplied universe (if it's not a flux parameter, then use the nominal universe)
    const auto integratedFlux = (
        group == "flux"
            ? scalingData.pFluxReweightor->GetIntegratedFluxVariation(paramName, universeIndex)
            : scalingData.pFluxReweightor->GetIntegratedNominalFlux()
    );

    if(group == "misc" && paramName == "sidebandWeights") // todo remove
    {
        std::cout<<"selected:";
        for(unsigned int i=0; i<selected.GetRows(); i++) std::cout<<" "<<selected.At(i,0);
        std::cout<<"\nbackgrounds:";
        for(unsigned int i=0; i<backgrounds.GetRows(); i++) std::cout<<" "<<backgrounds.At(i,0);
        std::cout<<"\nintegratedFlux: "<<integratedFlux<<std::endl;
    }

    // Get the cross-section
    return this->GetCrossSection(selected, backgrounds, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::CrossSection::GetBNBDataCrossSectionForUnisim(const std::string &group, const std::string &paramName, const ScalingData &scalingData) const
{
    // Get the number of selected events in BNB data
    const auto selected = CrossSectionHelper::GetMatrixFromHist(m_pBNBData_selected_reco);

    // Get the number of predicted backgrounds in the supplied universe
    const auto backgrounds = CrossSectionHelper::GetMatrixFromHist(m_background_selected_reco_unisims.at(group).at(paramName));

    // Get the integrated flux in the nominal universe
    const auto integratedFlux = scalingData.pFluxReweightor->GetIntegratedNominalFlux();

    // Get the cross-section
    return this->GetCrossSection(selected, backgrounds, integratedFlux, scalingData.exposurePOT, scalingData.nTargets);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetDistributionParams(const std::function<std::shared_ptr<ubsmear::UBMatrix>(const unsigned int)> &func, const unsigned int nUniverses, const std::shared_ptr<ubsmear::UBMatrix> &pNominal) const
{
    if (nUniverses == 0)
        throw std::logic_error("CrossSection::GetDisributionParams - No universes supplied");

    // Insist the nominal input is a column vector
    if (pNominal->GetColumns() != 1)
        throw std::logic_error("CrossSection::GetDisributionParams - Input nominal is not a column vector");

    // Get the number of bins
    const auto nBins = pNominal->GetRows();

    // Setup some empty matrices to hold the mean vector and error matrices
    auto meanSum = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, 1);
    auto errorMatrixSum = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);

    // Loop over the universes
    unsigned int nValidUniverses = 0u;
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        //// BEGIN DEBUG
        AnalysisHelper::PrintLoadingBar(iUni, nUniverses);
        //// END DEBUG

        // Get the value in this universe
        const auto pUniverse = func(iUni);


        // Insist that we get a valid quanitiy in this universe
        if (!pUniverse)
        {
            std::cout<<"Debug GetDistributionParams iUni: "<<iUni<<" skipped."<<std::endl; //todo remove this cout
            continue;
        }

        // Count the number of valid universes
        nValidUniverses++;

        // Loop over the bins
        // ATTN for speed we don't expcitly check here if the universe has the same number of bins as the nominal. Instead we rely on the
        // internal range checking of the ubsmear::UBMatrix class
        for (unsigned int iBin = 0; iBin < nBins; ++iBin)
        {
            // Add up the universes
            meanSum.SetElement(iBin, 0, meanSum.At(iBin, 0) + pUniverse->At(iBin, 0));
            // Loop over the bins again
            // ATTN the error matrix is symmetric so here we only set one half of the off-diagonals within the universe loop and then copy
            // them over to the other half afterwards. This is done for performance reasons as it halves the number of inserts required
            const auto diffI = pUniverse->At(iBin, 0) - pNominal->At(iBin, 0);
            for (unsigned int jBin = 0; jBin <= iBin; ++jBin)
            {
                const auto diffJ = pUniverse->At(jBin, 0) - pNominal->At(jBin, 0);
                if((!std::isfinite(diffI)) || (!std::isfinite(diffJ)))
                {
                    std::cout<<"Debug GetDistributionParams uni jBin "<<jBin<<": "<<pUniverse->At(jBin, 0)<<" - nom: "<< pNominal->At(jBin, 0)<<" - uni iBin "<<iBin<<": "<<pUniverse->At(iBin, 0)<<" - nom: "<< pNominal->At(iBin, 0)<<" - diffI"<<diffI<<" - diffJ"<<diffJ<<std::endl;
                }
                // Add up the error matrix elements
                errorMatrixSum.SetElement(iBin, jBin, errorMatrixSum.At(iBin, jBin) + diffI*diffJ);
            }
        }
    }

    
    // Scale the sums by the number of universes
    if (nValidUniverses == 0)
        throw std::logic_error("CrossSection::GetDisributionParams - Desired quantity was invalid in all universes");

    const auto scaleFactor = 1.f / static_cast<float>(nValidUniverses);
    const auto mean = meanSum * scaleFactor;
    const auto errorMatrix = errorMatrixSum * scaleFactor;

    std::cout<<"Debug GetDistributionParams errorMatrix:";
    errorMatrix.Print();

    std::cout<<"Debug GetDistributionParams mean:";
    for(unsigned int i=0; i< mean.GetRows(); i++) 
        std::cout<<" "<<mean.At(i,0);    
    std::cout<<std::endl; //todo remove these cout

    // Get the bias vector
    const auto pBias = std::make_shared<ubsmear::UBMatrix>(mean - *pNominal);

    std::cout<<"Debug GetDistributionParams bias:";
    for(unsigned int i=0; i< pBias->GetRows(); i++) 
        std::cout<<" "<<pBias->At(i,0);

    // Get the covariance matrix
    auto pCovarianceMatrix = std::make_shared<ubsmear::UBMatrix>(ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins));
    for (unsigned int iBin = 0; iBin < nBins; ++iBin)
    {
        for (unsigned int jBin = 0; jBin <= iBin; ++jBin)
        {
            // ATTN here we use the fact that V_ij = E_ij - b_i*b_j
            //   - E_ij = error matrix element
            //   - V_ij = covariance matrix element
            //   - b_i = bias vector element
            //
            // This is done instead of calculating the covariance directly as this only requires one pass through the universe loop.
            const auto covariance = errorMatrix.At(iBin, jBin) - pBias->At(iBin, 0)*pBias->At(jBin, 0);

            pCovarianceMatrix->SetElement(iBin, jBin, covariance);
            pCovarianceMatrix->SetElement(jBin, iBin, covariance);
        }
    }

    std::cout<<"Debug GetDistributionParams covarianceMatrix:";
    pCovarianceMatrix->Print();

    return { pBias, pCovarianceMatrix };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetDistributionParamsUnisim(const std::shared_ptr<ubsmear::UBMatrix> &pVaried, const std::shared_ptr<ubsmear::UBMatrix> &pCentralValue, const std::shared_ptr<ubsmear::UBMatrix> &pNominal) const
{
    // Check the input matrices are valid
    if (!pVaried || !pCentralValue || !pNominal)
        throw std::invalid_argument("CrossSection::GetDistributionParamsUnisim - One or more of the inputs are not valid");
    
    // Check ths input matrices have the desired dimensions
    if (pVaried->GetColumns() != 1 || pCentralValue->GetColumns() != 1 || pNominal->GetColumns() != 1)
        throw std::invalid_argument("CrossSection::GetDistributionParamsUnisim - One or more of the inputs are not a column vector");

    const auto nBins = pVaried->GetRows();
    if (pCentralValue->GetRows() != nBins || pNominal->GetRows() != nBins)
        throw std::invalid_argument("CrossSection::GetDistributionParamsUnisim - Input column vectors have different sizes");

    // Get the fractional bias of the cross-section away from the central value
    const auto fracBias = ubsmear::ElementWiseOperation(
        (*pVaried - *pCentralValue), *pCentralValue,
        [](const float numerator, const float denominator)
        {
            // ATTN here we decide what to do if there's a bin of the central-value that's zero
            // In this case we don't have any information, so here we set the fractional bias to 0.
            if (std::abs(denominator) <= std::numeric_limits<float>::epsilon())
                return 0.f;

            return numerator / denominator;
        }
    );

    // Apply the fractional bias to the nominal cross-section
    const auto bias = ubsmear::ElementWiseOperation(
        fracBias, *pNominal,
        [](const float lhs, const float rhs) { return lhs * rhs; }
    );

    // Return the result (use a zero matrix for the covariance)
    return {
        std::make_shared<ubsmear::UBMatrix>(bias),
        std::make_shared<ubsmear::UBMatrix>(ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins))
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetDistributionParamsNormalisation(const std::shared_ptr<ubsmear::UBMatrix> &pNominal, const float fracUncertainty) const
{
    // Check the input matrix are valid
    if (!pNominal)
        throw std::invalid_argument("CrossSection::GetDistributionParamsNormalisation - The input nominal is not valid");

    // Check ths input matrix has the desired dimensions
    if (pNominal->GetColumns() != 1)
        throw std::invalid_argument("CrossSection::GetDistributionParamsNormalisation - The input is not a column vector");

    const auto nBins = pNominal->GetRows();

    // Scale the nominal value by the fractional uncertainty to get the bias, and use a zero matrix for the covariance
    return {
        std::make_shared<ubsmear::UBMatrix>((*pNominal) * fracUncertainty),
        std::make_shared<ubsmear::UBMatrix>(ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins))
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetBNBDataCrossSectionDistributionParams(const std::string &group, const std::string &paramName, const ScalingData &scalingData) const
{
    std::cout << "\n DEBUG - Processing parameter: " << group << ", " << paramName << "\n" << std::endl;

    // Get the number of universes for this systematic parameter
    // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
    const auto nUniverses = m_signal_true_multisims.at(group).at(paramName).size();

    // Get the nominal cross-section
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>(this->GetBNBDataCrossSection(scalingData));

    // Get the distirbution parameters
    return this->GetDistributionParams([&](const unsigned int &iUni)
    {
        // This is the function we call to get the value in each universe
        return std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSectionInUniverse(group, paramName, iUni, scalingData) );

    }, nUniverses, pXSecNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetBNBDataCrossSectionDistributionParamsUnisim(const std::string &group, const std::string &paramName, const std::string &cvName, const ScalingData &scalingData) const
{
    std::cout << "DEBUG - Processing parameter: " << group << ", " << paramName << " (with central value: " << cvName << ")" << std::endl;

    // Get the cross-section using the unisim sample, in the central value sample and in the nominal universe
    const auto pXSec = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSectionForUnisim(group, paramName, scalingData) );
    const auto pXSecCV = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSectionForUnisim(group, cvName, scalingData) );
    const auto pXSecNom = std::make_shared<ubsmear::UBMatrix>( this->GetBNBDataCrossSection(scalingData) );

    // Get the parameters
    return this->GetDistributionParamsUnisim(pXSec, pXSecCV, pXSecNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------s

std::unordered_map<std::string, CrossSectionHelper::SystTH1FMap> CrossSectionHelper::CrossSection::GetSelectedBackgroundRecoMap() const
{
    return m_background_selected_reco_multisims;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::unordered_map<std::string, CrossSectionHelper::SystTH2FMap> CrossSectionHelper::CrossSection::GetSelectedSignalRecoTruthMap() const
{
    return m_signal_selected_recoTrue_multisims;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixInUniverse(const std::string &group, const std::string &paramName, const unsigned int universeIndex) const
{
    return this->GetSmearingMatrix(m_signal_true_multisims.at(group).at(paramName).at(universeIndex), m_signal_selected_recoTrue_multisims.at(group).at(paramName).at(universeIndex));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::CrossSection::GetSmearingMatrixForUnisim(const std::string &group, const std::string &paramName) const
{
    return this->GetSmearingMatrix(m_signal_true_unisims.at(group).at(paramName), m_signal_selected_recoTrue_unisims.at(group).at(paramName));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetSmearingMatrixDistributionParams(const std::string &group, const std::string &paramName) const
{
    std::cout << "DEBUG - Smearing matrix. Processing parameter: " << group << ", " << paramName << std::endl;

    // Get the number of universes for this systematic parameter
    // ATTN this choice of m_signal_true_multisims here is arbitrary, any of the mutlisims maps would do
    const auto nUniverses = m_signal_true_multisims.at(group).at(paramName).size();

    // Get the nominal smearing matrix
    const auto pSmearingNom = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixNominal());
    if (!pSmearingNom)
    {
        std::cout<<"CrossSection::GetSmearingMatrixDistributionParams - The nominal smearing matrix can't be calculated"<<std::endl;
        throw std::logic_error("CrossSection::GetSmearingMatrixDistributionParams - The nominal smearing matrix can't be calculated");
    }
    // Get the distirbution parameters
    return this->GetDistributionParams([&](const unsigned int &iUni)
    {
        // This is the function we call to get the value in each universe
        auto flattened = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixInUniverse(group, paramName, iUni));
        return flattened;

    }, nUniverses, pSmearingNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystBiasCovariancePair CrossSectionHelper::CrossSection::GetSmearingMatrixDistributionParamsUnisim(const std::string &group, const std::string &paramName, const std::string &cvName) const
{
    std::cout << "DEBUG - Smearing matrix. Processing parameter Unisim: " << group << ", " << paramName << " (with central value: " << cvName << ")" << std::endl;

    // Get the smearing matrix using the unisim sample, in the central value sample and in the nominal universe
    const auto pSmearing = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixForUnisim(group, paramName));
    std::cout << "DEBUG - Smearing matrix. Processing parameter Unisim Point 1"<<std::endl;
    const auto pSmearingCV = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixForUnisim(group, cvName));
    std::cout << "DEBUG - Smearing matrix. Processing parameter Unisim Point 2"<<std::endl;
    const auto pSmearingNom = CrossSectionHelper::FlattenMatrix(this->GetSmearingMatrixNominal());

    std::cout << "DEBUG - Smearing matrix. Processing parameter Unisim Point 3"<<std::endl;
    // Get the parameters
    return this->GetDistributionParamsUnisim(pSmearing, pSmearingCV, pSmearingNom);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH1F> CrossSectionHelper::GetTH1F(const std::vector<float> &binEdges)
{
    return std::make_shared<TH1F>(("xSecHist_" + std::to_string(m_histCount++)).c_str(), "", binEdges.size() - 1, binEdges.data());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TH2F> CrossSectionHelper::GetTH2F(const std::vector<float> &binEdges, const std::vector<float> &trueBinEdges)
{
    const auto binEdges2 = trueBinEdges.empty() ? binEdges : trueBinEdges;
    return std::make_shared<TH2F>(("xSecHist_" + std::to_string(m_histCount++)).c_str(), "", binEdges.size() - 1, binEdges.data(), binEdges2.size() - 1, binEdges2.data()); //todo check order is correct
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystTH1FMap CrossSectionHelper::GetSystTH1FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions)
{
    SystTH1FMap map;

    for (const auto &[paramName, nUniverses] : dimensions)
    {
        std::vector< std::shared_ptr<TH1F> > universes;
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            universes.push_back(CrossSectionHelper::GetTH1F(binEdges));
        }
        map.emplace(paramName, universes);
    }
    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystUnisimTH1FMap CrossSectionHelper::GetSystUnisimTH1FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions)
{
    SystUnisimTH1FMap map;

    // Setup a vector to keep track of the used central-values sample names
    std::vector<std::string> allCVNames;

    for (const auto &[paramName, cvName] : dimensions)
    {
        // Include all parameter names
        map.emplace(paramName, CrossSectionHelper::GetTH1F(binEdges));

        // If we haven't seen this central value sample name before, then add it too!
        if (std::find(allCVNames.begin(), allCVNames.end(), cvName) == allCVNames.end())
        {
            if (map.find(cvName) != map.end())
                throw std::invalid_argument("CrossSectionHelper::GetSystUnisimTH1FMap - Name \"" + cvName + "\" is describing a parameter and a central-value sample");

            map.emplace(cvName, CrossSectionHelper::GetTH1F(binEdges));
            allCVNames.push_back(cvName);
        }
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystTH2FMap CrossSectionHelper::GetSystTH2FMap(const std::vector<float> &binEdges, const SystDimensionsMap &dimensions, const std::vector<float> &trueBinEdges)
{
    SystTH2FMap map;
    const auto binEdges2 = trueBinEdges.empty() ? binEdges : trueBinEdges;
    for (const auto &[paramName, nUniverses] : dimensions)
    {
        std::vector< std::shared_ptr<TH2F> > universes;
        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            universes.push_back(CrossSectionHelper::GetTH2F(binEdges, binEdges2));
        }

        map.emplace(paramName, universes);
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystUnisimTH2FMap CrossSectionHelper::GetSystUnisimTH2FMap(const std::vector<float> &binEdges, const SystUnisimDimensionsMap &dimensions, const std::vector<float> &trueBinEdges)
{
    SystUnisimTH2FMap map;

    // Setup a vector to keep track of the used central-values sample names
    std::vector<std::string> allCVNames;
    const auto binEdges2 = trueBinEdges.empty() ? binEdges : trueBinEdges;

    for (const auto &[paramName, cvName] : dimensions)
    {
        // Include all parameter names
        map.emplace(paramName, CrossSectionHelper::GetTH2F(binEdges, binEdges2));

        // If we haven't seen this central value sample name before, then add it too!
        if (std::find(allCVNames.begin(), allCVNames.end(), cvName) == allCVNames.end())
        {
            if (map.find(cvName) != map.end())
                throw std::invalid_argument("CrossSectionHelper::GetSystUnisimTH2FMap - Name \"" + cvName + "\" is describes a parameter and a central-value sample");

            map.emplace(cvName, CrossSectionHelper::GetTH2F(binEdges, binEdges2));
            allCVNames.push_back(cvName);
        }
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::GetUnitWeightsMap(const SystDimensionsMap &dimensions)
{
    // Make a new map for the output
    SystFloatMap map;

    // Loop over all parameters in the input dimensions map
    for (const auto &[paramName, nUniverses] : dimensions)
    {
        // Add a weight of 1.f for each universe
        map.emplace(paramName, std::vector<float>(nUniverses, 1.f));
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::ScaleWeightsMap(const SystFloatMap &weightsMap, const float &divisor, const std::unordered_map<std::string, bool> scalingMap)
{
    // Check the input divisor isn't zero
    if (std::abs(divisor) <= std::numeric_limits<float>::epsilon())
        throw std::invalid_argument("CrossSectionHelper::ScaleWeightsMap - The input divisor is zero");

    auto map = weightsMap;
    for (auto &[paramName, weights] : map)
    {
        if(scalingMap.empty() || scalingMap.at(paramName))
        {
            for (auto &weight : weights)
            {
                weight /= divisor;
            }
        }
    }

    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions)
{
    if (!truth.systParamNames.IsSet() || !truth.systParamFirstValueIndex.IsSet() || !truth.systParamValues.IsSet())
        throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Systematic parameters are not set in the input truth information");

    const auto &systParamNames = truth.systParamNames();
    const auto &systParamFirstValueIndex = truth.systParamFirstValueIndex();
    const auto &systParamValues = truth.systParamValues();

    // Get the total number of parameters available
    const auto nParameters = systParamNames.size();

    // Make a new map for the output
    SystFloatMap map;

    // Loop over all parameters in the input dimensions map
    for (const auto &[paramName, nUniverses] : dimensions)
    {
        // Find the desired name
        const auto iter = std::find(systParamNames.begin(), systParamNames.end(), paramName);
        if (iter == systParamNames.end())
        {
            std::cout<<"CrossSectionHelper::GetWeightsMap - Unknown parameter: "<<paramName<<std::endl;
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Unknown parameter: " + paramName);
        }
        // Get the index of the requested parameter
        const unsigned int index = std::distance(systParamNames.begin(), iter);

        // Get the first and last value in the weights vector
        const unsigned int firstValueIndex = systParamFirstValueIndex.at(index);
        const unsigned int lastValueIndex = ((index == nParameters - 1u) ? systParamValues.size() : systParamFirstValueIndex.at(index + 1u));

        // Pick out the weights corresponding the the desired parameter
        const std::vector<float> weights(std::next(systParamValues.begin(), firstValueIndex), std::next(systParamValues.begin(), lastValueIndex));

        // Check that we have the right number of weights
        if (weights.size() != nUniverses)
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Number of weights for parameter " + paramName + " doesn't match the input dimensions");

        // Setup the output map entry
        auto &outputWeights = map[paramName];

        // Fill the output weights
        for (const auto &weight : weights)
        {
            // ATTN sometimes we have non-physical weights in the input. Here we decide what to do with them!

            // If the weight is non-negative and non-infinite then it's okay!
            if (weight >= 0.f && std::abs(std::abs(weight) - std::numeric_limits<float>::max()) > std::numeric_limits<float>::epsilon() && std::isfinite(weight))
            {
                outputWeights.push_back(weight);
                continue;
            }

            // Force negative weights to be zero and infinite weights to be unity
            outputWeights.push_back((weight < 0.f) ? 0.f : 1.f);
        }
    }
    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

Double_t CrossSectionHelper::CrossSection::GetSidebandScaling(const float trueSidebandValue, const std::vector<Double_t> &cc0piNominalConstraintParam, const std::vector<float> &binEdges) const
{
    for (unsigned int b=0; b < binEdges.size()-1; b++)//Todo: Do this more efficiently
    {
        if (trueSidebandValue>=binEdges[b] && trueSidebandValue<binEdges[b+1])
        {
            const auto paramNominal = cc0piNominalConstraintParam[b];//std::max(, 0.01);//cc0piNominalConstraintParam[b];
            if(paramNominal<0.0)
            {
                std::cout<<"ExtractXSecs - Negative cc0piNominalConstraintParam: "<<paramNominal<<std::endl;
                throw std::logic_error("ExtractXSecs - Negative cc0piNominalConstraintParam.");
            }
            return paramNominal;
        }
    }
    std::cout<<"GetSidebandScaling - trueSidebandValue "<< trueSidebandValue <<" outside of bins with range: "<<binEdges[0]<<" - "<<binEdges[binEdges.size()-1]<<std::endl;
    throw std::logic_error("GetSidebandScaling - trueSidebandValue outside of bin range.");
    return -1;//Todo improve code
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::CrossSection::GetSidebandUniverseScaling(const CrossSectionHelper::SystFloatMap weights, const float trueSidebandValue, const std::map<std::string, std::vector<std::pair<std::vector<Double_t>,std::vector<Double_t>>>> &cc0piUniverseConstraints, const std::vector<float> &binEdges) const
{
    auto weightsScaled = weights;
    for (unsigned int b=0; b < binEdges.size()-1; b++)//Todo: Do this more efficiently
    {
        if (trueSidebandValue>=binEdges[b] && trueSidebandValue<binEdges[b+1])
        {
            // const auto paramNominal = cc0piNominalConstraintParam[b];//std::max(, 0.01);//cc0piNominalConstraintParam[b];
            // if(paramNominal<0.0)
            // {
            //     std::cout<<"ExtractXSecs - Negative cc0piNominalConstraintParam: "<<paramNominal<<std::endl;
            //     throw std::logic_error("ExtractXSecs - Negative cc0piNominalConstraintParam.");
            // }

            for (auto &[paramName,weightVector] : weightsScaled)
            {
                const auto cc0piUniverseConstraintParam = cc0piUniverseConstraints.at(paramName);
                if(cc0piUniverseConstraintParam.size()!= weightVector.size())
                {
                    std::cout<<"GetSidebandUniverseScaling - Wrong number of parameters.."<<std::endl;
                    throw std::logic_error("GetSidebandUniverseScaling - Wrong number of parameters.");
                }
                
                for (unsigned int u=0; u<weightVector.size(); u++)
                {
                    // The smearing matrix in a universe could contain zeros. In these cases no parameters can be computed and all are set to -1. Use only nominal values instead.
                    if(cc0piUniverseConstraintParam.at(u).first.at(b)>=0)
                    {
                        weightVector[u] *= cc0piUniverseConstraintParam.at(u).first.at(b);///paramNominal;
                    }
                }
            }
            return weightsScaled;
        }
    }
    std::cout<<"GetSidebandUniverseScaling - trueSidebandValue "<< trueSidebandValue <<" outside of bins with range: "<<binEdges[0]<<" - "<<binEdges[binEdges.size()-1]<<std::endl;
    throw std::logic_error("GetSidebandUniverseScaling - trueSidebandValue outside of bin range.");//Todo improve code
}

// -----------------------------------------------------------------------------------------------------------------------------------------
std::vector<float> CrossSectionHelper::CrossSection::GetSidebandParameterWeights(const float trueSidebandValue, const std::vector<Double_t> &cc0piNominalConstraintParam, const std::vector<Double_t> &cc0piNominalConstraintParamError, const std::vector<float> &binEdges) const
{
    for (unsigned int b=0; b < binEdges.size()-1; b++)//Todo: Do this more efficiently
    {
        if (trueSidebandValue>=binEdges[b] && trueSidebandValue<binEdges[b+1])
        {
            const auto paramNominal = cc0piNominalConstraintParam[b];//std::max(, 0.01);//cc0piNominalConstraintParam[b];
            const auto paramNominalError = cc0piNominalConstraintParamError[b];
            // std::cout<<"GetSidebandParameterWeights: "<<trueSidebandValue<<"("<<paramNominal<<"+-"<<paramNominalError<<"), "<<std::endl;
            if(paramNominal<0.0)
            {
                // std::cout<<"CrossSection::GetSidebandParameterWeights - Negative value for paramNominal. This means fit failed in universe and weight has been set to -1. Check results!"<<std::endl;
                throw std::logic_error("CrossSection::GetSidebandParameterWeights - Negative paramNominal.");
            } else if(paramNominal<std::numeric_limits<float>::epsilon())
            {
                std::cout<<"CrossSection::GetSidebandParameterWeights - cc0piNominalConstraintParam too small: "<<paramNominal<<std::endl;
                // throw std::logic_error("CrossSection::GetSidebandParameterWeights - cc0piNominalConstraintParam too small.");
            }

            const auto sidebandWeights = (
                paramNominal>std::numeric_limits<float>::epsilon()
                    ? GenerateNormalWeights(m_systParams.nBootstrapUniverses, paramNominalError, paramNominal, 1.f)//1.f/paramNominal) // todo: verify this is the correct way
                    : std::vector<float>(m_systParams.nBootstrapUniverses, 1.0));
            // const auto sidebandWeights = GenerateNormalWeights(m_systParams.nBootstrapUniverses, paramNominalError); // todo: verify this is the correct way (no /paramNominal)
            // std::cout<<"GetSidebandParameterWeights - return"<<std::endl; 
            return sidebandWeights;
        }
    }
    std::cout<<"GetSidebandParameterWeights - trueSidebandValue "<< trueSidebandValue <<" outside of bins with range: "<<binEdges[0]<<" - "<<binEdges[binEdges.size()-1]<<std::endl;
    throw std::logic_error("GetSidebandParameterWeights - trueSidebandValue outside of bin range.");
    return std::vector<float>(m_systParams.nBootstrapUniverses, 1.0);//Todo improve code
}

// -----------------------------------------------------------------------------------------------------------------------------------------

CrossSectionHelper::SystFloatMap CrossSectionHelper::GetWeightsMap(const Event::Truth &truth, const SystDimensionsMap &dimensions, const SystMutuallyExclusiveDimensionsMap &mutuallyExclusiveDimensions)
{
    // Make a new map for the output
    SystFloatMap map;

    // Make new dimensions map to populate with the parameters that are available in the input truth information directly
    SystDimensionsMap availableParamDimensions;

    const auto &systParamNames = truth.systParamNames();

    // Find which parameters (if any) are built from mutually exclusive parameters
    for (const auto &[paramName, nUniverses] : dimensions)
    {
        // Check if this parameter name is available in the input truth information, if so just skip it for later
        if (std::find(systParamNames.begin(), systParamNames.end(), paramName) != systParamNames.end())
        {
            availableParamDimensions.emplace(paramName, nUniverses);
            continue;
        }

        // Insist that parameters not available in the input truth information are listed in the mutually exlusive dimensions
        if (mutuallyExclusiveDimensions.find(paramName) == mutuallyExclusiveDimensions.end())
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Parameter \"" + paramName + "\" isn't listed in the input truth information or in the supplied mutually exclusive parameter names");

        // We have a mutually exclusive parameter, so get it's weights
        const auto &[parameters, nUniversesCheck] = mutuallyExclusiveDimensions.at(paramName);

        if (nUniverses != nUniversesCheck)
            throw std::invalid_argument("CrossSectionHelper::GetWeightsMap - Inconsistent number of universes specified for parameter: \"" + paramName + "\"");

        const auto weights = CrossSectionHelper::GetMutuallyExclusiveWeights(truth, parameters, nUniverses);

        // Add the result to the output map
        map.emplace(paramName, weights);
    }

    // Now get the weights from the rest of the parameters
    auto availableParamWeights = CrossSectionHelper::GetWeightsMap(truth, availableParamDimensions);

    // Combine this together and return it!
    map.insert(availableParamWeights.begin(), availableParamWeights.end());
    return map;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::GetMutuallyExclusiveWeights(const Event::Truth &truth, const std::vector<std::string> &parameters, const unsigned int nUniverses)
{
    // Build a dimensions map using the same number of universes for each parameter
    SystDimensionsMap dimensions;
    for (const auto &paramName : parameters)
        dimensions.emplace(paramName, nUniverses);

    // Get the weights map for these parameters individually
    const auto weightsMap = CrossSectionHelper::GetWeightsMap(truth, dimensions);

    // Extract the mutually exclusive weights
    std::vector<float> weights;
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        // Find the parameter name for which the weight is not exactly one
        float weight = 1.f;
        bool foundNonUnitWeight = false;

        for (const auto &paramName : parameters)
        {
            const auto universeWeight = weightsMap.at(paramName).at(iUni);
            if (std::abs(universeWeight - 1.f) > std::numeric_limits<float>::epsilon())
            {
                // This weight is something other than one!

                // If we find more than one weight that's not exactly one - then the parameters aren't mutually exclusive!
                if (foundNonUnitWeight)
                    throw std::logic_error("CrossSectionHelper::GetMutuallyExclusiveWeights - Input weights aren't mutually exclusive");

                // Store this weight
                weight = universeWeight;
                foundNonUnitWeight = true;
            }
        }

        // ATTN it's still posible that all universe weights are exactly one. This is okay - just use 1.f as the weight
        weights.push_back(weight);
    }

    return weights;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> CrossSectionHelper::GenerateBootstrapWeights(const unsigned int nUniverses)//, const std::string seedString="")
{
    std::vector<float> weights;
    std::poisson_distribution<int> poisson(1.f);

    // std::seed_seq seed(seedString.begin(),seedString.end());
    // std::default_random_engine generator;//(seed);
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        // weights.push_back(static_cast<float>(poisson(m_generator)));
        const auto weight = static_cast<float>(poisson(m_generator));
        // const auto weight = static_cast<float>(poisson(generator));

        weights.push_back(weight);
    }

    return weights;
}

std::vector<float> CrossSectionHelper::GenerateNormalWeights(const unsigned int nUniverses, const float sigma, const float mean, const float globalScaling)
{
    std::vector<float> weights;
    std::normal_distribution<float> normal{mean,sigma};
    for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
    {
        // weights.push_back(static_cast<float>(poisson(m_generator)));
        auto weight = static_cast<float>(normal(m_generator));
        weight = std::max(weight * globalScaling, 0.f);

        weights.push_back(weight);
    }

    return weights;
}


// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::FillSystTH1FMap(const float value, const float nominalWeight, const SystFloatMap &weights, SystTH1FMap &map)
{
    // Validate the dimensions of the input SystMaps
    CrossSectionHelper::ValidateSystMap(weights, CrossSectionHelper::GetSystMapDimensions(map));

    // Fill the histograms in each systematic universe
    for (const auto &[paramName, universeWeights] : weights)
    {
        const auto nUniverses = universeWeights.size();
        auto &histVector = map.at(paramName);

        // ATTN this should never happen because of the above ValidateSystMap check, but including here for sanity
        if (histVector.size() != nUniverses)
            throw std::logic_error("CrossSectionHelper::FillSystTH1FMap - Input weights map and SystTH1FMap have different dimensions");

        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            // ATTN we use the product of the nominal weight and the universe weight
            const auto weight = nominalWeight * universeWeights.at(iUni);

            histVector.at(iUni)->Fill(value, weight);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void CrossSectionHelper::FillSystTH2FMap(const float xValue, const float yValue, const float nominalWeight, const SystFloatMap &weights, SystTH2FMap &map)
{
    // Validate the dimensions of the input SystMaps
    CrossSectionHelper::ValidateSystMap(weights, CrossSectionHelper::GetSystMapDimensions(map));

    // Fill the histograms in each systematic universe
    for (const auto &[paramName, universeWeights] : weights)
    {
        const auto nUniverses = universeWeights.size();
        auto &histVector = map.at(paramName);

        // ATTN this should never happen because of the above ValidateSystMap check, but including here for sanity
        if (histVector.size() != nUniverses)
            throw std::logic_error("CrossSectionHelper::FillSystTH2FMap - Input weights map and SystTH2Map have different dimensions");

        for (unsigned int iUni = 0; iUni < nUniverses; ++iUni)
        {
            // ATTN we use the product of the nominal weight and the universe weight
            const auto weight = nominalWeight * universeWeights.at(iUni);

            histVector.at(iUni)->Fill(xValue, yValue, weight);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::GetMatrixFromHist(const std::shared_ptr<TH1F> &pHist)
{
    std::vector<float> elements;

    const unsigned int nBins = pHist->GetNbinsX();
    for (unsigned int iBin = 1; iBin <= nBins; ++iBin)
    {
        elements.push_back(pHist->GetBinContent(iBin));
    }

    return ubsmear::UBMatrix(elements, nBins, 1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::GetMatrixFromHist(const std::shared_ptr<TH2F> &pHist)
{
    std::vector< std::vector<float> > elements;

    const unsigned int nBinsX = pHist->GetNbinsX();
    const unsigned int nBinsY = pHist->GetNbinsY();

    for (unsigned int iBinX = 1; iBinX <= nBinsX; ++iBinX)
    {
        // Make a new row
        elements.emplace_back();
        auto &row = elements.back();

        for (unsigned int iBinY = 1; iBinY <= nBinsY; ++iBinY)
        {
            row.push_back(pHist->GetBinContent(iBinX, iBinY));
        }
    }

    return ubsmear::UBMatrix(elements);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<ubsmear::UBMatrix> CrossSectionHelper::FlattenMatrix(const std::shared_ptr<ubsmear::UBMatrix> &pMatrix)
{
    if (!pMatrix)
    {
        throw std::invalid_argument("CrossSectionHelper::FlattenMatrix - Input matrix pointer is null");
    }

    auto flattened0 = ubsmear::UBSmearingHelper::Flatten(*pMatrix);
    auto flattened = std::make_shared<ubsmear::UBMatrix>(flattened0);
    return flattened;
    // return std::make_shared<ubsmear::UBMatrix>(ubsmear::UBSmearingHelper::Flatten(*pMatrix));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::tuple< std::vector<float>, bool, bool > CrossSectionHelper::GetExtendedBinEdges(const float min, const float max, const std::vector<float> &binEdges)
{
    // Ensure that we have at least one input bin
    if (binEdges.size() < 2)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Fewer than two input bin edges were supplied");

    // Ensure that the input bin edges are sorted
    if (!std::is_sorted(binEdges.begin(), binEdges.end()))
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Input bin edges are not sorted into ascending order");

    // Get the first and last bin edges supplied
    const auto firstEdge = binEdges.front();
    const auto lastEdge = binEdges.back();

    // Check if the first/last edge is already at the min/max
    const auto isFirstEdgeAtMin = (std::abs(firstEdge - min) <= std::numeric_limits<float>::epsilon());
    const auto isLastEdgeAtMax = (std::abs(lastEdge - max) <= std::numeric_limits<float>::epsilon());

    // If the first/last edge isn't at min/max, then make sure that it's not below/above
    if (!isFirstEdgeAtMin && firstEdge < min)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Lowest input bin edge is smaller than the supplied minimum value");

    if (!isLastEdgeAtMax && lastEdge > max)
        throw std::invalid_argument("CrossSectionHelper::GetExtendedBinEdges - Uppermost input bin edge is larget than the supplied maximum value");

    // Determine if we need to extend to an extra underflow/overflow bin edge
    const auto hasUnderflow = !isFirstEdgeAtMin;
    const auto hasOverflow = !isLastEdgeAtMax;

    // Setup the output bin edges
    std::vector<float> extendedBinEdges;

    // Add the underflow bin if required
    if (hasUnderflow)
        extendedBinEdges.push_back(min);

    // Add the supplied bins
    extendedBinEdges.insert(extendedBinEdges.end(), binEdges.begin(), binEdges.end());

    // Add the overflow bin if required
    if (hasOverflow)
        extendedBinEdges.push_back(max);

    return {extendedBinEdges, hasUnderflow, hasOverflow};
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> CrossSectionHelper::GetNominalFluxHistNames(const std::vector<int> &nuPdgsSignal, const std::map<int, std::string> &nuPdgToHistName, const std::string &nomHistPattern)
{
    std::vector<std::string> histNames;

    for (const auto &nuPdg : nuPdgsSignal)
    {
        const auto iter = nuPdgToHistName.find(nuPdg);
        if (iter == nuPdgToHistName.end())
            throw std::logic_error("CrossSectionHelper::GetNominalFluxHistNames - Can't find name corresponding to PDG code " + std::to_string(nuPdg));

        const auto &nuName = iter->second;
        const auto histName = std::regex_replace(nomHistPattern, std::regex("NEUTRINO"), nuName);
        histNames.push_back(histName);
    }

    return histNames;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::pair< std::vector<float>, std::vector<float> > CrossSectionHelper::ReadNominalFlux(const std::string &fileName, const std::vector<std::string> &histNames, const float pot)
{
    if (histNames.empty())
        throw std::invalid_argument("CrossSectionHelper::ReadNominalFlux - No histogram names supplied");

    // Open the flux file for reading
    TFile *pFluxFile = new TFile(fileName.c_str(), "READ");
    if (!pFluxFile->IsOpen())
        throw std::invalid_argument("CrossSectionHelper::ReadNominalFlux - Can't open flux file: " + fileName);

    // Get the scaling factor to go from the event rate in the flux file, to the flux itself
    // Here we get the flux in neutrinos/POT/bin/cm2 by scaling the event rate (in the samples) down by POT and the cross-sectional area of the active volume
    // Here we also scale up the fluxes by 1e10 for the sake of comparison just so we are working with reasonable numbers
    const float unitsScaling = 1e10;
    const auto fluxScaleFactor = unitsScaling / (pot * (GeometryHelper::highX - GeometryHelper::lowX) * (GeometryHelper::highY - GeometryHelper::lowY));
    bool isFirstHist = true;
    std::vector<float> fluxBinEdges, fluxBinValuesNominal;

    for (const auto &histName : histNames)
    {
        // Get the nominal flux
        const auto pFluxHist = static_cast<TH1F *>(pFluxFile->Get(histName.c_str()));
        if (!pFluxHist)
            throw std::logic_error("CrossSectionHelper::ReadNominalFlux - Input file doesn't contain histogram with name: " + histName);

        // Get the number of bins
        const unsigned int nFluxBins = pFluxHist->GetNbinsX();
        if (!isFirstHist && nFluxBins != (fluxBinEdges.size() - 1))
            throw std::logic_error("CrossSectionHelper::ReadNominalFlux - Supplied flux histograms have an inconsistent number of bins");

        // Get the flux bin edges and content
        for (unsigned int iBin = 1; iBin <= nFluxBins; ++iBin)
        {
            const auto flux = pFluxHist->GetBinContent(iBin) * fluxScaleFactor;
            const float lowEdge = pFluxHist->GetBinLowEdge(iBin);
            const float binWidth = pFluxHist->GetBinWidth(iBin);

            if (isFirstHist)
            {
                fluxBinEdges.push_back(lowEdge);

                // Add the upper edge of the last bin
                if (iBin == nFluxBins)
                {
                    fluxBinEdges.push_back(lowEdge + binWidth);
                }

                fluxBinValuesNominal.push_back(flux);
            }
            else
            {
                // If this isn't the first flux histogram, then add to the existing flux
                fluxBinValuesNominal.at(iBin - 1) += flux;

                // Check that the bin edges are consistent
                if (std::abs(fluxBinEdges.at(iBin - 1) - lowEdge) > std::numeric_limits<float>::epsilon())
                    throw std::logic_error("CrossSectionHelper::ReadNominalFlux - Supplied flux histograms have an inconsistent binning (low edge)");
            }
        }

        isFirstHist = false;
    }

    return {fluxBinEdges, fluxBinValuesNominal};
}

// -----------------------------------------------------------------------------------------------------------------------------------------

ubsmear::UBMatrix CrossSectionHelper::GetErrorMatrix(const ubsmear::UBMatrix &biasVector, const ubsmear::UBMatrix &covarianceMatrix)
{
    // Make sure the dimensions of the input are valid
    if (biasVector.GetColumns() != 1)
    {
        std::cout<<"CrossSectionHelper::GetErrorMatrix - Input bias vector isn't a column vector"<<std::endl;
        throw std::invalid_argument("CrossSectionHelper::GetErrorMatrix - Input bias vector isn't a column vector");
    }

    if (!ubsmear::UBMatrixHelper::IsSquare(covarianceMatrix))
    {
        std::cout<<"CrossSectionHelper::GetErrorMatrix - Input covariance matrix isn't square"<<std::endl;
        throw std::invalid_argument("CrossSectionHelper::GetErrorMatrix - Input covariance matrix isn't square");
    }

    const auto nBins = covarianceMatrix.GetRows();
    if (biasVector.GetRows() != nBins)
    {
        std::cout<<"CrossSectionHelper::GetErrorMatrix - Input bias vector and covariance matrix have incompatible dimensions"<<std::endl;
        throw std::invalid_argument("CrossSectionHelper::GetErrorMatrix - Input bias vector and covariance matrix have incompatible dimensions");
    }

    // Make a new matrix for the output
    auto errorMatrix = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
    for (unsigned int iRow = 0; iRow < nBins; ++iRow)
    {
        for (unsigned int iCol = 0; iCol < nBins; ++iCol)
        {
            errorMatrix.SetElement(iRow, iCol, (biasVector.At(iRow, 0) * biasVector.At(iCol, 0)) + covarianceMatrix.At(iRow, iCol));
        }
    }

    return errorMatrix;
}

} // namespace ubcc1pi
