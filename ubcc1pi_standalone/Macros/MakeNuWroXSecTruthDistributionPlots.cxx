/**
 *  @file  ubcc1pi_standalone/Macros/MakeNuWroXSecTruthDistributionPlots.cxx
 *
 *  @brief The implementation file of the MakeNuWroXSecTruthDistributionPlots macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"

#include "ubsmear.h"

#include <TStyle.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeNuWroXSecTruthDistributionPlots(const Config &config)
{

    if(!config.global.useEfficiencyCorrection)
    {
        std::cout << "############################" << std::endl;
        std::cout << "WARNING: Efficiency effects are being ignored in the forward-folding matrix.\nThis is NOT correct and should only be used for illustrative plots showing efficiency effects!" << std::endl;
        std::cout << "############################" << std::endl;
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get the binning of each cross-section
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::map<std::string, ubsmear::UBXSecMeta> metadataMap;

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
        std::cout << "Getting bins for: " << name << std::endl;

        // Get the bin edges from the input configuration
        const auto &[extendedBinEdges, hasUnderflow, hasOverflow] = CrossSectionHelper::GetExtendedBinEdges(binning.min, binning.max, binning.binEdges);
        metadataMap.emplace(name, ubsmear::UBXSecMeta(extendedBinEdges, hasUnderflow, hasOverflow, scaleByBinWidth));
    }

    // Add the dummy metadata for the total cross-section (see ExtractXSecs for more details)
    metadataMap.emplace("total", ubsmear::UBXSecMeta({-1.f, 1.f}, false, false, false));

    const std::map<std::string, CrossSectionHelper::SystDimensionsMap> systDimensionMapBNB(
        {
            {"flux", config.extractXSecs.fluxDimensions},
            {"xsec", config.extractXSecs.xsecDimensions},
            {"reint", config.extractXSecs.reintDimensions},
            {"misc", {
                {"bootstrap", config.extractXSecs.nBootstrapUniverses},
                {"sidebandWeights", config.extractXSecs.nBootstrapUniverses},
                {"dirt", 2},
                {"POT", 0}
            }}
        });

    const std::map<std::string, CrossSectionHelper::SystDimensionsMap> systDimensionMapFakeData(
        {
            {"xsec", config.extractXSecs.xsecDimensions},
            {"misc", {
                {"bootstrap", config.extractXSecs.nBootstrapUniverses},
                {"sidebandWeights", config.extractXSecs.nBootstrapUniverses},
            }}
        });

    std::vector<std::string> dataTypeList;
    if(config.global.useNuWroAsData) dataTypeList.push_back("NuWro");
    if(config.global.useBNBAsData) dataTypeList.push_back("BNB");
    if(config.global.useGenieAsData) dataTypeList.push_back("Genie");

    // String that is used to label plots and tables when efficiency correction is ignored
    const auto infix = std::string(config.global.useEfficiencyCorrection ? "" : "_notEfficiencyCorrected");

    for(const auto dataTypeName: dataTypeList)
    {
        // -------------------------------------------------------------------------------------------------------------------------------------
        // Setup a table to hold the chi2 values between the data and prediction for each cross-section
        // -------------------------------------------------------------------------------------------------------------------------------------
        FormattingHelper::Table tableScaled({"Selection", "Cross-section", "nBins", "", "Chi2", "DoF", "Chi2/DoF", "p-value"});
        FormattingHelper::Table tableUnscaled({"Selection", "Cross-section", "nBins", "", "Chi2", "DoF", "Chi2/DoF", "p-value"});

        // -------------------------------------------------------------------------------------------------------------------------------------
        // Loop over all possible cross-sections
        // -------------------------------------------------------------------------------------------------------------------------------------
        for (const auto &[selectionName, enabledMap] : config.extractXSecs.crossSectionIsEnabled)
        {
            for (const auto &[xsecName, isEnabled] : enabledMap)
            {
                std::cout<<"DEBug Point 0"<<std::endl;
                // Skip cross-sections that aren't enabled
                if (!isEnabled)
                    continue;

                std::cout << selectionName << " - " << xsecName << std::endl;

                // Define a lambda function to get a matrix from a file for brevity
                const auto getMatrix = [&, selectionName = selectionName, xsecName = xsecName](const std::string &identifier) -> ubsmear::UBMatrix {
                    return ubsmear::UBFileHelper::ReadMatrix("xsec_" + selectionName + "_" + xsecName + "_" + identifier + ".txt");
                };

                // Define a lambda function to get a matrix from a file an trim any overflow / underflow bins
                const auto &metadata = metadataMap.at(xsecName);
                const auto getTrimmedMatrix = [&] (const std::string &identifier) -> ubsmear::UBMatrix {
                    return ubsmear::UBSmearingHelper::TrimUnderOverflowBins(getMatrix(identifier), metadata);
                };

                const auto extendedBinEdges = metadata.GetBinEdges();
                std::vector<float> binEdges;
                for (unsigned int iBin = 0; iBin < metadata.GetNBins(); ++iBin)
                {
                    // Skip underflow/overflow bins
                    if (metadata.IsUnderOverflowBin(iBin))
                        continue;

                    // If this is the first bin then add the lower edge
                    if (binEdges.empty())
                    {
                        if (metadata.IsScaledByBinWidth())
                        {
                            binEdges.push_back(extendedBinEdges.at(iBin));
                        }
                        else
                        {
                            // If we don't scale by bin width, then just use zero as the first bin edge
                            binEdges.push_back(0.f);
                        }
                    }

                    // Add the upper bin edge
                    if (metadata.IsScaledByBinWidth())
                    {
                        binEdges.push_back(extendedBinEdges.at(iBin + 1));
                    }
                    else
                    {
                        // If we don't scale by bin width, then just use unit width bins
                        binEdges.push_back(binEdges.back() + 1.f);
                    }
                }

                const std::string prefix = "xsecPlots_" + selectionName + "_" + xsecName;
                auto pPredictionHistScaled = std::make_shared<TH1F>((prefix + "_prediction_"+dataTypeName+"Scaled").c_str(), "", binEdges.size() - 1, binEdges.data());

                // Now get the predicted cross-section and it's error matrix
                const auto predictionScaled = getMatrix("prediction_"+dataTypeName+"Scaled");
                auto pNuWroPredictionHist = std::make_shared<TH1F>((prefix + "_predictionNuWro").c_str(), "", binEdges.size() - 1, binEdges.data());

                // Fill the bins
                // ATTN here we only show the diagonals of the error matrices
                float minY = +std::numeric_limits<float>::max();
                float maxY = -std::numeric_limits<float>::max();
                for (unsigned int iBin = 1; iBin <= binEdges.size() - 1; ++iBin)
                {
                    // Set the values of the prediction
                    const auto predictionValueScaled = predictionScaled.At(iBin - 1, 0);

                    pPredictionHistScaled->SetBinContent(iBin, predictionValueScaled);

                    // For the proton multiplicity plot, use explicit bin labels
                    if (xsecName == "nProtons")
                    {
                        for (auto &pHist : {pPredictionHistScaled})
                        {
                            pHist->GetXaxis()->SetBinLabel(1, "0");
                            pHist->GetXaxis()->SetBinLabel(2, "1");
                            pHist->GetXaxis()->SetBinLabel(3, ">1");
                        }
                    }

                    // Get the limiting values
                    minY = std::min(minY, predictionValueScaled);
                    maxY = std::max(maxY, predictionValueScaled);
                }

                // Set the y-range
                std::cout << "minY = " << minY << ", maxY = " << maxY << std::endl;
                const auto padding = maxY * 0.1;
                maxY += padding;
                minY -= padding;
                minY = std::max(minY, 0.f);
                minY = 0.f; // Remove this line to get a dynamic lower y-range


                pPredictionHistScaled->GetYaxis()->SetRangeUser(minY, maxY);

                // Set the colours of the histograms
                PlottingHelper::SetLineStyle(pPredictionHistScaled, PlottingHelper::Secondary);

                // Make the plot!
                auto pCanvas = PlottingHelper::GetCanvas();

                gStyle->SetEndErrorSize(4);

                // Draw the smeared prediction
                pPredictionHistScaled->Draw("hist");

                PlottingHelper::SaveCanvas(pCanvas, prefix + "_data-vs-smearedPrediction-truthDistribution_" + dataTypeName);
            }
        }

        // // Save the table
        // tableScaled.WriteToFile("XSecPlots_goodnessOfFitStatistics_"+dataTypeName+"Scaled"+ infix + ".md");
        // tableUnscaled.WriteToFile("XSecPlots_goodnessOfFitStatistics_"+dataTypeName+"Unscaled"+ infix + ".md");
    }
}

} // namespace ubcc1pi_macros
