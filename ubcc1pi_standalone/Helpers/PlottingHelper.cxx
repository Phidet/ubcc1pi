/**
 *  @file  ubcc1pi_standalone/Helpers/PlottingHelper.cxx
 *
 *  @brief The implementation file of the plotting helper class
 */

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include <stdexcept>

#include <THStack.h>
#include <TLine.h>
#include <TStyle.h>

namespace ubcc1pi
{

PlottingHelper::MultiPlot::MultiPlot(const std::string &xLabel, const std::string &yLabel, unsigned int nBins, float min, float max, bool drawErrors, const bool useAxisTitles) :
    MultiPlot(xLabel, yLabel, PlottingHelper::GenerateUniformBinEdges(nBins, min, max), drawErrors, useAxisTitles)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::MultiPlot::MultiPlot(const std::string &xLabel, const std::string &yLabel, const std::vector<float> &binEdges, bool drawErrors, const bool useAxisTitles) :
    /// @cond Doxygen can't handle this initilizer list
    m_xLabel(xLabel),
    m_nBins(binEdges.size() - 1),
    m_min(binEdges.front()),
    m_max(binEdges.back()),
    m_binEdges(binEdges),
    m_id(++m_lastId),
    m_cloneCount(0),
    m_drawErrors(drawErrors)
    /// @endcond
{
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Make a unique name for this plot to avoid collisionsi
        const auto nameStr = "ubcc1pi_multiPlot_" + std::to_string(m_id) + "_" + std::to_string(static_cast<int>(style));
        const auto name = nameStr.c_str();

        auto pHist = std::make_shared<TH1F>(name, "", m_nBins, m_binEdges.data());

        pHist->Sumw2();
        PlottingHelper::SetLineStyle(pHist.get(), style);
        if (useAxisTitles){
            pHist->GetXaxis()->SetTitle(xLabel.c_str());
            pHist->GetYaxis()->SetTitle(yLabel.c_str());
        }
        m_plotToHistMap.emplace(style, pHist);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SetIntegerBinLabels()
{
    // Check that the bin edges lie on an integer and get those integers
    std::vector<int> intBinEdges;
    for (const auto &edge : m_binEdges)
    {
        if (std::abs(edge - std::round(edge)) > std::numeric_limits<float>::epsilon())
            throw std::logic_error("MultiPlot::SetIntegerBinLabels - found bin edge of: " + std::to_string(edge) + ", which isn't an integer");

        intBinEdges.push_back(static_cast<int>(std::round(edge)));
    }

    // Get the bin labels
    std::vector<std::string> labels;
    for (unsigned int iBin = 0; iBin < m_nBins; ++iBin)
    {
        const auto lower = intBinEdges.at(iBin);
        const auto upper = intBinEdges.at(iBin + 1);
        const auto width = upper - lower;

        // If this bin contains one integer, then just label it with that value
        if (width == 1)
        {
            labels.push_back(std::to_string(lower));
            continue;
        }

        // If the bin contained multiple integers then use the range
        labels.push_back(std::to_string(lower) + " -> " + std::to_string(upper - 1));
    }

    // Set the bin labels
    this->SetBinLabels(labels);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SetBinLabels(const std::vector<std::string> &labels)
{
    if (labels.size() != m_nBins)
        throw std::invalid_argument("MultiPlot::SetBinLabels - there are " + std::to_string(m_nBins) + " bins, but only " + std::to_string(labels.size()) + " labels were provided");

    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto &pHist = m_plotToHistMap.at(style);
        auto pXAxis = pHist->GetXaxis();

        for (unsigned int iBin = 0; iBin < m_nBins; ++iBin)
        {
            // ATTN root indexes it's bins from 1
            pXAxis->SetBinLabel(iBin + 1, labels.at(iBin).c_str());
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::AddCutLine(const float value)
{
    if (value < m_min || value > m_max)
        throw std::invalid_argument("MultiPlot::AddCutLine - supplied value: " + std::to_string(value) + " is out of range: " + std::to_string(m_min) + " -> " + std::to_string(m_max));

    m_cutValues.push_back(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::Fill(const float value, const PlotStyle &plotStyle, const float weight)
{
    auto iter = m_plotToHistMap.find(plotStyle);
    if (iter == m_plotToHistMap.end())
        throw std::invalid_argument("Input plot style is unknown");

    if (value < m_min || value > m_max)
        return;

    // Fill the histogram
    iter->second->Fill(value, weight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::GetHistogramClones(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap)
{
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        const auto pHist = m_plotToHistMap.at(style);
        const auto pHistClone = static_cast<TH1F *>(pHist->Clone());
        const auto cloneNameStr = std::string(pHist->GetName()) + "_clone_" + std::to_string(m_cloneCount);
        const auto cloneName = cloneNameStr.c_str();
        pHistClone->SetName(cloneName);

        pHistClone->GetXaxis()->SetTitle(pHist->GetXaxis()->GetTitle());
        pHistClone->GetYaxis()->SetTitle(pHist->GetYaxis()->GetTitle());

        plotToHistCloneMap.emplace(style, pHistClone);
    }

    m_cloneCount++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::ScaleHistograms(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap, const bool scaleByBinWidth) const
{
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto pHist = plotToHistCloneMap.at(style);
        const auto option = scaleByBinWidth ? "width" : "";
        pHist->Scale(1.f, option);
        pHist->Scale(1.f / pHist->Integral());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SetHistogramYRanges(const bool useLogY, const unsigned int minEntriesToDraw, std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap) const
{
    float yMin = std::numeric_limits<float>::max();
    float yMax = -std::numeric_limits<float>::max();

    bool shouldScale = false;
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        const auto pHist = plotToHistCloneMap.at(style);

        if (pHist->GetEntries() < minEntriesToDraw || pHist->GetEntries() == 0u)
            continue;

        yMax = std::max(yMax, static_cast<float>(pHist->GetMaximum()));
        shouldScale = true;

        for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pHist->GetNbinsX()); ++iBin)
        {
            const float binContent = pHist->GetBinContent(iBin);
            if (useLogY && binContent < std::numeric_limits<float>::epsilon())
                continue;

            yMin = std::min(yMin, binContent);
        }
    }

    if (!shouldScale)
        return;

    // Add some padding to the top of the histogram
    yMax += (yMax - yMin) * 0.05;

    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto pHist = plotToHistCloneMap.at(style);
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SaveAs(const std::string &fileName, const bool useLogX, const bool scaleByBinWidth, const unsigned int minEntriesToDraw, const bool useLogY)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    if (useLogX)
    {
        pCanvas->SetLogx();
    }

    if (useLogY)
    {
        pCanvas->SetLogy();
    }

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<PlotStyle, TH1F* > plotToHistCloneMap;
    this->GetHistogramClones(plotToHistCloneMap);
    this->ScaleHistograms(plotToHistCloneMap, scaleByBinWidth);
    this->SetHistogramYRanges(useLogY, minEntriesToDraw, plotToHistCloneMap);

    // Draw the error bands if required
    bool isFirst = true;
    if (m_drawErrors)
    {
        for (const auto &style : PlottingHelper::AllPlotStyles)
        {
            // The points already have error bars
            if (PlottingHelper::ShouldUsePoints(style))
                continue;

            // Don't draw BNB data
            if (style == BNBData)
                continue;

            auto pHist = plotToHistCloneMap.at(style);

            if (pHist->GetEntries() < minEntriesToDraw || pHist->GetEntries() == 0u)
                continue;

            // Draw a clone of the histogram so we can safely change it's style
            auto pHistClone = static_cast<TH1F *>(pHist->Clone());
            const auto col = pHistClone->GetLineColor();
            pHistClone->SetFillStyle(1001);
            pHistClone->SetLineColorAlpha(col, 0.f);
            pHistClone->SetFillColorAlpha(col, 0.3f);

            pHistClone->Draw(isFirst ? "e2" : "e2 same");
            isFirst = false;
        }
    }

    // Draw the histogram itself
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Don't draw BNB data
        if (style == BNBData)
            continue;

        auto pHist = plotToHistCloneMap.at(style);

        if (pHist->GetEntries() < minEntriesToDraw || pHist->GetEntries() == 0u)
            continue;


        const bool usePoints = PlottingHelper::ShouldUsePoints(style);
        const TString tstyle = usePoints ? "e1" : "hist";
        pHist->Draw(isFirst ? tstyle : tstyle + " same");
        isFirst = false;
    }

    PlottingHelper::SaveCanvas(pCanvas, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SaveAsStacked(const std::string &fileName, const bool useLogX, const bool scaleByBinWidth, const bool useLogY, const bool useAxisTitles)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    if (useLogX)
    {
        pCanvas->SetLogx();
    }

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<PlotStyle, TH1F* > plotToHistCloneMap;
    this->GetHistogramClones(plotToHistCloneMap);

    // Scale by bin width
    if (scaleByBinWidth)
    {
        for (auto &entry : plotToHistCloneMap)
            entry.second->Scale(1.f, "width");
    }

    // Work out if we have BNB data to plot, and if we do then get the minimum and maximum Y coordinates
    auto bnbDataHistIter = plotToHistCloneMap.find(BNBData);
    const bool hasBNBData = bnbDataHistIter != plotToHistCloneMap.end();

    auto yMax = hasBNBData ? static_cast<float>(bnbDataHistIter->second->GetMaximum()) : -std::numeric_limits<float>::max();

    // For log-y plots we also want to know the smallest non-zero bin entry
    float yMinNonZero = std::numeric_limits<float>::max();

    // Sum the non BNB data histograms to get the "MC" total
    const auto nameTotalStr = "ubcc1pi_multiPlot_" + std::to_string(m_id) + "_total";
    const auto nameTotal = nameTotalStr.c_str();
    auto pHistTotal = std::make_shared<TH1F>(nameTotal, "", m_nBins, m_binEdges.data());
    pHistTotal->Sumw2();

    bool isFirst = true;
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Don't add BNB data to the total, this will be drawn seperately
        if (style == BNBData)
            continue;

        auto pHist = plotToHistCloneMap.at(style);

        if (pHist->GetEntries() == 0)
            continue;

        if (isFirst)
        {
            // If we are using alphanumeric bin labels, then we need to copy those labels to the total histogram or root will complain
            const auto pAxis = pHist->GetXaxis();
            if (pAxis->IsAlphanumeric())
            {
                for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pAxis->GetNbins()); ++iBin)
                    pHistTotal->GetXaxis()->SetBinLabel(iBin, pAxis->GetBinLabel(iBin));
            }

            isFirst = false;

            if (useAxisTitles){
                pHistTotal->GetXaxis()->SetTitle(pHist->GetXaxis()->GetTitle());
                pHistTotal->GetYaxis()->SetTitle(pHist->GetYaxis()->GetTitle());
            }
        }

        pHistTotal->Add(pHist);
        yMax = std::max(yMax, static_cast<float>(pHist->GetMaximum()));

        // For log-y plots we also want to know the smallest non-zero bin entry
        for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pHist->GetNbinsX()); ++iBin)
        {
            const float content = pHistTotal->GetBinContent(iBin);
            if (content <= std::numeric_limits<float>::epsilon())
                continue;

            yMinNonZero = std::min(yMinNonZero, content);
        }
    }

    // Get the maximum and minimum Y coordinates
    yMax = std::max(yMax, static_cast<float>(pHistTotal->GetMaximum()));
    yMax *= 1.05;

    // Draw the stacked histogram
    const auto nameStackStr = "ubcc1pi_plotPlot_" + std::to_string(m_id) + "_stack";
    const auto nameStack = nameStackStr.c_str();
    std::string histTitle = string(";")+pHistTotal->GetXaxis()->GetTitle()+string(";")+pHistTotal->GetYaxis()->GetTitle();
    if (!useAxisTitles) histTitle = std::string("");
    auto pHistStack = std::make_shared<THStack>(nameStack, histTitle.c_str());

    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Don't draw BNB data
        if (style == BNBData)
            continue;

        auto pHist = plotToHistCloneMap.at(style);

        if (pHist->GetEntries() == 0)
            continue;

        const auto col = pHist->GetLineColor();
        pHist->SetLineColorAlpha(col, 0.f);
        pHist->SetFillStyle(1001);
        pHist->SetFillColor(col);
        pHist->SetMarkerSize(0);
        pHist->SetMarkerStyle(0);


        pHistStack->Add(pHist);
    }


    if (useLogY)
    {
        // Add some padding at the bottom
        pHistStack->SetMinimum(yMinNonZero * 0.9);
        pCanvas->SetLogy();
    }
    else
    {
        pHistStack->SetMinimum(0.f);
    }

    pHistStack->SetMaximum(yMax / (1+gStyle->GetHistTopMargin()));
    pHistStack->Draw("hist");

    // Draw the error bands on the stack if required
    if (m_drawErrors)
    {
        pHistTotal->SetFillStyle(1001);
        pHistTotal->SetLineColorAlpha(kBlack, 0.f);
        pHistTotal->SetFillColorAlpha(kBlack, 0.3f);
        pHistTotal->SetMarkerSize(0);
        pHistTotal->SetMarkerStyle(0);
        pHistTotal->Draw("e2 same");
    }

    // Draw the BNB data
    if (hasBNBData)
    {
        bnbDataHistIter->second->Draw("e1 same");
    }

    // Draw the cut lines
    for (const auto &cutValue : m_cutValues)
    {
        TLine *pLine = new TLine(cutValue, 0.f, cutValue, yMax);
        pLine->SetLineWidth(2);
        pLine->Draw();
    }

    PlottingHelper::SaveCanvas(pCanvas, fileName);

    // If we don't have BNB data then we are done!
    if (!hasBNBData)
        return;

    // Now make the ratio plot
    const auto ratioPlotHeightDefault = 270u;
    auto pCanvasRatio = PlottingHelper::GetCanvas(960, 270);

    if (useLogX)
    {
        // ATTN the 10^N labels for a log-x plot are shifted downwards so we need a bit more space to fit them in
        // Get the length of the top and bottom margins
        const float topMarginLength = pCanvasRatio->GetTopMargin() * ratioPlotHeightDefault;
        const float bottomMarginLength = pCanvasRatio->GetBottomMargin() * ratioPlotHeightDefault;
        const auto remainingLength = ratioPlotHeightDefault - topMarginLength - bottomMarginLength;

        // Set the new desired bottom margin length
        const float newBottomMarginLength = bottomMarginLength * 1.5f;

        // Work out how tall we need to make the canvas to account for this
        const unsigned int newRatioPlotHeight = std::round(topMarginLength + remainingLength + newBottomMarginLength);
        const float newBottomMarginFrac = newBottomMarginLength / static_cast<float>(newRatioPlotHeight);

        // Setup a new canvas with this height
        pCanvasRatio = PlottingHelper::GetCanvas(960, newRatioPlotHeight);
        pCanvasRatio->SetBottomMargin(newBottomMarginFrac);
        pCanvasRatio->SetLogx();
    }

    auto pHistBNBData = static_cast<TH1F *>(bnbDataHistIter->second->Clone());
    if (useAxisTitles) pHistBNBData->GetYaxis()->SetTitle("Beam on / (Overlay + Beam off)");
    pHistBNBData->Divide(pHistTotal.get());

    // Set the y-axis range
    float minRatio = 0.75f;
    float maxRatio = 1.25f;
    for (unsigned int i = 1; i <= static_cast<unsigned int>(pHistBNBData->GetNbinsX()); ++i)
    {
        minRatio = std::min(minRatio, static_cast<float>(pHistBNBData->GetBinContent(i) - pHistBNBData->GetBinError(i)));
        maxRatio = std::max(maxRatio, static_cast<float>(pHistBNBData->GetBinContent(i) + pHistBNBData->GetBinError(i)));
    }
    const auto padding = (maxRatio - minRatio) * 0.05f;
    const auto ratioYMin = std::max(0.f, minRatio - padding);
    const auto ratioYMax = maxRatio + padding;
    pHistBNBData->GetYaxis()->SetRangeUser(ratioYMin, ratioYMax);
    pHistBNBData->GetYaxis()->SetNdivisions(7, 4, 0, true);

    // Draw
    pHistBNBData->Draw("e1");
    pCanvasRatio->Update();

    // Add the lines at 0.8, 1.0 and 1.2
    const auto ratioXMin = pHistBNBData->GetXaxis()->GetXmin();
    const auto ratioXMax = pHistBNBData->GetXaxis()->GetXmax();
    TLine *l=new TLine(ratioXMin,1.0,ratioXMax,1.0);
    TLine *lPlus=new TLine(ratioXMin,1.2,ratioXMax,1.2);
    TLine *lMinus=new TLine(ratioXMin,0.8,ratioXMax,0.8);
    l->SetLineColor(kBlue);
    lPlus->SetLineColor(kBlue);
    lMinus->SetLineColor(kBlue);
    lPlus->SetLineStyle(2);
    lMinus->SetLineStyle(2);
    l->Draw();
    lPlus->Draw();
    lMinus->Draw();

    // Draw the cut lines
    for (const auto &cutValue : m_cutValues)
    {
        TLine *pLine = new TLine(cutValue, ratioYMin, cutValue, ratioYMax);
        pLine->SetLineWidth(2);
        pLine->Draw();
    }

    PlottingHelper::SaveCanvas(pCanvasRatio, fileName + "_ratio");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::EfficiencyPlot::EfficiencyPlot(const std::string &xLabel, unsigned int nBins, float min, float max, const std::vector<string> &cuts, bool drawErrors, const bool useAxisTitles) :
    EfficiencyPlot(xLabel, PlottingHelper::GenerateUniformBinEdges(nBins, min, max), cuts, drawErrors)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::EfficiencyPlot::EfficiencyPlot(const std::string &xLabel, const std::vector<float> &binEdges, const std::vector<string> &cuts, bool drawErrors, const bool useAxisTitles) :
    /// @cond Doxygen can't handle this initilizer list
    m_xLabel(xLabel),
    m_binEdges(binEdges),
    m_nBins(binEdges.size() - 1),
    m_min(binEdges.front()),
    m_max(binEdges.back()),
    m_cuts(cuts),
    m_drawErrors(drawErrors),
    m_id(++m_lastId)
    /// @endcond
{
    if (m_cuts.empty())
        throw std::logic_error("EfficiencyPlot::EfficiencyPlot - There are no cuts set!");

    // Make histograms for the cut
    for (unsigned int iCut = 0; iCut < m_cuts.size(); ++iCut)
    {
        const auto &cut = m_cuts.at(iCut);

        if (std::count(m_cuts.begin(), m_cuts.end(), cut) != 1)
            throw std::invalid_argument("EfficiencyPlot::EfficiencyPlot - Repeated cut: \"" + cut + "\"");

        // Define the plot names for root
        const std::string nameStr = "ubcc1pi_efficiencyPlot_" + std::to_string(m_id) + "_" + cut;
        const auto nameNumerator = nameStr + "_numerator";
        const auto nameDenominator = nameStr + "_denominator";

        // Make the pair of plots, one for the numerator one for denominator
        auto plotPair = std::pair< std::shared_ptr<TH1F>, std::shared_ptr<TH1F> >(
            std::make_shared<TH1F>(nameNumerator.c_str(), "", m_nBins, m_binEdges.data()),
            std::make_shared<TH1F>(nameDenominator.c_str(), "", m_nBins, m_binEdges.data()));

        // Setup the plots
        for (auto &pHist : {plotPair.first, plotPair.second})
        {
            pHist->Sumw2();
            if (useAxisTitles){
                pHist->GetXaxis()->SetTitle(xLabel.c_str());
                pHist->GetYaxis()->SetTitle("Efficiency");
            }
        }

        m_cutToPlotsMap.emplace(cut, plotPair);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::AddCutLine(const float value)
{
    if (value < m_min || value > m_max)
        throw std::invalid_argument("EfficiencyPlot::AddCutLine - supplied value: " + std::to_string(value) + " is out of range: " + std::to_string(m_min) + " -> " + std::to_string(m_max));

    m_cutValues.push_back(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::AddEvent(const float value, const float weight, const std::string &cut, const bool passedCut)
{
    auto plotPairIter = m_cutToPlotsMap.find(cut);

    if (plotPairIter == m_cutToPlotsMap.end())
        throw std::invalid_argument("EfficiencyPlot::AddEvent - Unknown cut: \"" + cut + "\"");

    auto &pHistNumerator = plotPairIter->second.first;
    auto &pHistDenominator = plotPairIter->second.second;

    if (passedCut)
        pHistNumerator->Fill(value, weight);

    pHistDenominator->Fill(value, weight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::SaveAs(const std::vector<std::string> &cuts, const std::vector<int> &colors, const std::string &fileName, const bool useAxisTitles)
{
    if (cuts.size() != colors.size())
        throw std::invalid_argument("EfficiencyPlot::SaveAs - number of cuts and colors doesn't match");

    for (const auto &cut : cuts)
    {
        if (std::find(m_cuts.begin(), m_cuts.end(), cut) == m_cuts.end())
            throw std::invalid_argument("EfficiencyPlot::SaveAs - Unknown cut: \"" + cut + "\"");
    }

    auto pCanvas = PlottingHelper::GetCanvas();

    // Save the raw histogram - just use the denominator of the first cut
    auto pHistRaw = static_cast<TH1F *>(m_cutToPlotsMap.at(m_cuts.front()).second->Clone());
    if (useAxisTitles) pHistRaw->GetYaxis()->SetTitle("Number of events");
    PlottingHelper::SetLineStyle(pHistRaw, Default);

    if (m_drawErrors)
    {
        auto pHistClone = static_cast<TH1F *>(pHistRaw->Clone());
        const auto deafultCol = PlottingHelper::GetColor(Default);
        pHistClone->SetFillStyle(1001);
        pHistClone->SetLineColorAlpha(deafultCol, 0.f);
        pHistClone->SetFillColorAlpha(deafultCol, 0.3f);

        pHistClone->Draw("e2");
    }

    // ATTN the raw histogram only works for uniform binning!
    pHistRaw->Draw(m_drawErrors ? "hist same" : "hist");
    PlottingHelper::SaveCanvas(pCanvas, fileName + "_raw");

    // Get the efficiency histograms
    std::vector<TH1F *> efficiencyHists;
    std::vector<TH1F *> efficiencyErrorHists;
    for (unsigned int i = 0; i < cuts.size(); ++i)
    {
        const auto &cut = cuts.at(i);
        const auto col = colors.at(i);

        const auto &plotPair = m_cutToPlotsMap.at(cut);
        const auto &pHistNumerator = plotPair.first;
        const auto &pHistDenominator = plotPair.second;

        // Clone the numerator and divide by the denominator to get the efficiency
        auto pHistEfficiency = static_cast<TH1F *>(pHistNumerator->Clone());
        pHistEfficiency->Divide(pHistDenominator.get());
        PlottingHelper::SetLineStyle(pHistEfficiency, col);
        efficiencyHists.push_back(pHistEfficiency);

        if (m_drawErrors)
        {
            auto pHistClone = static_cast<TH1F *>(pHistEfficiency->Clone());
            pHistClone->SetFillStyle(1001);
            pHistClone->SetLineColorAlpha(col, 0.f);
            pHistClone->SetFillColorAlpha(col, 0.3f);

            // We now need to manually find the uncertainty because root assumes the numerator and denominator are independent Poisson
            // distributed variables - here we use the binomial distribution to account for their dependence.
            for (unsigned int iBin = 1; iBin <= static_cast<unsigned int>(pHistClone->GetNbinsX()); ++iBin)
            {
                const auto numerator = pHistNumerator->GetBinContent(iBin);
                const auto denominator = pHistDenominator->GetBinContent(iBin);
                const auto error = (denominator <= 0.f) ? 0.f : AnalysisHelper::GetEfficiencyUncertainty(numerator, denominator);
                pHistClone->SetBinError(iBin, error);
            }

            efficiencyErrorHists.push_back(pHistClone);
        }
    }

    // Get the y-range
    float yMin = std::numeric_limits<float>::max();
    float yMax = -std::numeric_limits<float>::max();

    for (const auto &pHist : efficiencyHists)
    {
        yMin = std::min(yMin, static_cast<float>(pHist->GetMinimum()));
        yMax = std::max(yMax, static_cast<float>(pHist->GetMaximum()));
    }

    // Add some padding to the top of the histogram
    yMax += (yMax - yMin) * 0.05;

    //// TEST
    yMin = 0.f;
    yMax = 1.05f;
    //// END TEST

    // Make the cut lines
    std::vector< shared_ptr<TLine> > lines;
    for (const auto &cutValue : m_cutValues)
    {
        lines.emplace_back(new TLine(cutValue, yMin, cutValue, yMax));
        auto &pLine = lines.back();
        pLine->SetLineWidth(2);
    }

    // Draw the individual efficiency plots for each cut
    for (unsigned int i = 0; i < cuts.size(); ++i)
    {
        const auto &cut = cuts.at(i);

        if (m_drawErrors)
        {
            auto pHistErr = efficiencyErrorHists.at(i);
            pHistErr->GetYaxis()->SetRangeUser(yMin, yMax);
            pHistErr->Draw("e2");
        }

        auto pHist = efficiencyHists.at(i);
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(!m_drawErrors ? "hist" : "hist same");

        // Draw the cut value lines
        for (const auto &pLine : lines)
            pLine->Draw();

        PlottingHelper::SaveCanvas(pCanvas, fileName + "_" + std::to_string(i) + "_" + cut);
    }

    // Draw the error histograms all together
    bool isFirst = true;
    for (const auto &pHist : efficiencyErrorHists)
    {
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(isFirst ? "e2" : "e2 same");
        isFirst = false;
    }

    // Draw the histogram lines all together
    for (const auto &pHist : efficiencyHists)
    {
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(isFirst ? "hist" : "hist same");
        isFirst = false;
    }

    // Draw the cut value lines
    for (const auto &pLine : lines)
        pLine->Draw();

    PlottingHelper::SaveCanvas(pCanvas, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::SaveAs(const std::vector<std::string> &cuts, const std::vector<PlotStyle> &styles, const std::string &fileName, const bool useAxisTitles)
{
    // Convert the styles to colors
    std::vector<int> colors;
    for (const auto &style : styles)
        colors.push_back(PlottingHelper::GetColor(style));

    this->SaveAs(cuts, colors, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::SaveAs(const std::string &fileName)
{
    // Get the colors - one for each plot
    const auto palette = PlottingHelper::GetColorVector();
    std::vector<int> colors;
    for (unsigned int i = 0; i < m_cuts.size(); ++i)
        colors.push_back(palette.at(i % palette.size()));

    this->SaveAs(m_cuts, colors, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::PlotStyle PlottingHelper::GetPlotStyle(const Event::Reco::Particle &particle, const AnalysisHelper::SampleType &sampleType, const std::vector<Event::Truth::Particle> &truthParticles, const bool usePoints, const bool useAbsPdg, const bool peleeMode /* = false*/)
{
    const auto externalType = usePoints ? ExternalPoints : External;
    const auto dirtType = usePoints ? DirtPoints : Dirt;

    if (sampleType == AnalysisHelper::DataBNB)
        return BNBData;

    if (sampleType == AnalysisHelper::Dirt)
        return dirtType;

    if (sampleType == AnalysisHelper::DataEXT)
        return externalType;

    if ((peleeMode && (!particle.pdgBacktracked.IsSet() || particle.pdgBacktracked() == -2147483647)) 
        || (!peleeMode && (!particle.hasMatchedMCParticle.IsSet() || !particle.hasMatchedMCParticle())))
    {
        return externalType;
    }

    if (truthParticles.empty())
        return externalType;

    // Get the true pdg code applying the visibility conditions
    auto truePdgCode = -std::numeric_limits<int>::max();
    bool isGolden = false;
    try
    {
        const auto truthParticle = AnalysisHelper::GetBestMatchedTruthParticle(particle, truthParticles);
        truePdgCode = useAbsPdg ? std::abs(truthParticle.pdgCode()) : truthParticle.pdgCode();
        isGolden = AnalysisHelper::IsGolden(truthParticle);
    }
    catch (const std::logic_error &)
    {
        return externalType;
    }

    switch (truePdgCode)
    {
        case 13:
            return usePoints ? MuonPoints : Muon;
        case 2212:
            return usePoints ? ProtonPoints : Proton;
        case 211:
            return (isGolden ? (usePoints ? GoldenPionPoints : GoldenPion) : (usePoints ? NonGoldenPionPoints : NonGoldenPion));
        case -211:
            return usePoints ? PiMinusPoints : PiMinus;
        case 11:
            return usePoints ? ElectronPoints : Electron;
        // ATTN Pi0s are included as primary particles in the truth information, even though they aren't "reconstruction targets" i.e. we
        // actually reconstruct the photons from their decay. So the backtracker will match a reconstructed particle that really represents
        // a photon (from a pi0 decay) to the pi0 itself. Here we choose to plot these reco particles as photons instead of pi0s because
        // it's a better representation of what they actually represent. This is a (good) approximation, because it's possible that the pi0
        // decays via another mode (pi0 -> 2 gamma has a branching ratio of 98.8%). In the case of a Dalitz decay (pi0 -> gamma e+ e-) we
        // might reconstruct the elecron/positron and here label it as a photon.
        case 111:
        case 22:
            return usePoints ? PhotonPoints : Photon;
        default:
            return usePoints ? OtherPoints : Other;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::PlotStyle PlottingHelper::GetPlotStyle(const AnalysisHelper::SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const bool useAbsPdg, const bool treatNuWroLikeData /* = false */)
{
    if (sampleType == AnalysisHelper::Dirt)
        return Dirt;

    if (sampleType == AnalysisHelper::DataEXT)
        return External;

    if (sampleType == AnalysisHelper::DataBNB || (treatNuWroLikeData && sampleType == AnalysisHelper::NuWro))
        return BNBData;

    if (!AnalysisHelper::IsFiducial(pEvent->truth.nuVertex()))
        return NonFiducial;

    // Count the particles by PDG code
    std::vector<int> foundPdgs;
    std::unordered_map<int, unsigned int> pdgCodeCountMap;
    AnalysisHelper::GetPdgCodeCountMap(pEvent->truth.particles, useAbsPdg, foundPdgs, pdgCodeCountMap);

    // Quick utility lambda
    const auto CountByPdg = [&](const int pdg) {
        if (std::find(foundPdgs.begin(), foundPdgs.end(), pdg) == foundPdgs.end())
            return 0u;

        return pdgCodeCountMap.at(pdg);
    };

    const auto hasMuon = CountByPdg(13) > 0;
    const auto hasElectron = CountByPdg(11) > 0;

    if (!hasMuon && hasElectron)
        return Nue;

    if (!hasMuon && !hasElectron)
        return NC;

    const auto nPiPlus = CountByPdg(211);
    const auto nPiMinus = CountByPdg(-211);
    const auto nPiZero = CountByPdg(111);

    const auto nPion = nPiPlus + nPiMinus + nPiZero;
    const auto nCharedPion = nPiPlus + nPiMinus;

    if (nPion == 0)
        return NumuCC0Pi;

    if (nPiZero == 1 && nCharedPion == 0)
        return NumuCC1PiZero;

    const auto hasGoldenPion = (AnalysisHelper::CountGoldenParticlesWithPdgCode(pEvent->truth.particles, 211, useAbsPdg) != 0);

    if (nCharedPion == 1 && nPiZero == 0)
    {
        if (hasGoldenPion)
            return NumuCC1PiChargedGolden;

        return NumuCC1PiChargedNonGolden;
    }

    return NumuCCOther;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::PlotStyle PlottingHelper::GetPlotStyle(const int interactionMode)
{
    switch (interactionMode)
    {
        case 0:
            return QE;
        case 1:
            return Res;
        case 2:
            return DIS;
        case 3:
            return Coh;
        case 10:
            return MEC;
        default:
            return OtherInteraction;
    }
}

PlottingHelper::PlotStyle PlottingHelper::GetPlotStyle2(const int interactionMode)
{
    switch (interactionMode)
    {
        case 0:
            return BNBData;
        case 1:
            return NumuCC0Pi;
        default:
            return Other;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TCanvas> PlottingHelper::GetCanvas(const unsigned int width, const unsigned int height)
{
    const auto nameStr = "canvas_" + std::to_string(++m_lastCanvasId);
    const auto name = nameStr.c_str();

    std::shared_ptr<TCanvas> pCanvas = std::make_shared<TCanvas>(name, "", width, height);
    pCanvas->cd();

    return pCanvas;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::SaveCanvas(std::shared_ptr<TCanvas> &pCanvas, const std::string &fileName)
{
    for (const auto &ext : {"png", "pdf", "C"})
    {
        pCanvas->SaveAs((fileName + "." + ext).c_str());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> PlottingHelper::GenerateUniformBinEdges(const unsigned int nBins, const float min, const float max)
{
    if (min >= max)
        throw std::invalid_argument("PlottingHelper::GenerateUniformBinEdges - The range: " + std::to_string(min) + " -> " + std::to_string(max) + ", is invalid");

    std::vector<float> outVect;

    for (unsigned int i = 0; i <= nBins; ++i)
        outVect.emplace_back(min + (max - min) * static_cast<float>(i) / static_cast<float>(nBins));

    return outVect;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<float> PlottingHelper::GenerateLogBinEdges(const unsigned int nBins, const float min, const float max)
{
    if (std::min(min, max) <= 0.f)
        throw std::invalid_argument("PlottingHelper::GenerateUniformBinEdges - The range: " + std::to_string(min) + " -> " + std::to_string(max) + ", includes region <= 0");

    // Get uniformly spaced bins in log-space
    const auto logMin = std::log10(min);
    const auto logMax = std::log10(max);
    auto logBinEdges = PlottingHelper::GenerateUniformBinEdges(nBins, logMin, logMax);

    // Transform these edges into non-log space
    std::vector<float> outVect;
    for (const auto &logBinEdge : logBinEdges)
    {
        outVect.push_back(std::pow(10, logBinEdge));
    }

    return outVect;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
/*
void PlottingHelper::SaveBiasVector(const std::shared_ptr<TH1F> &vector, const std::shared_ptr<TH1F> &crossSection, const bool hasUnderflow, const bool hasOverflow, const std::string &namePrefix)
{
    // Check we have sensible binning
    const unsigned int nBins = crossSection->GetNbinsX();

    const auto nUnderOverflowBins = (hasUnderflow ? 1u : 0u) + (hasOverflow ? 1u : 0u);
    if (nBins <= nUnderOverflowBins)
        throw std::invalid_argument("PlottingHelper::SaveBiasVector - Input cross section doesn't have any non-underflow/overflow bins!");

    const auto nAnalysisBins = nBins - nUnderOverflowBins;

    const unsigned int biasVectorBins = vector->GetNbinsX();
    if (biasVectorBins != nAnalysisBins)
    {
        throw std::invalid_argument("PlottingHelper::SaveBiasVector - Input bias vector has " + std::to_string(biasVectorBins) +
            " bins, but input cross-section has " + std::to_string(nBins) + " bins, of which " + std::to_string(nUnderOverflowBins) +
            " are under/overflow bins - this doesn't add up!");
    }

    // Make the plots for the bias & fractional bias vectors
    auto biasVector = std::make_shared<TH1F>(("bias_" + std::to_string(m_lastPlotId++)).c_str(), "", nAnalysisBins, 0, nAnalysisBins);
    auto fracBiasVector = std::make_shared<TH1F>(("bias_" + std::to_string(m_lastPlotId++)).c_str(), "", nAnalysisBins, 0, nAnalysisBins);

    float maxBias = -std::numeric_limits<float>::max();
    float maxFracBias = -std::numeric_limits<float>::max();

    for (unsigned int iBin = 1u; iBin <= nAnalysisBins; ++iBin)
    {
        // Set the bin labels of each plot (here we enumerate from zero)
        for (auto &pHist : {biasVector, fracBiasVector})
            pHist->GetXaxis()->SetBinLabel(iBin, std::to_string(iBin - 1).c_str());

        // Get the cross section in this bin
        const auto xSec = crossSection->GetBinContent(iBin + (hasUnderflow ? 1u : 0u));
        if (xSec <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("PlottingHelper::SaveBiasVector - Cross section in analysis bin: " + std::to_string(iBin) + " is <= 0");

        // Get the bias in this bin
        const float bias = vector->GetBinContent(iBin);
        biasVector->SetBinContent(iBin, bias);

        // Get the fractional bias
        const float fracBias = bias / xSec;
        fracBiasVector->SetBinContent(iBin, fracBias);

        // Store the maximum bias
        const auto padding = 1.05f;
        maxBias = std::max(maxBias, std::abs(bias) * padding);
        maxFracBias = std::max(maxFracBias, std::abs(fracBias) * padding);
    }

    // Make the plots
    gStyle->SetHistMinimumZero(true);
    auto pCanvas = PlottingHelper::GetCanvas();

    // Set the y-limits
    biasVector->GetYaxis()->SetRangeUser(-maxBias, +maxBias);
    fracBiasVector->GetYaxis()->SetRangeUser(-maxFracBias, +maxFracBias);


    biasVector->SetFillColor(PlottingHelper::GetColor(Default));
    fracBiasVector->SetFillColor(PlottingHelper::GetColor(Default));

    biasVector->Draw("bar");
    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_bias");

    fracBiasVector->Draw("bar");
    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_fracBias");
    gStyle->SetHistMinimumZero(false);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::SaveSmearingMatrixBiasVector(const std::shared_ptr<TH1F> &vector, const bool hasUnderflow, const bool hasOverflow, const std::string &namePrefix)
{
    // Get the number of bins of the original cross-section (this smearing matrix should have the square of this number of bins)
    const unsigned int nBins = vector->GetNbinsX();
    const auto nBinsSqrt = static_cast<unsigned int>(std::round(std::pow(static_cast<float>(nBins), 0.5f)));
    if (nBinsSqrt*nBinsSqrt != nBins)
        throw std::invalid_argument("PlottingHelper::SaveSmearingMatrixBiasVector - Number of bins for input bias vector isn't a perfect square");

    const auto nUnderOverflowBins = (hasUnderflow ? 1u : 0u) + (hasOverflow ? 1u : 0u);
    if (nBinsSqrt <= nUnderOverflowBins)
        throw std::invalid_argument("PlottingHelper::SaveSmearingMatrixBiasVector - Input bias vector doesn't have any non-underflow/overflow bins!");

    // Setup a new histogram so we can style it without modifying the input
    auto biasVector = std::make_shared<TH1F>(("smearingBias_" + std::to_string(m_lastPlotId++)).c_str(), "", nBins, 0, nBins);

    float maxBias = -std::numeric_limits<float>::max();
    for (unsigned int iBin = 1u; iBin <= nBins; ++iBin)
    {
        // Get the corresponding reco-true bins
        const unsigned int iTrue = (iBin-1) / nBinsSqrt;
        const unsigned int iReco = (iBin-1) % nBinsSqrt;

        // Set the bin labels
        const auto isTrueUnderflow = (iTrue == 0 && hasUnderflow);
        const auto isRecoUnderflow = (iReco == 0 && hasUnderflow);

        const auto isTrueOverflow = (iTrue == (nBinsSqrt-1) && hasOverflow);
        const auto isRecoOverflow = (iReco == (nBinsSqrt-1) && hasOverflow);

        const std::string trueBinName = isTrueUnderflow ? "UF" : (isTrueOverflow ? "OF" : std::to_string(iTrue - (hasUnderflow ? 1u : 0u)));
        const std::string recoBinName = isRecoUnderflow ? "UF" : (isRecoOverflow ? "OF" : std::to_string(iReco - (hasUnderflow ? 1u : 0u)));

        biasVector->GetXaxis()->SetBinLabel(iBin, ("T-" + trueBinName + " R-" + recoBinName).c_str());

        // Store this bias value
        const float bias = vector->GetBinContent(iBin);

        // Choose a nice y-axis range
        const auto padding = 1.05f;
        maxBias = std::max(maxBias, std::abs(bias) * padding);

        biasVector->SetBinContent(iBin, bias);
    }

    biasVector->GetYaxis()->SetRangeUser(-maxBias, +maxBias);
    biasVector->GetXaxis()->LabelsOption("v");
    biasVector->SetFillColor(PlottingHelper::GetColor(Default));

    // Draw the histogram
    gStyle->SetHistMinimumZero(true);
    auto pCanvas = PlottingHelper::GetCanvas();
    pCanvas->SetBottomMargin(0.15f);
    biasVector->Draw("bar");

    // Add lines to show the original binning
    std::vector<std::shared_ptr<TLine> > lines;
    for (unsigned int iBin = 0; iBin <= nBinsSqrt; ++iBin)
    {
        const auto binEdge = iBin * nBinsSqrt;
        lines.emplace_back(std::make_shared<TLine>(binEdge, -maxBias, binEdge, maxBias));
        auto &pLine = lines.back();
        pLine->SetLineWidth(2);
        pLine->Draw();
    }

    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_smearing_bias");
    gStyle->SetHistMinimumZero(false);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::SaveCovarianceMatrix(const std::shared_ptr<TH2F> &matrix, const std::shared_ptr<TH1F> &crossSection, const bool hasUnderflow, const bool hasOverflow, const std::string &namePrefix)
{
    // Check we have sensible binning
    const unsigned int nBins = crossSection->GetNbinsX();

    const auto nUnderOverflowBins = (hasUnderflow ? 1u : 0u) + (hasOverflow ? 1u : 0u);
    if (nBins <= nUnderOverflowBins)
        throw std::invalid_argument("PlottingHelper::SaveCovarianceMatrix - Input cross section doesn't have any non-underflow/overflow bins!");

    const auto nAnalysisBins = nBins - nUnderOverflowBins;

    const unsigned int matrixBinsX = matrix->GetNbinsX();
    const unsigned int matrixBinsY = matrix->GetNbinsY();

    if (matrixBinsX != matrixBinsY)
        throw std::invalid_argument("PlottingHelper::SaveCovarianceMatrix - Input covariance matrix is not square!");

    if (matrixBinsX != nAnalysisBins)
    {
        throw std::invalid_argument("PlottingHelper::SaveCovarianceMatrix - Input covariance matrix has " + std::to_string(matrixBinsX) +
            " bins, but input cross-section has " + std::to_string(nBins) + " bins, of which " + std::to_string(nUnderOverflowBins) +
            " are under/overflow bins - this doesn't add up!");
    }

    // Make the fractional covariance and linear correlation matricies
    auto covarianceMatrix = std::make_shared<TH2F>(("matrix_" + std::to_string(m_lastPlotId++)).c_str(), "", nAnalysisBins, 0, nAnalysisBins, nAnalysisBins, 0, nAnalysisBins);
    auto fracCovarianceMatrix = std::make_shared<TH2F>(("matrix_" + std::to_string(m_lastPlotId++)).c_str(), "", nAnalysisBins, 0, nAnalysisBins, nAnalysisBins, 0, nAnalysisBins);
    auto correlationMatrix = std::make_shared<TH2F>(("matrix_" + std::to_string(m_lastPlotId++)).c_str(), "", nAnalysisBins, 0, nAnalysisBins, nAnalysisBins, 0, nAnalysisBins);

    // Set the number of decimal places & text size to print in each bin
    gStyle->SetPaintTextFormat("2.2f");
    covarianceMatrix->SetMarkerSize(2.2);
    fracCovarianceMatrix->SetMarkerSize(2.2);
    correlationMatrix->SetMarkerSize(2.2);

    for (unsigned int iBin = 1u; iBin <= nAnalysisBins; ++iBin)
    {
        // Set the bin labels of each plot (here we enumerate from zero)
        for (auto &pHist : {covarianceMatrix, fracCovarianceMatrix, correlationMatrix})
        {
            for (auto &pAxis : {pHist->GetXaxis(), pHist->GetYaxis()})
            {
                pAxis->SetBinLabel(iBin, std::to_string(iBin - 1).c_str());
            }
        }

        // Get the cross section in bin i
        const auto xSecI = crossSection->GetBinContent(iBin + (hasUnderflow ? 1u : 0u));
        if (xSecI <= std::numeric_limits<float>::epsilon())
            throw std::logic_error("PlottingHelper::SaveCovarianceMatrix - Cross section in analysis bin: " + std::to_string(iBin) + " is <= 0");

        for (unsigned int jBin = 1u; jBin <= nAnalysisBins; ++jBin)
        {
            // Get the cross section in bin j
            const auto xSecJ = crossSection->GetBinContent(jBin + (hasUnderflow ? 1u : 0u));
            if (xSecJ <= std::numeric_limits<float>::epsilon())
                throw std::logic_error("PlottingHelper::SaveCovarianceMatrix - Cross section in analysis bin: " + std::to_string(jBin) + " is <= 0");

            // Set the covariance
            const auto covariance = matrix->GetBinContent(iBin, jBin);
            covarianceMatrix->SetBinContent(iBin, jBin, covariance);

            // Set the fractional covariance
            const auto fracCovariance = covariance / (xSecI * xSecJ);
            fracCovarianceMatrix->SetBinContent(iBin, jBin, fracCovariance);

            // Get the standard deviation in bin i
            const auto varI = matrix->GetBinContent(iBin, iBin);
            if (varI <= std::numeric_limits<float>::epsilon())
            {
                std::cout << "PlottingHelper::SaveCovarianceMatrix - WARNING - Variance in bin: " << iBin << " is " << varI << ". Setting to correlation matrix element to dummy value"<< std::endl;
                correlationMatrix->SetBinContent(iBin, jBin, -std::numeric_limits<float>::max());
                continue;
            }

            const auto stdI = std::pow(varI, 0.5);

            // Get the standard deviation in bin j
            const auto varJ = matrix->GetBinContent(jBin, jBin);
            if (varJ <= std::numeric_limits<float>::epsilon())
            {
                std::cout << "PlottingHelper::SaveCovarianceMatrix - WARNING - Variance in bin: " << jBin << " is " << varJ << ". Setting to correlation matrix element to dummy value"<< std::endl;
                correlationMatrix->SetBinContent(iBin, jBin, -std::numeric_limits<float>::max());
                continue;
            }

            const auto stdJ = std::pow(varJ, 0.5);

            // Set the linear correlation
            const auto correlation = covariance / (stdI * stdJ);
            correlationMatrix->SetBinContent(iBin, jBin, correlation);
        }
    }

    // Make the plots
    auto pCanvas = PlottingHelper::GetCanvas(960, 960);

    covarianceMatrix->Draw("colz text");
    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_covariance");

    fracCovarianceMatrix->Draw("colz text");
    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_fracCovariance");

    correlationMatrix->Draw("colz text");
    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_linearCorrelation");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::SaveSmearingMatrixCovarianceMatrix(const std::shared_ptr<TH2F> &matrix, const bool hasUnderflow, const bool hasOverflow, const std::string &namePrefix)
{
    // Check we have sensible binning
    const unsigned int matrixBinsX = matrix->GetNbinsX();
    const unsigned int matrixBinsY = matrix->GetNbinsY();

    if (matrixBinsX != matrixBinsY)
        throw std::invalid_argument("PlottingHelper::SaveSmearingMatrixCovarianceMatrix - Input covariance matrix is not square!");

    // Get the number of bins of the original cross-section (this smearing matrix should have the square of this number of bins)
    const auto nBins = matrixBinsX;
    const auto nBinsSqrt = static_cast<unsigned int>(std::round(std::pow(static_cast<float>(nBins), 0.5f)));
    if (nBinsSqrt*nBinsSqrt != nBins)
        throw std::invalid_argument("PlottingHelper::SaveSmearingMatrixCovarianceMatrix - Number of bins for input bias vector isn't a perfect square");

    const auto nUnderOverflowBins = (hasUnderflow ? 1u : 0u) + (hasOverflow ? 1u : 0u);
    if (nBinsSqrt <= nUnderOverflowBins)
        throw std::invalid_argument("PlottingHelper::SaveSmearingMatrixCovarianceMatrix - Input cross section doesn't have any non-underflow/overflow bins!");

    // Setup a new histogram so we can style it without modifying the input
    auto covarianceMatrix = std::make_shared<TH2F>(("matrix_" + std::to_string(m_lastPlotId++)).c_str(), "", nBins, 0, nBins, nBins, 0, nBins);

    // Loop over the bins
    for (unsigned int iBin = 1u; iBin <= nBins; ++iBin)
    {
        // Get the corresponding reco-true bins
        const unsigned int iTrue = (iBin-1) / nBinsSqrt;
        const unsigned int iReco = (iBin-1) % nBinsSqrt;

        // Set the bin labels
        const auto isTrueUnderflow = (iTrue == 0 && hasUnderflow);
        const auto isRecoUnderflow = (iReco == 0 && hasUnderflow);

        const auto isTrueOverflow = (iTrue == (nBinsSqrt-1) && hasOverflow);
        const auto isRecoOverflow = (iReco == (nBinsSqrt-1) && hasOverflow);

        const std::string trueBinName = isTrueUnderflow ? "UF" : (isTrueOverflow ? "OF" : std::to_string(iTrue - (hasUnderflow ? 1u : 0u)));
        const std::string recoBinName = isRecoUnderflow ? "UF" : (isRecoOverflow ? "OF" : std::to_string(iReco - (hasUnderflow ? 1u : 0u)));

        const std::string binLabel = "T-" + trueBinName + " R-" + recoBinName;
        covarianceMatrix->GetXaxis()->SetBinLabel(iBin, binLabel.c_str());
        covarianceMatrix->GetYaxis()->SetBinLabel(iBin, binLabel.c_str());

        // Copy the content of the bins
        for (unsigned int jBin = 1u; jBin <= nBins; ++jBin)
        {
            const auto binContent = matrix->GetBinContent(iBin, jBin);
            covarianceMatrix->SetBinContent(iBin, jBin, binContent);
        }
    }

    // Draw the histogram
    auto pCanvas = PlottingHelper::GetCanvas(960, 960);
    covarianceMatrix->GetXaxis()->LabelsOption("v");
    covarianceMatrix->GetYaxis()->LabelsOption("h");
    covarianceMatrix->Draw("colz");

    // Add lines to show the original binning
    std::vector<std::shared_ptr<TLine> > lines;
    for (unsigned int iBin = 0; iBin <= nBinsSqrt; ++iBin)
    {
        const auto binEdge = iBin * nBinsSqrt;
        lines.emplace_back(std::make_shared<TLine>(0, binEdge, nBins, binEdge));
        auto &pLineX = lines.back();
        pLineX->SetLineWidth(2);
        pLineX->Draw();

        lines.emplace_back(std::make_shared<TLine>(binEdge, 0, binEdge, nBins));
        auto &pLineY = lines.back();
        pLineY->SetLineWidth(2);
        pLineY->Draw();
    }

    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_smearing_covariance");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::SaveCrossSection(const std::shared_ptr<TH1F> &crossSection, const std::shared_ptr<TH1F> &statUncertainty, const std::shared_ptr<TH2F> &totalSystCovarianceMatrix, const bool hasUnderflow, const bool hasOverflow, const std::string &namePrefix)
{
    // Check we have sensible binning
    const unsigned int nBins = crossSection->GetNbinsX();

    const auto nUnderOverflowBins = (hasUnderflow ? 1u : 0u) + (hasOverflow ? 1u : 0u);
    if (nBins <= nUnderOverflowBins)
        throw std::invalid_argument("PlottingHelper::SaveCrossSection - Input cross section doesn't have any non-underflow/overflow bins!");

    const auto nAnalysisBins = nBins - nUnderOverflowBins;

    const unsigned int statUncertaintyBins = statUncertainty->GetNbinsX();
    if (statUncertaintyBins != nAnalysisBins)
    {
        throw std::invalid_argument("PlottingHelper::SaveCrossSection - Input cross-section stat uncertainty has " + std::to_string(nAnalysisBins) +
            " bins, but input cross-section has " + std::to_string(nBins) + " bins, of which " + std::to_string(nUnderOverflowBins) +
            " are under/overflow bins - this doesn't add up!");
    }

    const unsigned int matrixBinsX = totalSystCovarianceMatrix->GetNbinsX();
    const unsigned int matrixBinsY = totalSystCovarianceMatrix->GetNbinsY();

    if (matrixBinsX != matrixBinsY)
        throw std::invalid_argument("PlottingHelper::SaveCrossSection - Input covariance matrix is not square!");

    if (matrixBinsX != nAnalysisBins)
    {
        throw std::invalid_argument("PlottingHelper::SaveCrossSection - Input covariance matrix has " + std::to_string(matrixBinsX) +
            " bins, but input cross-section has " + std::to_string(nBins) + " bins, of which " + std::to_string(nUnderOverflowBins) +
            " are under/overflow bins - this doesn't add up!");
    }

    // Get the bin edges from the cross-section histogram
    std::vector<float> binEdges;
    const auto firstBin = 1u + (hasUnderflow ? 1u : 0u);
    const auto lastBin = nBins - (hasOverflow ? 1u : 0u);
    for (unsigned int iBin = firstBin; iBin <= lastBin; ++iBin)
    {
        binEdges.push_back(crossSection->GetBinLowEdge(iBin));
    }
    binEdges.push_back(crossSection->GetBinLowEdge(lastBin) + crossSection->GetBinWidth(lastBin));

    // Perform a sanity check just for debugging purposes
    if (binEdges.size() - 1 != nAnalysisBins)
        throw std::logic_error("PlottingHelper::SaveCrossSection - Sanity check failed! Something went wrong getting the bin edges");

    // Make new histograms for drawing. Here we make 2 histograms with the same values, but one has the stat-only uncertainty and the other
    // has the total uncertainty
    auto xsecStatOnly = std::make_shared<TH1F>(("xsecWithErr_" + std::to_string(m_lastPlotId++)).c_str(), "", nAnalysisBins, binEdges.data());
    auto xsecTotalErr = std::make_shared<TH1F>(("xsecWithErr_" + std::to_string(m_lastPlotId++)).c_str(), "", nAnalysisBins, binEdges.data());

    float maxValue = -std::numeric_limits<float>::max();
    float minValue = +std::numeric_limits<float>::max();

    for (unsigned int iBin = 1u; iBin <= nAnalysisBins; ++iBin)
    {
        // Get the cross-section value and the stat & systematic uncertaint
        // NB. Here we just use the diagonal of the covarianc matrix for the plot
        const float xsec = crossSection->GetBinContent(iBin + (hasUnderflow ? 1u : 0u));
        const float statErr = statUncertainty->GetBinContent(iBin);
        const float systErr2 = totalSystCovarianceMatrix->GetBinContent(iBin, iBin);
        const float totalErr = std::pow(statErr*statErr + systErr2, 0.5f);

        maxValue = std::max(maxValue, xsec + totalErr);
        minValue = std::min(minValue, xsec - totalErr);

        xsecStatOnly->SetBinContent(iBin, xsec);
        xsecStatOnly->SetBinError(iBin, statErr);

        xsecTotalErr->SetBinContent(iBin, xsec);
        xsecTotalErr->SetBinError(iBin, totalErr);
    }

    // Make the plots
    auto pCanvas = PlottingHelper::GetCanvas();
    gStyle->SetEndErrorSize(4);

    // Set the axis ranges
    const auto padding = (maxValue - minValue) * 0.05;
    minValue -= padding;
    maxValue += padding;
    minValue = std::max(0.f, minValue);
    xsecStatOnly->GetYaxis()->SetRangeUser(minValue, maxValue);
    xsecTotalErr->GetYaxis()->SetRangeUser(minValue, maxValue);

    // Draw the data points as error bars
    PlottingHelper::SetLineStyle(xsecStatOnly, Default);
    PlottingHelper::SetLineStyle(xsecTotalErr, Default);

    xsecStatOnly->Draw("e1");
    xsecTotalErr->Draw("e1 same");
    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_xsecWithErrors");
    gStyle->SetEndErrorSize(2);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::SaveSmearingMatrix(const std::shared_ptr<TH2F> &matrix, const bool hasUnderflow, const bool hasOverflow, const std::string &namePrefix)
{
    // Check we have sensible binning
    const unsigned int matrixBinsX = matrix->GetNbinsX();
    const unsigned int matrixBinsY = matrix->GetNbinsY();

    if (matrixBinsX != matrixBinsY)
        throw std::invalid_argument("PlottingHelper::SaveSmearingMatrix - Input smearing matrix is not square!");

    const auto nBins = matrixBinsX;
    const auto nUnderOverflowBins = (hasUnderflow ? 1u : 0u) + (hasOverflow ? 1u : 0u);
    if (nBins <= nUnderOverflowBins)
        throw std::invalid_argument("PlottingHelper::SaveSmearingMatrix - Input cross section doesn't have any non-underflow/overflow bins!");

    // Make a new histogram so we can style it without modifying the original
    auto smearingMatrix = std::make_shared<TH2F>(("matrix_" + std::to_string(m_lastPlotId++)).c_str(), "", nBins, 0, nBins, nBins, 0, nBins);
    for (unsigned int iTrue = 1u; iTrue <= nBins; ++iTrue)
    {
        // Set the bin labels
        const auto isUnderflow = (iTrue == 1 && hasUnderflow);
        const auto isOverflow = (iTrue == nBins && hasOverflow);
        const std::string binName = isUnderflow ? "UF" : (isOverflow ? "OF" : std::to_string(iTrue - 1u - (hasUnderflow ? 1u : 0u)));

        smearingMatrix->GetXaxis()->SetBinLabel(iTrue, ("T-" + binName).c_str());
        smearingMatrix->GetYaxis()->SetBinLabel(iTrue, ("R-" + binName).c_str());

        // Copy the bin contents
        for (unsigned int iReco = 1u; iReco <= nBins; ++iReco)
        {
            const auto binContent = matrix->GetBinContent(iTrue, iReco);
            smearingMatrix->SetBinContent(iTrue, iReco, binContent);
        }
    }

    // Draw the histogram
    auto pCanvas = PlottingHelper::GetCanvas(960, 960);
    gStyle->SetPaintTextFormat("2.2f");
    smearingMatrix->SetMarkerSize(1.2);
    smearingMatrix->GetXaxis()->LabelsOption("v");
    smearingMatrix->GetYaxis()->LabelsOption("h");
    smearingMatrix->Draw("colz text");

    PlottingHelper::SaveCanvas(pCanvas, namePrefix + "_smearing");
}
*/
// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::PlotErrorMatrix(const ubsmear::UBMatrix &errorMatrix, const std::string &fileName, const ubsmear::UBXSecMeta &metadata, const bool useDefaultPalette, const bool useSymmetricZRange, const float zRangeMax)
{
    if (!ubsmear::UBMatrixHelper::IsSquare(errorMatrix))
        throw std::invalid_argument("PlottingHelper::PlotErrorMatrix - Input error matrix isn't square");

    const auto nBins = errorMatrix.GetRows();

    // Check we have a valid number of bins
    const auto nAnalysisBins = metadata.GetNBins() - (metadata.HasUnderflow() ? 1 : 0) - (metadata.HasOverflow() ? 1 : 0);
    const auto nSmearingBins = metadata.GetNBins()*metadata.GetNBins();
    if (nBins != metadata.GetNBins() && nBins != nAnalysisBins && nBins != nSmearingBins)
        throw std::invalid_argument("PlottingHelper::PlotErrorMatrix - Number of bins of error matrix doesn't match input metadata");

    const auto isSmearingMatrix = (nBins == nSmearingBins);
    const auto showNumbers = !isSmearingMatrix;
    const auto nBinsSqrt = static_cast<unsigned int>(std::round(std::pow(nBins, 0.5)));

    // When we have lots of bin labels it can be hard to read, here we set the maximum number of labels to print
    const auto maxLabels = 30u;

    // If we have less than the maximum number of labels to print, then just print them all (bin size 1)
    // Otherwise, set the largest possible group size so we don't exceed max labels
    const unsigned int labelGroupSize = std::ceil(static_cast<float>(nBins) / static_cast<float>(maxLabels));

    // Setup a vector of lines to add if this is a smearing matrix error matrix
    std::vector< shared_ptr<TLine> > lines;
    if (isSmearingMatrix)
    {
        for (unsigned int iGroup = 1; iGroup < nBinsSqrt; ++iGroup)
        {
            const auto linePos = iGroup * nBinsSqrt;

            // Add the horizontal line for this group
            lines.emplace_back(new TLine(0, linePos, nBins, linePos));

            // Add the vertical line for this group
            lines.emplace_back(new TLine(linePos, 0, linePos, nBins));
        }
    }

    // Make the histogram for this matrix
    float maxBinContent = 0.f;
    auto pHist = std::make_shared<TH2F>(("errorMatrix_" + std::to_string(m_lastPlotId++)).c_str(), "", nBins, 0, nBins, nBins, 0, nBins);
    for (unsigned int iRow = 0; iRow < nBins; ++iRow)
    {
        // Set the bin labels
        std::string binLabel = "";
        if (iRow % labelGroupSize == 0)
        {
            if (isSmearingMatrix)
            {
                const unsigned int iReco = iRow / nBinsSqrt;
                const unsigned int iTrue = iRow % nBinsSqrt;

                const auto isUnderflowReco = metadata.HasUnderflow() && (iReco == 0);
                const auto isOverflowReco = metadata.HasOverflow() && (iReco == metadata.GetNBins() - 1);
                const auto recoBinName = isUnderflowReco ? "UF" : (isOverflowReco ? "OF" : std::to_string(iReco - (metadata.HasUnderflow() ? 1 : 0)));

                const auto isUnderflowTrue = metadata.HasUnderflow() && (iTrue == 0);
                const auto isOverflowTrue = metadata.HasOverflow() && (iTrue == metadata.GetNBins() - 1);
                const auto trueBinName = isUnderflowTrue ? "UF" : (isOverflowTrue ? "OF" : std::to_string(iTrue - (metadata.HasUnderflow() ? 1 : 0)));

                binLabel = "R-" + recoBinName + " / T-" + trueBinName;
            }
            else
            {
                const auto isUnderflow = metadata.HasUnderflow() && (iRow == 0);
                const auto isOverflow = metadata.HasOverflow() && (iRow == metadata.GetNBins() - 1);
                binLabel = isUnderflow ? "UF" : (isOverflow ? "OF" : std::to_string(iRow - (metadata.HasUnderflow() ? 1 : 0)));
            }
        }

        // ATTN we +1 to the indices as ROOT enumerates from 1
        pHist->GetXaxis()->SetBinLabel(iRow + 1, binLabel.c_str());
        pHist->GetYaxis()->SetBinLabel(iRow + 1, binLabel.c_str());

        for (unsigned int iCol = 0; iCol < nBins; ++iCol)
        {
            float binContent = errorMatrix.At(iRow, iCol);

            // If the bin is empty then it will be drawn as as a blank square, here we set it to small value to force ROOT to fill the square
            if (std::abs(binContent) < std::numeric_limits<float>::epsilon())
                binContent = std::numeric_limits<float>::epsilon();

            pHist->SetBinContent(iCol + 1, iRow + 1, binContent);

            maxBinContent = std::max(maxBinContent, std::abs(binContent));
        }
    }

    // If desired, used the supplied maximum bin range
    if (zRangeMax > 0.f && useSymmetricZRange)
    {
        maxBinContent = zRangeMax;
    }

    // Plot the histogram
    auto pCanvas = PlottingHelper::GetCanvas(960, 960);
    pCanvas->SetTopMargin(0.10f);
    pCanvas->SetBottomMargin(0.20f);
    pCanvas->SetLeftMargin(0.17f);
    pCanvas->SetRightMargin(0.12f);

    // To make the zero-point clear, we choose an alternate colour mapping
    // Make some vectors that store the R, G and B values at different "stops" between 0 -> 1
    if (!useDefaultPalette)
    {
        std::vector<double> reds =   {  22./255., 255./255., 255./255. };
        std::vector<double> greens = { 103./255., 255./255., 149./255. };
        std::vector<double> blues =  { 255./255., 255./255.,   0./255. };
        std::vector<double> stops =  { 0.       , 0.5      , 1.        };
        const unsigned int nColourPoints = 50u;

        // Get the colour palette
        auto firstIndex = TColor::CreateGradientColorTable(stops.size(), stops.data(), reds.data(), greens.data(), blues.data(), nColourPoints);
        std::vector<int> palette;
        for (unsigned int i = 0; i < nColourPoints; ++i)
        {
            palette.push_back(firstIndex + i);
        }

        // Use this new palette
        gStyle->SetPalette(nColourPoints, palette.data());
    }

    // Set the range so its equal in the positive and negative directions
    if (useSymmetricZRange)
    {
        pHist->GetZaxis()->SetRangeUser(-maxBinContent, +maxBinContent);
    }

    // Remove the ticks
    pHist->GetXaxis()->SetTickLength(0.f);
    pHist->GetYaxis()->SetTickLength(0.f);

    // Rotate the labels on the x-axis
    if (isSmearingMatrix)
    {
        pHist->GetXaxis()->LabelsOption("v");
    }

    if (showNumbers)
    {
        // Decide on the marker size
        const float markerSizeMin = 1.f;
        const float markerSizeMax = 2.5f;
        const unsigned int nBinsMarkerMin = 1u;
        const unsigned int nBinsMarkerMax = 15u;

        //   - At nBins <= nBinsMarkerMin we want to be at markerSizeMax
        //   - At nBins >= nBinsMarkerMax we want to be at markerSizeMin
        //   - Between nBinsMarkerMin <= nBins <= nBinsMarkerThreshold, marker size should scale with 1/nBins
        //   - Choose form: markerSize = A / nBins + B
        const auto markerA = (markerSizeMax - markerSizeMin) * (nBinsMarkerMin * nBinsMarkerMax) / static_cast<float>(nBinsMarkerMax - nBinsMarkerMin);
        const auto markerB = markerSizeMax - (markerA / static_cast<float>(nBinsMarkerMin));
        const auto markerSize = (markerA / static_cast<float>(nBins)) + markerB;

        gStyle->SetPaintTextFormat("2.2f");
        pHist->SetMarkerSize(markerSize);
        pHist->Draw("colz text");
    }
    else
    {
        pHist->Draw("colz");
    }

    // Draw the lines
    for (const auto &pLine : lines)
    {
        pLine->Draw();
    }

    // Save the canvas
    PlottingHelper::SaveCanvas(pCanvas, fileName);

    // Set the palette back to default
    gStyle->SetPalette();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::PlotFractionalErrorMatrix(const ubsmear::UBMatrix &errorMatrix, const ubsmear::UBMatrix &quantityVector, const std::string &fileName, const ubsmear::UBXSecMeta &metadata, const bool useDefaultPalette, const bool useSymmetricZRange)
{
    // Check that the dimensions of the inputs are compatible
    if (quantityVector.GetColumns() != 1)
        throw std::invalid_argument("PlottingHelper::PlotErrorMatrix - Input quantity vector isn't a column vector");

    if (!ubsmear::UBMatrixHelper::IsSquare(errorMatrix))
        throw std::invalid_argument("PlottingHelper::PlotErrorMatrix - Input error matrix isn't square");

    const auto nBins = errorMatrix.GetRows();
    if (quantityVector.GetRows() != nBins)
        throw std::invalid_argument("PlottingHelper::PlotErrorMatrix - Input error matrix and quantity vetor have incompatible dimensions");

    // Scale down the error matrix to get the fractional error matrix
    auto fracErrorMatrix = ubsmear::UBMatrixHelper::GetZeroMatrix(nBins, nBins);
    for (unsigned int iRow = 0; iRow < nBins; ++iRow)
    {
        const auto quantityRow = quantityVector.At(iRow, 0);

        // Check if the quantity is zero
        if (std::abs(quantityRow) <= std::numeric_limits<float>::epsilon())
            continue;

        for (unsigned int iCol = 0; iCol < nBins; ++iCol)
        {
            const auto quantityCol = quantityVector.At(iCol, 0);

            // Check if the quantity is zero
            if (std::abs(quantityCol) <= std::numeric_limits<float>::epsilon())
                continue;

            const auto fracError = errorMatrix.At(iRow, iCol) / (quantityRow * quantityCol);
            fracErrorMatrix.SetElement(iRow, iCol, fracError);
        }
    }

    // Plot the fractional error matrix
    PlottingHelper::PlotErrorMatrix(fracErrorMatrix, fileName, metadata, useDefaultPalette, useSymmetricZRange);

    // Plot the fractional error matrix with a fixed scale
    PlottingHelper::PlotErrorMatrix(fracErrorMatrix, fileName + "_fixed", metadata, useDefaultPalette, useSymmetricZRange, 1.f);
}

} // namespace ubcc1pi
