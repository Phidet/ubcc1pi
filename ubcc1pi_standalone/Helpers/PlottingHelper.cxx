#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include <stdexcept>

#include <THStack.h>

namespace ubcc1pi
{

PlottingHelper::MultiPlot::MultiPlot(const std::string &xLabel, const std::string &yLabel, unsigned int nBins, float min, float max, bool drawErrors) :
    m_xLabel(xLabel),
    m_nBins(nBins),
    m_min(min),
    m_max(max),
    m_id(++m_lastId),
    m_cloneCount(0),
    m_drawErrors(drawErrors)
{
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Make a unique name for this plot to avoid collisionsi
        const auto nameStr = "ubcc1pi_multiPlot_" + std::to_string(m_id) + "_" + std::to_string(static_cast<int>(style));
        const auto name = nameStr.c_str();

        auto pHist = std::make_shared<TH1F>(name, "", nBins, min, max);
        
        pHist->Sumw2();
        PlottingHelper::SetLineStyle(pHist.get(), style);
        pHist->GetXaxis()->SetTitle(xLabel.c_str());
        pHist->GetYaxis()->SetTitle(yLabel.c_str());
        
        m_plotToHistMap.emplace(style, pHist);
    }
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

        plotToHistCloneMap.emplace(style, pHistClone);
    }

    m_cloneCount++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::ScaleHistograms(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap) const
{
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto pHist = plotToHistCloneMap.at(style);
        pHist->Scale(1.f / pHist->Integral());
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SetHistogramYRanges(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap) const
{
    float yMin = std::numeric_limits<float>::max();
    float yMax = -std::numeric_limits<float>::max();
    
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        const auto pHist = plotToHistCloneMap.at(style);
        yMin = std::min(yMin, static_cast<float>(pHist->GetMinimum()));
        yMax = std::max(yMax, static_cast<float>(pHist->GetMaximum()));
    }
    
    // Add some padding to the top of the histogram
    yMax += (yMax - yMin) * 0.05;
    
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        auto pHist = plotToHistCloneMap.at(style);
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SaveAs(const std::string &fileName)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<PlotStyle, TH1F* > plotToHistCloneMap;
    this->GetHistogramClones(plotToHistCloneMap);
    this->ScaleHistograms(plotToHistCloneMap);
    this->SetHistogramYRanges(plotToHistCloneMap);

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
    
            if (pHist->GetEntries() == 0)
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

        if (pHist->GetEntries() == 0)
            continue;
        
        const bool usePoints = PlottingHelper::ShouldUsePoints(style);
        const TString tstyle = usePoints ? "e1" : "hist";
        pHist->Draw(isFirst ? tstyle : tstyle + "same");
        isFirst = false;
    }

    PlottingHelper::SaveCanvas(pCanvas, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::MultiPlot::SaveAsStacked(const std::string &fileName)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    // Clone the histogtams and scale them (we clone to the original hisograms can be subsequently filled)
    std::unordered_map<PlotStyle, TH1F* > plotToHistCloneMap;
    this->GetHistogramClones(plotToHistCloneMap);

    // Work out if we have BNB data to plot, and if we do then get the minimum and maximum Y coordinates
    auto bnbDataHistIter = plotToHistCloneMap.find(BNBData);
    const bool hasBNBData = bnbDataHistIter != plotToHistCloneMap.end();
    auto yMin = hasBNBData ? static_cast<float>(bnbDataHistIter->second->GetMinimum()) : std::numeric_limits<float>::max();
    auto yMax = hasBNBData ? static_cast<float>(bnbDataHistIter->second->GetMaximum()) : -std::numeric_limits<float>::max();

    // Sum the non BNB data histograms to get the "MC" total
    const auto nameTotalStr = "ubcc1pi_multiPlot_" + std::to_string(m_id) + "_total";
    const auto nameTotal = nameTotalStr.c_str();
    auto pHistTotal = std::make_shared<TH1F>(nameTotal, "", m_nBins, m_min, m_max);
    pHistTotal->Sumw2();

    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        // Don't add BNB data to the total, this will be drawn seperately
        if (style == BNBData)
            continue;

        auto pHist = plotToHistCloneMap.at(style);
        
        std::cout << "Histogram: " << pHist->GetName() << " - nEntries = " << pHist->GetEntries() << std::endl;

        if (pHist->GetEntries() == 0)
            continue;

        pHistTotal->Add(pHist);
    }

    // Get the maximum and minimum Y coordinates
    yMin = std::min(yMin, static_cast<float>(pHistTotal->GetMinimum()));
    yMax = std::min(yMax, static_cast<float>(pHistTotal->GetMaximum()));
    yMax += (yMax - yMin) * 0.05;

    // Draw the stacked histogram
    const auto nameStackStr = "ubcc1pi_plotPlot_" + std::to_string(m_id) + "_stack";
    const auto nameStack = nameStackStr.c_str();
    auto pHistStack = std::make_shared<THStack>(nameStack, "");
    
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
       
        pHistStack->Add(pHist);
    }

    pHistStack->Draw("hist");
    pHistStack->GetYaxis()->SetRangeUser(yMin, yMax);
    pHistStack->Draw("hist");
           
    // Draw the error bands on the stack if required
    if (m_drawErrors)
    {
        pHistTotal->SetFillStyle(1001);
        pHistTotal->SetLineColorAlpha(kBlack, 0.f);
        pHistTotal->SetFillColorAlpha(kBlack, 0.3f);
        pHistTotal->Draw("e2 same");
    }

    // Draw the BNB data
    if (hasBNBData)
    {
        bnbDataHistIter->second->Draw("e1 same");
    }

    PlottingHelper::SaveCanvas(pCanvas, fileName);

    // Now make the ratio plot
    auto pCanvasRatio = PlottingHelper::GetCanvas(960, 270);
    auto pHistBNBData = bnbDataHistIter->second;
    pHistBNBData->GetYaxis()->SetTitle("Beam on / (Overlay + Beam off)");
    pHistBNBData->Divide(pHistTotal.get());
    pHistBNBData->Draw("e1");
    PlottingHelper::SaveCanvas(pCanvasRatio, fileName + "_ratio");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

PlottingHelper::EfficiencyPlot::EfficiencyPlot(const std::string &xLabel, unsigned int nBins, float min, float max, const std::vector<string> &cuts, bool drawErrors) :
    m_xLabel(xLabel),
    m_nBins(nBins),
    m_min(min),
    m_max(max),
    m_cuts(cuts),
    m_drawErrors(drawErrors),
    m_id(++m_lastId)
{
    // Get all colours
    const auto palette = PlottingHelper::GetColorVector();
    
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
            std::make_shared<TH1F>(nameNumerator.c_str(), "", nBins, min, max),
            std::make_shared<TH1F>(nameDenominator.c_str(), "", nBins, min, max));
        
        // Pick the color from the palette
        const auto col = palette.at(iCut % palette.size());

        // Setup the plots
        for (auto &pHist : {plotPair.first, plotPair.second})
        {
            pHist->Sumw2();
            PlottingHelper::SetLineStyle(pHist.get(), col);
            pHist->GetXaxis()->SetTitle(xLabel.c_str());
            pHist->GetYaxis()->SetTitle("Efficiency");
        }
        
        m_cutToPlotsMap.emplace(cut, plotPair);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
       
void PlottingHelper::EfficiencyPlot::AddEvent(const float value, const std::string &cut, const bool passedCut)
{
    auto plotPairIter = m_cutToPlotsMap.find(cut);

    if (plotPairIter == m_cutToPlotsMap.end())
        throw std::invalid_argument("EfficiencyPlot::AddEvent - Unknown cut: \"" + cut + "\"");

    auto &pHistNumerator = plotPairIter->second.first;
    auto &pHistDenominator = plotPairIter->second.second;

    if (passedCut)
        pHistNumerator->Fill(value);

    pHistDenominator->Fill(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void PlottingHelper::EfficiencyPlot::SaveAs(const std::string &fileName)
{
    auto pCanvas = PlottingHelper::GetCanvas();

    // Save the raw histogram - just use the denominator of the first cut
    auto pHistRaw = static_cast<TH1F *>(m_cutToPlotsMap.at(m_cuts.front()).second->Clone());
    pHistRaw->GetYaxis()->SetTitle("Number of events");
    pHistRaw->SetLineColor(kBlack);

    if (m_drawErrors)
    {
        auto pHistClone = static_cast<TH1F *>(pHistRaw->Clone());
        pHistClone->SetFillStyle(1001);
        pHistClone->SetLineColorAlpha(kBlack, 0.f);
        pHistClone->SetFillColorAlpha(kBlack, 0.3f);

        pHistClone->Draw("e2");
    }
    
    pHistRaw->Draw(m_drawErrors ? "hist same" : "hist"); 
    PlottingHelper::SaveCanvas(pCanvas, fileName + "_raw");

    // Get the efficiency histograms
    std::vector<TH1F *> efficiencyHists;
    std::vector<TH1F *> efficiencyErrorHists;
    for (const auto &cut : m_cuts)
    {
        const auto &plotPair = m_cutToPlotsMap.at(cut);
        const auto &pHistNumerator = plotPair.first;
        const auto &pHistDenominator = plotPair.second;

        // Clone the numerator and divide by the denominator to get the efficiency
        auto pHistEfficiency = static_cast<TH1F *>(pHistNumerator->Clone());
        pHistEfficiency->Divide(pHistDenominator.get());

        efficiencyHists.push_back(pHistEfficiency);
    
        if (m_drawErrors)
        {
            auto pHistClone = static_cast<TH1F *>(pHistEfficiency->Clone());
            const auto col = pHistClone->GetLineColor();
            pHistClone->SetFillStyle(1001);
            pHistClone->SetLineColorAlpha(col, 0.f);
            pHistClone->SetFillColorAlpha(col, 0.3f);

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

    // Draw the error histograms
    bool isFirst = true;
    for (const auto &pHist : efficiencyErrorHists)
    {
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(isFirst ? "e2" : "e2 same");
        isFirst = false;
    }
    
    // Draw the histogram lines
    for (const auto &pHist : efficiencyHists)
    {
        pHist->GetYaxis()->SetRangeUser(yMin, yMax);
        pHist->Draw(isFirst ? "hist" : "hist same");
        isFirst = false;
    }
    
    PlottingHelper::SaveCanvas(pCanvas, fileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------
        
PlottingHelper::PlotStyle PlottingHelper::GetPlotStyle(const Event::Reco::Particle &particle, const std::vector<Event::Truth::Particle> &truthParticles, const bool usePoints)
{
    const auto externalType = usePoints ? ExternalPoints : External;

    if (!particle.hasMatchedMCParticle.IsSet())
        return externalType;

    if (!particle.hasMatchedMCParticle())
        return externalType;

    if (truthParticles.empty())
        return externalType;

    // Get the true pdg code applying the visibility conditions
    auto truePdgCode = -std::numeric_limits<int>::max();
    bool isGolden = false;
    try
    {
        const auto truthParticle = AnalysisHelper::GetBestMatchedTruthParticle(particle, truthParticles);
        truePdgCode = truthParticle.pdgCode();
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
        case 22:
            return usePoints ? PhotonPoints : Photon;
        default:
            return usePoints ? OtherPoints : Other;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TCanvas> PlottingHelper::GetCanvas(const unsigned int width, const unsigned int height)
{
    const auto nameStr = "canvas_" + std::to_string(++m_lastCanvasId);
    const auto name = nameStr.c_str();

    std::shared_ptr<TCanvas> pCanvas = std::make_shared<TCanvas>(name, "", width, height);

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

} // namespace ubcc1pi
