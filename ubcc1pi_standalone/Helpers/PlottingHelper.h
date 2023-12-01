/**
 *  @file  ubcc1pi_standalone/Helpers/PlottingHelper.h
 *
 *  @brief The header file for the plotting helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_PLOTTING_HELPER
#define UBCC1PI_STANDALONE_HELPERS_PLOTTING_HELPER

#include <string>
#include <memory>
#include <unordered_map>

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

#include "ubsmear.h"

#include "ubcc1pi_standalone/Interface/Event.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The plotting helper class
 */
class PlottingHelper
{
    public:
        /**
         *  @brief  The plot style enumeration
         */
        enum PlotStyle
        {
            // Basic types
            Primary,
            Secondary,
            Tertiary,
            Quaternary,
            Quinary,
            Senary,
            Septenary,
            Octonary,

            // Particle types
            Muon,
            MuonPoints,
            Proton,
            ProtonPoints,
            GoldenPion,
            GoldenPionPoints,
            NonGoldenPion,
            NonGoldenPionPoints,
            PiMinus,
            PiMinusPoints,
            Electron,
            ElectronPoints,
            Photon,
            PhotonPoints,

            // Event types
            NonFiducial,
            NumuCC0Pi,
            NumuCC1PiChargedGolden,
            NumuCC1PiChargedNonGolden,
            NumuCC1PiZero,
            NumuCCOther,
            Nue,
            NC,

            // Common types
            External,
            ExternalPoints,
            Dirt,
            DirtPoints,
            BNBData,
            Other,
            OtherPoints,
            Default,

            // Interaction Modes
            Coh,
            QE,
            MEC,
            DIS,
            Res,
            OtherInteraction,
        };

        /**
         *  @brief  The vector of all plot styles
         */
        static const std::vector<PlotStyle> AllPlotStyles;

        /**
         *  @brief  If we should draw the input style with points (or a line)
         *
         *  @param  style the input style
         *
         *  @return boolean, true if points should be used
         */
        static bool ShouldUsePoints(const PlotStyle &style);

        /**
         *  @brief  The plot class that manages multiple related histograms
         */
        class MultiPlot
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  yLabel the y-label of the histogram
                 *  @param  nBins the number of bins
                 *  @param  min the minimum value
                 *  @param  max the maximum value
                 *  @param  drawErrors whether to draw the error bands
                 *  @param  useAxisTitles whether to draw axis lables and titles
                 */
                MultiPlot(const std::string &xLabel, const std::string &yLabel, unsigned int nBins, float min, float max, bool drawErrors = true, const bool useAxisTitles = false);

                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  yLabel the y-label of the histogram
                 *  @param  binEdges the bin edges (for variable binning)
                 *  @param  drawErrors whether to draw the error bands
                 *  @param  useAxisTitles whether to draw axis lables and titles
                 */
                MultiPlot(const std::string &xLabel, const std::string &yLabel, const std::vector<float> &binEdges, bool drawErrors = true, const bool useAxisTitles = false);

                /**
                 *  @brief  If the plot is to be filled with integer values - label each bin with the integer that lands in that bin
                 */
                void SetIntegerBinLabels();

                /**
                 *  @brief  Set the bin labels
                 *
                 *  @param  labels the bin labels
                 */
                void SetBinLabels(const std::vector<std::string> &labels);

                /**
                 *  @brief  Add a cut line at a given value
                 *
                 *  @param  value the value
                 */
                void AddCutLine(const float value);

                /**
                 *  @brief  Fill the histogram with the given value and plot style
                 *
                 *  @param  value the value to fill
                 *  @param  plotStyle the style of the hisogram to fill
                 *  @param  weight the weight
                 */
                void Fill(const float value, const PlotStyle &plotStyle, const float weight = 1.f);

                /**
                 *  @brief  Draw and save the plot
                 *
                 *  @param  fileName the output file name (don't include an extension)
                 *  @param  useLogX if we should use a log scale on the X-axis
                 *  @param  scaleByBinWidth if we should scale by bin width
                 *  @param  minEntriesToDraw the minimum number of entries required before we will draw a given histogram
                 *  @param  useLogY if we should use a log scale on the Y-axis
                 */
                void SaveAs(const std::string &fileName, const bool useLogX = false, const bool scaleByBinWidth = false, const unsigned int minEntriesToDraw = 0u, const bool useLogY = false);

                /**
                 *  @brief  Draw and save the plot as a stacked histogram
                 *
                 *  @param  fileName the output file name (don't include an extension)
                 *  @param  useLogX if we should use a log scale on the X-axis
                 *  @param  scaleByBinWidth if we should scale by bin width
                 *  @param  useLogY if we should use a log scale on the Y-axis
                 *  @param  useAxisTitles whether to draw axis lables and titles
                 */
                void SaveAsStacked(const std::string &fileName, const bool useLogX = false, const bool scaleByBinWidth = false, const bool useLogY = false, const bool useAxisTitles = false);

            private:
                /**
                 *  @brief  Get clones of the histograms so they can be safely manipulated
                 *
                 *  @param  plotToHistCloneMap the output map of cloned histograms
                 */
                void GetHistogramClones(std::unordered_map<PlotStyle, TH1F* > &plotToHistCloneMap);

                /**
                 *  @brief  Normalise the histograms by their number of entries
                 *
                 *  @param  plotToHistCloneMap the histograms to normalise
                 *  @param  scaleByBinWidth if we should scale by the bin width when normalising
                 */
                void ScaleHistograms(std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap, const bool scaleByBinWidth) const;

                /**
                 *  @brief  Set the Y-range of the histograms so all fit on the canvas
                 *
                 *  @param  useLogY if we are going to use log Y-scale
                 *  @param  minEntriesToDraw the minimum number of entries required to draw a histogram (won't be considered when finding the ranges)
                 *  @param  plotToHistCloneMap the histograms to modify
                 */
                void SetHistogramYRanges(const bool useLogY, const unsigned int minEntriesToDraw, std::unordered_map<PlotStyle, TH1F*> &plotToHistCloneMap) const;

                std::string        m_xLabel;     ///< The x-label of the histogram
                unsigned int       m_nBins;      ///< The number of bins
                float              m_min;        ///< The minimum value
                float              m_max;        ///< The maximum value
                std::vector<float> m_binEdges;   ///< The edges of the bins (for variable binning)
                unsigned int       m_id;         ///< The ID of this plot plot
                unsigned int       m_cloneCount; ///< A count of the number of clones to avoid name collisions
                bool               m_drawErrors; ///< Whether to draw the error bands
                std::vector<float> m_cutValues;  ///< The values of cuts to draw

                std::unordered_map<PlotStyle, std::shared_ptr<TH1F> > m_plotToHistMap; ///< The mapping from plot style to hist

                static unsigned int m_lastId; ///< The last plotPlot ID that was set
        };

        /**
         *  @brief  Efficiency plot class
         */
        class EfficiencyPlot
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  nBins the number of bins
                 *  @param  min the minimum value
                 *  @param  max the maximum value
                 *  @param  cuts the names of all possible cuts in order
                 *  @param  drawErrors whether to draw the error bands
                 *  @param  useAxisTitles whether to draw axis lables and titles
                 */
                EfficiencyPlot(const std::string &xLabel, unsigned int nBins, float min, float max, const std::vector<string> &cuts, bool drawErrors = true, const bool useAxisTitles = false);

                /**
                 *  @brief  Constructor
                 *
                 *  @param  xLabel the x-label of the histogram
                 *  @param  binEdges the bin edges
                 *  @param  cuts the names of all possible cuts in order
                 *  @param  drawErrors whether to draw the error bands
                 *  @param  useAxisTitles whether to draw axis lables and titles
                 */
                EfficiencyPlot(const std::string &xLabel, const std::vector<float> &binEdges, const std::vector<string> &cuts, bool drawErrors = true, const bool useAxisTitles = false);

                /**
                 *  @brief  Add an event to the plot for the given cut
                 *
                 *  @param  value the value at which to add the event
                 *  @param  weight the event weight
                 *  @param  cut the cut
                 *  @param  passedCut if the cut was passed
                 */
                void AddEvent(const float value, const float weight, const std::string &cut, const bool passedCut);

                /**
                 *  @brief  Add a cut line at a given value
                 *
                 *  @param  value the value
                 */
                void AddCutLine(const float value);

                /**
                 *  @brief  Draw and save the plots for the specified cuts, and styles
                 *
                 *  @param  cuts the cuts to plots
                 *  @param  styles the styles to uses per cut
                 *  @param  fileName the output file name without an extension
                 *  @param  useAxisTitles whether to draw axis lables and titles
                 */
                void SaveAs(const std::vector<std::string> &cuts, const std::vector<PlotStyle> &styles, const std::string &fileName, const bool useAxisTitles = false);

                /**
                 *  @brief  Draw and save the plots for the specified cuts, and styles
                 *
                 *  @param  cuts the cuts to plots
                 *  @param  colors the colors to uses per cut
                 *  @param  fileName the output file name without an extension
                 *  @param  useAxisTitles whether to draw axis lables and titles
                 */
                void SaveAs(const std::vector<std::string> &cuts, const std::vector<int> &colors, const std::string &fileName, const bool useAxisTitles = false);

                /**
                 *  @brief  Draw and save the plots for all cuts
                 *
                 *  @param  fileName the file name without an extension
                 */
                void SaveAs(const std::string &fileName);

            private:


                /**
                 *  @brief  Mapping from a cut name to a pair of histograms, first being the numerator, second being the denominator
                 */
                typedef std::unordered_map<std::string, std::pair< std::shared_ptr<TH1F>, std::shared_ptr<TH1F> > > CutToPlotsMap;

                std::string              m_xLabel;        ///< The x label
                std::vector<float>       m_binEdges;      ///< The bin edges
                unsigned int             m_nBins;         ///< The number of bins
                float                    m_min;           ///< The minimum histogram value
                float                    m_max;           ///< The maximum histogram value
                std::vector<std::string> m_cuts;          ///< The cuts
                bool                     m_drawErrors;    ///< If we should draw error bands
                std::vector<float>       m_cutValues;     ///< The values of cuts to draw

                CutToPlotsMap            m_cutToPlotsMap; ///< The mapping from cut name to the numerator and denominator histograms
                unsigned int             m_id;            ///< The ID of the plot
                static unsigned int      m_lastId;        ///< The last plot ID that was set - avoids ROOT name collisions
        };

        /**
         *  @brief  Get the particle style type for the input reco particle using it's truth match information
         *
         *  @param  particle the input reco particle
         *  @param  sampleType the input sample type
         *  @param  truthParticles the input list of all truth particles
         *  @param  usePoints if we should use datapoints instead of a line
         *  @param  useAbsPdg if we should use the absolute PDG codes
         *  @param  peleeMode if we are in pelee mode
         *
         *  @return the particle style
         */
        static PlotStyle GetPlotStyle(const Event::Reco::Particle &particle, const AnalysisHelper::SampleType &sampleType, const std::vector<Event::Truth::Particle> &truthParticles, const bool usePoints = false, const bool useAbsPdg = false, const bool peleeMode = false);

        /**
         *  @brief  Get a color for a given style
         *
         *  @param  plotStyle the input style
         *
         *  @return the color
         */
        static int GetColor(const PlotStyle plotStyle);

        /**
         *  @brief  Get the complete list of unique colors used
         *
         *  @return the colors
         */
        static std::vector<int> GetColorVector();

        /**
         *  @brief  Get the plot style of an event
         *
         *  @param  sampleType the input sample type
         *  @param  pEvent the input event
         *  @param  useAbsPdg if we should use the absolute PDG codes
         *
         *  @return the plot style
         */
        static PlotStyle GetPlotStyle(const AnalysisHelper::SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const bool useAbsPdg, const bool treatNuWroLikeData = false);

        /**
         *  @brief  Get the plot style of an interaction
         *
         *  @param  interactionMode the interaction integer value
         *
         *  @return the plot style
         */
        static PlotStyle GetPlotStyle(const int interactionMode);

        /**
         *  @brief  Get the plot style of an interaction
         *
         *  @param  interactionMode the interaction integer value
         *
         *  @return the plot style
         */
        static PlotStyle GetPlotStyle2(const int interactionMode);

        /**
         *  @brief  Set the line style for a given plot type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  plotStyle the plot style
         */
        template<typename T>
        static void SetLineStyle(T *pObject, const PlotStyle plotStyle);

        /**
         *  @brief  Set the line style for a given color
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  col the color
         */
        template<typename T>
        static void SetLineStyle(T *pObject, const int col);

        /**
         *  @brief  Set the line style for a given plot type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  plotStyle the style
         */
        template<typename T>
        static void SetLineStyle(std::shared_ptr<T> &pObject, const PlotStyle plotStyle);

        /**
         *  @brief  Set the line style for a given plot type
         *
         *  @tparam T the ROOT object class
         *  @param  pObject the ROOT object
         *  @param  col the color
         */
        template<typename T>
        static void SetLineStyle(std::shared_ptr<T> &pObject, const int col);

        /**
         *  @brief  Get a canvas to draw on
         *
         *  @param  width the canvas width (px)
         *  @param  height the canvas height (px)
         *
         *  @return the canvas
         */
        static std::shared_ptr<TCanvas> GetCanvas(const unsigned int width = 960, const unsigned int height = 540);

        /**
         *  @brief  Save the canvas
         *
         *  @param  canvas the canvas to save
         *  @param  fileName the output file name (no extension)
         */
        static void SaveCanvas(std::shared_ptr<TCanvas> &canvas, const std::string &fileName);

        /**
         *  @brief  Generate uniform bin edges in the specified range
         *
         *  @param  nBins the number of bins to generate (nEdges = nBins + 1)
         *  @param  min the minimum value of the bins
         *  @param  max the maximum value of hte bins
         *
         *  @return the output bin edges
         */
        static std::vector<float> GenerateUniformBinEdges(const unsigned int nBins, const float min, const float max);

        /**
         *  @brief  Generate log bin edges in the specified range
         *
         *  @param  nBins the number of bins to generate (nEdges = nBins + 1)
         *  @param  min the minimum value of the bins
         *  @param  max the maximum value of hte bins
         *
         *  @return the output bin edges
         */
        static std::vector<float> GenerateLogBinEdges(const unsigned int nBins, const float min, const float max);

        /**
        *  @brief  Plot an input error matrix
        *
        *  @param  errorMatrix the input error matrix
        *  @param  fileName the name of the file to save
        *  @param  metadata the metadata describing the binning
        *  @param  useDefaultPalette if the default palette should be used for the z-scale, if not then a temperature palette is used instead
        *  @param  useSymmetricZRange if the z-range should be symmetric around zero
        *  @param  zRangeMax the maximum z-value to use on the scale (only used if useSymmetricZRange is true and zRangeMax is non-negative)
        */
        static void PlotErrorMatrix(const ubsmear::UBMatrix &errorMatrix, const std::string &fileName, const ubsmear::UBXSecMeta &metadata, const bool useDefaultPalette = false, const bool useSymmetricZRange = true, const float zRangeMax = -1.f);

        /**
        *  @brief  Plot an input fractional error matrix
        *
        *  @param  errorMatrix the input error matrix
        *  @param  quantityVector the input quantity to which the error matrix applied
        *  @param  fileName the name of the file to save
        *  @param  metadata the metadata describing the binning
        *  @param  useDefaultPalette if the default palette should be used for the z-scale, if not then a temperature palette is used instead
        *  @param  useSymmetricZRange if the z-range should be symmetric around zero
        */
        static void PlotFractionalErrorMatrix(const ubsmear::UBMatrix &errorMatrix, const ubsmear::UBMatrix &quantityVector, const std::string &fileName, const ubsmear::UBXSecMeta &metadata, const bool useDefaultPalette = false, const bool useSymmetricZRange = true);

    private:

        static unsigned int m_lastCanvasId;  ///< The last canvas ID
        static unsigned int m_lastPlotId;  ///< The last plot ID
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<PlottingHelper::PlotStyle> PlottingHelper::AllPlotStyles = {
    Primary,
    Secondary,
    Tertiary,
    Quaternary,
    Quinary,
    Senary,
    Septenary,
    Octonary,
    External,
    ExternalPoints,
    Dirt,
    DirtPoints,
    NonFiducial,
    Other,
    OtherPoints,
    PiMinus,
    PiMinusPoints,
    Electron,
    ElectronPoints,
    Photon,
    PhotonPoints,
    Proton,
    ProtonPoints,
    Muon,
    MuonPoints,
    NonGoldenPion,
    NonGoldenPionPoints,
    GoldenPion,
    GoldenPionPoints,
    Nue,
    NC,
    NumuCCOther,
    NumuCC0Pi,
    NumuCC1PiZero,
    NumuCC1PiChargedGolden,
    NumuCC1PiChargedNonGolden,
    BNBData,
    Default,
    Coh,
    QE,
    MEC,
    DIS,
    Res,
    OtherInteraction
};

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int PlottingHelper::MultiPlot::m_lastId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int PlottingHelper::EfficiencyPlot::m_lastId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int PlottingHelper::m_lastCanvasId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int PlottingHelper::m_lastPlotId = 0;

// -----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline void PlottingHelper::SetLineStyle(std::shared_ptr<T> &pObject, const PlotStyle plotStyle)
{
    PlottingHelper::SetLineStyle(pObject.get(), plotStyle);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline void PlottingHelper::SetLineStyle(std::shared_ptr<T> &pObject, const int col)
{
    PlottingHelper::SetLineStyle(pObject.get(), col);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int PlottingHelper::GetColor(const PlotStyle plotStyle)
{
    auto col = static_cast<int>(kBlack);

    switch (plotStyle)
    {
        case Primary:
        case BNBData:
        case Default:
            col = kBlack;
            break;

        case Secondary:
        case Muon:
        case MuonPoints:
        case NumuCC0Pi:
            col = kAzure - 2;
            break;

        case Tertiary:
        case Proton:
        case ProtonPoints:
        case NumuCC1PiZero:
        case Res:
            col = kOrange - 3;
            break;

        case Quaternary:
        case GoldenPion:
        case GoldenPionPoints:
        case NumuCC1PiChargedGolden:
        case DIS:
            col = kGreen + 1;
            break;

        case Quinary:
        case NonGoldenPion:
        case NonGoldenPionPoints:
        case NumuCC1PiChargedNonGolden:
        case MEC:
            col = kMagenta + 1;
            break;

        case Senary:
        case PiMinus:
        case PiMinusPoints:
        case NonFiducial:
        case QE:
            col = kRed - 4;
            break;

        case Septenary:
        case Electron:
        case ElectronPoints:
        case NumuCCOther:
        case Coh:
            col = kCyan + 1;
            break;

        case Octonary:
        case Photon:
        case PhotonPoints:
        case Nue:
        case OtherInteraction:
            col = kYellow + 1;
            break;

        case NC:
            col = kGreen + 3;
            break;

        case External:
        case ExternalPoints:
            col = kGray + 3;
            break;

        case Dirt:
        case DirtPoints:
            col = kMagenta - 5;
            break;

        case Other:
        case OtherPoints:
            col = kGray;
            break;

        default: break;
    }

    return col;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> PlottingHelper::GetColorVector()
{
    std::vector<int> colors;

    // Use the ordered colors first
    colors.push_back(PlottingHelper::GetColor(Primary));
    colors.push_back(PlottingHelper::GetColor(Secondary));
    colors.push_back(PlottingHelper::GetColor(Tertiary));
    colors.push_back(PlottingHelper::GetColor(Quaternary));
    colors.push_back(PlottingHelper::GetColor(Quinary));
    colors.push_back(PlottingHelper::GetColor(Senary));
    colors.push_back(PlottingHelper::GetColor(Septenary));
    colors.push_back(PlottingHelper::GetColor(Octonary));

    // Then add the rest
    for (const auto &style : PlottingHelper::AllPlotStyles)
    {
        const auto col = PlottingHelper::GetColor(style);

        if (std::find(colors.begin(), colors.end(), col) == colors.end())
            colors.push_back(col);
    }

    return colors;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline void PlottingHelper::SetLineStyle(T *pObject, const int col)
{
    pObject->SetLineWidth(2);
    pObject->SetLineColor(col);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline void PlottingHelper::SetLineStyle(T *pObject, const PlotStyle plotStyle)
{
    const auto col = PlottingHelper::GetColor(plotStyle);

    pObject->SetLineWidth(2);
    pObject->SetLineColor(col);

    if (PlottingHelper::ShouldUsePoints(plotStyle))
    {
        pObject->SetMarkerColor(col);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool PlottingHelper::ShouldUsePoints(const PlotStyle &style)
{
    if (style == ExternalPoints ||
        style == DirtPoints ||
        style == MuonPoints ||
        style == ProtonPoints ||
        style == GoldenPionPoints ||
        style == NonGoldenPionPoints ||
        style == PiMinusPoints ||
        style == ElectronPoints ||
        style == PhotonPoints ||
        style == OtherPoints)
        return true;

    return false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

} // namespace ubcc1pi

#endif
