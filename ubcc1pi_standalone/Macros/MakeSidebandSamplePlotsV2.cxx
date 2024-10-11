/**
 *  @file  ubcc1pi_standalone/Macros/MakeSidebandSamplePlotsV2.cxx
 *
 *  @brief The implementation file of the MakeSidebandSamplePlotsV2 macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include <fstream> // Todo: not use txt files

#include <TH2F.h>
#include <TMath.h>
#include <TLegend.h>
#include <THStack.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void DrawStackedHistograms(const std::map<std::string, std::shared_ptr<TH1F>> &TrueCC0pi_histograms, const std::string &histName, const std::string &postfix, const bool normalize = false)
{
    // std::cout<<"DEBUG Point X0 "<<histName<<std::endl;
    // Create a copy of the map
    std::map<std::string, std::shared_ptr<TH1F>> histogramsCopy;
    for (const auto &pair : TrueCC0pi_histograms)
    {
        histogramsCopy[pair.first] = std::shared_ptr<TH1F>(static_cast<TH1F*>(pair.second->Clone()));
    }
    // std::cout<<"DEBUG Point X0.1"<<std::endl;

    auto pCanvas = PlottingHelper::GetCanvas();
    const auto yAxisLabel = normalize ? "Area normalized number of events / bin width" : "Number of events / bin width";
    THStack *hs = new THStack(("hs: " + histName).c_str(), ("Stacked Histograms: " + histName+";;"+yAxisLabel).c_str());

    // std::cout<<"DEBUG Point X1"<<std::endl;
    if(normalize)
    {
        const auto cc1piIntegral = histogramsCopy.at("SelCC1pi_truth")->Integral();
        const auto cc0piIntegral = histogramsCopy.at("SelCC0pi_truth")->Integral();
        const double epsilon = 1e-8;
        if(cc1piIntegral<epsilon || cc0piIntegral<epsilon) throw std::runtime_error("Cannot normalize histograms as one of the integrals is zero");
        histogramsCopy.at("SelCC1pi_truth")->Sumw2();
        histogramsCopy.at("SelCC1pi_truth")->Scale(1.0/cc1piIntegral);
        histogramsCopy.at("SelCC0pi_truth_truep")->Sumw2();
        histogramsCopy.at("SelCC0pi_truth_truep")->Scale(1.0/cc0piIntegral);
        histogramsCopy.at("SelCC0pi_truth_truemu")->Sumw2();
        histogramsCopy.at("SelCC0pi_truth_truemu")->Scale(1.0/cc0piIntegral);
    }

    // std::cout<<"DEBUG Point X2"<<std::endl;
    // Get the histograms
    auto hist1 = histogramsCopy.at("SelCC0pi_truth_truep").get();
    auto hist2 = histogramsCopy.at("SelCC0pi_truth_truemu").get();
    hist1->SetFillColor(kGreen);
    hist2->SetFillColor(kGreen+2);
    // std::cout<<"DEBUG Point X2.1"<<std::endl;
    hs->Add(hist1);
    hs->Add(hist2);
    // Get the maximum value of the histograms
    double maxVal1 = hs->GetMaximum();
    double maxVal2 = histogramsCopy.at("SelCC1pi_truth")->GetMaximum();
    // std::cout<<"DEBUG Point X2.2"<<std::endl;
    // Set the y-axis range to accommodate the maximum value
    hs->SetMaximum(std::max(maxVal1, maxVal2) * 1.1);
    hs->Draw("hist");

    // std::cout<<"DEBUG Point X3"<<std::endl;
    if(!normalize)
    {
        // Create a new TH1F that is a copy of hist1
        TH1F *histSum = dynamic_cast<TH1F*>(hist1->Clone("histSum"));
        // Add hist2 to histSum
        histSum->Add(hist2);
        histSum->SetLineWidth( 3 );
        histSum->SetFillColor(kGray + 2);
        histSum->SetFillStyle(3002);
        histSum->Draw("same E2");
    }

    histogramsCopy.at("SelCC1pi_truth")->SetLineColor(kBlue);
    histogramsCopy.at("SelCC1pi_truth")->SetLineWidth( 3 );
    const auto drawStyle = normalize ? "hist same" : "E hist same";
    
    // // std::cout<<"DEBUG Point X4"<<std::endl;
    histogramsCopy.at("SelCC1pi_truth")->Draw(drawStyle);
    // TLegend *legend = new TLegend(0.4,0.2);
    // legend->SetHeader("True CC0pi"); // Add this line
    // legend->AddEntry(hist1, "Sideband selection (pion candidate is true proton)", "f");
    // legend->AddEntry(hist2, "Sideband selection (pion candidate is true muon)", "f");
    // legend->AddEntry(histogramsCopy.at("SelCC1pi_truth").get(), "Main selection", "f");
    // legend->Draw();

    // std::cout<<"DEBUG Point X5"<<std::endl;
    const auto normString = normalize ? "_areanorm" : "";

    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_"+histName+normString+"_"+postfix);
    // std::cout<<"DEBUG Point X6"<<std::endl;
}

void DrawHistograms(const std::map<std::string, std::shared_ptr<TH1F>> &TrueCC0pi_histograms, const std::string &histName, const std::string &postfix, const std::string type = "reco", const bool normalize = false)
{
    // std::cout<<"DEBUG Point E0 "<<histName<<std::endl;
    // Create a copy of the map
    std::map<std::string, std::shared_ptr<TH1F>> histogramsCopy;
    for (const auto &pair : TrueCC0pi_histograms)
    {
        histogramsCopy[pair.first] = std::shared_ptr<TH1F>(static_cast<TH1F*>(pair.second->Clone()));
    }
    // std::cout<<"DEBUG Point E0.1"<<std::endl;

    auto pCanvas = PlottingHelper::GetCanvas();
    const auto yAxisLabel = normalize ? "Area normalized number of events / bin width" : "Number of events / bin width";
    // THStack *hs = new THStack(("hs: " + histName).c_str(), ("Stacked Histograms: " + histName+";;"+yAxisLabel).c_str());

    // std::cout<<"DEBUG Point E1"<<std::endl;
    if(normalize)
    {
        const auto cc1piIntegral = histogramsCopy.at(("SelCC1pi_"+type).c_str())->Integral();
        const auto cc0piIntegral = histogramsCopy.at(("SelCC0pi_"+type).c_str())->Integral();
        const double epsilon = 1e-8;
        if(cc1piIntegral<epsilon || cc0piIntegral<epsilon) throw std::runtime_error("Cannot normalize histograms as one of the integrals is zero");
        histogramsCopy.at(("SelCC1pi_"+type).c_str())->Sumw2();
        histogramsCopy.at(("SelCC1pi_"+type).c_str())->Scale(1.0/cc1piIntegral);
        histogramsCopy.at(("SelCC0pi_"+type).c_str())->Sumw2();
        histogramsCopy.at(("SelCC0pi_"+type).c_str())->Scale(1.0/cc0piIntegral);
    }

    // std::cout<<"DEBUG Point E2"<<std::endl;
    // Get the histograms
    // auto hist1 = histogramsCopy.at("SelCC0pi_reco_truep").get();
    // auto hist2 = histogramsCopy.at("SelCC0pi_reco_truemu").get();
    // hist1->SetFillColor(kGreen);
    // hist2->SetFillColor(kGreen+2);
    // std::cout<<"DEBUG Point E2.1"<<std::endl;
    // hs->Add(hist1);
    // hs->Add(hist2);
    // Get the maximum value of the histograms
    double maxVal1 = histogramsCopy.at(("SelCC0pi_"+type).c_str())->GetMaximum();
    double maxVal2 = histogramsCopy.at(("SelCC1pi_"+type).c_str())->GetMaximum();
    // std::cout<<"DEBUG Point E2.2"<<std::endl;
    // Set the y-axis range to accommodate the maximum value
    histogramsCopy.at(("SelCC0pi_"+type).c_str())->SetMaximum(std::max(maxVal1, maxVal2) * 1.1);
    histogramsCopy.at(("SelCC0pi_"+type).c_str())->Draw("hist");

    // std::cout<<"DEBUG Point E3"<<std::endl;
    if(!normalize)
    {
        // Create a new TH1F that is a copy of hist1
        // TH1F *histSum = dynamic_cast<TH1F*>(hist1->Clone("histSum"));
        // // Add hist2 to histSum
        // histSum->Add(hist2);
        // histSum->SetLineWidth( 3 );
        histogramsCopy.at(("SelCC0pi_"+type).c_str())->SetFillColor(kGray + 2);
        histogramsCopy.at(("SelCC0pi_"+type).c_str())->SetFillStyle(3002);
        histogramsCopy.at(("SelCC0pi_"+type).c_str())->Draw("same E2");
    }

    histogramsCopy.at(("SelCC1pi_"+type).c_str())->SetLineColor(kBlue);
    histogramsCopy.at(("SelCC1pi_"+type).c_str())->SetLineWidth( 3 );
    const auto drawStyle = normalize ? "hist same" : "E hist same";
    
    // std::cout<<"DEBUG Point E4"<<std::endl;
    histogramsCopy.at(("SelCC1pi_"+type).c_str())->Draw(drawStyle);
    // TLegend *legend = new TLegend(0.4,0.2);
    // legend->SetHeader("True CC0pi"); // Add this line
    // legend->AddEntry(histogramsCopy.at(("SelCC0pi_"+type).c_str()).get(), "Sideband selection", "f");
    // legend->AddEntry(histogramsCopy.at(("SelCC1pi_"+type).c_str()).get(), "Main selection", "f");
    // legend->Draw();

    // std::cout<<"DEBUG Point E5"<<std::endl;
    const auto normString = normalize ? "_areanorm" : "";

    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_"+type+"_"+histName+normString+"_"+postfix);
    // std::cout<<"DEBUG Point E6"<<std::endl;
}

void MakeSidebandSamplePlotsV2(const Config &config)
{
    // ROOT::EnableImplicitMT(2);


    ExtractionHelper::AnalysisValueMap getValueCC1Pi;
    ExtractionHelper::AnalysisValueMap getValueCC0Pi;
    ExtractionHelper::PopulateAnalysisValueMap(getValueCC1Pi, false);
    ExtractionHelper::PopulateAnalysisValueMap(getValueCC0Pi, true);

    // A map from each cross-section to the limits of the phase-space that should be considered. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;
    for (const auto &[name, binning, scaleByBinWidth] : std::vector< std::tuple<std::string, Config::Global::Binning, bool> > {

        // The names of the cross-section kinematic parameters, and their binning information.
        // The third (boolean) parameter indicates if the cross-section bins should be scaled by their width
        { "muonCosTheta", config.global.muonCosTheta, true },
        { "muonPhi", config.global.muonPhi, true },
        { "muonMomentum", config.global.muonMomentum, true },

        { "pionCosTheta", config.global.pionCosTheta, true },
        { "pionPhi", config.global.pionPhi, true  },
        { "pionMomentum", config.global.pionMomentum, true },

        { "muonPionAngle", config.global.muonPionAngle, true },
        { "nProtons", config.global.nProtons, false }
    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));
    }

    // auto selection = SelectionHelper::GetCC0piSelectionModifiedPeLEE(-0.48f, 0.12f, true);
    const auto doNotCheckPionCandidateIsTrueProton = true; // If true, do not require that the reco pion candidate's backtracked true particle is a proton
    const auto protonBDTThresholdHigh = 0.10f;
    const auto protonBDTThresholdLow = -0.06f;
    const auto leadingProtonViaBDTResponse = true;
    const auto pProtonBDT = leadingProtonViaBDTResponse ? std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames) : nullptr;

    const std::string postfix = leadingProtonViaBDTResponse ? "ViaProtonBDT" : "ViaProtonMom";
    // auto selection = leadingProtonViaBDTResponse ? SelectionHelper::GetCC0piSelectionModifiedAgainPeLEE(-0.48f, protonBDTThresholdHigh, true) : SelectionHelper::GetCC0piSelectionModifiedPeLEE(-0.48f, protonBDTThresholdHigh, true);
    auto selection = leadingProtonViaBDTResponse ? SelectionHelper::GetCC0piSelectionModifiedAgainPeLEE(-0.4f, protonBDTThresholdHigh, true) : SelectionHelper::GetCC0piSelectionModifiedPeLEE(-0.4f, protonBDTThresholdHigh, true);

    const auto color = leadingProtonViaBDTResponse ? kRed : kGreen;

    // std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-??f, barelyResemblingProtonValue=??f\n.........................................."<<std::endl;
    // auto selection = SelectionHelper::GetCC0piSelectionModified(-0.48f, 0.12f);

    const auto cuts = selection.GetCuts();
    auto CC1piSelection = SelectionHelper::GetDefaultSelection2(true);
    const auto CC1piCuts = CC1piSelection.GetCuts();

    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    // Get the muon BDT
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

    // Get the proton BDT
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);

    // Get the golden pion BDT
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    BDTHelper::BDT goldenpionBDT("goldenPion",goldenPionFeatureNames);

    //
    // Setup the plots
    //
    const std::string yLabel = "Number of events / bin width";
    std::vector<PlottingHelper::MultiPlot> muonMomentumPlots, muonCosThetaPlots, muonPhiPlots;
    std::vector<PlottingHelper::MultiPlot> protonMomentumPlots, protonCosThetaPlots, protonPhiPlots;
    for (const auto &cut : cuts)
    {
        // muonMomentumPlots.emplace_back("Muon momentum / GeV", yLabel, config.global.muonMomentum.binEdges, true, config.global.axisTitles);
        // muonCosThetaPlots.emplace_back("Muon cos(theta)", yLabel, config.global.muonCosTheta.binEdges, true, config.global.axisTitles);
        // muonPhiPlots.emplace_back("Muon phi / rad", yLabel, config.global.muonPhi.binEdges, true, config.global.axisTitles);
        //
        // protonMomentumPlots.emplace_back("Proton momentum / GeV", yLabel, config.global.pionMomentum.binEdges, true, config.global.axisTitles);
        // protonCosThetaPlots.emplace_back("Proton cos(theta)", yLabel, config.global.pionCosTheta.binEdges, true, config.global.axisTitles);
        // protonPhiPlots.emplace_back("Proton phi / rad", yLabel, config.global.pionPhi.binEdges, true, config.global.axisTitles);

        muonMomentumPlots.emplace_back("Muon momentum / GeV", yLabel, 10, 0, 1.5, true, config.global.axisTitles);
        muonCosThetaPlots.emplace_back("Muon cos(theta)", yLabel, 10, -1, 1, true, config.global.axisTitles);
        muonPhiPlots.emplace_back("Muon phi / rad", yLabel, 15, -TMath::Pi(), TMath::Pi(), true, config.global.axisTitles);

        protonMomentumPlots.emplace_back("Leading proton momentum / GeV", yLabel, 10, 0, 1, true, config.global.axisTitles);
        protonCosThetaPlots.emplace_back("Leading proton cos(theta)", yLabel, 10, -1, 1, true, config.global.axisTitles);
        protonPhiPlots.emplace_back("Proton phi / rad", yLabel, 15, -TMath::Pi(), TMath::Pi(), true, config.global.axisTitles);
    }

    // Kinematic plots for true CC0pi1p events, in CC1pi and CC0pi selections
    std::map<std::string, std::shared_ptr<TH1F>> TrueCC0pi_ProtonMultiplicity, TrueCC0pi_MuonMomentum, TrueCC0pi_MuonCosTheta, TrueCC0pi_MuonPhi, TrueCC0pi_ProtonMomentum, TrueCC0pi_ProtonCosTheta, TrueCC0pi_ProtonPhi;

    // TrueCC0pi_MuonMomentum.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), config.global.muonMomentum.binEdges.size()-1,config.global.muonMomentum.binEdges.data()));
    // TrueCC0pi_MuonCosTheta.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), config.global.muonCosTheta.binEdges.size()-1,config.global.muonCosTheta.binEdges.data()));
    // TrueCC0pi_MuonPhi.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), config.global.muonPhi.binEdges.size()-1,config.global.muonPhi.binEdges.data()));
    // TrueCC0pi_ProtonMomentum.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), config.global.pionMomentum.binEdges.size()-1,config.global.pionMomentum.binEdges.data()));
    // TrueCC0pi_ProtonCosTheta.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), config.global.pionCosTheta.binEdges.size()-1,config.global.pionCosTheta.binEdges.data()));
    // TrueCC0pi_ProtonPhi.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), config.global.pionPhi.binEdges.size()-1,config.global.pionPhi.binEdges.data()));
    //
    // TrueCC0pi_MuonMomentum.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), config.global.muonMomentum.binEdges.size()-1,config.global.muonMomentum.binEdges.data()));
    // TrueCC0pi_MuonCosTheta.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), config.global.muonCosTheta.binEdges.size()-1,config.global.muonCosTheta.binEdges.data()));
    // TrueCC0pi_MuonPhi.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), config.global.muonPhi.binEdges.size()-1,config.global.muonPhi.binEdges.data()));
    // TrueCC0pi_ProtonMomentum.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), config.global.pionMomentum.binEdges.size()-1,config.global.pionMomentum.binEdges.data()));
    // TrueCC0pi_ProtonCosTheta.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), config.global.pionCosTheta.binEdges.size()-1,config.global.pionCosTheta.binEdges.data()));
    // TrueCC0pi_ProtonPhi.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), config.global.pionPhi.binEdges.size()-1,config.global.pionPhi.binEdges.data()));

    // ------------- Truth

    TrueCC0pi_ProtonMultiplicity.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_ProtonMultiplicity",(string(";Proton Multiplicity;")+yLabel).c_str(), 4, 0, 4));
    TrueCC0pi_MuonMomentum.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC1pi_truth", new TH1F("TrueCC0pi_SelCC1pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));

    // std::cout<<"DEBUG Point Y"<<std::endl;
    TrueCC0pi_ProtonMultiplicity.emplace("SelCC0pi_truth_truep", new TH1F("TrueCC0pi_SelCC0pi_ProtonMultiplicity_truep",(string(";Proton Multiplicity;")+yLabel).c_str(), 4, 0, 4));
    TrueCC0pi_MuonMomentum.emplace("SelCC0pi_truth_truep", new TH1F("TrueCC0pi_SelCC0pi_MuonMomentum_truep",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC0pi_truth_truep", new TH1F("TrueCC0pi_SelCC0pi_MuonCosTheta_truep",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC0pi_truth_truep", new TH1F("TrueCC0pi_SelCC0pi_MuonPhi_truep",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC0pi_truth_truep", new TH1F("TrueCC0pi_SelCC0pi_ProtonMomentum_truep",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC0pi_truth_truep", new TH1F("TrueCC0pi_SelCC0pi_ProtonCosTheta_truep",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC0pi_truth_truep", new TH1F("TrueCC0pi_SelCC0pi_ProtonPhi_truep",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));

    // std::cout<<"DEBUG Point Y1"<<std::endl;

    TrueCC0pi_ProtonMultiplicity.emplace("SelCC0pi_truth_truemu", new TH1F("TrueCC0pi_SelCC0pi_ProtonMultiplicity_truemu",(string(";Proton Multiplicity;")+yLabel).c_str(), 4, 0, 4));
    TrueCC0pi_MuonMomentum.emplace("SelCC0pi_truth_truemu", new TH1F("TrueCC0pi_SelCC0pi_MuonMomentum_truemu",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC0pi_truth_truemu", new TH1F("TrueCC0pi_SelCC0pi_MuonCosTheta_truemu",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC0pi_truth_truemu", new TH1F("TrueCC0pi_SelCC0pi_MuonPhi_truemu",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC0pi_truth_truemu", new TH1F("TrueCC0pi_SelCC0pi_ProtonMomentum_truemu",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC0pi_truth_truemu", new TH1F("TrueCC0pi_SelCC0pi_ProtonCosTheta_truemu",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC0pi_truth_truemu", new TH1F("TrueCC0pi_SelCC0pi_ProtonPhi_truemu",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    // std::cout<<"DEBUG Point Y2"<<std::endl;

    TrueCC0pi_ProtonMultiplicity.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_ProtonMultiplicity",(string(";Proton Multiplicity;")+yLabel).c_str(), 4, 0, 4));
    TrueCC0pi_MuonMomentum.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_MuonMomentum",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_MuonCosTheta",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_MuonPhi",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_ProtonMomentum",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_ProtonCosTheta",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC0pi_truth", new TH1F("TrueCC0pi_SelCC0pi_ProtonPhi",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));

    // ------------- Reco

    TrueCC0pi_ProtonMultiplicity.emplace("SelCC1pi_reco", new TH1F("TrueCC0pi_SelCC1pi_ProtonMultiplicity_reco",(string(";Proton Multiplicity;")+yLabel).c_str(), 4, 0, 4));
    TrueCC0pi_MuonMomentum.emplace("SelCC1pi_reco", new TH1F("TrueCC0pi_SelCC1pi_MuonMomentum_reco",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC1pi_reco", new TH1F("TrueCC0pi_SelCC1pi_MuonCosTheta_reco",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC1pi_reco", new TH1F("TrueCC0pi_SelCC1pi_MuonPhi_reco",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC1pi_reco", new TH1F("TrueCC0pi_SelCC1pi_ProtonMomentum_reco",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC1pi_reco", new TH1F("TrueCC0pi_SelCC1pi_ProtonCosTheta_reco",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC1pi_reco", new TH1F("TrueCC0pi_SelCC1pi_ProtonPhi_reco",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));

    TrueCC0pi_ProtonMultiplicity.emplace("SelCC0pi_reco", new TH1F("TrueCC0pi_SelCC0pi_ProtonMultiplicity_reco",(string(";Proton Multiplicity;")+yLabel).c_str(), 4, 0, 4));
    TrueCC0pi_MuonMomentum.emplace("SelCC0pi_reco", new TH1F("TrueCC0pi_SelCC0pi_MuonMomentum_reco",(string(";Muon Momentum / GeV;")+yLabel).c_str(), 10, 0, 1.5));
    TrueCC0pi_MuonCosTheta.emplace("SelCC0pi_reco", new TH1F("TrueCC0pi_SelCC0pi_MuonCosTheta_reco",(string(";Muon cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_MuonPhi.emplace("SelCC0pi_reco", new TH1F("TrueCC0pi_SelCC0pi_MuonPhi_reco",(string(";Muon phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentum.emplace("SelCC0pi_reco", new TH1F("TrueCC0pi_SelCC0pi_ProtonMomentum_reco",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5));
    TrueCC0pi_ProtonCosTheta.emplace("SelCC0pi_reco", new TH1F("TrueCC0pi_SelCC0pi_ProtonCosTheta_reco",(string(";Proton cos(theta);")+yLabel).c_str(), 10, -1, 1));
    TrueCC0pi_ProtonPhi.emplace("SelCC0pi_reco", new TH1F("TrueCC0pi_SelCC0pi_ProtonPhi_reco",(string(";Proton phi / rad;")+yLabel).c_str(), 15, -TMath::Pi(), TMath::Pi()));


    std::map<std::string, std::shared_ptr<TH2F>>  TrueCC0pi_ProtonMomentumCosTheta, TrueCC0pi_ProtonCosThetaPhi, TrueCC0pi_ProtonMomentumPhi;

    TrueCC0pi_ProtonMomentumCosTheta.emplace("SelCC1pi_truth", new TH2F("TrueCC0pi_SelCC1pi_ProtonMomentumCosTheta",";True Proton Momentum / GeV;True Proton cos(theta) / rad", 30,0,1.5,20, -1, 1));
    TrueCC0pi_ProtonCosThetaPhi.emplace("SelCC1pi_truth", new TH2F("TrueCC0pi_SelCC1pi_ProtonCosThetaPhi",";True Proton cos(theta);True Proton phi / rad", 20,-1,1.5,15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentumPhi.emplace("SelCC1pi_truth", new TH2F("TrueCC0pi_SelCC1pi_ProtonMomentumPhi",";True Proton Momentum / GeV;True Proton phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi()));

    TrueCC0pi_ProtonMomentumCosTheta.emplace("SelCC0pi_truth", new TH2F("TrueCC0pi_SelCC0pi_ProtonMomentumCosTheta",";True Proton Momentum / GeV;True Proton cos(theta) / rad", 30,0,1.5,20, -1, 1));
    TrueCC0pi_ProtonCosThetaPhi.emplace("SelCC0pi_truth", new TH2F("TrueCC0pi_SelCC0pi_ProtonCosThetaPhi",";True Proton cos(theta);True Proton phi / rad", 20,-1,1.5,15, -TMath::Pi(), TMath::Pi()));
    TrueCC0pi_ProtonMomentumPhi.emplace("SelCC0pi_truth", new TH2F("TrueCC0pi_SelCC0pi_ProtonMomentumPhi",";True Proton Momentum / GeV;True Proton phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi()));

    // Plot 2D proton kinematics vs PID variables
    // Proton BDT score
    // Muon BDT score
    // Golden pion BDT score
    // logBragg_pToMIP
    // logBragg_piToMIP
    // truncMeandEdx
    // protonForwrd likelihood
    // muonForward likelihood
    // nDescendents
    // nSpacePointsNearEnd
    // wiggliness
    // trackScore
    std::vector<TH2F*> TrueCC0pi_SelCC0pi_ProtonPIDMomentum, TrueCC0pi_SelCC0pi_ProtonPIDCosTheta, TrueCC0pi_SelCC0pi_ProtonPIDPhi;
    std::vector<TH2F*> TrueCC0pi_SelCC1pi_ProtonPIDMomentum, TrueCC0pi_SelCC1pi_ProtonPIDCosTheta, TrueCC0pi_SelCC1pi_ProtonPIDPhi;
    for (size_t i_feat=0; i_feat<featureNames.size()+3; i_feat++){
        std::string xlabel;
        std::vector<float> binning;

        if (i_feat < featureNames.size()){
            if (featureNames.at(i_feat) == "logBragg_pToMIP"){
                xlabel = "log(L_p / L_MIP)";
                binning = {30,-9,7};
            }
            if (featureNames.at(i_feat) == "logBragg_piToMIP"){
                xlabel = "log(L_pi / L_MIP)";
                binning = {30,-4,6};
            }
            if (featureNames.at(i_feat) == "truncMeandEdx"){
                xlabel = "Truncated Mean dEdx";
                binning = {20,0,10};
            }
            if (featureNames.at(i_feat) == "protonForward"){
                xlabel = "Proton forward likelihood";
                binning = {30,0.42,0.62};
            }
            if (featureNames.at(i_feat) == "muonForward"){
                xlabel = "Muon forward likelihood";
                binning = {30,0.35,0.65};
            }
            if (featureNames.at(i_feat) == "nDescendents"){
                xlabel = "Number of descendents";
                binning = {4,0,4};
            }
            if (featureNames.at(i_feat) == "nSpacePointsNearEnd"){
                xlabel = "Number of spacepoints near track end";
                binning = {45,0,90};
            }
            if (featureNames.at(i_feat) == "wiggliness"){
                xlabel = "Wiggliness";
                binning = {30,0,0.1};
            }
            if (featureNames.at(i_feat) == "trackScore"){
                xlabel = "Track Score";
                binning = {20,0,1};
            }
        }
        else if (i_feat == featureNames.size()){ // Proton BDT score
            xlabel = "Proton BDT Score";
            binning = {20,-1,1};
        }
        else if (i_feat == featureNames.size()+1){ // Muon BDT score
            xlabel = "Muon BDT Score";
            binning = {20,-1,1};
        }
        else if (i_feat == featureNames.size()+2){ // Golden pion BDT score
            xlabel = "Golden pion BDT Score";
            binning = {20,-1,1};
        }


        TrueCC0pi_SelCC0pi_ProtonPIDMomentum.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC0pi_ProtonMomentumPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Leading Proton Candidate ")+xlabel+std::string(";True Proton Momentum / GeV")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,0,1.5));
        TrueCC0pi_SelCC0pi_ProtonPIDCosTheta.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC0pi_ProtonCosThetaPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Leading Proton Candidate ")+xlabel+std::string(";True Proton cos(theta)")).c_str(),binning.at(0),binning.at(1),binning.at(2),20,-1,1));
        TrueCC0pi_SelCC0pi_ProtonPIDPhi.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC0pi_ProtonPhiPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Leading Proton Candidate ")+xlabel+std::string(";True Proton Phi / rad.")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,-TMath::Pi(),TMath::Pi()));

        TrueCC0pi_SelCC1pi_ProtonPIDMomentum.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC1pi_ProtonMomentumPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Pion Candidate ")+xlabel+std::string(";True Momentum / GeV")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,0,1.5));
        TrueCC0pi_SelCC1pi_ProtonPIDCosTheta.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC1pi_ProtonCosThetaPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Pion Candidate ")+xlabel+std::string(";True cos(theta)")).c_str(),binning.at(0),binning.at(1),binning.at(2),20,-1,1));
        TrueCC0pi_SelCC1pi_ProtonPIDPhi.emplace_back(new TH2F(std::string(std::string("TrueCC0pi_SelCC1pi_ProtonPhiPID")+std::to_string(i_feat)).c_str(),std::string(std::string(";Pion Candidate ")+xlabel+std::string(";True Phi / rad.")).c_str(),binning.at(0),binning.at(1),binning.at(2),30,-TMath::Pi(),TMath::Pi()));
    }

    std::vector<std::vector<float>> TrueCC0pi_CC1piSelected_PIDConfusionMatrix { { 0, 0 }, { 0, 0 }, { 0, 0 } };

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get list of good runs
    // -------------------------------------------------------------------------------------------------------------------------------------
    // std::vector<int> goodRuns;
    // std::ifstream file(config.global.goodRunListFile);
    // int run;
    // while (file >> run) {
    //     goodRuns.push_back(run);
    // }

    // Define a function scoped to this macro to determine if a cut was passed
    auto passedCut = [&](const std::string &cut, const std::vector<std::string> &cutsPassed) -> bool {
        return (std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end());
    };

    // Loop over the events
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        std::cout << "Reading input file: " << filePath << " normalisation: " << normalisation << std::endl;
    
        if(!useThisFile) continue;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        if(!isOverlay) continue;
        // const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        // const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        // const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        // const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        // const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);
        // const auto isMC = (sampleType != AnalysisHelper::DataBNB) && (sampleType != AnalysisHelper::DataEXT);
        const auto isMC = true;
        // if(sampleType != AnalysisHelper::Overlay && sampleType != AnalysisHelper::Dirt && sampleType != AnalysisHelper::DataBNB && sampleType != AnalysisHelper::DataEXT) continue;

        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        if (isMC) readerPeLEE.EnableSystematicBranches(); // Todo: Is this correct/optimal?
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

        // std::cout<<"\n##############\nOnly counting 10\% of event!\n##############"<<std::endl;
        for (unsigned int i = 0; i < nEvents; ++i)
        {
            // std::cout<<"DEBUG Point Y0"<<std::endl;
            readerPeLEE.LoadEvent(i);

            // std::cout<<"DEBUG Point Y0.1"<<std::endl;

            // const auto run = pEventPeLEE->metadata.run();
            // const auto isGoodRun = (isDataBNB || isDataEXT) ? std::find(goodRuns.begin(), goodRuns.end(), run) != goodRuns.end() : true; // Apply good runs cuts to data
            // // if(!isGoodRun) continue;
            // if(!isGoodRun)
            // {
            //     continue;
            // }

            // std::cout<<"DEBUG Point Y1"<<std::endl;
            Event event(*pEventPeLEE, false);// true or false decides wether to cut generation!=2 particles
            // std::cout<<"DEBUG Point Y1.1"<<std::endl;
            const auto pEvent = std::make_shared<Event>(event);
            // std::cout << "DEBUG Point Y1.2" << std::endl;

            // Fill kinematic comparison and efficiency plots for true CC0pi events that pass the full selection
            // First, check if event is true CC0pi
            const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
            if (isTrueCC0Pi)
            {
                AnalysisHelper::PrintLoadingBar(i, nEvents);
                const auto recoParticles = pEvent->reco.particles;
                // std::cout << "DEBUG Point Y1.3" << std::endl;

                // Run the event selection and store which cuts are passed
                // std::cout << "DEBUG Point Y1.4" << std::endl;
                const auto &[isSelectedCC0pi, cutsPassedCC0pi, assignedPdgCodes] = selection.Execute(pEvent);
                const auto isSelectedCC0piGeneric = SelectionHelper::IsCutPassed(cutsPassedCC0pi, config.global.lastCutGeneric);

                // std::cout << "DEBUG Point Y1.5" << std::endl;
                if(!passedCut("topologicalScoreCC", cutsPassedCC0pi)) continue; // replacement for passesCCInclusive
                // std::cout << "DEBUG Point Y1.6" << std::endl;

                // Get the truth and reco analysis data
                const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

                // Get the proton reco data
                float protonMomentum = -std::numeric_limits<float>::max();
                float protonCosTheta = -std::numeric_limits<float>::max();
                float protonPhi = -std::numeric_limits<float>::max();
                float pbdtResponse = -std::numeric_limits<float>::max();
                float mbdtResponse = -std::numeric_limits<float>::max();
                float gbdtResponse = -std::numeric_limits<float>::max();
                std::vector<float> allfeatures;
                std::vector<float> pfeatures;
                std::vector<float> mfeatures;
                std::vector<float> gfeatures;

                // std::cout<<"DEBUG Point Y2"<<std::endl;

                const auto protonIndex = leadingProtonViaBDTResponse ? SelectionHelper::GetLeadingProtonCandidateIndexByBDTScore(recoParticles, assignedPdgCodes, protonBDTThresholdHigh, protonBDTThresholdLow, pProtonBDT) : SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, assignedPdgCodes);
                // std::cout<<"DEBUG Point Y2.01"<<std::endl;

                if (protonIndex!=std::numeric_limits<unsigned int>::max()){
                    // std::cout<<"DEBUG Point Y2.02"<<std::endl;
                    const auto proton = recoParticles.at(protonIndex);
                    const auto protonDir = TVector3(proton.directionX(), proton.directionY(), proton.directionZ()).Unit();
                    // std::cout<<"DEBUG Point Y2.03"<<std::endl;
                    protonCosTheta = protonDir.Z();
                    protonPhi = std::atan2(protonDir.Y(), protonDir.X());
                    protonMomentum = AnalysisHelper::GetProtonMomentumFromRange(proton.range());
                    // std::cout<<"DEBUG Point Y2.04"<<std::endl;
                    const auto hasallFeatures = BDTHelper::GetBDTFeatures(proton, featureNames, allfeatures);

                    const auto haspFeatures = BDTHelper::GetBDTFeatures(proton, protonFeatureNames, pfeatures);
                    // std::cout<<"DEBUG Point Y2.05"<<std::endl;

                    if (haspFeatures){
                        // std::cout<<"DEBUG Point Y2.1"<<std::endl;
                        pbdtResponse = protonBDT.GetResponse(pfeatures);
                        // std::cout<<"DEBUG Point Y2.2"<<std::endl;
                    }

                    const auto hasmFeatures = BDTHelper::GetBDTFeatures(proton, muonFeatureNames, mfeatures);
                    if (hasmFeatures){
                        // std::cout<<"DEBUG Point Y2.3"<<std::endl;
                        mbdtResponse = muonBDT.GetResponse(mfeatures);
                        // std::cout<<"DEBUG Point Y2.4"<<std::endl;
                    }

                    const auto hasgFeatures = BDTHelper::GetBDTFeatures(proton, goldenPionFeatureNames, gfeatures);
                    if (hasgFeatures){
                        // std::cout<<"DEBUG Point Y2.5"<<std::endl;
                        gbdtResponse = goldenpionBDT.GetResponse(gfeatures);
                        // std::cout<<"DEBUG Point Y2.6"<<std::endl;
                    }
                    // std::cout<<"DEBUG Point Y2.7"<<std::endl;
                }

                // Get the features we want to plot (truth, not reco)
                const auto truemuonidx = AnalysisHelper::GetTrueMuonIndex(pEvent->truth, config.global.useAbsPdg);
                const auto truemu = pEvent->truth.particles.at(truemuonidx);
                const auto TrueMuDir = TVector3(truemu.momentumX(), truemu.momentumY(), truemu.momentumZ()).Unit();
                const auto TrueMuCosTheta = TrueMuDir.Z();
                const auto TrueMuPhi = std::atan2(TrueMuDir.Y(),TrueMuDir.X());

                // Todo: Remove this 
                const auto truthDataCC0Pi = AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
                const auto truthCC0PiMuonMomentum = truthDataCC0Pi.muonMomentum;
                const auto truthCC0PiMuonCosTheta = truthDataCC0Pi.muonCosTheta;
                const auto truthCC0PiMuonPhi = truthDataCC0Pi.muonPhi;
                const auto smallValue = 1e-4f;
                if(std::abs(truemu.momentum() - truthCC0PiMuonMomentum) > smallValue)
                    throw std::logic_error("truemu.momentum() != truthCC0PiMuonMomentum");
                if(std::abs(TrueMuCosTheta - truthCC0PiMuonCosTheta) > smallValue)
                    throw std::logic_error("TrueMuCosTheta != truthCC0PiMuonCosTheta");
                if(std::abs(TrueMuPhi - truthCC0PiMuonPhi) > smallValue)
                    throw std::logic_error("TrueMuPhi != truthCC0PiMuonPhi");

                // If event passes CC0pi selection, fill plots
                if (isSelectedCC0piGeneric){

                    const auto recoDataCC0Pi = AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, isSelectedCC0pi);

                    // auto passesPhaseSpaceRecoCC0Pi = false;
                    // if (isSelectedCC0piGeneric)
                    // {
                    // Start by assuming the event passes the phase-space cuts
                    auto passesPhaseSpaceRecoCC0Pi = true;

                    // Check the value of the kinematic quantities are within the phase-space limits
                    for (const auto &[name, minMax] : phaseSpaceMap)
                    {
                        const auto &[min, max] = minMax;
                        const auto value = getValueCC0Pi.at(name)(recoDataCC0Pi);

                        std::cout<<"CC0pi phase-space - name: " << name << " - passes: " << ((value < min || value > max) ? "false" : "true") << " " << min << "<" << value << "<" << max << std::endl;

                        if (value < min || value > max)
                        {
                            passesPhaseSpaceRecoCC0Pi = false;
                            break;
                        }
                    }
                    
                    if(!passesPhaseSpaceRecoCC0Pi)
                    {
                        std::cout<<"Skipping event that passed cc0pi generic selection but did not pass phase space cuts"<<std::endl;
                        continue;
                    }

                    TrueCC0pi_ProtonMultiplicity.at("SelCC0pi_reco")->Fill(recoDataCC0Pi.nProtons-1, weight); // Do not count the pion candidate as a proton here
                    TrueCC0pi_MuonMomentum.at("SelCC0pi_reco")->Fill(recoDataCC0Pi.muonMomentum, weight);
                    TrueCC0pi_MuonCosTheta.at("SelCC0pi_reco")->Fill(recoDataCC0Pi.muonCosTheta, weight);
                    TrueCC0pi_MuonPhi.at("SelCC0pi_reco")->Fill(recoDataCC0Pi.muonPhi, weight);
                    TrueCC0pi_ProtonMomentum.at("SelCC0pi_reco")->Fill(recoDataCC0Pi.protonMomentum, weight);
                    TrueCC0pi_ProtonCosTheta.at("SelCC0pi_reco")->Fill(recoDataCC0Pi.protonCosTheta, weight);
                    TrueCC0pi_ProtonPhi.at("SelCC0pi_reco")->Fill(recoDataCC0Pi.protonPhi, weight);

                    // For proton kinematics, use true leading proton. If no true leading proton is found, don't fill any kinematics
                    // const auto leadpidx = AnalysisHelper::GetTrueLeadingProtonIndex(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
                    // std::cout<<"DEBUG Point Y5"<<std::endl;

                    const auto leadingRecoProton = recoParticles.at(protonIndex);

                    unsigned int leadpidx = std::numeric_limits<unsigned int>::max(); // todo check if this is the right way to do this !!!!!!!!!!!!!!!!!
                    try{
                        leadpidx = AnalysisHelper::GetBestMatchedTruthParticleIndex(leadingRecoProton,pEvent->truth.particles);
                    }
                    catch(const std::logic_error &){
                        std::cout<<"try failed for leadpidx\n"<<std::endl;
                        // throw std::logic_error("try failed for leadpidx");
                        // continue;
                    }

                    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
                    // const auto nProtons = AnalysisHelper::CountParticlesAboveMomentumThreshold(visibleParticles, 2212, config.global.useAbsPdg, config.global.protonMomentumThreshold);
                    const auto nProtons = std::min(AnalysisHelper::CountParticlesAboveMomentumThreshold(visibleParticles, 2212, config.global.useAbsPdg, config.global.protonMomentumThreshold)-1, 2u); //-1 as one proton is treated as the pion
                    // std::cout<<"DEBUG Point Y5.1"<<std::endl;

                    TrueCC0pi_ProtonMultiplicity.at("SelCC0pi_truth")->Fill(nProtons, weight);
                    TrueCC0pi_MuonMomentum.at("SelCC0pi_truth")->Fill(truemu.momentum(), weight);
                    TrueCC0pi_MuonCosTheta.at("SelCC0pi_truth")->Fill(TrueMuCosTheta, weight);
                    TrueCC0pi_MuonPhi.at("SelCC0pi_truth")->Fill(TrueMuPhi, weight);

                    if(leadpidx!=std::numeric_limits<unsigned int>::max())
                    {
                        const auto trueleadp = pEvent->truth.particles.at(leadpidx);
                        const auto dir = TVector3(trueleadp.momentumX(), trueleadp.momentumY(), trueleadp.momentumZ()).Unit();
                        const auto cosTheta = dir.Z();
                        const auto phi = std::atan2(dir.Y(),dir.X());
                        if (doNotCheckPionCandidateIsTrueProton || trueleadp.pdgCode() == 2212)
                        {
                            const auto category = trueleadp.pdgCode() == 2212 ? "SelCC0pi_truth_truep" : "SelCC0pi_truth_truemu";
                            TrueCC0pi_ProtonMultiplicity.at(category)->Fill(nProtons, weight);
                            TrueCC0pi_MuonMomentum.at(category)->Fill(truemu.momentum(), weight);
                            TrueCC0pi_MuonCosTheta.at(category)->Fill(TrueMuCosTheta, weight);
                            TrueCC0pi_MuonPhi.at(category)->Fill(TrueMuPhi, weight);

                            TrueCC0pi_ProtonMomentum.at(category)->Fill(trueleadp.momentum(), weight);
                            TrueCC0pi_ProtonCosTheta.at(category)->Fill(cosTheta, weight);
                            TrueCC0pi_ProtonPhi.at(category)->Fill(phi, weight);

                            TrueCC0pi_ProtonMomentum.at("SelCC0pi_truth")->Fill(trueleadp.momentum(), weight);
                            TrueCC0pi_ProtonCosTheta.at("SelCC0pi_truth")->Fill(cosTheta, weight);
                            TrueCC0pi_ProtonPhi.at("SelCC0pi_truth")->Fill(phi, weight);

                            TrueCC0pi_ProtonMomentumCosTheta.at("SelCC0pi_truth")->Fill(trueleadp.momentum(), cosTheta, weight);
                            TrueCC0pi_ProtonCosThetaPhi.at("SelCC0pi_truth")->Fill(cosTheta, phi, weight);
                            TrueCC0pi_ProtonMomentumPhi.at("SelCC0pi_truth")->Fill(trueleadp.momentum(), phi, weight);
                        }

                        // std::cout<<"DEBUG Point Y5.2"<<std::endl;

                        for (size_t i_feat=0; i_feat<allfeatures.size()+3;i_feat++){
                            float fillx = -std::numeric_limits<float>::max();

                            if (i_feat < allfeatures.size())
                                fillx = allfeatures.at(i_feat);
                            else if (i_feat == allfeatures.size()) // Proton BDT score
                                fillx = pbdtResponse;
                            else if (i_feat == allfeatures.size()+1) // Muon BDT score
                                fillx = mbdtResponse;
                            else if (i_feat == allfeatures.size()+2) // Golden pion BDT score
                                fillx = gbdtResponse;

                            TrueCC0pi_SelCC0pi_ProtonPIDMomentum.at(i_feat)->Fill(fillx, trueleadp.momentum(), weight);
                            TrueCC0pi_SelCC0pi_ProtonPIDCosTheta.at(i_feat)->Fill(fillx, cosTheta, weight);
                            TrueCC0pi_SelCC0pi_ProtonPIDPhi.at(i_feat)->Fill(fillx, phi, weight);
                        }
                        // }
                    }
                }
                // std::cout<<"DEBUG Point Y6"<<std::endl;

                // Find out if the event passes the generic CC1pi selection
                const auto &[isSelectedCC1pi, cutsPassedCC1pi, assignedPdgCodesCC1pi] = CC1piSelection.Execute(pEvent);
                const auto &isSelectedCC1piGeneric = SelectionHelper::IsCutPassed(cutsPassedCC1pi, config.global.lastCutGeneric);
                // std::cout<<"DEBUG Point Y6.1"<<std::endl;

                // If event passes CC1pi+ selection, fill plots
                if (isSelectedCC1piGeneric){
                    // For proton kinematics, use proton that is selected as the pi+ candidate. If the pi+ candidate is a true proton, fill muon and proton kinematic plots
                    const auto recoPion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1pi, 211));
                    const auto recoMuon = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1pi, 13));

                    const auto recoDataCC1Pi = AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodesCC1pi, isSelectedCC1pi);

                    // auto passesPhaseSpaceRecoCC1Pi = false;
                    // Start by assuming the event passes the phase-space cuts
                    auto passesPhaseSpaceRecoCC1Pi = true;

                    // Check the value of the kinematic quantities are within the phase-space limits
                    for (const auto &[name, minMax] : phaseSpaceMap)
                    {
                        const auto &[min, max] = minMax;
                        const auto value = getValueCC1Pi.at(name)(recoDataCC1Pi);

                        std::cout<<"CC1pi phase-space - name: " << name << " - passes: " << ((value < min || value > max) ? "false" : "true") << " " << min << "<" << value << "<" << max << std::endl;
                        if (value < min || value > max)
                        {
                            passesPhaseSpaceRecoCC1Pi = false;
                            break;
                        }
                    }

                    if(!passesPhaseSpaceRecoCC1Pi)
                    {
                        std::cout<<"Skipping event that passed cc1pi generic selection but did not pass phase space cuts"<<std::endl;
                        continue;
                    }

                    TrueCC0pi_ProtonMultiplicity.at("SelCC1pi_reco")->Fill(recoDataCC1Pi.nProtons, weight);
                    TrueCC0pi_MuonMomentum.at("SelCC1pi_reco")->Fill(recoDataCC1Pi.muonMomentum, weight);
                    TrueCC0pi_MuonCosTheta.at("SelCC1pi_reco")->Fill(recoDataCC1Pi.muonCosTheta, weight);
                    TrueCC0pi_MuonPhi.at("SelCC1pi_reco")->Fill(recoDataCC1Pi.muonPhi, weight);
                    TrueCC0pi_ProtonMomentum.at("SelCC1pi_reco")->Fill(recoDataCC1Pi.pionMomentum, weight);
                    TrueCC0pi_ProtonCosTheta.at("SelCC1pi_reco")->Fill(recoDataCC1Pi.pionCosTheta, weight);
                    TrueCC0pi_ProtonPhi.at("SelCC1pi_reco")->Fill(recoDataCC1Pi.pionPhi, weight);

                    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
                    // const auto nProtons = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 2212, config.global.useAbsPdg);
                    const auto nProtons = std::min(AnalysisHelper::CountParticlesAboveMomentumThreshold(visibleParticles, 2212, config.global.useAbsPdg, config.global.protonMomentumThreshold)-1, 2u); //-1 as one proton is treated as the pion
                    TrueCC0pi_ProtonMultiplicity.at("SelCC1pi_truth")->Fill(nProtons, weight);
                    TrueCC0pi_MuonMomentum.at("SelCC1pi_truth")->Fill(truemu.momentum(), weight);
                    TrueCC0pi_MuonCosTheta.at("SelCC1pi_truth")->Fill(TrueMuCosTheta, weight);
                    TrueCC0pi_MuonPhi.at("SelCC1pi_truth")->Fill(TrueMuPhi, weight);

                    Event::Truth::Particle pionMatch;
                    Event::Truth::Particle muonMatch;
                    Event::Truth::Particle leadingProtonMatch;

                    try{
                        pionMatch = AnalysisHelper::GetBestMatchedTruthParticle(recoPion,pEvent->truth.particles);
                    }
                    catch(const std::logic_error &){
                        std::cout<<"try failed for pionMatch\n"<<std::endl;
                        continue;
                    }
                    switch(pionMatch.pdgCode())
                    {
                        case 13:
                            TrueCC0pi_CC1piSelected_PIDConfusionMatrix[0][0]+=weight;
                            break;
                        case 2212:
                            TrueCC0pi_CC1piSelected_PIDConfusionMatrix[0][1]+=weight;
                    }

                    try{
                        muonMatch = AnalysisHelper::GetBestMatchedTruthParticle(recoMuon,pEvent->truth.particles);
                        switch(muonMatch.pdgCode())
                        {
                            case 13:
                                TrueCC0pi_CC1piSelected_PIDConfusionMatrix[1][0]+=weight;
                                break;
                            case 2212:
                                TrueCC0pi_CC1piSelected_PIDConfusionMatrix[1][1]+=weight;
                        }
                    }
                    catch(const std::logic_error &){
                        std::cout<<"try failed for muonMatch\n"<<std::endl;
                    }

                    if (protonIndex!=std::numeric_limits<unsigned int>::max())
                    {
                        const auto recoLeadingPion = recoParticles.at(protonIndex);
                        Event::Truth::Particle protonMatch;
                        try{
                            leadingProtonMatch = AnalysisHelper::GetBestMatchedTruthParticle(recoLeadingPion,pEvent->truth.particles);
                            switch(leadingProtonMatch.pdgCode())
                            {
                                case 13:
                                    TrueCC0pi_CC1piSelected_PIDConfusionMatrix[2][0]+=weight;
                                    break;
                                case 2212:
                                    TrueCC0pi_CC1piSelected_PIDConfusionMatrix[2][1]+=weight;
                            }
                        }
                        catch(const std::logic_error &){
                            std::cout<<"try failed for leadingProtonMatch\n"<<std::endl;
                            // continue;
                        }
                    }

                    // std::cout<<"DEBUG Point Y7"<<std::endl;
                    if (doNotCheckPionCandidateIsTrueProton || pionMatch.pdgCode()==2212)
                    {
                        const auto dir = TVector3(pionMatch.momentumX(), pionMatch.momentumY(), pionMatch.momentumZ()).Unit();
                        const auto cosTheta = dir.Z();
                        const auto phi = std::atan2(dir.Y(),dir.X());

                        TrueCC0pi_ProtonMomentum.at("SelCC1pi_truth")->Fill(pionMatch.momentum(), weight);
                        TrueCC0pi_ProtonCosTheta.at("SelCC1pi_truth")->Fill(cosTheta, weight);
                        TrueCC0pi_ProtonPhi.at("SelCC1pi_truth")->Fill(phi, weight);

                        TrueCC0pi_ProtonMomentumCosTheta.at("SelCC1pi_truth")->Fill(pionMatch.momentum(),cosTheta, weight);
                        TrueCC0pi_ProtonCosThetaPhi.at("SelCC1pi_truth")->Fill(cosTheta,phi, weight);
                        TrueCC0pi_ProtonMomentumPhi.at("SelCC1pi_truth")->Fill(pionMatch.momentum(),phi, weight);

                        std::vector<float> allfeatures;
                        std::vector<float> pfeatures;
                        std::vector<float> mfeatures;
                        std::vector<float> gfeatures;

                        const auto hasallFeatures = BDTHelper::GetBDTFeatures(recoPion, featureNames, allfeatures);

                        const auto haspFeatures = BDTHelper::GetBDTFeatures(recoPion, protonFeatureNames, pfeatures);
                        if (haspFeatures){
                            // std::cout<<"DEBUG Point Y7.1"<<std::endl;
                            pbdtResponse = protonBDT.GetResponse(pfeatures);
                            // std::cout<<"DEBUG Point Y7.2"<<std::endl;
                        }

                        const auto hasmFeatures = BDTHelper::GetBDTFeatures(recoPion, muonFeatureNames, mfeatures);
                        if (hasmFeatures){
                            // std::cout<<"DEBUG Point Y7.3"<<std::endl;
                            mbdtResponse = muonBDT.GetResponse(mfeatures);
                            // std::cout<<"DEBUG Point Y7.4"<<std::endl;
                        }

                        const auto hasgFeatures = BDTHelper::GetBDTFeatures(recoPion, goldenPionFeatureNames, gfeatures);
                        if (hasgFeatures){
                            // std::cout<<"DEBUG Point Y7.5"<<std::endl;
                            gbdtResponse = goldenpionBDT.GetResponse(gfeatures);
                            // std::cout<<"DEBUG Point Y7.6"<<std::endl;
                        }

                        for (size_t i_feat=0; i_feat<allfeatures.size()+3;i_feat++){
                            float fillx = -std::numeric_limits<float>::max();

                            if (i_feat < allfeatures.size())
                                 fillx = allfeatures.at(i_feat);
                            else if (i_feat == allfeatures.size()) // Proton BDT score
                                fillx = pbdtResponse;
                            else if (i_feat == allfeatures.size()+1) // Muon BDT score
                                fillx = mbdtResponse;
                            else if (i_feat == allfeatures.size()+2) // Golden pion BDT score
                                fillx = gbdtResponse;

                            TrueCC0pi_SelCC1pi_ProtonPIDMomentum.at(i_feat)->Fill(fillx,pionMatch.momentum(), weight);
                            TrueCC0pi_SelCC1pi_ProtonPIDCosTheta.at(i_feat)->Fill(fillx,cosTheta, weight);
                            TrueCC0pi_SelCC1pi_ProtonPIDPhi.at(i_feat)->Fill(fillx,phi, weight);
                        }
                    }
                    // std::cout<<"DEBUG Point Y8"<<std::endl;
                }
            } // end of true numu CC0pi

        }
    }

    // std::cout<<"DEBUG Point Y9"<<std::endl;
    for (unsigned int iCut = 0; iCut < cuts.size(); ++iCut)
    {
        const auto &cut = cuts.at(iCut);

        const std::string suffix = std::to_string(iCut) + "_" + cut + "_" + postfix;
        muonMomentumPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonMomentum_" + suffix, false, config.global.scaleByBinWidth, false, config.global.axisTitles);
        muonCosThetaPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonCosTheta_" + suffix, false, config.global.scaleByBinWidth, false, config.global.axisTitles);
        muonPhiPlots.at(iCut).SaveAsStacked("reco_cc1pi_muonPhi_" + suffix, false, config.global.scaleByBinWidth, false, config.global.axisTitles);

        protonMomentumPlots.at(iCut).SaveAsStacked("reco_cc1pi_protonMomentum_" + suffix, false, config.global.scaleByBinWidth, false, config.global.axisTitles);
        protonCosThetaPlots.at(iCut).SaveAsStacked("reco_cc1pi_protonCosTheta_" + suffix, false, config.global.scaleByBinWidth, false, config.global.axisTitles);
        protonPhiPlots.at(iCut).SaveAsStacked("reco_cc1pi_protonPhi_" + suffix, false, config.global.scaleByBinWidth, false, config.global.axisTitles);
    }

    // std::cout<<"DEBUG Point Y10"<<std::endl;
    // Save kinematic comparison plots for true CC0pi


    DrawStackedHistograms(TrueCC0pi_ProtonMultiplicity, "ProtonMultiplicity", postfix);
    DrawStackedHistograms(TrueCC0pi_MuonMomentum, "MuonMomentum", postfix);
    DrawStackedHistograms(TrueCC0pi_MuonCosTheta, "MuonCosTheta", postfix);
    DrawStackedHistograms(TrueCC0pi_MuonPhi, "MuonPhi", postfix);
    DrawStackedHistograms(TrueCC0pi_ProtonMomentum, "ProtonMomentum", postfix);
    DrawStackedHistograms(TrueCC0pi_ProtonCosTheta, "ProtonCosTheta", postfix);
    DrawStackedHistograms(TrueCC0pi_ProtonPhi, "ProtonPhi", postfix);

    DrawStackedHistograms(TrueCC0pi_ProtonMultiplicity, "ProtonMultiplicity", postfix, true);
    DrawStackedHistograms(TrueCC0pi_MuonMomentum, "MuonMomentum", postfix, true);
    DrawStackedHistograms(TrueCC0pi_MuonCosTheta, "MuonCosTheta", postfix, true);
    DrawStackedHistograms(TrueCC0pi_MuonPhi, "MuonPhi", postfix, true);
    DrawStackedHistograms(TrueCC0pi_ProtonMomentum, "ProtonMomentum", postfix, true);
    DrawStackedHistograms(TrueCC0pi_ProtonCosTheta, "ProtonCosTheta", postfix, true);
    DrawStackedHistograms(TrueCC0pi_ProtonPhi, "ProtonPhi", postfix, true);

    const auto norm = false;

    DrawHistograms(TrueCC0pi_ProtonMultiplicity, "ProtonMultiplicity", postfix, "reco", norm);
    DrawHistograms(TrueCC0pi_MuonMomentum, "MuonMomentum", postfix, "reco", norm);
    DrawHistograms(TrueCC0pi_MuonCosTheta, "MuonCosTheta", postfix, "reco", norm);
    DrawHistograms(TrueCC0pi_MuonPhi, "MuonPhi", postfix, "reco", norm);
    DrawHistograms(TrueCC0pi_ProtonMomentum, "ProtonMomentum", postfix, "reco", norm);
    DrawHistograms(TrueCC0pi_ProtonCosTheta, "ProtonCosTheta", postfix, "reco", norm);
    DrawHistograms(TrueCC0pi_ProtonPhi, "ProtonPhi", postfix, "reco", norm);

    DrawHistograms(TrueCC0pi_ProtonMultiplicity, "ProtonMultiplicity", postfix, "truth", norm);
    DrawHistograms(TrueCC0pi_MuonMomentum, "MuonMomentum", postfix, "truth", norm);
    DrawHistograms(TrueCC0pi_MuonCosTheta, "MuonCosTheta", postfix, "truth", norm);
    DrawHistograms(TrueCC0pi_MuonPhi, "MuonPhi", postfix, "truth", norm);
    DrawHistograms(TrueCC0pi_ProtonMomentum, "ProtonMomentum", postfix, "truth", norm);
    DrawHistograms(TrueCC0pi_ProtonCosTheta, "ProtonCosTheta", postfix, "truth", norm);
    DrawHistograms(TrueCC0pi_ProtonPhi, "ProtonPhi", postfix, "truth", norm);


    // // Now normalise plots so we can compare shapes
    auto pCanvas = PlottingHelper::GetCanvas();
    TrueCC0pi_ProtonMomentumCosTheta.at("SelCC0pi_truth")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumCosTheta_CC0pi_"+postfix);
    TrueCC0pi_ProtonCosThetaPhi.at("SelCC0pi_truth")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonCosThetaPhi_CC0pi_"+postfix);
    TrueCC0pi_ProtonMomentumPhi.at("SelCC0pi_truth")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumPhi_CC0pi_"+postfix);

    TrueCC0pi_ProtonMomentumCosTheta.at("SelCC1pi_truth")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumCosTheta_CC1pi_"+postfix);
    TrueCC0pi_ProtonCosThetaPhi.at("SelCC1pi_truth")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonCosThetaPhi_CC1pi_"+postfix);
    TrueCC0pi_ProtonMomentumPhi.at("SelCC1pi_truth")->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonMomentumPhi_CC1pi_"+postfix);

    pCanvas->SetLogz();
    // TrueCC0pi_ProtonBDTscoreMomentum.at("SelCC0pi_truth")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonMomentum_CC0pi");
    // TrueCC0pi_ProtonBDTscoreCosTheta.at("SelCC0pi_truth")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonCosTheta_CC0pi");
    // TrueCC0pi_ProtonBDTscorePhi.at("SelCC0pi_truth")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonPhi_CC0pi");
    //
    // TrueCC0pi_ProtonBDTscoreMomentum.at("SelCC1pi_truth")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonMomentum_CC1pi");
    // TrueCC0pi_ProtonBDTscoreCosTheta.at("SelCC1pi_truth")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonCosTheta_CC1pi");
    // TrueCC0pi_ProtonBDTscorePhi.at("SelCC1pi_truth")->Draw("colz");
    // PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_BDTScoreProtonPhi_CC1pi");
    // std::cout<<"DEBUG Point Y11"<<std::endl;
    for (unsigned int i_hist=0; i_hist<TrueCC0pi_SelCC0pi_ProtonPIDMomentum.size(); i_hist++){
        TrueCC0pi_SelCC0pi_ProtonPIDMomentum.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDMomentum_"+std::to_string(i_hist)+"_CC0pi_"+postfix);
        TrueCC0pi_SelCC0pi_ProtonPIDCosTheta.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDCosTheta_"+std::to_string(i_hist)+"_CC0pi_"+postfix);
        TrueCC0pi_SelCC0pi_ProtonPIDPhi.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDPhi_"+std::to_string(i_hist)+"_CC0pi_"+postfix);

        TrueCC0pi_SelCC1pi_ProtonPIDMomentum.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDMomentum_"+std::to_string(i_hist)+"_CC1pi_"+postfix);
        TrueCC0pi_SelCC1pi_ProtonPIDCosTheta.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDCosTheta_"+std::to_string(i_hist)+"_CC1pi_"+postfix);
        TrueCC0pi_SelCC1pi_ProtonPIDPhi.at(i_hist)->Draw("colz");
        PlottingHelper::SaveCanvas(pCanvas,"SidebandComparisons_ProtonPIDPhi_"+std::to_string(i_hist)+"_CC1pi_"+postfix);
    }

    // std::cout<<"DEBUG Point Y12"<<std::endl;
    std::ofstream outfile;
    outfile.open("TrueCC0pi_CC1piSelected_PIDConfusionMatrix.txt");
    outfile<<"Format: recoPion recoMuon recoLeadingProton\ntrueMuon \n trueProton\nValues:\n";
    for(const auto row : TrueCC0pi_CC1piSelected_PIDConfusionMatrix)
    {
        for (const auto elem: row)
        {
            outfile<<elem<<" ";
        }
        outfile<<"\n";
    }
    outfile.close();

    // std::cout<<"DEBUG Point Y13"<<std::endl;
    std::cout<<"---------All Done---------"<<std::endl;
}

} // namespace ubcc1pi_macros
