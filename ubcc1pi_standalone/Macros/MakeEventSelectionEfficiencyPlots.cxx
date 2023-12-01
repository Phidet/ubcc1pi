/**
 *  @file  ubcc1pi_standalone/Macros/MakeEventSelectionEfficiencyPlots.cxx
 *
 *  @brief The implementation file of the MakeEventSelectionEfficiencyPlots macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void MakeEventSelectionEfficiencyPlots(const Config &config)
{
    // Get the selection
    auto selection = SelectionHelper::GetDefaultSelection(true);
    const auto allCuts = selection.GetCuts();

    std::cout << "Making plots for cuts:" << std::endl;
    for (const auto &cut : allCuts)
        std::cout << " - " << cut << std::endl;

    // Set up the plots
    const auto drawErrors = config.efficiencyPlots.drawErrors;

    // True efficiencies
    auto plot_nuEnergy = PlottingHelper::EfficiencyPlot("Neutrino energy / GeV", 40u, 0.4f, 2.5f, allCuts, drawErrors);
    auto plot_nProtons = PlottingHelper::EfficiencyPlot("Proton multiplicity", 5u, 0, 5, allCuts, drawErrors);
    auto plot_muMomentum = PlottingHelper::EfficiencyPlot("True muon momentum / GeV", 40u, 0.f, 1.5f, allCuts, drawErrors);
    auto plot_muCosTheta = PlottingHelper::EfficiencyPlot("True muon cos(theta)", 40u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_muPhi = PlottingHelper::EfficiencyPlot("True muon phi / rad", 40u, -3.142f, 3.142f, allCuts, drawErrors);
    auto plot_muPiAngle = PlottingHelper::EfficiencyPlot("True muon-pion opening angle", 40u, 0.f, 3.142f, allCuts, drawErrors);
    auto plot_piMomentum = PlottingHelper::EfficiencyPlot("True pion momentum / GeV", 40u, 0.f, 1.0f, allCuts, drawErrors);
    auto plot_piCosTheta = PlottingHelper::EfficiencyPlot("True pion cos(theta)", 40u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_piPhi = PlottingHelper::EfficiencyPlot("True pion phi / rad", 40u, -3.142, 3.142f, allCuts, drawErrors);
    auto plot_piMomentumGolden = PlottingHelper::EfficiencyPlot("True golden pion momentum / GeV", 40u, 0.f, 0.4f, allCuts, drawErrors);
    auto plot_piCosThetaGolden = PlottingHelper::EfficiencyPlot("True golden pion cos(theta)", 40u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_piPhiGolden = PlottingHelper::EfficiencyPlot("True golden pion phi / rad", 40u, -3.142, 3.142f, allCuts, drawErrors);

    // Efficiencies wrt CC inclusive
    auto plot_ccinc_nuEnergy = PlottingHelper::EfficiencyPlot("Neutrino energy / GeV", 30u, 0.4, 2.5f, allCuts, drawErrors);
    auto plot_ccinc_nProtons = PlottingHelper::EfficiencyPlot("Proton multiplicity", 5u, 0, 5, allCuts, drawErrors);
    auto plot_ccinc_muMomentum = PlottingHelper::EfficiencyPlot("True muon momentum / GeV", 30u, 0.f, 1.5f, allCuts, drawErrors);
    auto plot_ccinc_muCosTheta = PlottingHelper::EfficiencyPlot("True muon cos(theta)", 30u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_ccinc_muPhi = PlottingHelper::EfficiencyPlot("True muon phi / rad", 30u, -3.142f, 3.142f, allCuts, drawErrors);
    auto plot_ccinc_muPiAngle = PlottingHelper::EfficiencyPlot("True muon-pion opening angle", 30u, 0.f, 3.142f, allCuts, drawErrors);
    auto plot_ccinc_piMomentum = PlottingHelper::EfficiencyPlot("True pion momentum / GeV", 30u, 0.f, 1.0f, allCuts, drawErrors);
    auto plot_ccinc_piCosTheta = PlottingHelper::EfficiencyPlot("True pion cos(theta)", 30u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_ccinc_piPhi = PlottingHelper::EfficiencyPlot("True pion phi / rad", 30u, -3.142, 3.142f, allCuts, drawErrors);
    auto plot_ccinc_piMomentumGolden = PlottingHelper::EfficiencyPlot("True golden pion momentum / GeV", 30u, 0.f, 0.4f, allCuts, drawErrors);
    auto plot_ccinc_piCosThetaGolden = PlottingHelper::EfficiencyPlot("True golden pion cos(theta)", 30u, -1.f, 1.0f, allCuts, drawErrors);
    auto plot_ccinc_piPhiGolden = PlottingHelper::EfficiencyPlot("True golden pion phi / rad", 30u, -3.142, 3.142f, allCuts, drawErrors);

    // 2D efficiencies
    // Only do these for true efficiency (i.e. not wrt CCinclusive), and final selection
    TH2F *Generated_AllPionMomentumCosTheta = new TH2F("Generated_AllPionMomentumCosTheta",";True #pi^{+} Momentum / GeV;True #pi^{+} cos(theta) / rad", 30,0,1.5,20, -1, 1);
    TH2F *Generated_AllPionCosThetaPhi = new TH2F("Generated_AllPionCosThetaPhi",";True #pi^{+} cos(theta);True #pi^{+} phi / rad",20,-1,1,15, -TMath::Pi(), TMath::Pi());
    TH2F *Generated_AllPionMomentumPhi = new TH2F("Generated_AllPionMomentumPhi",";True #pi^{+} Momentum / GeV;True #pi^{+} phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi());

    TH2F *Generated_GoldenPionMomentumCosTheta = new TH2F("Generated_GoldenPionMomentumCosTheta",";True #pi^{+} Momentum / GeV;True #pi^{+} cos(theta) / rad", 30,0,1.5,20, -1, 1);
    TH2F *Generated_GoldenPionCosThetaPhi = new TH2F("Generated_GoldenPionCosThetaPhi",";True #pi^{+} cos(theta);True #pi^{+} phi / rad",20,-1,1,15, -TMath::Pi(), TMath::Pi());
    TH2F *Generated_GoldenPionMomentumPhi = new TH2F("Generated_GoldenPionMomentumPhi",";True #pi^{+} Momentum / GeV;True #pi^{+} phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi());

    TH2F *GenericSel_AllPionMomentumCosTheta = new TH2F("GenericSel_AllPionMomentumCosTheta",";True #pi^{+} Momentum / GeV;True #pi^{+} cos(theta) / rad", 30,0,1.5,20, -1, 1);
    TH2F *GenericSel_AllPionCosThetaPhi = new TH2F("GenericSel_AllPionCosThetaPhi",";True #pi^{+} cos(theta);True #pi^{+} phi / rad",20,-1,1,15, -TMath::Pi(), TMath::Pi());
    TH2F *GenericSel_AllPionMomentumPhi = new TH2F("GenericSel_AllPionMomentumPhi",";True #pi^{+} Momentum / GeV;True #pi^{+} phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi());

    TH2F *GenericSel_GoldenPionMomentumCosTheta = new TH2F("GenericSel_GoldenPionMomentumCosTheta",";True #pi^{+} Momentum / GeV;True #pi^{+} cos(theta) / rad", 30,0,1.5,20, -1, 1);
    TH2F *GenericSel_GoldenPionCosThetaPhi = new TH2F("GenericSel_GoldenPionCosThetaPhi",";True #pi^{+} cos(theta);True #pi^{+} phi / rad",20,-1,1,15, -TMath::Pi(), TMath::Pi());
    TH2F *GenericSel_GoldenPionMomentumPhi = new TH2F("GenericSel_GoldenPionMomentumPhi",";True #pi^{+} Momentum / GeV;True #pi^{+} phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi());

    TH2F *GoldenSel_AllPionMomentumCosTheta = new TH2F("GoldenSel_AllPionMomentumCosTheta",";True #pi^{+} Momentum / GeV;True #pi^{+} cos(theta) / rad", 30,0,1.5,20, -1, 1);
    TH2F *GoldenSel_AllPionCosThetaPhi = new TH2F("GoldenSel_AllPionCosThetaPhi",";True #pi^{+} cos(theta);True #pi^{+} phi / rad",20,-1,1,15, -TMath::Pi(), TMath::Pi());
    TH2F *GoldenSel_AllPionMomentumPhi = new TH2F("GoldenSel_AllPionMomentumPhi",";True #pi^{+} Momentum / GeV;True #pi^{+} phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi());

    TH2F *GoldenSel_GoldenPionMomentumCosTheta = new TH2F("GoldenSel_GoldenPionMomentumCosTheta",";True #pi^{+} Momentum / GeV;True #pi^{+} cos(theta) / rad", 30,0,1.5,20, -1, 1);
    TH2F *GoldenSel_GoldenPionCosThetaPhi = new TH2F("GoldenSel_GoldenPionCosThetaPhi",";True #pi^{+} cos(theta);True #pi^{+} phi / rad",20,-1,1,15, -TMath::Pi(), TMath::Pi());
    TH2F *GoldenSel_GoldenPionMomentumPhi = new TH2F("GoldenSel_GoldenPionMomentumPhi",";True #pi^{+} Momentum / GeV;True #pi^{+} phi / rad", 20,0,1.5,15, -TMath::Pi(), TMath::Pi());

    for (const auto &[run, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        if(!useThisFile) continue;
        if(sampleType != AnalysisHelper::Overlay) continue;
        const auto isMC = (sampleType != AnalysisHelper::DataEXT && sampleType != AnalysisHelper::DataBNB);
        // Read the input file
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        if (isMC) readerPeLEE.EnableSystematicBranches(); // Todo: Is this correct/optimal?
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

        // Run the selection
        std::cout << "### Only processing 5\% of events ###" << std::endl;
        for (unsigned int i = 0; i < nEvents/20; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            readerPeLEE.LoadEvent(i);
            Event event(*pEventPeLEE);
            const auto pEvent = std::make_shared<Event>(event);


            // Only care about signal events
            if (!AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg))
                continue;

            // Check which event selection cuts are passed by this event
            const auto &[isSelected, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);

            // Get the event weight
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);//*normalisation;

            // Get the features we want to plot
            const auto analysisData = AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);

            // Fill the plots at each step
            for (unsigned int i = 0; i < allCuts.size(); ++i)
            {
                const auto &cut = allCuts.at(i);
                const auto passedCut = (std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end());

                plot_nuEnergy.AddEvent(pEvent->truth.nuEnergy(), weight, cut, passedCut);
                plot_nProtons.AddEvent(analysisData.nProtons, weight, cut, passedCut);
                plot_muMomentum.AddEvent(analysisData.muonMomentum, weight, cut, passedCut);
                plot_muCosTheta.AddEvent(analysisData.muonCosTheta, weight, cut, passedCut);
                plot_muPhi.AddEvent(analysisData.muonPhi, weight, cut, passedCut);
                plot_muPiAngle.AddEvent(analysisData.muonPionAngle, weight, cut, passedCut);
                plot_piMomentum.AddEvent(analysisData.pionMomentum, weight, cut, passedCut);
                plot_piCosTheta.AddEvent(analysisData.pionCosTheta, weight, cut, passedCut);
                plot_piPhi.AddEvent(analysisData.pionPhi, weight, cut, passedCut);

                if (analysisData.hasGoldenPion)
                {
                    plot_piMomentumGolden.AddEvent(analysisData.pionMomentum, weight, cut, passedCut);
                    plot_piCosThetaGolden.AddEvent(analysisData.pionCosTheta, weight, cut, passedCut);
                    plot_piPhiGolden.AddEvent(analysisData.pionPhi, weight, cut, passedCut);
                }

                if (pEvent->reco.passesCCInclusive())
                {
                    plot_ccinc_nuEnergy.AddEvent(pEvent->truth.nuEnergy(), weight, cut, passedCut);
                    plot_ccinc_nProtons.AddEvent(analysisData.nProtons, weight, cut, passedCut);
                    plot_ccinc_muMomentum.AddEvent(analysisData.muonMomentum, weight, cut, passedCut);
                    plot_ccinc_muCosTheta.AddEvent(analysisData.muonCosTheta, weight, cut, passedCut);
                    plot_ccinc_muPhi.AddEvent(analysisData.muonPhi, weight, cut, passedCut);
                    plot_ccinc_muPiAngle.AddEvent(analysisData.muonPionAngle, weight, cut, passedCut);
                    plot_ccinc_piMomentum.AddEvent(analysisData.pionMomentum, weight, cut, passedCut);
                    plot_ccinc_piCosTheta.AddEvent(analysisData.pionCosTheta, weight, cut, passedCut);
                    plot_ccinc_piPhi.AddEvent(analysisData.pionPhi, weight, cut, passedCut);

                    if (analysisData.hasGoldenPion)
                    {
                        plot_ccinc_piMomentumGolden.AddEvent(analysisData.pionMomentum, weight, cut, passedCut);
                        plot_ccinc_piCosThetaGolden.AddEvent(analysisData.pionCosTheta, weight, cut, passedCut);
                        plot_ccinc_piPhiGolden.AddEvent(analysisData.pionPhi, weight, cut, passedCut);
                    }
                }
            }

            const auto passedGenericSelection = (std::find(cutsPassed.begin(), cutsPassed.end(), config.global.lastCutGeneric) != cutsPassed.end());

            // Fill 2D plots
            Generated_AllPionMomentumCosTheta->Fill(analysisData.pionMomentum, analysisData.pionCosTheta, weight);
            Generated_AllPionCosThetaPhi->Fill(analysisData.pionCosTheta, analysisData.pionPhi, weight);
            Generated_AllPionMomentumPhi->Fill(analysisData.pionMomentum, analysisData.pionPhi, weight);

            if (passedGenericSelection){
                GenericSel_AllPionMomentumCosTheta->Fill(analysisData.pionMomentum, analysisData.pionCosTheta, weight);
                GenericSel_AllPionCosThetaPhi->Fill(analysisData.pionCosTheta, analysisData.pionPhi, weight);
                GenericSel_AllPionMomentumPhi->Fill(analysisData.pionMomentum, analysisData.pionPhi, weight);
            }

            if (isSelected){
                GoldenSel_AllPionMomentumCosTheta->Fill(analysisData.pionMomentum, analysisData.pionCosTheta, weight);
                GoldenSel_AllPionCosThetaPhi->Fill(analysisData.pionCosTheta, analysisData.pionPhi, weight);
                GoldenSel_AllPionMomentumPhi->Fill(analysisData.pionMomentum, analysisData.pionPhi, weight);
            }


            // Fill efficiency plots for the case that there is a true golden pion
            if (analysisData.hasGoldenPion){
                Generated_GoldenPionMomentumCosTheta->Fill(analysisData.pionMomentum, analysisData.pionCosTheta, weight);
                Generated_GoldenPionCosThetaPhi->Fill(analysisData.pionCosTheta, analysisData.pionPhi, weight);
                Generated_GoldenPionMomentumPhi->Fill(analysisData.pionMomentum, analysisData.pionPhi, weight);

                if (passedGenericSelection){
                    GenericSel_GoldenPionMomentumCosTheta->Fill(analysisData.pionMomentum, analysisData.pionCosTheta, weight);
                    GenericSel_GoldenPionCosThetaPhi->Fill(analysisData.pionCosTheta, analysisData.pionPhi, weight);
                    GenericSel_GoldenPionMomentumPhi->Fill(analysisData.pionMomentum, analysisData.pionPhi, weight);
                }

                if (isSelected){
                    GoldenSel_GoldenPionMomentumCosTheta->Fill(analysisData.pionMomentum, analysisData.pionCosTheta, weight);
                    GoldenSel_GoldenPionCosThetaPhi->Fill(analysisData.pionCosTheta, analysisData.pionPhi, weight);
                    GoldenSel_GoldenPionMomentumPhi->Fill(analysisData.pionMomentum, analysisData.pionPhi, weight);
                }
            }

        }
    }

    plot_nuEnergy.SaveAs("efficiency_nuEnergy");
    plot_nProtons.SaveAs("efficiency_nProtons");
    plot_muMomentum.SaveAs("efficiency_muMomentum");
    plot_muCosTheta.SaveAs("efficiency_muCosTheta");
    plot_muPhi.SaveAs("efficiency_muPhi");
    plot_muPiAngle.SaveAs("efficiency_muPiAngle");
    plot_piMomentum.SaveAs("efficiency_piMomentum");
    plot_piCosTheta.SaveAs("efficiency_piCosTheta");
    plot_piPhi.SaveAs("efficiency_piPhi");
    plot_piMomentumGolden.SaveAs("efficiency_piMomentumGolden");
    plot_piCosThetaGolden.SaveAs("efficiency_piCosThetaGolden");
    plot_piPhiGolden.SaveAs("efficiency_piPhiGolden");

    plot_ccinc_nuEnergy.SaveAs("efficiency_ccinc_nuEnergy");
    plot_ccinc_nProtons.SaveAs("efficiency_ccinc_nProtons");
    plot_ccinc_muMomentum.SaveAs("efficiency_ccinc_muMomentum");
    plot_ccinc_muCosTheta.SaveAs("efficiency_ccinc_muCosTheta");
    plot_ccinc_muPhi.SaveAs("efficiency_ccinc_muPhi");
    plot_ccinc_muPiAngle.SaveAs("efficiency_ccinc_muPiAngle");
    plot_ccinc_piMomentum.SaveAs("efficiency_ccinc_piMomentum");
    plot_ccinc_piCosTheta.SaveAs("efficiency_ccinc_piCosTheta");
    plot_ccinc_piPhi.SaveAs("efficiency_ccinc_piPhi");
    plot_ccinc_piMomentumGolden.SaveAs("efficiency_ccinc_piMomentumGolden");
    plot_ccinc_piCosThetaGolden.SaveAs("efficiency_ccinc_piCosThetaGolden");
    plot_ccinc_piPhiGolden.SaveAs("efficiency_ccinc_piPhiGolden");

    // Save the bottom line selections
    std::vector<PlottingHelper::PlotStyle> styles({PlottingHelper::Primary, PlottingHelper::Secondary});
    std::vector<std::string> cuts({config.global.lastCutGeneric, allCuts.back()});
    plot_ccinc_nuEnergy.SaveAs(cuts, styles, "efficiency_ccinc_nuEnergy_bottomLine");
    plot_ccinc_nProtons.SaveAs(cuts, styles, "efficiency_ccinc_nProtons_bottomLine");
    plot_ccinc_muMomentum.SaveAs(cuts, styles, "efficiency_ccinc_muMomentum_bottomLine");
    plot_ccinc_muCosTheta.SaveAs(cuts, styles, "efficiency_ccinc_muCosTheta_bottomLine");
    plot_ccinc_muPhi.SaveAs(cuts, styles, "efficiency_ccinc_muPhi_bottomLine");
    plot_ccinc_muPiAngle.SaveAs(cuts, styles, "efficiency_ccinc_muPiAngle_bottomLine");
    plot_ccinc_piMomentum.SaveAs(cuts, styles, "efficiency_ccinc_piMomentum_bottomLine");
    plot_ccinc_piCosTheta.SaveAs(cuts, styles, "efficiency_ccinc_piCosTheta_bottomLine");
    plot_ccinc_piPhi.SaveAs(cuts, styles, "efficiency_ccinc_piPhi_bottomLine");
    plot_ccinc_piMomentumGolden.SaveAs(cuts, styles, "efficiency_ccinc_piMomentumGolden_bottomLine");
    plot_ccinc_piCosThetaGolden.SaveAs(cuts, styles, "efficiency_ccinc_piCosThetaGolden_bottomLine");
    plot_ccinc_piPhiGolden.SaveAs(cuts, styles, "efficiency_ccinc_piPhiGolden_bottomLine");

    auto pCanvas = PlottingHelper::GetCanvas();

    GenericSel_AllPionMomentumCosTheta->Sumw2();
    GenericSel_AllPionMomentumCosTheta->Divide(Generated_AllPionMomentumCosTheta);
    GenericSel_AllPionMomentumCosTheta->GetZaxis()->SetRangeUser(0,0.3);
    GenericSel_AllPionMomentumCosTheta->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GenericSel_AllPionMomentumCosTheta");
    GenericSel_AllPionCosThetaPhi->Sumw2();
    GenericSel_AllPionCosThetaPhi->Divide(Generated_AllPionCosThetaPhi);
    GenericSel_AllPionCosThetaPhi->GetZaxis()->SetRangeUser(0,0.3);
    GenericSel_AllPionCosThetaPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GenericSel_AllPionCosThetaPhi");
    GenericSel_AllPionMomentumPhi->Sumw2();
    GenericSel_AllPionMomentumPhi->Divide(Generated_AllPionMomentumPhi);
    GenericSel_AllPionMomentumPhi->GetZaxis()->SetRangeUser(0,0.3);
    GenericSel_AllPionMomentumPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GenericSel_AllPionMomentumPhi");
    GenericSel_GoldenPionMomentumCosTheta->Sumw2();
    GenericSel_GoldenPionMomentumCosTheta->Divide(Generated_GoldenPionMomentumCosTheta);
    GenericSel_GoldenPionMomentumCosTheta->GetZaxis()->SetRangeUser(0,0.3);
    GenericSel_GoldenPionMomentumCosTheta->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GenericSel_GoldenPionMomentumCosTheta");
    GenericSel_GoldenPionCosThetaPhi->Sumw2();
    GenericSel_GoldenPionCosThetaPhi->Divide(Generated_GoldenPionCosThetaPhi);
    GenericSel_GoldenPionCosThetaPhi->GetZaxis()->SetRangeUser(0,0.3);
    GenericSel_GoldenPionCosThetaPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GenericSel_GoldenPionCosThetaPhi");
    GenericSel_GoldenPionMomentumPhi->Sumw2();
    GenericSel_GoldenPionMomentumPhi->Divide(Generated_GoldenPionMomentumPhi);
    GenericSel_GoldenPionMomentumPhi->GetZaxis()->SetRangeUser(0,0.3);
    GenericSel_GoldenPionMomentumPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GenericSel_GoldenPionMomentumPhi");

    GoldenSel_AllPionMomentumCosTheta->Sumw2();
    GoldenSel_AllPionMomentumCosTheta->Divide(Generated_AllPionMomentumCosTheta);
    GoldenSel_AllPionMomentumCosTheta->GetZaxis()->SetRangeUser(0,0.3);
    GoldenSel_AllPionMomentumCosTheta->Draw("colz");
    GoldenSel_AllPionCosThetaPhi->Sumw2();
    PlottingHelper::SaveCanvas(pCanvas,"GoldenSel_AllPionMomentumCosTheta");
    GoldenSel_AllPionCosThetaPhi->Divide(Generated_AllPionCosThetaPhi);
    GoldenSel_AllPionCosThetaPhi->GetZaxis()->SetRangeUser(0,0.3);
    GoldenSel_AllPionCosThetaPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GoldenSel_AllPionCosThetaPhi");
    GoldenSel_AllPionMomentumPhi->Sumw2();
    GoldenSel_AllPionMomentumPhi->Divide(Generated_AllPionMomentumPhi);
    GoldenSel_AllPionMomentumPhi->GetZaxis()->SetRangeUser(0,0.3);
    GoldenSel_AllPionMomentumPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GoldenSel_AllPionMomentumPhi");
    GoldenSel_GoldenPionMomentumCosTheta->Sumw2();
    GoldenSel_GoldenPionMomentumCosTheta->Divide(Generated_GoldenPionMomentumCosTheta);
    GoldenSel_GoldenPionMomentumCosTheta->GetZaxis()->SetRangeUser(0,0.3);
    GoldenSel_GoldenPionMomentumCosTheta->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GoldenSel_GoldenPionMomentumCosTheta");
    GoldenSel_GoldenPionCosThetaPhi->Sumw2();
    GoldenSel_GoldenPionCosThetaPhi->Divide(Generated_GoldenPionCosThetaPhi);
    GoldenSel_GoldenPionCosThetaPhi->GetZaxis()->SetRangeUser(0,0.3);
    GoldenSel_GoldenPionCosThetaPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GoldenSel_GoldenPionCosThetaPhi");
    GoldenSel_GoldenPionMomentumPhi->Sumw2();
    GoldenSel_GoldenPionMomentumPhi->Divide(Generated_GoldenPionMomentumPhi);
    GoldenSel_GoldenPionMomentumPhi->GetZaxis()->SetRangeUser(0,0.3);
    GoldenSel_GoldenPionMomentumPhi->Draw("colz");
    PlottingHelper::SaveCanvas(pCanvas,"GoldenSel_GoldenPionMomentumPhi");
}

} // namespace ubcc1pi_macros
