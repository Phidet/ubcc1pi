/**
 *  @file  ubcc1pi_standalone/Macros/PlotMuonConfusion.cxx
 *
 *  @brief The implementation file of the PlotMuonConfusion macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"

#include <TStyle.h>
#include <TColor.h>
#include <TText.h>
#include <TPaletteAxis.h>
#include <TLatex.h>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PlotMuonConfusion(const Config &config)
{
    const auto binMax = 99.0; 
    // Create two vectors with bin edges
    // std::vector<double> muonMomentumBinEdges = {0.15, 0.23, 0.32, 0.45, 0.66, 1.5, binMax};
    // std::vector<double> pionMomentumBinEdges = {0.0, 0.1, 0.16, 0.19, 0.22, 0.6, binMax};
    // const cool addPercentagesToCells = true;

    std::vector<double> muonMomentumBinEdges = {0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75 ,0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, binMax};
    std::vector<double> pionMomentumBinEdges = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 ,0.4, 0.45, 0.5, 0.55, 0.6, 0.65, binMax};
    const bool addPercentagesToCells = false;

    // Create a new 2D ROOT histogram
    TH2D *hist_correctPID = new TH2D("hist_correctPID", "Right muon candidate;X;Y", 
                        muonMomentumBinEdges.size() - 1, &muonMomentumBinEdges[0], 
                        pionMomentumBinEdges.size() - 1, &pionMomentumBinEdges[0]);

    TH2D *hist_wrongPID = new TH2D("hist_wrongPID", "Wrong muon candidate;X;Y", 
                        muonMomentumBinEdges.size() - 1, &muonMomentumBinEdges[0], 
                        pionMomentumBinEdges.size() - 1, &pionMomentumBinEdges[0]);

    // Check if the sizes of input and output files are identical
    if (config.input.files.size() != config.output.files.size())
        throw std::runtime_error("Input and output file sizes are not identical");

    // Loop over the input files
    for (unsigned int f = 0; f < config.input.files.size(); ++f)
    {
        const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] = config.input.files.at(f);
        const auto &[fileRunOut, normalisationOut, sampleTypeOut, useThisFileOut, filePathOut] = config.output.files.at(f);

        if(!useThisFileOut) continue;
        if(sampleTypeOut != AnalysisHelper::Overlay) continue;

        std::cout << "Processing file - " << filePath << " and " << filePathOut << std::endl;
        const auto isMC = true;

        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        readerPeLEE.EnableSystematicBranches();
        FileReader<EventXSec, SubrunXSec> readerXSec(filePathOut, isMC);
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        const auto pEventXSec = readerXSec.GetBoundEventAddress();
        // std::cout<<"DEBUG PlotMuonConfusion 3"<<std::endl;

        // Loop over the events
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto nEventsXSec = readerXSec.GetNumberOfEvents();
        if(nEvents != nEventsXSec) throw std::runtime_error("Number of events in input and output files are not identical");

        // std::cout << "### Only processing 50\% of events ###" << std::endl;
        for (unsigned int i = 0; i < nEvents; i++) //nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            // std::cout<<"DEBUG PlotMuonConfusion 4"<<std::endl;
            readerXSec.LoadEvent(i);
            // std::cout<<"DEBUG PlotMuonConfusion 5"<<std::endl;
            const auto isCC1PiSignal = pEventXSec->truth.cc1pi_signal();
            const auto isSelectedGenericCC1pi = pEventXSec->reco.cc1pi_selected_generic();
            // std::cout<<"DEBUG PlotMuonConfusion 6"<<std::endl;
            if(!isCC1PiSignal || !isSelectedGenericCC1pi) continue; // Only consider selected signal events for these confusion plots

            AnalysisHelper::PrintLoadingBar(i, nEvents);
            readerPeLEE.LoadEvent(i);
            Event event(*pEventPeLEE, true); // true: cut out particles with generation > 2
            const auto pEvent = std::make_shared<Event>(event);

            // const auto isCC1Pi = pEventXSec->truth.true_cc1pi();
            const auto trueMuonMomentum = pEventXSec->truth.cc1pi_truth_muonMomentum();
            const auto truePionMomentum = pEventXSec->truth.cc1pi_truth_pionMomentum();
            const auto recoMuonCosTheta = pEventXSec->reco.cc1pi_reco_muonCosTheta();


             // Match muon candidate values in the xsec_analyzer ntuple to one of the reco particles in the PelEE ntuple
            int matchedIndex = -1;
            const auto& recoParticles = pEvent->reco.particles;
            auto oldDiff = std::numeric_limits<double>::max();
            for (unsigned int p = 0; p < recoParticles.size(); ++p) {
                const auto recoParticle = recoParticles.at(p);
                if(!recoParticle.directionX.IsSet() || !recoParticle.directionY.IsSet() || !recoParticle.directionZ.IsSet())
                    continue;
                const auto dir = TVector3(recoParticle.directionX(), recoParticle.directionY(), recoParticle.directionZ()).Unit();
                const auto cosTheta = dir.Z();
                const auto diff = std::abs(cosTheta - recoMuonCosTheta); 
                if ( diff < 1e-6) { // use some small threshold to compare floating point numbers
                    if(diff < oldDiff)
                    {
                        if(oldDiff != std::numeric_limits<double>::max()) std::cout << "WARNING: Multiple reco particles with similar cosTheta" << std::endl;
                        matchedIndex = p;
                        oldDiff = diff;
                    }
                }
            } // end of loop over reco particles

            if (matchedIndex == -1) {
                 throw std::runtime_error("No matching particle found");
            }

            // Match reco particle to truth particle
            const auto truthParticles = pEvent->truth.particles;
            if(!recoParticles.at(matchedIndex).pdgBacktracked.IsSet())
            {
                std::cout << "WARNING: pdgBacktracked not set" << std::endl;
                continue;
            }
            const auto truthParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticles.at(matchedIndex), truthParticles);
            const auto muonCandiateIsTrueMuon = std::abs(truthParticles.at(truthParticleIndex).pdgCode()) == 13;
            // std::cout << "muonCandiateIsTrueMuon: " << muonCandiateIsTrueMuon << " with pdgCode: " << truthParticles.at(truthParticleIndex).pdgCode();
            // std::cout << " and trueMuonMom: " << trueMuonMomentum << " and truePionMom: " << truePionMomentum << std::endl;

            // Get the nominal event weight, scaled by the sample normalisation
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;

            // Fill the histogram with the true muon and pion momentum as x and y
            if(muonCandiateIsTrueMuon)
            {
                hist_correctPID->Fill(trueMuonMomentum, truePionMomentum, weight);
            }
            else
            {
                hist_wrongPID->Fill(trueMuonMomentum, truePionMomentum, weight);
            }
        } // end of loop over events
    } // end of loop over input files


    TH2D *hist_totalPID = new TH2D("hist_totalPID", "Total muon candidate;X;Y",
                        muonMomentumBinEdges.size() - 1, &muonMomentumBinEdges[0],
                        pionMomentumBinEdges.size() - 1, &pionMomentumBinEdges[0]);

    hist_totalPID->Add(hist_correctPID, hist_wrongPID);
    
    TH2D *hist_ratioPID = new TH2D("hist_ratioPID", "Muon candidate ratio;X;Y",
                            muonMomentumBinEdges.size() - 1, &muonMomentumBinEdges[0],
                            pionMomentumBinEdges.size() - 1, &pionMomentumBinEdges[0]);

    hist_ratioPID->Divide(hist_correctPID, hist_totalPID, 1, 1, "B");

    // Copy contents of hist_ratioPID into histogram with all bins equally sized
    int numBinsX = hist_ratioPID->GetNbinsX();
    int numBinsY = hist_ratioPID->GetNbinsY();
    TH2D *hist_equalBins = new TH2D("hist_equalBins", "Correctly identified muon fraction for selected signal events;X;Y",
                                    numBinsX, 0, numBinsX,
                                    numBinsY, 0, numBinsY);

    for (int i = 1; i <= numBinsX; i++)
    {
        for (int j = 1; j <= numBinsY; j++)
        {
            double content = hist_ratioPID->GetBinContent(i, j);
            hist_equalBins->SetBinContent(i, j, content);
        }
    }

    // // Print out the values in each bin of hist_correctPID:
    // std::cout << "Values in hist_correctPID:" << std::endl;
    // for (int i = 1; i <= hist_correctPID->GetNbinsX(); i++)
    // {
    //     for (int j = 1; j <= hist_correctPID->GetNbinsY(); j++)
    //     {
    //         std::cout << hist_correctPID->GetBinContent(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // // Print out the values in each bin of hist_wrongPID:
    // std::cout << "Values in hist_wrongPID:" << std::endl;
    // for (int i = 1; i <= hist_wrongPID->GetNbinsX(); i++)
    // {
    //     for (int j = 1; j <= hist_wrongPID->GetNbinsY(); j++)
    //     {
    //         std::cout << hist_wrongPID->GetBinContent(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // // Print out the values in each bin of hist_ratioPID:
    // std::cout << "Values in hist_ratioPID:" << std::endl;
    // for (int i = 1; i <= hist_ratioPID->GetNbinsX(); i++)
    // {
    //     for (int j = 1; j <= hist_ratioPID->GetNbinsY(); j++)
    //     {
    //         std::cout << hist_ratioPID->GetBinContent(i, j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // // Print out the bin edges for hist_ratioPID:
    // std::cout << "Bin edges for hist_ratioPID:" << std::endl;
    // for (int i = 1; i <= hist_ratioPID->GetNbinsX(); i++)
    // {
    //     std::cout << hist_ratioPID->GetXaxis()->GetBinLowEdge(i) << ", ";
    // }
    // std::cout << std::endl;
    // for (int i = 1; i <= hist_ratioPID->GetNbinsY(); i++)
    // {
    //     std::cout << hist_ratioPID->GetYaxis()->GetBinLowEdge(i) << ", ";
    // }
    // std::cout << std::endl;

    

    // Set the color palette to inverse mint
    gStyle->SetPalette(kMint);
    TColor::InvertPalette();

    // Create a new canvas
    TCanvas *canvas = new TCanvas("canvas", "Muon Confusion Plot", 800, 600);
    // Increase the right margin
    canvas->SetRightMargin(0.15);

    // Plot the 2D histogram as colz
    hist_equalBins->Draw("colz");

    // Remove the x axis ticks
    hist_equalBins->GetXaxis()->SetTickLength(0);

    // Remove the y axis ticks
    hist_equalBins->GetYaxis()->SetTickLength(0);

    // Decrease the font size of the x axis labels
    hist_equalBins->GetXaxis()->SetLabelSize(0.05);

    // Decrease the font size of the y axis labels
    hist_equalBins->GetYaxis()->SetLabelSize(0.05);

    // Decrease the font size of the z axis labels
    // hist_equalBins->GetZaxis()->SetLabelSize(0.05);

    // Label the x axis
    hist_equalBins->GetXaxis()->SetTitle("True muon momentum");

    // Label the y axis
    hist_equalBins->GetYaxis()->SetTitle("True pion momentum");

    // Set the plot title
    hist_equalBins->SetTitleSize(0.05);
    hist_equalBins->SetTitle("Correctly identified muon fraction for selected signal events");

    hist_equalBins->GetXaxis()->SetLabelOffset(999); // Hide x axis labels
    hist_equalBins->GetYaxis()->SetLabelOffset(999); // Hide y axis labels
    hist_equalBins->GetZaxis()->SetLabelOffset(0); // Hide z axis labels


    // Increase the offset of the y axis title
    hist_equalBins->GetYaxis()->SetTitleOffset(1.5);

    // Set the x axis ticks
    for (unsigned int i = 0; i < muonMomentumBinEdges.size() - 1; i++)
    {
        TText *text = new TText(i, -0.15, Form("%.2f", muonMomentumBinEdges[i]));
        text->SetTextAlign(22); // Center alignment
        text->SetTextSize(0.03); // Adjust text size as needed
        text->Draw();
    }

    // Set the y axis ticks
    for (unsigned int i = 0; i < pionMomentumBinEdges.size() - 1; i++)
    {
        TText *text = new TText(-0.4, i, Form("%.2f", pionMomentumBinEdges[i]));
        text->SetTextAlign(22); // Center alignment
        text->SetTextSize(0.03); // Adjust text size as needed
        text->Draw();
    }

    if(addPercentagesToCells)
    {
        // Write the percentage value of each bin into the cell on the plot
        for (int i = 1; i <= hist_equalBins->GetNbinsX(); i++)
        {
            for (int j = 1; j <= hist_equalBins->GetNbinsY(); j++)
            {
                double content = hist_equalBins->GetBinContent(i, j);
                double percentage = content * 100.0;
                TLatex *latex = new TLatex(hist_equalBins->GetXaxis()->GetBinCenter(i), hist_equalBins->GetYaxis()->GetBinCenter(j), Form("#splitline{%.1f%%}{(%.1f)}", percentage, hist_totalPID->GetBinContent(i, j)));
                latex->SetTextAlign(22); // Center alignment
                latex->SetTextSize(0.03); // Adjust text size as needed
                latex->Draw();
            }
        }
    }
    
    // Get the palette
    TPaletteAxis *palette = (TPaletteAxis*)hist_equalBins->GetListOfFunctions()->FindObject("palette");

    // Check if the palette exists
    if (palette) {
        // Decrease the font size of the color legend labels
        palette->SetLabelSize(0.03);
    }

    // Update the canvas
    canvas->Update();

    // Save the canvas as a PDF
    canvas->SaveAs("muonEfficiency.pdf");

} // end of macro
} // namespace ubcc1pi_macros
