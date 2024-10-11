/**
 *  @file  ubcc1pi_standalone/Macros/TestSidebandSelectionsPeLEEVsUbcc1pi.cxx
 *
 *  @brief The implementation file of the TestSidebandSelectionsPeLEEVsUbcc1pi macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
// #include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
// #include "ubcc1pi_standalone/Helpers/FittingHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
// #include "ubcc1pi_standalone/ubsmear/inc/ubsmear/Helpers/UBSmearingHelper.h"

// Boost libraries
// #include "binary_iarchive.hpp"
#include "binary_oarchive.hpp"
#include "binary_object.hpp"
#include "map.hpp"
#include "vector.hpp"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void TestSidebandSelectionsPeLEEVsUbcc1pi(const Config &config)
{
    std::cout<<"..........................................\nUSING Modified CC0pi Selection: muonLikeProtonValue=-0.48f, barelyResemblingProtonValue=0.12f\n.........................................."<<std::endl;
    auto sidebandSelection = SelectionHelper::GetCC0piSelectionModified(-0.48f, 0.12f);

    ExtractionHelper::InputFileList inputData;

    inputData.push_back(std::make_tuple(AnalysisHelper::Overlay,"",std::string("/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root"),1.f));
    inputData.push_back(std::make_tuple(AnalysisHelper::Overlay,"",std::string("/uboone/data/users/jdetje/ubcc1piVSpelee/ubcc1pi/ubcc1piAnalysis_0_4k.root"),1.f));

    // Loop over the files
    FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(std::get<2>(inputData.at(0)), true);
    FileReader<Event, Subrun> reader(std::get<2>(inputData.at(1)), true);

    if(readerPeLEE.GetNumberOfEvents()!=reader.GetNumberOfEvents())
    {
        std::cout<<"ERROR: Number of events in PeLEE and ubcc1pi files are different - PeLEE: "<<readerPeLEE.GetNumberOfEvents()<<" vs ubcc1pi: "<<reader.GetNumberOfEvents()<<std::endl;
        throw std::logic_error("Unequal number of events in PeLEE and ubcc1pi files!");
    }

    // const auto isPeLEE    = (sampleType == AnalysisHelper::Overlay);
    // const auto isUbcc1pi  = (sampleType == AnalysisHelper::DataBNB);

    auto pEvent = reader.GetBoundEventAddress();
    auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

    // Loop over the events in the file
    if(readerPeLEE.GetNumberOfEvents() != reader.GetNumberOfEvents())
        throw std::logic_error("Number of events in PeLEE and ubcc1pi files are different!");

    const auto nEvents = readerPeLEE.GetNumberOfEvents();

    // const auto style = isPeLEE ? PlottingHelper::ExternalPoints : PlottingHelper::External;//BNBData;

    auto nIdenticalGenericSelection = 0u;
    auto nIdenticalGoldenSelection = 0u;
    auto nDifferentCuts = 0u;
    auto nDifferentPDGs = 0u;
    for (unsigned int i = 0; i < nEvents; ++i)
    {
        readerPeLEE.LoadEvent(i);
        reader.LoadEvent(i);

        // Create a shared pointer to the event object
        Event event(*pEventPeLEE, true);
        const auto pEventPeLEEAsUbCC1Pi = std::make_shared<Event>(event);

        // ******************************************************************************************
        // Check truth
        // ******************************************************************************************

        // Determine if this is truly a CC0Pi event
        // const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
        const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
        const auto isTrueCC0PiPeLEE = AnalysisHelper::IsTrueCC0Pi(pEventPeLEEAsUbCC1Pi, config.global.useAbsPdg, config.global.protonMomentumThreshold);
        if(isTrueCC0Pi != isTrueCC0PiPeLEE)
        {
            std::cout<<"ERROR: isTrueCC0Pi != isTrueCC0PiPeLEE"<<std::endl;
            throw std::logic_error("isTrueCC0Pi != isTrueCC0PiPeLEE");
        }

        // ******************************************************************************************
        // Check the selection
        // ******************************************************************************************
        const auto &[passedGoldenSelectionUbcc1pi, cutsPassedUbcc1pi, assignedPdgCodesUbcc1pi] = sidebandSelection.Execute(pEvent);
        const auto passedGenericSelectionUbcc1pi = SelectionHelper::IsCutPassed(cutsPassedUbcc1pi, config.global.lastCutGeneric);
        const auto &[passedGoldenSelectionPelee, cutsPassedPelee, assignedPdgCodesPelee] = sidebandSelection.Execute(pEventPeLEEAsUbCC1Pi);
        const auto passedGenericSelectionPelee = SelectionHelper::IsCutPassed(cutsPassedPelee, config.global.lastCutGeneric);

        for (unsigned int c = 0; c < cutsPassedUbcc1pi.size(); c++) {
            const auto cutUbcc1pi = cutsPassedUbcc1pi.at(c);
            if (std::find(cutsPassedPelee.begin(), cutsPassedPelee.end(), cutUbcc1pi) == cutsPassedPelee.end())
            {
                std::cout<<"! Event "<<i<<" ubcc1pi cut: "<<cutUbcc1pi<<" was not passed in pelee."<<std::endl;
                nDifferentCuts++;
                break;
            }
        }

        for (unsigned int c = 0; c < cutsPassedPelee.size(); c++) {
            const auto cutPelee = cutsPassedPelee.at(c);
            if (std::find(cutsPassedUbcc1pi.begin(), cutsPassedUbcc1pi.end(), cutPelee) == cutsPassedUbcc1pi.end())
            {
                std::cout<<"! Event "<<i<<" pelee cut: "<<cutPelee<<" was not passed in ubcc1[i]."<<std::endl;
                nDifferentCuts++;
                break;
            }
        }

        if(passedGenericSelectionPelee == passedGenericSelectionUbcc1pi)
        {
            nIdenticalGenericSelection++;
            // std::cout << "✓ Event " << i << " generic selections are identical (both " << passedGenericSelectionPelee  << ")" << std::endl;

            if(assignedPdgCodesUbcc1pi != assignedPdgCodesPelee && passedGenericSelectionPelee)
            {
                nDifferentPDGs++;
                std::cout<< "! Event "<<i<< " assigned pdg codes are different."<<std::endl;
                if(assignedPdgCodesUbcc1pi.size() == assignedPdgCodesPelee.size())
                {
                    for(unsigned int p=0; p<assignedPdgCodesUbcc1pi.size(); p++)
                    {
                        std::cout<<"    PeLEE: "<<assignedPdgCodesPelee.at(p)<<" vs UBCC1Pi: "<<assignedPdgCodesUbcc1pi.at(p)<<std::endl;
                    }
                }
                else
                {
                    std::cout<<"    Different number of assigned PDG codes: PeLEE: "<<assignedPdgCodesPelee.size()<<" vs UBCC1Pi: "<<assignedPdgCodesUbcc1pi.size()<<std::endl;
                }
            }

            if(passedGoldenSelectionPelee == passedGoldenSelectionUbcc1pi)
            {
                nIdenticalGoldenSelection++;
                // std::cout << "\t✓✓ Event " << i << " golden selections are identical (both " << passedGoldenSelectionPelee  << ")" << std::endl;
            }
            else
            {
                std::cout<< "✓✗ Event "<<i<< " golden selections are different: PeLEE: "<< passedGoldenSelectionPelee <<" vs UBCC1Pi: "<< passedGoldenSelectionUbcc1pi <<std::endl;

                std::cout<<"    passed PeLEE cuts:";
                for(const auto &cut : cutsPassedPelee) std::cout<<" "<<cut;
                std::cout<<std::endl;

                std::cout<<"    passed UBCC1Pi cuts:";
                for(const auto &cut : cutsPassedUbcc1pi) std::cout<<" "<<cut;
                std::cout<<std::endl;
            }
        }
        else
        {
            std::cout<< "✗ Event "<<i<< " generic selections are different: PeLEE: "<< passedGenericSelectionPelee <<" vs UBCC1Pi: "<< passedGenericSelectionUbcc1pi <<std::endl;

            std::cout<<"    passed PeLEE cuts:";
            for(const auto &cut : cutsPassedPelee) std::cout<<" "<<cut;
            std::cout<<std::endl;

            std::cout<<"    passed UBCC1Pi cuts:";
            for(const auto &cut : cutsPassedUbcc1pi) std::cout<<" "<<cut;
            std::cout<<std::endl;
        }
    }
    std::cout<<"Generic selection: "<<nIdenticalGenericSelection<<"/"<<nEvents<<std::endl;
    std::cout<<"Golden selection: "<<nIdenticalGoldenSelection<<"/"<<nEvents<<std::endl;
    std::cout<<"Different cuts: "<<nDifferentCuts<<std::endl;
    std::cout<<"Different PDGs: "<<nDifferentPDGs<<std::endl;

    // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // // auto pFileIn = new TFile(fileName, "READ");
    // // auto pFileIn = std::make_shared<TFile>(fileName.c_str(), "READ"); // std::unique_ptr<TFile>(myFile(TFile::Open(fileName, "READ")));
    // // auto pTreeEvents = (TTree*)pFileIn->Get("events");
    // // auto pTreeSubruns = (TTree*)pFileIn->Get("subruns");
    // // // pFileIn->Close();
    // // auto pFileOut = std::make_shared<TFile>("/uboone/data/users/jdetje/ubcc1pi/ubcc1pi_stage2/test.root", "RECREATE");
    // // pFileOut->WriteObject(pTreeEvents, "events");
    // // pFileOut->WriteObject(pTreeSubruns, "subruns");
    // // auto pTree = (TTree*)pFileOut->Get("events");
    // // pFileIn->Close();

    // // std::filesystem::copy(fileName, "/uboone/data/users/jdetje/ubcc1pi/ubcc1pi_stage2/test.root", std::filesystem::copy_options::overwrite_existing); // copy file

    // const std::string path = "/uboone/data/users/jdetje/ubcc1pi/ubcc1pi_stage2/test.root";
    //          // If there is no accessible file, copy first; avoids copying when using this macro consecutively
    //  // ATTN: This does not check whether file is corrupted or not; if in doubt, clear output directory before running
    // // auto pFileOut = std::make_shared<TFile>(path.c_str(), "UPDATE");
    //         // if(!pFileOut || pFileOut->IsZombie())//
    // if(access( path.c_str(), F_OK ) == -1)
    // {
    //     // pFileOut->Close();
    //     auto pFileIn = std::make_shared<TFile>(fileName.c_str(), "READ"); // std::unique_ptr<TFile>(myFile(TFile::Open(fileName, "READ")));
    //     pFileIn->Cp(path.c_str());
    //     pFileIn->Close();
    //     // pFileOut = std::make_shared<TFile>(path.c_str(), "UPDATE");
    // }


    // auto pFileOut = std::make_shared<TFile>(path.c_str(), "UPDATE");
    // auto pTreeFriend = (TTree*)pFileOut->Get("treeFriend");
    // auto pTree = (TTree*)pFileOut->Get("events");
    // if(!pTreeFriend)
    // {
    //     pTree->RemoveFriend(pTreeFriend);
    //     pTreeFriend->Delete();
    // }

    // // bool o_isSelectedGenericCC0Pi, o_isSelectedGenericCC1Pi, o_isSelectedGoldenCC0Pi, o_isSelectedGoldenCC1Pi;
    // // pTree->Branch("isSelectedGenericCC0Pi", &o_isSelectedGenericCC0Pi);
    // // pTree->Branch("isSelectedGenericCC1Pi", &o_isSelectedGenericCC1Pi);
    // // pTree->Branch("isSelectedGoldenCC0Pi", &o_isSelectedGoldenCC0Pi);
    // // pTree->Branch("isSelectedGoldenCC1Pi", &o_isSelectedGoldenCC1Pi);

    // // bool o_isTrueCC0Pi, o_isTrueCC1Pi, o_isSignalCC0Pi, o_isSignalCC1Pi;
    // // if(isOverlay || isDetVar || isNuWro) // Only MC
    // // {
    // //     pTree->Branch("isTrueCC0Pi", &o_isTrueCC0Pi);
    // //     pTree->Branch("isTrueCC1Pi", &o_isTrueCC1Pi);
    // //     pTree->Branch("isSignalCC0Pi", &o_isSignalCC0Pi);
    // //     pTree->Branch("isSignalCC1Pi", &o_isSignalCC1Pi);
    // // }

    // std::vector<Bool_t> isSelectedGenericCC0PiVect(nEvents, false); //todo remove UBCC1Pi value -  only for debugging with >100% of events!!!!!!!!!!!!!!!!!!
    // std::vector<Bool_t> isSelectedGenericCC1PiVect(nEvents, false);
    // std::vector<Bool_t> isSelectedGoldenCC0PiVect(nEvents, false);
    // std::vector<Bool_t> isSelectedGoldenCC1PiVect(nEvents, false);
    // std::vector<Bool_t> isTrueCC0PiVect(nEvents, false);
    // std::vector<Bool_t> isTrueCC1PiVect(nEvents, false);
    // std::vector<Bool_t> isSignalCC0PiVect(nEvents, false);
    // std::vector<Bool_t> isSignalCC1PiVect(nEvents, false);

    // for(unsigned h = 2*nEvents/3; h < nEvents; ++h)
    // {
    //     isSelectedGenericCC0PiVect.at(h) = true;
    //     isSignalCC1PiVect.at(h) = true;
    // }

    // auto b1 = pTree->Branch("isSignalCC1Pi", &isSignalCC1PiVect, "isSignalCC1Pi/O");


    // for (unsigned int i = 0; i < 5000/*nEvents*/; ++i) // todo change back to all events!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // {
    //     AnalysisHelper::PrintLoadingBar(i, int(nEvents));
    //     reader.LoadEvent(i);

    //     // ############################################################
    //     // ############### Run the sideband selection #################
    //     // ############################################################
    //     const auto &[passedGoldenSidebandSelection, cutsPassed, assignedPdgCodes] = sidebandSelection.Execute(pEvent);

    //     const auto passedGenericSidebandSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

    //     // o_isSelectedGenericCC0Pi = passedGenericSidebandSelection;
    //     // o_isSelectedGoldenCC0Pi =  passedGoldenSidebandSelection;
    //     isSelectedGenericCC0PiVect.at(i) = passedGenericSidebandSelection;
    //     isSelectedGoldenCC0PiVect.at(i) = passedGoldenSidebandSelection;

    //     // Get the reco analysis data (if available, otherwise set to dummy values)
    //     const auto recoSidebandData = (
    //         passedGenericSidebandSelection
    //             ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodes, passedGoldenSidebandSelection)
    //             : AnalysisHelper::GetDummyAnalysisData()
    //     );

    //     // Here we apply reco-level phase-space restrictions
    //     // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
    //     // min/max values supplied in the binning. If so, then reject the event.
    //     auto passesSidebandPhaseSpaceReco = false;
    //     if (passedGenericSidebandSelection)
    //     {
    //         // Start by assuming the event passes the phase-space cuts
    //         passesSidebandPhaseSpaceReco = true;

    //         // Check the value of the kinematic quantities are within the phase-space limits
    //         for (const auto &[name, minMax] : phaseSpaceMap)
    //         {
    //             const auto &[min, max] = minMax;
    //             const auto value = getSidebandValue.at(name)(recoSidebandData);

    //             if (value < min || value > max)
    //             {
    //                 passesSidebandPhaseSpaceReco = false;
    //                 break;
    //             }
    //         }
    //     }
    //     const auto isSelectedSidebandGeneric = passedGenericSidebandSelection && passesSidebandPhaseSpaceReco;
    //     const auto isSelectedSidebandGolden = passedGoldenSidebandSelection && passesSidebandPhaseSpaceReco;


    //     // ############################################################
    //     // #################### Run the main selection #####################
    //     // ############################################################
    //     const auto &[passedGoldenSelection, cutsPassed, assignedPdgCodes] = selection.Execute(pEvent);
    //     const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassed, config.global.lastCutGeneric);

    //     // o_isSelectedGenericCC1Pi = passedGenericSelection;
    //     // o_isSelectedGoldenCC1Pi =  passedGoldenSelection;
    //     isSelectedGenericCC1PiVect.at(i) = passedGenericSelection;
    //     isSelectedGoldenCC1PiVect.at(i) = passedGoldenSelection;

    //     // Get the reco analysis data (if available, otherwise set to dummy values)
    //     const auto recoData = (
    //         passedGenericSelection
    //             ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodes, passedGoldenSelection)
    //             : AnalysisHelper::GetDummyAnalysisData()
    //     );

    //     // Here we apply reco-level phase-space restrictions
    //     // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
    //     // min/max values supplied in the binning. If so, then reject the event.
    //     auto passesPhaseSpaceReco = false;
    //     if (passedGenericSelection)
    //     {
    //         // Start by assuming the event passes the phase-space cuts
    //         passesPhaseSpaceReco = true;

    //         // Check the value of the kinematic quantities are within the phase-space limits
    //         for (const auto &[name, minMax] : phaseSpaceMap)
    //         {
    //             const auto &[min, max] = minMax;
    //             const auto value = getValue.at(name)(recoData);

    //             if (value < min || value > max)
    //             {
    //                 passesPhaseSpaceReco = false;
    //                 break;
    //             }
    //         }
    //     }
    //     const auto isSelectedGolden = passedGoldenSelection && passesPhaseSpaceReco;
    //     const auto isSelectedGeneric = passedGenericSelection && passesPhaseSpaceReco;


    //     if (isDataBNB) continue; // For BNB data that's all we need to do!

    //     // Determine if this is truly a CC0Pi event
    //     const auto isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
    //     // o_isTrueCC0Pi = isTrueCC0Pi;
    //     isTrueCC0PiVect.at(i) = isTrueCC0Pi;

    //     // Get the truth analysis data (if available, otherwise set to dummy values)
    //     const auto truthSidebandData = (
    //         (isTrueCC0Pi)
    //             ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
    //             : AnalysisHelper::GetDummyAnalysisData()
    //     );

    //     // Here we apply truth-level phase-space restrictions
    //     // For all true CC0Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
    //     // event is not classed as "signal"
    //     bool passesSidebandPhaseSpaceTruth = false;
    //     if (isTrueCC0Pi)
    //     {
    //         // Start by assuming the event passes the phase-space cuts
    //         passesSidebandPhaseSpaceTruth = true;

    //         // Check the value of the kinematic quantities are within the phase-space limits
    //         for (const auto &[name, minMax] : phaseSpaceMap)
    //         {
    //             // if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection // todo check this
    //             const auto &[min, max] = minMax;
    //             const auto value = getSidebandValue.at(name)(truthSidebandData);

    //             if (value < min || value > max)
    //             {
    //                 passesSidebandPhaseSpaceTruth = false;
    //                 break;
    //             }
    //         }
    //     }
    //     const auto isCC0PiSignal = isTrueCC0Pi && passesSidebandPhaseSpaceTruth;
    //     // o_isSignalCC0Pi = isCC0PiSignal;
    //     isSignalCC0PiVect.at(i) = isCC0PiSignal;


    //     const auto isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);
    //     // o_isTrueCC1Pi = isTrueCC1Pi;
    //     isTrueCC1PiVect.at(i) = isTrueCC1Pi;

    //     // Get the truth analysis data (if available, otherwise set to dummy values)
    //     const auto truthData = (
    //         isTrueCC1Pi
    //             ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
    //             : AnalysisHelper::GetDummyAnalysisData()
    //     );

    //     // Here we apply truth-level phase-space restrictions
    //     // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
    //     // event is not classed as "signal"
    //     bool passesPhaseSpaceTruth = false;
    //     if (isTrueCC1Pi)
    //     {
    //         // Start by assuming the event passes the phase-space cuts
    //         passesPhaseSpaceTruth = true;

    //         // Check the value of the kinematic quantities are within the phase-space limits
    //         for (const auto &[name, minMax] : phaseSpaceMap)
    //         {
    //             const auto &[min, max] = minMax;
    //             const auto value = getValue.at(name)(truthData);

    //             if (value < min || value > max)
    //             {
    //                 passesPhaseSpaceTruth = false;
    //                 break;
    //             }
    //         }
    //     }

    //     const auto isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruth;
    //     // o_isSignalCC1Pi = isCC1PiSignal;
    //     // isSignalCC1PiVect.at(i) = isCC1PiSignal;
    //     isSignalCC1PiVect.at(i) = true; //todo remove debugging only

    //     // pTree->Fill();
    // } // End of event-level iteration

    // Float_t isSelectedGenericCC0Pi = 44.14f;

    // // auto pB0 = pTree->GetBranch("isSelectedGenericCC0Pi");
    // // if(pB0)
    // // {
    // //     auto pLeaf = pTree->GetLeaf("isSelectedGenericCC0Pi");
    // //     pTree->GetListOfBranches()->Remove(pB0);
    // //     pTree->GetListOfLeaves()->Remove(pLeaf);
    // //     pTree->Write();
    // // }

    // // const auto pB1 = pTree->Branch("isSelectedGenericCC0Pi", &isSelectedGenericCC0Pi, "isSelectedGenericCC0Pi/F");
    // TTree treeFriend("treeFriend", "ubcc1pi");
    // const auto pB1 = treeFriend.Branch("isSelectedGenericCC0Pi", &isSelectedGenericCC0Pi, "isSelectedGenericCC0Pi/F");

    // // // Disable everything...
    // // pTree->SetBranchStatus("*", true);
    // // // ...but the branch we need
    // // pTree->SetBranchStatus("isSelectedGenericCC0Pi", true);
    // // pTree->SetBranchAddress("isSelectedGenericCC0Pi", &isSelectedGenericCC0Pi);
    // Long64_t nentries = pTree->GetEntries();
    // for (Long64_t k = 0; k < nentries; k++) {
    //     if(k%100==0) std::cout<<(1.0*k)/nentries<<std::endl;
    //     // pTree->GetEntry(k);
    //     if(k<nentries/3) isSelectedGenericCC0Pi = 41.3f; //isSelectedGenericCC0PiVect.at(k);
    //     else isSelectedGenericCC0Pi = 41.5f;
    //     pB1->Fill(); // Use Branch->BackFill() with newer ROOT versions
    // }

    // treeFriend.Write("", TObject::kOverwrite);
    // pTree->AddFriend(&treeFriend);
    // pTree->Write("", TObject::kOverwrite);
    // std::cout<<std::endl;

    // // // pTree->Branch("isSelectedGenericCC0Pi", &isSelectedGenericCC0PiVect, "isSelectedGenericCC0Pi/O")->Fill();
    // // pTree->Branch("isSelectedGenericCC1Pi", &isSelectedGenericCC1PiVect, "isSelectedGenericCC1Pi/O")->Fill();
    // // pTree->Branch("isSelectedGoldenCC0Pi", &isSelectedGoldenCC0PiVect, "isSelectedGoldenCC0Pi/O")->Fill();
    // // pTree->Branch("isSelectedGoldenCC1Pi", &isSelectedGoldenCC1PiVect, "isSelectedGoldenCC1Pi/O")->Fill();


    // // for(int l = 0; l < 50; l++)
    // // {
    // //     std::cout<<" "<<isTrueCC0PiVect.at(l);
    // // }
    // // std::cout<<"\n--------------------"<<std::endl;
    // // for(int l = 0; l < 50; l++)
    // // {
    // //     std::cout<<" "<<isTrueCC0PiVect.at(nEvents-1-l);
    // // }
    // // std::cout<<std::endl;

    // // if(isOverlay || isDetVar || isNuWro) // Only MC
    // // {
    // //     pTree->Branch("isTrueCC0Pi", &isTrueCC0PiVect, "isTrueCC0Pi/O")->Fill();
    // //     pTree->Branch("isTrueCC1Pi", &isTrueCC1PiVect, "isTrueCC1Pi/O")->Fill();
    // //     pTree->Branch("isSignalCC0Pi", &isSignalCC0PiVect, "isSignalCC0Pi/O")->Fill();
    // //     pTree->Branch("isSignalCC1Pi", &isSignalCC1PiVect, "isSignalCC1Pi/O")->Fill();
    // //     // b1->Fill();
    // // }

    // // pTree->Fill();

    // //////////////////////////////////////////////////////////////
    // // Add values to input ntuples and save them as new files
    // //////////////////////////////////////////////////////////////
    // pFileOut->Write("", TObject::kOverwrite);
    // pFileOut->Close();
    // break; //TODO Remove after first tests !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    // truth_particle_pdgCode_plot.SaveAs("truth_particle_pdgCode");
    // truth_particle_endX_plot.SaveAs("truth_particle_endX");
    // truth_particle_endY_plot.SaveAs("truth_particle_endY");
    // truth_particle_endZ_plot.SaveAs("truth_particle_endZ");
    // truth_particle_momentumX_plot.SaveAs("truth_particle_momentumX");
    // truth_particle_momentumY_plot.SaveAs("truth_particle_momentumY");
    // truth_particle_momentumZ_plot.SaveAs("truth_particle_momentumZ");
    // // truth_particle_momentum_plot
    // truth_particle_energy_plot.SaveAs("truth_particle_energy");

    std::cout<<"-----------------Done-----------------"<<std::endl;


}

} // namespace ubcc1pi_macros
