/**
 *  @file  ubcc1pi_standalone/Macros/OptimizeSidebandCuts.cxx
 *
 *  @brief The implementation file of the OptimizeSidebandCuts macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"
// #include <fstream> // Todo: don't use txt files
#include <iomanip>

#include <TH2F.h>
#include <TMath.h>

float muonCutValue(const int k) { return  -0.60f+k*0.05f; }
float protonCutValue(const int k) { return 0.05f+k*0.025f; }
// float muonCutValue(const int k) { return  -0.3f-k*0.01f; }
// float protonCutValue(const int k) { return 0.4f+k*0.01f; }

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void OptimizeSidebandCuts(const Config &config)
{
    ROOT::EnableImplicitMT(2);

    // Define a function scoped to this macro to determine if a cut was passed
    auto passedCut = [&](const std::string &cut, const std::vector<std::string> &cutsPassed) -> bool {
        return (std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end());
    };

    // Get the selections
    const std::string yLabel = "Number of particles (norm. to bin width)";
    auto selectionPart1 = SelectionHelper::GetCC0piSelectionModifiedPeLEEPart1(true);

    std::map<std::pair<int, int>, TH1F> TrueCC0pi_ProtonMomentum_SelCC0pi, TrueCC0pi_nProtons_SelCC0pi;
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<3; j++)
        {
            const auto coords = std::make_pair(i, j);
            auto histProtonMomentum = TH1F(("TrueCC0pi_SelCC0pi_ProtonMomentum"+std::to_string(i)+"+"+std::to_string(j)).c_str(),(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5);
            auto histnProtons = TH1F(("TrueCC0pi_SelCC0pi_nProtons"+std::to_string(i)+"+"+std::to_string(j)).c_str(),(string(";Proton Multiplicity (-1);")+yLabel).c_str(), 3, 0, 3);
            TrueCC0pi_ProtonMomentum_SelCC0pi.emplace(coords, histProtonMomentum);
            TrueCC0pi_nProtons_SelCC0pi.emplace(coords, histnProtons);
        }
    }

    auto TrueCC0pi_ProtonMomentum_SelCC1pi = std::make_shared<TH1F>("TrueCC0pi_ProtonMomentum_SelCC1pi",(string(";Proton Momentum / GeV;")+yLabel).c_str(), 15, 0, 1.5);
    auto TrueCC0pi_nProtons_SelCC1pi = std::make_shared<TH1F>("TrueCC0pi_nProtons_SelCC1pi",(string(";Proton Multiplicity (-1);")+yLabel).c_str(), 3, 0, 3);

    auto CC1piSelection = SelectionHelper::GetDefaultSelection2(true);

    const auto featureNames = BDTHelper::ParticleBDTFeatureNames;
    // Get the muon BDT
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);

    // Get the proton BDT
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);

    // Get the golden pion BDT
    const auto goldenPionFeatureNames = BDTHelper::GoldenPionBDTFeatureNames;
    BDTHelper::BDT goldenpionBDT("goldenPion", goldenPionFeatureNames);

    auto summedWeights = 0.f;
    // Loop over the events
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        if(!useThisFile) continue;
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        if(!isOverlay) continue;
        const auto isMC = true;

        std::cout << "Reading input file: " << filePath << std::endl;
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        std::cout << "    -> Opened file successfully" << std::endl;
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();

        std::cout<<"\n##############\nOnly counting 5\% of event!\n##############"<<std::endl;
        for (unsigned int e = 0; e < nEvents/20; ++e)
        {
            AnalysisHelper::PrintLoadingBar(e, nEvents);
            readerPeLEE.LoadEvent(e);

            Event event(*pEventPeLEE, true);// true or false decides wether to cut generation!=2 particles
            const auto pEvent = std::make_shared<Event>(event);

            const auto isTrueCC0Pi = AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);
            if (!isTrueCC0Pi) continue;

            // Get the truth and reco analysis data
            // const auto plotStyle = PlottingHelper::GetPlotStyle(sampleType, pEvent, config.global.useAbsPdg);

            const auto recoParticles = pEvent->reco.particles;
            const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent) * normalisation;
            // // Get the muon reco data
            // const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, muonFeatureNames, muonBDT);
            // const auto muon = recoParticles.at(muonIndex);
            // const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
            // const float muonCosTheta = muonDir.Z();
            // const float muonPhi = std::atan2(muonDir.Y(), muonDir.X());
            // const float muonMomentum = AnalysisHelper::GetMuonMomentum(muon);

            // // Get the proton reco data
            // float protonMomentum = -std::numeric_limits<float>::max();
            // float protonCosTheta = -std::numeric_limits<float>::max();
            // float protonPhi = -std::numeric_limits<float>::max();

            const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
            const auto nProtons = std::min(AnalysisHelper::CountParticlesAboveMomentumThreshold(visibleParticles, 2212, config.global.useAbsPdg, config.global.protonMomentumThreshold)-1, 2u); //-1 as one proton is treated as the pion

            // std::cout<<"DEBUG OptimiseSidebandCuts Point 0"<<std::endl;
            const auto &[isSelectedPart1, cutsPassedPart1, assignedPdgCodesPart1] = selectionPart1.Execute(pEvent);
            // std::cout<<"DEBUG OptimiseSidebandCuts Point 1"<<std::endl;
            if(isSelectedPart1)
            {
                // std::cout<<"DEBUG OptimiseSidebandCuts Point 1.1"<<std::endl;
                const auto recoMuonIndex = AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesPart1, 13);
                // std::cout<<"DEBUG OptimiseSidebandCuts Point 1.2"<<std::endl;
                auto selectionPart2 = SelectionHelper::GetCC0piSelectionModifiedPeLEEPart2(recoMuonIndex, 0.f, 0.f);
                // std::cout<<"DEBUG OptimiseSidebandCuts Point 1.3"<<std::endl;
                // const auto protonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(recoParticles, assignedPdgCodesPart1);

                // Default value in case no proton can be found
                unsigned int protonIndex = std::numeric_limits<unsigned int>::max();
                float leadingProtonMom = 0;
                for (unsigned int i = 0; i < recoParticles.size(); ++i) // Loop copied and modified from SelectionHelper::GetLeadingProtonCandidateIndex
                {
                    if (i == recoMuonIndex)
                        continue;

                    // Now check if this is the leading proton (i.e. the highest-reconstructed-momentum proton in the event)
                    const auto proton = recoParticles.at(i);
                    if (!AnalysisHelper::HasTrackFit(proton))
                        continue;
                    float protonmom = AnalysisHelper::GetProtonMomentumFromRange(proton.range());
                    if (protonmom > leadingProtonMom){
                        leadingProtonMom = protonmom;
                        protonIndex = i;
                    }
                }

                if(protonIndex == std::numeric_limits<unsigned int>::max())
                    throw std::logic_error("protonIndex == std::numeric_limits<unsigned int>::max()");

                const auto leadingRecoProton = recoParticles.at(protonIndex);
                // std::cout<<"DEBUG OptimiseSidebandCuts Point 1.3.1"<<std::endl;
                unsigned int leadpidx = std::numeric_limits<unsigned int>::max(); // todo check if this is the right way to do this !!!!!!!!!!!!!!!!!
                try{
                    leadpidx = AnalysisHelper::GetBestMatchedTruthParticleIndex(leadingRecoProton, pEvent->truth.particles);
                }
                catch(const std::logic_error &){
                    std::cout<<"try failed for leadpidx\n"<<std::endl;
                    // throw std::logic_error("try failed for leadpidx");
                    // continue;
                }
                // auto leadpidx = AnalysisHelper::GetTrueLeadingProtonIndex(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold);
                // std::cout<<"DEBUG OptimiseSidebandCuts Point 1.4"<<std::endl;
                if (leadpidx!=std::numeric_limits<unsigned int>::max())
                {
                    // std::cout<<"DEBUG OptimiseSidebandCuts Point 2"<<std::endl;
                    summedWeights += weight;
                    for(int i=0; i<4; i++)
                    {
                        for(int j=0; j<3; j++)
                        {
                            // // std::cout<<"DEBUG OptimiseSidebandCuts Point 3"<<std::endl;
                            // Run the event selection and store which cuts are passed
                            const auto coords = std::make_pair(i, j);
                            selectionPart2.SetCutValue("muonLikeProton", muonCutValue(i));
                            selectionPart2.SetCutValue("barelyResemblingProton", protonCutValue(j));
                            // std::cout<<"DEBUG OptimiseSidebandCuts Point 4"<<std::endl;
                            const auto &[selectedGoldenCC0pi, cutsPassedPart2, assignedPdgCodes] = selectionPart2.Execute(pEvent);                            
                            // std::cout<<"DEBUG OptimiseSidebandCuts Point 5"<<std::endl;
                            // Due to the split of the selction into two parts, the last cut of the generic selection here is different than in the other macros
                            const auto selectedGenericCC0pi = SelectionHelper::IsCutPassed(cutsPassedPart2, "openingAngle");
                            // std::cout<<"DEBUG OptimiseSidebandCuts Point 5.1"<<std::endl;
                            // If event passes CC0pi selection, fill plots
                            if (selectedGenericCC0pi)
                            {
                                // std::cout<<"DEBUG OptimiseSidebandCuts Point 6"<<std::endl;
                                // For proton kinematics, use true leading proton. If no true leading proton is found, don't fill any kinematics
                                const auto trueleadp = pEvent->truth.particles.at(leadpidx);
                                TrueCC0pi_ProtonMomentum_SelCC0pi[coords].Fill(trueleadp.momentum(), weight);
                                TrueCC0pi_nProtons_SelCC0pi[coords].Fill(nProtons, weight);
                            }
                            // // std::cout<<"DEBUG OptimiseSidebandCuts Point 7"<<std::endl;
                        }
                    }
                }
            }
            // std::cout<<"DEBUG OptimiseSidebandCuts Point 8"<<std::endl;

            if(passedCut("topologicalScoreCC", cutsPassedPart1)) // Not very optimised but at least check that the event previously passed the ccinc selection for cc0pi
            {
                // // std::cout<<"DEBUG OptimiseSidebandCuts Point 9"<<std::endl;
                // Find out if the event passes the generic CC1pi selection
                const auto &[isSelectedCC1pi, cutsPassedCC1pi, assignedPdgCodesCC1pi] = CC1piSelection.Execute(pEvent);
                const auto &isSelectedCC1piGeneric = (std::find(cutsPassedCC1pi.begin(),cutsPassedCC1pi.end(), config.global.lastCutGeneric) != cutsPassedCC1pi.end());

                // If event passes CC1pi+ selection, fill plots
                if (isSelectedCC1piGeneric){
                    // // std::cout<<"DEBUG OptimiseSidebandCuts Point 10"<<std::endl;
                    // For proton kinematics, use proton that is selected as the pi+ candidate. If the pi+ candidate is a true proton, fill muon and proton kinematic plots
                    const auto recoPion = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1pi, 211));
                    // const auto recoMuon = recoParticles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1pi, 13));
                    Event::Truth::Particle pionMatch;

                    try{
                        pionMatch = AnalysisHelper::GetBestMatchedTruthParticle(recoPion,pEvent->truth.particles);
                    }
                    catch(const std::logic_error &){
                        continue;
                    }
                    // // std::cout<<"DEBUG OptimiseSidebandCuts Point 11"<<std::endl;

                    // if (pionMatch.pdgCode()==2212){
                    TrueCC0pi_ProtonMomentum_SelCC1pi->Fill(pionMatch.momentum(), weight);
                    TrueCC0pi_nProtons_SelCC1pi->Fill(nProtons, weight); // The number of protons here does not include the missidentified proton that is a pion candidate 
                    // }
                    // // std::cout<<"DEBUG OptimiseSidebandCuts Point 12"<<std::endl;
                }
            }
        }
    }

    std::cout<<"-----summedWeights: "<<summedWeights<<std::endl;

    // TrueCC0pi_ProtonMomentum_SelCC1pi->Sumw2();
    const auto selCC1piIntegralProtonMomentum = TrueCC0pi_ProtonMomentum_SelCC1pi->Integral();
    const auto selCC1piIntegralnProtons = TrueCC0pi_nProtons_SelCC1pi->Integral();
    TrueCC0pi_ProtonMomentum_SelCC1pi->Scale(1.0/selCC1piIntegralProtonMomentum);
    TrueCC0pi_nProtons_SelCC1pi->Scale(1.0/selCC1piIntegralnProtons);
    std::cout<<"cc0piIntegralProtonMomentum: "<<selCC1piIntegralProtonMomentum<<std::endl;
    std::cout<<"cc0piIntegralnProtons: "<<selCC1piIntegralnProtons<<std::endl;
    std::map<std::pair<int, int>, float> selcc0piIntegralProtonMomentum, selcc0piIntegralnProtons;

    std::cout<<"Results (proton momentum):"<<std::endl;
    std::vector<float> values(6, std::numeric_limits<float>::max());
    std::vector<float> lowestValues(6, std::numeric_limits<float>::max());
    std::vector<std::pair<int, int>> bestCoords(6, std::make_pair(-1,-1));
    // std::pair<int, int> bestCoordsProtonMomentumOnly, bestCoordsEfficiencyWeightedProtonMomentumOnly, bestCoords, bestCoordsEfficiencyWeighted;
    std::cout<<"\t"<<std::scientific<<std::setprecision(2);
    for(int j=0; j<3; j++)
        std::cout<<"j:"<<protonCutValue(j)<<"\t";
    std::cout<<std::endl;
    for(int i=0; i<4; i++)
    {
        std::cout<<"i:"<<muonCutValue(i)<<"\t";
        for(int j=0; j<3; j++)
        {
            const auto coords = std::make_pair(i, j);

            const auto cc0piIntegralProtonMomentum = TrueCC0pi_ProtonMomentum_SelCC0pi[coords].Integral();
            const auto cc0piIntegralnProtons = TrueCC0pi_nProtons_SelCC0pi[coords].Integral();
            selcc0piIntegralProtonMomentum.emplace(coords, cc0piIntegralProtonMomentum);
            selcc0piIntegralnProtons.emplace(coords, cc0piIntegralnProtons);
            // std::cout<<"cc0piIntegralProtonMomentum - i: "<<i<<" - j: "<<j<<" "<<cc0piIntegralProtonMomentum<<std::endl;
            TrueCC0pi_ProtonMomentum_SelCC0pi[coords].Scale(1.0/cc0piIntegralProtonMomentum);
            TrueCC0pi_ProtonMomentum_SelCC0pi[coords].GetYaxis()->SetRangeUser(0,0.23);
            TrueCC0pi_nProtons_SelCC0pi[coords].Scale(1.0/cc0piIntegralnProtons);
            TrueCC0pi_nProtons_SelCC0pi[coords].GetYaxis()->SetRangeUser(0,0.8);

            const auto efficiencyWeightProtonMomentum = selCC1piIntegralProtonMomentum/cc0piIntegralProtonMomentum;
            const auto efficiencyWeightnProtons = selCC1piIntegralnProtons/cc0piIntegralnProtons;

            const auto subtractedProtonMomentum = *TrueCC0pi_ProtonMomentum_SelCC1pi - TrueCC0pi_ProtonMomentum_SelCC0pi[coords];
            const auto subtractednProtons = *TrueCC0pi_nProtons_SelCC1pi - TrueCC0pi_nProtons_SelCC0pi[coords];
            const auto normalisedDifferenceProtonMomentum = subtractedProtonMomentum/TrueCC0pi_ProtonMomentum_SelCC0pi[coords];
            const auto normalisedDifferencenProtons = subtractednProtons/TrueCC0pi_nProtons_SelCC0pi[coords];
            const auto squaredProtonMomentum = normalisedDifferenceProtonMomentum*normalisedDifferenceProtonMomentum;
            const auto squarednProtons = normalisedDifferencenProtons*normalisedDifferencenProtons;
            const auto squaredProtonMomentumIntegral = squaredProtonMomentum.Integral()/squaredProtonMomentum.GetNbinsX();
            const auto squarednProtonsIntegral = squarednProtons.Integral()/squarednProtons.GetNbinsX();
            // const auto value = squared.Integral()/squared.GetNbinsX();

            std::cout<<"i: "<<i<<"j: "<<j<<" eff.W.nP: "<<efficiencyWeightnProtons<<" sqrdPMomentum: "<<squaredProtonMomentumIntegral<<" sqrdnP: "<<squarednProtonsIntegral<<std::endl;
            // std::cout<<"subtractedProtonMomentum: "<<subtractedProtonMomentum.Integral()<<"  subtractednProtons: "<<subtractednProtons.Integral()<<std::endl;
            // std::cout<<"normalisedDifferenceProtonMomentum: "<<normalisedDifferenceProtonMomentum.Integral()<<"  normalisedDifferencenProtons: "<<normalisedDifferencenProtons.Integral()<<std::endl;
            // std::cout<<"squaredProtonMomentum: "<<squaredProtonMomentumIntegral<<" squarednProtons: "<<squarednProtonsIntegral<<std::endl;

            values[0] = squaredProtonMomentumIntegral;
            values[1] = efficiencyWeightProtonMomentum*squaredProtonMomentumIntegral*squaredProtonMomentumIntegral;
            values[2] = squaredProtonMomentumIntegral*squaredProtonMomentumIntegral*squarednProtonsIntegral;
            values[3] = efficiencyWeightProtonMomentum*squaredProtonMomentumIntegral*squarednProtonsIntegral*squaredProtonMomentumIntegral;
            values[4] = efficiencyWeightProtonMomentum*std::sqrt(squarednProtonsIntegral)*squaredProtonMomentumIntegral;
            values[5] = std::sqrt(efficiencyWeightProtonMomentum*squarednProtonsIntegral)*squaredProtonMomentumIntegral;
            for(int m=0; m<6; m++)
            {
                if(values[m] < lowestValues[m])
                {
                    lowestValues[m] = values[m];
                    bestCoords[m] = coords;
                }
            }
        }
        std::cout<<std::endl;
    }
    auto i = bestCoords[0].first;
    auto j = bestCoords[0].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum considered): "<<lowestValues[0]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[1].first;
    j = bestCoords[1].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + number of selected events considered): "<<lowestValues[1]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[2].first;
    j = bestCoords[2].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + proton multiplicity considered): "<<lowestValues[2]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[3].first;
    j = bestCoords[3].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + proton multiplicity + number of selected events considered): "<<lowestValues[3]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[4].first;
    j = bestCoords[4].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum + sqrt(proton multiplicity) + number of selected events considered): "<<lowestValues[4]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;

    i = bestCoords[5].first;
    j = bestCoords[5].second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Lowest value (proton momentum * sqrt(proton multiplicity * number of selected events considered)): "<<lowestValues[5]<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;


    std::cout<<"# events:"<<std::endl;
    auto highestValue = -std::numeric_limits<float>::max();
    std::pair<int, int> bestCoordsNumEvents;
    std::cout<<"           \t";
    for(int j=0; j<3; j++)
        std::cout<<"j:"<<protonCutValue(j)<<"\t";
    std::cout<<std::endl;
    for(int i=0; i<4; i++)
    {
        std::cout<<"i:"<<muonCutValue(i)<<"\t";
        for(int j=0; j<3; j++)
        {
            const auto coords = std::make_pair(i, j);
            const auto value = selcc0piIntegralProtonMomentum[coords];
            if(value>highestValue)
            {
                highestValue = value;
                bestCoordsNumEvents = coords;
            }
            std::cout<<value<<"\t";
        }
        std::cout<<std::endl;
    }
    i = bestCoordsNumEvents.first;
    j = bestCoordsNumEvents.second;
    std::cout<<"----------------------------------------------------"<<std::endl;
    std::cout<<"Largest number of events (CC0pi selection): "<<highestValue<<" at "<<i<<"(muonCutValue: "<<muonCutValue(i)<<"),"<<j<<"(protonCutValue: "<<protonCutValue(j)<<")"<<std::endl;
    std::cout<<"----------------------------------------------------"<<std::endl;



    auto pCanvas = PlottingHelper::GetCanvas();

    PlottingHelper::SetLineStyle(TrueCC0pi_ProtonMomentum_SelCC1pi, PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[0]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[1]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[2]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[3]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[4]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[5]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(TrueCC0pi_nProtons_SelCC1pi, PlottingHelper::Primary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[0]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[1]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[2]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[3]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[4]], PlottingHelper::Secondary);
    PlottingHelper::SetLineStyle(&TrueCC0pi_nProtons_SelCC0pi[bestCoords[5]], PlottingHelper::Secondary);

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[0]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_raw_ProtonMomentumOnly");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[0]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_raw_ProtonMomentumOnly");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[1]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted_ProtonMomentumOnly");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[1]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted_ProtonMomentumOnly");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[2]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_raw");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[2]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_raw");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[3]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[3]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[4]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted_sqrt1");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[4]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted_sqrt1");

    TrueCC0pi_ProtonMomentum_SelCC0pi[bestCoords[5]].Draw("hist");
    TrueCC0pi_ProtonMomentum_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_ProtonMomentum_areanorm_weighted_sqrt2");

    TrueCC0pi_nProtons_SelCC0pi[bestCoords[5]].Draw("hist");
    TrueCC0pi_nProtons_SelCC1pi->Draw("hist same");
    PlottingHelper::SaveCanvas(pCanvas,"OptimizeSidebandCuts_nProtons_areanorm_weighted_sqrt2");
}

} // namespace ubcc1pi_macros
