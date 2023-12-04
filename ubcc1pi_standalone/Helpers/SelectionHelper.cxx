/**
 *  @file  ubcc1pi_standalone/Helpers/SelectionHelper.cxx
 *
 *  @brief The implementation file for the selection helper class
 */

#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

namespace ubcc1pi
{

SelectionHelper::EventSelection::Cut::Cut() :
    m_name(""),
    m_hasValue(false),
    m_value(-std::numeric_limits<float>::max())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::Cut::Cut(const std::string &name) :
    m_name(name),
    m_hasValue(false),
    m_value(-std::numeric_limits<float>::max())
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::Cut::Cut(const std::string &name, const float value) :
    m_name(name),
    m_hasValue(true),
    m_value(value)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string SelectionHelper::EventSelection::Cut::GetName() const
{
    if (m_name.empty())
        throw std::logic_error("Cut::GetName - Name has not been set");

    return m_name;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::Cut::HasValue() const
{
    return m_hasValue;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::Cut::GetValue() const
{
    if (!m_hasValue)
        throw std::logic_error("Cut::GetValue - Can't get value of cut: \"" + m_name + "\" - it has no value");

    return m_value;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::Cut::SetValue(const float value)
{
    if (!m_hasValue)
        throw std::logic_error("Cut::SetValue - Can't set value of cut: \"" + m_name + "\" - it has no value");

    m_value = value;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::CutTracker::CutTracker(const std::vector<Cut> &cuts, const unsigned int nRecoParticles) :
    m_cuts(cuts),
    m_assignedPdgs(nRecoParticles, 0) // ATTN here we use a PDG code of zero to mean "unassigned"
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::CutTracker::GetCutValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("CutTracker::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->GetValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> SelectionHelper::EventSelection::CutTracker::GetCutsPassed() const
{
    return m_cutsPassed;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<int> SelectionHelper::EventSelection::CutTracker::GetAssignedPDGCodes() const
{
    return m_assignedPdgs;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::CutTracker::MarkCutAsPassed(const std::string &name)
{
    // Count the number of cuts that have already passed
    const auto nCutsPassed = m_cutsPassed.size();

    // Check we haven't passed all cuts already
    if (nCutsPassed == m_cuts.size())
        throw std::logic_error("CutTracker::MarkCutAsPassed - All cuts are already marked as passed");

    // Check that the cut we are marking as "passed" matches up with the next cut in the vector
    if (name != m_cuts.at(nCutsPassed).GetName())
        throw std::logic_error("CutTracker::MarkCutAsPassed - The cut \"" + name + "\" is not the next cut in the sequence");

    m_cutsPassed.push_back(name);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::CutTracker::AssignPDGCode(const unsigned int recoParticleIndex, const int pdg)
{
    if (recoParticleIndex >= m_assignedPdgs.size())
        throw std::out_of_range("CutTracker::AssignPDGCode - The input recoParticleIndex is out of range");

    // Set the PDG code
    m_assignedPdgs.at(recoParticleIndex) = pdg;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::EventSelection(const std::vector<Cut> &cuts, BDTMap &bdtMap, const SelectionLogic &logic) :
    m_cuts(cuts),
    m_bdtMap(bdtMap),
    m_logic(logic)
{
    // Make sure none of the input cut names are not repeated
    for (const auto &cut : m_cuts)
    {
        const auto cutName = cut.GetName();

        if (std::count_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) { return x.GetName() == cutName; }) != 1)
            throw std::invalid_argument("EventSelection::EventSelection - Repeated cut name: \"" + cutName + "\"");
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> SelectionHelper::EventSelection::GetCuts() const
{
    std::vector<std::string> cutNames;
    for (const auto &cut : m_cuts)
        cutNames.push_back(cut.GetName());

    return cutNames;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::EventSelection::CutHasValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::CutHasValue - Unknown cut: \"" + name + "\"");

    return iter->HasValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float SelectionHelper::EventSelection::GetCutValue(const std::string &name) const
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->GetValue();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SelectionHelper::EventSelection::SetCutValue(const std::string &name, const float value)
{
    // Find the cut by name
    const auto iter = std::find_if(m_cuts.begin(), m_cuts.end(), [&](const auto &x) {
        return x.GetName() == name;
    });

    if (iter == m_cuts.end())
        throw std::invalid_argument("EventSelection::GetCutValue - Unknown cut: \"" + name + "\"");

    return iter->SetValue(value);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection::SelectionResult SelectionHelper::EventSelection::Execute(const std::shared_ptr<Event> &pEvent)
{
    // Setup the cut tracker
    CutTracker cutTracker(m_cuts, pEvent->reco.particles.size());

    // Pass the event, BDTs and cut tracker through to the selection logic and run the selection
    const auto passed = m_logic(pEvent, m_bdtMap, cutTracker);

    // Return the selection result
    return {
        passed,
        cutTracker.GetCutsPassed(),
        cutTracker.GetAssignedPDGCodes()
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetSelection(const std::string &name)
{
    if (name == "CCInclusive")//todo switch
        return SelectionHelper::GetCCInclusiveSelection();
    else if (name == "Default")
        return SelectionHelper::GetDefaultSelection();
    else if (name == "CC0pi")
        return SelectionHelper::GetCC0piSelection();
    else if (name == "CC0piModified")
        return SelectionHelper::GetCC0piSelectionModified();
    else if (name == "CC0piModifiedPeLEE")
        return SelectionHelper::GetCC0piSelectionModifiedPeLEE();
    else throw std::invalid_argument("Cannot find selection: \"" + name + "\"");
};


// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCCInclusiveSelection()
{
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"passesCCInclusive"}
    };

    // We don't need any BDTs so just use an empty map
    EventSelection::BDTMap bdtMap;

    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // We must pass the CC inclusive selection
        if (!pEvent->reco.passesCCInclusive())
            return false;

        // Mark the cut "passesCCInclusive" as passed
        cutTracker.MarkCutAsPassed("passesCCInclusive");

        // Identify the muon candidate
        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (particle.isCCInclusiveMuonCandidate())
            {
                // Assign the muon candidate a PDG code of 13
                cutTracker.AssignPDGCode(i, 13);
            }
        }

        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetOriginalSelection()
{
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"2NonProtons", -0.06f},
        {"pionHasValiddEdx", 1.0f},
        {"pionNotInGap"},
        {"muonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f},
        {"likelyGoldenPion", -0.03f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)},
        {"goldenPion", std::make_shared<BDTHelper::BDT>("goldenPion", BDTHelper::GoldenPionBDTFeatureNames)}
    };

    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // ----------------------------------------------------------------------------------
        // passesCCInclusive
        // ----------------------------------------------------------------------------------

        // Insist the event passes the CC inclusive preselection
        if (!pEvent->reco.passesCCInclusive())
            return false;

        // Mark the cut "passesCCInclusive" as passed
        cutTracker.MarkCutAsPassed("passesCCInclusive");

        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");
        //std::cout<<"CC1pi: min2Tracks"<<std::endl;

        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;

        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");
        //std::cout<<"CC1pi: max1Uncontained"<<std::endl;

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // ----------------------------------------------------------------------------------
        // 2NonProtons
        // ----------------------------------------------------------------------------------

        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThreshold = cutTracker.GetCutValue("2NonProtons");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        std::vector<unsigned int> pionIndices;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
            {
                continue;
            }

            // Assume particles without a track fit are just small protons
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("GetOriginalSelection - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            if (!hasFeatures)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThreshold)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nPions + nMuons != recoParticles.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen!");

        // Insist that we exacly one pion (i.e we have have 2 non-protons)
        if (nPions != 1)
            return false;

        // Mark the cut "2NonProtons" as passed
        cutTracker.MarkCutAsPassed("2NonProtons");
        //std::cout<<"CC1pi: 2NonProtons"<<std::endl;

        // ----------------------------------------------------------------------------------
        // pionHasValiddEdx
        // ----------------------------------------------------------------------------------

        const auto piondEdxThreshold = cutTracker.GetCutValue("pionHasValiddEdx");

        // Get the pion reco particle
        const auto pionIndex = pionIndices.front();
        const auto &pion = recoParticles.at(pionIndex);

        // Insist that the pion has a valid deEdx measurement (i.e. truncated mean dE/dx is above threshold)
        if (pion.truncatedMeandEdx() < piondEdxThreshold)
            return false;

        // Mark the cut "pionHasValiddEdx" as passed
        cutTracker.MarkCutAsPassed("pionHasValiddEdx");
        //std::cout<<"CC1pi: pionHasValiddEdx"<<std::endl;

        // ----------------------------------------------------------------------------------
        // pionNotInGap
        // ----------------------------------------------------------------------------------

        // Sanity check that our muon and pion are not the same particle
        if (muonIndex == pionIndex)
            throw std::logic_error("DefaultSelection - The muon and the pion candidates are the same particle. This shouldn't happen!");

        // Insist that the pion has at least one hit in each view (i.e. not in a gap)
        if (pion.nHitsU() == 0 || pion.nHitsV() == 0 || pion.nHitsW() == 0)
            return false;

        // Mark the cut "pionNotInGap" as passed
        cutTracker.MarkCutAsPassed("pionNotInGap");
        //std::cout<<"CC1pi: pionNotInGap"<<std::endl;

        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");
        //std::cout<<"CC1pi: muonNotInGap"<<std::endl;

        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------

        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");

        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(pionDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");
        //std::cout<<"CC1pi: openingAngle"<<std::endl;

        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");
        //std::cout<<"CC1pi: topologicalScore"<<std::endl;
        
        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------

        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
                return false;
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");
        //std::cout<<"CC1pi: startNearVertex"<<std::endl;

        // ----------------------------------------------------------------------------------
        // likelyGoldenPion
        // ----------------------------------------------------------------------------------

        // Get the likely golden pion cut value
        const auto goldenPionBDTThreshold = cutTracker.GetCutValue("likelyGoldenPion");

        // Get the golden pion BDT
        auto &pGoldenPionBDT = bdtMap.at("goldenPion");

        // Get the features of the pion
        std::vector<float> features;
        if (!BDTHelper::GetBDTFeatures(pion, BDTHelper::GoldenPionBDTFeatureNames, features))
            throw std::logic_error("DefaultSelection - Can't get golden pion BDT features of pion candidate");

        // Insist that the BDT response is greater than the cut value to identify the pion as a golden pion
        const auto bdtResponse = pGoldenPionBDT->GetResponse(features);
        if (bdtResponse <= goldenPionBDTThreshold)
            return false;

        // Mark the cut "likelyGoldenPion" as passed
        cutTracker.MarkCutAsPassed("likelyGoldenPion");
        //std::cout<<"CC1pi: likelyGoldenPion"<<std::endl;

        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetDefaultSelection(const bool skipOldCuts /*= false*/)
{
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        // ########### CC inclusive selection starting from sliceID cuts; all conditions need to be met by at least one particle ###########
        // See docID: 25589-v2 for details

        // track score > 0.85
        {"particleTrackScore", 0.85f},
        // vertex distance < 4 cm
        {"particleVertexDistance", 4.0f},
        // generation = 2 (direct neutrino daughter)
        {"particleGeneration"},
        // track length > 20 cm
        {"particleTrackLength", 20.0f}, // ATTN: This analysis and pelee use SCE corrected lengths numucc does not
        // Proton χ2 > 60
        {"particleProtonChi2", 60.0f},
        // Muon χ2 < 30
        {"particleMuonChi2", 30.0f},
        // Proton χ2/Muon χ2 > 7
        {"particleProtonChi2OverMuonChi2", 7.0f},

        // ########### CC inclusive selection continued; conditions now apply to events ###########
        // Pandora PDG code = 14 (reconstructed object with the most hits is track-like).
        {"pandoraNuPDGIsNumu"},
        // Start vertex of all daughters is at least 10 cm from all borders
        {"daughterVerticesContained", 10.0f},
        // Neutrino vertex is in the fiducial volume.
        {"nuVertexFiducial"},
        // Flash χ2 < 10 OR topological score > 0.25
        {"topologicalOrFlashMatch"},
        // Topological score > 0.06
        {"topologicalScoreCC", 0.06f},

        // {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"2NonProtons", -0.06f},
        {"pionHasValiddEdx", 1.0f},
        {"pionNotInGap"},
        {"muonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f},
        {"likelyGoldenPion", -0.03f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)},
        {"goldenPion", std::make_shared<BDTHelper::BDT>("goldenPion", BDTHelper::GoldenPionBDTFeatureNames)}
    };

    // Define the actual selection logic
    const auto logic = [skipOldCuts](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        const auto &recoParticles = pEvent->reco.particles; // Use all neutrino downstream particles 
        std::vector<bool> selectedParticles(recoParticles.size(), true);
        
        // ----------------------------------------------------------------------------------
        // particleTrackScore
        // ----------------------------------------------------------------------------------
        
        // Get the track score cut value
        const auto minTrackScore = cutTracker.GetCutValue("particleTrackScore");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(recoParticles.at(i)))
            {
                selectedParticles.at(i) = false;
                continue;
            }

            // Skip particles that didn't pass the CC inclusive selection
            // if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).trackScore.IsSet() || recoParticles.at(i).trackScore() < minTrackScore) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackScore" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackScore");

        // ----------------------------------------------------------------------------------
        // particleVertexDistance
        // ----------------------------------------------------------------------------------

        // Get the vertex distance cut value
        const auto maxParticleVertexDist = cutTracker.GetCutValue("particleVertexDistance");
        const auto maxParticleVertexDist2 = maxParticleVertexDist*maxParticleVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertexCC1 = pEvent->reco.nuVertexNoSCC();
        for (unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(recoParticles.at(i).vertexX(), recoParticles.at(i).vertexY(), recoParticles.at(i).vertexZ());
            const auto vertexDist2 = (start - recoVertexCC1).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxParticleVertexDist2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleVertexDistance" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleVertexDistance");

        // ----------------------------------------------------------------------------------
        // particleGeneration
        // ----------------------------------------------------------------------------------

        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).generation() != 2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleGeneration" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleGeneration");

        // ----------------------------------------------------------------------------------
        // particleTrackLength
        // ----------------------------------------------------------------------------------

        // Get the track length cut value
        const auto minTrackLength = cutTracker.GetCutValue("particleTrackLength");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() < minTrackLength) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackLength" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackLength");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 cut value
        const auto minProtonChi2 = cutTracker.GetCutValue("particleProtonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardProtonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardProtonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).chi2ForwardProtonW.IsSet() || recoParticles.at(i).chi2ForwardProtonW() < minProtonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2");

        // ----------------------------------------------------------------------------------
        // particleMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the muon chi2 cut value
        const auto maxMuonChi2 = cutTracker.GetCutValue("particleMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardMuonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardMuonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardMuonW() > maxMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleMuonChi2");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2OverMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 over muon chi2 cut value
        const auto minProtonChi2OverMuonChi2 = cutTracker.GetCutValue("particleProtonChi2OverMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardProtonW()/recoParticles.at(i).chi2ForwardMuonW() < minProtonChi2OverMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2OverMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2OverMuonChi2");

        // ----------------------------------------------------------------------------------
        // Not a cut: In case there are multiple muon candidates, we select the longest one 
        // ----------------------------------------------------------------------------------

        // Get the longest particle
        auto ccInclusiveMuonCandidateIndex = 0u;
        auto longestParticleLength = 0.0f;
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() > longestParticleLength)
            {
                longestParticleLength = recoParticles.at(i).range();
                ccInclusiveMuonCandidateIndex = i;
            }
        }

        // ----------------------------------------------------------------------------------
        // pandoraNuPDGIsNumu
        // ----------------------------------------------------------------------------------

        if(!skipOldCuts)
        {
            // Get the pandora neutrino PDG
            const auto pandoraNuPDGIsNumu = pEvent->reco.nuPdgCode();

            // Insist that the pandora PDG is a muon neutrino
            if (pandoraNuPDGIsNumu != 14)
                return false;
        }
        // else
        // {
        //     std::cout<<"Warning: skipping pandoraNuPDGIsNumu cut"<<std::endl;
        // }

        // Mark the cut "pandoraNuPDGIsNumu" as passed
        cutTracker.MarkCutAsPassed("pandoraNuPDGIsNumu");

        // ----------------------------------------------------------------------------------
        // daughterVerticesContained
        // ----------------------------------------------------------------------------------

        // Get the start near vertex cut value
        const auto margin = cutTracker.GetCutValue("daughterVerticesContained");

        // Insist that all particles with a fitted track start near the vertex
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            // if (!(particle.startX.IsSet() && particle.startY.IsSet() && particle.startZ.IsSet()))
            //     continue;
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                continue;
            }

            // Get the distance between the particle's start position and the vertex
            if(!particle.vertexX.IsSet() || !particle.vertexY.IsSet() || !particle.vertexZ.IsSet())
            {
                return false;
            }
            

            const TVector3 vertex(particle.vertexX(), particle.vertexY(), particle.vertexZ());
            
            // Insist that all vertices are contained within the margin
            if (!AnalysisHelper::IsPointWithinMargins(vertex, margin, margin, margin, margin, margin, margin))
                return false;
        }

        // Mark the cut "daughterVerticesContained" as passed
        cutTracker.MarkCutAsPassed("daughterVerticesContained");

        // ----------------------------------------------------------------------------------
        // nuVertexFiducial
        // ----------------------------------------------------------------------------------

        // Get the vertex position
        const auto recoVertexCC = pEvent->reco.nuVertexNoSCC();
        if(!AnalysisHelper::IsFiducial(recoVertexCC))
            return false;
        
        // Mark the cut "nuVertexFiducial" as passed
        cutTracker.MarkCutAsPassed("nuVertexFiducial");

        // ----------------------------------------------------------------------------------
        // topologicalOrFlashMatch
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC1 = 0.25f;
        const auto maxFlashMatchScore = 10.0f;

        // Insist that the topological score is above the cut value or the flash match score is below the cut value
        // Here: when neither is the case, the event is rejected 
        // if(pEvent->reco.flashChi2.IsSet()) std::cout<<"flashChi2: "<<pEvent->reco.flashChi2.IsSet()<<std::endl;
        // else std::cout<<"flashChi2 not set"<<std::endl;

        if(!skipOldCuts)
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1 && pEvent->reco.flashChi2.IsSet() && pEvent->reco.flashChi2() >= maxFlashMatchScore)
                return false;
        }
        else
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1)
                return false;
        }

        // Mark the cut "topologicalOrFlashMatch" as passed
        cutTracker.MarkCutAsPassed("topologicalOrFlashMatch");

        // ----------------------------------------------------------------------------------
        // topologicalScoreCC
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC2 = cutTracker.GetCutValue("topologicalScoreCC");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC2)
            return false;

        // Mark the cut "topologicalScoreCC" as passed
        cutTracker.MarkCutAsPassed("topologicalScoreCC");

        // ##################################################################################
        // Start of pion selection
        // ##################################################################################

        // // ----------------------------------------------------------------------------------
        // // passesCCInclusive
        // // ----------------------------------------------------------------------------------

        // // Insist the event passes the CC inclusive preselection
        // if (!pEvent->reco.passesCCInclusive())
        //     return false;

        // // Mark the cut "passesCCInclusive" as passed
        // cutTracker.MarkCutAsPassed("passesCCInclusive");

        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");
        //std::cout<<"CC1pi: min2Tracks"<<std::endl;

        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;


        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");
        //std::cout<<"CC1pi: max1Uncontained"<<std::endl;

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT, ccInclusiveMuonCandidateIndex);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // ----------------------------------------------------------------------------------
        // 2NonProtons
        // ----------------------------------------------------------------------------------

        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThreshold = cutTracker.GetCutValue("2NonProtons");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        std::vector<unsigned int> pionIndices;

        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
            {
                continue;
            }

            // Assume particles without a track fit are just small protons
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("DefaultSelection - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            if (!hasFeatures)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThreshold)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nPions + nMuons != recoParticles.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen!");

        // Insist that we exacly one pion (i.e we have have 2 non-protons)
        if (nPions != 1)
            return false;

        // Mark the cut "2NonProtons" as passed
        cutTracker.MarkCutAsPassed("2NonProtons");
        //std::cout<<"CC1pi: 2NonProtons"<<std::endl;

        // ----------------------------------------------------------------------------------
        // pionHasValiddEdx
        // ----------------------------------------------------------------------------------

        const auto piondEdxThreshold = cutTracker.GetCutValue("pionHasValiddEdx");

        // Get the pion reco particle
        const auto pionIndex = pionIndices.front();
        const auto &pion = recoParticles.at(pionIndex);

        // Insist that the pion has a valid deEdx measurement (i.e. truncated mean dE/dx is above threshold)
        if (pion.truncatedMeandEdx() < piondEdxThreshold)
            return false;

        // Mark the cut "pionHasValiddEdx" as passed
        cutTracker.MarkCutAsPassed("pionHasValiddEdx");
        //std::cout<<"CC1pi: pionHasValiddEdx"<<std::endl;

        // ----------------------------------------------------------------------------------
        // pionNotInGap
        // ----------------------------------------------------------------------------------

        // Sanity check that our muon and pion are not the same particle
        if (muonIndex == pionIndex)
            throw std::logic_error("DefaultSelection - The muon and the pion candidates are the same particle. This shouldn't happen!");

        // Insist that the pion has at least one hit in each view (i.e. not in a gap)
        if (pion.nHitsU() == 0 || pion.nHitsV() == 0 || pion.nHitsW() == 0)
            return false;

        // Mark the cut "pionNotInGap" as passed
        cutTracker.MarkCutAsPassed("pionNotInGap");
        //std::cout<<"CC1pi: pionNotInGap"<<std::endl;

        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");
        //std::cout<<"CC1pi: muonNotInGap"<<std::endl;

        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------

        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");

        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(pionDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");
        //std::cout<<"CC1pi: openingAngle"<<std::endl;

        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");
        //std::cout<<"CC1pi: topologicalScore"<<std::endl;
        
        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------

        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                continue;
            }

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
            {
                return false;
            }
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");
        //std::cout<<"CC1pi: startNearVertex"<<std::endl;

        // ----------------------------------------------------------------------------------
        // likelyGoldenPion
        // ----------------------------------------------------------------------------------

        // Get the likely golden pion cut value
        const auto goldenPionBDTThreshold = cutTracker.GetCutValue("likelyGoldenPion");

        // Get the golden pion BDT
        auto &pGoldenPionBDT = bdtMap.at("goldenPion");

        // Get the features of the pion
        std::vector<float> features;
        if (!BDTHelper::GetBDTFeatures(pion, BDTHelper::GoldenPionBDTFeatureNames, features))
            throw std::logic_error("DefaultSelection - Can't get golden pion BDT features of pion candidate");

        // Insist that the BDT response is greater than the cut value to identify the pion as a golden pion
        const auto bdtResponse = pGoldenPionBDT->GetResponse(features);
        if (bdtResponse <= goldenPionBDTThreshold)
            return false;

        // Mark the cut "likelyGoldenPion" as passed
        cutTracker.MarkCutAsPassed("likelyGoldenPion");
        //std::cout<<"CC1pi: likelyGoldenPion"<<std::endl;

        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetDefaultSelection2(const bool skipOldCuts /*= false*/)
{
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        // ########### CC inclusive selection starting from sliceID cuts; all conditions need to be met by at least one particle ###########
        // See docID: 25589-v2 for details

        // track score > 0.85
        {"particleTrackScore", 0.85f},
        // vertex distance < 4 cm
        {"particleVertexDistance", 4.0f},
        // generation = 2 (direct neutrino daughter)
        {"particleGeneration"},
        // track length > 20 cm
        {"particleTrackLength", 20.0f}, // ATTN: This analysis and pelee use SCE corrected lengths numucc does not
        // Proton χ2 > 60
        {"particleProtonChi2", 60.0f},
        // Muon χ2 < 30
        {"particleMuonChi2", 30.0f},
        // Proton χ2/Muon χ2 > 7
        {"particleProtonChi2OverMuonChi2", 7.0f},

        // ########### CC inclusive selection continued; conditions now apply to events ###########
        // Pandora PDG code = 14 (reconstructed object with the most hits is track-like).
        {"pandoraNuPDGIsNumu"},
        // Start vertex of all daughters is at least 10 cm from all borders
        {"daughterVerticesContained", 10.0f},
        // Neutrino vertex is in the fiducial volume.
        {"nuVertexFiducial"},
        // Flash χ2 < 10 OR topological score > 0.25
        {"topologicalOrFlashMatch"},
        // Topological score > 0.06
        {"topologicalScoreCC", 0.06f},

        // {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"2NonProtons", -0.06f},
        {"pionHasValiddEdx", 1.0f},
        {"pionNotInGap"},
        {"muonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f},
        {"likelyGoldenPion", -0.03f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)},
        {"goldenPion", std::make_shared<BDTHelper::BDT>("goldenPion", BDTHelper::GoldenPionBDTFeatureNames)}
    };

    // Define the actual selection logic
    const auto logic = [skipOldCuts](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        const auto &recoParticles = pEvent->reco.particles; // Use all neutrino downstream particles 
        std::vector<bool> selectedParticles(recoParticles.size(), true);
        // ----------------------------------------------------------------------------------
        // particleTrackScore
        // ----------------------------------------------------------------------------------
        
        // Get the track score cut value
        const auto minTrackScore = cutTracker.GetCutValue("particleTrackScore");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(recoParticles.at(i)))
            {
                selectedParticles.at(i) = false;
                continue;
            }

            // Skip particles that didn't pass the CC inclusive selection
            if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).trackScore.IsSet() || recoParticles.at(i).trackScore() < minTrackScore) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackScore" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackScore");

        // ----------------------------------------------------------------------------------
        // particleVertexDistance
        // ----------------------------------------------------------------------------------

        // Get the vertex distance cut value
        const auto maxParticleVertexDist = cutTracker.GetCutValue("particleVertexDistance");
        const auto maxParticleVertexDist2 = maxParticleVertexDist*maxParticleVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertexCC1 = pEvent->reco.nuVertexNoSCC();
        for (unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(recoParticles.at(i).vertexX(), recoParticles.at(i).vertexY(), recoParticles.at(i).vertexZ());
            const auto vertexDist2 = (start - recoVertexCC1).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxParticleVertexDist2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleVertexDistance" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleVertexDistance");

        // ----------------------------------------------------------------------------------
        // particleGeneration
        // ----------------------------------------------------------------------------------
        std::vector<unsigned int> daughterRecoParticleIndices{};
        for(unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            // if(!selectedParticles.at(i)) continue; // This saves some time but here we want to iterate over all particles to also crete daughterRecoParticleIndices
            if(recoParticles.at(i).generation() == 2)
            {
                daughterRecoParticleIndices.push_back(i);
            }
            else
            {
                selectedParticles.at(i) = false;
            }
        }
        // Mark the cut "particleGeneration" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleGeneration");

        // ----------------------------------------------------------------------------------
        // particleTrackLength
        // ----------------------------------------------------------------------------------

        // Get the track length cut value
        const auto minTrackLength = cutTracker.GetCutValue("particleTrackLength");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() < minTrackLength) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackLength" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackLength");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 cut value
        const auto minProtonChi2 = cutTracker.GetCutValue("particleProtonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardProtonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardProtonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).chi2ForwardProtonW.IsSet() || recoParticles.at(i).chi2ForwardProtonW() < minProtonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2");

        // ----------------------------------------------------------------------------------
        // particleMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the muon chi2 cut value
        const auto maxMuonChi2 = cutTracker.GetCutValue("particleMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardMuonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardMuonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardMuonW() > maxMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleMuonChi2");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2OverMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 over muon chi2 cut value
        const auto minProtonChi2OverMuonChi2 = cutTracker.GetCutValue("particleProtonChi2OverMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardProtonW()/recoParticles.at(i).chi2ForwardMuonW() < minProtonChi2OverMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2OverMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2OverMuonChi2");

        // ----------------------------------------------------------------------------------
        // Not a cut: In case there are multiple muon candidates, we select the longest one 
        // ----------------------------------------------------------------------------------

        // Get the longest particle
        auto ccInclusiveMuonCandidateIndex = 0u;
        auto longestParticleLength = 0.0f;
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() > longestParticleLength)
            {
                longestParticleLength = recoParticles.at(i).range();
                ccInclusiveMuonCandidateIndex = i;
            }
        }

        // ----------------------------------------------------------------------------------
        // pandoraNuPDGIsNumu
        // ----------------------------------------------------------------------------------

        if(!skipOldCuts)
        {
            // Get the pandora neutrino PDG
            const auto pandoraNuPDGIsNumu = pEvent->reco.nuPdgCode();

            // Insist that the pandora PDG is a muon neutrino
            if (pandoraNuPDGIsNumu != 14)
                return false;
        }
        // else
        // {
        //     std::cout<<"Warning: skipping pandoraNuPDGIsNumu cut"<<std::endl;
        // }

        // Mark the cut "pandoraNuPDGIsNumu" as passed
        cutTracker.MarkCutAsPassed("pandoraNuPDGIsNumu");

        // ----------------------------------------------------------------------------------
        // daughterVerticesContained
        // ----------------------------------------------------------------------------------

        // Get the start near vertex cut value
        const auto margin = cutTracker.GetCutValue("daughterVerticesContained");

        // Insist that all particles with a fitted track start near the vertex
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            // if (!(particle.startX.IsSet() && particle.startY.IsSet() && particle.startZ.IsSet()))
            //     continue;
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                continue;
            }

            // Get the distance between the particle's start position and the vertex
            if(!particle.vertexX.IsSet() || !particle.vertexY.IsSet() || !particle.vertexZ.IsSet())
            {
                return false;
            }
            

            const TVector3 vertex(particle.vertexX(), particle.vertexY(), particle.vertexZ());
            
            // Insist that all vertices are contained within the margin
            if (!AnalysisHelper::IsPointWithinMargins(vertex, margin, margin, margin, margin, margin, margin))
            {
                return false;
            }
        }

        // Mark the cut "daughterVerticesContained" as passed
        cutTracker.MarkCutAsPassed("daughterVerticesContained");

        // ----------------------------------------------------------------------------------
        // nuVertexFiducial
        // ----------------------------------------------------------------------------------

        // Get the vertex position
        const auto recoVertexCC = pEvent->reco.nuVertexNoSCC();
        if(!AnalysisHelper::IsFiducial(recoVertexCC))
            return false;
        
        // Mark the cut "nuVertexFiducial" as passed
        cutTracker.MarkCutAsPassed("nuVertexFiducial");

        // ----------------------------------------------------------------------------------
        // topologicalOrFlashMatch
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC1 = 0.25f;
        const auto maxFlashMatchScore = 10.0f;

        // Insist that the topological score is above the cut value or the flash match score is below the cut value
        // Here: when neither is the case, the event is rejected 
        // if(pEvent->reco.flashChi2.IsSet()) std::cout<<"flashChi2: "<<pEvent->reco.flashChi2.IsSet()<<std::endl;
        // else std::cout<<"flashChi2 not set"<<std::endl;

        if(!skipOldCuts)
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1 && pEvent->reco.flashChi2.IsSet() && pEvent->reco.flashChi2() >= maxFlashMatchScore)
                return false;
        }
        else
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1)
                return false;
        }
        // else
        // {
        //     std::cout<<"Warning: skipping topologicalOrFlashMatch cut"<<std::endl;
        // }

        // Mark the cut "topologicalOrFlashMatch" as passed
        cutTracker.MarkCutAsPassed("topologicalOrFlashMatch");

        // ----------------------------------------------------------------------------------
        // topologicalScoreCC
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC2 = cutTracker.GetCutValue("topologicalScoreCC");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC2)
            return false;

        // Mark the cut "topologicalScoreCC" as passed
        cutTracker.MarkCutAsPassed("topologicalScoreCC");

        // ##################################################################################
        // Start of pion selection
        // ##################################################################################

        // // ----------------------------------------------------------------------------------
        // // passesCCInclusive
        // // ----------------------------------------------------------------------------------

        // // Insist the event passes the CC inclusive preselection
        // if (!pEvent->reco.passesCCInclusive())
        //     return false;

        // // Mark the cut "passesCCInclusive" as passed
        // cutTracker.MarkCutAsPassed("passesCCInclusive");

        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG min2Tracks"<<std::endl;

        // auto pionSectionRecoParticles = recoParticles;
        // for (int i = pionSectionRecoParticles.size() - 1; i >= 0; --i) // Loop backwards so we can erase elements
        // {
        //     const auto &particle = pionSectionRecoParticles.at(i);
        //     if (particle.generation() != 2)
        //     {
        //         pionSectionRecoParticles.erase(pionSectionRecoParticles.begin() + i);
        //     }
        // }

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        for (const auto i:daughterRecoParticleIndices)//unsigned int i = 0; i < pionSectionRecoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");
        //std::cout<<"CC1pi: min2Tracks"<<std::endl;

        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG max1Uncontained"<<std::endl;

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;


        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");
        //std::cout<<"CC1pi: max1Uncontained"<<std::endl;

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT, ccInclusiveMuonCandidateIndex);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // ----------------------------------------------------------------------------------
        // 2NonProtons
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG 2NonProtons"<<std::endl;

        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThreshold = cutTracker.GetCutValue("2NonProtons");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        std::vector<unsigned int> pionIndices;

        for (const auto i:daughterRecoParticleIndices)//unsigned int i = 0; i < pionSectionRecoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
            {
                continue;
            }

            // Assume particles without a track fit are just small protons
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("GetDefaultSelection2 - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            if (!hasFeatures)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThreshold)
            {
                nProtons++;
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nPions + nMuons != daughterRecoParticleIndices.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen! # protons: " + std::to_string(nProtons) + ", # pions: " + std::to_string(nPions) + ", # muons: " + std::to_string(nMuons) + ", # particles: " + std::to_string(daughterRecoParticleIndices.size())); 

        // Insist that we exacly one pion (i.e we have have 2 non-protons)
        if (nPions != 1)
            return false;

        // Mark the cut "2NonProtons" as passed
        cutTracker.MarkCutAsPassed("2NonProtons");
        //std::cout<<"CC1pi: 2NonProtons"<<std::endl;

        // ----------------------------------------------------------------------------------
        // pionHasValiddEdx
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG pionHasValiddEdx"<<std::endl;

        const auto piondEdxThreshold = cutTracker.GetCutValue("pionHasValiddEdx");

        // Get the pion reco particle
        const auto pionIndex = pionIndices.front();
        const auto &pion = recoParticles.at(pionIndex);

        // Insist that the pion has a valid deEdx measurement (i.e. truncated mean dE/dx is above threshold)
        if (pion.truncatedMeandEdx() < piondEdxThreshold)
            return false;

        // Mark the cut "pionHasValiddEdx" as passed
        cutTracker.MarkCutAsPassed("pionHasValiddEdx");
        //std::cout<<"CC1pi: pionHasValiddEdx"<<std::endl;

        // ----------------------------------------------------------------------------------
        // pionNotInGap
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG pionNotInGap"<<std::endl;

        // Sanity check that our muon and pion are not the same particle
        if (muonIndex == pionIndex)
            throw std::logic_error("DefaultSelection - The muon and the pion candidates are the same particle. This shouldn't happen!");

        // Insist that the pion has at least one hit in each view (i.e. not in a gap)
        if (pion.nHitsU() == 0 || pion.nHitsV() == 0 || pion.nHitsW() == 0)
            return false;

        // Mark the cut "pionNotInGap" as passed
        cutTracker.MarkCutAsPassed("pionNotInGap");
        //std::cout<<"CC1pi: pionNotInGap"<<std::endl;

        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG muonNotInGap"<<std::endl;

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");
        //std::cout<<"CC1pi: muonNotInGap"<<std::endl;

        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG openingAngle"<<std::endl;

        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");
        
        // if( !(muon.directionX.IsSet() && muon.directionY.IsSet() && muon.directionZ.IsSet() && pion.directionX.IsSet() && pion.directionY.IsSet() && pion.directionZ.IsSet()) ) return false;
        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(pionDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");
        //std::cout<<"CC1pi: openingAngle"<<std::endl;

        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG topologicalScore"<<std::endl;

        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");
        //std::cout<<"CC1pi: topologicalScore"<<std::endl;
        
        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG startNearVertex"<<std::endl;

        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto i:daughterRecoParticleIndices)
        {
            const auto &particle = recoParticles.at(i);
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                continue;
            }

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
            {
                return false;
            }
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");
        //std::cout<<"CC1pi: startNearVertex"<<std::endl;

        // ----------------------------------------------------------------------------------
        // likelyGoldenPion
        // ----------------------------------------------------------------------------------
        // std::cout<<"DEBUG likelyGoldenPion"<<std::endl;

        // Get the likely golden pion cut value
        const auto goldenPionBDTThreshold = cutTracker.GetCutValue("likelyGoldenPion");

        // Get the golden pion BDT
        auto &pGoldenPionBDT = bdtMap.at("goldenPion");

        // Get the features of the pion
        std::vector<float> features;
        if (!BDTHelper::GetBDTFeatures(pion, BDTHelper::GoldenPionBDTFeatureNames, features))
            throw std::logic_error("DefaultSelection - Can't get golden pion BDT features of pion candidate");

        // Insist that the BDT response is greater than the cut value to identify the pion as a golden pion
        const auto bdtResponse = pGoldenPionBDT->GetResponse(features);
        if (bdtResponse <= goldenPionBDTThreshold)
            return false;

        // Mark the cut "likelyGoldenPion" as passed
        cutTracker.MarkCutAsPassed("likelyGoldenPion");
        //std::cout<<"CC1pi: likelyGoldenPion"<<std::endl;

        // We passed all cuts!
        return true;
    };

    // std::cout<<"DEBUG done with Default2 logic"<<std::endl;
    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCCInclusiveSelection2(const bool skipOldCuts /*= false*/)
{
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        // ########### CC inclusive selection starting from sliceID cuts; all conditions need to be met by at least one particle ###########
        // See docID: 25589-v2 for details

        // track score > 0.85
        {"particleTrackScore", 0.85f},
        // vertex distance < 4 cm
        {"particleVertexDistance", 4.0f},
        // generation = 2 (direct neutrino daughter)
        {"particleGeneration"},
        // track length > 20 cm
        {"particleTrackLength", 20.0f}, // ATTN: This analysis and pelee use SCE corrected lengths numucc does not
        // Proton χ2 > 60
        {"particleProtonChi2", 60.0f},
        // Muon χ2 < 30
        {"particleMuonChi2", 30.0f},
        // Proton χ2/Muon χ2 > 7
        {"particleProtonChi2OverMuonChi2", 7.0f},

        // ########### CC inclusive selection continued; conditions now apply to events ###########
        // Pandora PDG code = 14 (reconstructed object with the most hits is track-like).
        {"pandoraNuPDGIsNumu"},
        // Start vertex of all daughters is at least 10 cm from all borders
        {"daughterVerticesContained", 10.0f},
        // Neutrino vertex is in the fiducial volume.
        {"nuVertexFiducial"},
        // Flash χ2 < 10 OR topological score > 0.25
        {"topologicalOrFlashMatch"},
        // Topological score > 0.06
        {"topologicalScoreCC", 0.06f},
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {};

    // Define the actual selection logic
    const auto logic = [skipOldCuts](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        const auto &recoParticles = pEvent->reco.particles; // Use all neutrino downstream particles 
        std::vector<bool> selectedParticles(recoParticles.size(), true);
        // ----------------------------------------------------------------------------------
        // particleTrackScore
        // ----------------------------------------------------------------------------------
        
        // Get the track score cut value
        const auto minTrackScore = cutTracker.GetCutValue("particleTrackScore");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(recoParticles.at(i)))
            {
                selectedParticles.at(i) = false;
                continue;
            }

            // Skip particles that didn't pass the CC inclusive selection
            if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).trackScore.IsSet() || recoParticles.at(i).trackScore() < minTrackScore) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackScore" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackScore");

        // ----------------------------------------------------------------------------------
        // particleVertexDistance
        // ----------------------------------------------------------------------------------

        // Get the vertex distance cut value
        const auto maxParticleVertexDist = cutTracker.GetCutValue("particleVertexDistance");
        const auto maxParticleVertexDist2 = maxParticleVertexDist*maxParticleVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertexCC1 = pEvent->reco.nuVertexNoSCC();
        for (unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(recoParticles.at(i).vertexX(), recoParticles.at(i).vertexY(), recoParticles.at(i).vertexZ());
            const auto vertexDist2 = (start - recoVertexCC1).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxParticleVertexDist2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleVertexDistance" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleVertexDistance");

        // ----------------------------------------------------------------------------------
        // particleGeneration
        // ----------------------------------------------------------------------------------
        std::vector<unsigned int> daughterRecoParticleIndices{};
        for(unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            // if(!selectedParticles.at(i)) continue; // This saves some time but here we want to iterate over all particles to also crete daughterRecoParticleIndices
            if(recoParticles.at(i).generation() == 2)
            {
                daughterRecoParticleIndices.push_back(i);
            }
            else
            {
                selectedParticles.at(i) = false;
            }
        }
        // Mark the cut "particleGeneration" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleGeneration");

        // ----------------------------------------------------------------------------------
        // particleTrackLength
        // ----------------------------------------------------------------------------------

        // Get the track length cut value
        const auto minTrackLength = cutTracker.GetCutValue("particleTrackLength");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() < minTrackLength) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackLength" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackLength");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 cut value
        const auto minProtonChi2 = cutTracker.GetCutValue("particleProtonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardProtonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardProtonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).chi2ForwardProtonW.IsSet() || recoParticles.at(i).chi2ForwardProtonW() < minProtonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2");

        // ----------------------------------------------------------------------------------
        // particleMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the muon chi2 cut value
        const auto maxMuonChi2 = cutTracker.GetCutValue("particleMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardMuonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardMuonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardMuonW() > maxMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleMuonChi2");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2OverMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 over muon chi2 cut value
        const auto minProtonChi2OverMuonChi2 = cutTracker.GetCutValue("particleProtonChi2OverMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardProtonW()/recoParticles.at(i).chi2ForwardMuonW() < minProtonChi2OverMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2OverMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2OverMuonChi2");

        // ----------------------------------------------------------------------------------
        // Not a cut: In case there are multiple muon candidates, we select the longest one 
        // ----------------------------------------------------------------------------------

        // Get the longest particle
        auto ccInclusiveMuonCandidateIndex = 0u;
        auto longestParticleLength = 0.0f;
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() > longestParticleLength)
            {
                longestParticleLength = recoParticles.at(i).range();
                ccInclusiveMuonCandidateIndex = i;
            }
        }

        // ----------------------------------------------------------------------------------
        // pandoraNuPDGIsNumu
        // ----------------------------------------------------------------------------------

        if(!skipOldCuts)
        {
            // Get the pandora neutrino PDG
            const auto pandoraNuPDGIsNumu = pEvent->reco.nuPdgCode();

            // Insist that the pandora PDG is a muon neutrino
            if (pandoraNuPDGIsNumu != 14)
                return false;
        }
        // else
        // {
        //     std::cout<<"Warning: skipping pandoraNuPDGIsNumu cut"<<std::endl;
        // }

        // Mark the cut "pandoraNuPDGIsNumu" as passed
        cutTracker.MarkCutAsPassed("pandoraNuPDGIsNumu");

        // ----------------------------------------------------------------------------------
        // daughterVerticesContained
        // ----------------------------------------------------------------------------------

        // Get the start near vertex cut value
        const auto margin = cutTracker.GetCutValue("daughterVerticesContained");

        // Insist that all particles with a fitted track start near the vertex
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            // if (!(particle.startX.IsSet() && particle.startY.IsSet() && particle.startZ.IsSet()))
            //     continue;
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                continue;
            }

            // Get the distance between the particle's start position and the vertex
            if(!particle.vertexX.IsSet() || !particle.vertexY.IsSet() || !particle.vertexZ.IsSet())
            {
                return false;
            }
            

            const TVector3 vertex(particle.vertexX(), particle.vertexY(), particle.vertexZ());
            
            // Insist that all vertices are contained within the margin
            if (!AnalysisHelper::IsPointWithinMargins(vertex, margin, margin, margin, margin, margin, margin))
            {
                return false;
            }
        }

        // Mark the cut "daughterVerticesContained" as passed
        cutTracker.MarkCutAsPassed("daughterVerticesContained");

        // ----------------------------------------------------------------------------------
        // nuVertexFiducial
        // ----------------------------------------------------------------------------------

        // Get the vertex position
        const auto recoVertexCC = pEvent->reco.nuVertexNoSCC();
        if(!AnalysisHelper::IsFiducial(recoVertexCC))
            return false;
        
        // Mark the cut "nuVertexFiducial" as passed
        cutTracker.MarkCutAsPassed("nuVertexFiducial");

        // ----------------------------------------------------------------------------------
        // topologicalOrFlashMatch
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC1 = 0.25f;
        const auto maxFlashMatchScore = 10.0f;

        // Insist that the topological score is above the cut value or the flash match score is below the cut value
        // Here: when neither is the case, the event is rejected 
        // if(pEvent->reco.flashChi2.IsSet()) std::cout<<"flashChi2: "<<pEvent->reco.flashChi2.IsSet()<<std::endl;
        // else std::cout<<"flashChi2 not set"<<std::endl;

        if(!skipOldCuts)
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1 && pEvent->reco.flashChi2.IsSet() && pEvent->reco.flashChi2() >= maxFlashMatchScore)
                return false;
        }
        else
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1)
                return false;
        }

        // Mark the cut "topologicalOrFlashMatch" as passed
        cutTracker.MarkCutAsPassed("topologicalOrFlashMatch");

        // ----------------------------------------------------------------------------------
        // topologicalScoreCC
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC2 = cutTracker.GetCutValue("topologicalScoreCC");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC2)
            return false;

        // Mark the cut "topologicalScoreCC" as passed
        cutTracker.MarkCutAsPassed("topologicalScoreCC");

        // We passed all cuts!
        return true;
    };

    // std::cout<<"DEBUG done with Default2 logic"<<std::endl;
    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCC0piSelectionModifiedPeLEE(const float muonLikeProtonValue /*= -0.4f*/, const float barelyResemblingProtonValue /*= 0.1f*/, const bool skipOldCuts /*= false*/)
{
    // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 0\n";
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        // ########### CC inclusive selection starting from sliceID cuts; all conditions need to be met by at least one particle ###########
        // See docID: 25589-v2 for details

        // track score > 0.85
        {"particleTrackScore", 0.85f},
        // vertex distance < 4 cm
        {"particleVertexDistance", 4.0f},
        // generation = 2 (direct neutrino daughter)
        {"particleGeneration"},
        // track length > 20 cm
        {"particleTrackLength", 20.0f}, // ATTN: This analysis and pelee use SCE corrected lengths numucc does not
        // Proton χ2 > 60
        {"particleProtonChi2", 60.0f},
        // Muon χ2 < 30
        {"particleMuonChi2", 30.0f},
        // Proton χ2/Muon χ2 > 7
        {"particleProtonChi2OverMuonChi2", 7.0f},

        // ########### CC inclusive selection continued; conditions now apply to events ###########
        // Pandora PDG code = 14 (reconstructed object with the most hits is track-like).
        {"pandoraNuPDGIsNumu"},
        // Start vertex of all daughters is at least 10 cm from all borders
        {"daughterVerticesContained", 10.0f},
        // Neutrino vertex is in the fiducial volume.
        {"nuVertexFiducial"},
        // Flash χ2 < 10 OR topological score > 0.25
        {"topologicalOrFlashMatch"},
        // Topological score > 0.06
        {"topologicalScoreCC", 0.06f},

        // {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"1NonProton", -0.06f},
        {"atLeast1Proton"},
        {"muonLikeProton", muonLikeProtonValue}, //todo: why?
        {"barelyResemblingProton", barelyResemblingProtonValue},
        {"protonHasValiddEdx", 1.0f},
        {"muonNotInGap"},
        {"protonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f},
        {"likelyGoldenProton", -0.03f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)},
        {"goldenPion", std::make_shared<BDTHelper::BDT>("goldenPion", BDTHelper::GoldenPionBDTFeatureNames)}
    };

    // Define the actual selection logic
    const auto logic = [skipOldCuts](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        const auto &recoParticles = pEvent->reco.particles; // Use all neutrino downstream particles
        std::vector<bool> selectedParticles(recoParticles.size(), true);
        // ----------------------------------------------------------------------------------
        // particleTrackScore
        // ----------------------------------------------------------------------------------

        // Get the track score cut value
        const auto minTrackScore = cutTracker.GetCutValue("particleTrackScore");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(recoParticles.at(i)))
            {
                selectedParticles.at(i) = false;
                continue;
            }

            // Skip particles that didn't pass the CC inclusive selection
            if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).trackScore.IsSet() || recoParticles.at(i).trackScore() < minTrackScore) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackScore" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackScore");

        // ----------------------------------------------------------------------------------
        // particleVertexDistance
        // ----------------------------------------------------------------------------------

        // Get the vertex distance cut value
        const auto maxParticleVertexDist = cutTracker.GetCutValue("particleVertexDistance");
        const auto maxParticleVertexDist2 = maxParticleVertexDist*maxParticleVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertexCC1 = pEvent->reco.nuVertexNoSCC();
        for (unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(recoParticles.at(i).vertexX(), recoParticles.at(i).vertexY(), recoParticles.at(i).vertexZ());
            const auto vertexDist2 = (start - recoVertexCC1).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxParticleVertexDist2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleVertexDistance" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleVertexDistance");

        // ----------------------------------------------------------------------------------
        // particleGeneration
        // ----------------------------------------------------------------------------------
        std::vector<unsigned int> daughterRecoParticleIndices{};
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            // if(!selectedParticles.at(i)) continue; // This saves some time but here we want to iterate over all particles to also crete daughterRecoParticleIndices
            if(recoParticles.at(i).generation() == 2)
            {
                daughterRecoParticleIndices.push_back(i);
            }
            else
            {
                selectedParticles.at(i) = false;
            }
        }
        // Mark the cut "particleGeneration" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleGeneration");

        // ----------------------------------------------------------------------------------
        // particleTrackLength
        // ----------------------------------------------------------------------------------

        // Get the track length cut value
        const auto minTrackLength = cutTracker.GetCutValue("particleTrackLength");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() < minTrackLength) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleTrackLength" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleTrackLength");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 cut value
        const auto minProtonChi2 = cutTracker.GetCutValue("particleProtonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardProtonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardProtonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(!recoParticles.at(i).chi2ForwardProtonW.IsSet() || recoParticles.at(i).chi2ForwardProtonW() < minProtonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2");

        // ----------------------------------------------------------------------------------
        // particleMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the muon chi2 cut value
        const auto maxMuonChi2 = cutTracker.GetCutValue("particleMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // if(recoParticles.at(i).chi2ForwardMuonW.IsSet())std::cout<<recoParticles.at(i).chi2ForwardMuonW();
            // else std::cout<<"not set";
            // std::cout<<"("<<selectedParticles.at(i)<<")"<<std::endl;
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardMuonW() > maxMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleMuonChi2");

        // ----------------------------------------------------------------------------------
        // particleProtonChi2OverMuonChi2
        // ----------------------------------------------------------------------------------

        // Get the proton chi2 over muon chi2 cut value
        const auto minProtonChi2OverMuonChi2 = cutTracker.GetCutValue("particleProtonChi2OverMuonChi2");
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).chi2ForwardProtonW()/recoParticles.at(i).chi2ForwardMuonW() < minProtonChi2OverMuonChi2) selectedParticles.at(i) = false;
        }
        // Mark the cut "particleProtonChi2OverMuonChi2" as passed as long as there is at least one particle that passes the cut
        if(!std::any_of(selectedParticles.begin(), selectedParticles.end(), [](bool value) {return value;}))
            return false;
        cutTracker.MarkCutAsPassed("particleProtonChi2OverMuonChi2");

        // ----------------------------------------------------------------------------------
        // Not a cut: In case there are multiple muon candidates, we select the longest one
        // ----------------------------------------------------------------------------------

        // Get the longest particle
        auto ccInclusiveMuonCandidateIndex = 0u;
        auto longestParticleLength = 0.0f;
        for(unsigned int i = 0; i< recoParticles.size(); ++i)
        {
            // Skip particles that didn't pass the previous cuts
            if(!selectedParticles.at(i)) continue;
            if(recoParticles.at(i).range() > longestParticleLength)
            {
                longestParticleLength = recoParticles.at(i).range();
                ccInclusiveMuonCandidateIndex = i;
            }
        }

        // ----------------------------------------------------------------------------------
        // pandoraNuPDGIsNumu
        // ----------------------------------------------------------------------------------

        if(!skipOldCuts)
        {
            // Get the pandora neutrino PDG
            const auto pandoraNuPDGIsNumu = pEvent->reco.nuPdgCode();

            // Insist that the pandora PDG is a muon neutrino
            if (pandoraNuPDGIsNumu != 14)
                return false;
        }
        // else
        // {
        //     std::cout<<"Warning: skipping pandoraNuPDGIsNumu cut"<<std::endl;
        // }

        // Mark the cut "pandoraNuPDGIsNumu" as passed
        cutTracker.MarkCutAsPassed("pandoraNuPDGIsNumu");

        // ----------------------------------------------------------------------------------
        // daughterVerticesContained
        // ----------------------------------------------------------------------------------

        // Get the start near vertex cut value
        const auto margin = cutTracker.GetCutValue("daughterVerticesContained");

        // Insist that all particles with a fitted track start near the vertex
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            // if (!(particle.startX.IsSet() && particle.startY.IsSet() && particle.startZ.IsSet()))
            //     continue;
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                continue;
            }

            // Get the distance between the particle's start position and the vertex
            if(!particle.vertexX.IsSet() || !particle.vertexY.IsSet() || !particle.vertexZ.IsSet())
            {
                return false;
            }


            const TVector3 vertex(particle.vertexX(), particle.vertexY(), particle.vertexZ());

            // Insist that all vertices are contained within the margin
            if (!AnalysisHelper::IsPointWithinMargins(vertex, margin, margin, margin, margin, margin, margin))
            {
                return false;
            }
        }

        // Mark the cut "daughterVerticesContained" as passed
        cutTracker.MarkCutAsPassed("daughterVerticesContained");

        // ----------------------------------------------------------------------------------
        // nuVertexFiducial
        // ----------------------------------------------------------------------------------

        // Get the vertex position
        const auto recoVertexCC = pEvent->reco.nuVertexNoSCC();
        if(!AnalysisHelper::IsFiducial(recoVertexCC))
            return false;

        // Mark the cut "nuVertexFiducial" as passed
        cutTracker.MarkCutAsPassed("nuVertexFiducial");

        // ----------------------------------------------------------------------------------
        // topologicalOrFlashMatch
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC1 = 0.25f;
        const auto maxFlashMatchScore = 10.0f;

        // Insist that the topological score is above the cut value or the flash match score is below the cut value
        // Here: when neither is the case, the event is rejected
        // if(pEvent->reco.flashChi2.IsSet()) std::cout<<"flashChi2: "<<pEvent->reco.flashChi2.IsSet()<<std::endl;
        // else std::cout<<"flashChi2 not set"<<std::endl;

        if(!skipOldCuts)
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1 && pEvent->reco.flashChi2.IsSet() && pEvent->reco.flashChi2() >= maxFlashMatchScore)
                return false;
        }
        else
        {
            if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC1)
                return false;
        }

        // Mark the cut "topologicalOrFlashMatch" as passed
        cutTracker.MarkCutAsPassed("topologicalOrFlashMatch");

        // ----------------------------------------------------------------------------------
        // topologicalScoreCC
        // ----------------------------------------------------------------------------------

        // Get the topological score cut value
        const auto minTopologicalScoreCC2 = cutTracker.GetCutValue("topologicalScoreCC");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScoreCC2)
            return false;

        // Mark the cut "topologicalScoreCC" as passed
        cutTracker.MarkCutAsPassed("topologicalScoreCC");

        // ##################################################################################
        // Start of pion selection
        // ##################################################################################

        // // ----------------------------------------------------------------------------------
        // // passesCCInclusive
        // // ----------------------------------------------------------------------------------

        // // Insist the event passes the CC inclusive preselection
        // if (!pEvent->reco.passesCCInclusive())
        //     return false;

        // // Mark the cut "passesCCInclusive" as passed
        // cutTracker.MarkCutAsPassed("passesCCInclusive");
        // //std::cout<<"CC0pi: passesCCInclusive"<<std::endl;

        // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 3\n";
        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        for (const auto i:daughterRecoParticleIndices)//unsigned int i = 0; i < pionSectionRecoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");
        //std::cout<<"CC0pi: min2Tracks"<<std::endl;

        // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 4\n";
        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;


        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");
        //std::cout<<"CC0pi: max1Uncontained"<<std::endl;

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT, ccInclusiveMuonCandidateIndex);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 5\n";
        // ----------------------------------------------------------------------------------
        // 1NonProton
        // ----------------------------------------------------------------------------------
        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThresholdLow = cutTracker.GetCutValue("1NonProton");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        unsigned int nBadProtons = 0;
        unsigned int leadingProtonIndex = std::numeric_limits<unsigned int>::max();
        float leadingProtonMom = 0;
        std::vector<unsigned int> protonIndices;
        std::vector<unsigned int> pionIndices;

        for (const auto i:daughterRecoParticleIndices)//unsigned int i = 0; i < pionSectionRecoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
                continue;

            // Assume particles without a track fit are just small protons
            // Don't check for leading proton because without a track fit we can't get momentum
            // Don't "count" them as protons because we need at least one good proton for this selection
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nBadProtons ++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("GetCC0piSelectionModifiedPeLEE - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            // Don't "count" it as a proton or check if it's a leading proton because we need at least one good proton for this selection
            if (!hasFeatures)
            {
                nBadProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThresholdLow)// && bdtResponse <= protonBDTThresholdHigh) //todo: why? protonBDTThresholdHigh
            {
                nProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                float protonmom = AnalysisHelper::GetProtonMomentumFromRange(particle.range());
                if (protonmom > leadingProtonMom){
                    leadingProtonMom = protonmom;
                    leadingProtonIndex = i;
                }

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nBadProtons + nPions + nMuons != daughterRecoParticleIndices.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen!");

        // Insist that we exacly zero pions (i.e we have only 1 non-proton, the muon)
        if (nPions != 0)
            return false;

        // Mark the cut "1NonProton" as passed
        cutTracker.MarkCutAsPassed("1NonProton");
        //std::cout<<"CC0pi: 1NonProton"<<std::endl;

        // Also require that there is at least 1 identified proton
        if (nProtons == 0)
            return false;
        cutTracker.MarkCutAsPassed("atLeast1Proton");
        //std::cout<<"CC0pi: atLeast1Proton"<<std::endl;

        // ----------------------------------------------------------------------------------
        // muonLikeProton
        // ----------------------------------------------------------------------------------
        const auto muonBDTThreshold = cutTracker.GetCutValue("muonLikeProton");

        // Leading proton must have a muon-like muon BDT score
        // Get the leading proton reco particle
        const auto &leadingproton = recoParticles.at(leadingProtonIndex);

        // Get run the muon BDT
        std::vector<float> muonFeatures;
        const auto hasFeaturesMuonBDT = BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::MuonBDTFeatureNames, muonFeatures);
        if (!hasFeaturesMuonBDT)
        {
            return false;
        }

        // Insist that the BDT response is greater than the cut value
        const auto bdtResponseMu = pMuonBDT->GetResponse(muonFeatures);

        if (bdtResponseMu < muonBDTThreshold)
        {
            return false;
        }

        // Mark the cut "muonLikeProton" as passed
        cutTracker.MarkCutAsPassed("muonLikeProton");
        //std::cout<<"CC0pi: muonLikeProton"<<std::endl;


        // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 7\n";
        // ----------------------------------------------------------------------------------
        // barelyResemblingProton
        // ----------------------------------------------------------------------------------
        const auto protonBDTThresholdHigh = cutTracker.GetCutValue("barelyResemblingProton");

        // Leading proton must have a low (yet proton-like) proton BDT score
        // Get the leading proton reco particle
        // const auto &leadingproton = recoParticles.at(leadingProtonIndex);

        // Get run the proton BDT
        std::vector<float> protonFeatures;
        const auto hasFeaturesProtonBDT = BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::ProtonBDTFeatureNames, protonFeatures);
        if (!hasFeaturesProtonBDT)
        {
            return false;
        }

        // Insist that the BDT response is greater than the cut value
        const auto bdtResponseProton = pProtonBDT->GetResponse(protonFeatures);

        if (bdtResponseProton > protonBDTThresholdHigh)
        {
            return false;
        }

        // Mark the cut "barelyResemblingProton" as passed
        cutTracker.MarkCutAsPassed("barelyResemblingProton");

        // // Get the leading proton reco particle
        // const auto &leadingproton = recoParticles.at(leadingProtonIndex);
        // const auto muonBDTThreshold = cutTracker.GetCutValue("muonLikeProton");
        // const auto protonBDTThresholdHigh = cutTracker.GetCutValue("barelyResemblingProton");
        // bool suitableProton = false;
        // for(const auto &protonIndex: protonIndices)
        // {
        //     const auto proton = recoParticles.at(protonIndex);
        //     std::vector<float> muonFeatures, protonFeatures;
        //     const auto hasFeaturesMuonBDT = BDTHelper::GetBDTFeatures(proton, BDTHelper::MuonBDTFeatureNames, muonFeatures);
        //     const auto hasFeaturesProtonBDT = BDTHelper::GetBDTFeatures(proton, BDTHelper::ProtonBDTFeatureNames, protonFeatures);
        //     if (!hasFeaturesMuonBDT || !hasFeaturesProtonBDT)
        //         continue;

        //     if(pMuonBDT->GetResponse(muonFeatures) < muonBDTThreshold && pProtonBDT->GetResponse(protonFeatures) < protonBDTThresholdHigh)
        //     {
        //         suitableProton = true;
        //         break;
        //     }
        // }

        // if (!suitableProton)
        //     return false;

        // cutTracker.MarkCutAsPassed("muonLikeProton");
        // cutTracker.MarkCutAsPassed("barelyResemblingProton");
        // ----------------------------------------------------------------------------------
        // protonHasValiddEdx
        // ----------------------------------------------------------------------------------

        const auto protondEdxThreshold = cutTracker.GetCutValue("protonHasValiddEdx");

        // Insist that the leading proton has a valid deEdx measurement (i.e. truncated mean dE/dx is above threshold)
        if (leadingproton.truncatedMeandEdx() < protondEdxThreshold)
            return false;

        // Mark the cut "protonHasValiddEdx" as passed
        cutTracker.MarkCutAsPassed("protonHasValiddEdx");
        //std::cout<<"CC0pi: protonHasValiddEdx"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 8\n";
        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");
        //std::cout<<"CC0pi: muonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 9\n";
        // ----------------------------------------------------------------------------------
        // protonNotInGap
        // ----------------------------------------------------------------------------------

        // Sanity check that our muon and leading proton are not the same particle
        if (muonIndex == leadingProtonIndex){
            throw std::logic_error("CC0piSelection - The muon and the proton candidates are the same particle. This shouldn't happen!");
        }

        // Sanity check that we have a lehading proton
        if (leadingProtonIndex == std::numeric_limits<unsigned int>::max()){
            throw std::logic_error("CC0piSelection - No leading proton candidate found. This shouldn't happen!");
        }

        // Insist that the leading proton has at least one hit in each view (i.e. not in a gap)
        if (leadingproton.nHitsU() == 0 || leadingproton.nHitsV() == 0 || leadingproton.nHitsW() == 0)
            return false;

        // Mark the cut "protonNotInGap" as passed
        cutTracker.MarkCutAsPassed("protonNotInGap");
        //std::cout<<"CC0pi: protonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 10\n";
        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------
        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");

        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto protonDir = TVector3(leadingproton.directionX(), leadingproton.directionY(), leadingproton.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(protonDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");
        //std::cout<<"CC0pi: openingAngle"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 11\n";
        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------
        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");
        //std::cout<<"CC0pi: topologicalScore"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 12\n";
        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------
        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto i:daughterRecoParticleIndices)
        {
            const auto &particle = recoParticles.at(i);
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
                return false;
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");
        //std::cout<<"CC0pi: startNearVertex"<<std::endl;

        // ----------------------------------------------------------------------------------
        // likelyGoldenProton
        // ----------------------------------------------------------------------------------

        // Get the likely golden pion cut value
        const auto goldenPionBDTThreshold = cutTracker.GetCutValue("likelyGoldenProton");

        // Get the golden pion BDT
        auto &pGoldenPionBDT = bdtMap.at("goldenPion");

        // Get the features of the pion
        std::vector<float> features;
        if (!BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::GoldenPionBDTFeatureNames, features))
            throw std::logic_error("DefaultSelection - Can't get golden pion BDT features of pion candidate");

        // Insist that the BDT response is greater than the cut value to identify the pion as a golden pion
        const auto bdtResponse = pGoldenPionBDT->GetResponse(features);
        if (bdtResponse <= goldenPionBDTThreshold)
            return false;

        // Mark the cut "likelyGoldenProton" as passed
        cutTracker.MarkCutAsPassed("likelyGoldenProton");

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 13\n";
        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCC0piSelection()
{
    // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 0\n";
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"1NonProton", -0.06f},
        {"atLeast1Proton", 0.1f},
        {"muonLikeProton", -0.4f}, //todo: why?
        {"protonHasValiddEdx", 1.0f},
        {"muonNotInGap"},
        {"protonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)}
    };
    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 2\n";
        // ----------------------------------------------------------------------------------
        // passesCCInclusive
        // ----------------------------------------------------------------------------------
        // Insist the event passes the CC inclusive preselection
        if (!pEvent->reco.passesCCInclusive())
            return false;

        // Mark the cut "passesCCInclusive" as passed
        cutTracker.MarkCutAsPassed("passesCCInclusive");
        //std::cout<<"CC0pi: passesCCInclusive"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 3\n";
        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");
        //std::cout<<"CC0pi: min2Tracks"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 4\n";
        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;

        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");
        //std::cout<<"CC0pi: max1Uncontained"<<std::endl;

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 5\n";
        // ----------------------------------------------------------------------------------
        // 1NonProton
        // ----------------------------------------------------------------------------------
        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThresholdLow = cutTracker.GetCutValue("1NonProton");
        const auto protonBDTThresholdHigh = cutTracker.GetCutValue("atLeast1Proton"); //todo: why

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        unsigned int nBadProtons = 0;
        unsigned int leadingProtonIndex = std::numeric_limits<unsigned int>::max();
        float leadingProtonMom = 0;
        std::vector<unsigned int> protonIndices;
        std::vector<unsigned int> pionIndices;

        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
                continue;

            // Assume particles without a track fit are just small protons
            // Don't check for leading proton because without a track fit we can't get momentum
            // Don't "count" them as protons because we need at least one good proton for this selection
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nBadProtons ++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("GetCC0piSelection - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            // Don't "count" it as a proton or check if it's a leading proton because we need at least one good proton for this selection
            if (!hasFeatures)
            {
                nBadProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThresholdLow && bdtResponse <= protonBDTThresholdHigh) //todo: why? protonBDTThresholdHigh
            {
                nProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                float protonmom = AnalysisHelper::GetProtonMomentumFromRange(particle.range());
                if (protonmom > leadingProtonMom){
                    leadingProtonMom = protonmom;
                    leadingProtonIndex = i;
                }

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nBadProtons + nPions + nMuons != recoParticles.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen!");

        // Insist that we exacly zero pions (i.e we have only 1 non-proton, the muon)
        if (nPions != 0)
            return false;

        // Mark the cut "1NonProton" as passed
        cutTracker.MarkCutAsPassed("1NonProton");
        //std::cout<<"CC0pi: 1NonProton"<<std::endl;

        // Also require that there is at least 1 identified proton
        if (nProtons == 0)
            return false;
        cutTracker.MarkCutAsPassed("atLeast1Proton");
        //std::cout<<"CC0pi: atLeast1Proton"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 6\n";
        // ----------------------------------------------------------------------------------
        // muonLikeProton
        // ----------------------------------------------------------------------------------
        const auto muonBDTThreshold = cutTracker.GetCutValue("muonLikeProton");

        // Leading proton must have a muon-like muon BDT score
        // Get the leading proton reco particle
        const auto &leadingproton = recoParticles.at(leadingProtonIndex);

        // Get run the muon BDT
        std::vector<float> features;
        const auto hasFeatures = BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::MuonBDTFeatureNames, features);
        if (!hasFeatures)
        {
            return false;
        }

        // Insist that the BDT response is greater than the cut value
        const auto bdtResponsemu = pMuonBDT->GetResponse(features);

        if (bdtResponsemu < muonBDTThreshold)
        {
            return false;
        }

        // Mark the cut "muonLikeProton" as passed
        cutTracker.MarkCutAsPassed("muonLikeProton");
        //std::cout<<"CC0pi: muonLikeProton"<<std::endl;


        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 7\n";
        // ----------------------------------------------------------------------------------
        // protonHasValiddEdx
        // ----------------------------------------------------------------------------------

        const auto protondEdxThreshold = cutTracker.GetCutValue("protonHasValiddEdx");

        // Insist that the leading proton has a valid deEdx measurement (i.e. truncated mean dE/dx is above threshold)
        if (leadingproton.truncatedMeandEdx() < protondEdxThreshold)
            return false;

        // Mark the cut "protonHasValiddEdx" as passed
        cutTracker.MarkCutAsPassed("protonHasValiddEdx");
        //std::cout<<"CC0pi: protonHasValiddEdx"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 8\n";
        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");
        //std::cout<<"CC0pi: muonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 9\n";
        // ----------------------------------------------------------------------------------
        // protonNotInGap
        // ----------------------------------------------------------------------------------

        // Sanity check that our muon and leading proton are not the same particle
        if (muonIndex == leadingProtonIndex){
            throw std::logic_error("CC0piSelection - The muon and the proton candidates are the same particle. This shouldn't happen!");
        }

        // Sanity check that we have a lehading proton
        if (leadingProtonIndex == std::numeric_limits<unsigned int>::max()){
            throw std::logic_error("CC0piSelection - No leading proton candidate found. This shouldn't happen!");
        }

        // Insist that the leading proton has at least one hit in each view (i.e. not in a gap)
        if (leadingproton.nHitsU() == 0 || leadingproton.nHitsV() == 0 || leadingproton.nHitsW() == 0)
            return false;

        // Mark the cut "protonNotInGap" as passed
        cutTracker.MarkCutAsPassed("protonNotInGap");
        //std::cout<<"CC0pi: protonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 10\n";
        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------
        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");

        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto protonDir = TVector3(leadingproton.directionX(), leadingproton.directionY(), leadingproton.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(protonDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");
        //std::cout<<"CC0pi: openingAngle"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 11\n";
        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------
        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");
        //std::cout<<"CC0pi: topologicalScore"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 12\n";
        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------
        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
                return false;
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");
        //std::cout<<"CC0pi: startNearVertex"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelection - Point 13\n";
        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCC0piSelectionModified(const float muonLikeProtonValue /*= -0.4f*/, const float barelyResemblingProtonValue /*= 0.1f*/)
{
    // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 0\n";
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"1NonProton", -0.06f},
        {"atLeast1Proton"},
        {"muonLikeProton", muonLikeProtonValue}, //todo: why?
        {"barelyResemblingProton", barelyResemblingProtonValue},
        {"protonHasValiddEdx", 1.0f},
        {"muonNotInGap"},
        {"protonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f},
        {"likelyGoldenProton", -0.03f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)},
        {"goldenPion", std::make_shared<BDTHelper::BDT>("goldenPion", BDTHelper::GoldenPionBDTFeatureNames)}
    };
    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 2\n";
        // ----------------------------------------------------------------------------------
        // passesCCInclusive
        // ----------------------------------------------------------------------------------
        // Insist the event passes the CC inclusive preselection
        if (!pEvent->reco.passesCCInclusive())
            return false;

        // Mark the cut "passesCCInclusive" as passed
        cutTracker.MarkCutAsPassed("passesCCInclusive");
        //std::cout<<"CC0pi: passesCCInclusive"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 3\n";
        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");
        //std::cout<<"CC0pi: min2Tracks"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 4\n";
        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;

        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");
        //std::cout<<"CC0pi: max1Uncontained"<<std::endl;

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 5\n";
        // ----------------------------------------------------------------------------------
        // 1NonProton
        // ----------------------------------------------------------------------------------
        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThresholdLow = cutTracker.GetCutValue("1NonProton");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        unsigned int nBadProtons = 0;
        unsigned int leadingProtonIndex = std::numeric_limits<unsigned int>::max();
        float leadingProtonMom = 0;
        std::vector<unsigned int> protonIndices;
        std::vector<unsigned int> pionIndices;

        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
                continue;

            // Assume particles without a track fit are just small protons
            // Don't check for leading proton because without a track fit we can't get momentum
            // Don't "count" them as protons because we need at least one good proton for this selection
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nBadProtons ++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("GetCC0piSelectionModified - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            // Don't "count" it as a proton or check if it's a leading proton because we need at least one good proton for this selection
            if (!hasFeatures)
            {
                nBadProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThresholdLow)// && bdtResponse <= protonBDTThresholdHigh) //todo: why? protonBDTThresholdHigh
            {
                nProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                float protonmom = AnalysisHelper::GetProtonMomentumFromRange(particle.range());
                if (protonmom > leadingProtonMom){
                    leadingProtonMom = protonmom;
                    leadingProtonIndex = i;
                }

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nBadProtons + nPions + nMuons != recoParticles.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen!");

        // Insist that we exacly zero pions (i.e we have only 1 non-proton, the muon)
        if (nPions != 0)
            return false;

        // Mark the cut "1NonProton" as passed
        cutTracker.MarkCutAsPassed("1NonProton");
        //std::cout<<"CC0pi: 1NonProton"<<std::endl;

        // Also require that there is at least 1 identified proton
        if (nProtons == 0)
            return false;
        cutTracker.MarkCutAsPassed("atLeast1Proton");
        //std::cout<<"CC0pi: atLeast1Proton"<<std::endl;

        // ----------------------------------------------------------------------------------
        // muonLikeProton
        // ----------------------------------------------------------------------------------
        const auto muonBDTThreshold = cutTracker.GetCutValue("muonLikeProton");

        // Leading proton must have a muon-like muon BDT score
        // Get the leading proton reco particle
        const auto &leadingproton = recoParticles.at(leadingProtonIndex);

        // Get run the muon BDT
        std::vector<float> muonFeatures;
        const auto hasFeaturesMuonBDT = BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::MuonBDTFeatureNames, muonFeatures);
        if (!hasFeaturesMuonBDT)
        {
            return false;
        }

        // Insist that the BDT response is greater than the cut value
        const auto bdtResponseMu = pMuonBDT->GetResponse(muonFeatures);

        if (bdtResponseMu < muonBDTThreshold)
        {
            return false;
        }

        // Mark the cut "muonLikeProton" as passed
        cutTracker.MarkCutAsPassed("muonLikeProton");
        //std::cout<<"CC0pi: muonLikeProton"<<std::endl;


        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 7\n";
        // ----------------------------------------------------------------------------------
        // barelyResemblingProton
        // ----------------------------------------------------------------------------------
        const auto protonBDTThresholdHigh = cutTracker.GetCutValue("barelyResemblingProton");

        // Leading proton must have a low (yet proton-like) proton BDT score
        // Get the leading proton reco particle
        // const auto &leadingproton = recoParticles.at(leadingProtonIndex);

        // Get run the proton BDT
        std::vector<float> protonFeatures;
        const auto hasFeaturesProtonBDT = BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::ProtonBDTFeatureNames, protonFeatures);
        if (!hasFeaturesProtonBDT)
        {
            return false;
        }

        // Insist that the BDT response is greater than the cut value
        const auto bdtResponseProton = pProtonBDT->GetResponse(protonFeatures);

        if (bdtResponseProton > protonBDTThresholdHigh)
        {
            return false;
        }

        // Mark the cut "barelyResemblingProton" as passed
        cutTracker.MarkCutAsPassed("barelyResemblingProton");

        // // Get the leading proton reco particle
        // const auto &leadingproton = recoParticles.at(leadingProtonIndex);
        // const auto muonBDTThreshold = cutTracker.GetCutValue("muonLikeProton");
        // const auto protonBDTThresholdHigh = cutTracker.GetCutValue("barelyResemblingProton");
        // bool suitableProton = false; 
        // for(const auto &protonIndex: protonIndices)
        // {
        //     const auto proton = recoParticles.at(protonIndex);
        //     std::vector<float> muonFeatures, protonFeatures;
        //     const auto hasFeaturesMuonBDT = BDTHelper::GetBDTFeatures(proton, BDTHelper::MuonBDTFeatureNames, muonFeatures);
        //     const auto hasFeaturesProtonBDT = BDTHelper::GetBDTFeatures(proton, BDTHelper::ProtonBDTFeatureNames, protonFeatures);
        //     if (!hasFeaturesMuonBDT || !hasFeaturesProtonBDT)
        //         continue;

        //     if(pMuonBDT->GetResponse(muonFeatures) < muonBDTThreshold && pProtonBDT->GetResponse(protonFeatures) < protonBDTThresholdHigh)
        //     {
        //         suitableProton = true;
        //         break;
        //     }
        // }

        // if (!suitableProton)
        //     return false;

        // cutTracker.MarkCutAsPassed("muonLikeProton");
        // cutTracker.MarkCutAsPassed("barelyResemblingProton");
        // ----------------------------------------------------------------------------------
        // protonHasValiddEdx
        // ----------------------------------------------------------------------------------

        const auto protondEdxThreshold = cutTracker.GetCutValue("protonHasValiddEdx");

        // Insist that the leading proton has a valid deEdx measurement (i.e. truncated mean dE/dx is above threshold)
        if (leadingproton.truncatedMeandEdx() < protondEdxThreshold)
            return false;

        // Mark the cut "protonHasValiddEdx" as passed
        cutTracker.MarkCutAsPassed("protonHasValiddEdx");
        //std::cout<<"CC0pi: protonHasValiddEdx"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 8\n";
        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");
        //std::cout<<"CC0pi: muonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 9\n";
        // ----------------------------------------------------------------------------------
        // protonNotInGap
        // ----------------------------------------------------------------------------------

        // Sanity check that our muon and leading proton are not the same particle
        if (muonIndex == leadingProtonIndex){
            throw std::logic_error("CC0piSelection - The muon and the proton candidates are the same particle. This shouldn't happen!");
        }

        // Sanity check that we have a lehading proton
        if (leadingProtonIndex == std::numeric_limits<unsigned int>::max()){
            throw std::logic_error("CC0piSelection - No leading proton candidate found. This shouldn't happen!");
        }

        // Insist that the leading proton has at least one hit in each view (i.e. not in a gap)
        if (leadingproton.nHitsU() == 0 || leadingproton.nHitsV() == 0 || leadingproton.nHitsW() == 0)
            return false;

        // Mark the cut "protonNotInGap" as passed
        cutTracker.MarkCutAsPassed("protonNotInGap");
        //std::cout<<"CC0pi: protonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 10\n";
        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------
        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");

        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto protonDir = TVector3(leadingproton.directionX(), leadingproton.directionY(), leadingproton.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(protonDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");
        //std::cout<<"CC0pi: openingAngle"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 11\n";
        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------
        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");
        //std::cout<<"CC0pi: topologicalScore"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 12\n";
        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------
        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
                return false;
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");
        //std::cout<<"CC0pi: startNearVertex"<<std::endl;

        // ----------------------------------------------------------------------------------
        // likelyGoldenProton
        // ----------------------------------------------------------------------------------

        // Get the likely golden pion cut value
        const auto goldenPionBDTThreshold = cutTracker.GetCutValue("likelyGoldenProton");

        // Get the golden pion BDT
        auto &pGoldenPionBDT = bdtMap.at("goldenPion");

        // Get the features of the pion
        std::vector<float> features;
        if (!BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::GoldenPionBDTFeatureNames, features))
            throw std::logic_error("DefaultSelection - Can't get golden pion BDT features of pion candidate");

        // Insist that the BDT response is greater than the cut value to identify the pion as a golden pion
        const auto bdtResponse = pGoldenPionBDT->GetResponse(features);
        if (bdtResponse <= goldenPionBDTThreshold)
            return false;

        // Mark the cut "likelyGoldenProton" as passed
        cutTracker.MarkCutAsPassed("likelyGoldenProton");

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 13\n";
        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCC0piSelectionModifiedPart1()
{
    // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 0\n";
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"passesCCInclusive"},
        {"min2Tracks"},
        {"max1Uncontained"},
        {"1NonProton", -0.06f},
        {"atLeast1Proton"},
        // {"muonLikeProton", muonLikeProtonValue}, //todo: why?
        // {"barelyResemblingProton", barelyResemblingProtonValue},
        {"protonHasValiddEdx", 1.0f},
        {"muonNotInGap"},
        {"protonNotInGap"},
        {"openingAngle", 2.65f},
        {"topologicalScore", 0.67f},
        {"startNearVertex", 9.5f}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)}
    };
    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 2\n";
        // ----------------------------------------------------------------------------------
        // passesCCInclusive
        // ----------------------------------------------------------------------------------
        // Insist the event passes the CC inclusive preselection
        if (!pEvent->reco.passesCCInclusive())
            return false;

        // Mark the cut "passesCCInclusive" as passed
        cutTracker.MarkCutAsPassed("passesCCInclusive");
        //std::cout<<"CC0pi: passesCCInclusive"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 3\n";
        // ----------------------------------------------------------------------------------
        // min2Tracks
        // ----------------------------------------------------------------------------------

        // Count the particles with a track fit and check if they are contained
        unsigned int nTrackParticles = 0u;
        unsigned int nUncontainedParticles = 0u;

        const auto &recoParticles = pEvent->reco.particles;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            nTrackParticles++;

            if (!AnalysisHelper::IsContained(particle))
                nUncontainedParticles++;
        }

        // Insist that we have at least two tracks
        if (nTrackParticles < 2)
            return false;

        // Mark the cut "min2Tracks" as passed
        cutTracker.MarkCutAsPassed("min2Tracks");
        //std::cout<<"CC0pi: min2Tracks"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 4\n";
        // ----------------------------------------------------------------------------------
        // max1Uncontained
        // ----------------------------------------------------------------------------------

        // Insist that at most one particle is uncontained
        if (nUncontainedParticles > 1)
            return false;

        // Mark the cut "max1Uncontained" as passed
        cutTracker.MarkCutAsPassed("max1Uncontained");
        //std::cout<<"CC0pi: max1Uncontained"<<std::endl;

        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 5\n";
        // ----------------------------------------------------------------------------------
        // 1NonProton
        // ----------------------------------------------------------------------------------
        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThresholdLow = cutTracker.GetCutValue("1NonProton");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");

        // Keep track of the number of protons and pions we have identifies
        unsigned int nProtons = 0;
        unsigned int nBadProtons = 0;
        unsigned int leadingProtonIndex = std::numeric_limits<unsigned int>::max();
        float leadingProtonMom = 0;
        std::vector<unsigned int> protonIndices;
        std::vector<unsigned int> pionIndices;

        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
                continue;

            // Assume particles without a track fit are just small protons
            // Don't check for leading proton because without a track fit we can't get momentum
            // Don't "count" them as protons because we need at least one good proton for this selection
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                nBadProtons ++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("GetCC0piSelectionModifiedPart1 - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            // Don't "count" it as a proton or check if it's a leading proton because we need at least one good proton for this selection
            if (!hasFeatures)
            {
                nBadProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThresholdLow)// && bdtResponse <= protonBDTThresholdHigh) //todo: why? protonBDTThresholdHigh
            {
                nProtons++;
                protonIndices.push_back(i);
                cutTracker.AssignPDGCode(i, 2212);

                float protonmom = AnalysisHelper::GetProtonMomentumFromRange(particle.range());
                if (protonmom > leadingProtonMom){
                    leadingProtonMom = protonmom;
                    leadingProtonIndex = i;
                }

                continue;
            }

            // If we've got to this point then we haven't identified the particle as a muon or a proton.
            // Instead let's identify the particle as a pion
            pionIndices.push_back(i);
            cutTracker.AssignPDGCode(i, 211);
        }

        // Sanity check that we have identified every particle
        const auto nMuons = 1u;
        const auto nPions = pionIndices.size();
        if (nProtons + nBadProtons + nPions + nMuons != recoParticles.size())
            throw std::logic_error("DefaultSelection - Identified the wrong number of particles. This shouldn't happen!");

        // Insist that we exacly zero pions (i.e we have only 1 non-proton, the muon)
        if (nPions != 0)
            return false;

        // Mark the cut "1NonProton" as passed
        cutTracker.MarkCutAsPassed("1NonProton");
        //std::cout<<"CC0pi: 1NonProton"<<std::endl;

        // Also require that there is at least 1 identified proton
        if (nProtons == 0)
            return false;
        cutTracker.MarkCutAsPassed("atLeast1Proton");
        //std::cout<<"CC0pi: atLeast1Proton"<<std::endl;


        // ----------------------------------------------------------------------------------
        // protonHasValiddEdx
        // ----------------------------------------------------------------------------------

        const auto protondEdxThreshold = cutTracker.GetCutValue("protonHasValiddEdx");
        const auto &leadingproton = recoParticles.at(leadingProtonIndex);
        // Insist that the leading proton has a valid deEdx measurement (i.e. truncated mean dE/dx is above threshold)
        if (leadingproton.truncatedMeandEdx() < protondEdxThreshold)
            return false;

        // Mark the cut "protonHasValiddEdx" as passed
        cutTracker.MarkCutAsPassed("protonHasValiddEdx");
        //std::cout<<"CC0pi: protonHasValiddEdx"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 8\n";
        // ----------------------------------------------------------------------------------
        // muonNotInGap
        // ----------------------------------------------------------------------------------

        // Get the muon reco particle
        const auto &muon = recoParticles.at(muonIndex);

        // Insist that the muon has at least one hit in each view (i.e. not in a gap)
        if (muon.nHitsU() == 0 || muon.nHitsV() == 0 || muon.nHitsW() == 0)
            return false;

        // Mark the cut "muonNotInGap" as passed
        cutTracker.MarkCutAsPassed("muonNotInGap");
        //std::cout<<"CC0pi: muonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 9\n";
        // ----------------------------------------------------------------------------------
        // protonNotInGap
        // ----------------------------------------------------------------------------------

        // Sanity check that our muon and leading proton are not the same particle
        if (muonIndex == leadingProtonIndex){
            throw std::logic_error("CC0piSelection - The muon and the proton candidates are the same particle. This shouldn't happen!");
        }

        // Sanity check that we have a lehading proton
        if (leadingProtonIndex == std::numeric_limits<unsigned int>::max()){
            throw std::logic_error("CC0piSelection - No leading proton candidate found. This shouldn't happen!");
        }

        // Insist that the leading proton has at least one hit in each view (i.e. not in a gap)
        if (leadingproton.nHitsU() == 0 || leadingproton.nHitsV() == 0 || leadingproton.nHitsW() == 0)
            return false;

        // Mark the cut "protonNotInGap" as passed
        cutTracker.MarkCutAsPassed("protonNotInGap");
        //std::cout<<"CC0pi: protonNotInGap"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 10\n";
        // ----------------------------------------------------------------------------------
        // openingAngle
        // ----------------------------------------------------------------------------------
        // Get the opening angle cut value
        const auto maxOpeningAngle = cutTracker.GetCutValue("openingAngle");

        // Get the opening angle between the muon and pion
        const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
        const auto protonDir = TVector3(leadingproton.directionX(), leadingproton.directionY(), leadingproton.directionZ()).Unit();
        const auto openingAngle = muonDir.Angle(protonDir);

        // Insist that the opening angle isn't too wide
        if (openingAngle >= maxOpeningAngle)
            return false;

        // Mark the cut "openingAngle" as passed
        cutTracker.MarkCutAsPassed("openingAngle");
        //std::cout<<"CC0pi: openingAngle"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 11\n";
        // ----------------------------------------------------------------------------------
        // topologicalScore
        // ----------------------------------------------------------------------------------
        // Get the topological score cut value
        const auto minTopologicalScore = cutTracker.GetCutValue("topologicalScore");

        // Insist that the topological score is above the cut value
        if (pEvent->reco.selectedTopologicalScore() <= minTopologicalScore)
            return false;

        // Mark the cut "topologicalScore" as passed
        cutTracker.MarkCutAsPassed("topologicalScore");
        //std::cout<<"CC0pi: topologicalScore"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 12\n";
        // ----------------------------------------------------------------------------------
        // startNearVertex
        // ----------------------------------------------------------------------------------
        // Get the start near vertex cut value
        const auto maxVertexDist = cutTracker.GetCutValue("startNearVertex");
        const auto maxVertexDist2 = maxVertexDist*maxVertexDist;

        // Insist that all particles with a fitted track start near the vertex
        const auto recoVertex = pEvent->reco.nuVertex();
        for (const auto &particle : recoParticles)
        {
            // Skip particles without a track fit
            if (!AnalysisHelper::HasTrackFit(particle))
                continue;

            // Get the distance between the particle's start position and the vertex
            const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
            const auto vertexDist2 = (start - recoVertex).Mag2();

            // Insist that this isn't too large
            if (vertexDist2 > maxVertexDist2)
                return false;
        }

        // Mark the cut "startNearVertex" as passed
        cutTracker.MarkCutAsPassed("startNearVertex");
        //std::cout<<"CC0pi: startNearVertex"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 13\n";
        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

SelectionHelper::EventSelection SelectionHelper::GetCC0piSelectionModifiedPart2(const float muonLikeProtonValue /*= -0.4f*/, const float barelyResemblingProtonValue /*= 0.1f*/)
{
    // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point -1"<<std::endl;
    // // std::cout<<"SelectionHelper::GetCC0piSelectionModified - Point 0\n";
    // Define the cuts
    const std::vector<EventSelection::Cut> cuts = {
        {"1NonProton", -0.06f},// only for cut value
        {"muonLikeProton", muonLikeProtonValue},
        {"barelyResemblingProton", barelyResemblingProtonValue}
    };

    // Load up the BDTs and store them in a map
    EventSelection::BDTMap bdtMap = {
        {"muon",       std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames)},
        {"proton",     std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames)}
    };
    // Define the actual selection logic
    const auto logic = [](const auto &pEvent, auto &bdtMap, auto &cutTracker)
    {
        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 0"<<std::endl;
        const auto &recoParticles = pEvent->reco.particles;
        // Identify the muon candidate
        auto &pMuonBDT = bdtMap.at("muon");
        const auto muonIndex = SelectionHelper::GetMuonCandidateIndex(recoParticles, BDTHelper::MuonBDTFeatureNames, *pMuonBDT);

        // Assign the muon candidate a muon PDG code
        cutTracker.AssignPDGCode(muonIndex, 13);

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 1"<<std::endl;
        // ----------------------------------------------------------------------------------
        // 1NonProton
        // ----------------------------------------------------------------------------------
        // Identify the rest of the particles using the proton BDT
        const auto protonBDTThresholdLow = cutTracker.GetCutValue("1NonProton");

        // Get the proton BDT from the map
        auto &pProtonBDT = bdtMap.at("proton");
        unsigned int leadingProtonIndex = std::numeric_limits<unsigned int>::max();
        float leadingProtonMom = 0;
        std::vector<unsigned int> protonIndices;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 2"<<std::endl;
        for (unsigned int i = 0; i < recoParticles.size(); ++i)
        {
            const auto &particle = recoParticles.at(i);

            // Skip the muon candidate as we've already identified it
            if (i == muonIndex)
                continue;

            // Assume particles without a track fit are just small protons
            // Don't check for leading proton because without a track fit we can't get momentum
            // Don't "count" them as protons because we need at least one good proton for this selection
            if (!AnalysisHelper::HasTrackFit(particle))
            {
                protonIndices.push_back(i);
                continue;
            }

            // The particle should be contained (as only the muon candidate is allowed to escape). But for sanity, let's check
            if (!AnalysisHelper::IsContained(particle))
                throw std::logic_error("GetCC0piSelectionModifiedPart2 - Found an uncontained particle that isn't the muon. This shouldn't happen!");

            // Get run the proton BDT
            std::vector<float> features;
            const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, BDTHelper::ProtonBDTFeatureNames, features);

            // If one or more of the BDT features are missing, then assume the particle is a proton
            // Don't "count" it as a proton or check if it's a leading proton because we need at least one good proton for this selection
            if (!hasFeatures)
            {
                continue;
            }

            // Insist that the BDT response is greater than the cut value to identify the particle as a proton
            const auto bdtResponse = pProtonBDT->GetResponse(features);
            if (bdtResponse >= protonBDTThresholdLow)// && bdtResponse <= protonBDTThresholdHigh) //todo: why? protonBDTThresholdHigh
            {
                protonIndices.push_back(i);
                float protonmom = AnalysisHelper::GetProtonMomentumFromRange(particle.range());
                if (protonmom > leadingProtonMom){
                    leadingProtonMom = protonmom;
                    leadingProtonIndex = i;
                }
                continue;
            }
        }
        cutTracker.MarkCutAsPassed("1NonProton");

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 3"<<std::endl;
        // ----------------------------------------------------------------------------------
        // muonLikeProton
        // ----------------------------------------------------------------------------------
        const auto muonBDTThreshold = cutTracker.GetCutValue("muonLikeProton");
        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 3.01"<<std::endl;
        // Leading proton must have a muon-like muon BDT score
        // Get the leading proton reco particle
        const auto &leadingproton = recoParticles.at(leadingProtonIndex);

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 3.1"<<std::endl;
        // Get run the muon BDT
        std::vector<float> muonFeatures;
        const auto hasFeaturesMuonBDT = BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::MuonBDTFeatureNames, muonFeatures);
        if (!hasFeaturesMuonBDT)
        {
            return false;
        }

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 3.2"<<std::endl;
        // Insist that the BDT response is greater than the cut value
        const auto bdtResponseMu = pMuonBDT->GetResponse(muonFeatures);

        if (bdtResponseMu < muonBDTThreshold)
        {
            return false;
        }

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 3.3"<<std::endl;
        // Mark the cut "muonLikeProton" as passed
        cutTracker.MarkCutAsPassed("muonLikeProton");
        //std::cout<<"CC0pi: muonLikeProton"<<std::endl;

        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 4"<<std::endl;
        // ----------------------------------------------------------------------------------
        // barelyResemblingProton
        // ----------------------------------------------------------------------------------
        const auto protonBDTThresholdHigh = cutTracker.GetCutValue("barelyResemblingProton");

        // Leading proton must have a low (yet proton-like) proton BDT score
        // Get the leading proton reco particle
        // const auto &leadingproton = recoParticles.at(leadingProtonIndex);

        // Get run the proton BDT
        std::vector<float> protonFeatures;
        const auto hasFeaturesProtonBDT = BDTHelper::GetBDTFeatures(leadingproton, BDTHelper::ProtonBDTFeatureNames, protonFeatures);
        if (!hasFeaturesProtonBDT)
        {
            return false;
        }

        // Insist that the BDT response is greater than the cut value
        const auto bdtResponseProton = pProtonBDT->GetResponse(protonFeatures);

        if (bdtResponseProton > protonBDTThresholdHigh)
        {
            return false;
        }

        // Mark the cut "barelyResemblingProton" as passed
        cutTracker.MarkCutAsPassed("barelyResemblingProton");
        // // std::cout<<"SelectionHelper::GetCC0piSelectionModifiedPart2 - Point 5"<<std::endl;
        // We passed all cuts!
        return true;
    };

    return EventSelection(cuts, bdtMap, logic);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SelectionHelper::IsCutPassed(const std::vector<string> &cutsPassed, const std::string &cut)
{
    return std::find(cutsPassed.begin(), cutsPassed.end(), cut) != cutsPassed.end();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int SelectionHelper::GetMuonCandidateIndex(const std::vector<Event::Reco::Particle> &particles, const std::vector<std::string> &featureNames, BDTHelper::BDT &muonBDT, const int ccInclusiveMuonCandidateIndex /*= -1*/)
{
    bool foundCCInclusiveMuon = false;

    unsigned int ccInclusiveMuonIndex = std::numeric_limits<unsigned int>::max();
    unsigned int muonIndex = std::numeric_limits<unsigned int>::max();

    unsigned int nUncontainedParticles = 0;

    if(ccInclusiveMuonCandidateIndex >= 0)
    {
        ccInclusiveMuonIndex = ccInclusiveMuonCandidateIndex;
        foundCCInclusiveMuon = true;
    }

    for (unsigned int index = 0; index < particles.size(); ++index)
    {
        const auto &particle = particles.at(index);
        if(particle.generation() != 2) continue; // This is needed because the PeLEE ntuples comtain all particles

        // When no cc-inclusive muon candidate is explicitly provided (ccInclusiveMuonCandidateIndex), get it via the cc-inclusive muon candidate flag provided by the numucc filter
        if (ccInclusiveMuonCandidateIndex < 0 && particle.isCCInclusiveMuonCandidate())
        {
            if (foundCCInclusiveMuon)
                throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found multiple CC inclusive muon candidates");

            foundCCInclusiveMuon = true;
            ccInclusiveMuonIndex = index;
        }

        if (!AnalysisHelper::HasTrackFit(particle))
            continue;

        // For now make the last escaping particle the muon candidate, if there are multiple we will fall back on the CC inclusive candidate
        if (!AnalysisHelper::IsContained(particle))
        {
            muonIndex = index;
            nUncontainedParticles++;
        }
    }

    if (!foundCCInclusiveMuon)
        throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found no CC inclusive muon candidate");

    if (nUncontainedParticles > 1)
        return ccInclusiveMuonIndex;

    if (nUncontainedParticles == 1)
        return muonIndex;

    // If we are here all particles are contained, choose the muon using the BDT
    float maxMuonBDTResponse = -std::numeric_limits<float>::max();
    bool foundMuon = false;

    for (unsigned int index = 0; index < particles.size(); ++index)
    {
        const auto &particle = particles.at(index);
        if(particle.generation() != 2) continue; // This is needed because the PeLEE ntuples comtain all particles

        if (!AnalysisHelper::HasTrackFit(particle))
            continue;

        if (!AnalysisHelper::IsContained(particle))
            throw std::logic_error("SelectionHelper::GetMuonCandidateIndex - found escaping particle when not expecting to!");

        std::vector<float> features;
        const auto hasFeatures = BDTHelper::GetBDTFeatures(particle, featureNames, features);

        if (!hasFeatures)
            continue;

        const auto muonBDTResponse = muonBDT.GetResponse(features);

        if (muonBDTResponse < maxMuonBDTResponse)
            continue;

        maxMuonBDTResponse = muonBDTResponse;
        muonIndex = index;
        foundMuon = true;
    }

    // If no muon can be found, then default to the CC inclusive candidate
    if (!foundMuon)
        return ccInclusiveMuonIndex;

    return muonIndex;
}
// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int SelectionHelper::GetLeadingProtonCandidateIndex(const std::vector<Event::Reco::Particle> &particles, std::vector<int> const &assignedPdgCodes)
{
  // Default value in case no proton can be found
  unsigned int protonIndex = std::numeric_limits<unsigned int>::max();

  float leadingProtonMom = 0;

  for (unsigned int i = 0; i < assignedPdgCodes.size(); ++i)
  {
      if (assignedPdgCodes.at(i) != 2212)
          continue;

      // Now check if this is the leading proton (i.e. the highest-reconstructed-momentum proton in the event)
      const auto proton = particles.at(i);
      if (!AnalysisHelper::HasTrackFit(proton))
          continue;
      float protonmom = AnalysisHelper::GetProtonMomentumFromRange(proton.range());
      if (protonmom > leadingProtonMom){
        leadingProtonMom = protonmom;
        protonIndex = i;
      }
  }

  return protonIndex;
}

} // namespace ubcc1pi
