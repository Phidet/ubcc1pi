/**
 *  @file  ubcc1pi_standalone/Helpers/SelectionHelper.h
 *
 *  @brief The header file for the selection helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_SELECTION_HELPER
#define UBCC1PI_STANDALONE_HELPERS_SELECTION_HELPER

#include "ubcc1pi_standalone/Interface/Event.h"

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"

#include <memory>
#include <functional>

namespace ubcc1pi
{

/**
 *  @brief  The selection helper class
 */
class SelectionHelper
{
    public:
        /**
         *  @brief  The event selection class
         */
        class EventSelection
        {
            public:

                /**
                *  @brief  The event selection cut class
                */
                class Cut
                {
                    public:
                        /**
                        *  @brief  Default constructor
                        */
                        Cut();

                        /**
                        *  @brief  Construct a cut with a name
                        *
                        *  @param  name the cut name
                        */
                        Cut(const std::string &name);

                        /**
                        *  @brief  Construct a cut with a name and a value
                        *
                        *  @param  name the cut name
                        *  @param  value the cut value
                        */
                        Cut(const std::string &name, const float value);

                        /**
                        *  @brief  Get the cut name
                        *
                        *  @return the name of the cut
                        */
                        std::string GetName() const;

                        /**
                        *  @brief  Get if this cut has a value
                        *
                        *  @return boolean, true if the cut has a value
                        */
                        bool HasValue() const;

                        /**
                        *  @brief  Get the cut value (if available)
                        *
                        *  @return the cut value
                        */
                        float GetValue() const;

                        /**
                        *  @brief  Set the cut value
                        *
                        *  @param  value the cut value
                        */
                        void SetValue(const float value);

                    private:

                        std::string m_name;     ///< The name of the cut
                        bool        m_hasValue; ///< If the cut has a value
                        float       m_value;    ///< The cut value (if it exists)
                };

                /**
                *  @brief  The cut tracker class
                */
                class CutTracker
                {
                    public:
                        /**
                        *  @brief  Constructor
                        *
                        *  @param  cuts the cuts that we are tracking
                        *  @param  nRecoParticles the number of reconstructed particles in the event
                        */
                        CutTracker(const std::vector<Cut> &cuts, const unsigned int nRecoParticles);

                        /**
                        *  @brief  Get the value of one of the cuts by name
                        *
                        *  @param  name the name of the cut
                        *
                        *  @return the value of the cut
                        */
                        float GetCutValue(const std::string &name) const;

                        /**
                        *  @brief  Get the cuts that have been passed
                        *
                        *  @return the cuts passed
                        */
                        std::vector<std::string> GetCutsPassed() const;

                        /**
                        *  @brief  Get the PDG codes that have been assigned to the reco particles
                        *
                        *  @return the assigned PDG codes
                        */
                        std::vector<int> GetAssignedPDGCodes() const;

                        /**
                        *  @brief  Mark the cut with the supplied name as "passed"
                        *          Note that the cuts must be marked in the order they are declared
                        *
                        *  @param  name the name of the cut to mark as passed
                        */
                        void MarkCutAsPassed(const std::string &name);

                        /**
                        *  @brief  Assign a PDG code to a reco particle
                        *
                        *  @param  recoParticleIndex the index of the reco particle
                        *  @param  pdg the PDG code to assign
                        */
                        void AssignPDGCode(const unsigned int recoParticleIndex, const int pdg);

                    private:

                        std::vector<Cut>         m_cuts;         ///< The cuts
                        std::vector<std::string> m_cutsPassed;   ///< The cuts that have been marked as passed
                        std::vector<int>         m_assignedPdgs; ///< The PDG codes that have been assigned to the reco particles
                };

                /**
                *  @brief  The selection result type definition. A tuple of three entries:
                *            0. bool = if the event passed the last cut of the selection
                *            1. vector<string> = the names of the cuts that the selection passed (even if it didn't pass the last cut)
                *            2. vector<int> = one entry per reco particle, the PDG codes that were assigned by the selection
                */
                typedef std::tuple< bool, std::vector<std::string>, std::vector<int> > SelectionResult;

                /**
                *  @brief  A mapping from the name of a BDT to the BDT itself
                */
                typedef std::unordered_map<std::string, std::shared_ptr<BDTHelper::BDT> > BDTMap;

                /**
                *  @brief  The selection logic type definition. A function that return true/false if a given event passed some selection
                *          logic. The function takes as input:
                *            pEvent     an event object to be considered
                *            bdtMap     the BDTs that can be utilized by the selection
                *            cutTracker an object that the selection can use to mark which cuts are passed
                */
                typedef std::function<bool (const std::shared_ptr<Event> &, BDTMap &, CutTracker &)> SelectionLogic;

                /**
                *  @brief  Constructor
                *
                *  @param  cuts the ordered vector of cuts that will be applied
                *  @param  bdtMap the BDTs that the selection will use (can be empty)
                *  @param  logic the event selection logic
                */
                EventSelection(const std::vector<Cut> &cuts, BDTMap &bdtMap, const SelectionLogic &logic);

                /**
                *  @brief  Get the names of the event selection cuts
                *
                *  @return the cut names
                */
                std::vector<std::string> GetCuts() const;

                /**
                *  @brief  Get if a cut has a value
                *
                *  @param  name the name of the cut
                *
                *  @return boolean, true if the cut has a value
                */
                bool CutHasValue(const std::string &name) const;

                /**
                 *  @brief  Get the value of one of the cuts by name
                 *
                 *  @param  name the name of the cut
                 *
                 *  @return the value of the cut
                 */
                float GetCutValue(const std::string &name) const;

                /**
                 *  @brief  Set the value of one of the cuts by name
                 *
                 *  @param  name the name of the cut
                 *  @param  value the new value of the cut
                 */
                void SetCutValue(const std::string &name, const float value);

                /**
                *  @brief  Execute the selection on a given event
                *
                *  @param  pEvent the input event
                *
                *  @return the result of the selection
                */
                SelectionResult Execute(const std::shared_ptr<Event> &pEvent);

            private:

                std::vector<Cut>            m_cuts;    ///< The selection cuts
                BDTMap                      m_bdtMap;  ///< The BDTs
                SelectionLogic              m_logic;   ///< The event selection logic
        };

        /**
         *  @brief  Get a selection based on an input string. String can be "CCInclusive" for the CCInclusive selection, "Default" for the default selection, or "CC0pi" for the CC0pi selection
         *
         *  @return the event selection
         */
        static EventSelection GetSelection(const std::string &name);

        /**
         *  @brief  Get the CC inclusive selection
         *
         *  @return the event selection
         */
        static EventSelection GetCCInclusiveSelection();

        /**
         *  @brief  Get the 'original' CC1Pi event selection (relies on NuMuCC filter)
         *
         *  @return the event selection
         */
        static EventSelection GetOriginalSelection();

        /**
         *  @brief  Get the default CC1Pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetDefaultSelection(const bool skipOldCuts = false);

        /**
         *  @brief  Get the default CC1Pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetDefaultSelection2(const bool skipOldCuts = false);


        /**
         *  @brief  Get the cc inclusive numu event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCCInclusiveSelection2(const bool skipOldCuts = false);

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelection();

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModified(const float muonLikeProtonValue = -0.4f, const float barelyResemblingProtonValue = 0.1f);

        /**
         *  @brief  Get the original CC0pi event selection for PeLEE
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionOriginalPeLEE(const float muonLikeProtonValue = -0.4f, const float barelyResemblingProtonValue = 0.1f, const bool skipOldCuts = false);

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModifiedPeLEE(const float muonLikeProtonValue = -0.4f, const float barelyResemblingProtonValue = 0.1f, const bool skipOldCuts = false);

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModifiedPeLEENoCCInclusive(const float muonLikeProtonValue = -0.4f, const float barelyResemblingProtonValue = 0.1f, const bool skipOldCuts = false);

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModifiedAgainPeLEE(const float muonLikeProtonValue = -0.4f, const float barelyResemblingProtonValue = 0.1f, const bool skipOldCuts = false);

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModifiedPart1();

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModifiedPart2(const float muonLikeProtonValue = -0.4f, const float barelyResemblingProtonValue = 0.1f);

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModifiedPeLEEPart1(const bool skipOldCuts = false);

        /**
         *  @brief  Get the CC0pi event selection
         *
         *  @return the event selection
         */
        static EventSelection GetCC0piSelectionModifiedPeLEEPart2(const int ccInclusiveMuonCandidateIndex, const float muonLikeProtonValue = -0.4f, const float barelyResemblingProtonValue = 0.1f);


        /**
        *  @brief  Check if a given cut is listed in the input vector of cuts passed
        *
        *  @param  cutsPassed the cuts that have been passed
        *  @param  cut the cut to check
        *
        *  @return if the cut is in the cutsPassed vector
        */
        static bool IsCutPassed(const std::vector<string> &cutsPassed, const std::string &cut);

        /**
         *  @brief  Get the muon candidate index
         *
         *  @param  particles the input list of all reco particles
         *  @param  featureNames the input list of muon BDT feature names
         *  @param  muonBDT the muon BDT
         *  @param  ccInclusiveMuonCandidateIndex the override for the index of the muon candidate in the CC inclusive selection
         *
         *  @return the index of the muon candidate in the input list
         */
        static unsigned int GetMuonCandidateIndex(const std::vector<Event::Reco::Particle> &particles, const std::vector<std::string> &featureNames, BDTHelper::BDT &muonBDT, const int ccInclusiveMuonCandidateIndex = -1);

        /**
         *  @brief  Get the leading proton candidate index (works for the CC0pi1p sideband selection only)
         *
         *  @param  particles the input list of all reco particles
         *  @param  assignedPdgCodes the input list of proton BDT feature names
         *
         *  @return the index of the muon candidate in the input list
         */
        static unsigned int GetLeadingProtonCandidateIndex(const std::vector<Event::Reco::Particle> &particles, std::vector<int> const &assignedPdgCodes);

        /**
         *  @brief  Get the leading proton candidate index (works for the CC0pi1p sideband selection only)
         *
         *  @param  particles the input list of all reco particles
         *  @param  assignedPdgCodes the input list of proton BDT feature names
         *  @param  protonBDTThresholdHigh the upper proton BDT threshold
         *  @param  protonBDTThresholdLow the lower proton BDT threshold
         *
         *  @return the index of the muon candidate in the input list
         */
        static unsigned int GetLeadingProtonCandidateIndexByBDTScore(const std::vector<Event::Reco::Particle> &particles, std::vector<int> const &assignedPdgCodes, const float protonBDTThresholdHigh = 0.12f, const float protonBDTThresholdLow = -0.06, std::shared_ptr<BDTHelper::BDT> pProtonBDT = nullptr);

};

} // namespace ubcc1pi

#endif
