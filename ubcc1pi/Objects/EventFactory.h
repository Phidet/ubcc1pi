/**
 *  @file  ubcc1pi/Objects/EventFactory.h
 *
 *  @brief The header file for the event factory class
 */

#ifndef UBCC1PI_OBJECTS_EVENT_FACTORY
#define UBCC1PI_OBJECTS_EVENT_FACTORY

#include "ubcc1pi_standalone/Interface/Event.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"

#include "ubcc1pi/Helpers/CollectionHelper.h"
#include "ubcc1pi/Helpers/BacktrackHelper.h"
#include "ubcc1pi/Helpers/TruthHelper.h"
#include "ubcc1pi/Helpers/RecoHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The event factory class is used to populate a ubcc1pi::Event from an art::Event
 *          This implementation is separated from the Event class itself so that one can use the Event without any coupling to LArSoft
 */
class EventFactory
{
    public:

        /**
         *  @brief  The configuration required to build a ubcc1pi::Event
         */
        struct Config
        {
            /**
             *  @brief  If the input file has truth info
             */
            fhicl::Atom<bool> HasTruthInfo
            {
                fhicl::Name("HasTruthInfo"),
                fhicl::Comment("If the input events have truth info that we should look for")
            };

            /**
             *  @brief  If we should try to collect event weights
             */
            fhicl::Atom<bool> GetEventWeights
            {
                fhicl::Name("GetEventWeights"),
                fhicl::Comment("If we should try to collect event weights")
            };

            /**
             *  @brief  The MCEventWeight label
             */
            fhicl::Atom<art::InputTag> MCEventWeightLabel
            {
                fhicl::Name("MCEventWeightLabel"),
                fhicl::Comment("The label for the MC event weight producer")
            };

            /**
             *  @brief  The MCTruth label
             */
            fhicl::Atom<art::InputTag> MCTruthLabel
            {
                fhicl::Name("MCTruthLabel"),
                fhicl::Comment("The label for the neutrino MCTruth producer")
            };

            /**
             *  @brief  The MCParticle label
             */
            fhicl::Atom<art::InputTag> MCParticleLabel
            {
                fhicl::Name("MCParticleLabel"),
                fhicl::Comment("The label for the MCParticle producer")
            };

            /**
             *  @brief  The backtracker label
             */
            fhicl::Atom<art::InputTag> BacktrackerLabel
            {
                fhicl::Name("BacktrackerLabel"),
                fhicl::Comment("The label for the MCParticle to hit backtracker producer")
            };

            /**
             *  @brief  The hit label
             */
            fhicl::Atom<art::InputTag> HitLabel
            {
                fhicl::Name("HitLabel"),
                fhicl::Comment("The label for the Hit producer")
            };

            /**
             *  @brief  The slice label
             */
            fhicl::Atom<art::InputTag> SliceLabel
            {
                fhicl::Name("SliceLabel"),
                fhicl::Comment("The label for the Slice producer")
            };

            /**
             *  @brief  The PFParticle label
             */
            fhicl::Atom<art::InputTag> PFParticleLabel
            {
                fhicl::Name("PFParticleLabel"),
                fhicl::Comment("The label for the PFParticle producer")
            };

            /**
             *  @brief  The vertex label
             */
            fhicl::Atom<art::InputTag> VertexLabel
            {
                fhicl::Name("VertexLabel"),
                fhicl::Comment("The label for the Vertex producer")
            };

            /**
             *  @brief  The spacepoint label
             */
            fhicl::Atom<art::InputTag> SpacePointLabel
            {
                fhicl::Name("SpacePointLabel"),
                fhicl::Comment("The label for the SpacePoint producer")
            };

            /**
             *  @brief  The track label
             */
            fhicl::Atom<art::InputTag> TrackLabel
            {
                fhicl::Name("TrackLabel"),
                fhicl::Comment("The label for the Track producer")
            };

            /**
             *  @brief  The MCS fit result label
             */
            fhicl::Atom<art::InputTag> MCSFitResultLabel
            {
                fhicl::Name("MCSFitResultLabel"),
                fhicl::Comment("The label for the MCSFitResult producer")
            };

            /**
             *  @brief  The calorimetry label
             */
            fhicl::Atom<art::InputTag> CalorimetryLabel
            {
                fhicl::Name("CalorimetryLabel"),
                fhicl::Comment("The label for the Calorimetry producer")
            };

            /**
             *  @brief  The PID label
             */
            fhicl::Atom<art::InputTag> PIDLabel
            {
                fhicl::Name("PIDLabel"),
                fhicl::Comment("The label for the PID producer")
            };

            /**
             *  @brief  The CC inclusive label
             */
            fhicl::Atom<art::InputTag> CCInclusiveLabel
            {
                fhicl::Name("CCInclusiveLabel"),
                fhicl::Comment("The label for the CC inclusive T0 tag producer")
            };

            /**
             *  @brief  The track end space point distance (see fhicl comment)
             */
            fhicl::Atom<float> TrackEndSpacePointDistance
            {
                fhicl::Name("TrackEndSpacePointDistance"),
                fhicl::Comment("The distance within which a spacepoint must be from the end of a track to be counted in the nSpacePointsNearEnd feature")
            };

            /**
             *  @brief  The squared sin angle threshold (see fhicl comment)
             */
            fhicl::Atom<float> Sin2AngleThreshold
            {
                fhicl::Name("Sin2AngleThreshold"),
                fhicl::Comment("The squared sin of the angular threshold outside of which a track must point with respect to the collection plane wires to use information that plane")
            };

            /**
             *  @brief  The number of hits to skip (see fhicl comment)
             */
            fhicl::Atom<float> NHitsToSkip
            {
                fhicl::Name("NHitsToSkip"),
                fhicl::Comment("The number of hits to skip at the start of a track when calculating the truncated mean dEdx")
            };

            /**
             *  @brief  The length fraction (see fhicl comment)
             */
            fhicl::Atom<float> LengthFraction
            {
                fhicl::Name("LengthFraction"),
                fhicl::Comment("The fration of the track range to consider at the start of the track when calculating the truncated mean dEdx")
            };

            /**
             *  @brief  The PFP-flash matching label
             */
            fhicl::Atom<art::InputTag> FlashMatchLabel
            {
                fhicl::Name("FlashMatchLabel"),
                fhicl::Comment("The label for the Flash Matching producer")
            };

            /**
             *  @brief  The flash label (OpFlash data product)
             */
            fhicl::Atom<art::InputTag> FlashLabel
            {
                fhicl::Name("FlashLabel"),
                fhicl::Comment("The label for the Flash producer (OpFlash)")
            };
        };

        /**
         *  @brief  Populate the output event with information from the input event
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  pOutputEvent the output event to populate
         */
        static void PopulateEvent(const art::Event &event, const Config &config, Event *pOutputEvent);

    private:

        /**
         *  @brief  Populate the metadata of the event
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  metadata the output metadata
         */
        static void PopulateEventMetadata(const art::Event &event, const Config &config, Event::Metadata &metadata);

        /**
         *  @brief  Populate the truth information of the event
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  truth the output truth information
         *  @param  finalStateMCParticles the MCParticles that have been populated, outputted for use in reco-true matching
         */
        static void PopulateEventTruthInfo(const art::Event &event, const Config &config, Event::Truth &truth, MCParticleVector &finalStateMCParticles);

        /**
         *  @brief  Populate the truth information of the slices
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  truth the output truth information
         */
        static void PopulateEventTruthSliceInfo(const art::Event &event, const Config &config, Event::Truth &truth);

        /**
         *  @brief  Populate the truth particle information
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  mcParticle the input MCParticle
         *  @param  mcParticleToHitWeights the input mapping from MCParticles to hits along with the weight
         *  @param  mcParticleMap the hierarchy map for MCParticles
         *  @param  particle the output particle
         */
        static void PopulateEventTruthParticleInfo(const art::Event &event, const Config &config, const art::Ptr<simb::MCParticle> &mcParticle, const MCParticlesToHitWeights &mcParticleToHitWeights, const MCParticleMap &mcParticleMap, Event::Truth::Particle &particle);

        /**
         *  @brief  Populate the reco information of the event
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  finalStateMCParticles the final state MCParticles, used for truth matching
         *  @param  reco the output reco information
         */
        static void PopulateEventRecoInfo(const art::Event &event, const Config &config, const MCParticleVector &finalStateMCParticles, Event::Reco &reco);

        /**
         *  @brief  Populate the reco slice information
         *
         *  @param  event the input event
         *  @param  config the configuration options
         *  @param  reco the output reco information
         */
        static void PopulateEventRecoSliceInfo(const art::Event &event, const Config &config, Event::Reco &reco);

        /**
         *  @brief  Populate the reco particle information
         *
         *  @param  config the congiuration options
         *  @param  pfParticle the input PFParticle
         *  @param  pfParticleMap the PFParticle map
         *  @param  pfParticleToT0s the input mapping from PFParticles to T0s
         *  @param  pfParticleToMetadata the input mapping from PFParticle to metadata
         *  @param  pfParticleToHits the input mapping from PFParticles to hits
         *  @param  pfParticleToTracks the input mapping from PFParticles to tracks
         *  @param  trackToMCSFitResults the input mapping from tracks to MCS fit results
         *  @param  trackToPIDs the input mapping from tracks to PIDs
         *  @param  trackToCalorimetries the input mapping from tracks to Calorimetry objects
         *  @param  spacePoints the input spacepoints
         *  @param  finalStateMCParticles the final state MCParticles, used for truth matching
         *  @param  pBacktrackerData the backtracker data (can pass as null if truth information isn't available)
         *  @param  nuVertex the (not space-charge corrected) reconstructed neutrino vertex
         *  @param  particle the output particle
         */
        static void PopulateEventRecoParticleInfo(const Config &config, const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleMap &pfParticleMap, const PFParticleToT0s &pfParticleToT0s, const PFParticleToMetadata &pfParticleToMetadata, const PFParticleToHits &pfParticleToHits, const PFParticleToTracks &pfParticleToTracks, const TrackToMCSFitResults &trackToMCSFitResults, const TrackToPIDs &trackToPIDs, const TrackToCalorimetries &trackToCalorimetries, const SpacePointVector &spacePoints, const MCParticleVector &finalStateMCParticles, const std::shared_ptr<BacktrackHelper::BacktrackerData> &pBacktrackerData, const TVector3 &nuVertex, Event::Reco::Particle &particle);

        /**
         *  @brief  Check if the input particle is the CC inclusive muon candidate
         *
         *  @param  pfParticle the input PFParticle
         *  @param  pfParticleToT0s the input assocation from PFParticles to T0s
         *
         *  @return boolean, true if the particle is the muon candidate
         */
        static bool IsCCInclusiveMuonCandidate(const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleToT0s &pfParticleToT0s);

        /**
         *  @brief  Populate the reco particle pattern recognition information
         *
         *  @param  pfParticle the input PFParticle
         *  @param  pfParticleMap the PFParticle map
         *  @param  pfParticleToMetadata the input mapping from PFParticle to metadata
         *  @param  pfParticleToHits the input mapping from PFParticle to hits
         *  @param  particle the output particle
         */
        static void PopulateEventRecoParticlePatRecInfo(const art::Ptr<recob::PFParticle> &pfParticle, const PFParticleMap &pfParticleMap, const PFParticleToMetadata &pfParticleToMetadata, const PFParticleToHits &pfParticleToHits, Event::Reco::Particle &particle);

        /**
         *  @brief  Populate the reco particle track information
         *
         *  @param  config the configuration options
         *  @param  track the input track
         *  @param  spacePoints all spacepoins in the event
         *  @param  nuVertex the (not space-charge corrected) reconstructed neutrino vertex
         *  @param  particle the output particle
         */
        static void PopulateEventRecoParticleTrackInfo(const Config &config, const art::Ptr<recob::Track> &track, const SpacePointVector &spacePoints, const TVector3 &nuVertex, Event::Reco::Particle &particle);

        /**
         *  @brief  Populare the MCS fit result information
         *
         *  @param  mcsFitResult the input MCS fit result object
         *  @param  particle the output particle
         */
        static void PopulateEventRecoParticleMCSFitInfo(const art::Ptr<recob::MCSFitResult> &mcsFitResult, Event::Reco::Particle &particle);

        /**
         *  @brief  Populate the reco particle PID information
         *
         *  @param  config the configuration options
         *  @param  pid the input PID object
         *  @param  particle the output particle
         */
        static void PopulateEventRecoParticlePIDInfo(const Config &config, const art::Ptr<anab::ParticleID> &pid, Event::Reco::Particle &particle);

        /**
         *  @brief  Populate the reco particle calorimetry information
         *
         *  @param  config the configuration options
         *  @param  calos the input calorimetry objects
         *  @param  particle the output particle
         */
        static void PopulateEventRecoParticleCalorimetryInfo(const Config &config, const CalorimetryVector &calos, Event::Reco::Particle &particle);

        /**
         *  @brief  Set the value of the input member variable to the bragg likelihood from the PID with given parameters
         *
         *  @param  pid the input PID object
         *  @param  pdg the assumed pdg
         *  @param  view the view
         *  @param  dir if the fit is forward or backward
         *  @param  member the member to set
         */
        static void SetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const geo::View_t &view, const anab::kTrackDir &dir, Member<float> &member);

        /**
         *  @brief  Set the value of the input member variable to the bragg likelihood from the PID with given parameters using all planes
         *
         *  @param  pid the input PID object
         *  @param  pdg the assumed pdg
         *  @param  dir if the fit is forward or backward
         *  @param  yzAngle the angle in the YZ plane of the input track
         *  @param  sin2AngleThreshold threshold angle to collection wires to use colletion information
         *  @param  nHitsU the number of hits on the U plane
         *  @param  nHitsV the number of hits on the V plane
         *  @param  member the member to set
         */
        static void SetBraggLikelihood(const art::Ptr<anab::ParticleID> &pid, const int &pdg, const anab::kTrackDir &dir, const float yzAngle, const float sin2AngleThreshold, const unsigned int nHitsU, const unsigned int nHitsV, Member<float> &member);
};

} // namespace ubcc1pi

#endif
