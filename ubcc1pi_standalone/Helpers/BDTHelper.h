/**
 *  @file  ubcc1pi_standalone/Helpers/BDTHelper.h
 *
 *  @brief The header file for the bdt helper class
 */

#ifndef UBCC1PI_STANDALONE_HELPERS_BDT_HELPER
#define UBCC1PI_STANDALONE_HELPERS_BDT_HELPER

#include <string>
#include <vector>
#include <unordered_map>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include <TFile.h>

#include "ubcc1pi_standalone/Interface/Event.h"
#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The BDT helper class
 */
class BDTHelper
{
    public:

        typedef std::vector< std::pair<float, float> > BDTResponseWeightVector; ///< A vector of BDT responses and their corresponding weights

        /**
         *  @brief  A mapping with keys [missingFeature][particleType] and value of a vector of BDT responses as pair <bdtResponse, eventWeight>
         *          Used by the N-1 study
         */
        typedef std::map< std::string, std::map< PlottingHelper::PlotStyle, BDTResponseWeightVector > > BDTResponseMap;

        /**
         *  @brief  Base wrapper class for at root TMVA BDT
         */
        class BDTBase
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  bdtName the name of the BDT
                 *  @param  featureNames the vector of feature variable names
                 */
                BDTBase(const std::string &bdtName, const std::vector<std::string> &featureNames);

                // Destructor
                virtual ~BDTBase() = default;

                /**
                 *  @brief  Get the name
                 *
                 *  @return BDT name
                 */
                std::string GetName() const;

            protected:

                std::string               m_bdtName;      ///< The BDT name
                std::vector<std::string>  m_featureNames; ///< The feature variable names
        };

        /**
         *  @brief  Event shuffler class selects a random set of event indices for training and for testing
         */
        class EventShuffler
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  nEvents the total number of events to shuffle
                 *  @param  nTraining the number of events that should be used for training
                 */
                EventShuffler(const unsigned int nEvents, const unsigned int nTraining);

                /**
                 *  @brief  Determine if the input event should be used for training or testing
                 *
                 *  @param  eventIndex the input eventIndex
                 *
                 *  @return boolean, true if event should be used for training
                 */
                bool IsTrainingEvent(const unsigned int eventIndex);

            private:
                std::vector<bool> m_isTrainingEvent; ///< If the event with a given index should be used for training
        };

        /**
         *  @brief  BDT factory class used to train BDTs
         */
        class BDTFactory : public BDTBase
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  bdtName the name of the BDT
                 *  @param  featureNames the vector of feature variable names
                 */
                BDTFactory(const std::string &bdtName, const std::vector<std::string> &featureNames);

                /**
                 *  @brief  Destructor
                 */
                ~BDTFactory();

                /**
                 *  @brief  Add an entry
                 *
                 *  @param  features the features of the entry
                 *  @param  isSignal if the entry is signal or background
                 *  @param  isTraining if the entry is for training or testing
                 *  @param  weight the weight of the entry
                 */
                void AddEntry(const std::vector<float> &features, const bool isSignal, const bool isTraining, const float weight = 1.f);

                /**
                 *  @brief  Add a signal training entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddSignalTrainingEntry(const std::vector<float> &features, const float weight = 1.f);

                /**
                 *  @brief  Add a background training entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddBackgroundTrainingEntry(const std::vector<float> &features, const float weight = 1.f);

                /**
                 *  @brief  Add a signal test entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddSignalTestEntry(const std::vector<float> &features, const float weight = 1.f);

                /**
                 *  @brief  Add a background test entry
                 *
                 *  @param  features the features of the entry
                 *  @param  weight the weight of the entry
                 */
                void AddBackgroundTestEntry(const std::vector<float> &features, const float weight = 1.f);

                /**
                 *  @brief  Optimize the BDT parameters
                 */
                void OptimizeParameters();

                /**
                 *  @brief  Train and test the BDT using the entries added so far
                 */
                void TrainAndTest();

                /**
                 *  @brief  Get the feature importance for a trained BDT
                 * 
                 * @param  featureImportance the output feature importance histogram
                 * 
                 */
                void GetFeatureImportance(TH1F *featureImportance);

            private:

                /**
                 *  @brief  Book the BDT method
                 */
                void BookMethod();

                TFile             *m_pOutputFile; ///< The TMVA output file
                TMVA::Factory     *m_pFactory;    ///< The TMVA factory
                TMVA::DataLoader  *m_pDataLoader; ///< The data loader
                bool               m_isBooked;    ///< If the BDT has been booked

                unsigned int       m_nSignalTraining;     ///< The number of signal training events
                unsigned int       m_nSignalTest;         ///< The number of signal testing events
                unsigned int       m_nBackgroundTraining; ///< The number of background training events
                unsigned int       m_nBackgroundTest;     ///< The number of background testing events
        };

        /**
         *  @brief  The class for using a trained BDT
         */
        class BDT : public BDTBase
        {
            public:
                /**
                 *  @brief  Constructor
                 *
                 *  @param  bdtName the name of the BDT
                 *  @param  featureNames the vector of feature variable names
                 */
                BDT(const std::string &bdtName, const std::vector<std::string> &featureNames);

                // Destructor
                ~BDT() override;

                /**
                 *  @brief  Get the BDT response for a given set of features
                 *
                 *  @param  features the input vector of features
                 *
                 *  @return the BDT response
                 */
                float GetResponse(const std::vector<float> &features);

            private:
                TMVA::Reader      *m_pReader;     ///< The TMVA reader
                std::string        m_weightsFile; ///< The weights file
                std::vector<float> m_features;    ///< The BDT features
        };

        /**
         *  @brief  The particle BDT feature names
         */
        static const std::vector<std::string> ParticleBDTFeatureNames;

        /**
         *  @brief  The particle BDT feature names for the muon BDT
         */
        static const std::vector<std::string> MuonBDTFeatureNames;

        /**
         *  @brief  The particle BDT feature names for the proton BDT
         */
        static const std::vector<std::string> ProtonBDTFeatureNames;

        /**
         *  @brief  The particle BDT feature names for the golden pion BDT
         */
        static const std::vector<std::string> GoldenPionBDTFeatureNames;

        /**
         *  @brief  Get the features corresponding to the given names for the input reco particle
         *
         *  @param  recoParticle the input reco particle
         *  @param  featureNames the feature names to calculate
         *  @param  features the output features
         *  @param  shouldDebug if we should print debug statements
         *
         *  @return boolean, true if all features are available
         */
        static bool GetBDTFeatures(const Event::Reco::Particle &recoParticle, const std::vector<std::string> &featureNames, std::vector<float> &features, const bool shouldDebug = false);

        /**
         *  @brief  Get the ROC curve for a set of BDT responses for signal and background
         *
         *  @param  signalResponses the BDT responses of signal
         *  @param  backgroundResponses the BDT responses of background
         *  @param  nSamplePoints the number of sample points (values of the BDT response)
         *
         *  @return the ROC curve
         */
        static std::shared_ptr<TGraph> GetROCCurve(const BDTResponseWeightVector &signalResponses, const BDTResponseWeightVector &backgroundResponses, const unsigned int nSamplePoints = 100u);
        /**
         *  @brief  Get the ROC curve for a set of BDT responses for signal and background
         *
         *  @param  signalResponses the input BDT responses of signal
         *  @param  backgroundResponses the input BDT responses of background
         *  @param  nSamplePoints the number of sample points (values of the BDT response)
         *  @param  signalPassingRates the output signal passing rates
         *  @param  backgroundRejectionRates the output background rejection rates
         *  @param  signalPassingRateErrors the output uncertainties on the signal passing rates
         *  @param  backgroundRejectionRateErrors the output uncertainties onf the background rejection rates
         */
        static void GetROCCurve(const BDTResponseWeightVector &signalResponses, const BDTResponseWeightVector &backgroundResponses, const unsigned int nSamplePoints, std::vector<float> &signalPassingRates, std::vector<float> &backgroundRejectionRates, std::vector<float> &signalPassingRateErrors, std::vector<float> &backgroundRejectionRateErrors);

        /**
         *  @brief  Get the integral under the ROC curve
         *
         *  @param  signalPassingRates the input signal passing rates
         *  @param  backgroundRejectionRates the input background rejection rates
         *
         *  @return the ROC integral
         */
        static float GetROCIntegral(const std::vector<float> &signalPassingRates, const std::vector<float> &backgroundRejectionRates);

        /**
         *  @brief  Get the uncertainty on the ROC curve integral
         *
         *  @param  signalPassingRates the input signal passing rates
         *  @param  backgroundRejectionRates the input background rejection rates
         *  @param  signalPassingRateErrs the uncertainty on the input signal passing rates
         *  @param  backgroundRejectionRateErrs the uncertainty on the input background rejection rates
         *  @param  nUniverses the number of universes to use when calculating the uncertainty
         *
         *  @return The uncertainty on the ROC curve integral
         */
        static float GetROCIntegralError(const std::vector<float> &signalPassingRates, const std::vector<float> &backgroundRejectionRates, const std::vector<float> &signalPassingRateErrs, const std::vector<float> &backgroundRejectionRateErrs, const unsigned int nUniverses = 1000u);

};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<std::string> BDTHelper::ParticleBDTFeatureNames = {
    "logBragg_pToMIP",
    "logBragg_piToMIP",
    "truncMeandEdx",
    "protonForward",
    "muonForward",
    "nDescendents",
    "nSpacePointsNearEnd",
    "wiggliness",
    "trackScore"
};

// -----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<std::string> BDTHelper::MuonBDTFeatureNames = {
    "logBragg_pToMIP",
    "logBragg_piToMIP",
    "truncMeandEdx",
    "nDescendents",
    "wiggliness",
    "trackScore"
};

// -----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<std::string> BDTHelper::ProtonBDTFeatureNames = {
    "logBragg_pToMIP",
    "logBragg_piToMIP",
    "truncMeandEdx",
    "wiggliness",
    "trackScore"
};

// -----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<std::string> BDTHelper::GoldenPionBDTFeatureNames = {
    "logBragg_pToMIP",
    "logBragg_piToMIP",
    "truncMeandEdx",
    "nDescendents",
    "wiggliness",
    "trackScore"
};

} // namespace ubcc1pi

#endif
