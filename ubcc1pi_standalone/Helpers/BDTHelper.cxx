/**
 *  @file  ubcc1pi_standalone/Helpers/BDTHelper.cxx
 *
 *  @brief Implementation of the BDT helper class
 */

#include "ubcc1pi_standalone/Helpers/BDTHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include <TString.h>
#include <random>
#include <algorithm>

namespace ubcc1pi
{

BDTHelper::EventShuffler::EventShuffler(const unsigned int nEvents, const unsigned int nTraining)
{
    if (nTraining > nEvents)
        throw std::invalid_argument("BDTHelper::EventShuffler::EventShuffler - Number of training events greater than total number of events");

    // Setup the training and testing events
    const auto nTesting = nEvents - nTraining;
    m_isTrainingEvent.insert(m_isTrainingEvent.end(), nTraining, true);
    m_isTrainingEvent.insert(m_isTrainingEvent.end(), nTesting, false);

    // Randomize
    std::shuffle(m_isTrainingEvent.begin(), m_isTrainingEvent.end(), std::default_random_engine());
}
// -----------------------------------------------------------------------------------------------------------------------------------------

bool BDTHelper::EventShuffler::IsTrainingEvent(const unsigned int eventIndex)
{
    if (eventIndex >= m_isTrainingEvent.size())
        throw std::invalid_argument("BDTHelper::EventShuffler::IsTrainingEvent - Input eventIndex out of range");

    return m_isTrainingEvent.at(eventIndex);
}


// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

BDTHelper::BDTBase::BDTBase(const std::string &bdtName, const std::vector<std::string> &featureNames) :
    m_bdtName(bdtName),
    m_featureNames(featureNames)
{
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string BDTHelper::BDTBase::GetName() const
{
    return m_bdtName;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

BDTHelper::BDTFactory::BDTFactory(const std::string &bdtName, const std::vector<std::string> &featureNames) :
    BDTBase(bdtName, featureNames),
    m_pOutputFile(nullptr),
    m_pFactory(nullptr),
    m_isBooked(false),
    m_nSignalTraining(0u),
    m_nSignalTest(0u),
    m_nBackgroundTraining(0u),
    m_nBackgroundTest(0)
{
    // Setup the output file and make the TMVA factory
    m_pDataLoader = new TMVA::DataLoader((bdtName + "_dataset").c_str());

    // Add the features to the data loader
    for (const auto &feature : featureNames)
        m_pDataLoader->AddVariable(TString(feature));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

BDTHelper::BDTFactory::~BDTFactory()
{
    if (m_pOutputFile)
        m_pOutputFile->Close();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddEntry(const std::vector<float> &features, const bool isSignal, const bool isTraining, const float weight)
{
    if (isSignal && isTraining)
        this->AddSignalTrainingEntry(features, weight);

    if (!isSignal && isTraining)
        this->AddBackgroundTrainingEntry(features, weight);

    if (isSignal && !isTraining)
        this->AddSignalTestEntry(features, weight);

    if (!isSignal && !isTraining)
        this->AddBackgroundTestEntry(features, weight);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddSignalTrainingEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddSignalTrainingEntry - Wrong number of features");

    m_pDataLoader->AddSignalTrainingEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nSignalTraining++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddBackgroundTrainingEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddBackgroundTrainingEntry - Wrong number of features");

    m_pDataLoader->AddBackgroundTrainingEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nBackgroundTraining++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddSignalTestEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddSignalTestEntry - Wrong number of features");

    m_pDataLoader->AddSignalTestEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nSignalTest++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::AddBackgroundTestEntry(const std::vector<float> &features, const float weight)
{
    if (features.size() != m_featureNames.size())
        throw std::invalid_argument("BDTHelper::BDTFactory::AddBackgroundTestEntry - Wrong number of features");

    m_pDataLoader->AddBackgroundTestEvent(std::vector<double>(features.begin(), features.end()), weight);
    m_nBackgroundTest++;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::BookMethod()
{
    if (!m_isBooked)
    {
        m_pOutputFile = TFile::Open((m_bdtName + "_BDT.root").c_str(), "RECREATE" );
        m_pFactory = new TMVA::Factory(TString(m_bdtName), m_pOutputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

        m_pDataLoader->PrepareTrainingAndTestTree("0==0", m_nSignalTraining, m_nBackgroundTraining, m_nSignalTest, m_nBackgroundTest);
        m_pFactory->BookMethod(m_pDataLoader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=505:MinNodeSize=1.26436%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

        m_isBooked = true;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::OptimizeParameters()
{
    this->BookMethod();

    if (!m_pOutputFile || !m_pFactory)
        throw std::logic_error("BDTFactory::OptimizeParameters - Invalid output file or TMVA factory");

    m_pFactory->OptimizeAllMethodsForClassification();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::BDTFactory::TrainAndTest()
{
    this->BookMethod();

    if (!m_pOutputFile || !m_pFactory)
        throw std::logic_error("BDTFactory::OptimizeParameters - Invalid output file or TMVA factory");

    std::cout << "DEBUG - TrainAndTest - Training BDT" << std::endl;
    m_pFactory->TrainAllMethods();
    std::cout << "DEBUG - TrainAndTest - Testing BDT" << std::endl;
    m_pFactory->TestAllMethods();
    std::cout << "DEBUG - TrainAndTest - Evaluating BDT" << std::endl;
    m_pFactory->EvaluateAllMethods();
    std::cout << "DEBUG - TrainAndTest - Done" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

BDTHelper::BDT::BDT(const std::string &bdtName, const std::vector<std::string> &featureNames) :
    BDTBase(bdtName, featureNames),
    m_pReader(new TMVA::Reader("Silent")),
    m_weightsFile(bdtName + "_dataset/weights/" + bdtName + "_BDT.weights.xml"),
    m_features(featureNames.size(), -std::numeric_limits<float>::max())
{
    // Bind the feature vector to the reader
    for (unsigned int iFeature = 0; iFeature < featureNames.size(); ++iFeature)
        m_pReader->AddVariable(featureNames.at(iFeature).c_str(), &m_features.at(iFeature));

    // Book the MVA using the weights file
    m_pReader->BookMVA("BDT", m_weightsFile.c_str());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BDTHelper::BDT::GetResponse(const std::vector<float> &features)
{
    if (features.size() != m_features.size())
        throw std::invalid_argument("BDTHelper::BDT::GetResponse - Wrong number of features, expected " + std::to_string(m_features.size()) + " but passed " + std::to_string(features.size()));

    // Shallow copy the features
    for (unsigned int iFeature = 0; iFeature < m_features.size(); ++iFeature)
        m_features.at(iFeature) = features.at(iFeature);

    // Run the BDT
    return m_pReader->EvaluateMVA("BDT");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

bool BDTHelper::GetBDTFeatures(const Event::Reco::Particle &recoParticle, const std::vector<std::string> &featureNames, std::vector<float> &features, const bool shouldDebug)
{
    if (!features.empty())
        throw std::logic_error("BDTHelper::GetBDTFeatures - Input feature vector isn't empty");

    for (const auto &name : featureNames)
    {
        if (shouldDebug)
            std::cout << "DEBUG - Calculating: " << name << std::endl;

        if (name == "logBragg_pToMIP")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetLogLikelihoodRatio(recoParticle.likelihoodForwardProton, recoParticle.likelihoodMIP, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "logBragg_piToMIP")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetLogLikelihoodRatio(recoParticle.likelihoodForwardPion, recoParticle.likelihoodMIP, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "truncMeandEdx")
        {
            if (!recoParticle.truncatedMeandEdx.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(recoParticle.truncatedMeandEdx());
            continue;
        }

        if (name == "protonForward")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetSoftmax(recoParticle.likelihoodForwardProton, recoParticle.likelihoodBackwardProton, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "muonForward")
        {
            float feature = -std::numeric_limits<float>::max();
            if (!AnalysisHelper::GetSoftmax(recoParticle.likelihoodForwardMuon, recoParticle.likelihoodBackwardMuon, feature))
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(feature);
            continue;
        }

        if (name == "nDescendents")
        {
            if (!recoParticle.nDescendents.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(recoParticle.nDescendents());
            continue;
        }

        if (name == "nSpacePointsNearEnd")
        {
            if (!recoParticle.nSpacePointsNearEnd.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(recoParticle.nSpacePointsNearEnd());
            continue;
        }

        if (name == "wiggliness")
        {
            if (!recoParticle.wiggliness.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(recoParticle.wiggliness());
            continue;
        }

        if (name == "trackScore")
        {
            if (!recoParticle.trackScore.IsSet())
            {
                if (shouldDebug)
                    std::cout << "DEBUG - Can't calculate: " << name << std::endl;

                features.clear();
                return false;
            }

            features.push_back(recoParticle.trackScore());
            continue;
        }

        throw std::invalid_argument("BDTHelper::GetBDTFeatures - Unknown feature: \"" + name + "\"");
    }

    // This should never happen - just a sanity check
    if (features.size() != featureNames.size())
        throw std::logic_error("BDTHelper::GetBDTFeatures - Number of features calculated is incorrect!");

    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void BDTHelper::GetROCCurve(const BDTResponseWeightVector &signalResponses, const BDTResponseWeightVector &backgroundResponses, const unsigned int nSamplePoints, std::vector<float> &signalPassingRates, std::vector<float> &backgroundRejectionRates, std::vector<float> &signalPassingRateErrors, std::vector<float> &backgroundRejectionRateErrors)
{
    if (signalResponses.empty() || backgroundResponses.empty())
        throw std::invalid_argument("BDTHelper::GetROCCurve - input vector of responses is empty");

    if (nSamplePoints < 2)
        throw std::invalid_argument("BDTHelper::GetROCCurve - need at least 2 sample points");

    if (!signalPassingRates.empty() || !backgroundRejectionRates.empty() || !signalPassingRateErrors.empty() || !backgroundRejectionRateErrors.empty())
        throw std::invalid_argument("BDTHelper::GetROCCurve - output vectors need to be empty");

    // Define a comparison function for an element of a BDTResponseWeightVector - it returns true if the first response is less than the second
    const auto &comp = [](const auto &a, const auto &b) { return a.first < b.first; };

    // Sort the responses
    auto sortedSignalResponses = signalResponses;
    std::sort(sortedSignalResponses.begin(), sortedSignalResponses.end(), comp);

    auto sortedBackgroundResponses = backgroundResponses;
    std::sort(sortedBackgroundResponses.begin(), sortedBackgroundResponses.end(), comp);

    // Get the maximum and minimum BDT responses
    const auto epsilon = std::numeric_limits<float>::epsilon();
    const auto minX = std::min(sortedSignalResponses.front(), sortedBackgroundResponses.front()).first - epsilon;
    const auto maxX = std::max(sortedSignalResponses.back(), sortedBackgroundResponses.back()).first + epsilon;

    // Sample between minX and maxX
    std::vector<float> vectS, vectB;
    unsigned int iSignal = 0u;
    unsigned int iBackground = 0u;
    float signalWeight = 0.f;
    float backgroundWeight = 0.f;

    for (unsigned int iSample = 0; iSample < nSamplePoints; ++iSample)
    {
        const auto X = minX + (maxX - minX) * (iSample / static_cast<float>(nSamplePoints - 1));

        // Sum up the signal and background responses until we find one that is greater than the current cut value
        while (iSignal < sortedSignalResponses.size())
        {
            if (sortedSignalResponses.at(iSignal).first > X)
                break;

            signalWeight += sortedSignalResponses.at(iSignal).second;
            iSignal++;
        }

        while (iBackground < sortedBackgroundResponses.size())
        {
            if (sortedBackgroundResponses.at(iBackground).first > X)
                break;

            backgroundWeight += sortedBackgroundResponses.at(iBackground).second;
            iBackground++;
        }

        // Store the values at this point
        vectS.push_back(signalWeight);
        vectB.push_back(backgroundWeight);
    }

    // Sanity check that we have seen every BDT response
    if (iSignal != sortedSignalResponses.size())
        throw std::logic_error("BDTHelper::GetROCCurve - sanity check failed, haven't used all signal responses");

    if (iBackground != sortedBackgroundResponses.size())
        throw std::logic_error("BDTHelper::GetROCCurve - sanity check failed, haven't used all background responses");

    const auto signalWeightTotal = signalWeight;
    const auto backgroundWeightTotal = backgroundWeight;

    // Get the fractions used in the ROC curve
    for (unsigned int iSample = 0; iSample < nSamplePoints; ++iSample)
    {
        // The signal passing rate
        const auto fracS = 1.f - (vectS.at(iSample) / signalWeightTotal);
        const auto errS = AnalysisHelper::GetEfficiencyUncertainty(vectS.at(iSample), signalWeightTotal);

        // The background rejection rate
        const auto fracB = vectB.at(iSample) / backgroundWeightTotal;
        const auto errB = AnalysisHelper::GetEfficiencyUncertainty(vectB.at(iSample), backgroundWeightTotal);

        signalPassingRates.push_back(fracS);
        backgroundRejectionRates.push_back(fracB);

        signalPassingRateErrors.push_back(errS);
        backgroundRejectionRateErrors.push_back(errB);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BDTHelper::GetROCIntegral(const std::vector<float> &signalPassingRates, const std::vector<float> &backgroundRejectionRates)
{
    if (signalPassingRates.size() != backgroundRejectionRates.size())
        throw std::invalid_argument("BDTHelper::GetROCIntegral - sizes of input vectors doesn't match");

    const auto nSamples = signalPassingRates.size();
    if (nSamples < 2)
        throw std::invalid_argument("BDTHelper::GetROCIntegral - need at least 2 samples");

    // Sort the input data points by signal passing rate
    std::vector< std::pair<float, float> > data;
    for (unsigned int i = 0; i < nSamples; ++i)
    {
        data.emplace_back(signalPassingRates.at(i), backgroundRejectionRates.at(i));
    }
    std::sort(data.begin(), data.end(), [](const auto &a, const auto &b) { return a.first < b.first; });

    // Perform trapezium integration
    float integral = 0.f;
    for (unsigned int i = 1; i < nSamples; ++i)
    {
        const auto f1 = data.at(i).first;
        const auto f0 = data.at(i-1).first;
        const auto g1 = data.at(i).second;
        const auto g0 = data.at(i-1).second;

        const auto df = f1 - f0;
        const auto dg = g1 - g0;

        // Calculate the area of the rectangle between f0 -> f1 at height g0
        const auto rectArea = df * g0;

        // Calculate the area of the triangle between f0 -> f1 on top of this rectangle made between g0 -> g1
        const auto triArea = 0.5 * df * dg;

        // Add to the integral
        integral += rectArea + triArea;
    }

    return integral;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float BDTHelper::GetROCIntegralError(const std::vector<float> &signalPassingRates, const std::vector<float> &backgroundRejectionRates, const std::vector<float> &signalPassingRateErrs, const std::vector<float> &backgroundRejectionRateErrs, const unsigned nUniverses)
{
    const auto nSamples = signalPassingRates.size();
    if (backgroundRejectionRates.size() != nSamples || signalPassingRateErrs.size() != nSamples || backgroundRejectionRateErrs.size() != nSamples)
        throw std::invalid_argument("BDTHelper::GetROCIntegralError - sizes of input vectors doesn't match");

    if (nUniverses <= 1)
        throw std::invalid_argument("BDTHelper::GetROCIntegralError - need at least two universes");

    // Create normal distributions for each sample point in the ROC curve
    std::default_random_engine generator;
    std::vector< std::normal_distribution<float> > signalDistributions, backgroundDistributions;
    for (unsigned int i = 0; i < nSamples; ++i)
    {
        signalDistributions.emplace_back(signalPassingRates.at(i), signalPassingRateErrs.at(i));
        backgroundDistributions.emplace_back(backgroundRejectionRates.at(i), backgroundRejectionRateErrs.at(i));
    }

    // Here we use a MC technique to estimate the uncertainty
    std::vector<float> rocIntegrals;
    float summedROCIntegral = 0.f;

    for (unsigned int iUniv = 0; iUniv < nUniverses; ++iUniv)
    {
        // Generate signal passing rate and background rejection rate vectors with random 1 sigma variations
        std::vector<float> generatedSignalPassingRates, generatedBackgroundRejectionRates;

        for (unsigned int i = 0; i < nSamples; ++i)
        {
            generatedSignalPassingRates.push_back(signalDistributions.at(i)(generator));
            generatedBackgroundRejectionRates.push_back(backgroundDistributions.at(i)(generator));
        }

        // Now get the ROC integral for this universe
        const auto rocIntegral = BDTHelper::GetROCIntegral(generatedSignalPassingRates, generatedBackgroundRejectionRates);

        rocIntegrals.push_back(rocIntegral);
        summedROCIntegral += rocIntegral;
    }

    // Find the standard deviation of the generate ROC integrals
    const auto mean = summedROCIntegral / static_cast<float>(nUniverses);

    float summedSquareDifferences = 0.f;
    for (unsigned int iUniv = 0; iUniv < nUniverses; ++iUniv)
    {
        summedSquareDifferences += std::pow(rocIntegrals.at(iUniv) - mean, 2);
    }

    const auto variance = summedSquareDifferences / static_cast<float>(nUniverses - 1);
    const auto uncertainty = std::pow(variance, 0.5f);

    return uncertainty;
}

} // namespace ubcc1pi
