#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include "ubcc1pi_standalone/Helpers/GeometryHelper.h"

#include <stdexcept>
#include <algorithm>

namespace ubcc1pi
{
        
std::string AnalysisHelper::GetSampleTypeName(const SampleType &sampleType)
{
    switch (sampleType)
    {
        case DataBNB:
            return "Data BNB";
        case Overlay:
            return "Overlay";
        case DataEXT:
            return "Data EXT";
        case Dirt:
            return "Dirt";
        default: break;
    }
            
    throw std::invalid_argument("AnalysisHelper::GetSampleTypeName - unknown sample type");
}
    
// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::CountEvent(const std::string &tag, const SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const float weight)
{
    // Keep track of this tag if we haven't seen it before
    if (std::find(m_tags.begin(), m_tags.end(), tag) == m_tags.end())
        m_tags.push_back(tag);

    // Get the mapping from the classification string to the count/weight, making it if it doesn't exist for this tag & sampleType
    auto &stringToCountMap = m_eventCountMap[tag][sampleType];
    auto &stringToWeightMap = m_eventWeightMap[tag][sampleType];
        
    const auto classification = AnalysisHelper::GetClassificationString(pEvent);

    // Add to the count map
    auto countIter = stringToCountMap.find(classification);
    if (countIter == stringToCountMap.end())
    {
        stringToCountMap.emplace(classification, 1u);
    }
    else
    {
        countIter->second++;
    }

    // Add to the weight map
    auto weightIter = stringToWeightMap.find(classification);
    if (weightIter == stringToWeightMap.end())
    {
        stringToWeightMap.emplace(classification, weight);
    }
    else
    {
        weightIter->second += weight;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::PrintBreakdown(const unsigned int nBackgrounds) const
{
    this->PrintBreakdown(m_eventCountMap, nBackgrounds);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::PrintBreakdownWithWeights(const unsigned int nBackgrounds) const
{
    this->PrintBreakdown(m_eventWeightMap, nBackgrounds);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsTrueCC1Pi(const std::shared_ptr<Event> &pEvent)
{
    if (!pEvent->metadata.hasTruthInfo())
        throw std::invalid_argument("AnalysisHelper::IsTrueCC1Pi - Input event doesn't have truth information!");

    const auto truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return false;
    
    // Count the visible particles
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
    const auto nMu = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 13);
    const auto nProton = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 2212);
    const auto nPion = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 211);
    const auto nOther = visibleParticles.size() - (nMu + nProton + nPion);

    // Insist the CC1Pi topology
    return (nMu == 1 && nPion == 1 && nOther == 0);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
    
bool AnalysisHelper::PassesVisibilityThreshold(const Event::Truth::Particle &particle)
{
    // First only consider specific particle types that can lead to hits
    const auto absPDG = std::abs(particle.pdgCode());
    const auto isVisible = (absPDG == 11   || // Electron
                            absPDG == 13   || // Muon
                            absPDG == 2212 || // Proton
                            absPDG == 2112 || // Neutron
                            absPDG == 22   || // Photon
                            absPDG == 211  || // Charged pion
                            absPDG == 321);   // Charged kaon

    if (!isVisible)
        return false;

    // Apply the momentum threhsolds
    float thresholdMomentum = -std::numeric_limits<float>::max();
    switch (absPDG)
    {
        case 2112: // Neutron
            thresholdMomentum = std::numeric_limits<float>::max(); // ATTN infinite threshold = never let a neutron pass
            break;
        default: break;
    }

    return (particle.momentum() >= thresholdMomentum);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
std::vector<Event::Truth::Particle> AnalysisHelper::SelectVisibleParticles(const std::vector<Event::Truth::Particle> &particles)
{
    std::vector<Event::Truth::Particle> visibleParticles;

    for (const auto &particle : particles)
    {
        if (AnalysisHelper::PassesVisibilityThreshold(particle))
            visibleParticles.push_back(particle);
    }

    return visibleParticles;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
unsigned int AnalysisHelper::CountParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode)
{
    unsigned int count = 0;

    for (const auto &particle : particles)
    {
        if (particle.pdgCode() == pdgCode)
            count++;
    }

    return count;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::CountGoldenParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode)
{
    std::vector<Event::Truth::Particle> goldenParticles;
    for (const auto &particle : particles)
    {
        if (AnalysisHelper::IsGolden(particle))
            goldenParticles.push_back(particle);
    }

    return AnalysisHelper::CountParticlesWithPdgCode(goldenParticles, pdgCode);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
void AnalysisHelper::GetPdgCodeCountMap(const std::vector<Event::Truth::Particle> &particles, std::vector<int> &foundPdgs, std::unordered_map<int, unsigned int> &pdgCodeCountMap)
{
    if (!foundPdgs.empty())
        throw std::invalid_argument("AnalysisHelper::GetPdgCodeCountMap - input foundPdgs vector isn't empty");

    if (!pdgCodeCountMap.empty())
        throw std::invalid_argument("AnalysisHelper::GetPdgCodeCountMap - input pdgCodeCountMap isn't empty");

    // Get the mapping
    for (const auto &particle : particles)
    {
        const auto pdg = particle.pdgCode();
        auto iter = pdgCodeCountMap.find(pdg);

        if (iter == pdgCodeCountMap.end())
        {
            // Add a new entry
            pdgCodeCountMap.emplace(pdg, 1);
            foundPdgs.push_back(pdg);
        }
        else
        {
            // Increment the counter
            iter->second++;
        }
    }

    // Sort the vector for reproducibility. Here we sort numerically in increasing order, but place particles and antiparticles next to each other
    std::sort(foundPdgs.begin(), foundPdgs.end(), [](const int &a, const int &b) {

        const auto aAbs = std::abs(a);
        const auto bAbs = std::abs(b);

        return (aAbs == bAbs) ? a > b : aAbs < bAbs;
    });
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string AnalysisHelper::GetTopologyString(const std::vector<Event::Truth::Particle> &particles, const bool countProtonsInclusively)
{
    // Count the particles by PDG code
    std::vector<int> foundPdgs;
    std::unordered_map<int, unsigned int> pdgCodeCountMap;
    AnalysisHelper::GetPdgCodeCountMap(particles, foundPdgs, pdgCodeCountMap);

    unsigned int nOther = 0;
    std::string topology = "";

    for (const auto &pdg : foundPdgs)
    {
        const auto count = pdgCodeCountMap.at(pdg);
        const auto countStr = std::to_string(count);
    
        switch (pdg)
        {
            case 11:
                topology += countStr + " e-  ";
                break;
            case -11:
                topology += countStr + " e+  ";
                break;
            case 13:
                topology += countStr + " Mu-  ";
                break;
            case -13:
                topology += countStr + " Mu+  ";
                break;
            case 22:
                topology += countStr + " Gamma  ";
                break;
            case 211:
                topology += countStr + " Pi+  ";
                break;
            case -211:
                topology += countStr + " Pi-  ";
                break;
            case 321:
                topology += countStr + " K+  ";
                break;
            case -321:
                topology += countStr + " K-  ";
                break;
            case 2112:
                topology += countStr + " n  ";
                break;
            case 2212:
                if (!countProtonsInclusively)
                    topology += countStr + " p  ";
                break;
            default:
                nOther++;
                break;
        }
    }

    if (countProtonsInclusively)
        topology += "X p  ";

    if (nOther != 0)
        topology += std::to_string(nOther) + " other  ";


    return topology;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::string AnalysisHelper::GetClassificationString(const std::shared_ptr<Event> &pEvent, const bool countProtonsInclusively)
{
    // Check if data
    if (!pEvent->metadata.hasTruthInfo())
        return "D,  ";

    const auto truth = pEvent->truth;
    
    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return "NF, ";

    std::string classification = "";
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(truth.particles);
    
    // Signal or background
    const auto isTrueCC1Pi = AnalysisHelper::IsTrueCC1Pi(pEvent);
    classification += isTrueCC1Pi ? "S" : "B,  ";
    
    if (isTrueCC1Pi)
    {
        // Check if we have a golden pion
        const auto hasGoldenPion = (AnalysisHelper::CountGoldenParticlesWithPdgCode(visibleParticles, 211) != 0);

        if (hasGoldenPion)
            classification += " G,";
        else
            classification += ",  ";
    }

    // Add the topology classification
    classification += "  " + AnalysisHelper::GetTopologyString(visibleParticles, countProtonsInclusively);

    return classification;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsPointWithinMargins(const TVector3 &point, const float lowXMargin, const float highXMargin, const float lowYMargin, const float highYMargin, const float lowZMargin, const float highZMargin)
{
    return ((point.x() > GeometryHelper::lowX + lowXMargin)   && 
            (point.x() < GeometryHelper::highX - highXMargin) &&
            (point.y() > GeometryHelper::lowY + lowYMargin)   &&
            (point.y() < GeometryHelper::highY - highYMargin) &&
            (point.z() > GeometryHelper::lowZ + lowZMargin)   &&
            (point.z() < GeometryHelper::highZ - highZMargin) );
}
        
// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsFiducial(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 10.f, 10.f, 10.f, 10.f, 10.f, 50.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 5.f, 5.f, 5.f, 5.f, 5.f, 5.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const Event::Truth::Particle &particle)
{
    const auto start = TVector3(particle.startX(), particle.startY(), particle.startZ());
    const auto end = TVector3(particle.endX(), particle.endY(), particle.endZ());

    return (AnalysisHelper::IsContained(start) && AnalysisHelper::IsContained(end));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::HasTrackFit(const Event::Reco::Particle &particle)
{
    return (particle.startX.IsSet() && particle.startY.IsSet() && particle.startZ.IsSet() && particle.endX.IsSet() && particle.endY.IsSet() && particle.endZ.IsSet());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsContained(const Event::Reco::Particle &particle)
{
    if (!AnalysisHelper::HasTrackFit(particle))
        throw std::invalid_argument("AnalysisHelper::IsContained - input reco particle doesn't have a fitted track start-end points");

    const auto start = TVector3(particle.startX(), particle.startY(), particle.startZ());
    const auto end = TVector3(particle.endX(), particle.endY(), particle.endZ());

    return (AnalysisHelper::IsContained(start) && AnalysisHelper::IsContained(end));
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetBestMatchedTruthParticleIndex(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    if (!recoParticle.hasMatchedMCParticle.IsSet())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle is from an event without truth info");

    if (!recoParticle.hasMatchedMCParticle())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle has no matched truth particle");

    if (truthParticles.size() != recoParticle.truthMatchPurities().size())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - wrong number of input truth particles");

    if (truthParticles.empty())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndex - no truth particles supplied but reco particle has a match!");


    unsigned int bestMatchedParticleId = std::numeric_limits<unsigned int>::max();
    float bestMatchScore = -std::numeric_limits<float>::max();
    bool foundMatch = false;

    for (unsigned int i = 0; i < truthParticles.size(); ++i)
    {
        if (!AnalysisHelper::PassesVisibilityThreshold(truthParticles.at(i)))
            continue;

        const auto score = recoParticle.truthMatchPurities().at(i);
        if (score < bestMatchScore)
            continue;

        bestMatchScore = score;
        bestMatchedParticleId = i;
        foundMatch = true;
    }

    if (!foundMatch)
        throw std::logic_error("AnalysisHelper::GetBestMatchedTruthParticleIndex - input reco particle has no matched truth particle");
    
    return bestMatchedParticleId;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

Event::Truth::Particle AnalysisHelper::GetBestMatchedTruthParticle(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    return truthParticles.at(AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, truthParticles, applyVisibilityThreshold));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
            
bool AnalysisHelper::IsGolden(const Event::Truth::Particle &particle)
{
    return (particle.nElasticScatters() == 0 &&
            particle.nInelasticScatters() == 0 &&
            particle.isStopping() &&
            AnalysisHelper::IsContained(particle));
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
bool AnalysisHelper::GetLikelihoodRatio(const Member<float> &numerator, const Member<float> &denominator, float &ratio)
{
    ratio = -std::numeric_limits<float>::max();

    if (!numerator.IsSet() || !denominator.IsSet())
        return false;

    if (std::abs(denominator()) <= std::numeric_limits<float>::epsilon())
        return false;

    ratio = numerator() / denominator();
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
bool AnalysisHelper::GetSoftmax(const Member<float> &signal, const Member<float> &background, float &softmax)
{
    softmax = -std::numeric_limits<float>::max();

    if (!signal.IsSet() || !background.IsSet())
        return false;

    const auto denominator = std::exp(signal()) + std::exp(background());
    if (std::abs(denominator) <= std::numeric_limits<float>::epsilon())
        return false;

    softmax = std::exp(signal()) / denominator;
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::PrintLoadingBar(const unsigned int numerator, const unsigned int denominator)
{
    if (denominator == 0)
        return;

    const auto width = 85u;
    const auto loadedWidth = static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator) / static_cast<float>(denominator))));
    const auto unloadedWidth = static_cast<unsigned int>(width - loadedWidth);
    
    // Work out if it's worth printing 
    const bool shouldPrint = (numerator == 0) || (static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator - 1) / static_cast<float>(denominator)))) != loadedWidth);

    if (!shouldPrint)
        return;

    if (numerator != 0)
        std::cout << "\033[F\033[F";

    std::cout << numerator << " / " << denominator << std::endl;
    std::cout << "|" << std::string(loadedWidth, '=') << std::string(unloadedWidth, ' ') << "|" << std::endl;
}

} // namespace ubcc1pi
