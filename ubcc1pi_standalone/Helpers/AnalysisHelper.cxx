/**
 *  @file  ubcc1pi_standalone/Helpers/AnalysisHelper.cxx
 *
 *  @brief The implementation of the analysis helper class
 */

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

#include "ubcc1pi_standalone/Helpers/GeometryHelper.h"
#include "ubcc1pi_standalone/Helpers/FormattingHelper.h"
#include "ubcc1pi_standalone/Helpers/NormalisationHelper.h"

#include <stdexcept>
#include <algorithm>
#include <numeric>

#include <TGraph.h>
#include <TFitResult.h>

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
        case NuWro:
            return "NuWro";
        case DataEXT:
            return "Data EXT";
        case Dirt:
            return "Dirt";
        case DetectorVariation:
            return "Detector variation";
        default: break;
    }

    throw std::invalid_argument("AnalysisHelper::GetSampleTypeName - unknown sample type");
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::EventCounter::EventCounter()
{
    // The "all" is treated differently when printing out
    m_tags.push_back("all");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetWeight(const std::string &tag, const SampleType &sampleType, const std::string &classification) const
{
    float weight = 0.f;
    if (this->GetWeight(tag, sampleType, classification, weight))
        return weight;

    return 0.f;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::EventCounter::GetWeight(const std::string &tag, const SampleType &sampleType, const std::string &classification, float &weight) const
{
    const auto tagIter = m_eventWeightMap.find(tag);
    if (tagIter == m_eventWeightMap.end())
        return false;

    const auto tagMap = tagIter->second;

    const auto sampleIter = tagMap.find(sampleType);
    if (sampleIter == tagMap.end())
        return false;

    const auto sampleMap = sampleIter->second;

    const auto classificationIter = sampleMap.find(classification);
    if (classificationIter == sampleMap.end())
        return false;

    weight = classificationIter->second;
    return true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetTotalMCWeight(const std::string &tag) const
{
    // Get the total weight for this tag
    float totalWeight = 0.f;
    for (const auto &type : AnalysisHelper::AllSampleTypes)
    {
        if (type == DataBNB)
            continue;

        for (const auto &classification : m_classifications)
            totalWeight += this->GetWeight(tag, type, classification);
    }

    return totalWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetBNBDataWeight(const std::string &tag) const
{
    float weight = 0.f;

    for (const auto &classification : m_classifications)
        weight += this->GetWeight(tag, DataBNB, classification);

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::EventCounter::IsSignalClassification(const std::string &classification, const std::string &searchQuery) const
{
    if (classification.empty())
        return false;

    // Signal classifications always begin with 'S'
    if (classification.at(0) != 'S')
        return false;

    if (searchQuery.empty())
        return true;

    // Insist that there is also the additional query
    return (classification.find(searchQuery) != std::string::npos);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::string> AnalysisHelper::EventCounter::GetSignalClassifications(const std::string &searchQuery) const
{
    std::vector<std::string> signalClassifications;

    for (const auto &classification : m_classifications)
    {
        if (!this->IsSignalClassification(classification, searchQuery))
            continue;

        signalClassifications.push_back(classification);
    }

    return signalClassifications;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetSignalWeight(const std::string &tag, const std::string &searchQuery) const
{
    float weight = 0.f;

    for (const auto &classification : m_classifications)
    {
        if (!this->IsSignalClassification(classification, searchQuery))
            continue;

        weight += this->GetWeight(tag, Overlay, classification);
    }

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetBackgroundWeight(const std::string &tag, const std::string &searchQuery) const
{
    float weight = 0.f;

    for (const auto &classification : m_classifications)
    {
        if (this->IsSignalClassification(classification, searchQuery))
            continue;

        for (const auto &type : AnalysisHelper::AllSampleTypes)
        {
            if (type == DataBNB)
                continue;

            weight += this->GetWeight(tag, type, classification);
        }
    }

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetSignalEfficiency(const std::string &tag, const std::string &searchQuery) const
{
    const auto allWeight = this->GetSignalWeight("all", searchQuery);

    if (allWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    const auto weight = this->GetSignalWeight(tag, searchQuery);

    return weight / allWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetSignalPurity(const std::string &tag, const std::string &searchQuery) const
{
    if (tag == "all")
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetSignalPurity - Can't get purity for tag \"all\", as it won't be available for off-beam data");

    const auto totalWeight = this->GetTotalMCWeight(tag);

    if (totalWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    const auto weight = this->GetSignalWeight(tag, searchQuery);

    return weight / totalWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetEfficiency(const std::string &tag, const SampleType &sampleType, const std::string &classification) const
{
    if (sampleType != Overlay)
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetEfficiency - Can't get efficiency of data!");

    const auto allWeight = this->GetWeight("all", sampleType, classification);

    if (allWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::max();

    const auto weight = this->GetWeight(tag, sampleType, classification);

    return weight / allWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::EventCounter::GetPurity(const std::string &tag, const SampleType &sampleType, const std::string &classification) const
{
    if (sampleType == DataBNB)
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetPurity - Can't get purity of BNB data!");

    if (tag == "all")
        throw std::invalid_argument("AnalysisHelper::EventCounter::GetPurity - Can't get purity for tag \"all\", as it won't be available for off-beam data");

    const auto totalWeight = this->GetTotalMCWeight(tag);

    if (totalWeight <= std::numeric_limits<float>::epsilon())
        return -std::numeric_limits<float>::epsilon();

    const auto weight = this->GetWeight(tag, sampleType, classification);

    return weight / totalWeight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::CountEvent(const std::string &tag, const SampleType &sampleType, const std::shared_ptr<Event> &pEvent, const float weight, const bool selectCC0pi, const float protonMomentumThreshold)
{
    // Keep track of this tag if we haven't seen it before
    if (std::find(m_tags.begin(), m_tags.end(), tag) == m_tags.end())
        m_tags.push_back(tag);

    // Classify the event
    const auto useAbsPdg = true; // TODO make this configurable
    const auto countProtonsInclusively = true; // TODO make this configurable
    const auto classification = AnalysisHelper::GetClassificationString(pEvent, useAbsPdg, countProtonsInclusively, selectCC0pi, protonMomentumThreshold);

    // Keep track of this classification if we haven't seen it before
    if (std::find(m_classifications.begin(), m_classifications.end(), classification) == m_classifications.end())
        m_classifications.push_back(classification);

    std::sort(m_classifications.begin(), m_classifications.end());

    // Get the mapping from the classification string to the weight, making it if it doesn't exist for this tag and sampleType
    auto &stringToWeightMap = m_eventWeightMap[tag][sampleType];
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

std::vector<std::string> AnalysisHelper::EventCounter::GetTags() const
{
    return m_tags;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::PrintBreakdownSummary(const std::string &outputFileName) const
{
    FormattingHelper::Table table({"Tag", "", "Signal", "Background", "Efficiency", "Purity", "E*P", "Golden fraction", "", "BNB Data", "Data/MC ratio"});
    for (const auto &tag : m_tags)
    {
        table.AddEmptyRow();
        table.SetEntry("Tag", tag);

        table.SetEntry("Signal", this->GetSignalWeight(tag));
        table.SetEntry("Background", this->GetBackgroundWeight(tag));

        const auto efficiency = this->GetSignalEfficiency(tag);
        table.SetEntry("Efficiency", efficiency);

        const auto isAll = (tag == "all");
        if (tag == "all")
        {
            table.SetEntry("Purity", "?");
            table.SetEntry("E*P", "?");
        }
        else
        {
            const auto purity = this->GetSignalPurity(tag);

            table.SetEntry("Purity", purity);
            table.SetEntry("E*P", efficiency * purity);
        }

        table.SetEntry("Golden fraction", this->GetSignalWeight(tag, "S G") / this->GetSignalWeight(tag));

        const auto bnbDataWeight = this->GetBNBDataWeight(tag);
        const auto totalMCWeight = this->GetTotalMCWeight(tag);

        table.SetEntry("BNB Data", bnbDataWeight);

        if (totalMCWeight <= std::numeric_limits<float>::epsilon())
        {
            table.SetEntry("Data/MC ratio", "?");
        }
        else
        {
            table.SetEntry("Data/MC ratio", bnbDataWeight / totalMCWeight);
        }
    }

    table.WriteToFile(outputFileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::EventCounter::PrintBreakdownDetails(const std::string &outputFileName, const unsigned int nEntries) const
{
    FormattingHelper::Table table({"Tag", "", "Sample", "Classification", "", "Weight", "Efficiency", "Purity"});

    for (const auto &tag : m_tags)
    {
        // Print a separator row
        if (tag != m_tags.front())
            table.AddEmptyRow();

        // The total MC weight we have and haven't printed
        float totalMCWeight = 0.f;
        float totalMCWeightPrinted = 0.f;

        for (const auto &sample : AnalysisHelper::AllSampleTypes)
        {
            const auto sampleName = AnalysisHelper::GetSampleTypeName(sample);
            std::vector<std::pair<std::string, float> > classificationWeightVector;

            for (const auto &classification : m_classifications)
            {
                float weight = 0.f;
                if (!this->GetWeight(tag, sample, classification, weight))
                    continue;

                classificationWeightVector.emplace_back(classification, weight);
            }

            std::sort(classificationWeightVector.begin(), classificationWeightVector.end(), [&](const auto &a, const auto &b){
                const bool isASignal = this->IsSignalClassification(a.first);
                const bool isBSignal = this->IsSignalClassification(b.first);

                // Put signal before background
                if (isASignal != isBSignal)
                    return isASignal;

                // Order in numerically, largest numbers first
                return a.second > b.second;
            });

            // Print the entries
            unsigned int entriesPrinted = 0;
            for (const auto &entry : classificationWeightVector)
            {
                const auto &classification = entry.first;
                const auto weight = entry.second;

                if (sample != AnalysisHelper::DataBNB)
                    totalMCWeight += weight;

                if (entriesPrinted == nEntries)
                    continue;

                if (sample != AnalysisHelper::DataBNB)
                    totalMCWeightPrinted += weight;

                table.AddEmptyRow();
                table.SetEntry("Tag", tag);
                table.SetEntry("Sample", sampleName);
                table.SetEntry("Classification", classification);
                table.SetEntry("Weight", weight);

                if (sample == Overlay)
                {
                    table.SetEntry("Efficiency", this->GetEfficiency(tag, sample, classification));
                }

                if (sample != DataBNB && tag != "all")
                {
                    table.SetEntry("Purity", this->GetPurity(tag, sample, classification));
                }

                entriesPrinted++;
            }
        }

        // Print the other catagory
        const auto mcWeightNotPrinted = totalMCWeight - totalMCWeightPrinted;
        table.AddEmptyRow();
        table.SetEntry("Tag", tag);
        table.SetEntry("Classification", "Other");
        table.SetEntry("Weight", mcWeightNotPrinted);

        if (totalMCWeight > std::numeric_limits<float>::epsilon())
        {
            const auto otherPurity = mcWeightNotPrinted / totalMCWeight;
            table.SetEntry("Purity", otherPurity);
        }

    }

    table.WriteToFile(outputFileName);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetNominalEventWeight(const std::shared_ptr<Event> &pEvent)
{
    float weight = 1.f;

    const auto &truth = pEvent->truth;

    // ATTN unfortunately when the event weights are applied a small fraction of events are given an infinite weight, I assume due to a
    // failure in the event weighting code - for now we just skip infinite weights
    if (truth.splineEventWeight.IsSet() && std::abs(truth.splineEventWeight()) < std::numeric_limits<float>::max())
        weight *= truth.splineEventWeight();

    if (truth.genieTuneEventWeight.IsSet() && std::abs(truth.genieTuneEventWeight()) < std::numeric_limits<float>::max())
        weight *= truth.genieTuneEventWeight();

    // if(!truth.splineEventWeight.IsSet() || !truth.genieTuneEventWeight.IsSet())
    //     std::cout << "DEBUG: splineEventWeight = " << truth.splineEventWeight.IsSet() << ", genieTuneEventWeight = " << truth.genieTuneEventWeight.IsSet() << std::endl;

    return weight;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsTrueCCInclusive(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg)
{
    if (!pEvent->metadata.hasTruthInfo())
        throw std::invalid_argument("AnalysisHelper::IsTrueCCInclusive - Input event doesn't have truth information!");

    const auto truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return false;

    // Count the visible particles
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
    const auto nMu = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 13, useAbsPdg);

    // Insist the CCInclusive topology
    return (nMu == 1);

}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsTrueCC1Pi(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg)
{
    if (!pEvent->metadata.hasTruthInfo())
        throw std::invalid_argument("AnalysisHelper::IsTrueCC1Pi - Input event doesn't have truth information!");

    const auto& truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return false;

    // Count the visible particles
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(pEvent->truth.particles);
    const auto nMu = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 13, useAbsPdg);
    const auto nProton = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 2212, useAbsPdg);
    const auto nPion = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 211, useAbsPdg);
    const auto nOther = visibleParticles.size() - (nMu + nProton + nPion);

    // Insist on the CC1Pi topology
    return (nMu == 1 && nPion == 1 && nOther == 0);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::IsTrueCC0Pi(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg, const float protonMomentumThreshold)
{
    if (!pEvent->metadata.hasTruthInfo())
        throw std::invalid_argument("AnalysisHelper::IsTrueCC0Pi - Input event doesn't have truth information!");

    const auto& truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return false;

    // Count the visible particles
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(truth.particles);
    const auto nMu = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 13, useAbsPdg);
    const auto nProton = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 2212, useAbsPdg);
    const auto nPion = AnalysisHelper::CountParticlesWithPdgCode(visibleParticles, 211, useAbsPdg);
    const auto nProtonAboveMomentumThreshold = AnalysisHelper::CountParticlesAboveMomentumThreshold(visibleParticles, 2212, useAbsPdg, protonMomentumThreshold);
    const auto nOther = visibleParticles.size() - (nMu + nProton + nPion);

    // Insist on the CC0Pi topology
    return (nMu == 1 && nPion == 0 && nOther == 0 && nProtonAboveMomentumThreshold>=1);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::PassesVisibilityThreshold(const Event::Truth::Particle &particle)
{
    // First only consider specific particle types that can lead to hits
    const auto absPDG = std::abs(particle.pdgCode());

    return (absPDG == 11   || // Electron
            absPDG == 13   || // Muon
            absPDG == 2212 || // Proton
            absPDG == 22   || // Photon
            absPDG == 111  || // Pi0 (visible via decay products)
            absPDG == 211  || // Charged pion
            absPDG == 321);   // Charged kaon
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

unsigned int AnalysisHelper::CountParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode, const bool useAbsPdg)
{
    unsigned int count = 0;

    for (const auto &particle : particles)
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();
        if (pdg == pdgCode)
            count++;
    }

    return count;
}

unsigned int AnalysisHelper::CountParticlesAboveMomentumThreshold(const std::vector<Event::Truth::Particle> &particles, const int pdgCode, const bool useAbsPdg, const float momentumThreshold)
{
    unsigned int count = 0;

    for (const auto &particle : particles)
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();
        if (pdg == pdgCode && particle.momentum()>=momentumThreshold)
            count++;
    }

    return count;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::CountGoldenParticlesWithPdgCode(const std::vector<Event::Truth::Particle> &particles, const int pdgCode, const bool useAbsPdg)
{
    std::vector<Event::Truth::Particle> goldenParticles;
    for (const auto &particle : particles)
    {
        if (AnalysisHelper::IsGolden(particle))
            goldenParticles.push_back(particle);
    }

    return AnalysisHelper::CountParticlesWithPdgCode(goldenParticles, pdgCode, useAbsPdg);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::GetPdgCodeCountMap(const std::vector<Event::Truth::Particle> &particles, const bool useAbsPdg, std::vector<int> &foundPdgs, std::unordered_map<int, unsigned int> &pdgCodeCountMap)
{
    if (!foundPdgs.empty())
        throw std::invalid_argument("AnalysisHelper::GetPdgCodeCountMap - input foundPdgs vector isn't empty");

    if (!pdgCodeCountMap.empty())
        throw std::invalid_argument("AnalysisHelper::GetPdgCodeCountMap - input pdgCodeCountMap isn't empty");

    // Get the mapping
    for (const auto &particle : particles)
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();
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

std::string AnalysisHelper::GetTopologyString(const std::vector<Event::Truth::Particle> &particles, const bool useAbsPdg, const bool countProtonsInclusively)
{
    // Count the particles by PDG code
    std::vector<int> foundPdgs;
    std::unordered_map<int, unsigned int> pdgCodeCountMap;
    AnalysisHelper::GetPdgCodeCountMap(particles, useAbsPdg, foundPdgs, pdgCodeCountMap);

    unsigned int nOther = 0;
    unsigned int nProtons = 0;
    std::string topology = "";

    for (const auto &pdg : foundPdgs)
    {
        const auto count = pdgCodeCountMap.at(pdg);
        const auto countStr = std::to_string(count);

        switch (pdg)
        {
            case 11:
                topology += countStr + " e" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case -11:
                topology += countStr + " e" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case 13:
                topology += countStr + " Mu" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case -13:
                topology += countStr + " Mu" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case 22:
                topology += countStr + " Gamma  ";
                break;
            case 211:
                topology += countStr + " Pi" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case -211:
                topology += countStr + " Pi" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case 111:
                topology += countStr + " Pi0  ";
                break;
            case 321:
                topology += countStr + " K" + (useAbsPdg ? "" : "+") + "  ";
                break;
            case -321:
                topology += countStr + " K" + (useAbsPdg ? "" : "-") + "  ";
                break;
            case 2112:
                topology += countStr + " n  ";
                break;
            case 2212:
                nProtons = count;
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

std::string AnalysisHelper::GetClassificationString(const std::shared_ptr<Event> &pEvent, const bool useAbsPdg, const bool countProtonsInclusively, const bool selectCC0pi, const float protonMomentumThreshold)
{
    // Check if data
    if (!pEvent->metadata.hasTruthInfo())
        return "D,  ";

    const auto& truth = pEvent->truth;

    // Insist the true neutrino is fiducial
    if (!AnalysisHelper::IsFiducial(truth.nuVertex()))
        return "NF, ";

    std::string classification = "";
    const auto visibleParticles = AnalysisHelper::SelectVisibleParticles(truth.particles);

    // Signal or background
    const auto isTrueSignal = selectCC0pi ? AnalysisHelper::IsTrueCC0Pi(pEvent, useAbsPdg, protonMomentumThreshold) : AnalysisHelper::IsTrueCC1Pi(pEvent, useAbsPdg);
    classification += isTrueSignal ? "S" : "B,  ";

    if (isTrueSignal)
    {
        // Check if we have a golden pion
        const auto hasGoldenPion = (AnalysisHelper::CountGoldenParticlesWithPdgCode(visibleParticles, 211, useAbsPdg) != 0);

        if (hasGoldenPion)
            classification += " G,";
        else
            classification += ",  ";
    }

    // Add the topology classification
    classification += "  " + AnalysisHelper::GetTopologyString(visibleParticles, useAbsPdg, countProtonsInclusively);

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

float AnalysisHelper::GetFiducialVolume()
{
    // Get the width of the FV in each direction - here we take off the fiducial margins (see AnalysisHelper::IsFiducial)
    const auto widthX = (GeometryHelper::highX - GeometryHelper::lowX) - 10.f - 10.f;
    const auto widthY = (GeometryHelper::highY - GeometryHelper::lowY) - 10.f - 10.f;
    const auto widthZ = (GeometryHelper::highZ - GeometryHelper::lowZ) - 10.f - 50.f;

    return widthX * widthY * widthZ;
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

bool AnalysisHelper::IsInActiveVolume(const TVector3 &point)
{
    return AnalysisHelper::IsPointWithinMargins(point, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------
unsigned int AnalysisHelper::GetBestMatchedTruthParticleIndex(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    if(recoParticle.hasMatchedMCParticle.IsSet()) 
    {
        std::cout<<"DEBUG AnalysisHelper::GetBestMatchedTruthParticleIndex - reco particle has matched MC particle"<<std::endl;
        return GetBestMatchedTruthParticleIndexUbCC1pi(recoParticle, truthParticles, applyVisibilityThreshold);
    }
    else
    {
        return GetBestMatchedTruthParticleIndexPeLEE(recoParticle, truthParticles, applyVisibilityThreshold);
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    if (truthParticles.empty())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - no truth particles supplied but reco particle has a match!");

    std::string errorMessage = "AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - input reco particle has missing or invalid backtracked values: ";

    bool hasError = false;

    if (!recoParticle.momentumXBacktracked.IsSet()) {
        errorMessage += "momentumXBacktracked is not set, ";
        hasError = true;
    }

    if (!recoParticle.momentumYBacktracked.IsSet()) {
        errorMessage += "momentumYBacktracked is not set, ";
        hasError = true;
    }

    if (!recoParticle.momentumZBacktracked.IsSet()) {
        errorMessage += "momentumZBacktracked is not set, ";
        hasError = true;
    }

    if (!recoParticle.pdgBacktracked.IsSet()) {
        errorMessage += "pdgBacktracked is not set, ";
        hasError = true;
    } else if (recoParticle.pdgBacktracked() == -2147483647) {
        errorMessage += "pdgBacktracked is invalid, ";
        hasError = true;
    }

    if (hasError) {
        throw std::logic_error(errorMessage);
    }

    constexpr auto tolerance = 1e-6;


    unsigned int bestMatchedParticleId = std::numeric_limits<unsigned int>::max();
    const auto backtrackedPDG = recoParticle.pdgBacktracked();
    const auto momentumX = recoParticle.momentumXBacktracked();
    const auto momentumY = recoParticle.momentumYBacktracked();
    const auto momentumZ = recoParticle.momentumZBacktracked();
    const auto backtrackedMomentum = std::sqrt(momentumX * momentumX + momentumY * momentumY + momentumZ * momentumZ);

    for (unsigned int i = 0; i < truthParticles.size(); ++i)
    {
        if (applyVisibilityThreshold && !AnalysisHelper::PassesVisibilityThreshold(truthParticles.at(i)))
        {
            // std::cout<<"DEBUG AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - truth particle with pdg "<<truthParticles.at(i).pdgCode()<<" doesn't pass visibility threshold"<<std::endl;
            continue;
        }

        if (backtrackedPDG != truthParticles.at(i).pdgCode()) {
            std::cout << "DEBUG AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - truth particle pdg doesn't match. Backtracked: " << backtrackedPDG << " vs Truth particle: " << truthParticles.at(i).pdgCode() << std::endl;
            continue;
        }

        const auto diff = std::abs(backtrackedMomentum - truthParticles.at(i).momentum());
        if(diff < tolerance)
        {
            if(bestMatchedParticleId != std::numeric_limits<unsigned int>::max())
                throw std::logic_error("AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - input reco particle has more than one matched truth particle");
            bestMatchedParticleId = i;
        }
        else
        {
            // std::cout<<"DEBUG AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - truth particle momentum doesn't match. Backtracked: "<<backtrackedMomentum<<" vs Truth particle: "<<truthParticles.at(i).momentum()<<std::endl;
        }
    }

    if (bestMatchedParticleId == std::numeric_limits<unsigned int>::max())
        throw std::logic_error("AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - input reco particle with pdg: "+std::to_string(backtrackedPDG)+" has no matched truth particle");


    // std::cout<<"DEBUG AnalysisHelper::GetBestMatchedTruthParticleIndexPeLEE - successfully matched reco to truth particle"<<std::endl;
    return bestMatchedParticleId;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetBestMatchedTruthParticleIndexUbCC1pi(const Event::Reco::Particle &recoParticle, const std::vector<Event::Truth::Particle> &truthParticles, const bool applyVisibilityThreshold)
{
    if (!recoParticle.hasMatchedMCParticle.IsSet())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndexUbCC1pi - input reco particle is from an event without truth info");

    if (!recoParticle.hasMatchedMCParticle())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndexUbCC1pi - input reco particle has no matched truth particle");

    if (truthParticles.size() != recoParticle.truthMatchPurities().size())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndexUbCC1pi - wrong number of input truth particles");

    if (truthParticles.empty())
        throw std::invalid_argument("AnalysisHelper::GetBestMatchedTruthParticleIndexUbCC1pi - no truth particles supplied but reco particle has a match!");


    unsigned int bestMatchedParticleId = std::numeric_limits<unsigned int>::max();
    float bestMatchScore = -std::numeric_limits<float>::max();
    bool foundMatch = false;

    for (unsigned int i = 0; i < truthParticles.size(); ++i)
    {
        if (applyVisibilityThreshold && !AnalysisHelper::PassesVisibilityThreshold(truthParticles.at(i)))
            continue;

        const auto score = recoParticle.truthMatchPurities().at(i);
        if (score < bestMatchScore)
            continue;

        bestMatchScore = score;
        bestMatchedParticleId = i;
        foundMatch = true;
    }

    if (!foundMatch)
        throw std::logic_error("AnalysisHelper::GetBestMatchedTruthParticleIndexUbCC1pi - input reco particle has no matched truth particle");

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

std::shared_ptr<TF1> AnalysisHelper::GetRangeToMomentumFunction()
{
    // Work out the maximum straight-line range a particle can have, so we can set sensible limits on the fit
    const auto minRange = 0.f;
    const auto maxRange = std::pow( std::pow(GeometryHelper::highX - GeometryHelper::lowX, 2) +
                                    std::pow(GeometryHelper::highY - GeometryHelper::lowY, 2) +
                                    std::pow(GeometryHelper::highZ - GeometryHelper::lowZ, 2) , 0.5f);

    std::shared_ptr<TF1> pFunc(new TF1(("fitFunc_" + std::to_string(m_fitFunctionIndex++)).c_str(), "[0] + [1]*x - [2]*pow(x, -[3])", minRange, maxRange));
    //std::shared_ptr<TF1> pFunc(new TF1(("fitFunc_" + std::to_string(m_fitFunctionIndex++)).c_str(), "[0] + [1]*x", minRange, maxRange));
    return pFunc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TF1> AnalysisHelper::GetRangeToMomentumFunctionMuon()
{
    RangeToMomentumFitParameters params;
    params.a = 0.19581;
    params.b = 0.002107;
    params.c = 0.15896;
    params.d = 0.19468;

    auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();
    AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

    return pFunc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TF1> AnalysisHelper::GetRangeToMomentumFunctionPion()
{
    RangeToMomentumFitParameters params;
    params.a = 0.25798;
    params.b = 0.0024088;
    params.c = 0.18828;
    params.d = 0.11687;

    auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();
    AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

    return pFunc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

std::shared_ptr<TF1> AnalysisHelper::GetRangeToMomentumFunctionProton()
{
    RangeToMomentumFitParameters params;
    params.a = 14.96; // Fit full range: 14.96, Fit 0-30cm: 5.7364
    params.b = 0.0043489; // Fit full range: 0.0043489, Fit 0-30cm: 0.0074425
    params.c = 14.688; // Fit full range: 14.688, Fit 0-30cm: 5.4472
    params.d = 0.0053518; // Fit full range: 0.0053518, Fit 0-30cm: 0.0097336

    auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();
    AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

    return pFunc;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::SetRangeToMomentumFunctionParameters(const RangeToMomentumFitParameters &params, std::shared_ptr<TF1> &pFunc)
{
    pFunc->SetParameter(0, params.a);
    pFunc->SetParameter(1, params.b);
    pFunc->SetParameter(2, params.c);
    pFunc->SetParameter(3, params.d);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::GetRangeToMomentumFunctionParameters(const std::shared_ptr<TF1> &pFunc, RangeToMomentumFitParameters &params)
{
    params.a = pFunc->GetParameter(0);
    params.b = pFunc->GetParameter(1);
    params.c = pFunc->GetParameter(2);
    params.d = pFunc->GetParameter(3);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::GetRangeToMomentumFitParameters(const std::vector<float> &ranges, const std::vector<float> &momenta, RangeToMomentumFitParameters &params, const float limitLow, const float limitHigh)
{
    if (ranges.size() != momenta.size())
        throw std::invalid_argument("AnalysisHelper::GetRangeToMomentumFitParameters - input ranges and momenta are of different size");

    // In general for a given range, there can be a significant (and asymmetric) spread of momenta around the most probable value.
    // The idea here is to first bin the ranges, and get the most probable value at each range, then fit those most probable values

    // The number of data points to put in each bin
    const unsigned int nPointsPerBin = 100u;

    // The maximum allowed bin width
    const float maxBinWidth = 10.f;

    // The momentum resolution scale to use when finding the MPV
    const float momentumScale = 0.01f;

    // First turn the two input vectors into a vector of pairs so they can be sorted
    std::vector< std::pair<float, float> > dataPoints;
    for (unsigned int i = 0; i < ranges.size(); ++i)
    {
        dataPoints.emplace_back(ranges.at(i), momenta.at(i));
    }

    // Now sort the points by range
    std::sort(dataPoints.begin(), dataPoints.end(), [](const auto &a, const auto &b){ return a.first < b.first; });

    // Now organise the data points into bins, each bin corresponds to a set of nearby ranges
    std::vector< std::vector< std::pair<float, float> > > binnedPoints;
    for (unsigned int i = 0; i < dataPoints.size(); ++i)
    {
        // If we have enough points in the current bin, then add another
        if (i % nPointsPerBin == 0)
            binnedPoints.emplace_back();

        // Add this data point to the last bin
        binnedPoints.back().push_back(dataPoints.at(i));
    }

    // Next find the most probable value in each bin
    std::vector<float> fitRanges, fitMomenta;
    for (const auto &bin : binnedPoints)
    {
        // Ignore the last bin if it doesn't have enough points
        if (bin.size() < nPointsPerBin)
            continue;

        // Extract the ranges and momenta in this bin
        std::vector<float> rangesInBin, momentaInBin;
        for (const auto &point : bin)
        {
            rangesInBin.push_back(point.first);
            momentaInBin.push_back(point.second);
        }

        // Get the spread of ranges in the bin
        const auto minRange = *std::min_element(rangesInBin.begin(), rangesInBin.end());
        const auto maxRange = *std::max_element(rangesInBin.begin(), rangesInBin.end());
        const auto binWidth = maxRange - minRange;

        // Don't use bins that have a large spread of values
        if (binWidth > maxBinWidth)
            continue;

        // Get the mean of the ranges in the bin
        const auto meanRange = std::accumulate(rangesInBin.begin(), rangesInBin.end(), 0.f) / static_cast<float>(rangesInBin.size());

        // Get the minimum momenta in the bin
        const auto minMomentum = *std::min_element(momentaInBin.begin(), momentaInBin.end());

        // Organise the momenta into their own bins
        std::map<unsigned int, std::vector<float> > momentumBins;
        for (const auto &momentum : momentaInBin)
        {
            const unsigned int binIndex = std::floor((momentum - minMomentum) / momentumScale);
            momentumBins[binIndex].push_back(momentum);
        }

        // Sort the momentum bins by the number of entries
        std::vector< std::pair<unsigned int, std::vector<float> > > momentumBinsVector;
        for (const auto &momentumBin : momentumBins)
        {
            momentumBinsVector.emplace_back(momentumBin.first, momentumBin.second);
        }
        std::sort(momentumBinsVector.begin(), momentumBinsVector.end(), [](const auto &a, const auto &b){ return a.second.size() > b.second.size(); });

        // Find the mean momenta in the most probable bin
        const auto mostProbableMomenta = momentumBinsVector.front().second;
        const auto meanMomentum = std::accumulate(mostProbableMomenta.begin(), mostProbableMomenta.end(), 0.f) / static_cast<float>(mostProbableMomenta.size());

        // Save this as a data point to fit
        fitRanges.push_back(meanRange);
        fitMomenta.push_back(meanMomentum);
    }

    // Determine the range over which to fit
    auto minRange = *std::min_element(fitRanges.begin(), fitRanges.end());
    if (limitLow != -std::numeric_limits<float>::max()){
        minRange = limitLow;
    }
    auto maxRange = *std::max_element(fitRanges.begin(), fitRanges.end());
    if (limitHigh != -std::numeric_limits<float>::max()){
        maxRange = limitHigh;
    }


    // Make a TGraph to fit
    const auto nPoints = fitRanges.size();
    std::shared_ptr<TGraph> pGraph(new TGraph(nPoints, fitRanges.data(), fitMomenta.data()));

    // Get the function we want to fit
    auto pFunc = AnalysisHelper::GetRangeToMomentumFunction();

    // Set some initial values to help the fit converge quickly. For reference, I just eyeballed these numbers from a plot of muons
    params.a = 2e-1f;
    params.b = 2e-3f;
    params.c = 1e-1f;
    params.d = 2e-1f;
    AnalysisHelper::SetRangeToMomentumFunctionParameters(params, pFunc);

    // Run the fit
    //  - M = do it gooder (about as much information as the ROOT documentation gives you!)
    //  - N = don't draw the result
    //  - Q = quiet mode
    //  - S = return the fit result
    const auto fitResult = pGraph->Fit(pFunc->GetName(), "MNQS", "", minRange, maxRange);
    if (!fitResult->IsValid())
        std::cout << "AnalysisHelper::GetRangeToMomentumFitParameters - WARNING. Fit was deemed problematic, check the result!" << std::endl;

    // Store the fitted parameters
    AnalysisHelper::GetRangeToMomentumFunctionParameters(pFunc, params);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetPionMomentumFromRange(const float &range)
{
    auto pFunc = AnalysisHelper::GetRangeToMomentumFunctionPion();
    return pFunc->Eval(range);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetMuonMomentumFromRange(const float &range)
{
    auto pFunc = AnalysisHelper::GetRangeToMomentumFunctionMuon();
    return pFunc->Eval(range);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetMuonMomentumFromMCS(const Event::Reco::Particle &muon)
{
    if (muon.mcsMomentumForwardMuon.IsSet())
        return muon.mcsMomentumForwardMuon();

    if (muon.mcsMomentumBackwardMuon.IsSet())
        return muon.mcsMomentumBackwardMuon();

    // ATTN could throw here instead - would need to add event selection cut to remove ents in which muon momentum isn't available
    return -std::numeric_limits<float>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetMuonMomentum(const Event::Reco::Particle &muon)
{
    if (AnalysisHelper::IsContained(muon))
        return AnalysisHelper::GetMuonMomentumFromRange(muon.range());

    return AnalysisHelper::GetMuonMomentumFromMCS(muon);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetProtonMomentumFromRange(const float &range)
{
    auto pFunc = AnalysisHelper::GetRangeToMomentumFunctionProton();
    return pFunc->Eval(range);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetProtonMomentumFromRangeLarsoft(const float &range)
{
    // Reimplemented from larsoft: https://nusoft.fnal.gov/larsoft/doxsvn/html/TrackMomentumCalculator_8cxx_source.html

     /*Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
     website: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html

     CSDA values:
     double KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300,
     350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
     1500, 2000, 2500, 3000, 4000, 5000};

     double Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296,
     2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1,
     9.184E1, 1.138E2, 1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2,
     2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
     7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3};

     Functions below are obtained by fitting power and polynomial fits to
     KE_MeV vs Range (cm) graph. A better fit was obtained by splitting the
     graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
     polynomial of power 6 was a good fit

     Fit errors for future purposes:
     For power function fit: a=0.388873; and b=0.00347075
     Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3
     err=4.45542E-7; p4 err=4.16428E-10; p5 err=1.81679E-13;p6
     err=2.96958E-17;*/

     //*********For proton, the calculations are valid up to 3.022E3 cm range
     //corresponding to a Muon KE of 5 GeV**********//

     if (range < 0 || std::isnan(range)) {
         throw std::logic_error("AnalysisHelper::GetProtonMomentumFromRangeLarsoft - Invalid track range " + std::to_string(range));
         return -std::numeric_limits<float>::max();
     }

     float KE = -999;
     float Momentum = -std::numeric_limits<float>::max();
     float M = 938.272;
     if (range > 0 && range <= 80){
         KE = 29.9317 * std::pow(range, 0.586304);
     }
     else if (range > 80 && range <= 3.022E3){
         KE = 149.904 + (3.34146 * range)
         + (-0.00318856 * range * range)
         + (4.34587E-6 * range * range * range)
         + (-3.18146E-9 * range * range * range * range)
         + (1.17854E-12 * range * range * range * range * range)
         + (-1.71763E-16 * range * range * range * range * range * range);
     }

    if (KE >= 0)
    Momentum = std::sqrt((KE * KE) + (2 * M * KE));

    Momentum = Momentum / 1000;

    return Momentum;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetRecoAnalysisData(const Event::Reco &reco, const std::vector<int> &assignedPdgCodes, const bool passesGoldenPionSelection)
{
    AnalysisData data;

    // Sanity check
    const auto recoParticles = reco.particles;
    if (assignedPdgCodes.size() != recoParticles.size())
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisData - The assigned PDG codes is the wrong size");

    auto muonIndex = std::numeric_limits<unsigned int>::max();
    auto pionIndex = std::numeric_limits<unsigned int>::max();
    unsigned int nMuons = 0u;
    unsigned int nPions = 0u;
    unsigned int nProtons = 0u;
    unsigned int nOther = 0u;

    for (unsigned int index = 0; index < assignedPdgCodes.size(); ++index)
    {
        const auto recoPdg = assignedPdgCodes.at(index);
        switch (recoPdg)
        {
            case 13:
                muonIndex = index;
                nMuons++;
                break;
            case 211:
                pionIndex = index;
                nPions++;
                break;
            case 2212:
            {
                nProtons++;
                break;
            }
            default:
                nOther++;
                break;
        }
    }

    // Make sure the reco PDGs make sense
    if (nMuons != 1)
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisData - Reconstructed " + std::to_string(nMuons) + " muons");

    if (nPions != 1)
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisData - Reconstructed " + std::to_string(nMuons) + " pions");

    // Get the muon and pion reconstructed particles
    const auto muon = recoParticles.at(muonIndex);
    const auto pion = recoParticles.at(pionIndex);

    // Find the reconstructed variables
    const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
    data.muonCosTheta = muonDir.Z();
    data.muonPhi = std::atan2(muonDir.Y(), muonDir.X());
    data.muonMomentum = AnalysisHelper::GetMuonMomentum(muon);

    const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
    data.pionCosTheta = pionDir.Z();
    data.pionPhi = std::atan2(pionDir.Y(), pionDir.X());
    data.pionMomentum = AnalysisHelper::GetPionMomentumFromRange(pion.range());

    data.muonPionAngle = std::acos(muonDir.Dot(pionDir));
    data.nProtons = nProtons;
    data.hasGoldenPion = passesGoldenPionSelection;

    // Fill everything else with dummy values
    data.protonMomentum = -std::numeric_limits<float>::max();
    data.protonCosTheta = -std::numeric_limits<float>::max();
    data.protonPhi = -std::numeric_limits<float>::max();
    data.muonProtonAngle = -std::numeric_limits<float>::max();

    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetRecoAnalysisDataCC0Pi(const Event::Reco &reco, const std::vector<int> &assignedPdgCodes, const bool passesGoldenProtonSelection)
{
    AnalysisData data;


    // Sanity check
    const auto recoParticles = reco.particles;
    if (assignedPdgCodes.size() != recoParticles.size())
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisDataCC0pi - The assigned PDG codes is the wrong size");

    auto muonIndex = std::numeric_limits<unsigned int>::max();
    auto protonIndex = std::numeric_limits<unsigned int>::max();
    auto highestProtonMomentum = -std::numeric_limits<float>::max();
    unsigned int nMuons = 0u;
    unsigned int nPions = 0u;
    unsigned int nProtons = 0u;
    unsigned int nOther = 0u;

    for (unsigned int index = 0; index < assignedPdgCodes.size(); ++index)
    {
        const auto recoPdg = assignedPdgCodes.at(index);
        switch (recoPdg)
        {
            case 13:
                muonIndex = index;
                nMuons++;
                break;
            case 211:
                nPions++;
                break;
            case 2212:
            {
                nProtons++;
                const auto proton = recoParticles.at(index);

                if(proton.range.IsSet())
                {
                    auto range = proton.range();
                    auto protonMomentum = AnalysisHelper::GetPionMomentumFromRange(range); //ATTN: Use pion momentum here for reco values
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    if (protonMomentum>highestProtonMomentum)
                    {
                        protonIndex = index;
                        highestProtonMomentum=protonMomentum;
                    }
                }

                break;
            }
            default:
                nOther++;
                break;
        }
    }

    // Make sure the reco PDGs make sense
    if (nMuons != 1)
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisDataCC0pi - Reconstructed " + std::to_string(nMuons) + " muons for CC0pi");

    if (nPions != 0)
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisDataCC0pi - Reconstructed " + std::to_string(nPions) + " pions for CC0pi");

    if (nProtons == 0)
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisDataCC0pi - Reconstructed " + std::to_string(nProtons) + " protons for CC0pi");

    // Get the muon and pion reconstructed particles
    const auto muon = recoParticles.at(muonIndex);
    const auto proton = recoParticles.at(protonIndex);

    // Find the reconstructed variables
    const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
    data.muonCosTheta = muonDir.Z();
    data.muonPhi = std::atan2(muonDir.Y(), muonDir.X());
    data.muonMomentum = AnalysisHelper::GetMuonMomentum(muon);

    const auto protonDir = TVector3(proton.directionX(), proton.directionY(), proton.directionZ()).Unit();
    data.protonCosTheta = protonDir.Z();
    data.protonPhi = std::atan2(protonDir.Y(), protonDir.X());
    if(!proton.range.IsSet())
        throw std::logic_error("AnalysisHelper::GetRecoAnalysisDataCC0Pi - range not set for proton");
    data.protonMomentum = highestProtonMomentum; //AnalysisHelper::GetProtonMomentumFromRange(proton.range());

    data.muonProtonAngle = std::acos(muonDir.Dot(protonDir));
    data.nProtons = nProtons;
    data.hasGoldenProton = passesGoldenProtonSelection;

    // Fill everything else with dummy values
    data.pionMomentum = -std::numeric_limits<float>::max();
    data.pionCosTheta = -std::numeric_limits<float>::max();
    data.pionPhi = -std::numeric_limits<float>::max();
    data.muonPionAngle = -std::numeric_limits<float>::max();
    data.hasGoldenPion = false;


    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetParticleIndexWithPdg(const std::vector<int> &assignedPdgCodes, const int pdgCode)
{
    std::vector<unsigned int> indices;

    for (unsigned int i = 0; i < assignedPdgCodes.size(); ++i)
    {
        if (assignedPdgCodes.at(i) != pdgCode)
            continue;

        indices.push_back(i);
    }

    if (indices.size() != 1)
        throw std::logic_error("AnalysisHelper::GetParticleIndexWithPdg - Found " + std::to_string(indices.size()) + " particles with assigned PDG code: " + std::to_string(pdgCode) + ", expected 1");

    return indices.front();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetTruthAnalysisData(const Event::Truth &truth, const bool useAbsPdg, const float protonMomentumThreshold)
{
    AnalysisData data;

    bool foundPion = false;
    bool foundMuon = false;

    data.nProtons = 0u;
    TVector3 pionDir, muonDir;


    for (const auto &particle : AnalysisHelper::SelectVisibleParticles(truth.particles))
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();

        const auto dir = TVector3(particle.momentumX(), particle.momentumY(), particle.momentumZ()).Unit();
        const auto cosTheta = dir.Z();
        const auto phi = std::atan2(dir.Y(), dir.X());


        if (pdg == 211)
        {
            if (foundPion)
                throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Found multiple pions! Are you sure this is a signal event?");

            foundPion = true;

            data.pionMomentum = particle.momentum();
            data.pionCosTheta = cosTheta;
            data.pionPhi = phi;
            data.hasGoldenPion = AnalysisHelper::IsGolden(particle);

            pionDir = dir;

            continue;
        }

        if (pdg == 13)
        {
            if (foundMuon)
                throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Found multiple muons! Are you sure this is a signal event?");

            foundMuon = true;

            data.muonMomentum = particle.momentum();
            data.muonCosTheta = cosTheta;
            data.muonPhi = phi;

            muonDir = dir;

            continue;
        }

        if (pdg == 2212)
        {
            const auto passesMomentumThreshold = (particle.momentum() > protonMomentumThreshold);
            if (!passesMomentumThreshold)
                continue;

            data.nProtons++;
            continue;
        }

        throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Found particle with unexpected PDG code = " + std::to_string(particle.pdgCode()));
    }

    if (!foundPion || !foundMuon)
        throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Input event doesn't contain both a muon and a pion!");


    data.muonPionAngle = std::acos(muonDir.Dot(pionDir));

    // Fill everything else with dummy values
    data.protonMomentum = -std::numeric_limits<float>::max();
    data.protonCosTheta = -std::numeric_limits<float>::max();
    data.protonPhi = -std::numeric_limits<float>::max();
    data.muonProtonAngle = -std::numeric_limits<float>::max();
    data.hasGoldenProton = false; // leading proton can be golden but not relevant for regular analysis (only sideband)

    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetTruthAnalysisDataCC0Pi(const Event::Truth &truth/*, const Event::Reco &reco, const std::vector<int> &assignedPdgCodes*/, const bool useAbsPdg, const float protonMomentumThreshold)
{
    AnalysisData data;
    bool foundMuon = false;
    bool foundProton = false;
    auto highestProtonMomentum = -std::numeric_limits<float>::max();

    data.nProtons = 0u;
    TVector3 protonDir, muonDir;


    // // Get the index of the proton in the truth particles that matches the pion in the reco particles
    // if (assignedPdgCodes.size() != reco.particles.size())
    //     throw std::logic_error("AnalysisHelper::GetRecoAnalysisDataCC0pi - The assigned PDG codes is the wrong size");
    // const auto recoPionCandidate = reco.particles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodes, 211));
    // protonAsPionTruthIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoPionCandidate, truth.particles);

    for (const auto &particle : AnalysisHelper::SelectVisibleParticles(truth.particles))
    {
        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();

        const auto dir = TVector3(particle.momentumX(), particle.momentumY(), particle.momentumZ()).Unit();
        const auto cosTheta = dir.Z();
        const auto phi = std::atan2(dir.Y(), dir.X());

        if (pdg == 13)
        {
            if (foundMuon)
                throw std::logic_error("AnalysisHelper::GetTruthAnalysisDataCC0pi - Found multiple muons! Are you sure this is a CC0pi signal event?");

            foundMuon = true;

            data.muonMomentum = particle.momentum();
            data.muonCosTheta = cosTheta;
            data.muonPhi = phi;

            muonDir = dir;

            continue;
        }
        if (pdg == 2212)
        {
            const auto passesMomentumThreshold = (particle.momentum() > protonMomentumThreshold);
            if (!passesMomentumThreshold)
                continue;
            data.nProtons++;
            const auto isHighestEnergyProton = (particle.momentum() > highestProtonMomentum);
            if(!isHighestEnergyProton) // Use highest energy proton
                continue;
            foundProton = true;
            highestProtonMomentum = particle.momentum();
            data.protonMomentum = particle.momentum();
            data.protonCosTheta = cosTheta;
            data.protonPhi = phi;

            protonDir = dir;
            data.hasGoldenProton = AnalysisHelper::IsGolden(particle);

            continue;
        }

        throw std::logic_error("AnalysisHelper::GetTruthAnalysisDataCC0pi - Found particle with unexpected PDG code for CC0pi event = " + std::to_string(particle.pdgCode()));
    }

    if (!foundMuon | !foundProton)
        throw std::logic_error("AnalysisHelper::GetTruthAnalysisDataCC0pi - Input event doesn't contain a muon event!");

    data.muonProtonAngle = std::acos(muonDir.Dot(protonDir));

    // Fill everything else with dummy values
    data.pionMomentum = -std::numeric_limits<float>::max();
    data.pionCosTheta = -std::numeric_limits<float>::max();
    data.pionPhi = -std::numeric_limits<float>::max();
    data.muonPionAngle = -std::numeric_limits<float>::max();
    data.hasGoldenPion = false;

    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

AnalysisHelper::AnalysisData AnalysisHelper::GetDummyAnalysisData()
{
    AnalysisData data;

    data.muonMomentum = -std::numeric_limits<float>::max();
    data.muonCosTheta = -std::numeric_limits<float>::max();
    data.muonPhi = -std::numeric_limits<float>::max();
    data.pionMomentum = -std::numeric_limits<float>::max();
    data.pionCosTheta = -std::numeric_limits<float>::max();
    data.pionPhi = -std::numeric_limits<float>::max();
    data.muonPionAngle = -std::numeric_limits<float>::max();
    data.nProtons = std::numeric_limits<unsigned int>::max();
    data.protonMomentum = -std::numeric_limits<float>::max();
    data.protonCosTheta = -std::numeric_limits<float>::max();
    data.protonPhi = -std::numeric_limits<float>::max();
    data.muonProtonAngle = -std::numeric_limits<float>::max();
    data.hasGoldenPion = false;
    data.hasGoldenProton = false;

    return data;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool AnalysisHelper::GetLogLikelihoodRatio(const Member<float> &numerator, const Member<float> &denominator, float &ratio)
{
    ratio = -std::numeric_limits<float>::max();

    if (!numerator.IsSet() || !denominator.IsSet())
        return false;

    if (numerator() < 0.f)
        return false;

    if (denominator() <= std::numeric_limits<float>::epsilon())
        return false;

    ratio = std::log(numerator() / denominator());
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

float AnalysisHelper::GetCountUncertainty(const float &count)
{
    if (count < 0)
        throw std::logic_error("AnalysisHelper::GetCountUncertainty - count is < 0");

    return std::pow(count + 1, 0.5f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

float AnalysisHelper::GetEfficiencyUncertainty(const float &numerator, const float &denominator)
{
    if (denominator <= std::numeric_limits<float>::epsilon())
        throw std::logic_error("AnalysisHelper::GetEfficiencyUncertainty - denominator is <= 0");

    if (numerator > denominator)
        throw std::logic_error("AnalysisHelper::GetEfficiencyUncertainty - numerator > denomintaor");

    const auto n1 = numerator + 1;
    const auto n2 = numerator + 2;
    const auto d2 = denominator + 2;
    const auto d3 = denominator + 3;

    return std::pow((n1*n2)/(d2*d3) - (n1*n1)/(d2*d2), 0.5f);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::PrintLoadingBar(const unsigned int numerator, const unsigned int denominator, const bool fancy)
{
    if (denominator == 0)
        return;

    // Get the current stage of the loading bar
    const auto width = 85u;
    const auto loadedWidth = static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator) / static_cast<float>(denominator))));
    const auto unloadedWidth = static_cast<unsigned int>(width - loadedWidth);

    // Get the fraction of the way through the current stage
    const auto stageLength = static_cast<float>(denominator) / static_cast<float>(width);
    // const unsigned int tenth = std::floor(10 * ((static_cast<float>(numerator) - stageLength * std::floor(static_cast<float>(numerator) / stageLength)) / stageLength));
    const unsigned int tenth = std::floor(10 * std::fmod(static_cast<float>(numerator) / stageLength, 1.f));

    // Now do the same but for the last value of the numerator
    const auto lastLoadedWidth = static_cast<unsigned int>(std::floor(width * (static_cast<float>(numerator - 1) / static_cast<float>(denominator))));
    // const unsigned int lastTenth = std::floor(10 * ((static_cast<float>(numerator - 1) - stageLength * std::floor(static_cast<float>(numerator - 1) / stageLength)) / stageLength));
    const unsigned int lastTenth = std::floor(10 * std::fmod(static_cast<float>(numerator-1) / stageLength, 1.f));

    // Work out if it's worth printing
    const bool shouldPrint = (numerator == 0) || (lastLoadedWidth != loadedWidth) || (lastTenth != tenth);

    if (!shouldPrint)
        return;

    // Jump back 2 lines
    if (numerator != 0)
        std::cout << "\033[F\033[F";

    // Print the fraction and loading bar
    std::cout << numerator << " / " << denominator << std::endl;
    // std::cout << "|" << std::string(loadedWidth, '=') << tenth << std::string(unloadedWidth, ' ') << "|" << std::endl;
    std::cout << "|";
    if(fancy)
    {
        std::cout<<"\033[0;42m" << std::string(loadedWidth, '=') << "\033[0m" << "\033[1;" << 3 << tenth/2+1 << "m" << tenth << "\033[0m";
    }
    else
    {
        std::cout << std::string(loadedWidth, '=') << tenth;
    }
    std::cout << std::string(unloadedWidth, ' ') << "|" << std::endl;;

}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetTrueLeadingProtonIndex(const Event::Truth &truth, const bool useAbsPdg, const float protonMomentumThreshold)
{
    float leadingProtonMom = -std::numeric_limits<float>::max();
    unsigned int leadingProtonIndex = std::numeric_limits<unsigned int>::max();

    for (unsigned int i = 0; i < truth.particles.size(); ++i){
        auto particle = truth.particles.at(i);

        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();

        if (pdg == 2212){
            const auto passesMomentumThreshold = (particle.momentum() > protonMomentumThreshold);
            if (!passesMomentumThreshold)
                continue;

            const auto momentum = particle.momentum();
            if (momentum > leadingProtonMom){
                leadingProtonMom = momentum;
                leadingProtonIndex = i;
            }
        }
    }
    return leadingProtonIndex;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

unsigned int AnalysisHelper::GetTrueMuonIndex(const Event::Truth &truth, const bool useAbsPdg)
{
    bool foundMuon = false;
    unsigned int muonidx = std::numeric_limits<unsigned int>::max();

    for (unsigned int i = 0; i < truth.particles.size(); ++i){
        auto particle = truth.particles.at(i);

        const auto pdg = useAbsPdg ? std::abs(particle.pdgCode()) : particle.pdgCode();

        if (pdg == 13)
        {
            if (foundMuon)
                throw std::logic_error("AnalysisHelper::GetTruthAnalysisData - Found multiple muons! Are you sure this is a signal event?");

            foundMuon = true;

            muonidx = i;

            continue;
        }
    }

    return muonidx;
}

} // namespace ubcc1pi
