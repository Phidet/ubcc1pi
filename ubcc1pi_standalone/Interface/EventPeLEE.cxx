/**
 *  @file  ubcc1pi_standalone/Interface/Event.cxx
 *
 *  @brief The implementation of the event class
 */

#include "ubcc1pi_standalone/Interface/EventPeLEE.h"
#include <ctime> // DEBUG

namespace ubcc1pi
{

EventPeLEE::EventPeLEE(const bool hasTruthInfo) : hasTruthWeights(hasTruthInfo)
{
    PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", PELEE_MACRO_INIT_MEMBER_VECTOR)
    PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", PELEE_MACRO_INIT_MEMBER_VECTOR)
    PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", PELEE_MACRO_INIT_MEMBER_VECTOR)
}
// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPeLEE::Print() const
{
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "METADATA" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    PELEE_MACRO_EVENT_METADATA_MEMBERS("", metadata, PELEE_MACRO_PRINT_MEMBER)

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "TRUTH" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    PELEE_MACRO_EVENT_TRUTH_MEMBERS("", truth, PELEE_MACRO_PRINT_MEMBER)
    if(hasTruthWeights) {PELEE_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth, PELEE_MACRO_PRINT_MEMBER)}

    unsigned int i = 0;
    for (const auto &particle : truth.particles)
    {
        std::cout << "TRUTH PARTICLE "<< i << ": " << std::endl;
        PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", particle, PELEE_MACRO_PRINT_MEMBER)
        i++;
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "RECO" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    PELEE_MACRO_EVENT_RECO_MEMBERS("", reco, PELEE_MACRO_PRINT_MEMBER)

    unsigned int j = 0;
    for (const auto &particle : reco.particles)
    {
        std::cout << "RECO PARTICLE "<< j << ": " << std::endl;
        PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS("", particle, PELEE_MACRO_PRINT_MEMBER)
        j++;
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool EventPeLEE::HasTruthWeights() const
{
    return hasTruthWeights;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPeLEE::BindToOutputTree(TTree * pTree)
{
    PELEE_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    PELEE_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    if(hasTruthWeights) {PELEE_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth, PELEE_MACRO_BIND_OUTPUT_BRANCH)}
    PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", PELEE_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
    PELEE_MACRO_EVENT_RECO_MEMBERS(reco, reco, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", PELEE_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
    PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", PELEE_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPeLEE::BindToInputTree(TTree * pTree)
{
    PELEE_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, PELEE_MACRO_BIND_INPUT_BRANCH)
    PELEE_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, PELEE_MACRO_BIND_INPUT_BRANCH)
    if(hasTruthWeights) {PELEE_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth,  PELEE_MACRO_BIND_INPUT_BRANCH)}
    PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", PELEE_MACRO_BIND_INPUT_VECTOR_BRANCH)
    PELEE_MACRO_EVENT_RECO_MEMBERS(reco, reco, PELEE_MACRO_BIND_INPUT_BRANCH)
    PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", PELEE_MACRO_BIND_INPUT_VECTOR_BRANCH)
    // PELEE_MACRO_EVENT_BACKTRACKED_MEMBERS(backtracked, backtracked, PELEE_MACRO_BIND_INPUT_BRANCH)
    PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", PELEE_MACRO_BIND_INPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPeLEE::Reset()
{
    PELEE_MACRO_EVENT_METADATA_MEMBERS("", metadata, PELEE_MACRO_RESET_MEMBER)

    PELEE_MACRO_EVENT_TRUTH_MEMBERS("", truth, PELEE_MACRO_RESET_MEMBER)
    if(hasTruthWeights) {PELEE_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth, PELEE_MACRO_RESET_MEMBER)}
    PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", PELEE_MACRO_RESET_MEMBER_VECTOR)
    truth.particles.clear();

    PELEE_MACRO_EVENT_RECO_MEMBERS("", reco, PELEE_MACRO_RESET_MEMBER)
    PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", PELEE_MACRO_RESET_MEMBER_VECTOR)
    reco.particles.clear();

    PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", PELEE_MACRO_RESET_MEMBER_VECTOR)
    backtracked.particles.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPeLEE::PrepareForTreeFill()
{
    unsigned int i = 0;
    for (const auto &particle : truth.particles)
    {
            std::cout<<"pelee truth particle "<<i<<std::endl; i++;
            PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, particle, PELEE_MACRO_FILL_MEMBER_VECTOR)
    }

    i = 0;
    for (const auto &particle : reco.particles)
    {
            std::cout<<"pelee reco particle "<<i<<std::endl; i++;
            PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, particle, PELEE_MACRO_FILL_MEMBER_VECTOR)
    }

    i = 0;
    for (const auto &particle : backtracked.particles)
    {
            std::cout<<"pelee backtracked particle "<<i<<std::endl; i++;
            PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, particle, PELEE_MACRO_FILL_MEMBER_VECTOR)
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventPeLEE::PrepareAfterTreeRead()
{
    unsigned int nTruthParticles;
    PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, &nTruthParticles, PELEE_MACRO_GET_MEMBER_VECTOR_SIZE)

    truth.particles.resize(nTruthParticles);
    PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, truth.particles, PELEE_MACRO_READ_MEMBER_VECTOR)

    unsigned int nRecoParticles;
    PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, &nRecoParticles, PELEE_MACRO_GET_MEMBER_VECTOR_SIZE)

    reco.particles.resize(nRecoParticles);
    PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, reco.particles, PELEE_MACRO_READ_MEMBER_VECTOR)

    unsigned int nBacktrackedParticles;
    PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, &nBacktrackedParticles, PELEE_MACRO_GET_MEMBER_VECTOR_SIZE)

    backtracked.particles.resize(nBacktrackedParticles);
    PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, backtracked.particles, PELEE_MACRO_READ_MEMBER_VECTOR)
}

} // namespace ubcc1pi
