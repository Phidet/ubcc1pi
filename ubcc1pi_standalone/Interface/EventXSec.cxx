/**
 *  @file  ubcc1pi_standalone/Interface/Event.cxx
 *
 *  @brief The implementation of the event class
 */

#include "ubcc1pi_standalone/Interface/EventXSec.h"
#include <ctime> // DEBUG

namespace ubcc1pi
{

EventXSec::EventXSec(const bool hasTruthInfo) : hasTruthWeights(hasTruthInfo)
{
    // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", XSEC_MACRO_INIT_MEMBER_VECTOR)
    // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", XSEC_MACRO_INIT_MEMBER_VECTOR)
    // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", XSEC_MACRO_INIT_MEMBER_VECTOR)
}
// -----------------------------------------------------------------------------------------------------------------------------------------

void EventXSec::Print() const
{
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "METADATA" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    // XSEC_MACRO_EVENT_METADATA_MEMBERS("", metadata, PELEE_MACRO_PRINT_MEMBER)

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "TRUTH" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    XSEC_MACRO_EVENT_TRUTH_MEMBERS("", truth, PELEE_MACRO_PRINT_MEMBER)
    // if(hasTruthWeights) {XSEC_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth, XSEC_MACRO_PRINT_MEMBER)}

    // unsigned int i = 0;
    // for (const auto &particle : truth.particles)
    // {
    //     std::cout << "TRUTH PARTICLE "<< i << ": " << std::endl;
    //     // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", particle, XSEC_MACRO_PRINT_MEMBER)
    //     i++;
    // }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "RECO" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    XSEC_MACRO_EVENT_RECO_MEMBERS("", reco, PELEE_MACRO_PRINT_MEMBER)

    // unsigned int j = 0;
    // for (const auto &particle : reco.particles)
    // {
    //     std::cout << "RECO PARTICLE "<< j << ": " << std::endl;
    //     // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS("", particle, XSEC_MACRO_PRINT_MEMBER)
    //     j++;
    // }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool EventXSec::HasTruthWeights() const
{
    return hasTruthWeights;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventXSec::BindToOutputTree(TTree * pTree)
{
    // XSEC_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    XSEC_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    // if(hasTruthWeights) {XSEC_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth, XSEC_MACRO_BIND_OUTPUT_BRANCH)}
    // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", XSEC_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
    XSEC_MACRO_EVENT_RECO_MEMBERS(reco, reco, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", XSEC_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
    // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", XSEC_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventXSec::BindToInputTree(TTree * pTree)
{
    // XSEC_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, PELEE_MACRO_BIND_INPUT_BRANCH)
    std::cout << "DEBUG EventXSec::BindToInputTree() Point 0" << std::endl;
    XSEC_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, PELEE_MACRO_BIND_INPUT_BRANCH)
    std::cout << "DEBUG EventXSec::BindToInputTree() Point 1" << std::endl;
    // if(hasTruthWeights) {XSEC_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth,  XSEC_MACRO_BIND_INPUT_BRANCH)}
    // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", XSEC_MACRO_BIND_INPUT_VECTOR_BRANCH)
    std::cout << "DEBUG EventXSec::BindToInputTree() Point 2" << std::endl;
    XSEC_MACRO_EVENT_RECO_MEMBERS(reco, reco, PELEE_MACRO_BIND_INPUT_BRANCH)
    std::cout << "DEBUG EventXSec::BindToInputTree() Point 3" << std::endl;
    // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", XSEC_MACRO_BIND_INPUT_VECTOR_BRANCH)
    // XSEC_MACRO_EVENT_BACKTRACKED_MEMBERS(backtracked, backtracked, XSEC_MACRO_BIND_INPUT_BRANCH)
    // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", XSEC_MACRO_BIND_INPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventXSec::Reset()
{
    // XSEC_MACRO_EVENT_METADATA_MEMBERS("", metadata, PELEE_MACRO_RESET_MEMBER)

    XSEC_MACRO_EVENT_TRUTH_MEMBERS("", truth, PELEE_MACRO_RESET_MEMBER)
    // if(hasTruthWeights) {XSEC_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(truth, truth, XSEC_MACRO_RESET_MEMBER)}
    // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", XSEC_MACRO_RESET_MEMBER_VECTOR)
    // truth.particles.clear();

    XSEC_MACRO_EVENT_RECO_MEMBERS("", reco, PELEE_MACRO_RESET_MEMBER)
    // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", XSEC_MACRO_RESET_MEMBER_VECTOR)
    // reco.particles.clear();

    // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", XSEC_MACRO_RESET_MEMBER_VECTOR)
    // backtracked.particles.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventXSec::PrepareForTreeFill()
{
    std::cout << "EventXSec::PrepareForTreeFill() called" << std::endl;
    // unsigned int i = 0;
    // for (const auto &particle : truth.particles)
    // {
    //         std::cout<<"XSEC truth particle "<<i<<std::endl; i++;
    //         // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, particle, XSEC_MACRO_FILL_MEMBER_VECTOR)
    // }

    // i = 0;
    // for (const auto &particle : reco.particles)
    // {
    //         std::cout<<"XSEC reco particle "<<i<<std::endl; i++;
    //         // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, particle, XSEC_MACRO_FILL_MEMBER_VECTOR)
    // }

    // i = 0;
    // for (const auto &particle : backtracked.particles)
    // {
    //         std::cout<<"XSEC backtracked particle "<<i<<std::endl; i++;
    //         // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, particle, XSEC_MACRO_FILL_MEMBER_VECTOR)
    // }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void EventXSec::PrepareAfterTreeRead()
{
    // unsigned int nTruthParticles;
    // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, &nTruthParticles, XSEC_MACRO_GET_MEMBER_VECTOR_SIZE)

    // truth.particles.resize(nTruthParticles);
    // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, truth.particles, XSEC_MACRO_READ_MEMBER_VECTOR)

    // unsigned int nRecoParticles;
    // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, &nRecoParticles, XSEC_MACRO_GET_MEMBER_VECTOR_SIZE)

    // reco.particles.resize(nRecoParticles);
    // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, reco.particles, XSEC_MACRO_READ_MEMBER_VECTOR)

    // unsigned int nBacktrackedParticles;
    // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, &nBacktrackedParticles, XSEC_MACRO_GET_MEMBER_VECTOR_SIZE)

    // backtracked.particles.resize(nBacktrackedParticles);
    // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, backtracked.particles, XSEC_MACRO_READ_MEMBER_VECTOR)
}

} // namespace ubcc1pi
