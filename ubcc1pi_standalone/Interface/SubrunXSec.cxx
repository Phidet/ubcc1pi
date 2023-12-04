/**
 *  @file  ubcc1pi_standalone/Interface/Subrun.cxx
 *
 *  @brief The implementation of the subrun class
 */

#include "ubcc1pi_standalone/Interface/SubrunXSec.h"

namespace ubcc1pi
{

SubrunXSec::SubrunXSec(const bool hasTruthInfo) : hasTruthWeights(hasTruthInfo) {}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunXSec::Print() const
{
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "SUBRUN" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    // auto &subrun = *this;
    // XSEC_MACRO_SUBRUN_MEMBERS("", subrun, PELEE_MACRO_PRINT_MEMBER)
    // if(hasTruthWeights){XSEC_MACRO_SUBRUN_OPTIONAL_MEMBERS("", subrun, XSEC_MACRO_PRINT_MEMBER)}
}

// -----------------------------------------------------------------------------------------------------------------------------------------

bool SubrunXSec::HasTruthWeights() const
{
    return hasTruthWeights;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunXSec::BindToOutputTree(TTree * pTree)
{
    std::cout << "SubrunXSec::BindToOutputTree" << std::endl;
    // auto &subrun = *this;
    // XSEC_MACRO_SUBRUN_MEMBERS(subrun, subrun, PELEE_MACRO_BIND_OUTPUT_BRANCH)
    // if(hasTruthWeights){XSEC_MACRO_SUBRUN_OPTIONAL_MEMBERS(subrun, subrun, XSEC_MACRO_BIND_OUTPUT_BRANCH)}
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunXSec::BindToInputTree(TTree * pTree)
{
    std::cout << "SubrunXSec::BindToInputTree" << std::endl;
    // auto &subrun = *this;
    // XSEC_MACRO_SUBRUN_MEMBERS(subrun, subrun, PELEE_MACRO_BIND_INPUT_BRANCH)
    // if(hasTruthWeights){XSEC_MACRO_SUBRUN_OPTIONAL_MEMBERS(subrun, subrun, XSEC_MACRO_BIND_INPUT_BRANCH)}
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void SubrunXSec::Reset()
{
    std::cout << "SubrunXSec::Reset" << std::endl;
    // auto &subrun = *this;
    // XSEC_MACRO_SUBRUN_MEMBERS("", subrun, PELEE_MACRO_RESET_MEMBER)
    // if(hasTruthWeights){XSEC_MACRO_SUBRUN_OPTIONAL_MEMBERS("", subrun, XSEC_MACRO_RESET_MEMBER)}
}

} // namespace ubcc1pi
