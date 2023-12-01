/**
 *  @file  ubcc1pi_standalone/Interface/SubrunMembers.h
 *
 *  @brief The file defining the members of the subrun class
 */


#ifndef UBCC1PI_STANDALONE_INTERFACE_SUBRUNMEMBERS
#define UBCC1PI_STANDALONE_INTERFACE_SUBRUNMEMBERS

// In the definition of the members below the first two parameters (p, q) are used by later macros, the third parameter is a boolean which
// should be true if ROOT needs the address of a pointer when setting branch addresses (for non-standard types). The fourth is the type of
// the varaiable, and the fifth is the variable name itself

/** The subrun members */
#define UBCC1PI_MACRO_SUBRUN_MEMBERS(p, q, f)                                                                                              \
    f(p, q, false, int,         run)                                                                                                       \
    f(p, q, false, int,         subRun)                                                                                                    \
    f(p, q, false, bool,        hasTruthInfo)                                                                                              \
    f(p, q, false, float,       totalPOT)                                                                                                  \
    f(p, q, false, float,       totalGoodPOT)                                                                                              \
    f(p, q, false, int,         totalSpills)                                                                                               \
    f(p, q, false, int,         goodSpills)

/** The PeLEE subrun members */
#define PELEE_MACRO_SUBRUN_MEMBERS(p, q, f)                                                                                                \
    f(p, q, false, int,         run)                                                                                                       \
    f(p, q, false, int,         subRun)                                                                                                    \

#define PELEE_MACRO_SUBRUN_OPTIONAL_MEMBERS(p, q, f)                                                                                                    \
    f(p, q, false, float,       pot)

// /** The PeLEE subrun members */
// #define XSEC_MACRO_SUBRUN_MEMBERS(p, q, f)
//     f(p, q, false, int,         run)

#endif
