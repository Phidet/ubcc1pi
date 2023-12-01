/**
 *  @file  ubcc1pi_standalone/Interface/EventMembers.h
 *
 *  @brief The file defining the members of the event class
 *
 *  In this file we define all of the member variable that will go into the analysis files. The types and names of the variables are defined
 *  using pre-processor directives.
 *
 *  By defining the variables at the pre-processer level, all of the bookwork can be handled automatically. Each variable need only be
 *  defined in one place and will automatically:
 *      - Be added as a member variable of the ubcc1pi::Event class
 *      - Be bound to a branch if an analysis file is to be written or read using a consistent branch naming convention
 *      - Be assigned a default value and be checked such that you can't accidentally use a default variable at the analysis level
 *      - Show up as a compiler error if used in C++ but not defined here
 *
 *  The Event class has three categories
 *      - Metadata
 *          - Defined by the UBCC1PI_MACRO_EVENT_METADATA_MEMBERS macro
 *          - Accessible in C++ via: event.metadata.variable()
 *          - Accessible in the root tree via: metadata_variable
 *      - Truth
 *          - Defined by the UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS macro
 *          - Accessible in C++ via: event.truth.variable()
 *          - Accessible in the root tree via: truth_variable
 *          - Particle level information is defined by the UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS macro
 *              - Accessible in C++ via: event.truth.particles.at(i).variable()
 *              - Accessible in the root tree via: truth_particle_variable_vect
 *      - Reco
 *          - Defined by the UBCC1PI_MACRO_EVENT_RECO_MEMBERS macro
 *          - Accessible in C++ via: event.reco.variable()
 *          - Accessible in the root tree via: reco_variable
 *          - Particle level information is defined by the UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS macro
 *              - Accessible in C++ via: event.reco.particles.at(i).variable()
 *              - Accessible in the root tree via: reco_particle_variable_vect
 */


#ifndef UBCC1PI_STANDALONE_INTERFACE_EVENTMEMBERS
#define UBCC1PI_STANDALONE_INTERFACE_EVENTMEMBERS

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <TVector3.h>

// In the definition of the members below the first two parameters (p, q) are used by later macros, the third parameter is a boolean which
// should be true if ROOT needs the address of a pointer when setting branch addresses (for non-standard types). The fourth is the type of
// the varaiable, and the fifth is the variable name itself

/** The event metadata members */
#define UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(p, q, f)                                                                                      \
    f(p, q, false, int,         run)                                                                                                       \
    f(p, q, false, int,         subRun)                                                                                                    \
    f(p, q, false, int,         event)                                                                                                     \
    f(p, q, false, bool,        hasTruthInfo)

/** The event truth information members */
#define UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(p, q, f)                                                                                         \
    f(p, q, false, float,                      splineEventWeight)                                                                          \
    f(p, q, false, float,                      genieTuneEventWeight)                                                                       \
    f(p, q, true,  std::vector<std::string>,   systParamNames)                                                                             \
    f(p, q, true,  std::vector<int>,           systParamFirstValueIndex)                                                                   \
    f(p, q, true,  std::vector<float>,         systParamValues)                                                                            \
    f(p, q, false, bool,                       isCC)                                                                                       \
    f(p, q, false, int,                        interactionMode)                                                                            \
    f(p, q, true,  std::string,                interactionString)                                                                          \
    f(p, q, false, int,                        nuPdgCode)                                                                                  \
    f(p, q, false, float,                      nuEnergy)                                                                                   \
    f(p, q, true,  TVector3,                   nuVertex)                                                                                   \
    f(p, q, false, int,                        nFinalStates)                                                                               \
    f(p, q, true,  std::vector<float>,         slicePurities)                                                                              \
    f(p, q, true,  std::vector<float>,         sliceCompletenesses)


/** The event truth particle information members */
#define UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(p, q, f)                                                                                \
    f(p, q, false, int,                 pdgCode)                                                                                           \
    f(p, q, false, float,               startX)                                                                                            \
    f(p, q, false, float,               startY)                                                                                            \
    f(p, q, false, float,               startZ)                                                                                            \
    f(p, q, false, float,               endX)                                                                                              \
    f(p, q, false, float,               endY)                                                                                              \
    f(p, q, false, float,               endZ)                                                                                              \
    f(p, q, false, float,               momentumX)                                                                                         \
    f(p, q, false, float,               momentumY)                                                                                         \
    f(p, q, false, float,               momentumZ)                                                                                         \
    f(p, q, false, float,               momentum)                                                                                          \
    f(p, q, false, float,               energy)                                                                                            \
    f(p, q, false, float,               mass)                                                                                              \
    f(p, q, false, float,               hitWeightU)                                                                                        \
    f(p, q, false, float,               hitWeightV)                                                                                        \
    f(p, q, false, float,               hitWeightW)                                                                                        \
    f(p, q, false, int,                 nElasticScatters)                                                                                  \
    f(p, q, false, int,                 nInelasticScatters)                                                                                \
    f(p, q, false, std::vector<float>,  scatterCosThetas)                                                                                  \
    f(p, q, false, std::vector<float>,  scatterMomentumFracsLost)                                                                          \
    f(p, q, false, std::vector<int>,    scatterIsElastic)                                                                                  \
    f(p, q, false, float,               scatteredMomentum)                                                                                 \
    f(p, q, false, float,               endMomentum)                                                                                       \
    f(p, q, false, bool,                isStopping)                                                                                        \
    f(p, q, false, int,                 endState)                                                                                          \
    f(p, q, false, float,               endStateProductsHitWeightU)                                                                        \
    f(p, q, false, float,               endStateProductsHitWeightV)                                                                        \
    f(p, q, false, float,               endStateProductsHitWeightW)

/** The event reco information members */
#define UBCC1PI_MACRO_EVENT_RECO_MEMBERS(p, q, f)                                                                                          \
    f(p, q, false, bool,                passesCCInclusive)                                                                                 \
    f(p, q, false, int,                 nSlices)                                                                                           \
    f(p, q, false, bool,                hasSelectedSlice)                                                                                  \
    f(p, q, false, float,               selectedTopologicalScore)                                                                          \
    f(p, q, true,  std::vector<float>,  sliceTopologicalScores)                                                                            \
    f(p, q, true,  std::vector<bool>,   sliceIsSelectedAsNu)                                                                               \
    f(p, q, false, bool,                hasNeutrino)                                                                                       \
    f(p, q, false, int,                 nuPdgCode)                                                                                         \
    f(p, q, true,  TVector3,            nuVertex)                                                                                          \
    f(p, q, false, int,                 nFinalStates)                                                                                      \
    f(p, q, false, float,               flashChi2)                                                                                         \
    f(p, q, false, float,               flashTime)                                                                                         \
    f(p, q, false, float,               largestFlashPE)                                                                                    \
    f(p, q, false, float,               largestFlashTime)                                                                                  \
    f(p, q, false, float,               largestFlashTimeWidth)                                                                             \
    f(p, q, false, float,               largestFlashYCtr)                                                                                  \
    f(p, q, false, float,               largestFlashYWidth)                                                                                \
    f(p, q, false, float,               largestFlashZCtr)                                                                                  \
    f(p, q, false, float,               largestFlashZWidth)                                                                                \
    f(p, q, true,  TVector3,            nuVertexNoSCC) /* New PeLEE variable */

/** The event reco particle information members */
#define UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(p, q, f)                                                                                 \
    f(p, q, false, bool,                isCCInclusiveMuonCandidate)                                                                        \
    f(p, q, false, int,                 pdgCode)                                                                                           \
    f(p, q, false, int,                 nHitsU)                                                                                            \
    f(p, q, false, int,                 nHitsV)                                                                                            \
    f(p, q, false, int,                 nHitsW)                                                                                            \
    f(p, q, false, int,                 nDaughters)                                                                                        \
    f(p, q, false, int,                 nDescendents)                                                                                      \
    f(p, q, false, int,                 nDescendentHitsU)                                                                                  \
    f(p, q, false, int,                 nDescendentHitsV)                                                                                  \
    f(p, q, false, int,                 nDescendentHitsW)                                                                                  \
    f(p, q, false, int,                 nHitsInLargestDescendent)                                                                          \
    f(p, q, false, float,               trackScore)                                                                                        \
    f(p, q, false, float,               startX)                                                                                            \
    f(p, q, false, float,               startY)                                                                                            \
    f(p, q, false, float,               startZ)                                                                                            \
    f(p, q, false, float,               endX)                                                                                              \
    f(p, q, false, float,               endY)                                                                                              \
    f(p, q, false, float,               endZ)                                                                                              \
    f(p, q, false, float,               directionX)                                                                                        \
    f(p, q, false, float,               directionY)                                                                                        \
    f(p, q, false, float,               directionZ)                                                                                        \
    f(p, q, false, float,               yzAngle)                                                                                           \
    f(p, q, false, float,               xyAngle)                                                                                           \
    f(p, q, false, float,               xzAngle)                                                                                           \
    f(p, q, false, float,               length)                                                                                            \
    f(p, q, false, float,               range)                                                                                             \
    f(p, q, false, float,               transverseVertexDist)                                                                              \
    f(p, q, false, float,               longitudinalVertexDist)                                                                            \
    f(p, q, false, float,               mcsMomentumForwardMuon)                                                                            \
    f(p, q, false, float,               mcsMomentumUncertaintyForwardMuon)                                                                 \
    f(p, q, false, float,               mcsLogLikelihoodForwardMuon)                                                                       \
    f(p, q, false, float,               mcsMomentumBackwardMuon)                                                                           \
    f(p, q, false, float,               mcsMomentumUncertaintyBackwardMuon)                                                                \
    f(p, q, false, float,               mcsLogLikelihoodBackwardMuon)                                                                      \
    f(p, q, false, float,               wiggliness)                                                                                        \
    f(p, q, false, int,                 nSpacePointsNearEnd)                                                                               \
    f(p, q, false, float,               likelihoodForwardMuonU)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuonV)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuonW)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuon)                                                                             \
    f(p, q, false, float,               likelihoodBackwardMuonU)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuonV)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuonW)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuon)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionU)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionV)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionW)                                                                            \
    f(p, q, false, float,               likelihoodForwardPion)                                                                             \
    f(p, q, false, float,               likelihoodBackwardPionU)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPionV)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPionW)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPion)                                                                            \
    f(p, q, false, float,               likelihoodForwardProtonU)                                                                          \
    f(p, q, false, float,               likelihoodForwardProtonV)                                                                          \
    f(p, q, false, float,               likelihoodForwardProtonW)                                                                          \
    f(p, q, false, float,               likelihoodForwardProton)                                                                           \
    f(p, q, false, float,               likelihoodBackwardProtonU)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProtonV)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProtonW)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProton)                                                                          \
    f(p, q, false, float,               likelihoodMIPU)                                                                                    \
    f(p, q, false, float,               likelihoodMIPV)                                                                                    \
    f(p, q, false, float,               likelihoodMIPW)                                                                                    \
    f(p, q, false, float,               likelihoodMIP)                                                                                     \
    f(p, q, false, float,               truncatedMeandEdxU)                                                                                \
    f(p, q, false, float,               truncatedMeandEdxV)                                                                                \
    f(p, q, false, float,               truncatedMeandEdxW)                                                                                \
    f(p, q, false, float,               truncatedMeandEdx)                                                                                 \
    f(p, q, false, std::vector<float>,  truthMatchPurities)                                                                                \
    f(p, q, false, std::vector<float>,  truthMatchCompletenesses)                                                                          \
    f(p, q, false, bool,                hasMatchedMCParticle)                                                                              \
    f(p, q, false, int,                 bestMatchedMCParticleIndex)                                                                        \
    f(p, q, false, int,                 generation) /* New PeLEE variable */                                                               \
    f(p, q, false, float,               chi2ForwardMuonW) /* New PeLEE variable */                                                         \
    f(p, q, false, float,               chi2ForwardProtonW) /* New PeLEE variable */                                                       \
    f(p, q, false, float,               distance) /* New PeLEE variable; not SCE corrected (unlike transverse- and logitudinalVertexDist)*/\
    f(p, q, false, float,               vertexX) /* New variable !ONLY USED FOR UBCC1PI INPUT*/                                            \
    f(p, q, false, float,               vertexY) /* New variable !ONLY USED FOR UBCC1PI INPUT*/                                            \
    f(p, q, false, float,               vertexZ) /* New variable !ONLY USED FOR UBCC1PI INPUT*/                                            \
    f(p, q, false, int,                 pdgBacktracked) /* New PeLEE variable */                                                           \
    f(p, q, false, float,               purityBacktracked) /* New PeLEE variable */                                                        \
    f(p, q, false, float,               completenessBacktracked) /* New PeLEE variable */                                                  \
    f(p, q, false, float,               momentumXBacktracked) /* New PeLEE variable */                                                           \
    f(p, q, false, float,               momentumYBacktracked) /* New PeLEE variable */                                                           \
    f(p, q, false, float,               momentumZBacktracked) /* New PeLEE variable */                                                           \
    // f(p, q, false, double,              vertexDistanceToXBoundary) /* New PeLEE variable !ONLY USED FOR PELEE INPUT*/
    // f(p, q, false, double,              vertexDistanceToZBoundary) /* New PeLEE variable !ONLY USED FOR PELEE INPUT*/
    // f(p, q, false, double,              vertexDistanceToYBoundary) /* New PeLEE variable !ONLY USED FOR PELEE INPUT*/




//*****************************************************************************************************************************************
//**************************************************************************************************************************************
//***********************************************************************************************************************************
// PeLEE Mode:
//***********************************************************************************************************************************
//***************************************************************************************************************************************
//*****************************************************************************************************************************************
/** The event metadata members */
#define PELEE_MACRO_EVENT_METADATA_MEMBERS(p, q, f)                                                                                        \
    f(p, q, false, int,         run)                                                                                                       \
    f(p, q, false, int,         sub) /*subRun*/                                                                                            \
    f(p, q, false, int,         evt) /*event*/

/** The event truth information members */
#define COMMA , // This is needed for the std::map macro below since it uses commas to separate arguments
#define PELEE_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(p, q, f)                                                                                                   \
    f(p, q, false, float,                                           weightSpline) /*splineEventWeight*/                                                    \
    f(p, q, false, float,                                           weightTune) /*splineEventWeight ?*/                                                    \
    f(p, q, false, float,                                           weightSplineTimesTune) /*genieTuneEventWeight*/                                        \
    f(p, q, true,  std::map<std::string COMMA std::vector<double>>, weights) /*systParamNames and systParamValues*/                                        \
    f(p, q, true,  std::vector<unsigned short>,                     weightsFlux) /*cast from u short to int*/ /*systParamFirstValueIndex*/                 \
    f(p, q, true,  std::vector<unsigned short>,                     weightsGenie) /*cast from u short to int*/ /*systParamFirstValueIndex*/                \
    f(p, q, true,  std::vector<unsigned short>,                     weightsReint) /*cast from u short to int*/ /*systParamFirstValueIndex*/

#define PELEE_MACRO_EVENT_TRUTH_MEMBERS(p, q, f)                                                                                                           \
    f(p, q, false, int,                                             ccnc) /*bool isCC*/                                                                    \
    f(p, q, false, int,                                             interaction) /*interactionMode and interactionString*/                                 \
    f(p, q, false, int,                                             nu_pdg) /*nuPDGCode*/                                                                  \
    f(p, q, false, float,                                           nu_e) /*nuEnergy*/                                                                     \
    f(p, q, false, float,                                           true_nu_vtx_x) /*Part of nuVertex*/                                                    \
    f(p, q, false, float,                                           true_nu_vtx_y) /*Part of nuVertex*/                                                    \
    f(p, q, false, float,                                           true_nu_vtx_z) /*Part of nuVertex*/                                                    \
    f(p, q, true,  std::vector<float>,                              mc_purity) /*slicePurities*/                                                           \
    f(p, q, true,  std::vector<float>,                              mc_completeness) /*sliceCompletenesses*/
    /*f(p, q, false, int,                                           nFinalStates)*/

/** The event truth particle information members */
#define PELEE_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(p, q, f)                                                                                  \
    f(p, q, false, int,                 mc_pdg) /*pdgCode*/                                                                                \
    f(p, q, false, float,               mc_vx) /*startX*/                                                                                  \
    f(p, q, false, float,               mc_vy) /*startY*/                                                                                  \
    f(p, q, false, float,               mc_vz) /*startZ*/                                                                                  \
    f(p, q, false, float,               mc_endx) /*endX*/                                                                                  \
    f(p, q, false, float,               mc_endy) /*endY*/                                                                                  \
    f(p, q, false, float,               mc_endz) /*endZ*/                                                                                  \
    f(p, q, false, float,               mc_px) /*momentumX and part of momentum*/                                                          \
    f(p, q, false, float,               mc_py) /*momentumY and part of momentum*/                                                          \
    f(p, q, false, float,               mc_pz) /*momentumZ and part of momentum*/                                                          \
    f(p, q, false, float,               mc_E) /*energy*/                                                                                   \
    f(p, q, false, float,               mc_end_p) /*endMomentum and isStopping*/                                                           \
    f(p, q, false, int,                 mc_n_elastic)                                                                                      \
    f(p, q, false, int,                 mc_n_inelastic)                                                                                    \
    // f(p, q, false, bool,                mc_nu_truth) /*Used to dismiss some mcParticles*/
    /*f(p, q, false, float,               mass)                                                                                            \
    f(p, q, false, float,               hitWeightU)                                                                                        \
    f(p, q, false, float,               hitWeightV)                                                                                        \
    f(p, q, false, float,               hitWeightW)                                                                                        \
    f(p, q, false, std::vector<float>,  scatterCosThetas)                                                                                  \
    f(p, q, false, std::vector<float>,  scatterMomentumFracsLost)                                                                          \
    f(p, q, false, std::vector<int>,    scatterIsElastic)                                                                                  \
    f(p, q, false, float,               scatteredMomentum)                                                                                 \
    f(p, q, false, int,                 endState)                                                                                          \
    f(p, q, false, float,               endStateProductsHitWeightU)                                                                        \
    f(p, q, false, float,               endStateProductsHitWeightV)                                                                        \
    f(p, q, false, float,               endStateProductsHitWeightW)*/

/** The event reco information members */
#define PELEE_MACRO_EVENT_RECO_MEMBERS(p, q, f)                                                                                            \
    f(p, q, false, int,                 filter_ccinclusive) /*passesCCInclusivex*/                                                         \
    f(p, q, false, int,                 nslice) /*nSlices*/                                                                                \
    f(p, q, false, float,               topological_score) /*selectedTopologicalScore*/                                                    \
    f(p, q, false, int,                 slpdg) /*nuPdgCode and hasNeutrino*/                                                               \
    f(p, q, false, float,               reco_nu_vtx_x) /*Part of nuVertexNoSCC*/                                                           \
    f(p, q, false, float,               reco_nu_vtx_y) /*Part of nuVertexNoSCC*/                                                           \
    f(p, q, false, float,               reco_nu_vtx_z) /*Part of nuVertexNoSCC*/                                                           \
    f(p, q, false, float,               reco_nu_vtx_sce_x) /*Part of nuVertex*/                                                            \
    f(p, q, false, float,               reco_nu_vtx_sce_y) /*Part of nuVertex*/                                                            \
    f(p, q, false, float,               reco_nu_vtx_sce_z) /*Part of nuVertex*/                                                            \
    f(p, q, false, float,               nu_flashmatch_score) /*flashChi2*/
    /*f(p, q, false, bool,                hasSelectedSlice)                                                                                \
    f(p, q, true,  std::vector<float>,  sliceTopologicalScores)                                                                            \
    f(p, q, true,  std::vector<bool>,   sliceIsSelectedAsNu)                                                                               \
    f(p, q, false, bool,                hasNeutrino)                                                                                       \
    f(p, q, false, int,                 nFinalStates)                                                                                      \
    f(p, q, false, float,               flashTime)                                                                                         \
    f(p, q, false, float,               largestFlashPE)                                                                                    \
    f(p, q, false, float,               largestFlashTime)                                                                                  \
    f(p, q, false, float,               largestFlashTimeWidth)                                                                             \
    f(p, q, false, float,               largestFlashYCtr)                                                                                  \
    f(p, q, false, float,               largestFlashYWidth)                                                                                \
    f(p, q, false, float,               largestFlashZCtr)                                                                                  \
    f(p, q, false, float,               largestFlashZWidth)*/

/** The event reco particle information members */
#define PELEE_MACRO_EVENT_RECO_PARTICLE_MEMBERS(p, q, f)                                                                                   \
    f(p, q, false, int,                 pfpdg) /*pdgCode*/                                                                                 \
    f(p, q, false, int,                 pfnplanehits_U) /*nHitsU*/                                                                         \
    f(p, q, false, int,                 pfnplanehits_V) /*nHitsV*/                                                                         \
    f(p, q, false, int,                 pfnplanehits_Y) /*nHitsW*/                                                                         \
    f(p, q, false, int,                 pfp_trk_daughters_v) /*to compute nDaughters*/                                                     \
    f(p, q, false, int,                 pfp_shr_daughters_v) /*to compute nDaughters*/                                                     \
    f(p, q, false, int,                 pfp_n_descendents_v) /*nDescendents*/                                                              \
    f(p, q, false, float,               trk_score_v) /*trackScore*/                                                                        \
    f(p, q, false, int,                 pfp_generation_v) /*used to reject granddaughter entries & generation; ATTN: original is unsigned int */ \
    f(p, q, false, float,               trk_pid_chipr_v) /*chi2ForwardProtonW*/ /*todo: check this is squared*/                            \
    f(p, q, false, float,               trk_pid_chimu_v) /*chi2ForwardMuonW*/ /*todo: check this is squared*/                              \
    f(p, q, false, float,               trk_sce_start_x_v) /*startX*/                                                                      \
    f(p, q, false, float,               trk_sce_start_y_v) /*startY*/                                                                      \
    f(p, q, false, float,               trk_sce_start_z_v) /*startZ*/                                                                      \
    f(p, q, false, float,               trk_start_x_v) /*vertexX*/                                                                         \
    f(p, q, false, float,               trk_start_y_v) /*vertexY*/                                                                         \
    f(p, q, false, float,               trk_start_z_v) /*vertexZ*/                                                                         \
    f(p, q, false, float,               shr_start_x_v) /*vertexX*/                                                                         \
    f(p, q, false, float,               shr_start_y_v) /*vertexY*/                                                                         \
    f(p, q, false, float,               shr_start_z_v) /*vertexZ*/                                                                         \
    f(p, q, false, float,               shr_tkfit_start_x_v) /*vertexX*/                                                                   \
    f(p, q, false, float,               shr_tkfit_start_y_v) /*vertexY*/                                                                   \
    f(p, q, false, float,               shr_tkfit_start_z_v) /*vertexZ*/                                                                   \
    f(p, q, false, float,               pfp_vtx_x_v) /*vertexX*/                                                                             \
    f(p, q, false, float,               pfp_vtx_y_v) /*vertexY*/                                                                             \
    f(p, q, false, float,               pfp_vtx_z_v) /*vertexZ*/                                                                             \
    f(p, q, false, float,               trk_end_x_v) /*used for vertexX*/                                                                  \
    f(p, q, false, float,               trk_end_y_v) /*used for vertexY*/                                                                  \
    f(p, q, false, float,               trk_end_z_v) /*used for vertexZ*/                                                                  \
    f(p, q, false, float,               trk_sce_end_x_v) /*endX*/                                                                          \
    f(p, q, false, float,               trk_sce_end_y_v) /*endY*/                                                                          \
    f(p, q, false, float,               trk_sce_end_z_v) /*endZ*/                                                                          \
    f(p, q, false, float,               trk_dir_x_v) /*directionX*/                                                                        \
    f(p, q, false, float,               trk_dir_y_v) /*directionY*/                                                                        \
    f(p, q, false, float,               trk_dir_z_v) /*directionZ*/                                                                        \
    f(p, q, false, float,               trk_len_v) /*range not length*/                                                                    \
    f(p, q, false, float,               trk_avg_deflection_stdev_v) /*wiggliness*/                                                         \
    f(p, q, false, int,                 trk_end_spacepoints_v) /*nSpacePointsNearEnd*/                                                     \
    f(p, q, false, float,               trk_trunk_rr_dEdx_u_v) /*truncatedMeandEdxU*/                                                      \
    f(p, q, false, float,               trk_trunk_rr_dEdx_v_v) /*truncatedMeandEdxV*/                                                      \
    f(p, q, false, float,               trk_trunk_rr_dEdx_y_v) /*truncatedMeandEdxW*/                                                      \
    f(p, q, false, float,               trk_distance_v) /*distance*/                                                                       \
    f(p, q, false, float,               trk_bragg_p_v)                                                                                     \
    f(p, q, false, float,               trk_bragg_mu_v)                                                                                    \
    f(p, q, false, float,               trk_bragg_pion_v)                                                                                  \
    f(p, q, false, float,               trk_bragg_mip_v)                                                                                   \
    f(p, q, false, float,               trk_bragg_p_u_v)                                                                                   \
    f(p, q, false, float,               trk_bragg_mu_u_v)                                                                                  \
    f(p, q, false, float,               trk_bragg_pion_u_v)                                                                                \
    f(p, q, false, float,               trk_bragg_mip_u_v)                                                                                 \
    f(p, q, false, float,               trk_bragg_p_v_v)                                                                                   \
    f(p, q, false, float,               trk_bragg_mu_v_v)                                                                                  \
    f(p, q, false, float,               trk_bragg_pion_v_v)                                                                                \
    f(p, q, false, float,               trk_bragg_mip_v_v)                                                                                 \
    f(p, q, false, float,               trk_bragg_p_alt_dir_v)                                                                             \
    f(p, q, false, float,               trk_bragg_mu_alt_dir_v)                                                                            \
    f(p, q, false, float,               trk_bragg_pion_alt_dir_v)                                                                          \
    f(p, q, false, bool,                trk_bragg_p_fwd_preferred_v)                                                                       \
    f(p, q, false, bool,                trk_bragg_mu_fwd_preferred_v)                                                                      \
    f(p, q, false, bool,                trk_bragg_pion_fwd_preferred_v)                                                                    \
    f(p, q, false, float,               trk_bragg_p_alt_dir_u_v)                                                                           \
    f(p, q, false, float,               trk_bragg_mu_alt_dir_u_v)                                                                          \
    f(p, q, false, float,               trk_bragg_pion_alt_dir_u_v)                                                                        \
    f(p, q, false, bool,                trk_bragg_p_fwd_preferred_u_v)                                                                     \
    f(p, q, false, bool,                trk_bragg_mu_fwd_preferred_u_v)                                                                    \
    f(p, q, false, bool,                trk_bragg_pion_fwd_preferred_u_v)                                                                  \
    f(p, q, false, float,               trk_bragg_p_alt_dir_v_v)                                                                           \
    f(p, q, false, float,               trk_bragg_mu_alt_dir_v_v)                                                                          \
    f(p, q, false, float,               trk_bragg_pion_alt_dir_v_v)                                                                        \
    f(p, q, false, bool,                trk_bragg_p_fwd_preferred_v_v)                                                                     \
    f(p, q, false, bool,                trk_bragg_mu_fwd_preferred_v_v)                                                                    \
    f(p, q, false, bool,                trk_bragg_pion_fwd_preferred_v_v)                                                                  \
    f(p, q, false, int,                 trk_nhits_u_v) /*Used for truncatedMeandEdx*/                                                      \
    f(p, q, false, int,                 trk_nhits_v_v) /*Used for truncatedMeandEdx*/                                                      \
    f(p, q, false, int,                 trk_nhits_y_v) /*Used for truncatedMeandEdx*/
    
    // f(p, q, false, float,               trk_bragg_p_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_mu_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_pion_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_mip_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_p_u_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_mu_u_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_pion_u_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_mip_u_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_p_v_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_mu_v_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_pion_v_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, float,               trk_bragg_mip_v_v) /*ATTN: Not same as ubcc1pi - max of fwd and bwd*/
    // f(p, q, false, double,              dvtx_x_boundary) /*vertexDistanceToXBoundary; todo: Warning conversion from double to float*/
    // f(p, q, false, double,              dvtx_y_boundary) /*vertexDistanceToYBoundary; todo: Warning conversion from double to float*/
    // f(p, q, false, double,              dvtx_z_boundary) /*vertexDistanceToZBoundary; todo: Warning conversion from double to float*/


/** The event backtracked members */
// #define PELEE_MACRO_EVENT_BACKTRACKED_MEMBERS(p, q, f)

/** The event backtracked particle information members */
#define PELEE_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(p, q, f)                                                                            \
    f(p, q, false, int,                 backtracked_pdg) /*pdgBacktracked*/                                                                \
    f(p, q, false, float,               backtracked_purity) /*purityBacktracked*/                                                          \
    f(p, q, false, float,               backtracked_completeness) /*completenessBacktracked*/                                              \
    f(p, q, false, float,               backtracked_px) /*momentumXBacktracked*/                                                           \
    f(p, q, false, float,               backtracked_py) /*momentumYBacktracked*/                                                           \
    f(p, q, false, float,               backtracked_pz) /*momentumZBacktracked*/

    /*f(p, q, false, bool,                isCCInclusiveMuonCandidate)                                                                      \
    f(p, q, false, int,                 nDescendentHitsU)                                                                                  \
    f(p, q, false, int,                 nDescendentHitsV)                                                                                  \
    f(p, q, false, int,                 nDescendentHitsW)                                                                                  \
    f(p, q, false, int,                 nHitsInLargestDescendent)                                                                          \
    f(p, q, false, float,               yzAngle)                                                                                           \
    f(p, q, false, float,               xyAngle)                                                                                           \
    f(p, q, false, float,               xzAngle)                                                                                           \
    f(p, q, false, float,               range)                                                                                             \
    f(p, q, false, float,               transverseVertexDist)                                                                              \
    f(p, q, false, float,               longitudinalVertexDist)                                                                            \
    f(p, q, false, float,               mcsMomentumForwardMuon)                                                                            \
    f(p, q, false, float,               mcsMomentumUncertaintyForwardMuon)                                                                 \
    f(p, q, false, float,               mcsLogLikelihoodForwardMuon)                                                                       \
    f(p, q, false, float,               mcsMomentumBackwardMuon)                                                                           \
    f(p, q, false, float,               mcsMomentumUncertaintyBackwardMuon)                                                                \
    f(p, q, false, float,               mcsLogLikelihoodBackwardMuon)                                                                      \
    f(p, q, false, float,               likelihoodForwardMuonU)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuonV)                                                                            \
    f(p, q, false, float,               likelihoodForwardMuonW)                                                                            \
    f(p, q, false, float,               likelihoodBackwardMuonU)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuonV)                                                                           \
    f(p, q, false, float,               likelihoodBackwardMuonW)                                                                           \
    f(p, q, false, float,               likelihoodForwardPionU)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionV)                                                                            \
    f(p, q, false, float,               likelihoodForwardPionW)                                                                            \
    f(p, q, false, float,               likelihoodBackwardPionU)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPionV)                                                                           \
    f(p, q, false, float,               likelihoodBackwardPionW)                                                                           \
    f(p, q, false, float,               likelihoodForwardProtonU)                                                                          \
    f(p, q, false, float,               likelihoodForwardProtonV)                                                                          \
    f(p, q, false, float,               likelihoodForwardProtonW)                                                                          \
    f(p, q, false, float,               likelihoodBackwardProtonU)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProtonV)                                                                         \
    f(p, q, false, float,               likelihoodBackwardProtonW)                                                                         \
    f(p, q, false, float,               likelihoodMIPU)                                                                                    \
    f(p, q, false, float,               likelihoodMIPV)                                                                                    \
    f(p, q, false, float,               likelihoodMIPW)                                                                                    \
    f(p, q, false, float,               truncatedMeandEdx)                                                                                 \
    f(p, q, false, std::vector<float>,  truthMatchPurities)                                                                                \
    f(p, q, false, std::vector<float>,  truthMatchCompletenesses)                                                                          \
    f(p, q, false, bool,                hasMatchedMCParticle)                                                                              \
    f(p, q, false, int,                 bestMatchedMCParticleIndex)*/


//*****************************************************************************************************************************************
//**************************************************************************************************************************************
//***********************************************************************************************************************************
// XSec Mode:
//***********************************************************************************************************************************
//***************************************************************************************************************************************
//*****************************************************************************************************************************************
// /** The event metadata members */
// #define XSEC_MACRO_EVENT_METADATA_MEMBERS(p, q, f)
//     f(p, q, false, int,         run)

// /** The event truth information members */
// #define COMMA , // This is needed for the std::map macro below since it uses commas to separate arguments
// #define XSEC_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS(p, q, f)
//     f(p, q, false, float,                                           weightSpline) /*splineEventWeight*/

#define XSEC_MACRO_EVENT_TRUTH_MEMBERS(p, q, f) \
    f(p, q, false, int,                        run) \
    f(p, q, false, int,                        category) \
    f(p, q, false, bool,                       true_cc0pi) \
    f(p, q, false, bool,                       true_cc1pi) \
    f(p, q, false, bool,                       cc0pi_signal) \
    f(p, q, false, bool,                       cc1pi_signal) \
    f(p, q, false, float,                      cc1pi_truth_muonMomentum) \
    f(p, q, false, float,                      cc1pi_truth_muonCosTheta) \
    f(p, q, false, float,                      cc1pi_truth_muonPhi) \
    f(p, q, false, float,                      cc1pi_truth_pionMomentum) \
    f(p, q, false, float,                      cc1pi_truth_pionCosTheta) \
    f(p, q, false, float,                      cc1pi_truth_pionPhi) \
    f(p, q, false, float,                      cc1pi_truth_muonPionAngle) \
    f(p, q, false, int,                        cc1pi_truth_nProtons) \
    f(p, q, false, float,                      cc0pi_truth_muonMomentum) \
    f(p, q, false, float,                      cc0pi_truth_muonCosTheta) \
    f(p, q, false, float,                      cc0pi_truth_muonPhi) \
    f(p, q, false, float,                      cc0pi_truth_pionMomentum) \
    f(p, q, false, float,                      cc0pi_truth_pionCosTheta) \
    f(p, q, false, float,                      cc0pi_truth_pionPhi) \
    f(p, q, false, float,                      cc0pi_truth_muonPionAngle) \
    f(p, q, false, int,                        cc0pi_truth_nProtons)

// /** The event truth particle information members */
// #define XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(p, q, f)
//     f(p, q, false, int,                 mc_pdg) /*pdgCode*/

/** The event reco information members */
#define XSEC_MACRO_EVENT_RECO_MEMBERS(p, q, f) \
    f(p, q, false, bool,                       cc0pi_selected_generic) \
    f(p, q, false, bool,                       cc0pi_selected_golden) \
    f(p, q, false, bool,                       cc1pi_selected_generic) \
    f(p, q, false, bool,                       cc1pi_selected_golden) \
    f(p, q, false, float,                      cc1pi_reco_muonMomentum) \
    f(p, q, false, float,                      cc1pi_reco_muonCosTheta) \
    f(p, q, false, float,                      cc1pi_reco_muonPhi) \
    f(p, q, false, float,                      cc1pi_reco_pionMomentum) \
    f(p, q, false, float,                      cc1pi_reco_pionCosTheta) \
    f(p, q, false, float,                      cc1pi_reco_pionPhi) \
    f(p, q, false, float,                      cc1pi_reco_muonPionAngle) \
    f(p, q, false, float,                      cc1pi_reco_nProtons) \
    f(p, q, false, float,                      cc0pi_reco_muonMomentum) \
    f(p, q, false, float,                      cc0pi_reco_muonCosTheta) \
    f(p, q, false, float,                      cc0pi_reco_muonPhi) \
    f(p, q, false, float,                      cc0pi_reco_pionMomentum) \
    f(p, q, false, float,                      cc0pi_reco_pionCosTheta) \
    f(p, q, false, float,                      cc0pi_reco_pionPhi) \
    f(p, q, false, float,                      cc0pi_reco_muonPionAngle) \
    f(p, q, false, float,                      cc0pi_reco_nProtons) \
    f(p, q, false, bool,                       passed_particleTrackScore) \
    f(p, q, false, bool,                       passed_particleVertexDistance) \
    f(p, q, false, bool,                       passed_particleGeneration) \
    f(p, q, false, bool,                       passed_particleTrackLength) \
    f(p, q, false, bool,                       passed_particleProtonChi2) \
    f(p, q, false, bool,                       passed_particleMuonChi2) \
    f(p, q, false, bool,                       passed_particleProtonChi2OverMuonChi2) \
    f(p, q, false, bool,                       passed_pandoraNuPDGIsNumu) \
    f(p, q, false, bool,                       passed_daughterVerticesContained) \
    f(p, q, false, bool,                       passed_nuVertexFiducial) \
    f(p, q, false, bool,                       passed_topologicalOrFlashMatch) \
    f(p, q, false, bool,                       passed_min2Tracks) \
    f(p, q, false, bool,                       passed_max1Uncontained) \
    f(p, q, false, bool,                       passed_2NonProtons) \
    f(p, q, false, bool,                       passed_pionHasValiddEdx) \
    f(p, q, false, bool,                       passed_pionNotInGap) \
    f(p, q, false, bool,                       passed_muonNotInGap) \
    f(p, q, false, bool,                       passed_openingAngle) \
    f(p, q, false, bool,                       passed_topologicalScore) \
    f(p, q, false, bool,                       passed_startNearVertex) \
    f(p, q, false, bool,                       passed_likelyGoldenPion)


// /** The event reco particle information members */
// #define XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(p, q, f)
//     f(p, q, false, int,                 pfpdg) /*pdgCode*/






// pTrackScoreCutValues
// pVertexDistanceCutValues
// pGenerationCutValues
// pTrackLengthCutValues
// pProtonChi2CutValues
// pMuonChi2CutValues
// event_cutValue_nuPdgCode
// event_cutValue_topologicalScore
// event_cutValue_flashChi2
// event_cutValue_pionTruncatedMeandEdx
// event_cutValue_openingAngle
// event_cutValue_maxVertexDist
// event_cutValue_goldenPionBDT

// pBestMatchedTruthPDGs
// mc_nu_pdg
// spline_weight
// tuned_cv_weight


#endif
