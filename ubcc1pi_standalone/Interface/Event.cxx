/**
 *  @file  ubcc1pi_standalone/Interface/Event.cxx
 *
 *  @brief The implementation of the event class
 */

#include "ubcc1pi_standalone/Interface/Event.h"
// #include "ubcc1pi_standalone/Interface/EventPeLEE.h"
#include <ctime> // DEBUG

namespace ubcc1pi
{

Event::Event(const bool hasTruthInfo)
{
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_INIT_MEMBER_VECTOR)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_INIT_MEMBER_VECTOR)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

Event::Event(const EventPeLEE& eventPeLEE, const bool excludeGranddaughterParticles /*= false*/)
{
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_INIT_MEMBER_VECTOR)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_INIT_MEMBER_VECTOR)

    const auto hasTruthInfo = eventPeLEE.HasTruthWeights();
    const int intMinValue = -2147483648; // hardcoded to match root type
    const float floatMinValue = -340282346638528859811704183484516925440.000000f;
    // const double doubleMinValue = -1.7976931348623158e+308;

    // PELEE_TO_UBCC1PI_MACRO_EVENT_METADATA(eventPeLEE.metadata, metadata, PELEE_TO_UBCC1PI_MEMBER_CONVERSION)
    // PELEE_TO_UBCC1PI_MACRO_EVENT_TRUTH(eventPeLEE.truth, truth, PELEE_TO_UBCC1PI_MEMBER_CONVERSION)

    // ************* EVENT METADATA *************
    metadata.run.Set(eventPeLEE.metadata.run());
    metadata.subRun.Set(eventPeLEE.metadata.sub());
    metadata.event.Set(eventPeLEE.metadata.evt());
    metadata.hasTruthInfo.Set(hasTruthInfo);

    if(hasTruthInfo)
    {
        // ************* EVENT TRUTH *************
        if(eventPeLEE.truth.weightSpline() != floatMinValue) truth.splineEventWeight.Set(eventPeLEE.truth.weightSpline());
        if(eventPeLEE.truth.weightTune() != floatMinValue) truth.genieTuneEventWeight.Set(eventPeLEE.truth.weightTune());
        // std::cout<<"DEBUG Event.cxx: weightSpline: "<<eventPeLEE.truth.weightSpline()<<" weightTune: "<<eventPeLEE.truth.weightTune()<<std::endl;

        std::vector<std::string> systParamNames;
        std::vector<int> systParamFirstValueIndex;
        std::vector<float> systParamValues;
        for (const auto &[paramName, values]: eventPeLEE.truth.weights())
        {
            systParamFirstValueIndex.push_back(systParamValues.size());
            systParamNames.push_back(paramName);
            for (const auto &value: values)
            {
                systParamValues.push_back(value);
            }
        }
        truth.systParamNames.Set(systParamNames);
        truth.systParamFirstValueIndex.Set(systParamFirstValueIndex); // Might need to use weightsFlux, weightsGenie, weightsReint at some point to sort out the indices
        truth.systParamValues.Set(systParamValues);
        truth.isCC.Set(eventPeLEE.truth.ccnc() == 0);
        truth.interactionMode.Set(eventPeLEE.truth.interaction());
        // truth.interactionString = DebugHelper::GetInteractionString(eventPeLEE.truth.interaction(), true); // To do: Needs to be implemented in DebugHelper
        if(eventPeLEE.truth.nu_pdg() != intMinValue) truth.nuPdgCode.Set(eventPeLEE.truth.nu_pdg());
        if(eventPeLEE.truth.nu_e() != floatMinValue) truth.nuEnergy.Set(eventPeLEE.truth.nu_e());
        if(eventPeLEE.truth.true_nu_vtx_x() != floatMinValue && eventPeLEE.truth.true_nu_vtx_y() != floatMinValue && eventPeLEE.truth.true_nu_vtx_z() != floatMinValue) truth.nuVertex.Set(TVector3(eventPeLEE.truth.true_nu_vtx_x(), eventPeLEE.truth.true_nu_vtx_y(), eventPeLEE.truth.true_nu_vtx_z()));
        // truth.nFinalStates;
        // truth.slicePurities.Set(eventPeLEE.truth.mc_purity());
        // truth.sliceCompletenesses.Set(eventPeLEE.truth.mc_completeness());

        // ************* EVENT TRUTH PARTICLE *************
        unsigned int nTruthParticles = eventPeLEE.truth.particles.size();

        unsigned int nTruthParticlesToSkip = 0u;
        // for (unsigned int i = 0; i < nTruthParticles; ++i) {
        //     if (!eventPeLEE.truth.particles.at(i).mc_nu_truth())
        //     {
        //         nTruthParticlesToSkip++; // Count mcParticles that do not originate from right mcTruth
        //     }
        // }

        truth.particles.resize(nTruthParticles - nTruthParticlesToSkip);
        auto truthTargetIndex = 0u;
        for(unsigned int i = 0; i<nTruthParticles; ++i)
        {
            const auto pParticlePeLEE = &eventPeLEE.truth.particles.at(i);
            // if (!pParticlePeLEE->mc_nu_truth())
            // {
            //     continue; // Skip mcParticles that do not originate from right mcTruth
            // }
            const auto pParticle = &truth.particles.at(truthTargetIndex);
            truthTargetIndex++;

            if(pParticlePeLEE->mc_pdg() != intMinValue) pParticle->pdgCode.Set(pParticlePeLEE->mc_pdg());
            if(pParticlePeLEE->mc_vx() != floatMinValue) pParticle->startX.Set(pParticlePeLEE->mc_vx());
            if(pParticlePeLEE->mc_vy() != floatMinValue) pParticle->startY.Set(pParticlePeLEE->mc_vy());
            if(pParticlePeLEE->mc_vz() != floatMinValue) pParticle->startZ.Set(pParticlePeLEE->mc_vz());
            if(pParticlePeLEE->mc_endx() != floatMinValue) pParticle->endX.Set(pParticlePeLEE->mc_endx());
            if(pParticlePeLEE->mc_endy() != floatMinValue) pParticle->endY.Set(pParticlePeLEE->mc_endy());
            if(pParticlePeLEE->mc_endz() != floatMinValue) pParticle->endZ.Set(pParticlePeLEE->mc_endz());
            if(pParticlePeLEE->mc_px() != floatMinValue) pParticle->momentumX.Set(pParticlePeLEE->mc_px());
            if(pParticlePeLEE->mc_py() != floatMinValue) pParticle->momentumY.Set(pParticlePeLEE->mc_py());
            if(pParticlePeLEE->mc_pz() != floatMinValue) pParticle->momentumZ.Set(pParticlePeLEE->mc_pz());
            const auto momentumVector = std::hypot(pParticlePeLEE->mc_px(), pParticlePeLEE->mc_py(), pParticlePeLEE->mc_pz());
            if(pParticlePeLEE->mc_px() != floatMinValue && pParticlePeLEE->mc_py() != floatMinValue && pParticlePeLEE->mc_pz() != floatMinValue) pParticle->momentum.Set(momentumVector);
            if(pParticlePeLEE->mc_E() != floatMinValue) pParticle->energy.Set(pParticlePeLEE->mc_E());
            if(pParticlePeLEE->mc_end_p() != floatMinValue) pParticle->endMomentum.Set(pParticlePeLEE->mc_end_p());
            if(pParticlePeLEE->mc_end_p() != floatMinValue) pParticle->isStopping.Set(pParticlePeLEE->mc_end_p() <= std::numeric_limits<float>::epsilon());
            if(pParticlePeLEE->mc_n_elastic() != intMinValue) pParticle->nElasticScatters.Set(pParticlePeLEE->mc_n_elastic());
            if(pParticlePeLEE->mc_n_inelastic() != intMinValue) pParticle->nInelasticScatters.Set(pParticlePeLEE->mc_n_inelastic());

            // pParticle->mass.Set(pParticlePeLEE->());
            // pParticle->hitWeightU.Set(pParticlePeLEE->());
            // pParticle->hitWeightV.Set(pParticlePeLEE->());
            // pParticle->hitWeightW.Set(pParticlePeLEE->());
            // pParticle->scatterCosThetas.Set(pParticlePeLEE->());
            // pParticle->scatterMomentumFracsLost.Set(pParticlePeLEE->());
            // pParticle->scatterIsElastic.Set(pParticlePeLEE->());
            // pParticle->scatteredMomentum.Set(pParticlePeLEE->());
            // pParticle->endState.Set(pParticlePeLEE->());
            // pParticle->endStateProductsHitWeightU.Set(pParticlePeLEE->());
            // pParticle->endStateProductsHitWeightV.Set(pParticlePeLEE->());
            // pParticle->endStateProductsHitWeightW.Set(pParticlePeLEE->());
        }
    }
    

    // ************* EVENT RECO *************
    reco.passesCCInclusive.Set(eventPeLEE.reco.filter_ccinclusive());
    // if(eventPeLEE.reco.nslice() != intMinValue) reco.nSlices.Set(eventPeLEE.reco.nslice());
    // reco.hasSelectedSlice.Set(eventPeLEE.reco..());
    if(eventPeLEE.reco.topological_score() != floatMinValue) reco.selectedTopologicalScore.Set(eventPeLEE.reco.topological_score());
    // reco.sliceTopologicalScores.Set(eventPeLEE.reco..());
    // reco.sliceIsSelectedAsNu.Set(eventPeLEE.reco..());
    // reco.hasNeutrino.Set(eventPeLEE.reco..());
    if(eventPeLEE.reco.slpdg() != intMinValue) reco.nuPdgCode.Set(eventPeLEE.reco.slpdg());
    if(eventPeLEE.reco.reco_nu_vtx_sce_x() != floatMinValue && eventPeLEE.reco.reco_nu_vtx_sce_y() != floatMinValue && eventPeLEE.reco.reco_nu_vtx_sce_z() != floatMinValue) reco.nuVertex.Set(TVector3(eventPeLEE.reco.reco_nu_vtx_sce_x(), eventPeLEE.reco.reco_nu_vtx_sce_y(), eventPeLEE.reco.reco_nu_vtx_sce_z()));
    // reco.nFinalStates.Set(eventPeLEE.reco..());
    if(eventPeLEE.reco.nu_flashmatch_score() < 1e6f && eventPeLEE.reco.nu_flashmatch_score() != floatMinValue) reco.flashChi2.Set(eventPeLEE.reco.nu_flashmatch_score()); // The not-set value is a bit unusual
    // reco.flashTime.Set(eventPeLEE.reco..());
    // reco.largestFlashPE.Set(eventPeLEE.reco..());
    // reco.largestFlashTime.Set(eventPeLEE.reco..());
    // reco.largestFlashTimeWidth.Set(eventPeLEE.reco..());
    // reco.largestFlashYCtr.Set(eventPeLEE.reco..());
    // reco.largestFlashYWidth.Set(eventPeLEE.reco..());
    // reco.largestFlashZCtr.Set(eventPeLEE.reco..());
    // reco.largestFlashZWidth.Set(eventPeLEE.reco..());
    if(eventPeLEE.reco.reco_nu_vtx_x() != floatMinValue && eventPeLEE.reco.reco_nu_vtx_y() != floatMinValue && eventPeLEE.reco.reco_nu_vtx_z() != floatMinValue) reco.nuVertexNoSCC.Set(TVector3(eventPeLEE.reco.reco_nu_vtx_x(), eventPeLEE.reco.reco_nu_vtx_y(), eventPeLEE.reco.reco_nu_vtx_z()));


    // ************* EVENT RECO PARTICLE *************
    unsigned int nRecoParticles = eventPeLEE.reco.particles.size();

    unsigned int nRecoParticlesToSkip = 0u;
    if(excludeGranddaughterParticles)
    {
        for (unsigned int i = 0; i < nRecoParticles; ++i) {
            if (eventPeLEE.reco.particles.at(i).pfp_generation_v() > 2) {
                nRecoParticlesToSkip++; // Skip granddaughter and higher order descendants
            }
        }
    }
    reco.particles.resize(nRecoParticles - nRecoParticlesToSkip);

    const auto nBacktrackedParticles = eventPeLEE.backtracked.particles.size();
    const auto hasBacktrackedInfo = nBacktrackedParticles > 0u;
    if(hasBacktrackedInfo && nBacktrackedParticles != nRecoParticles) throw std::logic_error("EventPeLEE::EventPeLEE - The number of reco particles does not match the number of backtracked particles");

    auto recoTargetIndex = 0u;
    for (unsigned int i = 0; i < nRecoParticles; ++i)
    {
        const auto pParticlePeLEE = &eventPeLEE.reco.particles.at(i);
        if (excludeGranddaughterParticles && pParticlePeLEE->pfp_generation_v() > 2) continue; // Skip granddaughter and higher order descendants
        const auto pParticle = &reco.particles.at(recoTargetIndex);
        recoTargetIndex++;

        // pParticle->isCCInclusiveMuonCandidate.Set(pParticlePeLEE->);
        if(pParticlePeLEE->pfpdg() != intMinValue) pParticle->pdgCode.Set(pParticlePeLEE->pfpdg());
        if(pParticlePeLEE->pfnplanehits_U() != intMinValue) pParticle->nHitsU.Set(pParticlePeLEE->pfnplanehits_U());
        if(pParticlePeLEE->pfnplanehits_V() != intMinValue) pParticle->nHitsV.Set(pParticlePeLEE->pfnplanehits_V());
        if(pParticlePeLEE->pfnplanehits_Y() != intMinValue) pParticle->nHitsW.Set(pParticlePeLEE->pfnplanehits_Y());
        if(pParticlePeLEE->pfp_trk_daughters_v() != intMinValue && pParticlePeLEE->pfp_shr_daughters_v() != intMinValue) pParticle->nDaughters.Set(pParticlePeLEE->pfp_trk_daughters_v() + pParticlePeLEE->pfp_shr_daughters_v());
        if(pParticlePeLEE->pfp_n_descendents_v() != intMinValue) pParticle->nDescendents.Set(pParticlePeLEE->pfp_n_descendents_v());
        // pParticle->nDescendentHitsU.Set(pParticlePeLEE->);
        // pParticle->nDescendentHitsV.Set(pParticlePeLEE->);
        // pParticle->nDescendentHitsW.Set(pParticlePeLEE->);
        // pParticle->nHitsInLargestDescendent.Set(pParticlePeLEE->);
        if(pParticlePeLEE->trk_score_v() >= 0.f) pParticle->trackScore.Set(pParticlePeLEE->trk_score_v());
        if(pParticlePeLEE->trk_sce_start_x_v() != floatMinValue) pParticle->startX.Set(pParticlePeLEE->trk_sce_start_x_v());
        if(pParticlePeLEE->trk_sce_start_y_v() != floatMinValue) pParticle->startY.Set(pParticlePeLEE->trk_sce_start_y_v());
        if(pParticlePeLEE->trk_sce_start_z_v() != floatMinValue) pParticle->startZ.Set(pParticlePeLEE->trk_sce_start_z_v());
        if(pParticlePeLEE->trk_sce_end_x_v() != floatMinValue) pParticle->endX.Set(pParticlePeLEE->trk_sce_end_x_v());
        if(pParticlePeLEE->trk_sce_end_y_v() != floatMinValue) pParticle->endY.Set(pParticlePeLEE->trk_sce_end_y_v());
        if(pParticlePeLEE->trk_sce_end_z_v() != floatMinValue) pParticle->endZ.Set(pParticlePeLEE->trk_sce_end_z_v());
        if(pParticlePeLEE->trk_dir_x_v() != floatMinValue) pParticle->directionX.Set(pParticlePeLEE->trk_dir_x_v());
        if(pParticlePeLEE->trk_dir_y_v() != floatMinValue) pParticle->directionY.Set(pParticlePeLEE->trk_dir_y_v());
        if(pParticlePeLEE->trk_dir_z_v() != floatMinValue) pParticle->directionZ.Set(pParticlePeLEE->trk_dir_z_v());
        // else std::cout<<"DEBUG Event.cxx directionZ: "<<pParticlePeLEE->trk_dir_z_v()<<std::endl;
        if(pParticlePeLEE->trk_dir_z_v() != floatMinValue && pParticlePeLEE->trk_dir_y_v() != floatMinValue) pParticle->yzAngle.Set(std::atan2(pParticlePeLEE->trk_dir_z_v(), pParticlePeLEE->trk_dir_y_v()));
        if(pParticlePeLEE->trk_dir_y_v()!= floatMinValue && pParticlePeLEE->trk_dir_x_v() != floatMinValue) pParticle->xyAngle.Set(std::atan2(pParticlePeLEE->trk_dir_y_v(), pParticlePeLEE->trk_dir_x_v()));
        if(pParticlePeLEE->trk_dir_z_v()!= floatMinValue && pParticlePeLEE->trk_dir_x_v() != floatMinValue) pParticle->xzAngle.Set(std::atan2(pParticlePeLEE->trk_dir_z_v(), pParticlePeLEE->trk_dir_x_v()));
        // if(pParticlePeLEE->trk_len_v() != floatMinValue) pParticle->length.Set(pParticlePeLEE->trk_len_v());
        if(pParticlePeLEE->trk_len_v() != floatMinValue) pParticle->range.Set(pParticlePeLEE->trk_len_v());

        // pParticle->transverseVertexDist.Set(pParticlePeLEE->);
        // pParticle->longitudinalVertexDist.Set(pParticlePeLEE->);
        // pParticle->mcsMomentumForwardMuon.Set(pParticlePeLEE->);
        // pParticle->mcsMomentumUncertaintyForwardMuon.Set(pParticlePeLEE->);
        // pParticle->mcsLogLikelihoodForwardMuon.Set(pParticlePeLEE->);
        // pParticle->mcsMomentumBackwardMuon.Set(pParticlePeLEE->);
        // pParticle->mcsMomentumUncertaintyBackwardMuon.Set(pParticlePeLEE->);
        // pParticle->mcsLogLikelihoodBackwardMuon.Set(pParticlePeLEE->);
        if(pParticlePeLEE->trk_avg_deflection_stdev_v() != floatMinValue) pParticle->wiggliness.Set(pParticlePeLEE->trk_avg_deflection_stdev_v());
        if(pParticlePeLEE->trk_end_spacepoints_v() != intMinValue) pParticle->nSpacePointsNearEnd.Set(pParticlePeLEE->trk_end_spacepoints_v());


        const auto protonLikelihoodFwdW = pParticlePeLEE->trk_bragg_p_fwd_preferred_v() ? pParticlePeLEE->trk_bragg_p_v() : pParticlePeLEE->trk_bragg_p_alt_dir_v();
        if(protonLikelihoodFwdW > -1.f) pParticle->likelihoodForwardProtonW.Set(protonLikelihoodFwdW);
        const auto protonLikelihoodBwdW = !pParticlePeLEE->trk_bragg_p_fwd_preferred_v() ? pParticlePeLEE->trk_bragg_p_v() : pParticlePeLEE->trk_bragg_p_alt_dir_v();
        if(protonLikelihoodBwdW > -1.f) pParticle->likelihoodBackwardProtonW.Set(protonLikelihoodBwdW);
        const auto protonLikelihoodFwdU = pParticlePeLEE->trk_bragg_p_fwd_preferred_u_v() ? pParticlePeLEE->trk_bragg_p_u_v() : pParticlePeLEE->trk_bragg_p_alt_dir_u_v();
        if(protonLikelihoodFwdU > -1.f) pParticle->likelihoodForwardProtonU.Set(protonLikelihoodFwdU);
        const auto protonLikelihoodBwdU = !pParticlePeLEE->trk_bragg_p_fwd_preferred_u_v() ? pParticlePeLEE->trk_bragg_p_u_v() : pParticlePeLEE->trk_bragg_p_alt_dir_u_v();
        if(protonLikelihoodBwdU > -1.f) pParticle->likelihoodBackwardProtonU.Set(protonLikelihoodBwdU);
        const auto protonLikelihoodFwdV = pParticlePeLEE->trk_bragg_p_fwd_preferred_v_v() ? pParticlePeLEE->trk_bragg_p_v_v() : pParticlePeLEE->trk_bragg_p_alt_dir_v_v();
        if(protonLikelihoodFwdV > -1.f) pParticle->likelihoodForwardProtonV.Set(protonLikelihoodFwdV);
        const auto protonLikelihoodBwdV = !pParticlePeLEE->trk_bragg_p_fwd_preferred_v_v() ? pParticlePeLEE->trk_bragg_p_v_v() : pParticlePeLEE->trk_bragg_p_alt_dir_v_v();
        if(protonLikelihoodBwdV > -1.f) pParticle->likelihoodBackwardProtonV.Set(protonLikelihoodBwdV);

        const auto pionLikelihoodFwdW = pParticlePeLEE->trk_bragg_pion_fwd_preferred_v() ? pParticlePeLEE->trk_bragg_pion_v() : pParticlePeLEE->trk_bragg_pion_alt_dir_v();
        if(pionLikelihoodFwdW > -1.f) pParticle->likelihoodForwardPionW.Set(pionLikelihoodFwdW);
        const auto pionLikelihoodBwdW = !pParticlePeLEE->trk_bragg_pion_fwd_preferred_v() ? pParticlePeLEE->trk_bragg_pion_v() : pParticlePeLEE->trk_bragg_pion_alt_dir_v();
        if(pionLikelihoodBwdW > -1.f) pParticle->likelihoodBackwardPionW.Set(pionLikelihoodBwdW);
        const auto pionLikelihoodFwdU = pParticlePeLEE->trk_bragg_pion_fwd_preferred_u_v() ? pParticlePeLEE->trk_bragg_pion_u_v() : pParticlePeLEE->trk_bragg_pion_alt_dir_u_v();
        if(pionLikelihoodFwdU > -1.f) pParticle->likelihoodForwardPionU.Set(pionLikelihoodFwdU);
        const auto pionLikelihoodBwdU = !pParticlePeLEE->trk_bragg_pion_fwd_preferred_u_v() ? pParticlePeLEE->trk_bragg_pion_u_v() : pParticlePeLEE->trk_bragg_pion_alt_dir_u_v();
        if(pionLikelihoodBwdU > -1.f) pParticle->likelihoodBackwardPionU.Set(pionLikelihoodBwdU);
        const auto pionLikelihoodFwdV = pParticlePeLEE->trk_bragg_pion_fwd_preferred_v_v() ? pParticlePeLEE->trk_bragg_pion_v_v() : pParticlePeLEE->trk_bragg_pion_alt_dir_v_v();
        if(pionLikelihoodFwdV > -1.f) pParticle->likelihoodForwardPionV.Set(pionLikelihoodFwdV);
        const auto pionLikelihoodBwdV = !pParticlePeLEE->trk_bragg_pion_fwd_preferred_v_v() ? pParticlePeLEE->trk_bragg_pion_v_v() : pParticlePeLEE->trk_bragg_pion_alt_dir_v_v();
        if(pionLikelihoodBwdV > -1.f) pParticle->likelihoodBackwardPionV.Set(pionLikelihoodBwdV);

        const auto muonLikelihoodFwdW = pParticlePeLEE->trk_bragg_mu_fwd_preferred_v() ? pParticlePeLEE->trk_bragg_mu_v() : pParticlePeLEE->trk_bragg_mu_alt_dir_v();
        if(muonLikelihoodFwdW > -1.f) pParticle->likelihoodForwardMuonW.Set(muonLikelihoodFwdW);
        const auto muonLikelihoodBwdW = !pParticlePeLEE->trk_bragg_mu_fwd_preferred_v() ? pParticlePeLEE->trk_bragg_mu_v() : pParticlePeLEE->trk_bragg_mu_alt_dir_v();
        if(muonLikelihoodBwdW > -1.f) pParticle->likelihoodBackwardMuonW.Set(muonLikelihoodBwdW);
        const auto muonLikelihoodFwdU = pParticlePeLEE->trk_bragg_mu_fwd_preferred_u_v() ? pParticlePeLEE->trk_bragg_mu_u_v() : pParticlePeLEE->trk_bragg_mu_alt_dir_u_v();
        if(muonLikelihoodFwdU > -1.f) pParticle->likelihoodForwardMuonU.Set(muonLikelihoodFwdU);
        const auto muonLikelihoodBwdU = !pParticlePeLEE->trk_bragg_mu_fwd_preferred_u_v() ? pParticlePeLEE->trk_bragg_mu_u_v() : pParticlePeLEE->trk_bragg_mu_alt_dir_u_v();
        if(muonLikelihoodBwdU > -1.f) pParticle->likelihoodBackwardMuonU.Set(muonLikelihoodBwdU);
        const auto muonLikelihoodFwdV = pParticlePeLEE->trk_bragg_mu_fwd_preferred_v_v() ? pParticlePeLEE->trk_bragg_mu_v_v() : pParticlePeLEE->trk_bragg_mu_alt_dir_v_v();
        if(muonLikelihoodFwdV > -1.f) pParticle->likelihoodForwardMuonV.Set(muonLikelihoodFwdV);
        const auto muonLikelihoodBwdV = !pParticlePeLEE->trk_bragg_mu_fwd_preferred_v_v() ? pParticlePeLEE->trk_bragg_mu_v_v() : pParticlePeLEE->trk_bragg_mu_alt_dir_v_v();
        if(muonLikelihoodBwdV > -1.f) pParticle->likelihoodBackwardMuonV.Set(muonLikelihoodBwdV);

        const auto mipLikelihoodW = pParticlePeLEE->trk_bragg_mip_v();
        if(mipLikelihoodW > -1.f) pParticle->likelihoodMIPW.Set(mipLikelihoodW);
        const auto mipLikelihoodU = pParticlePeLEE->trk_bragg_mip_u_v();
        if(mipLikelihoodU > -1.f) pParticle->likelihoodMIPU.Set(mipLikelihoodU);
        const auto mipLikelihoodV = pParticlePeLEE->trk_bragg_mip_v_v();
        if(mipLikelihoodV > -1.f) pParticle->likelihoodMIPV.Set(mipLikelihoodV);

        // Values match default ubcc1pi ntuples
        const auto sin2AngleThreshold = 0.175f;
        const auto GetMultiplaneBraggLikelihood = [pParticle, sin2AngleThreshold](const float &likelihoodU, const float &likelihoodV, const float &likelihoodW)
        {
            if(!pParticle->yzAngle.IsSet())
                return -std::numeric_limits<float>::max();

            const auto piBy3 = std::acos(0.5f);
            const bool isTrackAlongWWire = (std::pow(std::sin(pParticle->yzAngle()), 2) < sin2AngleThreshold);
            const auto hasW = (likelihoodW > -1.f && !isTrackAlongWWire);

            if (hasW)
                return likelihoodW;

            // Otherwise the track is along the W wire direction, so just average the other two planes weighted by the number of degrees of freedom
            const bool isTrackAlongUWire = (std::pow(std::sin(pParticle->yzAngle() - piBy3), 2) < sin2AngleThreshold);
            const auto hasU = (likelihoodU > -1.f && !isTrackAlongUWire);

            const bool isTrackAlongVWire = (std::pow(std::sin(pParticle->yzAngle() + piBy3), 2) < sin2AngleThreshold);
            const auto hasV = (likelihoodV > -1.f && !isTrackAlongVWire);

            const auto nHitsU = pParticle->nHitsU();
            const auto nHitsV = pParticle->nHitsV();

            const auto dofU = hasU ? static_cast<float>(nHitsU) : 0.f;
            const auto dofV = hasV ? static_cast<float>(nHitsV) : 0.f;
            const auto dofUV = dofU + dofV;

            const auto weightU = hasU ? (likelihoodU * dofU) : 0.f;
            const auto weightV = hasV ? (likelihoodV * dofV) : 0.f;

            if ((!hasU && !hasV) || (nHitsU == 0 && nHitsV == 0) || dofUV <= std::numeric_limits<float>::epsilon())
                return -std::numeric_limits<float>::max();

            const auto likelihood = (weightU + weightV) / dofUV;
            if(likelihood > -1.f) return likelihood;
            else return -std::numeric_limits<float>::max();
        };

        const auto protonLikelihoodFwdUVW = GetMultiplaneBraggLikelihood(protonLikelihoodFwdU, protonLikelihoodFwdV, protonLikelihoodFwdW);
        const auto protonLikelihoodBwdUVW = GetMultiplaneBraggLikelihood(protonLikelihoodBwdU, protonLikelihoodBwdV, protonLikelihoodBwdW);
        const auto pionLikelihoodFwdUVW = GetMultiplaneBraggLikelihood(pionLikelihoodFwdU, pionLikelihoodFwdV, pionLikelihoodFwdW);
        const auto pionLikelihoodBwdUVW = GetMultiplaneBraggLikelihood(pionLikelihoodBwdU, pionLikelihoodBwdV, pionLikelihoodBwdW);
        const auto muonLikelihoodFwdUVW = GetMultiplaneBraggLikelihood(muonLikelihoodFwdU, muonLikelihoodFwdV, muonLikelihoodFwdW);
        const auto muonLikelihoodBwdUVW = GetMultiplaneBraggLikelihood(muonLikelihoodBwdU, muonLikelihoodBwdV, muonLikelihoodBwdW);
        const auto mipLikelihoodUVW = GetMultiplaneBraggLikelihood(mipLikelihoodU, mipLikelihoodV, mipLikelihoodW);

        if(protonLikelihoodFwdUVW != floatMinValue) pParticle->likelihoodForwardProton.Set(protonLikelihoodFwdUVW);
        if(protonLikelihoodBwdUVW != floatMinValue) pParticle->likelihoodBackwardProton.Set(protonLikelihoodBwdUVW);
        if(pionLikelihoodFwdUVW != floatMinValue) pParticle->likelihoodForwardPion.Set(pionLikelihoodFwdUVW);
        if(pionLikelihoodBwdUVW != floatMinValue) pParticle->likelihoodBackwardPion.Set(pionLikelihoodBwdUVW);
        if(muonLikelihoodFwdUVW != floatMinValue) pParticle->likelihoodForwardMuon.Set(muonLikelihoodFwdUVW);
        if(muonLikelihoodBwdUVW != floatMinValue) pParticle->likelihoodBackwardMuon.Set(muonLikelihoodBwdUVW);
        if(mipLikelihoodUVW != floatMinValue) pParticle->likelihoodMIP.Set(mipLikelihoodUVW);

        if(pParticlePeLEE->trk_trunk_rr_dEdx_u_v() != floatMinValue) pParticle->truncatedMeandEdxU.Set(pParticlePeLEE->trk_trunk_rr_dEdx_u_v());
        if(pParticlePeLEE->trk_trunk_rr_dEdx_v_v() != floatMinValue) pParticle->truncatedMeandEdxV.Set(pParticlePeLEE->trk_trunk_rr_dEdx_v_v());
        if(pParticlePeLEE->trk_trunk_rr_dEdx_y_v() != floatMinValue) pParticle->truncatedMeandEdxW.Set(pParticlePeLEE->trk_trunk_rr_dEdx_y_v());

        if(pParticle->yzAngle.IsSet())
        {
            // Values match default ubcc1pi ntuples
            const auto lengthFraction = 0.333333333f;
            const auto nHitsToSkip = 3;

            const bool isTrackAlongWWire = (std::pow(std::sin(pParticle->yzAngle()), 2) < sin2AngleThreshold);
            if (!isTrackAlongWWire && pParticle->truncatedMeandEdxW.IsSet())
            {
                pParticle->truncatedMeandEdx.Set(pParticle->truncatedMeandEdxW());
            }
            else
            {
                const auto uWeight = std::max(0.f, (pParticlePeLEE->trk_nhits_u_v() /*dedxPerHit.size()*/ * lengthFraction) - nHitsToSkip);
                const auto vWeight = std::max(0.f, (pParticlePeLEE->trk_nhits_v_v() /*dedxPerHit.size()*/ * lengthFraction) - nHitsToSkip);
                const auto hasU = pParticle->truncatedMeandEdxU.IsSet() && uWeight > 0.f;
                const auto hasV = pParticle->truncatedMeandEdxV.IsSet() && vWeight > 0.f;

                if (hasU || hasV)
                {
                    float truncatedMeandEdx = 0.f;
                    if (hasU) truncatedMeandEdx += pParticle->truncatedMeandEdxU() * uWeight;
                    if (hasV) truncatedMeandEdx += pParticle->truncatedMeandEdxV() * vWeight;

                    truncatedMeandEdx /= (uWeight + vWeight);
                    pParticle->truncatedMeandEdx.Set(truncatedMeandEdx);
                }
            }
        }

        // pParticle->truncatedMeandEdx.Set(pParticlePeLEE->);
        // pParticle->truthMatchPurities.Set(pParticlePeLEE->);
        // pParticle->truthMatchCompletenesses.Set(pParticlePeLEE->);
        // pParticle->hasMatchedMCParticle.Set(pParticlePeLEE->);
        // pParticle->bestMatchedMCParticleIndex.Set(pParticlePeLEE->);
        if(pParticlePeLEE->pfp_generation_v() > 0) pParticle->generation.Set(pParticlePeLEE->pfp_generation_v());
        if(pParticlePeLEE->trk_pid_chipr_v() >= 0) pParticle->chi2ForwardProtonW.Set(pParticlePeLEE->trk_pid_chipr_v());
        if(pParticlePeLEE->trk_pid_chimu_v() >= 0) pParticle->chi2ForwardMuonW.Set(pParticlePeLEE->trk_pid_chimu_v());
        if(pParticlePeLEE->trk_distance_v() != floatMinValue) pParticle->distance.Set(pParticlePeLEE->trk_distance_v());

        // if(pParticlePeLEE->dvtx_x_boundary() != doubleMinValue) pParticle->vertexDistanceToXBoundary.Set(pParticlePeLEE->dvtx_x_boundary());
        // if(pParticlePeLEE->dvtx_y_boundary() != doubleMinValue) pParticle->vertexDistanceToYBoundary.Set(pParticlePeLEE->dvtx_y_boundary());
        // if(pParticlePeLEE->dvtx_z_boundary() != doubleMinValue) pParticle->vertexDistanceToZBoundary.Set(pParticlePeLEE->dvtx_z_boundary());

        if(pParticlePeLEE->pfp_vtx_x_v() != floatMinValue) pParticle->vertexX.Set(pParticlePeLEE->pfp_vtx_x_v());
        if(pParticlePeLEE->pfp_vtx_y_v() != floatMinValue) pParticle->vertexY.Set(pParticlePeLEE->pfp_vtx_y_v());
        if(pParticlePeLEE->pfp_vtx_z_v() != floatMinValue) pParticle->vertexZ.Set(pParticlePeLEE->pfp_vtx_z_v());

        if(hasBacktrackedInfo)
        {
            const auto pBacktrackedParticlePeLEE = &eventPeLEE.backtracked.particles.at(i);
            if(pBacktrackedParticlePeLEE->backtracked_pdg() != 0) pParticle->pdgBacktracked.Set(pBacktrackedParticlePeLEE->backtracked_pdg());
            if(pBacktrackedParticlePeLEE->backtracked_purity() != floatMinValue) pParticle->purityBacktracked.Set(pBacktrackedParticlePeLEE->backtracked_purity());
            if(pBacktrackedParticlePeLEE->backtracked_completeness() != floatMinValue) pParticle->completenessBacktracked.Set(pBacktrackedParticlePeLEE->backtracked_completeness());
            if(pBacktrackedParticlePeLEE->backtracked_px() != floatMinValue) pParticle->momentumXBacktracked.Set(pBacktrackedParticlePeLEE->backtracked_px());
            if(pBacktrackedParticlePeLEE->backtracked_py() != floatMinValue) pParticle->momentumYBacktracked.Set(pBacktrackedParticlePeLEE->backtracked_py());
            if(pBacktrackedParticlePeLEE->backtracked_pz() != floatMinValue) pParticle->momentumZBacktracked.Set(pBacktrackedParticlePeLEE->backtracked_pz());
        }

        // if(pParticlePeLEE->pfpdg() == 11)
        // {
        //     if(pParticlePeLEE->pfp_vtx_x_v() != floatMinValue) pParticle->vertexX.Set(pParticlePeLEE->pfp_vtx_x_v());
        //     if(pParticlePeLEE->pfp_vtx_y_v() != floatMinValue) pParticle->vertexY.Set(pParticlePeLEE->pfp_vtx_y_v());
        //     if(pParticlePeLEE->pfp_vtx_z_v() != floatMinValue) pParticle->vertexZ.Set(pParticlePeLEE->pfp_vtx_z_v());

        //     // if(pParticlePeLEE->shr_start_x_v() != floatMinValue) pParticle->vertexX.Set(pParticlePeLEE->shr_start_x_v());
        //     // else if(pParticlePeLEE->trk_start_x_v() != floatMinValue) pParticle->vertexX.Set(pParticlePeLEE->trk_start_x_v());

        //     // if(pParticlePeLEE->shr_start_y_v() != floatMinValue) pParticle->vertexY.Set(pParticlePeLEE->shr_start_y_v());
        //     // else if(pParticlePeLEE->trk_start_y_v() != floatMinValue) pParticle->vertexY.Set(pParticlePeLEE->trk_start_y_v());

        //     // if(pParticlePeLEE->shr_start_z_v() != floatMinValue) pParticle->vertexZ.Set(pParticlePeLEE->shr_start_z_v());
        //     // else if(pParticlePeLEE->trk_start_z_v() != floatMinValue) pParticle->vertexZ.Set(pParticlePeLEE->trk_start_z_v());

        //     if(pParticlePeLEE->shr_start_x_v() != floatMinValue) std::cout<<"shr x: "<<pParticlePeLEE->shr_start_x_v()<<", ";
        //     if(pParticlePeLEE->shr_start_y_v() != floatMinValue) std::cout<<"shr y: "<<pParticlePeLEE->shr_start_y_v()<<", ";
        //     if(pParticlePeLEE->shr_start_z_v() != floatMinValue) std::cout<<"shr z: "<<pParticlePeLEE->shr_start_z_v()<<", ";
        //     if(pParticlePeLEE->shr_tkfit_start_x_v() != floatMinValue) std::cout<<"shr_tkfit x: "<<pParticlePeLEE->shr_tkfit_start_x_v()<<", ";
        //     if(pParticlePeLEE->shr_tkfit_start_y_v() != floatMinValue) std::cout<<"shr_tkfit y: "<<pParticlePeLEE->shr_tkfit_start_y_v()<<", ";
        //     if(pParticlePeLEE->shr_tkfit_start_z_v() != floatMinValue) std::cout<<"shr_tkfit z: "<<pParticlePeLEE->shr_tkfit_start_z_v()<<", ";
        //     if(pParticlePeLEE->trk_start_x_v() != floatMinValue) std::cout<<"trk x: "<<pParticlePeLEE->trk_start_x_v()<<", ";
        //     if(pParticlePeLEE->trk_start_y_v() != floatMinValue) std::cout<<"trk y: "<<pParticlePeLEE->trk_start_y_v()<<", ";
        //     if(pParticlePeLEE->trk_start_z_v() != floatMinValue) std::cout<<"trk z: "<<pParticlePeLEE->trk_start_z_v()<<", ";
        //     std::cout<<std::endl;
        // }
        // else if(pParticlePeLEE->pfpdg() == 13)
        // {
        //     if(pParticlePeLEE->trk_start_x_v() != floatMinValue) pParticle->vertexX.Set(pParticlePeLEE->trk_start_x_v());
        //     if(pParticlePeLEE->trk_start_y_v() != floatMinValue) pParticle->vertexY.Set(pParticlePeLEE->trk_start_y_v());
        //     if(pParticlePeLEE->trk_start_z_v() != floatMinValue) pParticle->vertexZ.Set(pParticlePeLEE->trk_start_z_v());
        // }
    }

}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::Print() const
{
    std::cout << std::string(80, '=') << std::endl;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "METADATA" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", metadata, UBCC1PI_MACRO_PRINT_MEMBER)

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "TRUTH" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", truth, UBCC1PI_MACRO_PRINT_MEMBER)

    for (const auto &particle : truth.particles)
    {
        std::cout << "TRUTH PARTICLE" << std::endl;
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", particle, UBCC1PI_MACRO_PRINT_MEMBER)
    }

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "RECO" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    UBCC1PI_MACRO_EVENT_RECO_MEMBERS("", reco, UBCC1PI_MACRO_PRINT_MEMBER)

    for (const auto &particle : reco.particles)
    {
        std::cout << "RECO PARTICLE" << std::endl;
        UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS("", particle, UBCC1PI_MACRO_PRINT_MEMBER)
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::BindToOutputTree(TTree * pTree)
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_MEMBERS(reco, reco, UBCC1PI_MACRO_BIND_OUTPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_BIND_OUTPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::BindToInputTree(TTree * pTree)
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS(metadata, metadata, UBCC1PI_MACRO_BIND_INPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS(truth, truth, UBCC1PI_MACRO_BIND_INPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_BIND_INPUT_VECTOR_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_MEMBERS(reco, reco, UBCC1PI_MACRO_BIND_INPUT_BRANCH)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_BIND_INPUT_VECTOR_BRANCH)
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::Reset()
{
    UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", metadata, UBCC1PI_MACRO_RESET_MEMBER)

    UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", truth, UBCC1PI_MACRO_RESET_MEMBER)
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_RESET_MEMBER_VECTOR)
    truth.particles.clear();

    UBCC1PI_MACRO_EVENT_RECO_MEMBERS("", reco, UBCC1PI_MACRO_RESET_MEMBER)
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_RESET_MEMBER_VECTOR)
    reco.particles.clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::PrepareForTreeFill()
{
    for (const auto &particle : truth.particles)
    {
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, particle, UBCC1PI_MACRO_FILL_MEMBER_VECTOR)
    }

    for (const auto &particle : reco.particles)
    {
        UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, particle, UBCC1PI_MACRO_FILL_MEMBER_VECTOR)
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

void Event::PrepareAfterTreeRead()
{
    unsigned int nTruthParticles;
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, &nTruthParticles, UBCC1PI_MACRO_GET_MEMBER_VECTOR_SIZE)

    truth.particles.resize(nTruthParticles);
    UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, truth.particles, UBCC1PI_MACRO_READ_MEMBER_VECTOR)

    unsigned int nRecoParticles;
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, &nRecoParticles, UBCC1PI_MACRO_GET_MEMBER_VECTOR_SIZE)

    reco.particles.resize(nRecoParticles);
    UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, reco.particles, UBCC1PI_MACRO_READ_MEMBER_VECTOR)
}


} // namespace ubcc1pi
