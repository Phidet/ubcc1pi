/**
 *  @file  ubcc1pi_standalone/Objects/TreeWriter.cxx
 *
 *  @brief  Code to facilitate setting output TTree branch addresses
 *
 */

#include "ubcc1pi_standalone/Objects/TreeWriter.h"

namespace ubcc1pi
{

TreeWriter::TreeWriter(const std::vector<std::string> inputFiles, const std::string outputFile)
  {
    // std::cout<<"DEBUG TreeWriter Point 0"<<std::endl;
    // OUTPUT TTREE
    // Make an output TTree for plotting (one entry per event)
    m_pOutFile = new TFile(outputFile.c_str(), "recreate");
    // std::cout<<"DEBUG TreeWriter Point 1"<<std::endl;
    m_pOutTree = new TTree( "stv_tree", "Ubcc1pi analysis tree" ); // Todo: maybe this can be changed to a custom name at some point when xsec_analyzer allows it
    // std::cout<<"DEBUG TreeWriter Point 2"<<std::endl;

    m_pEventChain = new TChain("nuselection/NeutrinoSelectionFilter");
    // std::cout<<"DEBUG TreeWriter Point 3"<<std::endl;
    m_pSubrunChain = new TChain("nuselection/SubRun");
    // std::cout<<"DEBUG TreeWriter Point 4"<<std::endl;
    for (const auto& fileName : inputFiles)
    {
        m_pEventChain->Add(fileName.c_str());
        m_pSubrunChain->Add(fileName.c_str());
    }
    // std::cout<<"DEBUG TreeWriter Point 5"<<std::endl;

    // Get the total POT from the subruns TTree. Save it in the output
    // TFile as a TParameter<float>. Real data doesn't have this TTree,
    // so check that it exists first.
    float pot;
    auto totalExposurePOT = 0.f;
    auto hasPotBranch = (m_pSubrunChain->GetBranch("pot") != nullptr);
    if (hasPotBranch) {
        m_pSubrunChain->SetBranchAddress("pot", &pot);
        for (int se = 0; se < m_pSubrunChain->GetEntries(); ++se) {
            m_pSubrunChain->GetEntry(se);
            totalExposurePOT += pot;
        }
    }
    // std::cout<<"DEBUG TreeWriter Point 6"<<std::endl;

    m_pTotalExposurePOTParam = new TParameter<float>("summed_pot", totalExposurePOT);
    m_pTotalExposurePOTParam->Write();
    // std::cout<<"DEBUG TreeWriter Point 7"<<std::endl;
  }

// Create destructor
TreeWriter::~TreeWriter()
{
    // std::cout<<"DEBUG ~TreeWriter Point -1"<<std::endl;
    m_pOutFile->cd();
    // std::cout<<"DEBUG ~TreeWriter Point 0"<<std::endl;
    m_pOutTree->Write();
    // std::cout<<"DEBUG ~TreeWriter Point 1"<<std::endl;
    m_pOutFile->Close();
    // std::cout<<"DEBUG ~TreeWriter Point 2"<<std::endl;
    // delete m_pOutTree;
    // // std::cout<<"DEBUG ~TreeWriter Point 2.1"<<std::endl;
    delete m_pOutFile;
    // std::cout<<"DEBUG ~TreeWriter Point 2.2"<<std::endl;
    delete m_pEventChain;
    // std::cout<<"DEBUG ~TreeWriter Point 2.2"<<std::endl;
    delete m_pSubrunChain;
    // std::cout<<"DEBUG ~TreeWriter Point 2.3"<<std::endl;
    delete m_pTotalExposurePOTParam;
    // std::cout<<"DEBUG ~TreeWriter Point 3"<<std::endl;
}

void TreeWriter::Fill()
{
    m_pOutTree->Fill();
}

// template <typename T>
// class TreeWriter::Pointer : public std::unique_ptr<T>
// {
//     public:
//         Pointer() : std::unique_ptr<T>( new T ) {}

//         T*& GetBarePtr()
//         {
//             m_pT = this->get();
//             return m_pT;
//         }

//     protected:
//         T* m_pT = nullptr;
// };

void TreeWriter::SetOutputBranchAddress(const std::string& branchName, void* address, const std::string& leafSpec)
{
    if ( !m_createdOutputBranches )
    {
        if ( leafSpec != "" ) m_pOutTree->Branch( branchName.c_str(), address, leafSpec.c_str() );
        else m_pOutTree->Branch( branchName.c_str(), address );
    }
    else m_pOutTree->SetBranchAddress( branchName.c_str(), address );
}

template <typename T>
void TreeWriter::SetObjectOutputBranchAddress(const std::string& branchName, const T*& address)
{
    if ( !m_createdOutputBranches ) m_pOutTree->Branch( branchName.c_str(), &address );
    else m_pOutTree->SetBranchAddress( branchName.c_str(), &address );
}

// template <typename T>
// void TreeWriter::SetObjectOutputBranchAddress(const std::string& branchName, T* address)
// {
//     if ( !m_createdOutputBranches ) m_pOutTree->Branch( branchName.c_str(), &address );
//     else m_pOutTree->SetBranchAddress( branchName.c_str(), &address );
// }

// template <typename T>
// void TreeWriter::SetObjectOutputBranchAddress(const std::string& branchName, TreeWriter::Pointer<T>& ptr)
// {
//   T*& address = ptr.getBarePtr();
//   TreeWriter::SetObjectOutputBranchAddress(branchName, address);
// }

void TreeWriter::CreateNoNewBranches()
{
    m_createdOutputBranches = true;
}

// Need to explicitly instantiate the template we use
template void TreeWriter::SetObjectOutputBranchAddress<std::vector<double>>(const std::string& branchName, const std::vector<double>*& address);
template void TreeWriter::SetObjectOutputBranchAddress<std::vector<int>>(const std::string& branchName, const std::vector<int>*& address);
template void TreeWriter::SetObjectOutputBranchAddress<std::vector<float>>(const std::string& branchName, const std::vector<float>*& address);


} // namespace ubcc1pi











//   // Signal definition flags
//   SetOutputBranchAddress("is_mc", &pEventPeLEE->is_mc_, "is_mc/O" );
//   SetOutputBranchAddress("mc_neutrino_is_numu", &pEventPeLEE->mc_neutrino_is_numu_, "mc_neutrino_is_numu/O" );
//   SetOutputBranchAddress("mc_vertex_in_FV", &pEventPeLEE->mc_vertex_in_FV_, "mc_vertex_in_FV/O" );
//   SetOutputBranchAddress("mc_muon_in_mom_range", &pEventPeLEE->mc_muon_in_mom_range_, "mc_muon_in_mom_range/O" );
//   SetOutputBranchAddress("mc_lead_p_in_mom_range", &pEventPeLEE->mc_lead_p_in_mom_range_, "mc_lead_p_in_mom_range/O" );
//   SetOutputBranchAddress("mc_no_fs_pi0", &pEventPeLEE->mc_no_fs_pi0_, "mc_no_fs_pi0/O" );
//   SetOutputBranchAddress("mc_no_charged_pi_above_threshold", &pEventPeLEE->mc_no_charged_pi_above_threshold_, "mc_no_charged_pi_above_threshold/O" );
//   SetOutputBranchAddress("mc_no_fs_mesons", &pEventPeLEE->mc_no_fs_mesons_, "mc_no_fs_mesons/O" );
//   SetOutputBranchAddress("mc_is_signal", &pEventPeLEE->mc_is_signal_, "mc_is_signal/O" );

//   // MC event category
//   SetOutputBranchAddress("category", &pEventPeLEE->category_, "category/I" );

//   // Event weights
//   SetOutputBranchAddress("spline_weight", &pEventPeLEE->spline_weight_, "spline_weight/F" );
//   SetOutputBranchAddress("tuned_cv_weight", &pEventPeLEE->tuned_cv_weight_, "tuned_cv_weight/F" );

//   // If MC weights are available, prepare to store them in the output TTree
//   if ( pEventPeLEE->mc_weights_map_ ) {

//     // Make separate branches for the various sets of systematic variation
//     // weights in the map
//     for ( auto& pair : *pEventPeLEE->mc_weights_map_ ) {

//       // Prepend "weight_" to the name of the vector of weights in the map
//       std::string weight_branch_name = "weight_" + pair.first;

//       // Store a pointer to the vector of weights (needed to set the branch
//       // address properly) in the temporary map of pointers
//       pEventPeLEE->mc_weights_ptr_map_[ weight_branch_name ] = &pair.second;

//       // Set the branch address for this vector of weights
//       SetObjectOutputBranchAddress< std::vector<double> >(weight_branch_name, pEventPeLEE->mc_weights_ptr_map_.at(weight_branch_name), m_createdOutputBranches);
//     }
//   }

//   // Backtracked neutrino purity and completeness
//   SetOutputBranchAddress("nu_completeness_from_pfp", &pEventPeLEE->nu_completeness_from_pfp_, "nu_completeness_from_pfp/F" );
//   SetOutputBranchAddress("nu_purity_from_pfp", &pEventPeLEE->nu_purity_from_pfp_, "nu_purity_from_pfp/F" );

//   // Number of neutrino slices identified by the SliceID
//   SetOutputBranchAddress("nslice", &pEventPeLEE->nslice_, "nslice/I" );

//   // CCNp0pi selection criteria
//   SetOutputBranchAddress("sel_nu_mu_cc", &pEventPeLEE->sel_nu_mu_cc_, "sel_nu_mu_cc/O" );
//   SetOutputBranchAddress("sel_reco_vertex_in_FV", &pEventPeLEE->sel_reco_vertex_in_FV_, "sel_reco_vertex_in_FV/O" );
//   SetOutputBranchAddress("sel_topo_cut_passed", &pEventPeLEE->sel_topo_cut_passed_, "sel_topo_cut_passed/O" );
//   SetOutputBranchAddress("sel_cosmic_ip_cut_passed", &pEventPeLEE->sel_cosmic_ip_cut_passed_, "sel_cosmic_ip_cut_passed/O" );
//   SetOutputBranchAddress("sel_pfp_starts_in_PCV", &pEventPeLEE->sel_pfp_starts_in_PCV_, "sel_pfp_starts_in_PCV/O" );
//   SetOutputBranchAddress("sel_no_reco_showers", &pEventPeLEE->sel_no_reco_showers_, "sel_no_reco_showers/O" );
//   SetOutputBranchAddress("sel_has_muon_candidate", &pEventPeLEE->sel_has_muon_candidate_, "sel_has_muon_candidate/O" );
//   SetOutputBranchAddress("sel_muon_contained", &pEventPeLEE->sel_muon_contained_, "sel_muon_contained/O" );
//   SetOutputBranchAddress("sel_muon_passed_mom_cuts", &pEventPeLEE->sel_muon_passed_mom_cuts_, "sel_muon_passed_mom_cuts/O" );
//   SetOutputBranchAddress("sel_muon_quality_ok", &pEventPeLEE->sel_muon_quality_ok_, "sel_muon_quality_ok/O" );
//   SetOutputBranchAddress("sel_has_p_candidate", &pEventPeLEE->sel_has_p_candidate_, "sel_has_p_candidate/O" );
//   SetOutputBranchAddress("sel_passed_proton_pid_cut", &pEventPeLEE->sel_passed_proton_pid_cut_, "sel_passed_proton_pid_cut/O" );
//   SetOutputBranchAddress("sel_protons_contained", &pEventPeLEE->sel_protons_contained_, "sel_protons_contained/O" );
//   SetOutputBranchAddress("sel_lead_p_passed_mom_cuts", &pEventPeLEE->sel_lead_p_passed_mom_cuts_, "sel_lead_p_passed_mom_cuts/O" );
//   SetOutputBranchAddress("sel_CCNp0pi", &pEventPeLEE->sel_CCNp0pi_, "sel_CCNp0pi/O" );

//   // Index for the muon candidate in the vectors of PFParticles
//   SetOutputBranchAddress("muon_candidate_idx", &pEventPeLEE->muon_candidate_idx_, "muon_candidate_idx/I" );

//   // Index for the leading proton candidate in the vectors of PFParticles
//   SetOutputBranchAddress("lead_p_candidate_idx", &pEventPeLEE->lead_p_candidate_idx_, "lead_p_candidate_idx/I" );

//   // Reco 3-momenta (muon, leading proton)
//   SetObjectOutputBranchAddress< TVector3 >("p3_mu", pEventPeLEE->p3_mu_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< TVector3 >("p3_lead_p", pEventPeLEE->p3_lead_p_, m_createdOutputBranches);

//   // Reco 3-momenta (all proton candidates, ordered from highest to lowest
//   // magnitude)
//   SetObjectOutputBranchAddress< std::vector<TVector3> >("p3_p_vec", pEventPeLEE->p3_p_vec_, m_createdOutputBranches);

//   // True 3-momenta (muon, leading proton)
//   SetObjectOutputBranchAddress< TVector3 >("mc_p3_mu", pEventPeLEE->mc_p3_mu_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< TVector3 >("mc_p3_lead_p", pEventPeLEE->mc_p3_lead_p_, m_createdOutputBranches);

//   // True 3-momenta (all protons, ordered from highest to lowest magnitude)
//   SetObjectOutputBranchAddress< std::vector<TVector3> >("mc_p3_p_vec", pEventPeLEE->mc_p3_p_vec_, m_createdOutputBranches);

// //   // Reco STVs
// //   SetOutputBranchAddress("delta_pT", &pEventPeLEE->delta_pT_, "delta_pT/F" );
// //   SetOutputBranchAddress("delta_phiT", &pEventPeLEE->delta_phiT_, "delta_phiT/F" );
// //   SetOutputBranchAddress("delta_alphaT", &pEventPeLEE->delta_alphaT_, "delta_alphaT/F" );
// //   SetOutputBranchAddress("delta_pL", &pEventPeLEE->delta_pL_, "delta_pL/F" );
// //   SetOutputBranchAddress("pn", &pEventPeLEE->pn_, "pn/F" );
// //   SetOutputBranchAddress("delta_pTx", &pEventPeLEE->delta_pTx_, "delta_pTx/F" );
// //   SetOutputBranchAddress("delta_pTy", &pEventPeLEE->delta_pTy_, "delta_pTy/F" );
// //   SetOutputBranchAddress("theta_mu_p", &pEventPeLEE->theta_mu_p_, "theta_mu_p/F" );

// //   // MC STVs (only filled for signal events)
// //   SetOutputBranchAddress("mc_delta_pT", &pEventPeLEE->mc_delta_pT_, "mc_delta_pT/F" );
// //   SetOutputBranchAddress("mc_delta_phiT", &pEventPeLEE->mc_delta_phiT_, "mc_delta_phiT/F" );
// //   SetOutputBranchAddress("mc_delta_alphaT", &pEventPeLEE->mc_delta_alphaT_, "mc_delta_alphaT/F" );
// //   SetOutputBranchAddress("mc_delta_pL", &pEventPeLEE->mc_delta_pL_, "mc_delta_pL/F" );
// //   SetOutputBranchAddress("mc_pn", &pEventPeLEE->mc_pn_, "mc_pn/F" );
// //   SetOutputBranchAddress("mc_delta_pTx", &pEventPeLEE->mc_delta_pTx_, "mc_delta_pTx/F" );
// //   SetOutputBranchAddress("mc_delta_pTy", &pEventPeLEE->mc_delta_pTy_, "mc_delta_pTy/F" );
// //   SetOutputBranchAddress("mc_theta_mu_p", &pEventPeLEE->mc_theta_mu_p_, "mc_theta_mu_p/F" );

//   // *** Branches copied directly from the input ***

//   // Cosmic rejection parameters for numu CC inclusive selection
//   SetOutputBranchAddress("topological_score", &pEventPeLEE->topological_score_, "topological_score/F" );
//   SetOutputBranchAddress("CosmicIP", &pEventPeLEE->cosmic_impact_parameter_, "CosmicIP/F" );

//   // Reconstructed neutrino vertex position
//   SetOutputBranchAddress("reco_nu_vtx_sce_x", &pEventPeLEE->nu_vx_, "reco_nu_vtx_sce_x/F" );
//   SetOutputBranchAddress("reco_nu_vtx_sce_y", &pEventPeLEE->nu_vy_, "reco_nu_vtx_sce_y/F" );
//   SetOutputBranchAddress("reco_nu_vtx_sce_z", &pEventPeLEE->nu_vz_, "reco_nu_vtx_sce_z/F" );

//   // MC truth information for the neutrino
//   SetOutputBranchAddress("mc_nu_pdg", &pEventPeLEE->mc_nu_pdg_, "mc_nu_pdg/I" );
//   SetOutputBranchAddress("mc_nu_vtx_x", &pEventPeLEE->mc_nu_vx_, "mc_nu_vtx_x/F" );
//   SetOutputBranchAddress("mc_nu_vtx_y", &pEventPeLEE->mc_nu_vy_, "mc_nu_vtx_y/F" );
//   SetOutputBranchAddress("mc_nu_vtx_z", &pEventPeLEE->mc_nu_vz_, "mc_nu_vtx_z/F" );
//   SetOutputBranchAddress("mc_nu_energy", &pEventPeLEE->mc_nu_energy_, "mc_nu_energy/F" );
//   SetOutputBranchAddress("mc_ccnc", &pEventPeLEE->mc_nu_ccnc_, "mc_ccnc/I" );
//   SetOutputBranchAddress("mc_interaction", &pEventPeLEE->mc_nu_interaction_type_, "mc_interaction/I" );

//   // PFParticle properties
//   SetObjectOutputBranchAddress< std::vector<unsigned int> >("pfp_generation_v", pEventPeLEE->pfp_generation_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<unsigned int> >("pfp_trk_daughters_v", pEventPeLEE->pfp_trk_daughters_count_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<unsigned int> >("pfp_shr_daughters_v", pEventPeLEE->pfp_shr_daughters_count_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_score_v", pEventPeLEE->pfp_track_score_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<int> >("pfpdg", pEventPeLEE->pfp_reco_pdg_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<int> >("pfnhits", pEventPeLEE->pfp_hits_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<int> >("pfnplanehits_U", pEventPeLEE->pfp_hitsU_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<int> >("pfnplanehits_V", pEventPeLEE->pfp_hitsV_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<int> >("pfnplanehits_Y", pEventPeLEE->pfp_hitsY_, m_createdOutputBranches);

//   // Backtracked PFParticle properties
//   SetObjectOutputBranchAddress< std::vector<int> >("backtracked_pdg", pEventPeLEE->pfp_true_pdg_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("backtracked_e", pEventPeLEE->pfp_true_E_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("backtracked_px", pEventPeLEE->pfp_true_px_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("backtracked_py", pEventPeLEE->pfp_true_py_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("backtracked_pz", pEventPeLEE->pfp_true_pz_, m_createdOutputBranches);

//   // Shower properties
//   // For some ntuples, reconstructed shower information is excluded.
//   // In such cases, skip writing these branches to the output TTree.
//   if ( pEventPeLEE->shower_startx_ ) {
//     SetObjectOutputBranchAddress< std::vector<float> >("shr_start_x_v", pEventPeLEE->shower_startx_, m_createdOutputBranches);
//     SetObjectOutputBranchAddress< std::vector<float> >("shr_start_y_v", pEventPeLEE->shower_starty_, m_createdOutputBranches);
//     SetObjectOutputBranchAddress< std::vector<float> >("shr_start_z_v", pEventPeLEE->shower_startz_, m_createdOutputBranches);
//     // Shower start distance from reco neutrino vertex (pre-calculated for
//     // convenience)
//     SetObjectOutputBranchAddress< std::vector<float> >("shr_dist_v", pEventPeLEE->shower_start_distance_, m_createdOutputBranches);
//   }

//   // Track properties
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_len_v", pEventPeLEE->track_length_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_sce_start_x_v", pEventPeLEE->track_startx_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_sce_start_y_v", pEventPeLEE->track_starty_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_sce_start_z_v", pEventPeLEE->track_startz_, m_createdOutputBranches);

//   // Track start distance from reco neutrino vertex (pre-calculated for
//   // convenience)
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_distance_v", pEventPeLEE->track_start_distance_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_sce_end_x_v", pEventPeLEE->track_endx_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_sce_end_y_v", pEventPeLEE->track_endy_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_sce_end_z_v", pEventPeLEE->track_endz_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_dir_x_v", pEventPeLEE->track_dirx_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_dir_y_v", pEventPeLEE->track_diry_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_dir_z_v", pEventPeLEE->track_dirz_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_energy_proton_v", pEventPeLEE->track_kinetic_energy_p_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_range_muon_mom_v", pEventPeLEE->track_range_mom_mu_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("trk_mcs_muon_mom_v", pEventPeLEE->track_mcs_mom_mu_, m_createdOutputBranches);

//   // Some ntuples exclude the old chi^2 proton PID score. Only include it in
//   // the output if it is available.
//   if ( pEventPeLEE->track_chi2_proton_ )
//   {
//     SetObjectOutputBranchAddress< std::vector<float> >("trk_pid_chipr_v", pEventPeLEE->track_chi2_proton_, m_createdOutputBranches);
//   }

// //   // Log-likelihood-based particle ID information
// //   SetObjectOutputBranchAddress< std::vector<float> >("trk_llr_pid_v", pEventPeLEE->track_llr_pid_, m_createdOutputBranches);
// //   SetObjectOutputBranchAddress< std::vector<float> >("trk_llr_pid_u_v", pEventPeLEE->track_llr_pid_U_, m_createdOutputBranches);
// //   SetObjectOutputBranchAddress< std::vector<float> >("trk_llr_pid_v_v", pEventPeLEE->track_llr_pid_V_, m_createdOutputBranches);
// //   SetObjectOutputBranchAddress< std::vector<float> >("trk_llr_pid_y_v", pEventPeLEE->track_llr_pid_Y_, m_createdOutputBranches);
// //   SetObjectOutputBranchAddress< std::vector<float> >("trk_llr_pid_score_v", pEventPeLEE->track_llr_pid_score_, m_createdOutputBranches);

//   // MC truth information for the final-state primary particles
//   SetObjectOutputBranchAddress< std::vector<int> >("mc_pdg", pEventPeLEE->mc_nu_daughter_pdg_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("mc_E", pEventPeLEE->mc_nu_daughter_energy_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("mc_px", pEventPeLEE->mc_nu_daughter_px_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("mc_py", pEventPeLEE->mc_nu_daughter_py_, m_createdOutputBranches);
//   SetObjectOutputBranchAddress< std::vector<float> >("mc_pz", pEventPeLEE->mc_nu_daughter_pz_, m_createdOutputBranches);