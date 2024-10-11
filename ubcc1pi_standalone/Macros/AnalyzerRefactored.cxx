/**
 *  @file  ubcc1pi_standalone/Macros/AnalyzerRefactored.cxx
 *
 *  @brief The implementation file of the AnalyzerRefactored macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"
#include "ubcc1pi_standalone/Objects/FileReader.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"
#include "ubcc1pi_standalone/Helpers/SelectionHelper.h"
#include "ubcc1pi_standalone/Helpers/ExtractionHelper.h"
#include "ubcc1pi_standalone/Objects/TreeWriter.h"


using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void AnalyzerRefactored(const Config &config)
{

    // -------------------------------------------------------------------------------------------------------------------------------------
    // A map from each cross-section to the limits of the phase-space that should be considered. The key is the
    // identifier for the kinematic quantity and the mapped value is a pair containing the limits [min, max]
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::map< std::string, std::pair<float, float> > phaseSpaceMap;
    for (const auto &[name, binning] : std::vector< std::tuple<std::string, Config::Global::Binning> > {
        { "muonCosTheta", config.global.muonCosTheta},
        { "muonPhi", config.global.muonPhi},
        { "muonMomentum", config.global.muonMomentum},
        { "pionCosTheta", config.global.pionCosTheta},
        { "pionPhi", config.global.pionPhi},
        { "pionMomentum", config.global.pionMomentum},
        { "muonPionAngle", config.global.muonPionAngle},
        { "nProtons", config.global.nProtons}
    })
    {
        // Add to the phase-space map
        phaseSpaceMap.emplace(name, std::pair<float, float>({binning.min, binning.max}));
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the relevent "getters" for each cross-section and for the CC0pi selection
    // -------------------------------------------------------------------------------------------------------------------------------------
    ExtractionHelper::AnalysisValueMap getValueCC1Pi;
    ExtractionHelper::AnalysisValueMap getValueCC0Pi;
    ExtractionHelper::PopulateAnalysisValueMap(getValueCC1Pi, false);
    ExtractionHelper::PopulateAnalysisValueMap(getValueCC0Pi, true);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the selections
    // The CC0pi selection is modified to not execute the BDT cuts on the proton kinematics
    // -------------------------------------------------------------------------------------------------------------------------------------
    auto selectionCC0Pi = SelectionHelper::GetCC0piSelectionModifiedPeLEE(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), true); // Do not cut on proton kinematics here
    auto selectionCC1Pi = SelectionHelper::GetDefaultSelection2(true); // todo decide on final selection !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // ExtractionHelper::InputFileList inputData;
    // typedef std::vector< std::tuple<AnalysisHelper::SampleType, std::string, std::string, float> > InputFileList;
    // inputData.emplace_back(AnalysisHelper::Overlay, "", "/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root", 1);
    // inputData.emplace_back(AnalysisHelper::Overlay, "", "/uboone/app/users/jdetje/searchingfornues/files/steps4/neutrinoselection_filt_upodated5.root", 1);
    // inputData.emplace_back(AnalysisHelper::DataBNB, "", "/uboone/data/users/jdetje/pelee_v08_00_00_70/bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2_E1_head5.root", 1);

    // Get the muon BDTs
    const auto muonFeatureNames = BDTHelper::MuonBDTFeatureNames;
    BDTHelper::BDT muonBDT("muon", muonFeatureNames);
    const auto protonFeatureNames = BDTHelper::ProtonBDTFeatureNames;
    BDTHelper::BDT protonBDT("proton", protonFeatureNames);

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Get list of good runs
    // -------------------------------------------------------------------------------------------------------------------------------------
    std::vector<int> goodRuns;
    std::ifstream file(config.global.goodRunListFile);
    int run;
    while (file >> run) {
        goodRuns.push_back(run);
    }

    // -------------------------------------------------------------------------------------------------------------------------------------
    // Setup the output file
    // -------------------------------------------------------------------------------------------------------------------------------------
    // std::vector<std::string> inputFiles;
    // inputFiles.push_back("/uboone/data/users/jdetje/ubcc1piVSpelee/pelee/neutrinoselection_filt_0_4k.root");
    // TreeWriter treeWriter(inputFiles, config.global.outputFile.c_str());

    // Loop over the files
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        // Check if we should use this file
        if(!useThisFile) continue;


        // -------------------------------------------------------------------------------------------------------------------------------------
        // Setup the output ntuple file
        // -------------------------------------------------------------------------------------------------------------------------------------
        // Find the last occurrence of '/'
        size_t lastSlash = filePath.find_last_of('/');
        // Find the last occurrence of '.'
        size_t lastDot = filePath.find_last_of('.');
        // Extract the substring between the last '/' and the last '.'fileList
        std::string fileName = filePath.substr(lastSlash + 1, lastDot - lastSlash - 1);
        std::cout << "Reading input file: " << fileName << std::endl;
        const auto outputFilePath = config.global.outputPath + fileName + "_ubcc1pi.root";
        TreeWriter treeWriter(std::vector<std::string>({filePath}), outputFilePath.c_str());


        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // File-level variables
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        const auto isOverlay = (sampleType == AnalysisHelper::Overlay);
        const auto isDirt    = (sampleType == AnalysisHelper::Dirt);
        const auto isNuWro   = (sampleType == AnalysisHelper::NuWro);
        const auto isDataBNB = (sampleType == AnalysisHelper::DataBNB);
        const auto isDetVar  = (sampleType == AnalysisHelper::DetectorVariation);
        const auto isDataEXT = (sampleType == AnalysisHelper::DataEXT);
        const auto isMC = (sampleType != AnalysisHelper::DataBNB) && (sampleType != AnalysisHelper::DataEXT);

        treeWriter.SetOutputBranchAddress("is_mc", (void*)&isMC, "is_mc/O");
        treeWriter.SetOutputBranchAddress("is_Overlay", (void*)&isOverlay, "is_Overlay/O");
        treeWriter.SetOutputBranchAddress("is_Dirt", (void*)&isDirt, "is_Dirt/O");
        treeWriter.SetOutputBranchAddress("is_NuWro", (void*)&isNuWro, "is_NuWro/O");
        treeWriter.SetOutputBranchAddress("is_DataBNB", (void*)&isDataBNB, "is_DataBNB/O");
        treeWriter.SetOutputBranchAddress("is_DetVar", (void*)&isDetVar, "is_DetVar/O");
        treeWriter.SetOutputBranchAddress("is_DataEXT", (void*)&isDataEXT, "is_DataEXT/O");
        treeWriter.SetOutputBranchAddress("run", (void*)&fileRun, "run/I");

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Reco variables
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bool passesGenericCC0Pi, passesGoldenCC0Pi, passesGenericCC1Pi, passesGoldenCC1Pi;
        treeWriter.SetOutputBranchAddress("cc0pi_passes_generic", (void*)&passesGenericCC0Pi, "cc0pi_passes_generic/O" );
        treeWriter.SetOutputBranchAddress("cc0pi_passes_golden", (void*)&passesGoldenCC0Pi, "cc0pi_passes_golden/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_passes_generic", (void*)&passesGenericCC1Pi, "cc1pi_passes_generic/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_passes_golden", (void*)&passesGoldenCC1Pi, "cc1pi_passes_golden/O" );

        // Includes phace space restriction
        bool isSelectedGenericCC0Pi, isSelectedGoldenCC0Pi, isSelectedGenericCC1Pi, isSelectedGoldenCC1Pi;
        treeWriter.SetOutputBranchAddress("cc0pi_selected_generic", (void*)&isSelectedGenericCC0Pi, "cc0pi_selected_generic/O" );
        treeWriter.SetOutputBranchAddress("cc0pi_selected_golden", (void*)&isSelectedGoldenCC0Pi, "cc0pi_selected_golden/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_selected_generic", (void*)&isSelectedGenericCC1Pi, "cc1pi_selected_generic/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_selected_golden", (void*)&isSelectedGoldenCC1Pi, "cc1pi_selected_golden/O" );

        float cc1pi_reco_muonMomentum, cc1pi_reco_muonCosTheta, cc1pi_reco_muonPhi, cc1pi_reco_pionMomentum, cc1pi_reco_pionCosTheta, cc1pi_reco_pionPhi, cc1pi_reco_muonPionAngle;
        int   cc1pi_reco_nProtons;
        treeWriter.SetOutputBranchAddress("cc1pi_reco_muonMomentum", (void*)&cc1pi_reco_muonMomentum, "cc1pi_reco_muonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_reco_muonCosTheta", (void*)&cc1pi_reco_muonCosTheta, "cc1pi_reco_muonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_reco_muonPhi", (void*)&cc1pi_reco_muonPhi, "cc1pi_reco_muonPhi/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_reco_pionMomentum", (void*)&cc1pi_reco_pionMomentum, "cc1pi_reco_pionMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_reco_pionCosTheta", (void*)&cc1pi_reco_pionCosTheta, "cc1pi_reco_pionCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_reco_pionPhi", (void*)&cc1pi_reco_pionPhi, "cc1pi_reco_pionPhi/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_reco_muonPionAngle", (void*)&cc1pi_reco_muonPionAngle, "cc1pi_reco_muonPionAngle/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_reco_nProtons", (void*)&cc1pi_reco_nProtons, "cc1pi_reco_nProtons/I" );

        bool cc1pi_recoMuon_IsContained, cc1pi_recoPion_IsContained;
        float cc1pi_recoMuon_TrackLength, cc1pi_recoMuon_protonBDTScore, cc1pi_recoMuon_muonBDTScore, cc1pi_recoPion_TrackLength, cc1pi_recoPion_protonBDTScore, cc1pi_recoPion_muonBDTScore;
        treeWriter.SetOutputBranchAddress("cc1pi_recoMuon_IsContained", (void*)&cc1pi_recoMuon_IsContained, "cc1pi_recoMuon_IsContained/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_recoMuon_TrackLength", (void*)&cc1pi_recoMuon_TrackLength, "cc1pi_recoMuon_TrackLength/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_recoMuon_protonBDTScore", (void*)&cc1pi_recoMuon_protonBDTScore, "cc1pi_recoMuon_protonBDTScore/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_recoMuon_muonBDTScore", (void*)&cc1pi_recoMuon_muonBDTScore, "cc1pi_recoMuon_muonBDTScore/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_recoPion_IsContained", (void*)&cc1pi_recoPion_IsContained, "cc1pi_recoPion_IsContained/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_recoPion_TrackLength", (void*)&cc1pi_recoPion_TrackLength, "cc1pi_recoPion_TrackLength/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_recoPion_protonBDTScore", (void*)&cc1pi_recoPion_protonBDTScore, "cc1pi_recoPion_protonBDTScore/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_recoPion_muonBDTScore", (void*)&cc1pi_recoPion_muonBDTScore, "cc1pi_recoPion_muonBDTScore/F" );

        int   cc1pi_backtracked_protonPDG;
        float cc1pi_backtracked_protonMomentum, cc1pi_backtracked_protonCosTheta, cc1pi_backtracked_protonPhi;
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_protonPDG", (void*)&cc1pi_backtracked_protonPDG, "cc1pi_backtracked_protonPDG/I" );
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_protonMomentum", (void*)&cc1pi_backtracked_protonMomentum, "cc1pi_backtracked_protonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_protonCosTheta", (void*)&cc1pi_backtracked_protonCosTheta, "cc1pi_backtracked_protonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_protonPhi", (void*)&cc1pi_backtracked_protonPhi, "cc1pi_backtracked_protonPhi/F" );

        int   cc1pi_backtracked_muonPDG;
        float cc1pi_backtracked_muonMomentum, cc1pi_backtracked_muonCosTheta, cc1pi_backtracked_muonPhi;
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_muonPDG", (void*)&cc1pi_backtracked_muonPDG, "cc1pi_backtracked_muonPDG/I" );
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_muonMomentum", (void*)&cc1pi_backtracked_muonMomentum, "cc1pi_backtracked_muonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_muonCosTheta", (void*)&cc1pi_backtracked_muonCosTheta, "cc1pi_backtracked_muonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_backtracked_muonPhi", (void*)&cc1pi_backtracked_muonPhi, "cc1pi_backtracked_muonPhi/F" );

        float cc0pi_reco_muonMomentum, cc0pi_reco_muonCosTheta, cc0pi_reco_muonPhi, cc0pi_reco_protonMomentum, cc0pi_reco_protonCosTheta, cc0pi_reco_protonPhi, cc0pi_reco_muonProtonAngle;
        int   cc0pi_reco_nProtons;
        treeWriter.SetOutputBranchAddress("cc0pi_reco_muonMomentum", (void*)&cc0pi_reco_muonMomentum, "cc0pi_reco_muonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_reco_muonCosTheta", (void*)&cc0pi_reco_muonCosTheta, "cc0pi_reco_muonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_reco_muonPhi", (void*)&cc0pi_reco_muonPhi, "cc0pi_reco_muonPhi/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_reco_protonMomentum", (void*)&cc0pi_reco_protonMomentum, "cc0pi_reco_protonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_reco_protonCosTheta", (void*)&cc0pi_reco_protonCosTheta, "cc0pi_reco_protonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_reco_protonPhi", (void*)&cc0pi_reco_protonPhi, "cc0pi_reco_protonPhi/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_reco_muonProtonAngle", (void*)&cc0pi_reco_muonProtonAngle, "cc0pi_reco_muonProtonAngle/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_reco_nProtons", (void*)&cc0pi_reco_nProtons, "cc0pi_reco_nProtons/I" );
        int cc0pi_backtracked_protonPDG;
        float cc0pi_backtracked_protonMomentum, cc0pi_backtracked_protonCosTheta, cc0pi_backtracked_protonPhi, cc0pi_leadingProton_protonBDTScore, cc0pi_leadingProton_muonBDTScore;
        treeWriter.SetOutputBranchAddress("cc0pi_backtracked_protonPDG", (void*)&cc0pi_backtracked_protonPDG, "cc0pi_backtracked_protonPDG/I" );
        treeWriter.SetOutputBranchAddress("cc0pi_backtracked_protonMomentum", (void*)&cc0pi_backtracked_protonMomentum, "cc0pi_backtracked_protonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_backtracked_protonCosTheta", (void*)&cc0pi_backtracked_protonCosTheta, "cc0pi_backtracked_protonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_backtracked_protonPhi", (void*)&cc0pi_backtracked_protonPhi, "cc0pi_backtracked_protonPhi/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_leadingProton_protonBDTScore", (void*)&cc0pi_leadingProton_protonBDTScore, "cc0pi_leadingProton_protonBDTScore/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_leadingProton_muonBDTScore", (void*)&cc0pi_leadingProton_muonBDTScore, "cc0pi_leadingProton_muonBDTScore/F" );

        bool passed_particleTrackScore, passed_particleVertexDistance, passed_particleGeneration, passed_particleTrackLength, passed_particleProtonChi2, passed_particleMuonChi2;
        bool passed_particleProtonChi2OverMuonChi2, passed_pandoraNuPDGIsNumu, passed_daughterVerticesContained, passed_nuVertexFiducial, passed_topologicalOrFlashMatch;
        bool passed_topologicalScoreCC, passed_min2Tracks, passed_max1Uncontained, passed_2NonProtons, passed_pionHasValiddEdx, passed_pionNotInGap, passed_muonNotInGap;
        bool passed_topologicalScore, passed_startNearVertex, passed_openingAngle, passed_likelyGoldenPion;
        treeWriter.SetOutputBranchAddress("passed_particleTrackScore", (void*)&passed_particleTrackScore, "passed_particleTrackScore/O" );
        treeWriter.SetOutputBranchAddress("passed_particleVertexDistance", (void*)&passed_particleVertexDistance, "passed_particleVertexDistance/O" );
        treeWriter.SetOutputBranchAddress("passed_particleGeneration", (void*)&passed_particleGeneration, "passed_particleGeneration/O" );
        treeWriter.SetOutputBranchAddress("passed_particleTrackLength", (void*)&passed_particleTrackLength, "passed_particleTrackLength/O" );
        treeWriter.SetOutputBranchAddress("passed_particleProtonChi2", (void*)&passed_particleProtonChi2, "passed_particleProtonChi2/O" );
        treeWriter.SetOutputBranchAddress("passed_particleMuonChi2", (void*)&passed_particleMuonChi2, "passed_particleMuonChi2/O" );
        treeWriter.SetOutputBranchAddress("passed_particleProtonChi2OverMuonChi2", (void*)&passed_particleProtonChi2OverMuonChi2, "passed_particleProtonChi2OverMuonChi2/O" );
        treeWriter.SetOutputBranchAddress("passed_pandoraNuPDGIsNumu", (void*)&passed_pandoraNuPDGIsNumu, "passed_pandoraNuPDGIsNumu/O" );
        treeWriter.SetOutputBranchAddress("passed_daughterVerticesContained", (void*)&passed_daughterVerticesContained, "passed_daughterVerticesContained/O" );
        treeWriter.SetOutputBranchAddress("passed_nuVertexFiducial", (void*)&passed_nuVertexFiducial, "passed_nuVertexFiducial/O" );
        treeWriter.SetOutputBranchAddress("passed_topologicalOrFlashMatch", (void*)&passed_topologicalOrFlashMatch, "passed_topologicalOrFlashMatch/O" );
        treeWriter.SetOutputBranchAddress("passed_topologicalScoreCC", (void*)&passed_topologicalScoreCC, "passed_topologicalScoreCC/O" );
        treeWriter.SetOutputBranchAddress("passed_min2Tracks", (void*)&passed_min2Tracks, "passed_min2Tracks/O" );
        treeWriter.SetOutputBranchAddress("passed_max1Uncontained", (void*)&passed_max1Uncontained, "passed_max1Uncontained/O" );
        treeWriter.SetOutputBranchAddress("passed_2NonProtons", (void*)&passed_2NonProtons, "passed_2NonProtons/O" );
        treeWriter.SetOutputBranchAddress("passed_pionHasValiddEdx", (void*)&passed_pionHasValiddEdx, "passed_pionHasValiddEdx/O" );
        treeWriter.SetOutputBranchAddress("passed_pionNotInGap", (void*)&passed_pionNotInGap, "passed_pionNotInGap/O" );
        treeWriter.SetOutputBranchAddress("passed_muonNotInGap", (void*)&passed_muonNotInGap, "passed_muonNotInGap/O" );
        treeWriter.SetOutputBranchAddress("passed_topologicalScore", (void*)&passed_topologicalScore, "passed_topologicalScore/O" );
        treeWriter.SetOutputBranchAddress("passed_startNearVertex", (void*)&passed_startNearVertex, "passed_startNearVertex/O" );
        treeWriter.SetOutputBranchAddress("passed_openingAngle", (void*)&passed_openingAngle, "passed_openingAngle/O" );
        treeWriter.SetOutputBranchAddress("passed_likelyGoldenPion", (void*)&passed_likelyGoldenPion, "passed_likelyGoldenPion/O" );

        std::vector<int> generationCutValues;
        std::vector<float> trackScoreCutValues, vertexDistanceCutValues, trackLengthCutValues, protonChi2CutValues, muonChi2CutValues;
        // int track_num;
        // treeWriter.SetOutputBranchAddress("track_num", (void*)&track_num, "track_num/I" ); // Needed for root to know the size of the arrays
        // treeWriter.SetOutputBranchAddress("particle_cutValue_trackScore", (void*)&trackScoreCutValues, "particle_cutValue_trackScore[track_num]/F" );
        // treeWriter.SetOutputBranchAddress("particle_cutValue_vertexDistance", (void*)&vertexDistanceCutValues, "particle_cutValue_vertexDistance[track_num]/F" );
        // treeWriter.SetOutputBranchAddress("particle_cutValue_generation", (void*)&generationCutValues, "particle_cutValue_generation[track_num]/I" );
        // treeWriter.SetOutputBranchAddress("particle_cutValue_trackLength", (void*)&trackLengthCutValues, "particle_cutValue_trackLength[track_num]/F" );
        // treeWriter.SetOutputBranchAddress("particle_cutValue_protonChi2", (void*)&protonChi2CutValues, "particle_cutValue_protonChi2[track_num]/F" );
        // treeWriter.SetOutputBranchAddress("particle_cutValue_muonChi2", (void*)&muonChi2CutValues, "particle_cutValue_muonChi2[track_num]/F" );

        const std::vector<int>* pGenerationCutValues(&generationCutValues);
        const std::vector<float>* pTrackScoreCutValues(&trackScoreCutValues);
        const std::vector<float>* pVertexDistanceCutValues(&vertexDistanceCutValues);
        const std::vector<float>* pTrackLengthCutValues(&trackLengthCutValues);
        const std::vector<float>* pProtonChi2CutValues(&protonChi2CutValues);
        const std::vector<float>* pMuonChi2CutValues(&muonChi2CutValues);

        treeWriter.SetObjectOutputBranchAddress<std::vector<float>>("particle_cutValue_trackScore", pTrackScoreCutValues);
        treeWriter.SetObjectOutputBranchAddress<std::vector<float>>("particle_cutValue_vertexDistance", pVertexDistanceCutValues);
        treeWriter.SetObjectOutputBranchAddress<std::vector<int>>("particle_cutValue_generation", pGenerationCutValues);
        treeWriter.SetObjectOutputBranchAddress<std::vector<float>>("particle_cutValue_trackLength", pTrackLengthCutValues);
        treeWriter.SetObjectOutputBranchAddress<std::vector<float>>("particle_cutValue_protonChi2", pProtonChi2CutValues);
        treeWriter.SetObjectOutputBranchAddress<std::vector<float>>("particle_cutValue_muonChi2", pMuonChi2CutValues);

        int   nuPdgCode, category;
        float topologicalScore, flashChi2, piontruncatedMeandEdx, openingAngle, maxVertexDist, goldenPionBDTResponse;
        treeWriter.SetOutputBranchAddress("event_cutValue_nuPdgCode", (void*)&nuPdgCode, "event_cutValue_nuPdgCode/I" );
        treeWriter.SetOutputBranchAddress("event_cutValue_topologicalScore", (void*)&topologicalScore, "event_cutValue_topologicalScore/F" );
        treeWriter.SetOutputBranchAddress("event_cutValue_flashChi2", (void*)&flashChi2, "event_cutValue_flashChi2/F" );
        treeWriter.SetOutputBranchAddress("event_cutValue_pionTruncatedMeandEdx", (void*)&piontruncatedMeandEdx, "event_cutValue_pionTruncatedMeandEdx/F" );
        treeWriter.SetOutputBranchAddress("event_cutValue_maxVertexDist", (void*)&maxVertexDist, "event_cutValue_maxVertexDist/F" );
        treeWriter.SetOutputBranchAddress("event_cutValue_openingAngle", (void*)&openingAngle, "event_cutValue_openingAngle/F" );
        treeWriter.SetOutputBranchAddress("event_cutValue_goldenPionBDT", (void*)&goldenPionBDTResponse, "event_cutValue_goldenPionBDT/F" );
        treeWriter.SetOutputBranchAddress("category", (void*)&category, "category/I" );

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Truth variables
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bool isTrueCC0Pi, isTrueCC1Pi, isCC0PiSignal, isCC1PiSignal;
        treeWriter.SetOutputBranchAddress("true_cc0pi", (void*)&isTrueCC0Pi, "true_cc0pi/O" );
        treeWriter.SetOutputBranchAddress("true_cc1pi", (void*)&isTrueCC1Pi, "true_cc1pi/O" );
        treeWriter.SetOutputBranchAddress("cc0pi_signal", (void*)&isCC0PiSignal, "cc0pi_signal/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_signal", (void*)&isCC1PiSignal, "cc1pi_signal/O" );

        // Calculated with true muon and true pion
        float truthCC1PiMuonMomentum, truthCC1PiMuonCosTheta, truthCC1PiMuonPhi, truthCC1PiPionMomentum, truthCC1PiPionCosTheta, truthCC1PiPionPhi, truthCC1PiMuonPionAngle;
        int   truthCC1PiNProtons;
        treeWriter.SetOutputBranchAddress("cc1pi_truth_muonMomentum", (void*)&truthCC1PiMuonMomentum, "cc1pi_truth_muonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truth_muonCosTheta", (void*)&truthCC1PiMuonCosTheta, "cc1pi_truth_muonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truth_muonPhi", (void*)&truthCC1PiMuonPhi, "cc1pi_truth_muonPhi/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truth_pionMomentum", (void*)&truthCC1PiPionMomentum, "cc1pi_truth_pionMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truth_pionCosTheta", (void*)&truthCC1PiPionCosTheta, "cc1pi_truth_pionCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truth_pionPhi", (void*)&truthCC1PiPionPhi, "cc1pi_truth_pionPhi/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truth_muonPionAngle", (void*)&truthCC1PiMuonPionAngle, "cc1pi_truth_muonPionAngle/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truth_nProtons", (void*)&truthCC1PiNProtons, "cc1pi_truth_nProtons/I" );

        bool cc1pi_truthMuon_IsContained, cc1pi_truthPion_IsContained;
        float cc1pi_truthMuon_TrackLength, cc1pi_truthPion_TrackLength;
        treeWriter.SetOutputBranchAddress("cc1pi_truthMuon_IsContained", (void*)&cc1pi_truthMuon_IsContained, "cc1pi_truthMuon_IsContained/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_truthMuon_TrackLength", (void*)&cc1pi_truthMuon_TrackLength, "cc1pi_truthMuon_TrackLength/F" );
        treeWriter.SetOutputBranchAddress("cc1pi_truthPion_IsContained", (void*)&cc1pi_truthPion_IsContained, "cc1pi_truthPion_IsContained/O" );
        treeWriter.SetOutputBranchAddress("cc1pi_truthPion_TrackLength", (void*)&cc1pi_truthPion_TrackLength, "cc1pi_truthPion_TrackLength/F" );

        bool truecc1pi_recoMuon_IsContained;
        float truecc1pi_recoMuon_TrackLength, truecc1pi_recoMuon_protonBDTScore, truecc1pi_recoMuon_muonBDTScore;
        treeWriter.SetOutputBranchAddress("truecc1pi_recoMuon_IsContained", (void*)&truecc1pi_recoMuon_IsContained, "truecc1pi_recoMuon_IsContained/O" );
        treeWriter.SetOutputBranchAddress("truecc1pi_recoMuon_TrackLength", (void*)&truecc1pi_recoMuon_TrackLength, "truecc1pi_recoMuon_TrackLength/F" );
        treeWriter.SetOutputBranchAddress("truecc1pi_recoMuon_protonBDTScore", (void*)&truecc1pi_recoMuon_protonBDTScore, "truecc1pi_recoMuon_protonBDTScore/F" );
        treeWriter.SetOutputBranchAddress("truecc1pi_recoMuon_muonBDTScore", (void*)&truecc1pi_recoMuon_muonBDTScore, "truecc1pi_recoMuon_muonBDTScore/F" );

        bool truecc1pi_recoPion_IsContained;
        float truecc1pi_recoPion_TrackLength, truecc1pi_recoPion_protonBDTScore, truecc1pi_recoPion_muonBDTScore;
        treeWriter.SetOutputBranchAddress("truecc1pi_recoPion_IsContained", (void*)&truecc1pi_recoPion_IsContained, "truecc1pi_recoPion_IsContained/O" );
        treeWriter.SetOutputBranchAddress("truecc1pi_recoPion_TrackLength", (void*)&truecc1pi_recoPion_TrackLength, "truecc1pi_recoPion_TrackLength/F" );
        treeWriter.SetOutputBranchAddress("truecc1pi_recoPion_protonBDTScore", (void*)&truecc1pi_recoPion_protonBDTScore, "truecc1pi_recoPion_protonBDTScore/F" );
        treeWriter.SetOutputBranchAddress("truecc1pi_recoPion_muonBDTScore", (void*)&truecc1pi_recoPion_muonBDTScore, "truecc1pi_recoPion_muonBDTScore/F" );


        // Calculated with true muon and true leading proton
        float truthCC0PiMuonMomentum, truthCC0PiMuonCosTheta, truthCC0PiMuonPhi, truthCC0PiProtonMomentum, truthCC0PiProtonCosTheta, truthCC0PiProtonPhi, truthCC0PiMuonProtonAngle;
        int   truthCC0PiNProtons;
        treeWriter.SetOutputBranchAddress("cc0pi_truth_muonMomentum", (void*)&truthCC0PiMuonMomentum, "cc0pi_truth_muonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_truth_muonCosTheta", (void*)&truthCC0PiMuonCosTheta, "cc0pi_truth_muonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_truth_muonPhi", (void*)&truthCC0PiMuonPhi, "cc0pi_truth_muonPhi/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_truth_protonMomentum", (void*)&truthCC0PiProtonMomentum, "cc0pi_truth_protonMomentum/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_truth_protonCosTheta", (void*)&truthCC0PiProtonCosTheta, "cc0pi_truth_protonCosTheta/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_truth_protonPhi", (void*)&truthCC0PiProtonPhi, "cc0pi_truth_protonPhi/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_truth_muonProtonAngle", (void*)&truthCC0PiMuonProtonAngle, "cc0pi_truth_muonProtonAngle/F" );
        treeWriter.SetOutputBranchAddress("cc0pi_truth_nProtons", (void*)&truthCC0PiNProtons, "cc0pi_truth_nProtons/I" );

        std::vector<int> bestMatchedTruthPDGs;
        // treeWriter.SetOutputBranchAddress("track_bestMatchedTruthPDG", (void*)&bestMatchedTruthPDGs, "track_bestMatchedTruthPDG[track_num]/I" );
        const std::vector<int>* pBestMatchedTruthPDGs(&bestMatchedTruthPDGs);
        treeWriter.SetObjectOutputBranchAddress<std::vector<int>>("track_bestMatchedTruthPDG", pBestMatchedTruthPDGs);

        int mcNuPdg;
        treeWriter.SetOutputBranchAddress("mc_nu_pdg", (void*)&mcNuPdg, "mc_nu_pdg/I");


        std::vector<int> recoParticle_backtracked_absTruePDG;
        std::vector<bool> recoParticle_isContained, recoParticle_backtracked_isContained;
        std::vector<double> recoParticle_TrackLength, recoParticle_backtracked_TrackLength, recoParticle_protonBDTScore, recoParticle_muonBDTScore;

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Weights
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        float splineWeight;
        float tunedCVWeight;
        std::map<std::string, const std::vector<double>*> weightPtrMap;
        treeWriter.SetOutputBranchAddress("spline_weight", (void*)&splineWeight, "spline_weight/F");
        treeWriter.SetOutputBranchAddress("tuned_cv_weight", (void*)&tunedCVWeight, "tuned_cv_weight/F");

        // Open the input file for reading and enable the branches with systematic event weights (if required)
        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, isMC);
        if (isMC) readerPeLEE.EnableSystematicBranches(); // Todo: Is this correct/optimal?


        // Loop over the events in the file
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        // std::cout << "### Only processing 10\% of events (out of "<<nEvents<<") ###" << std::endl;
        for (unsigned int i = 0; i < nEvents; i++) // Todo: Remove!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            AnalysisHelper::PrintLoadingBar(i, nEvents);
            readerPeLEE.LoadEvent(i);

            const auto run = pEventPeLEE->metadata.run();
            const auto isGoodRun = (isDataBNB || isDataEXT) ? std::find(goodRuns.begin(), goodRuns.end(), run) != goodRuns.end() : true; // Apply good runs cuts to data
            if(!isGoodRun)
            {
                continue;
            }
            // else
            // {
            //     std::cout << "DEBUG - good run: " << run << std::endl;
            // }

            // std::cout << "DEBUG Point Y8" << std::endl;
            Event event(*pEventPeLEE, false);// true or false decides wether to cut generation!=2 particles
            const auto pEvent = std::make_shared<Event>(event);
            // std::cout << "DEBUG Point Y9" << std::endl;


            // ####################################################################
            // #################### Reco CC0pi ####################
            // ####################################################################
            const auto &[passedGoldenSelectionCC0Pi, cutsPassedCC0Pi, assignedPdgCodesCC0Pi] = selectionCC0Pi.Execute(pEvent);
            const auto passedGenericSelectionCC0Pi = SelectionHelper::IsCutPassed(cutsPassedCC0Pi, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoDataCC0Pi = (
                passedGenericSelectionCC0Pi
                    ? AnalysisHelper::GetRecoAnalysisDataCC0Pi(pEvent->reco, assignedPdgCodesCC0Pi, passedGoldenSelectionCC0Pi)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            cc0pi_reco_muonMomentum = recoDataCC0Pi.muonMomentum;
            cc0pi_reco_muonCosTheta = recoDataCC0Pi.muonCosTheta;
            cc0pi_reco_muonPhi = recoDataCC0Pi.muonPhi;
            cc0pi_reco_protonMomentum = recoDataCC0Pi.protonMomentum;
            cc0pi_reco_protonCosTheta = recoDataCC0Pi.protonCosTheta;
            cc0pi_reco_protonPhi = recoDataCC0Pi.protonPhi;
            cc0pi_reco_muonProtonAngle = recoDataCC0Pi.muonProtonAngle;
            cc0pi_reco_nProtons = recoDataCC0Pi.nProtons;


            cc0pi_leadingProton_protonBDTScore = -std::numeric_limits<float>::max();
            cc0pi_leadingProton_muonBDTScore = -std::numeric_limits<float>::max();
            unsigned int truthMatchedLeadingProtonIndex = std::numeric_limits<unsigned int>::max();
            if(passedGenericSelectionCC0Pi)
            {
                const auto leadingRecoProtonIndex = SelectionHelper::GetLeadingProtonCandidateIndex(pEvent->reco.particles, assignedPdgCodesCC0Pi);
                // std::cout<<"DEBUG Point E0.1"<<std::endl;
                const auto leadingRecoProton = pEvent->reco.particles.at(leadingRecoProtonIndex);
                // std::cout<<"DEBUG Point E1"<<std::endl;

                std::vector<float> protonFeatures;
                const auto hasProtonFeatures = BDTHelper::GetBDTFeatures(leadingRecoProton, BDTHelper::ProtonBDTFeatureNames, protonFeatures);
                std::vector<float> muonFeatures;
                const auto hasMuonFeatures = BDTHelper::GetBDTFeatures(leadingRecoProton, BDTHelper::MuonBDTFeatureNames, muonFeatures);
                if(!hasProtonFeatures || !hasMuonFeatures)
                    throw std::logic_error("Failed to get BDT features for proton or muon. This should not happen for the leading proton.");

                cc0pi_leadingProton_protonBDTScore = protonBDT.GetResponse(protonFeatures);
                cc0pi_leadingProton_muonBDTScore = muonBDT.GetResponse(muonFeatures);

                if(isMC)
                {
                    try{
                        truthMatchedLeadingProtonIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(leadingRecoProton,pEvent->truth.particles);
                    }
                    catch(const std::logic_error &){
                        std::cout<<"try failed for truthMatchedLeadingProtonIndex\n"<<std::endl;
                    }
                }
            }
            // std::cout<<"DEBUG Point E2"<<std::endl;
            if(truthMatchedLeadingProtonIndex!=std::numeric_limits<unsigned int>::max())
            {
                const auto proton = pEvent->truth.particles.at(truthMatchedLeadingProtonIndex);
                const auto dir = TVector3(proton.momentumX(), proton.momentumY(), proton.momentumZ()).Unit();
                const auto cosTheta = dir.Z();
                const auto phi = std::atan2(dir.Y(),dir.X());
                cc0pi_backtracked_protonPDG = proton.pdgCode();
                cc0pi_backtracked_protonMomentum = proton.momentum();
                cc0pi_backtracked_protonCosTheta = cosTheta;
                cc0pi_backtracked_protonPhi = phi;
                // cc0pi_backtracked_muonProtonAngle
            }
            else
            {
                cc0pi_backtracked_protonPDG = -std::numeric_limits<int>::max();
                cc0pi_backtracked_protonMomentum = -std::numeric_limits<float>::max();
                cc0pi_backtracked_protonCosTheta = -std::numeric_limits<float>::max();
                cc0pi_backtracked_protonPhi = -std::numeric_limits<float>::max();
            }
            // std::cout<<"DEBUG Point E3"<<std::endl;

            // std::cout << "DEBUG Point Y10" << std::endl;
            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            auto passesPhaseSpaceRecoCC0Pi = false;
            if (passedGenericSelectionCC0Pi)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceRecoCC0Pi = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValueCC0Pi.at(name)(recoDataCC0Pi);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceRecoCC0Pi = false;
                        break;
                    }
                }
            }
            passesGenericCC0Pi = passedGenericSelectionCC0Pi;
            passesGoldenCC0Pi = passedGoldenSelectionCC0Pi;
            isSelectedGenericCC0Pi = passedGenericSelectionCC0Pi && passesPhaseSpaceRecoCC0Pi;
            isSelectedGoldenCC0Pi = passedGoldenSelectionCC0Pi && passesPhaseSpaceRecoCC0Pi;

            // std::cout << "DEBUG Point Y11" << std::endl;
            // ################################################################
            // #################### Reco CC1pi ####################
            // ################################################################
            const auto &[passedGoldenSelection, cutsPassedCC1Pi, assignedPdgCodesCC1Pi] = selectionCC1Pi.Execute(pEvent);
            const auto passedGenericSelection = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, config.global.lastCutGeneric);

            // Get the reco analysis data (if available, otherwise set to dummy values)
            const auto recoDataCC1Pi = (
                passedGenericSelection
                    ? AnalysisHelper::GetRecoAnalysisData(pEvent->reco, assignedPdgCodesCC1Pi, passedGoldenSelection)
                    : AnalysisHelper::GetDummyAnalysisData()
            );

            cc1pi_reco_muonMomentum = recoDataCC1Pi.muonMomentum;
            cc1pi_reco_muonCosTheta = recoDataCC1Pi.muonCosTheta;
            cc1pi_reco_muonPhi = recoDataCC1Pi.muonPhi;
            cc1pi_reco_pionMomentum = recoDataCC1Pi.pionMomentum;
            cc1pi_reco_pionCosTheta = recoDataCC1Pi.pionCosTheta;
            cc1pi_reco_pionPhi = recoDataCC1Pi.pionPhi;
            cc1pi_reco_muonPionAngle = recoDataCC1Pi.muonPionAngle;
            cc1pi_reco_nProtons = recoDataCC1Pi.nProtons;

            unsigned int truthMatchedProtonIndex = std::numeric_limits<unsigned int>::max();
            cc1pi_recoMuon_IsContained = false;
            cc1pi_recoMuon_TrackLength = -std::numeric_limits<float>::max();
            cc1pi_recoMuon_protonBDTScore = -std::numeric_limits<float>::max();
            cc1pi_recoMuon_muonBDTScore = -std::numeric_limits<float>::max();

            unsigned int truthMatchedMuonIndex = std::numeric_limits<unsigned int>::max();
            cc1pi_recoPion_IsContained = false;
            cc1pi_recoPion_TrackLength = -std::numeric_limits<float>::max();
            cc1pi_recoPion_protonBDTScore = -std::numeric_limits<float>::max();
            cc1pi_recoPion_muonBDTScore = -std::numeric_limits<float>::max();
            if(passedGenericSelection)
            {
                const auto recoMuonCandidate = pEvent->reco.particles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1Pi, 13));
                cc1pi_recoMuon_IsContained = AnalysisHelper::IsContained(recoMuonCandidate);
                cc1pi_recoMuon_TrackLength = recoMuonCandidate.range();

                std::vector<float> muonProtonBDTFeatures;
                const auto muonHasProtonBDTFeatures = BDTHelper::GetBDTFeatures(recoMuonCandidate, BDTHelper::ProtonBDTFeatureNames, muonProtonBDTFeatures);
                std::vector<float> muonMuonBDTFeatures;
                const auto muonHasMuonBDTFeatures = BDTHelper::GetBDTFeatures(recoMuonCandidate, BDTHelper::MuonBDTFeatureNames, muonMuonBDTFeatures);
                // if(!muonHasProtonBDTFeatures)
                //     throw std::logic_error("Failed to get proton BDT features for muon.");  
                // if(!muonHasMuonBDTFeatures)
                //     throw std::logic_error("Failed to get muon BDT features for muon.");

                if(muonHasProtonBDTFeatures) cc1pi_recoMuon_protonBDTScore = protonBDT.GetResponse(muonProtonBDTFeatures);
                if(muonHasMuonBDTFeatures) cc1pi_recoMuon_muonBDTScore = muonBDT.GetResponse(muonMuonBDTFeatures);

                const auto recoPionCandidate = pEvent->reco.particles.at(AnalysisHelper::GetParticleIndexWithPdg(assignedPdgCodesCC1Pi, 211));
                cc1pi_recoPion_IsContained = AnalysisHelper::IsContained(recoPionCandidate);
                cc1pi_recoPion_TrackLength = recoPionCandidate.range();

                std::vector<float> pionProtonBDTFeatures;
                const auto pionHasProtonBDTFeatures = BDTHelper::GetBDTFeatures(recoPionCandidate, BDTHelper::ProtonBDTFeatureNames, pionProtonBDTFeatures);
                std::vector<float> pionMuonBDTFeatures;
                const auto pionHasMuonBDTFeatures = BDTHelper::GetBDTFeatures(recoPionCandidate, BDTHelper::MuonBDTFeatureNames, pionMuonBDTFeatures);
                // if(!pionHasProtonBDTFeatures)
                //     throw std::logic_error("Failed to get proton BDT features for pion.");  
                // if(!pionHasMuonBDTFeatures)
                //     throw std::logic_error("Failed to get muon BDT features for pion.");
                if(pionHasProtonBDTFeatures) cc1pi_recoPion_protonBDTScore = protonBDT.GetResponse(pionProtonBDTFeatures);
                if(pionHasMuonBDTFeatures) cc1pi_recoPion_muonBDTScore = muonBDT.GetResponse(pionMuonBDTFeatures);

                if(isMC)
                {
                    try{
                        truthMatchedProtonIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoPionCandidate,pEvent->truth.particles);
                    }
                    catch(const std::logic_error &){
                        std::cout<<"try failed for recoPionCandidate\n"<<std::endl;
                    }

                    try{
                        truthMatchedMuonIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoMuonCandidate,pEvent->truth.particles);
                        // std::cout<<"DEBUG try succeeded for recoMuonCandidate \\/\n"<<std::endl;
                    }
                    catch(const std::logic_error &){
                        std::cout<<"try failed for recoMuonCandidate\n"<<std::endl;
                    }
                }
            }

            if(truthMatchedProtonIndex!=std::numeric_limits<unsigned int>::max())
            {
                const auto matchedTrueProton = pEvent->truth.particles.at(truthMatchedProtonIndex);
                const auto dir = TVector3(matchedTrueProton.momentumX(), matchedTrueProton.momentumY(), matchedTrueProton.momentumZ()).Unit();
                const auto cosTheta = dir.Z();
                const auto phi = std::atan2(dir.Y(),dir.X());
                cc1pi_backtracked_protonPDG = matchedTrueProton.pdgCode();
                cc1pi_backtracked_protonMomentum = matchedTrueProton.momentum();
                cc1pi_backtracked_protonCosTheta = cosTheta;
                cc1pi_backtracked_protonPhi = phi;
                // cc0pi_backtracked_muonProtonAngle   
            }
            else
            {
                cc1pi_backtracked_protonPDG = -std::numeric_limits<int>::max();
                cc1pi_backtracked_protonMomentum = -std::numeric_limits<float>::max();
                cc1pi_backtracked_protonCosTheta = -std::numeric_limits<float>::max();
                cc1pi_backtracked_protonPhi = -std::numeric_limits<float>::max();
            }

            if(truthMatchedMuonIndex!=std::numeric_limits<unsigned int>::max())
            {
                const auto matchedTrueMuon = pEvent->truth.particles.at(truthMatchedMuonIndex);
                const auto dir = TVector3(matchedTrueMuon.momentumX(), matchedTrueMuon.momentumY(), matchedTrueMuon.momentumZ()).Unit();
                const auto cosTheta = dir.Z();
                const auto phi = std::atan2(dir.Y(),dir.X());
                cc1pi_backtracked_muonPDG = matchedTrueMuon.pdgCode();
                cc1pi_backtracked_muonMomentum = matchedTrueMuon.momentum();
                cc1pi_backtracked_muonCosTheta = cosTheta;
                cc1pi_backtracked_muonPhi = phi;
                // cc0pi_backtracked_muonMuonAngle   
            }
            else
            {
                cc1pi_backtracked_muonPDG = -std::numeric_limits<int>::max();
                cc1pi_backtracked_muonMomentum = -std::numeric_limits<float>::max();
                cc1pi_backtracked_muonCosTheta = -std::numeric_limits<float>::max();
                cc1pi_backtracked_muonPhi = -std::numeric_limits<float>::max();
            }

            // Here we apply reco-level phase-space restrictions
            // For any event that passes the generic selection, get the value of the kinematic quantity and check if it is outside of the
            // min/max values supplied in the binning. If so, then reject the event.
            auto passesPhaseSpaceRecoCC1Pi = false;
            if (passedGenericSelection)
            {
                // Start by assuming the event passes the phase-space cuts
                passesPhaseSpaceRecoCC1Pi = true;

                // Check the value of the kinematic quantities are within the phase-space limits
                for (const auto &[name, minMax] : phaseSpaceMap)
                {
                    const auto &[min, max] = minMax;
                    const auto value = getValueCC1Pi.at(name)(recoDataCC1Pi);

                    if (value < min || value > max)
                    {
                        passesPhaseSpaceRecoCC1Pi = false;
                        break;
                    }
                }
            }
            passesGenericCC1Pi = passedGenericSelection;
            passesGoldenCC1Pi = passedGoldenSelection;
            isSelectedGenericCC1Pi = passedGenericSelection && passesPhaseSpaceRecoCC1Pi;
            isSelectedGoldenCC1Pi = passedGoldenSelection && passesPhaseSpaceRecoCC1Pi;


            // std::cout << "DEBUG Point Y12" << std::endl;
            // if(isSelectedGenericCC0Pi || isSelectedGenericCC1Pi) {std::cout << "\nDEBUG isSelectedGenericCC0Pi: " << isSelectedGenericCC0Pi << " isSelectedGenericCC1Pi: " << isSelectedGenericCC1Pi << "\n" << std::endl;}

            passed_particleTrackScore = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "particleTrackScore");
            passed_particleVertexDistance = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "particleVertexDistance");
            passed_particleGeneration = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "particleGeneration");
            passed_particleTrackLength = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "particleTrackLength");
            passed_particleProtonChi2 = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "particleProtonChi2");
            passed_particleMuonChi2 = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "particleMuonChi2");
            passed_particleProtonChi2OverMuonChi2 = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "particleProtonChi2OverMuonChi2");
            passed_pandoraNuPDGIsNumu = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "pandoraNuPDGIsNumu");
            passed_daughterVerticesContained = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "daughterVerticesContained");
            passed_nuVertexFiducial = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "nuVertexFiducial");
            passed_topologicalOrFlashMatch = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "topologicalOrFlashMatch");
            passed_topologicalScoreCC = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "topologicalScoreCC");
            passed_min2Tracks = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "min2Tracks");
            passed_max1Uncontained = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "max1Uncontained");
            passed_2NonProtons = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "2NonProtons");
            passed_pionHasValiddEdx = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "pionHasValiddEdx");
            passed_pionNotInGap = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "pionNotInGap");
            passed_muonNotInGap = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "muonNotInGap");
            passed_topologicalScore = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "topologicalScore");
            passed_startNearVertex = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "startNearVertex");
            passed_openingAngle = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "openingAngle");
            passed_likelyGoldenPion = SelectionHelper::IsCutPassed(cutsPassedCC1Pi, "likelyGoldenPion");


            // ################### Cut variables for CC1pi ###################
            const auto &recoParticles = pEvent->reco.particles;
            trackScoreCutValues.clear();
            vertexDistanceCutValues.clear();
            generationCutValues.clear();
            trackLengthCutValues.clear();
            protonChi2CutValues.clear();
            muonChi2CutValues.clear();
            bestMatchedTruthPDGs.clear();
            for(unsigned int i = 0; i< recoParticles.size(); ++i)
            {
                const auto hasTrackFit = AnalysisHelper::HasTrackFit(recoParticles.at(i));
                if(!hasTrackFit) continue;
                bestMatchedTruthPDGs.push_back(-std::numeric_limits<int>::max());// Used later in MC block
                // particleTrackScore
                const auto trackScore = recoParticles.at(i).trackScore.IsSet() ? recoParticles.at(i).trackScore() : -std::numeric_limits<float>::max();
                trackScoreCutValues.push_back(trackScore);
                // particleVertexDistance
                const auto recoVertex = pEvent->reco.nuVertexNoSCC();
                const TVector3 start(recoParticles.at(i).vertexX(), recoParticles.at(i).vertexY(), recoParticles.at(i).vertexZ());
                vertexDistanceCutValues.push_back(sqrt((start-recoVertex).Mag2()));
                // particleGeneration
                generationCutValues.push_back(recoParticles.at(i).generation());
                // particleTrackLength
                trackLengthCutValues.push_back(recoParticles.at(i).range());
                // particleProtonChi2
                const auto chi2ForwardProtonW = recoParticles.at(i).chi2ForwardProtonW.IsSet() ? recoParticles.at(i).chi2ForwardProtonW() : -std::numeric_limits<float>::max();
                protonChi2CutValues.push_back(chi2ForwardProtonW);
                // particleMuonChi2
                const auto chi2ForwardMuonW = recoParticles.at(i).chi2ForwardMuonW.IsSet() ? recoParticles.at(i).chi2ForwardMuonW() : -std::numeric_limits<float>::max();
                muonChi2CutValues.push_back(chi2ForwardMuonW);
                // particleProtonChi2OverMuonChi2 - Use previous values
            }
            // track_num = trackScoreCutValues.size();

            // pandoraNuPDGIsNumu
            nuPdgCode = pEvent->reco.nuPdgCode.IsSet() ? pEvent->reco.nuPdgCode() : -std::numeric_limits<int>::max();
            // daughterVerticesContained - Use passes cut
            // nuVertexFiducial - Use passes cut
            // topologicalOrFlashMatch
            topologicalScore = pEvent->reco.selectedTopologicalScore.IsSet() ? pEvent->reco.selectedTopologicalScore() : -std::numeric_limits<float>::max();
            flashChi2 = pEvent->reco.flashChi2.IsSet() ? pEvent->reco.flashChi2() : -std::numeric_limits<float>::max();
            // min2Tracks - Use passes cut
            // max1Uncontained - Use passes cut
            // 2NonProtons - Use passes cut
            // pionHasValiddEdx
            auto muonIndex = -std::numeric_limits<int>::max();
            auto pionIndex = -std::numeric_limits<int>::max();
            bool maxOneMuonOnePion = true;
            for (unsigned int i = 0; i < assignedPdgCodesCC1Pi.size(); ++i)
            {
                if (assignedPdgCodesCC1Pi.at(i) == 13)
                {
                    if(muonIndex == -std::numeric_limits<int>::max()) muonIndex = i;
                    else maxOneMuonOnePion = false;
                } else if (assignedPdgCodesCC1Pi.at(i) == 211)
                {
                    if(pionIndex == -std::numeric_limits<int>::max()) pionIndex = i;
                    else maxOneMuonOnePion = false;
                }
            }

            piontruncatedMeandEdx = -std::numeric_limits<float>::max();
            openingAngle = -std::numeric_limits<float>::max();
            maxVertexDist = -std::numeric_limits<float>::max();
            goldenPionBDTResponse = -std::numeric_limits<float>::max();
            // Save some cut variables that are applied after muon and pion have been identified
            // Ensure that there is only one muon and one pion
            // All the precise selection logic is saved in the passed_* variables
            if(maxOneMuonOnePion && muonIndex!=-std::numeric_limits<int>::max() && pionIndex!=-std::numeric_limits<int>::max())
            {
                const auto muon = recoParticles.at(muonIndex);
                const auto pion = recoParticles.at(pionIndex);
                piontruncatedMeandEdx = pion.truncatedMeandEdx.IsSet() ? pion.truncatedMeandEdx() : std::numeric_limits<float>::max();

                if(muon.directionX.IsSet() && muon.directionY.IsSet() && muon.directionZ.IsSet() && pion.directionX.IsSet() && pion.directionY.IsSet() && pion.directionZ.IsSet())
                {
                    const auto muonDir = TVector3(muon.directionX(), muon.directionY(), muon.directionZ()).Unit();
                    const auto pionDir = TVector3(pion.directionX(), pion.directionY(), pion.directionZ()).Unit();
                    openingAngle = muonDir.Angle(pionDir);
                }

                const auto recoVertex = pEvent->reco.nuVertex();
                for (const auto &particle : recoParticles)
                {
                    // Skip particles without a track fit
                    if (!AnalysisHelper::HasTrackFit(particle)) continue;

                    // Get the distance between the particle's start position and the vertex
                    const TVector3 start(particle.startX(), particle.startY(), particle.startZ());
                    const auto vertexDist = (start - recoVertex).Mag2();

                    // Insist that this isn't too large
                    if (vertexDist > maxVertexDist) maxVertexDist = vertexDist;
                }
                // const auto pMuonBDT = std::make_shared<BDTHelper::BDT>("muon", BDTHelper::MuonBDTFeatureNames);
                // const auto pProtonBDT = std::make_shared<BDTHelper::BDT>("proton", BDTHelper::ProtonBDTFeatureNames);
                // const auto pGoldenPionBDT = std::make_shared<BDTHelper::BDT>("goldenPion", BDTHelper::GoldenPionBDTFeatureNames);
                // std::cout<<"DEBUG P Point 3.5"<<std::endl;
                // std::vector<float> features;
                // std::cout<<"DEBUG P Point 3.6"<<std::endl;
                // const auto hasFeatures = BDTHelper::GetBDTFeatures(pion, BDTHelper::GoldenPionBDTFeatureNames, features);
                // std::cout<<"DEBUG P Point 3.7"<<std::endl;
                // goldenPionBDTResponse = hasFeatures ? pGoldenPionBDT->GetResponse(features) : -std::numeric_limits<float>::max();
                goldenPionBDTResponse = 0;
            }

            // pionNotInGap - Use passes cut
            // muonNotInGap - Use passes cut
            // openingAngle
            // topologicalScore - Reuse previous topological score
            // startNearVertex
            // likelyGoldenPion

            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // Event variable
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            // treeWriter.SetOutputBranchAddress("muon_candidate_idx", &pEventPeLEE->muon_candidate_idx_, "muon_candidate_idx/I" );
            // treeWriter.SetOutputBranchAddress("pion_candidate_idx", &pEventPeLEE->pion_candidate_idx_, "pion_candidate_idx/I" );
            // treeWriter.SetOutputBranchAddress("lead_p_candidate_idx", &pEventPeLEE->lead_p_candidate_idx_, "lead_p_candidate_idx/I" );

            isTrueCC1Pi = false;
            isTrueCC0Pi = false;
            isCC1PiSignal = false;
            isCC0PiSignal = false;

            truthCC1PiMuonMomentum = -std::numeric_limits<float>::max();
            truthCC1PiMuonCosTheta = -std::numeric_limits<float>::max();
            truthCC1PiMuonPhi = -std::numeric_limits<float>::max();
            truthCC1PiPionMomentum = -std::numeric_limits<float>::max();
            truthCC1PiPionCosTheta = -std::numeric_limits<float>::max();
            truthCC1PiPionPhi = -std::numeric_limits<float>::max();
            truthCC1PiMuonPionAngle = -std::numeric_limits<float>::max();
            truthCC1PiNProtons = -std::numeric_limits<int>::max();

            truthCC0PiMuonMomentum = -std::numeric_limits<float>::max();
            truthCC0PiMuonCosTheta = -std::numeric_limits<float>::max();
            truthCC0PiMuonPhi = -std::numeric_limits<float>::max();
            truthCC0PiProtonMomentum = -std::numeric_limits<float>::max();
            truthCC0PiProtonCosTheta = -std::numeric_limits<float>::max();
            truthCC0PiProtonPhi = -std::numeric_limits<float>::max();
            truthCC0PiMuonProtonAngle = -std::numeric_limits<float>::max();
            truthCC0PiNProtons = -std::numeric_limits<int>::max();

            mcNuPdg = -std::numeric_limits<int>::max();

            cc1pi_truthMuon_IsContained = false;
            cc1pi_truthMuon_TrackLength = -std::numeric_limits<float>::max();
            cc1pi_truthPion_IsContained = false;
            cc1pi_truthPion_TrackLength = -std::numeric_limits<float>::max();

            truecc1pi_recoMuon_IsContained = false;
            truecc1pi_recoMuon_TrackLength = -std::numeric_limits<float>::max();
            truecc1pi_recoMuon_protonBDTScore = -std::numeric_limits<float>::max();
            truecc1pi_recoMuon_muonBDTScore = -std::numeric_limits<float>::max();

            truecc1pi_recoPion_IsContained = false;
            truecc1pi_recoPion_TrackLength = -std::numeric_limits<float>::max();
            truecc1pi_recoPion_protonBDTScore = -std::numeric_limits<float>::max();
            truecc1pi_recoPion_muonBDTScore = -std::numeric_limits<float>::max();

            

            if (isMC) // For BNB data that's all we need to do!
            {
                // #####################################################
                // #################### Truth CC0pi ####################
                // #####################################################
                // Determine if this is truly a CC0Pi event
                isTrueCC0Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC0Pi(pEvent, config.global.useAbsPdg, config.global.protonMomentumThreshold);

                // Get the truth analysis data (if available, otherwise set to dummy values)
                const auto truthDataCC0Pi = (
                    (isTrueCC0Pi)
                        ? AnalysisHelper::GetTruthAnalysisDataCC0Pi(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                        : AnalysisHelper::GetDummyAnalysisData()
                );

                // Here we apply truth-level phase-space restrictions
                // For all true CC0Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
                // event is not classed as "signal"
                bool passesPhaseSpaceTruthCC0Pi = false;
                if (isTrueCC0Pi)
                {
                    // Start by assuming the event passes the phase-space cuts
                    passesPhaseSpaceTruthCC0Pi = true;

                    // Check the value of the kinematic quantities are within the phase-space limits
                    for (const auto &[name, minMax] : phaseSpaceMap)
                    {
                        // if(name == "pionMomentum") continue; // Not really compatible with the pion momentum in the CC1pi selection // todo check this
                        const auto &[min, max] = minMax;
                        const auto value = getValueCC0Pi.at(name)(truthDataCC0Pi);

                        if (value < min || value > max)
                        {
                            passesPhaseSpaceTruthCC0Pi = false;
                            break;
                        }
                    }
                }
                isCC0PiSignal = isTrueCC0Pi && passesPhaseSpaceTruthCC0Pi;

                // #####################################################
                // #################### Truth CC1pi ####################
                // #####################################################
                isTrueCC1Pi = (isOverlay || isDetVar || isNuWro) && AnalysisHelper::IsTrueCC1Pi(pEvent, config.global.useAbsPdg);

                // Get the truth analysis data (if available, otherwise set to dummy values)
                const auto truthDataCC1Pi = (
                    isTrueCC1Pi
                        ? AnalysisHelper::GetTruthAnalysisData(pEvent->truth, config.global.useAbsPdg, config.global.protonMomentumThreshold)
                        : AnalysisHelper::GetDummyAnalysisData()
                );

                // Here we apply truth-level phase-space restrictions
                // For all true CC1Pi events, we check if the values of each kinematic variable are within the supplied limits. If not then the
                // event is not classed as "signal"
                bool passesPhaseSpaceTruthCC1Pi = false;
                if (isTrueCC1Pi)
                {
                    // Start by assuming the event passes the phase-space cuts
                    passesPhaseSpaceTruthCC1Pi = true;

                    // Check the value of the kinematic quantities are within the phase-space limits
                    for (const auto &[name, minMax] : phaseSpaceMap)
                    {
                        const auto &[min, max] = minMax;
                        const auto value = getValueCC1Pi.at(name)(truthDataCC1Pi);

                        if (value < min || value > max)
                        {
                            passesPhaseSpaceTruthCC1Pi = false;
                            break;
                        }
                    }
                }
                isCC1PiSignal = isTrueCC1Pi && passesPhaseSpaceTruthCC1Pi;

                truthCC1PiMuonMomentum = truthDataCC1Pi.muonMomentum;
                truthCC1PiMuonCosTheta = truthDataCC1Pi.muonCosTheta;
                truthCC1PiMuonPhi = truthDataCC1Pi.muonPhi;
                truthCC1PiPionMomentum = truthDataCC1Pi.pionMomentum;
                truthCC1PiPionCosTheta = truthDataCC1Pi.pionCosTheta;
                truthCC1PiPionPhi = truthDataCC1Pi.pionPhi;
                truthCC1PiMuonPionAngle = truthDataCC1Pi.muonPionAngle;
                truthCC1PiNProtons = truthDataCC1Pi.nProtons;

                truthCC0PiMuonMomentum = truthDataCC0Pi.muonMomentum;
                truthCC0PiMuonCosTheta = truthDataCC0Pi.muonCosTheta;
                truthCC0PiMuonPhi = truthDataCC0Pi.muonPhi;
                truthCC0PiProtonMomentum = truthDataCC0Pi.protonMomentum;
                truthCC0PiProtonCosTheta = truthDataCC0Pi.protonCosTheta;
                truthCC0PiProtonPhi = truthDataCC0Pi.protonPhi;
                truthCC0PiMuonProtonAngle = truthDataCC0Pi.muonProtonAngle;
                truthCC0PiNProtons = truthDataCC0Pi.nProtons;

                mcNuPdg = pEvent->truth.nuPdgCode();

                if(isTrueCC1Pi)
                {
                    const auto trueMuonIndex = AnalysisHelper::GetTrueMuonIndex(pEvent->truth, config.global.useAbsPdg);
                    const auto muon = pEvent->truth.particles.at(trueMuonIndex);
                    cc1pi_truthMuon_IsContained = AnalysisHelper::IsContained(muon);
                    cc1pi_truthMuon_TrackLength = std::sqrt(std::pow(muon.endX()-muon.startX(),2)+std::pow(muon.endY()-muon.startY(),2)+std::pow(muon.endZ()-muon.startZ(),2));

                    const auto truePionIndex = AnalysisHelper::GetTruePionIndex(pEvent->truth, config.global.useAbsPdg);
                    const auto pion = pEvent->truth.particles.at(truePionIndex);
                    cc1pi_truthPion_IsContained = AnalysisHelper::IsContained(pion);
                    cc1pi_truthPion_TrackLength = std::sqrt(std::pow(pion.endX()-pion.startX(),2)+std::pow(pion.endY()-pion.startY(),2)+std::pow(pion.endZ()-pion.startZ(),2));

                    // Find the reco particle that best matche the true muon and pion 
                    auto recoMuonIndex = -std::numeric_limits<int>::max();
                    auto recoPionIndex = -std::numeric_limits<int>::max();
                    for(unsigned int i=0; i<pEvent->reco.particles.size(); i++)
                    {
                        const auto recoParticle = pEvent->reco.particles.at(i);
                        auto matchedTrueParticleIndex = std::numeric_limits<unsigned int>::max();
                        try{
                            matchedTrueParticleIndex = AnalysisHelper::GetBestMatchedTruthParticleIndex(recoParticle, pEvent->truth.particles);
                        }
                        catch(const std::logic_error &){
                            std::cout<<"try failed for recoParticle\n"<<std::endl;
                        }

                        if(matchedTrueParticleIndex!=std::numeric_limits<unsigned int>::max())
                        {
                            const auto truthParticle = pEvent->truth.particles.at(matchedTrueParticleIndex);
                            recoParticle_backtracked_absTruePDG.push_back(truthParticle.pdgCode());
                            recoParticle_backtracked_isContained.push_back(AnalysisHelper::IsContained(truthParticle));
                            const auto length = std::sqrt(std::pow(truthParticle.endX()-truthParticle.startX(),2)+std::pow(truthParticle.endY()-truthParticle.startY(),2)+std::pow(truthParticle.endZ()-truthParticle.startZ(),2));
                            recoParticle_backtracked_TrackLength.push_back(length);
                        }
                        else
                        {
                            recoParticle_backtracked_absTruePDG.push_back(-std::numeric_limits<int>::max());
                            recoParticle_backtracked_isContained.push_back(false);
                            recoParticle_backtracked_TrackLength.push_back(-std::numeric_limits<float>::max());
                        }
                        
                        if(AnalysisHelper::HasTrackFit(recoParticle))
                        {
                            recoParticle_isContained.push_back(AnalysisHelper::IsContained(recoParticle));
                            recoParticle_TrackLength.push_back(recoParticle.range());
                        }
                        else
                        {
                            recoParticle_isContained.push_back(false);
                            recoParticle_TrackLength.push_back(-std::numeric_limits<float>::max());
                        }
                        std::vector<float> protonBDTFeatures;
                        const auto hasProtonBDTFeatures = BDTHelper::GetBDTFeatures(recoParticle, BDTHelper::ProtonBDTFeatureNames, protonBDTFeatures);
                        std::vector<float> muonBDTFeatures;
                        const auto hasMuonBDTFeatures = BDTHelper::GetBDTFeatures(recoParticle, BDTHelper::MuonBDTFeatureNames, muonBDTFeatures);                        
                        if(hasProtonBDTFeatures) recoParticle_protonBDTScore.push_back(protonBDT.GetResponse(protonBDTFeatures));
                        else recoParticle_protonBDTScore.push_back(-std::numeric_limits<float>::max());
                        
                        if(hasMuonBDTFeatures) recoParticle_muonBDTScore.push_back(muonBDT.GetResponse(muonBDTFeatures));
                        else recoParticle_muonBDTScore.push_back(-std::numeric_limits<float>::max());

                        if(hasMuonBDTFeatures && !hasProtonBDTFeatures) throw std::logic_error("Muon BDT features found but not proton BDT features");
                    }
                }


                unsigned int trackCounter = 0; // Count non skipped particles
                for(unsigned int i = 0; i< recoParticles.size(); ++i)
                {
                    const auto hasTrackFit = AnalysisHelper::HasTrackFit(recoParticles.at(i));
                    if(!hasTrackFit) continue; // Needed for consistency with cuts
                    // const auto matchedTruthParticle = AnalysisHelper::GetBestMatchedTruthParticle(recoParticles.at(i), pEvent->truth.particles, true);
                    // const auto mcpdg = matchedTruthParticle.pdgCode();
                    if(pEvent->reco.particles.at(i).pdgBacktracked.IsSet())
                    {
                        const auto matchedPDG = pEvent->reco.particles.at(i).pdgBacktracked(); // todo verifiy and use translation layer
                        bestMatchedTruthPDGs.at(trackCounter) = matchedPDG;
                    }
                    trackCounter++;
                }
            }

            splineWeight = isMC ? pEventPeLEE->truth.weightSpline() : 1.f;
            tunedCVWeight = isMC ? pEventPeLEE->truth.weightTune() : 1.f;

            weightPtrMap.clear();
            if(isMC)
            {
                // General systematic weights
                for (auto &weights : *pEventPeLEE->truth.weights.GetAddress())
                {
                    weightPtrMap[weights.first] = &weights.second;
                    treeWriter.SetObjectOutputBranchAddress<std::vector<double>>("weight_" + weights.first, weightPtrMap.at(weights.first)); //todo find a way to also move this up
                }
            }

            const std::vector<int>* pRecoParticle_backtracked_absTruePDG = &recoParticle_backtracked_absTruePDG;
            treeWriter.SetObjectOutputBranchAddress<std::vector<int>>("trueCC1pi_recoParticle_backtracked_absTruePDG_vector", pRecoParticle_backtracked_absTruePDG);
            const std::vector<bool>* pRecoParticle_isContained = &recoParticle_isContained;
            treeWriter.SetObjectOutputBranchAddress<std::vector<bool>>("trueCC1pi_recoParticle_isContained_vector", pRecoParticle_isContained);
            const std::vector<bool>* pRecoParticle_backtracked_isContained = &recoParticle_backtracked_isContained;
            treeWriter.SetObjectOutputBranchAddress<std::vector<bool>>("trueCC1pi_recoParticle_backtracked_isContained_vector", pRecoParticle_backtracked_isContained);
            const std::vector<double>* pRecoParticle_TrackLength = &recoParticle_TrackLength;
            treeWriter.SetObjectOutputBranchAddress<std::vector<double>>("trueCC1pi_recoParticle_TrackLength_vector", pRecoParticle_TrackLength);
            const std::vector<double>* pRecoParticle_backtracked_TrackLength = &recoParticle_backtracked_TrackLength;
            treeWriter.SetObjectOutputBranchAddress<std::vector<double>>("trueCC1pi_recoParticle_backtracked_TrackLength_vector", pRecoParticle_backtracked_TrackLength);
            const std::vector<double>* pRecoParticle_protonBDTScore = &recoParticle_protonBDTScore;
            treeWriter.SetObjectOutputBranchAddress<std::vector<double>>("trueCC1pi_recoParticle_protonBDTScore_vector", pRecoParticle_protonBDTScore);
            const std::vector<double>* pRecoParticle_muonBDTScore = &recoParticle_muonBDTScore;
            treeWriter.SetObjectOutputBranchAddress<std::vector<double>>("trueCC1pi_recoParticle_muonBDTScore_vector", pRecoParticle_muonBDTScore);

            const auto type = PlottingHelper::GetPlotStyle(sampleType, pEvent, true);
            switch (type) {
                case PlottingHelper::PlotStyle::BNBData:
                    category = 0;
                    break;
                case PlottingHelper::PlotStyle::NumuCC1PiChargedNonGolden:
                    category = isCC1PiSignal ? 1 : 11;
                    break;
                case PlottingHelper::PlotStyle::NumuCC1PiChargedGolden:
                    category = isCC1PiSignal ? 2 : 11;
                    break;
                case PlottingHelper::PlotStyle::External:
                    category = 3;
                    break;
                case PlottingHelper::PlotStyle::Dirt:
                    category = 4;
                    break;
                case PlottingHelper::PlotStyle::NonFiducial:
                    category = 5;
                    break;
                case PlottingHelper::PlotStyle::Nue:
                    category = 6;
                    break;
                case PlottingHelper::PlotStyle::NC:
                    category = 7;
                    break;
                case PlottingHelper::PlotStyle::NumuCC0Pi:
                    category = isCC0PiSignal ? 12 : 8; // Here 8 is CC0pi events not selected by CC0pi selection
                    break;
                case PlottingHelper::PlotStyle::NumuCC1PiZero:
                    category = 9;
                    break;
                case PlottingHelper::PlotStyle::NumuCCOther:
                    category = 10;
                    break;
                // category 11: CC1pi events outside phase space   
                // category 12: Selected CC0pi events   
                default:
                    throw std::invalid_argument("Invalid type encountered");
            }

            treeWriter.CreateNoNewBranches(); // Only creates branches for the first run through the loop // todo get rid of this method
            treeWriter.Fill();

            recoParticle_backtracked_absTruePDG.clear();
            recoParticle_isContained.clear();
            recoParticle_backtracked_isContained.clear();
            recoParticle_TrackLength.clear();
            recoParticle_backtracked_TrackLength.clear();
            recoParticle_protonBDTScore.clear();
            recoParticle_muonBDTScore.clear();
        } // End of event-level iteration

    } // End of file-level iterration
    std::cout << "-----------------Done-----------------" << std::endl;


}

} // namespace ubcc1pi_macros
