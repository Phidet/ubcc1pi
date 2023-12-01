/**
 *  @file  ubcc1pi_standalone/Objects/Config.h
 *
 *  @brief The header file for the config structure
 */

#ifndef UBCC1PI_STANDALONE_OBJECTS_CONFIG
#define UBCC1PI_STANDALONE_OBJECTS_CONFIG

#include <string>
#include <unordered_map>

#include "ubcc1pi_standalone/Helpers/PlottingHelper.h"
#include "ubcc1pi_standalone/Helpers/CrossSectionHelper.h"
#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

namespace ubcc1pi
{

/**
 *  @brief  The config struture
 */
struct Config
{

    struct Inputs {
        const std::vector<float> overlayPOT {1.28278e+21, 1.01592e+21, 1.31355e+21, 2.48772e+20, 7.64399e+20, 4.64842e+20, 8.66958e+20, 9.70634e+20}; // run 1, 2, 3, 4a, 4b, 4c, 4d, 5
        const std::vector<float> nuWroPOT {3.08672e+20, 3.14546e+20, 3.154e+20, -1.f, -1.f, -1.f, -1.f, -1.f};
        const std::vector<float> dirtPOT {3.1606e+20, 9.50551e+20, 3.26609e+20, 1.10496e+20, 3.02251e+20, 1.77721e+20, 3.52244e+20, 3.57315e+20}; // Done
        const std::vector<float> dataEXTTriggers {65744587.0, 153928909.0, 207380078.0, 27596585.0, 86095801.0, 50162518.0, 74593111.0, 107494728.0};
        const std::vector<float> dataBNBTor860WCut {1.669e+20, 2.616e+20, 2.562e+20, 3.607e+19, 1.39e+20, 8.586e+19, -1.f, 1.23e+20};
        const std::vector<float> dataEXTE1DCNTWCut {36906271.0, 62557664.0, 61744410.0, 7822167.0, 31158202.0, 19052861.0, 10960320.0, 28136255.0};
        // todo remove dedx ????
        const std::vector<float> detVarPOT {1.3422e+21, 6.28536e+20, 9.95366e+20, 1.31618e+21, 1.3532e+21, 9.177e+20, 6.344e+20, 6.23808e+20, 1.26005e+21, 1.34671e+21, 1.11179e+21, 1.34253e+21}; // CV, CVextra, WireModdEdx, LYAtteNuation, LYDown, LYRayleigh, Recomb2, SCE, WireModThetaXZ, WireModThetaYZ, WireModX, WireModYZ

        /**
        *  @brief  Input files
        *
        * Specify the type and location of PeLEE ntuples
        * The first entry is the run number, the second is the sample type, the third is whether to use the file, the fourth is the file path
        *
        */
        std::vector< std::tuple<unsigned int, float, AnalysisHelper::SampleType, bool, std::string> > files = {
 
            // Beam off
            { /* 1  */ 1, dataEXTE1DCNTWCut.at(0)/dataEXTTriggers.at(0), AnalysisHelper::DataEXT, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run1.root"},
            { /* 2  */ 2, dataEXTE1DCNTWCut.at(1)/dataEXTTriggers.at(1), AnalysisHelper::DataEXT, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run2.root"},
            { /* 3  */ 3, dataEXTE1DCNTWCut.at(2)/dataEXTTriggers.at(2), AnalysisHelper::DataEXT, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/bnb_beam_off_peleeTuple_uboone_v08_00_00_70_run3.root"},
            { /* 4a */ 4, dataEXTE1DCNTWCut.at(3)/dataEXTTriggers.at(3), AnalysisHelper::DataEXT, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_off_peleeTuple_uboone_v08_00_00_73_run4a.root"},
            { /* 4b */ 5, dataEXTE1DCNTWCut.at(4)/dataEXTTriggers.at(4), AnalysisHelper::DataEXT, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/Run4b_bnb_beamOff_PeLee_ntuples_run4b_ana.root"},
            { /* 4c */ 6, dataEXTE1DCNTWCut.at(5)/dataEXTTriggers.at(5), AnalysisHelper::DataEXT, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_off_peleeTuple_uboone_v08_00_00_73_run4c.root"},
            { /* 4d */ 7, dataEXTE1DCNTWCut.at(6)/dataEXTTriggers.at(6), AnalysisHelper::DataEXT, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_off_peleeTuple_uboone_v08_00_00_73_run4d.root"},
            { /* 5  */ 8, dataEXTE1DCNTWCut.at(7)/dataEXTTriggers.at(7), AnalysisHelper::DataEXT, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_off_peleeTuple_uboone_v08_00_00_72_run5.root"},
 
            // Beam on
            { /* 1  */ 1, 1.f, AnalysisHelper::DataBNB, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run1_C1.root"},
            { /* 2  */ 2, 1.f, AnalysisHelper::DataBNB, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run2.root"},
            { /* 3  */ 3, 1.f, AnalysisHelper::DataBNB, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/bnb_beam_on_peleeTuple_uboone_v08_00_00_70_run3.root"},
            { /* 4a */ 4, 1.f, AnalysisHelper::DataBNB, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/Run4a_bnb_beamOn_PeLee_ntuples_run4a_ana.root"},
            { /* 4b */ 5, 1.f, AnalysisHelper::DataBNB, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_on_peleeTuple_uboone_v08_00_00_74_run4b.root"},
            { /* 4c */ 6, 1.f, AnalysisHelper::DataBNB, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_on_peleeTuple_uboone_v08_00_00_73_run4c.root"},
            { /* 4d */ 7, 1.f, AnalysisHelper::DataBNB, false, ""},
            { /* 5  */ 8, 1.f, AnalysisHelper::DataBNB, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/bnb_beam_on_peleeTuple_uboone_v08_00_00_75_run5.root"},

            // Overlay dirt
            { /* 1  */ 1, dataBNBTor860WCut.at(0)/dirtPOT.at(0), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run1_dirt.root"},
            { /* 2  */ 2, dataBNBTor860WCut.at(1)/dirtPOT.at(1), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/overlay_peleeTuple_uboone_v08_00_00_70_run2_dirt.root"},
            { /* 3  */ 3, dataBNBTor860WCut.at(2)/dirtPOT.at(2), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run3_dirt.root"},
            { /* 4a */ 4, dataBNBTor860WCut.at(3)/dirtPOT.at(3), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/Run4a_BNB_dirt_overlay_Pandora_Unified_Reco2_fixed_run4a_ana.root"},
            { /* 4b */ 5, dataBNBTor860WCut.at(4)/dirtPOT.at(4), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4b_dirt.root"},
            { /* 4c */ 6, dataBNBTor860WCut.at(5)/dirtPOT.at(5), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4c_dirt.root"},
            { /* 4d */ 7, dataBNBTor860WCut.at(6)/dirtPOT.at(6), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4d_dirt.root"},
            { /* 5  */ 8, dataBNBTor860WCut.at(7)/dirtPOT.at(7), AnalysisHelper::Dirt, true, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run5_dirt.root"},

            // Overlay MC
            { /* 1  */ 1, dataBNBTor860WCut.at(0)/overlayPOT.at(0), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu.root"},
            { /* 2  */ 2, dataBNBTor860WCut.at(1)/overlayPOT.at(1), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu.root"},
            { /* 3  */ 3, dataBNBTor860WCut.at(2)/overlayPOT.at(2), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu.root"},
            { /* 4a */ 4, dataBNBTor860WCut.at(3)/overlayPOT.at(3), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_73_run4a_nu.root"},
            { /* 4b */ 5, dataBNBTor860WCut.at(4)/overlayPOT.at(4), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_73_run4b_nu.root"},
            { /* 4c */ 6, dataBNBTor860WCut.at(5)/overlayPOT.at(5), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4c_nu.root"},
            { /* 4d */ 7, dataBNBTor860WCut.at(6)/overlayPOT.at(6), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run4d_nu.root"},
            { /* 5  */ 8, dataBNBTor860WCut.at(7)/overlayPOT.at(7), AnalysisHelper::Overlay, false, "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/uboone/persistent/users/jdetje/pelee_v08_00_00_70/overlay_peleeTuple_uboone_v08_00_00_70_run5_nu.root"},

            // DetVar
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(0), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_CV_virtualParentsFix.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(1), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_CVextra.root"},
            // {3, dataBNBTor860WCut.at(2)/detVarPOT.at(2), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_WireModdEdx_virtualParentsFix.root.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(3), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_LYAtteNuation.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(4), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_LYDown.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(5), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_LYRayleigh_virtualParentsFix.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(6), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_Recomb2.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(7), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_SCE.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(8), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_WireModThetaXZ.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(9), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_WireModThetaYZ.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(10), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_WireModX_virtualParentsFix.root"},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(11), AnalysisHelper::DetectorVariation, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/detvar_run3_peleeTuple_uboone_v08_00_00_73_WireModYZ.root"},

            // NuWro MC
            {1, dataBNBTor860WCut.at(0)/nuWroPOT.at(0), AnalysisHelper::NuWro, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/nuwro_fakedata_peleeTuple_uboone_run1.root"},
            {2, dataBNBTor860WCut.at(1)/nuWroPOT.at(1), AnalysisHelper::NuWro, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/nuwro_fakedata_peleeTuple_uboone_run2.root"},
            {3, dataBNBTor860WCut.at(2)/nuWroPOT.at(2), AnalysisHelper::NuWro, false, "/pnfs/uboone/persistent/users/jdetje/pelee_v08_00_00_70/combined/nuwro_fakedata_peleeTuple_uboone_run3.root"}
        };

    };
    Inputs input;

struct Outputs {
        const std::vector<float> overlayPOT {1.28278e+21, 1.01592e+21, 1.31355e+21, 2.48772e+20, 7.64399e+20, 4.64842e+20, 8.66958e+20, 9.70634e+20}; // run 1, 2, 3, 4a, 4b, 4c, 4d, 5
        const std::vector<float> nuWroPOT {3.08672e+20, 3.14546e+20, 3.154e+20, -1.f}; // run 1, 2, 3, 4c
        const std::vector<float> dirtPOT {3.1606e+20, 9.50551e+20, 3.26609e+20, 1.77721e+20, 3.52244e+20, 3.57315e+20}; // run 1, 2, 3, 4c, 4d, 5
        const std::vector<float> dataEXTTriggers {65744587.f, 153928909.f, 207380078.f, 11212608.f}; // run 1, 2, 3, 4c
        const std::vector<float> dataBNBTor860WCut {1.669e+20, 2.616e+20, 2.562e+20, 0, 0, 8.586e+19, 0, 0}; // run 1, 2, 3, 4a, 4b, 4c, 4d, 5
        const std::vector<float> dataEXTE1DCNTWCut {36906271.f, 62557664.f, 61744410.f, 19476010.f}; // run 1, 2, 3, 4c // todo run 4c very low. Is this right?
        // todo remove dedx ????
        const std::vector<float> detVarPOT {1.3422e+21, 6.28536e+20, 9.95366e+20, 1.31618e+21, 1.3532e+21, 9.177e+20, 6.344e+20, 6.23808e+20, 1.26005e+21, 1.34671e+21, 1.11179e+21, 1.34253e+21}; // CV, CVextra, WireModdEdx, LYAtteNuation, LYDown, LYRayleigh, Recomb2, SCE, WireModThetaXZ, WireModThetaYZ, WireModX, WireModYZ

        /**
        *  @brief  Input files
        *
        * Specify the type and location of PeLEE ntuples
        * The first entry is the run number, the second is the sample type, the third is whether to use the file, the fourth is the file path
        *
        */
        std::vector< std::tuple<unsigned int, float, AnalysisHelper::SampleType, bool, std::string> > files = {
            // Beam off
            {1, dataEXTE1DCNTWCut.at(0)/dataEXTTriggers.at(0), AnalysisHelper::DataEXT, false, ""},
            {2, dataEXTE1DCNTWCut.at(1)/dataEXTTriggers.at(1), AnalysisHelper::DataEXT, false, ""},
            {3, dataEXTE1DCNTWCut.at(2)/dataEXTTriggers.at(2), AnalysisHelper::DataEXT, false, ""},
            {4, dataEXTE1DCNTWCut.at(3)/dataEXTTriggers.at(3), AnalysisHelper::DataEXT, false, ""},
            // {6, dataEXTE1DCNTWCut.at(3)/dataEXTTriggers.at(3), AnalysisHelper::DataEXT, false, ""},
            // Beam on
            {1, 1.f, AnalysisHelper::DataBNB, false, ""},
            {2, 1.f, AnalysisHelper::DataBNB, false, ""},
            {3, 1.f, AnalysisHelper::DataBNB, false, ""},
            {4, 1.f, AnalysisHelper::DataBNB, false, ""},
            // Overlay dirt
            {1, dataBNBTor860WCut.at(0)/dirtPOT.at(0), AnalysisHelper::Dirt, false, ""},
            {2, dataBNBTor860WCut.at(1)/dirtPOT.at(1), AnalysisHelper::Dirt, false, ""},
            {3, dataBNBTor860WCut.at(2)/dirtPOT.at(2), AnalysisHelper::Dirt, false, ""},
            {4, dataBNBTor860WCut.at(3)/dirtPOT.at(3), AnalysisHelper::Dirt, false, ""},
            {5, dataBNBTor860WCut.at(3)/dirtPOT.at(4), AnalysisHelper::Dirt, false, ""},
            {6, dataBNBTor860WCut.at(3)/dirtPOT.at(5), AnalysisHelper::Dirt, false, ""},

            // Overlay MC
            { /* 1  */ 1, dataBNBTor860WCut.at(0)/overlayPOT.at(0), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_70_run1_nu_ubcc1pi.root"},
            { /* 2  */ 2, dataBNBTor860WCut.at(1)/overlayPOT.at(1), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_70_run2_nu_ubcc1pi.root"},
            { /* 3  */ 3, dataBNBTor860WCut.at(2)/overlayPOT.at(2), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_70_run3_nu_ubcc1pi.root"},
            { /* 4a */ 4, dataBNBTor860WCut.at(3)/overlayPOT.at(3), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_73_run4a_nu_ubcc1pi.root"},
            { /* 4b */ 5, dataBNBTor860WCut.at(4)/overlayPOT.at(4), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_73_run4b_nu_ubcc1pi.root"},
            { /* 4c */ 6, dataBNBTor860WCut.at(5)/overlayPOT.at(5), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_70_run4c_nu_ubcc1pi.root"},
            { /* 4d */ 7, dataBNBTor860WCut.at(6)/overlayPOT.at(6), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_70_run4d_nu_ubcc1pi.root"},
            { /* 5  */ 8, dataBNBTor860WCut.at(7)/overlayPOT.at(7), AnalysisHelper::Overlay, true, "/uboone/data/users/jdetje/ubcc1piOutput/100Percent_5/overlay_peleeTuple_uboone_v08_00_00_70_run5_nu_ubcc1pi.root"},

            // DetVar
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(0), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(1), AnalysisHelper::DetectorVariation, false, ""},
            // {3, dataBNBTor860WCut.at(2)/detVarPOT.at(2), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(3), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(4), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(5), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(6), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(7), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(8), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(9), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(10), AnalysisHelper::DetectorVariation, false, ""},
            {3, dataBNBTor860WCut.at(2)/detVarPOT.at(11), AnalysisHelper::DetectorVariation, false, ""},
            // NuWro MC
            {1, dataBNBTor860WCut.at(0)/nuWroPOT.at(0), AnalysisHelper::NuWro, false, ""},
            {2, dataBNBTor860WCut.at(1)/nuWroPOT.at(1), AnalysisHelper::NuWro, false, ""},
            {3, dataBNBTor860WCut.at(2)/nuWroPOT.at(2), AnalysisHelper::NuWro, false, ""}
        };

    };
    Outputs output;


    /**
     *  @brief  Input files structure
     *
     * Specify the location of special ubcc1pi analysis root files that the standalone analysis code uses as input
     * You can produce them from art root files by running the ubcc1pi::AnalysisFileWriter module
     *
     */
    // struct FileList
    // {
    //     // September overlays files contain systematic weights
    //     std::string overlaysFileName, dirtFileName, dataEXTFileName, dataBNBFileName, nuWroFileName;
    //     std::vector< std::pair<std::string, std::string> > detVarFiles;
    // };

    // FileList filesRun1;
    // filesRun1.overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_08Oct21.root"; ///< Overlays file name input
    // // std::string overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_08Oct21_test_head20.root";
    // filesRun1.dirtFileName     = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_dirt_08Oct21.root";     ///< Dirt file name input
    // // std::string dirtFileName     = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/dirt_final/49396507_11/ubcc1piAnalysis.root";
    // filesRun1.dataEXTFileName  = "/pnfs/uboone/persistent/users/kduffy/ubcc1pi/ubcc1pi_extbnb_run1_combined_5Jan2021.root";  ///< EXT data file name input
    // filesRun1.dataBNBFileName  = "/pnfs/uboone/persistent/users/kduffy/ubcc1pi/ubcc1pi_bnb_run1-C1_5Jan2021.root";  ///< BNB data file name input
    // // std::string dataBNBFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_08Oct21_test_tail20.root";
    // // std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run1_28Feb22.root";  ///< NuWro file name input
    // // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisFiltered.root";
    // // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisUnfiltered.root";
    // // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisUnfilteredDetVar.root";
    // // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisFilteredDetVar.root";
    // // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisFilteredDetVar_larger.root";
    // // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisUnfilteredDetVar_larger.root";
    // // std::string nuWroFileName = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run1_corrected_backup/ubcc1piAnalysisCombined.root";
    // filesRun1.nuWroFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run1_28Apr22_noCCProduction.root";
    // /**
    //  *  @brief  The detector variation files
    //  *
    //  *          The events in each detector variation sample (in the same run) are identical apart from the detector parameter that's
    //  *          been changed - this way we can limit the effects of statistical variations. As a result we also need to have a
    //  *          central-value (CV) sample in which nothing is changed that we can compare to. Here, only some of the variations are
    //  *          available for run 1. Instead we use run3b, for which the varaitions are available. In every case the variation is
    //  *          considered wrt to the relevant "CV" file, and the fractional difference is quantity we care about.
    //  */
    // filesRun1.detVarFiles = {
    //     // Run-1 files
    //     {"CV_SCE",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run1_19Oct.root"},
    //     {"CV",             "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run1_19Oct.root"}, // Same CV for all run 1 variations
    //     {"LYDown",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYDown_run1_19Oct.root"},
    //     {"LYRayleigh",     "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run1_19Oct.root"},
    //     {"SCE",            "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_SCE_run1_19Oct.root"},
    //     {"Recomb2",        "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_Recomb2_run1_19Oct.root"},
    //     {"WireModX",       "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModX_run1_19Oct.root"},
    //     {"WireModYZ",      "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run1_19Oct.root"},
    //     {"WireModThetaXZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run1_19Oct.root"},
    //     {"WireModThetaYZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run1_19Oct.root"}
    // };



    // FileList filesRun2;
    // // September overlays files contain systematic weights
    // filesRun2.overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_run2_07Nov21.root"; ///< Overlays file name inputs
    // filesRun2.dirtFileName     = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_dirt_run2_07Nov21.root";     ///< Dirt file name input
    // filesRun2.dataEXTFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_extbnb_run2_17Jan22.root";  ///< EXT data file name input
    // filesRun2.dataBNBFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_bnb_run2_08Nov21.root";  ///< BNB data file name input

    // // std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run2_13Dec21.root";
    // // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run2a/nuwro_run2a_test.root";
    // // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run2b/nuwro_run2b_test.root";
    // filesRun2.nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run2ab_06May22_noCCProduction.root";

    // std::vector< std::pair<std::string, std::string> > detVarFiles = {
    //     // Run-3 files
    //     {"CV_SCE",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run3ab_09Mar22.root"},
    //     {"CV",             "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run3ab_high_09Mar22.root"},
    //     {"LYDown",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYDown_run3ab_high_09Mar22.root"},
    //     {"LYRayleigh",     "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run3ab_high_09Mar22.root"},
    //     {"SCE",            "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_SCE_run3ab_09Mar22.root"},
    //     {"Recomb2",        "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_Recomb2_run3ab_09Mar22.root"},
    //     {"WireModX",       "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModX_run3ab_high_09Mar22.root"},
    //     {"WireModYZ",      "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run3ab_high_09Mar22.root"},
    //     {"WireModThetaXZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run3ab_high_09Mar22.root"},
    //     {"WireModThetaYZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run3ab_high_09Mar22.root"}
    // };

    // filesRun2.detVarFiles = detVarFiles;

    // FileList filesRun3;
    // // September overlays files contain systematic weights
    // filesRun3.overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_run3_07Nov21.root"; ///< Overlays file name inputs
    // filesRun3.dirtFileName     = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_dirt_run3_07Nov21.root";     ///< Dirt file name input
    // filesRun3.dataEXTFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_extbnb_run3_17Jan22.root";  ///< EXT data file name input
    // filesRun3.dataBNBFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_bnb_run3_17Jan22.root";  ///< BNB data file name input

    // // std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run3_14Dec21.root";
    // // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run3a/nuwro_run3a_test.root";
    // // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run3b/nuwro_run3b_test.root";
    // filesRun3.nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run3ab_06May22_noCCProduction.root";
    // filesRun3.detVarFiles = detVarFiles;










    struct FilesRun1
    {
        // September overlays files contain systematic weights
        std::string overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_08Oct21.root"; ///< Overlays file name input
        // std::string overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_08Oct21_test_head20.root";

        std::string dirtFileName     = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_dirt_08Oct21.root";     ///< Dirt file name input
        // std::string dirtFileName     = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/dirt_final/49396507_11/ubcc1piAnalysis.root";

        std::string dataEXTFileName  = "/pnfs/uboone/persistent/users/kduffy/ubcc1pi/ubcc1pi_extbnb_run1_combined_5Jan2021.root";  ///< EXT data file name input

        std::string dataBNBFileName  = "/pnfs/uboone/persistent/users/kduffy/ubcc1pi/ubcc1pi_bnb_run1-C1_5Jan2021.root";  ///< BNB data file name input
        // std::string dataBNBFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_08Oct21_test_tail20.root";

        // std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run1_28Feb22.root";  ///< NuWro file name input
        // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisFiltered.root";
        // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisUnfiltered.root";
        // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisUnfilteredDetVar.root";
        // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisFilteredDetVar.root";
        // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisFilteredDetVar_larger.root";
        // std::string nuWroFileName  = "/uboone/app/users/jdetje/cc1pi_handover_57_2/test_ccinc_21Apr/ubcc1piAnalysisUnfilteredDetVar_larger.root";
        // std::string nuWroFileName = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run1_corrected_backup/ubcc1piAnalysisCombined.root";
        std::string nuWroFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run1_28Apr22_noCCProduction.root";

        /**
         *  @brief  The detector variation files
         *
         *          The events in each detector variation sample (in the same run) are identical apart from the detector parameter that's
         *          been changed - this way we can limit the effects of statistical variations. As a result we also need to have a
         *          central-value (CV) sample in which nothing is changed that we can compare to. Here, only some of the variations are
         *          available for run 1. Instead we use run3b, for which the varaitions are available. In every case the variation is
         *          considered wrt to the relevant "CV" file, and the fractional difference is quantity we care about.
         */
        std::vector< std::pair<std::string, std::string> > detVarFiles = {
            // Run-1 files
            {"CV_SCE",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run1_19Oct.root"},
            {"CV",             "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run1_19Oct.root"}, // Same CV for all run 1 variations
            {"LYDown",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYDown_run1_19Oct.root"},
            {"LYRayleigh",     "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run1_19Oct.root"},
            {"SCE",            "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_SCE_run1_19Oct.root"},
            {"Recomb2",        "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_Recomb2_run1_19Oct.root"},
            {"WireModX",       "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModX_run1_19Oct.root"},
            {"WireModYZ",      "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run1_19Oct.root"},
            {"WireModThetaXZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run1_19Oct.root"},
            {"WireModThetaYZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run1_19Oct.root"}
        };

        // std::vector< std::pair<std::string, std::string> > detVarFiles = {
        //     // Run-1 files
        //     {"CVRun1",         "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_CV_run1.root"},
        //     {"LYDown",         "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_LYDown_run1.root"},
        //     {"LYRayleigh",     "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run1.root"},

        //     // Run-3b files
        //     {"CVRun3b",        "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_CV_run3b.root"},
        //     {"SCE",            "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_SCE_run3b.root"},
        //     {"Recomb2",        "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_Recomb2_run3b.root"},
        //     // {"WireModX",       "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModX_run3b.root"},
        //     // {"WireModYZ",      "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run3b.root"},
        //     // {"WireModThetaXZ", "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run3b.root"},
        //     // {"WireModThetaYZ", "/uboone/data/users/asmith/ubcc1pi/samples/oct2020/samples/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run3b.root"}
        // };
    };
    FilesRun1 filesRun1; ///< The input files

    struct FilesRun2
    {
        // September overlays files contain systematic weights
        std::string overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_run2_07Nov21.root"; ///< Overlays file name inputs
        std::string dirtFileName     = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_dirt_run2_07Nov21.root";     ///< Dirt file name input
        std::string dataEXTFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_extbnb_run2_17Jan22.root";  ///< EXT data file name input
        std::string dataBNBFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_bnb_run2_08Nov21.root";  ///< BNB data file name input

        // std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run2_13Dec21.root";
        // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run2a/nuwro_run2a_test.root";
        // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run2b/nuwro_run2b_test.root";
        std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run2ab_06May22_noCCProduction.root";

        /**
         *  @brief  The detector variation files
         *
         *          The events in each detector variation sample (in the same run) are identical apart from the detector parameter that's
         *          been changed - this way we can limit the effects of statistical variations. As a result we also need to have a
         *          central-value (CV) sample in which nothing is changed that we can compare to. Here, only some of the variations are
         *          available for run 1. Instead we use run3b, for which the varaitions are available. In every case the variation is
         *          considered wrt to the relevant "CV" file, and the fractional difference is quantity we care about.
         */
        std::vector< std::pair<std::string, std::string> > detVarFiles = {
            // Run-3 files
            {"CV_SCE",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run3ab_09Mar22.root"},
            {"CV",             "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run3ab_high_09Mar22.root"},
            {"LYDown",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYDown_run3ab_high_09Mar22.root"},
            {"LYRayleigh",     "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run3ab_high_09Mar22.root"},
            {"SCE",            "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_SCE_run3ab_09Mar22.root"},
            {"Recomb2",        "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_Recomb2_run3ab_09Mar22.root"},
            {"WireModX",       "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModX_run3ab_high_09Mar22.root"},
            {"WireModYZ",      "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run3ab_high_09Mar22.root"},
            {"WireModThetaXZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run3ab_high_09Mar22.root"},
            {"WireModThetaYZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run3ab_high_09Mar22.root"}
        };

    };
    FilesRun2 filesRun2; ///< The input files

    struct FilesRun3
    {
        // September overlays files contain systematic weights
        std::string overlaysFileName = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_run3_07Nov21.root"; ///< Overlays file name inputs
        std::string dirtFileName     = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_dirt_run3_07Nov21.root";     ///< Dirt file name input
        std::string dataEXTFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_extbnb_run3_17Jan22.root";  ///< EXT data file name input
        std::string dataBNBFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_bnb_run3_17Jan22.root";  ///< BNB data file name input

        // std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run3_14Dec21.root";
        // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run3a/nuwro_run3a_test.root";
        // std::string nuWroFileName  = "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/nuwro_run3b/nuwro_run3b_test.root";
        std::string nuWroFileName  = "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_nuwro_run3ab_06May22_noCCProduction.root";

        /**
         *  @brief  The detector variation files
         *
         *          The events in each detector variation sample (in the same run) are identical apart from the detector parameter that's
         *          been changed - this way we can limit the effects of statistical variations. As a result we also need to have a
         *          central-value (CV) sample in which nothing is changed that we can compare to. Here, only some of the variations are
         *          available for run 1. Instead we use run3b, for which the varaitions are available. In every case the variation is
         *          considered wrt to the relevant "CV" file, and the fractional difference is quantity we care about.
         */

        std::vector< std::pair<std::string, std::string> > detVarFiles = {
            // Run-3 files
            {"CV_SCE",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run3ab_09Mar22.root"},
            {"CV",             "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run3ab_high_09Mar22.root"},
            {"LYDown",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYDown_run3ab_high_09Mar22.root"},
            {"LYRayleigh",     "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run3ab_high_09Mar22.root"},
            {"SCE",            "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_SCE_run3ab_09Mar22.root"},
            {"Recomb2",        "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_Recomb2_run3ab_09Mar22.root"},
            {"WireModX",       "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModX_run3ab_high_09Mar22.root"},
            {"WireModYZ",      "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run3ab_high_09Mar22.root"},
            {"WireModThetaXZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run3ab_high_09Mar22.root"},
            {"WireModThetaYZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run3ab_high_09Mar22.root"}
        };

        // std::vector< std::pair<std::string, std::string> > detVarFiles = {
        //     // Run-3 files
        //     {"CV_SCE",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_SCE_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"CV",             "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_CV_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"LYDown",         "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYDown_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"LYRayleigh",     "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_LYRayleigh_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"SCE",            "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_SCE_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"Recomb2",        "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_Recomb2_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"WireModX",       "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModX_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"WireModYZ",      "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModYZ_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"WireModThetaXZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaXZ_run3a_11Aug22_corrected5.root"}, // <----changed
        //     {"WireModThetaYZ", "/uboone/data/users/jdetje/ubcc1pi/sep2020/ubcc1piAnalysis_overlays_DetVar_WireModThetaYZ_run3a_11Aug22_corrected5.root"} // <----changed
        // };

        // std::vector< std::pair<std::string, std::string> > detVarFiles = {
        //     {"CV_run3a", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_CV_run3a_low_corrected4/36473026_1/ubcc1piAnalysis.root"},
        //     {"CV_run3b", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_CV_run3b_high_corrected4/36473043_2/ubcc1piAnalysis.root"},
        //     {"LYDown_run3a", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_LYDown_run3a_low_corrected4/36473029_1/ubcc1piAnalysis.root"},
        //     {"LYDown_run3b", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_LYDown_run3b_high_corrected4/36473051_1/ubcc1piAnalysis.root"},
        //     {"LYRayleigh_run3a", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_LYRayleigh_run3a_low_corrected4/36473039_1/ubcc1piAnalysis.root"},
        //     {"LYRayleigh_run3b", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_LYRayleigh_run3b_high_corrected4/36473054_1/ubcc1piAnalysis.root"},
        //     {"WireModThetaXZ_run3a", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModThetaXZ_run3a_low_corrected4/36473045_1/ubcc1piAnalysis.root"},
        //     {"WireModThetaXZ_run3b", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModThetaXZ_run3b_high_corrected4/36473067_1/ubcc1piAnalysis.root"},
        //     {"WireModThetaYZ_run3a", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModThetaYZ_run3a_low_corrected4/36473048_1/ubcc1piAnalysis.root"},
        //     {"WireModThetaYZ_run3b", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModThetaYZ_run3b_high_corrected4/36473064_1/ubcc1piAnalysis.root"},
        //     {"WireModX_run3a", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModX_run3a_low_corrected4/36473033_1/ubcc1piAnalysis.root"},
        //     {"WireModX_run3b", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModX_run3b_high_corrected4/36473057_1/ubcc1piAnalysis.root"},
        //     {"WireModYZ_run3a", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModYZ_run3a_low_corrected4/36473036_1/ubcc1piAnalysis.root"},
        //     {"WireModYZ_run3b", "/pnfs/uboone/scratch/users/jdetje/ubcc1pi/overlays_DetVar_WireModYZ_run3b_high_corrected4/36473061_1/ubcc1piAnalysis.root"}
        // };

    };
    FilesRun3 filesRun3; ///< The input files


    // -------------------------------------------------------------------------------------------------------------------------------------

    struct NormList
    {
        float overlaysPOT, dirtPOT, dataEXTTriggers, dataBNBTor875WCut, dataBNBE1DCNTWCut, nuWroPOT;
        std::unordered_map<std::string, float> detVarPOTs;
    };

    /**
     *  @brief  The sample normalisations structure
     *
     * These values can be found by:
     *      - Running the ubcc1pi::CountPOT macro (for overlays, dirt and detector variation files)
     *      - Running the ubcc1pi::GetRunSubrunList macro (for data), and then running Zarko's tool:
     *          /uboone/app/users/zarko/getDataInfo.py -v2 --run-subrun-list runSubrunList.txt
     */
    struct NormsRun1 : NormList
    {
        float  overlaysPOT        = 1.195e+21;// TODO: REMOVE!! 6.42546e+19 POT for overlays_08Oct21_test_head20.root //1.195e+21;      ///< The total POT for the overlays MC
        // Andy S files: 1.22447e+21, Kirsty D files (Jan 21):  1.18578e+21, Philip D files (Oct 08) 1.195e+21
        float  dirtPOT            = 2.63523e+20;    ///< The total POT for the dirt MC
        // Andy S files: 2.85049e+20, Kirsty D files (Jan 21): 2.91414e+20, Philip D files (Oct 08) 2.63523e+20
        float  dataEXTTriggers    = 64322029.0;     ///< The EXT triggers for the EXT data
        // Andy S files: 62540367.0, Kirsty D files (Jan 21): 64322029.0
        float  dataBNBTor875WCut  = 1.532e+20;      ///< The POT measured by the 875m toroid (with quality cuts)
        // Andy S files: 1.455e+20, Kirsty D files (Jan 21): 1.532e+20
        float  dataBNBE1DCNTWCut  = 34076199.0;     ///< The BNB spills sent by the accelerator division (with quality cuts)
        // Andy S files: 32339256.0, Kirsty D files (Jan 21): 34076199.0

        float nuWroPOT            = 3.08669e+20;//3.0867e+20;      ///< The total POT for the nuWro MC
        // float  nuWroTor875WCut  = 1.627e+19;      ///< The POT measured by the 875m toroid (with quality cuts)
        // float  nuWroE1DCNTWCut  = 3578340.0;     ///< The BNB spills sent by the accelerator division (with quality cuts)

        /**
         *  @brief  The detector variation POTs
         */
        std::unordered_map<std::string, float> detVarPOTs = {
            {"CV",             6.13708e+20},    // Andy/Kirst version: 1.14339e+20
            {"CV_SCE",         6.13708e+20},
            {"LYDown",         6.14119e+20},    // Andy/Kirst version: 1.05031e+20
            {"LYRayleigh",     6.21255e+20},    // Andy/Kirst version: 1.06661e+20
            {"SCE",            6.11946e+20},    // Andy/Kirst version: 1.02517e+20
            {"Recomb2",        6.18821e+20},    // Andy/Kirst version: 1.00832e+20
            {"WireModX",       6.16338e+20},    // Andy/Kirst version: 1.09739e+20
            {"WireModYZ",      6.1371e+20},     // Andy/Kirst version: 1.10877e+20
            {"WireModThetaXZ", 6.18155e+20},    // Andy/Kirst version: 1.12906e+20
            {"WireModThetaYZ", 6.15852e+20}     // Andy/Kirst version: 1.09244e+20
        };

        /**
         *  @brief  The detector variation POTs
         */
        // std::unordered_map<std::string, float> detVarPOTs = {
        //     {"CVRun1",         1.14339e+20},
        //     {"LYDown",         1.05031e+20},
        //     {"LYRayleigh",     1.06661e+20},
        //     {"CVRun3b",        9.82298e+19},
        //     {"SCE",            1.02517e+20},
        //     {"Recomb2",        1.00832e+20},
        //     {"WireModX",       1.09739e+20},
        //     {"WireModYZ",      1.10877e+20},
        //     {"WireModThetaXZ", 1.12906e+20},
        //     {"WireModThetaYZ", 1.09244e+20}
        // };
    };
    NormsRun1 normsRun1; ///< The sample normalisations

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  The sample normalisations structure
     *
     * These values can be found by:
     *      - Running the ubcc1pi::CountPOT macro (for overlays, dirt and detector variation files)
     *      - Running the ubcc1pi::GetRunSubrunList macro (for data), and then running Zarko's tool:
     *          /uboone/app/users/zarko/getDataInfo.py -v2 --run-subrun-list runSubrunList.txt
     */
    struct NormsRun2 : NormList
    {
        float  overlaysPOT        = 4.38501e+20;
        float  dirtPOT            = 4.40468e+20;
        float  dataEXTTriggers    = 87761971.0;
        float  dataBNBTor875WCut  = 1.94e+20;
        float  dataBNBE1DCNTWCut  = 42667752.0;

        // float  nuWroTor875WCut  = 2.278e+19;      ///< The POT measured by the 875m toroid (with quality cuts)
        // float  nuWroE1DCNTWCut  = 5397028.0;     ///< The BNB spills sent by the accelerator division (with quality cuts)
        float nuWroPOT            = 3.14564e+20;      ///< The total POT for the nuWro MC

        /**
         *  @brief  The detector variation POTs
         */
        std::unordered_map<std::string, float> detVarPOTs = {
            {"CV",             1.35602e+21},//4.85185e+20},
            {"CV_SCE",         5.46746e+20},//4.85185e+20},
            {"LYDown",         1.36375e+21},//4.72026e+20},
            {"LYRayleigh",     1.35301e+21},//5.11749e+20},
            {"SCE",            6.7823e+20},//6.17584e+20},
            {"Recomb2",        6.96454e+20},//6.34396e+20},
            {"WireModX",       1.33697e+21},//6.265e+20},
            {"WireModYZ",      1.35337e+21},//6.24617e+20},
            {"WireModThetaXZ", 1.27112e+21},//6.16842e+20},
            {"WireModThetaYZ", 1.35387e+21}//6.16682e+20}
        };
    };
    NormsRun2 normsRun2; ///< The sample normalisations

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  The sample normalisations structure
     *
     * These values can be found by:
     *      - Running the ubcc1pi::CountPOT macro (for overlays, dirt and detector variation files)
     *      - Running the ubcc1pi::GetRunSubrunList macro (for data), and then running Zarko's tool:
     *          /uboone/app/users/zarko/getDataInfo.py -v2 --run-subrun-list runSubrunList.txt
     */
    struct NormsRun3 : NormList
    {
        float  overlaysPOT        = 8.92215e+20;
        float  dirtPOT            = 1.61137e+20;
        float  dataEXTTriggers    = 126185602.0;//130439550.0;
        float  dataBNBTor875WCut  = 2.051e+20;//1.996e+20;
        float  dataBNBE1DCNTWCut  = 49186781.0;//47886106.0;

        // float  nuWroTor875WCut  = 2.195e+19;      ///< The POT measured by the 875m toroid (with quality cuts)
        // float  nuWroE1DCNTWCut  = 5215049.0;     ///< The BNB spills sent by the accelerator division (with quality cuts)
        float nuWroPOT            = 3.15446e+20;      ///< The total POT for the nuWro MC
        /**
         *  @brief  The detector variation POTs
         */
        std::unordered_map<std::string, float> detVarPOTs = {
            {"CV",             1.35602e+21},//4.85185e+20},
            {"CV_SCE",         5.46746e+20},//4.85185e+20},
            {"LYDown",         1.36375e+21},//4.72026e+20},
            {"LYRayleigh",     1.35301e+21},//5.11749e+20},
            {"SCE",            6.7823e+20},//6.17584e+20},
            {"Recomb2",        6.96454e+20},//6.34396e+20},
            {"WireModX",       1.33697e+21},//6.265e+20},
            {"WireModYZ",      1.35337e+21},//6.24617e+20},
            {"WireModThetaXZ", 1.27112e+21},//6.16842e+20},
            {"WireModThetaYZ", 1.35387e+21}//6.16682e+20}
        };
    };
    NormsRun3 normsRun3; ///< The sample normalisations

    const std::vector<NormList> inputNormalisations = {normsRun1, normsRun2, normsRun3};

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  The flux structure
     */
    struct FluxRun
    {
        std::string  fileName             = "/uboone/data/users/asmith/ubcc1pi/samples/flux/MCC9_FluxHist_volTPCActive.root";  ///< The path to the file containing the flux in each systematic universe
        float        pot                  = 4997 * 5e8;                                                                        ///< The total number of protons-on-target simulated in the input flux file

        /**
        *  @brief  A mapping from neutrino PDG codes, to the names used in the flux file
        */
        std::map<int, std::string> nuPdgToHistName = {
            {+12, "nue"},
            {-12, "nuebar"},
            {+14, "numu"},
            {-14, "numubar"},
        };

        std::vector<int> nuPdgsSignal = {-14, +14}; ///< The neutrino PDG codes for the fluxes to use in the cross-section calculation

        std::string  nomHistPattern       = "hENEUTRINO_cv";                         ///< The pattern for the nominal flux historam names (NEUTRINO is replaced by one of the names in nuPdgToHistMap)
        std::string  variationDirPattern  = "NEUTRINO_ms_PARAMNAME";                 ///< The pattern for the directory corresponding to each systematic paramter (NEUTRINO and PARAMNAME are replaced)
        std::string  variationHistPattern = "hENEUTRINO_PARAMNAME_ms_UNIVERSEINDEX"; ///< The pattern for the flux histogram corresponding to each universe in a given directory (NEUTRINO, PARAMNAME and UNIVERSEINDEX are replace)
    };
    FluxRun fluxRun1; ///< The flux
    FluxRun flux; ///< The flux

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Global configuration options structure
     */
    struct Global
    {
        bool        useAbsPdg                   = true;              ///< If we should use absolute PDG codes (this makes pi+ == pi- in the signal definition)
        bool        countProtonsInclusively     = true;              ///< If we should count protons inclusively (as Xp), or exclusively as (0p, 1p, 2p, ...)
        std::string lastCutGeneric              = "startNearVertex"; ///< The last cut of the generic selection (remaining cuts are part of the golden selection)
        float       protonMomentumThreshold     = 0.3f;              ///< The minimum proton momentum to be counted [GeV]
        float       targetDensity               = 8.44191f;          ///< The number of target nuclei per unit volume - units e23 / cm^3
        std::string selection                   = "Default";         ///< Which selection to use (can be "CCInclusive","Default", or "CC0pi")
        bool        axisTitles                  = true;              ///< If we want to draw axis lables and titles on the plots (if false, they are not drawn so you can add your own later)
        bool        scaleByBinWidth             = true;
        bool        useCC0piConstraint          = true;              ///< If we should use the CC0pi selection to constrain CC1pi cross-section
        bool        useBNBAsData                = true;             ///< If we should run with real data
        bool        useNuWroAsData              = false;              ///< If we should run with NuWro as data
        bool        useGenieAsData              = false;             ///< If we should run with Genie as data
        bool        useDetVar                   = true;              ///< If we should run with detector variations
        bool        fitInSystematicUniverses    = true;              ///< If we should fit not only in nominal but also in systematic universes
        bool        useEfficiencyCorrection     = false;             ///< If we should use the efficiency-effect-free smearing matrix to use when creating the plots
        std::vector<unsigned int> runs          = {1,2,3};           ///< The runs to use in the analysis
        std::string outputPath                  = "/uboone/data/users/jdetje/ubcc1piOutput/"; ///< The location to save processed PeLEE files
        std::string goodRunListFile             = "/uboone/data/users/jdetje/goodruns/pass_r1-5.txt"; ///< The location of the good run list file
        /**
         *  @brief  The Binning structure
         *          The supplied binEdges define the "analysis bins" in which the data will be presented. The min/max values give the limits
         *          of the phase-space that's includede in the measurment. If min/max extend beyond the bin edges, then an
         *          underflow/overflow bin will be added.
         */
        struct Binning
        {
            float               min;      ///< Minimum possible value
            float               max;      ///< Maximum possible value
            std::vector<float>  binEdges; ///< The analysis bin edges
        };

        // Here we define the binning for each kinematic variable. We're initializing each Binning object using member initialization lists.
        // I.e. we supply the member variable of the Binning struct, in the order they are deined (min, max, binEdges).

        Binning muonCosTheta {
            -1.f,                                                                              // min
             1.f,                                                                              // max
            {-1.f, -0.27f, 0.29f, 0.46f, 0.58f, 0.67f, 0.77f, 0.82f, 0.88f, 0.93f, 0.97f, 1.f} // binEdges
        }; ///< The muonCosTheta binning

        Binning muonPhi {
            -3.142f,                                                     // min
             3.142f,                                                     // max
            PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f) // binEdges
        }; ///< The muonPhi binning

        Binning muonMomentum {
            0.15f,                                    // min
            std::numeric_limits<float>::max(),        // max
            {0.15f, 0.23f, 0.32f, 0.45f, 0.66f, 1.5f} // binEdges
        }; ///< The muonMomentum binning

        Binning pionCosTheta {
            -1.f,                                                // min
             1.f,                                                // max
            {-1.f, -0.47f, 0.f, 0.39f, 0.65f, 0.84f, 0.93f, 1.f} // binEdges
        }; ///< The pionCosTheta binning

        Binning pionPhi {
            -3.142f,                                                     // min
             3.142f,                                                     // max
            PlottingHelper::GenerateUniformBinEdges(10, -3.142f, 3.142f) // binEdges
        }; ///< The pionPhi binning

        Binning pionMomentum {
            0.f,                               // min
            std::numeric_limits<float>::max(), // max
            {0.1f, 0.16f, 0.19f, 0.22f, 0.6f}  // binEdges
        }; ///< The pionMomentum binning

        Binning muonPionAngle {
            0.f,                                                   // min
            2.65f,                                                 // max
            {0.f, 0.49f, 0.93f, 1.26f, 1.57f, 1.88f, 2.21f, 2.65f} // binEdges
        }; ///< The muonPionAngle binning

        Binning nProtons {
            0,                                           // min
            std::numeric_limits<float>::max(),           // max
            {0, 1, 2, std::numeric_limits<float>::max()} // binEdges
        }; ///< The nProtons binning

        Binning sidebandMuonMomentum {
            0.15f,                                    // min
            std::numeric_limits<float>::max(),        // max
            {0.15f, 0.23f, 0.32f, 0.45f, 0.66f, 1.5f, std::numeric_limits<float>::max()} // binEdges
        }; ///< The muonMomentum binning without the overflow bins

        Binning sidebandPionMomentum {
            0.f,                               // min
            std::numeric_limits<float>::max(), // max
            {0.f, 0.1f, 0.16f, 0.19f, 0.22f, 0.6f, std::numeric_limits<float>::max()}
        }; ///< The pionMomentum binning without the underflow/overflow bins

        Binning sidebandProtonMomentum {
            0.f,                               // min
            std::numeric_limits<float>::max(), // max
            // {0.f, 0.32f, 0.63f, 0.93f, 1.23f, 1.38f, std::numeric_limits<float>::max()}  // binEdges //todo find optimal binning - currently just eyballed values
            // {0.f, 0.44f, 0.55f, 0.69f, 0.81f, 1.01f, std::numeric_limits<float>::max()}
            {0.f, 0.42f, 0.47f, 0.54f, 0.61f, 0.76f, std::numeric_limits<float>::max()}
        }; ///< Proton momentum binning for the pionMomentum cross-section



        // Binning muonPhiSideband {
        //     -3.142f,                                                     // min
        //      3.142f,                                                     // max
        //     PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f) // binEdges
        // }; ///< The muonPhi binning

        // Binning nProtonsSideband {
        //     0,                                           // min
        //     std::numeric_limits<float>::max(),           // max
        //     {0, 1, std::numeric_limits<float>::max()} // binEdges
        // }; ///< The nProtons binning




        // Additional plots requested by EB
        // TODO: Values are only placeholders
        Binning protonCosTheta {
            -1.f,                                                // min
             1.f,                                                // max
            {-1.f, -0.47f, 0.f, 0.39f, 0.65f, 0.84f, 0.93f, 1.f} // binEdges
        };  ///< The protonCosTheta binning

        Binning protonPhi {
            -3.142f,                                                     // min
             3.142f,                                                     // max
            PlottingHelper::GenerateUniformBinEdges(15, -3.142f, 3.142f) // binEdges
        }; ///< The protonPhi binning

        Binning protonMomentum {
            0.1f,                                    // min
            1.8f,        // max
            {0.1f, 0.23f, 0.32f, 0.45f, 0.66f, 1.3f} // binEdges
        }; ///< The protonMomentum binning

        Binning protonPionAngle {
            -0.f,                                                   // min
            3.142f,//2.65f,                                                 // max
            {0.f, 0.49f, 0.93f, 1.26f, 1.57f, 1.88f, 2.21f, 2.65f} // binEdges
        }; ///< The protonPionAngle binning

        Binning protonMuonAngle {
            -0.f,                                                   // min
            3.142f,//2.65f,                                                 // max
            {0.f, 0.49f, 0.93f, 1.26f, 1.57f, 1.88f, 2.21f, 2.65f} // binEdges
        }; ///< The protonMuonAngle binning

    };
    Global global; ///< The global configuration options

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the CountPOT macro
     */
    struct CountPOT
    {
        bool  useOverlays           = false; ///< If we should count the POT for the overlays
        bool  useDirt               = false; ///< If we should count the POT for the dirt
        bool  useDetectorVariations = true; ///< If we should count the POT for the detector variations
        bool  useNuWro              = false; ///< If we should count the POT for NuWro
    };
    CountPOT countPOT; ///< The configuration options for the CountPOT macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the GetRunSubrunList macro
     */
    struct GetRunSubrunList
    {
        bool  useDataEXT = true; ///< If we should run on the EXT data
        bool  useDataBNB = true; ///< If we should run on the BNB data
        bool  useNuWro = true;   ///< If we should run on the NuWro files
    };
    GetRunSubrunList getRunSubrunList; ///< The configuration options for the GetRunSubrunList macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the MultiPlanePIDDemo macro
     */
    struct MultiPlanePIDDemo
    {
        float sin2AngleThreshold = 0.175; ///< The squared sin angular threshold between a particle and a wire in the YZ plane to use dEdx information
    };
    MultiPlanePIDDemo multiPlanePIDDemo; ///< The configuration options for the MultiPlanePIDDemo macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the PlotInputVariables macro
     */
    struct PlotInputVariables
    {
        bool plotBDTResponses = true; ///< If we should plot the responses of the trained BDTs, set to true if you haven't already trained the BDTs
    };
    PlotInputVariables plotInputVariables; ///< The configuration options for the PlotInputVariables macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the NMinusOneBDTStudy macro
     */
    struct NMinusOneBDTStudy
    {
        bool shouldTrainBDTs = true;                                          ///< If we should run the BDT training (if false then look for trained BDTs)
        PlottingHelper::PlotStyle signalType = PlottingHelper::Muon;    ///< The type of particle considered signal by the BDT in question
        std::vector< std::string > featureNames = {
            "logBragg_pToMIP",
            "logBragg_piToMIP",
            "truncMeandEdx",
            "protonForward",
            "muonForward",
            "nDescendents",
            "nSpacePointsNearEnd",
            "wiggliness",
            "trackScore"
        };                                                                    ///< The features to consider by the BDT
        unsigned int nSamplePoints = 300u;                                    ///< The number of sampling points to use when finding the ROC curves
    };
    NMinusOneBDTStudy nMinusOneBDTStudy; ///< The configuration options for the NMinusOneBDTStudy macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for TrainBDTs macro
     */
    struct TrainBDTs
    {
        float trainingFraction     = 0.5f;  ///< Fraction of the sample on which we should train the BDTs
        bool  onlyGoodTruthMatches = false; ///< If we should only train on reco particles with a completeness >50%
        bool  weightByCompleteness = true;  ///< If we should weight the training examples by the reco-truth match completeness
        bool  shouldOptimize       = false; ///< If we should optimize the BDT parameters (this can take a while)
        bool  shouldMakePlots      = true;  ///< If we should make plots of the BDT responses
    };
    TrainBDTs trainBDTs; ///< The configuration options for the TrainBDTs macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the MakeSelectedPIDTable macro
     */
    struct MakeSelectedPIDTable
    {
        bool  useGenericSelection = false;  ///< If we should use the generic selection (if false, we use full golden selection)
        bool  goldenPionIsSignal = false;   ///< If we should only treat events containing golden pions as signal
        bool  onlyLowMomentumPions = false; ///< If we should only treat events with low momentum pions as "signal" - to check the performance in a restricted region of phase-space
        float pionMomentumThreshold = 0.1f; ///< The threshold pion momentum below which we consider "signal" if onlyLowMomentumPions == true
    };
    MakeSelectedPIDTable makeSelectedPIDTable; ///< The configuration options for the MakeSelectedPIDTable macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the EfficiencyPlots macro
     */
    struct EfficiencyPlots
    {
        bool drawErrors = true; ///< If we should draw errors
    };
    EfficiencyPlots efficiencyPlots; ///< The configuration options for the EfficiencyPlots macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the MakeBinningPlots macro
     */
    struct MakeBinningPlots
    {
        bool useFineBinEdges = true; ///< If we break up the analsis bins into finer sub-bins
    };
    MakeBinningPlots makeBinningPlots; ///< The configuration options for the MakeBinningPlots macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the ExtractXSecs macro
     */
    struct ExtractXSecs
    {
        /**
        *  @brief  A mapping to identify which cross-section should be extracted. The first index is an identifier for the selection that's
        *          used (either "generic" or "golden"). The second index is an identifier for the cross-section itself (either "total" for
        *          the total cross-section, or the name of the kinematic parameter, e.g. "muonMomentum"). The mapped type is a boolean,
        *          which is true if the cross-section should be extracted, and false otherwise.
        */
        std::unordered_map<std::string, std::unordered_map<std::string, bool> > crossSectionIsEnabled = {
            {
                "generic", {
                    {"total",         true }, // B
                    {"muonCosTheta",  true }, // A
                    {"muonPhi",       true }, // B
                    {"muonMomentum",  false }, // C
                    {"pionCosTheta",  true }, // B
                    {"pionPhi",       true }, // B
                    {"pionMomentum",  false }, // C
                    {"muonPionAngle", true }, // B
                    {"nProtons",      true } // B
                }
            },
            {
                "golden", {
                    {"total",         false },
                    {"muonCosTheta",  false },
                    {"muonPhi",       false },
                    {"muonMomentum",  false },
                    {"pionCosTheta",  false },
                    {"pionPhi",       false },
                    {"pionMomentum",  false },
                    {"muonPionAngle", false },
                    {"nProtons",      false } // C
                }
            }
        };

        /**
        *  @brief  Boolean indicating if we should scale the GENIE weights so not to double count the genieTuneEventWeight
        *          By default we apply (splineEventWeight * genieTuneEventWeight) as the "nominal weight" to all events. Then we apply the
        *          multisim universe weights on top of the nominal weight. For the GENIE systematic parameters, the universe weights already
        *          include the genieTuneEventWeight, and so require special treatment to avoid double counting this weight. If this option
        *          is set to true, then we scale the GENIE universe weights down by genieTuneEventWeight to put them on the same footing as
        *          all other parameters. If this option is false then the GENIE weights recieve no special treatment.
        */
        bool scaleXSecWeights = true;

        unsigned int nBootstrapUniverses = 1000u; ///< The number of bootrap universes to generate for the MC stat uncertainty
        unsigned int nSidebandFitUniverses = 1000u; ///< The number of universes to use for the sideband-fit parameter uncertainties

        float potFracUncertainty = 0.02f; ///< The fractional uncertainty on the POT normalisation

        CrossSectionHelper::SystDimensionsMap fluxDimensions = {
            {"hadronProduction",          1000u},
            {"expskin_FluxUnisim",        1000u},
            {"horncurrent_FluxUnisim",    1000u},
            {"nucleoninexsec_FluxUnisim", 1000u},
            {"nucleonqexsec_FluxUnisim",  1000u},
            {"nucleontotxsec_FluxUnisim", 1000u},
            {"pioninexsec_FluxUnisim",    1000u},
            {"pionqexsec_FluxUnisim",     1000u},
            {"piontotxsec_FluxUnisim",    1000u}
        }; ///< A mapping from the flux parameter names to the number of universes

        CrossSectionHelper::SystDimensionsMap xsecDimensions = {
            {"All_UBGenie",             100u},
            {"AxFFCCQEshape_UBGenie",   2u},
            {"DecayAngMEC_UBGenie",     2u},
            {"Theta_Delta2Npi_UBGenie", 2u},
            {"VecFFCCQEshape_UBGenie",  2u},
            {"xsr_scc_Fa3_SCC",         10u},
            {"xsr_scc_Fv3_SCC",         10u}
        }; ///< A mapping from the cross-section parameter names to the number of universes

        std::unordered_map<std::string, bool> xsecUBGenieScalingMap = {
            {"All_UBGenie",             true},
            {"AxFFCCQEshape_UBGenie",   true},
            {"DecayAngMEC_UBGenie",     true},
            {"Theta_Delta2Npi_UBGenie", true},
            {"VecFFCCQEshape_UBGenie",  true},
            {"xsr_scc_Fa3_SCC",         false},
            {"xsr_scc_Fv3_SCC",         false}
        }; ///< Mapping from cross-section parameter names to boolean indicating if we should scale down the parameters by the genieTuneEventWeight ...
        // ... to avoid double counting the factor in the universe and nominal weights

        CrossSectionHelper::SystDimensionsMap reintDimensions = {
            {"reinteractions_piminus_Geant4", 1000u},
            {"reinteractions_piplus_Geant4",  1000u},
            {"reinteractions_proton_Geant4",  1000u}
        }; ///< A mapping from the reinteraction parameter names to the number of universes

        CrossSectionHelper::SystUnisimDimensionsMap detVarDimensions = {
            {"LYDown",         "CV"},
            {"LYRayleigh",     "CV"},
            {"SCE",            "CV_SCE"},
            {"Recomb2",        "CV_SCE"},
            {"WireModX",       "CV"},
            {"WireModYZ",      "CV"},
            {"WireModThetaXZ", "CV"},
            {"WireModThetaYZ", "CV"}
        }; ///< A mapping from the detector variation sample identifiers, to the identifiers for their relevant central-value sample

        // CrossSectionHelper::SystUnisimDimensionsMap detVarDimensions = {
        //     {"LYDown",         "CVRun1"},
        //     {"LYRayleigh",     "CVRun1"},
        //     {"SCE",            "CVRun3b"},
        //     {"Recomb2",        "CVRun3b"},
        //     // {"WireModX",       "CVRun3b"},
        //     // {"WireModYZ",      "CVRun3b"},
        //     // {"WireModThetaXZ", "CVRun3b"},
        //     // {"WireModThetaYZ", "CVRun3b"}
        // }; ///< A mapping from the detector variation sample identifiers, to the identifiers for their relevant central-value sample

        /**
        *  @brief  A mapping from a user-defined parameter name, to corresponding set of mutually exclusive parameter names and the number of universes
        *
        *          A set of parameters is mutually exclusive if for any given universe at most one parameter has a weight that's not equal
        *          to one. In this case, the set of parameters can be "combined" into one parameter (with the user-defined name). The value
        *          of the combined parameter is taken from whichever parameter in the set has been varied(or equivalently just the product
        *          of the weights in the input parameter set).
        *
        *          For example, the flux hadron production channel can be any one of kminus, kplus, kzero, piminus, piplus, but never more
        *          than one channel at once. As a result, a given event will have a (possibly) non-unit weight for the relevant channel, and
        *          a weight of one for all other channels. In this sense, there is a "correlation" between the weights for each channel.
        *          Instead of applying each channel as a separate systematic parameter, we apply their product as a combined parameter which
        *          accounts for their "correlations"
        */
        CrossSectionHelper::SystMutuallyExclusiveDimensionsMap mutuallyExclusiveDimensions = {
            {
                "hadronProduction",
                {
                    {
                        "kminus_PrimaryHadronNormalization",
                        "kplus_PrimaryHadronFeynmanScaling",
                        "kzero_PrimaryHadronSanfordWang",
                        "piminus_PrimaryHadronSWCentralSplineVariation",
                        "piplus_PrimaryHadronSWCentralSplineVariation"
                    }, 1000u
                }
            }
        };

    };
    ExtractXSecs extractXSecs; ///< The configuration options for the ExtractXSecs macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the PrintUncertaintiesSummary macro
     */
    struct PrintUncertaintiesSummary
    {
        bool  useGenericSelection = true;  ///< If we should use the generic selection (if false, we use golden selection)
    };
    PrintUncertaintiesSummary printUncertaintiesSummary; ///< The configuration options for the PrintUncertaintiesSummary macro

    // -------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Configuration for the MakeXSecPlots macro
     */
    struct MakeXSecPlots
    {
        unsigned int nUniverses = 10000u;                                 ///< The number of universes to use when propagating the smearing matrix uncertainties to the smeared prediction
        float        precision  = std::numeric_limits<float>::epsilon();  ///< The precision to use when finding eigenvalues & eigenvectors
    };
    MakeXSecPlots makeXSecPlots; ///< The configuration options for the MakeXSecPlots macro
};

} // namespace ubcc1pi

#endif
