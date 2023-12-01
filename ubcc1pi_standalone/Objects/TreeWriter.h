/**
 *  @file  ubcc1pi_standalone/Objects/TreeWriter.h
 *
 *  @brief The header file for the event file reader class
 *
 *  Adapted from https://github.com/sjgardiner/stv-analysis-new/blob/master/analyzer.C
 *  & https://github.com/sjgardiner/stv-analysis-new/blob/4658cb1d4163ebc4f0378a1027e2f91629eefe25/TreeUtils.hh
 *
 */

#ifndef UBCC1PI_STANDALONE_OBJECTS_TREE_WRITER
#define UBCC1PI_STANDALONE_OBJECTS_TREE_WRITER

// #include "ubcc1pi_standalone/Interface/Event.h"
// #include "ubcc1pi_standalone/Interface/EventPeLEE.h"
// #include "ubcc1pi_standalone/Interface/Subrun.h"
// #include "ubcc1pi_standalone/Interface/SubrunPeLEE.h"

// #include <memory>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TParameter.h>

namespace ubcc1pi
{

/**
 *  @brief  The File Reader class
 */
class TreeWriter
{
    public:

        /**
         *  @brief  Constructor
         *
         *  @param  outputFileDir the output file path
         */
        TreeWriter(const std::vector<std::string> inputFiles, const std::string outputFile);

        /**
         *  @brief  Destructor
         */
        ~TreeWriter();

        // /**
        //  *  @brief Custom pointer to help with filling TTree branches
        //  */
        // template <typename T>
        // class Pointer;

        // /**
        //  *  @brief  Set the event output branch addresses
        //  *
        //  *  @param  pEventPeLEE the event object
        //  *  @param  create whether to create the branches
        //  */
        // void SetEventOutputBranchAddresses(const std::shared_ptr<EventPeLEE> pEventPeLEE);

        /**
         *  @brief Stop creating branches when using SetOutputBranchAddress & SetObjectOutputBranchAddress
         */
        void CreateNoNewBranches();

        /**
         * @brief  Fill the tree
        */
        void Fill();

        /**
         * @brief  Set the output branch addresses
         *
         * @param  branchName the name of the branch
         * @param  address the address of the branch
         * @param  leafSpec the leaf specification
         */
        void SetOutputBranchAddress(const std::string& branchName, void* address, const std::string& leafSpec = "");

        /**
         * @brief  Set the output branch addresses for objects
         *
         * @param  branchName the name of the branch
         * @param  address the address of the branch
         */
        template <typename T>
        void SetObjectOutputBranchAddress(const std::string& branchName, const T*& address);

        // /**
        //  * @brief  Set the output branch addresses for objects
        //  *
        //  * @param  branchName the name of the branch
        //  * @param  address the address of the branch
        //  */
        // template <typename T>
        // void SetObjectOutputBranchAddress(const std::string& branchName, T* address);

    private:

        TChain* m_pEventChain;  ///< The event chain
        TChain* m_pSubrunChain; ///< The subrun chain
        TTree* m_pOutTree;  ///< The output event tree
        TFile* m_pOutFile;  ///< The output file
        TParameter<float>* m_pTotalExposurePOTParam; ///< The output event weight
        bool m_createdOutputBranches = false; ///< Whether the output branches have been created
};

} // namespace ubcc1pi

#endif
