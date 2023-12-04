/**
 *  @file  ubcc1pi_standalone/Interface/EventXSec.h
 *
 *  @brief The header file for the event class
 */

#ifndef UBCC1PI_STANDALONE_INTERFACE_EVENT_XSEC
#define UBCC1PI_STANDALONE_INTERFACE_EVENT_XSEC

#include "ubcc1pi_standalone/Interface/EventMembers.h"
#include "ubcc1pi_standalone/Interface/Member.h"

#include <iostream>
#include <stdexcept>
#include <limits>
#include <string>
#include <vector>
#include <memory>

#include <TTree.h>

namespace ubcc1pi
{

class FileWriter;
// template <typename T> class FileReader;
class EventFactory;

/**
 *  @brief  The event class
 */
class EventXSec
{
    public:

        /**
         *  @brief  Constructor
         *  @param  hasTruthInfo whether the event has truth info
         */
        EventXSec(const bool hasTruthInfo = false);

        /**
         *  @brief  Print the member variables to the terminal
         */
        void Print() const;

        /**
         *  @brief  Return whether the event contains truth weights
         */
        bool HasTruthWeights() const; 

        /**
         *  @brief  The metadata information structure
         */
        // struct Metadata
        // {
        //     XSEC_MACRO_EVENT_METADATA_MEMBERS("", "", PELEE_MACRO_DECLARE_MEMBER)
        // };

        /**
         *  @brief  The truth information structure
         */
        struct Truth
        {
            XSEC_MACRO_EVENT_TRUTH_MEMBERS("", "", PELEE_MACRO_DECLARE_MEMBER)
            // XSEC_MACRO_EVENT_TRUTH_OPTIONAL_MEMBERS("", "", XSEC_MACRO_DECLARE_MEMBER)

            // /**
            //  *  @brief  The truth particle information structure
            //  */
            // struct Particle
            // {
            //     // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", "", XSEC_MACRO_DECLARE_MEMBER)
            // };

            // std::vector<Particle> particles; ///< The truth particles
        };

        /**
         *  @brief  The reco information structure
         */
        struct Reco
        {
            XSEC_MACRO_EVENT_RECO_MEMBERS("", "", PELEE_MACRO_DECLARE_MEMBER)

            // /**
            //  *  @brief  The reco particle information structure
            //  */
            // struct Particle
            // {
            //     // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS("", "", XSEC_MACRO_DECLARE_MEMBER)
            // };

            // std::vector<Particle> particles; ///< The reco particles
        };

        // struct Backtracked
        // {
        //     // XSEC_MACRO_EVENT_BACKTRACKED_MEMBERS("", "", XSEC_MACRO_DECLARE_MEMBER)

        //     /**
        //      *  @brief  The backtracked particle information structure
        //      */
        //     struct Particle
        //     {
        //         // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS("", "", XSEC_MACRO_DECLARE_MEMBER)
        //     };

        //     std::vector<Particle> particles; ///< The backtracked particles
        // };

        // Metadata metadata;  ///< The metadata
        Truth truth;        ///< The truth information
        Reco  reco;         ///< The reco information
        // Backtracked backtracked; ///< The backtracked information

    private:

        friend FileWriter;              ///< The file writer class is a friend
        // friend FileReader<EventXSec>;  ///< The file reader class is a friend
        template<typename T, typename U> friend class FileReader;
        friend EventFactory;            ///< The event factory class is a friend
        bool hasTruthWeights;

        /**
         *  @brief  Bind this event to an output tree
         *
         *  @param  pTree the tree with which to bind
         */
        void BindToOutputTree(TTree * pTree);

        /**
         *  @brief  Bind this event to an input tree
         *
         *  @param  pTree the tree with which to bind
         */
        void BindToInputTree(TTree * pTree);

        /**
         *  @brief  Reset the member variables as if the event were new
         */
        void Reset();

        /**
         *  @brief  Fill the output vectors with the data from the particles ready for a fill
         */
        void PrepareForTreeFill();

        /**
         *  @brief  Fill the input particles with the info from the vectors in the tree
         */
        void PrepareAfterTreeRead();

        // // Here we define private member variables for each of the particle parameters as a vector so they can be read from the root file
        // XSEC_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", XSEC_MACRO_DECLARE_MEMBER_VECTOR)
        // XSEC_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", XSEC_MACRO_DECLARE_MEMBER_VECTOR)
        // XSEC_MACRO_EVENT_BACKTRACKED_PARTICLE_MEMBERS(backtracked_particle, "", XSEC_MACRO_DECLARE_MEMBER_VECTOR)
};

} // namespace ubcc1pi

#endif
