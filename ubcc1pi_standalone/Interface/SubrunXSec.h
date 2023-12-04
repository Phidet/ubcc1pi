/**
 *  @file  ubcc1pi_standalone/Interface/SubrunXSec.h
 *
 *  @brief The header file for the subrun XSec class
 */

#ifndef UBCC1PI_STANDALONE_INTERFACE_SUBRUN_XSEC
#define UBCC1PI_STANDALONE_INTERFACE_SUBRUN_XSEC

#include "ubcc1pi_standalone/Interface/SubrunMembers.h"
#include "ubcc1pi_standalone/Interface/Member.h"

#include <TTree.h>

namespace ubcc1pi
{

class FileWriter;
// template <typename T> class FileReader;
class SubrunFactory;

/**
 *  @brief  The subrun class
 */
class SubrunXSec
{
    public:

        /**
         *  @brief  Constructor
         *  @param  hasTruthInfo whether the event has truth info; needed to be compatible with EventXSec in FileReader
         */
        SubrunXSec(const bool hasTruthInfo = false);

        /**
         *  @brief  Print the member variables to the terminal
         */
        void Print() const;

        /**
         *  @brief  Return whether the event contains truth weights
         */
        bool HasTruthWeights() const; 

        // XSEC_MACRO_SUBRUN_MEMBERS("", "", PELEE_MACRO_DECLARE_MEMBER)
        // XSEC_MACRO_SUBRUN_OPTIONAL_MEMBERS("", "", XSEC_MACRO_DECLARE_MEMBER)
        

    private:

        friend FileWriter;      ///< The file writer class is a friend
        // friend FileReader<Event>;      ///< The file reader class is a friend
        template<typename T, typename U> friend class FileReader;
        friend SubrunFactory;   ///< The subrun factory class is a friend
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
};

} // namespace ubcc1pi

#endif
