/**
 *  @file  ubcc1pi_standalone/Interface/SubrunPeLEE.h
 *
 *  @brief The header file for the subrun PeLEE class
 */

#ifndef UBCC1PI_STANDALONE_INTERFACE_SUBRUN_PELEE
#define UBCC1PI_STANDALONE_INTERFACE_SUBRUN_PELEE

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
class SubrunPeLEE
{
    public:

        /**
         *  @brief  Print the member variables to the terminal
         */
        void Print() const;

        PELEE_MACRO_SUBRUN_MEMBERS("", "", PELEE_MACRO_DECLARE_MEMBER)

    private:

        friend FileWriter;      ///< The file writer class is a friend
        // friend FileReader<Event>;      ///< The file reader class is a friend
        template<typename T, typename U> friend class FileReader;
        friend SubrunFactory;   ///< The subrun factory class is a friend

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
