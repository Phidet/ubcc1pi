/**
 *  @file  ubcc1pi_standalone/Objects/FileReader.h
 *
 *  @brief The header file for the event file reader class
 */

#ifndef UBCC1PI_STANDALONE_OBJECTS_FILE_READER
#define UBCC1PI_STANDALONE_OBJECTS_FILE_READER

#include "ubcc1pi_standalone/Interface/Event.h"
#include "ubcc1pi_standalone/Interface/EventPeLEE.h"
#include "ubcc1pi_standalone/Interface/EventXSec.h"
#include "ubcc1pi_standalone/Interface/Subrun.h"
#include "ubcc1pi_standalone/Interface/SubrunPeLEE.h"
#include "ubcc1pi_standalone/Interface/SubrunXSec.h"

#include <memory>
#include <TFile.h>
#include <TTree.h>

namespace ubcc1pi
{

/**
 *  @brief  The File Reader class
 */
template <class T = Event, class U = Subrun> // T: Event or EventPeLEE class
class FileReader
{
    public:

        /**
         *  @brief  Constructor
         *
         *  @param  inputFile the path to the input file
         *  @param  hasTruthInfo whether the event has truth info
         */
        FileReader(const std::string &inputFile, const bool hasTruthInfo);

        /**
         *  @brief  Destructor
         */
        ~FileReader();

        /**
         *  @brief  Get the total number of events
         *
         *  @return number of events
         */
        unsigned int GetNumberOfEvents() const;

        /**
         *  @brief  Get the total number of subruns
         *
         *  @return number of subruns
         */
        unsigned int GetNumberOfSubruns() const;

        // TODO make these bound output addresses point to const objects

        /**
         *  @brief  Get the event bound to the output tree
         *
         *  @return the address of the bound event
         */
        std::shared_ptr<T> GetBoundEventAddress();

        /**
         *  @brief  Get the subrun bound to the output tree
         *
         *  @return the address of the bound subrun
         */
        std::shared_ptr<U> GetBoundSubrunAddress();

        /**
         *  @brief  Disable the branches related to systematic weights
         */
        void DisableSystematicBranches();

        /**
         *  @brief  Enable the branches related to systematic weights
         */
        void EnableSystematicBranches();

        /**
         *  @brief  Load the event with the supplied index
         *
         *  @param  eventIndex the event to load
         */
        void LoadEvent(const unsigned int eventIndex);

        /**
         *  @brief  Load the subrun with the supplied index
         *
         *  @param  subrunIndex the subrun to load
         */
        void LoadSubrun(const unsigned int subrunIndex);

    private:

        /**
         *  @brief  Bind the event member to the input tree
         */
        void BindEventToTree();

        /**
         *  @brief  Bind the subrun member to the input tree
         */
        void BindSubrunToTree();

        std::string m_inputFile;  ///< The input file name

        TFile      *m_pFile;       ///< The input file
        TTree      *m_pEventTree;  ///< The input event tree
        TTree      *m_pSubrunTree; ///< The input subrun tree

        std::shared_ptr<T>           m_pEvent;       ///< The input event
        std::shared_ptr<U>           m_pSubrun;      ///< The input subrun
};

} // namespace ubcc1pi

#endif
