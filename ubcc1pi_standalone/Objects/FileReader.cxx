/**
 *  @file  ubcc1pi_standalone/Objects/FileReader.cxx
 *
 *  @brief The implementation of the file reader class
 */

#include "ubcc1pi_standalone/Objects/FileReader.h"

namespace ubcc1pi
{

template <class T, class U>
FileReader<T, U>::FileReader(const std::string &inputFile, const bool hasTruthInfo)
{
    std::cout << "DEBUG: Entering FileReader constructor" << std::endl;

    m_inputFile = inputFile;
    std::cout << "DEBUG: Input file set to: " << m_inputFile << std::endl;

    m_pFile = TFile::Open(m_inputFile.c_str());
    std::cout << "DEBUG: File opened: " << m_inputFile << std::endl;

    std::string eventTreeName = "events";
    if constexpr (std::is_same_v<T, EventPeLEE>) eventTreeName = "nuselection/NeutrinoSelectionFilter";
    else if constexpr (std::is_same_v<T, EventXSec>) eventTreeName = "stv_tree";
    std::cout << "DEBUG: Event tree name set to: " << eventTreeName << std::endl;

    m_pEventTree = static_cast<TTree*>(m_pFile->Get(eventTreeName.c_str()));
    std::cout << "DEBUG: Event tree obtained from file" << std::endl;

    std::string subrunTreeName = "subruns";
    if constexpr (std::is_same_v<U, SubrunPeLEE>) subrunTreeName = "nuselection/SubRun";
    else if constexpr (std::is_same_v<U, SubrunXSec>) subrunTreeName = "stv_tree/nuselection/SubRun";
    std::cout << "DEBUG: Subrun tree name set to: " << subrunTreeName << std::endl;

    m_pSubrunTree = static_cast<TTree*>(m_pFile->Get(subrunTreeName.c_str()));
    std::cout << "DEBUG: Subrun tree obtained from file" << std::endl;

    std::cout << "DEBUG: Event and Subrun trees initialized" << std::endl;

    m_pEvent = std::make_shared<T>(hasTruthInfo);
    m_pSubrun = std::make_shared<U>(hasTruthInfo);

    std::cout << "DEBUG: Event and Subrun objects created" << std::endl;

    this->BindEventToTree();
    this->BindSubrunToTree();

    std::cout << "DEBUG: Event and Subrun objects bound to trees" << std::endl;

    if (this->GetNumberOfEvents() == 0)
        std::cout << "ubcc1pi::FileReader: Warning, input file has no entries" << std::endl;

    std::cout << "DEBUG: Number of events: " << this->GetNumberOfEvents() << std::endl;

    if (this->GetNumberOfSubruns() == 0)
        std::cout << "ubcc1pi::FileReader: Warning, input file has no subruns" << std::endl;

    std::cout << "DEBUG: Number of subruns: " << this->GetNumberOfSubruns() << std::endl;

    // By default, disable the branches that hold the event weights as these are large and only needed in a few cases
    this->DisableSystematicBranches();

    std::cout << "DEBUG: Exiting FileReader constructor" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
FileReader<T, U>::~FileReader()
{
    // ATTN this function does the memory management for the pointers belonging to this class
    m_pFile->Close();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::BindEventToTree()
{
    m_pEvent->BindToInputTree(m_pEventTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::BindSubrunToTree()
{
    m_pSubrun->BindToInputTree(m_pSubrunTree);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
std::shared_ptr<T> FileReader<T, U>::GetBoundEventAddress()
{
    return m_pEvent;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
std::shared_ptr<U> FileReader<T, U>::GetBoundSubrunAddress()
{
    return m_pSubrun;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::DisableSystematicBranches()
{
    std::cout << "DEBUG: Disabling systematic branches" << std::endl;
    std::string weightBranches = "truth_systParam*";
    if constexpr (std::is_same_v<T, EventPeLEE>) weightBranches = "weights*";
    if constexpr (!std::is_same_v<T, EventXSec>) m_pEventTree->SetBranchStatus(weightBranches.c_str(), false);
    std::cout << "DEBUG: Disabling systematic branches" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::EnableSystematicBranches()
{
    std::cout << "DEBUG: Enabling systematic branches" << std::endl;
    std::string weightBranches = "truth_systParam*";
    if constexpr (std::is_same_v<T, EventPeLEE>) weightBranches = "weights*";
    if constexpr (!std::is_same_v<T, EventXSec>) m_pEventTree->SetBranchStatus(weightBranches.c_str(), true);
    std::cout << "DEBUG: Enabling systematic branches" << std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
unsigned int FileReader<T, U>::GetNumberOfEvents() const
{
    return m_pEventTree->GetEntries();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
unsigned int FileReader<T, U>::GetNumberOfSubruns() const
{
    return std::is_same_v<T, EventXSec> ? 0 : m_pSubrunTree->GetEntries();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::LoadEvent(const unsigned int eventIndex)
{
    // std::cout<<"DEBUG LoadEvent 0"<<std::endl;
    const auto nEvents = this->GetNumberOfEvents();
    // std::cout<<"DEBUG LoadEvent 1"<<std::endl;
    if (eventIndex >= nEvents)
        throw std::range_error("Can't load event, input eventIndex is out of bounds");

    // std::cout<<"DEBUG LoadEvent 2"<<std::endl;
    m_pEvent->Reset();
    // std::cout<<"DEBUG LoadEvent 3"<<std::endl;
    m_pEventTree->GetEntry(eventIndex);
    // std::cout<<"DEBUG LoadEvent 4"<<std::endl;
    m_pEvent->PrepareAfterTreeRead();
    // std::cout<<"DEBUG LoadEvent 5"<<std::endl;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <class T, class U>
void FileReader<T, U>::LoadSubrun(const unsigned int subrunIndex)
{
    const auto nSubruns = this->GetNumberOfSubruns();
    if (subrunIndex >= nSubruns)
        throw std::range_error("Can't load subrun, input subrunIndex is out of bounds");

    m_pSubrun->Reset();
    m_pSubrunTree->GetEntry(subrunIndex);
}

// -----------------------------------------------------------------------------------------------------------------------------------------

// Explicit instantiations
template class FileReader<Event, Subrun>;
template class FileReader<EventPeLEE, SubrunPeLEE>;
template class FileReader<EventXSec, SubrunXSec>;

} // namespace ubcc1pi
