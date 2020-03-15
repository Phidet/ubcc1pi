/**
 *  @file  ubcc1pi_standalone/Interface/Event.h
 *
 *  @brief The header file for the event class
 */

#ifndef UBCC1PI_STANDALONE_INTERFACE_EVENT
#define UBCC1PI_STANDALONE_INTERFACE_EVENT

#include "ubcc1pi_standalone/Interface/EventMembers.h"

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
class FileReader;
class EventFactory;

class Event
{
    public:

        /**
         *  @brief  Constructor
         */
        Event();

        /**
         *  @brief  Print the member variables to the terminal
         */
        void Print() const;

        /**
         *  @brief  The class for member variable which will throw if not set
         */
        template <typename T>
        class Member
        {
            public:
                /**
                 *  @brief  Default constructor
                 */
                Member();
                
                /**
                 *  @brief  Returns if this member has been set
                 *
                 *  @return true if set, false otherwise
                 */
                bool IsSet() const;

                /**
                 *  @brief  Gets the value of the member if set and throws otherwise
                 *
                 *  @return the value
                 */
                const T Get() const;

                /**
                 *  @brief  Overload the () operator as a shorthand Get()
                 */
                const T operator()() const;

            private:

                // Allow the event class to access private members of this class, all other classes will have to go through the public
                // accessor functions Get() and Set()
                friend class Event;
                friend class EventFactory;
                
                /**
                 *  @brief  Set the value of the member
                 *
                 *  @param  value the value to set
                 */
                void Set(const T &value);
                
                /**
                 *  @brief  Reset the member value to default
                 */
                void Reset();
                
 
                /**
                 *  @brief  Set the input value to default
                 *
                 *  @param  value the value to set to default
                 */
                void SetDefault();
                
                /**
                 *  @brief  Convert the member variable to a string
                 *
                 *  @return the string
                 */
                std::string ToString() const;

                std::shared_ptr<T>   m_pValue;        ///< The value of the member
                T                   *m_pAddress;      ///< The address owned by the shared_ptr needed by root... sigh
                bool                 m_isSet;         ///< If the value has been set
        };

        /**
         *  @brief  The metadata information structure
         */
        struct Metadata
        {
            UBCC1PI_MACRO_EVENT_METADATA_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)
        };

        /**
         *  @brief  The truth information structure
         */
        struct Truth
        {
            UBCC1PI_MACRO_EVENT_TRUTH_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)

            /**
             *  @brief  The truth particle information structure
             */
            struct Particle
            {
                UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)
            };

            std::vector<Particle> particles; ///< The truth particles
        };
        
        /**
         *  @brief  The reco information structure
         */
        struct Reco
        {
            UBCC1PI_MACRO_EVENT_RECO_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)

            /**
             *  @brief  The reco particle information structure
             */
            struct Particle
            {
                UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS("", "", UBCC1PI_MACRO_DECLARE_MEMBER)
            };

            std::vector<Particle> particles; ///< The reco particles
        };

        Metadata metadata;  ///< The metadata
        Truth truth;        ///< The truth information
        Reco  reco;         ///< The reco information

    private:

        friend FileWriter;
        friend FileReader;
        friend EventFactory;

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

        // Here we define private member variables for each of the particle parameters as a vector so they can be read from the root file
        UBCC1PI_MACRO_EVENT_TRUTH_PARTICLE_MEMBERS(truth_particle, "", UBCC1PI_MACRO_DECLARE_MEMBER_VECTOR)
        UBCC1PI_MACRO_EVENT_RECO_PARTICLE_MEMBERS(reco_particle, "", UBCC1PI_MACRO_DECLARE_MEMBER_VECTOR)
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
Event::Member<T>::Member() :
    m_pValue(std::make_shared<T>()),
    m_pAddress(m_pValue.get())
{
    this->Reset();
}   

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline bool Event::Member<T>::IsSet() const
{
    return m_isSet;
}

// -----------------------------------------------------------------------------------------------------------------------------------------
        
template <typename T>
inline void Event::Member<T>::Set(const T &value)
{
    *m_pValue = value;
    m_isSet = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const T Event::Member<T>::Get() const
{
    if (this->IsSet())
        return *m_pValue;

    throw std::logic_error("You are trying to access a member value of ubcc1pi::Event that hasn't been set");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>        
inline const T Event::Member<T>::operator()() const
{
    return this->Get();
}
        
// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void Event::Member<T>::Reset()
{
    this->SetDefault();
    m_isSet = false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member<bool>::SetDefault()
{
    *m_pValue =  false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member<int>::SetDefault()
{
    *m_pValue =  -std::numeric_limits<int>::max();
}
// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member<float>::SetDefault()
{
    *m_pValue =  -std::numeric_limits<float>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member<std::string>::SetDefault()
{
    *m_pValue =  "";
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member<TVector3>::SetDefault()
{
    m_pValue->SetXYZ(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member< std::vector<float> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member< std::vector<bool> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Event::Member< std::vector<int> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member<bool>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return this->Get() ? "true" : "false";
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member<int>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return std::to_string(this->Get());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member<float>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return std::to_string(this->Get());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member<std::string>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return this->Get();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member<TVector3>::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    return ("(" + std::to_string(value.X()) + ", " + std::to_string(value.Y()) + ", " + std::to_string(value.Z()) + ")");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member< std::vector<float> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";
    for (const auto &entry : value)
        str += std::to_string(entry) + "  ";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member< std::vector<bool> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";
    for (const auto &entry : value)
        str += std::string(entry ? "true" : "false") + "  ";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline std::string Event::Member< std::vector<int> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";
    for (const auto &entry : value)
        str += std::to_string(entry) + "  ";

    return str;
}

} // namespace ubcc1pi

#endif
