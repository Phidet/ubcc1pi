/**
 *  @file  ubcc1pi_standalone/Interface/Member.h
 *
 *  @brief The header file for the member class
 */

#ifndef UBCC1PI_STANDALONE_INTERFACE_MEMBER
#define UBCC1PI_STANDALONE_INTERFACE_MEMBER

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <TVector3.h>

namespace ubcc1pi
{

class Event;
class EventPeLEE;
class EventXSec;
class Subrun;
class SubrunPeLEE;
class SubrunXSec;
class EventFactory;
class SubrunFactory;

/**
 *  @brief  The class for member variable which will throw if not set
 */
template <typename T>
class Member
{
    public:
        /**
         *  @brief  Default constructor
         *
         *  @param  name the name of the member
         */
        Member(const std::string &name = "");

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
         *  @brief  Gets the address of the value of the member if set and throws otherwise
         *
         *  @return the value
         */
        const T* GetAddress() const;

        /**
         *  @brief  Overload the () operator as a shorthand Get()
         */
        const T operator()() const;

        /**
         *  @brief  Set the value of the member
         *
         *  @param  value the value to set
         */
        void Set(const T &value);

    private:

        friend class Event;
        friend class EventPeLEE;
        friend class EventXSec;
        friend class EventFactory;

        friend class Subrun;
        friend class SubrunPeLEE;
        friend class SubrunXSec;
        friend class SubrunFactory;

        /**
         *  @brief  Reset the member value to default
         */
        void Reset();


        /**
         *  @brief  Set the input value to default
         */
        void SetDefault();

        /**
         *  @brief  Convert the member variable to a string
         *
         *  @return the string
         */
        std::string ToString() const;

        std::shared_ptr<T>   m_pValue;         ///< The value of the member
        T*                   m_pAddress;       ///< The address owned by the shared_ptr needed by root... sigh
        bool                 m_isSet;          ///< If the value has been set
        static unsigned int  m_maxVectorPrint; ///< The maximum number of entries of a vector to print
        std::string          m_name;           ///< The name of the member
};

// -----------------------------------------------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline unsigned int Member<T>::m_maxVectorPrint = 3;

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
Member<T>::Member(const std::string &name) :
    m_pValue(std::make_shared<T>()),
    m_pAddress(m_pValue.get()),
    m_name(name)
{
    this->Reset();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline bool Member<T>::IsSet() const
{
    return m_isSet;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void Member<T>::Set(const T &value)
{
    *m_pValue = value;
    m_isSet = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <>
inline void Member<unsigned short>::Set(const unsigned short &value)
{
    *m_pValue = static_cast<int>(value);
    m_isSet = true;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const T Member<T>::Get() const
{
    if (this->IsSet())
        return *m_pValue;

    throw std::logic_error("Member::Get() - You are trying to access a member value: \"" + m_name + "\" that hasn't been set");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const T* Member<T>::GetAddress() const
{
    if (this->IsSet())
        return m_pAddress;

    throw std::logic_error("Member::GetAddress() - You are trying to access a member value: \"" + m_name + "\" that hasn't been set");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const T Member<T>::operator()() const
{
    return this->Get();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void Member<T>::Reset()
{
    this->SetDefault();
    m_isSet = false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member<bool>::SetDefault()
{
    *m_pValue =  false;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member<int>::SetDefault()
{
    *m_pValue =  -std::numeric_limits<int>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member<unsigned short>::SetDefault()
{
    *m_pValue = 0u;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member<float>::SetDefault()
{
    *m_pValue =  -std::numeric_limits<float>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member<double>::SetDefault()
{
    *m_pValue =  -std::numeric_limits<double>::max();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member<std::string>::SetDefault()
{
    *m_pValue =  "";
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member<TVector3>::SetDefault()
{
    m_pValue->SetXYZ(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member< std::vector<float> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member< std::vector<double> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member< std::vector<bool> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member< std::vector<int> >::SetDefault()
{
    m_pValue->clear();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member< std::vector<unsigned short> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member< std::map<std::string, std::vector<double>> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/** @brief doxygen wants this :( */
/**
 *  @brief  Set the member variable to default
 */
template <>
inline void Member< std::vector<std::string> >::SetDefault()
{
    m_pValue->clear();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member<bool>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return this->Get() ? "true" : "false";
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member<int>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return std::to_string(this->Get());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member<float>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return std::to_string(this->Get());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member<double>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return std::to_string(this->Get());
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member<std::string>::ToString() const
{
    if (!this->IsSet())
        return "?";

    return this->Get();
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member<TVector3>::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    return ("(" + std::to_string(value.X()) + ", " + std::to_string(value.Y()) + ", " + std::to_string(value.Z()) + ")");
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member< std::vector<float> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";

    unsigned int i = 0;
    for (const auto &entry : value)
    {
        str += std::to_string(entry) + "  ";

        if (++i >= m_maxVectorPrint)
            break;
    }

    if (i < value.size())
        str += "...";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member< std::vector<double> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";

    unsigned int i = 0;
    for (const auto &entry : value)
    {
        str += std::to_string(entry) + "  ";

        if (++i >= m_maxVectorPrint)
            break;
    }

    if (i < value.size())
        str += "...";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member< std::vector<bool> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";

    unsigned int i = 0;
    for (const auto &entry : value)
    {
        str += std::string(entry ? "true" : "false") + "  ";

        if (++i >= m_maxVectorPrint)
            break;
    }

    if (i < value.size())
        str += "...";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member< std::vector<int> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";

    unsigned int i = 0;
    for (const auto &entry : value)
    {
        str += std::to_string(entry) + "  ";

        if (++i >= m_maxVectorPrint)
            break;
    }

    if (i < value.size())
        str += "...";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member< std::vector<unsigned short> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";

    unsigned int i = 0;
    for (const auto &entry : value)
    {
        str += std::to_string(entry) + "  ";

        if (++i >= m_maxVectorPrint)
            break;
    }

    if (i < value.size())
        str += "...";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member< std::vector<std::string> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";

    unsigned int i = 0;
    for (const auto &entry : value)
    {
        str += entry + "  ";

        if (++i >= m_maxVectorPrint)
            break;
    }

    if (i < value.size())
        str += "...";

    return str;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get the member variable as a string
 *
 *  @return the member variable string
 */
template <>
inline std::string Member< std::map<std::string, std::vector<double>> >::ToString() const
{
    if (!this->IsSet())
        return "?";

    const auto value = this->Get();
    std::string str = "[" + std::to_string(value.size()) + "]  ";

    unsigned int i = 0;
    for (const auto &[key, entry] : value)
    {
        str += key + " : [";
        for ( unsigned int j = 0; j < std::min(m_maxVectorPrint, (unsigned int)entry.size()); j++ )
        {
            str += std::to_string(entry[j]) + "  ";
        }
        str += "]";

        if (++i >= m_maxVectorPrint)
            break;
    }

    if (i < value.size())
        str += "...";

    return str;
}

// =========================================================================================================================================

// The macros that come below are all envoked with one of the lists of member variables (either for an event or a subrun). These lists give
// the three parameters (r, t, n) for each member variable, and the options (p, q) are supplied at the point the macros are called.
//     p = a prefix to the variable name in the list (e.g. "metadata" or "truth_particle")
//     q = an object that owns the member variable (e.g. the metadata object)
//     r = if root needs to use the address of a pointer when setting branch addresses (used to output non-standard types)
//     t = the type of the member variable (e.g. int)
//     n = the name of the member variable (e.g. run)

/** Define a macro that declares a member variable and an associated boolean to check if that variable has been set */
#define UBCC1PI_MACRO_DECLARE_MEMBER(p, q, r, t, n)                                                                                        \
    Member<t> n = Member<t>(#n);

/** Define a macro that declares a member vector, used to bind particles to vectors in trees */
#define UBCC1PI_MACRO_DECLARE_MEMBER_VECTOR(p, q, r, t, n)                                                                                 \
    std::shared_ptr< std::vector<t> > p##_##n##_vect;                                                                                      \
    std::shared_ptr< std::vector<bool> > p##_##n##_isSet_vect;                                                                             \
    std::vector<t> * p##_##n##_inputVect;                                                                                                  \
    std::vector<bool> * p##_##n##_isSet_inputVect;

/** Define a macro that initializes the member vectors */
#define UBCC1PI_MACRO_INIT_MEMBER_VECTOR(p, q, r, t, n)                                                                                    \
    p##_##n##_vect = std::make_shared< std::vector<t> >();                                                                                 \
    p##_##n##_isSet_vect = std::make_shared< std::vector<bool> >();                                                                        \
    p##_##n##_inputVect = nullptr;                                                                                                         \
    p##_##n##_isSet_inputVect = nullptr;

/** Define a macro that resets a member variable */
#define UBCC1PI_MACRO_RESET_MEMBER(p, q, r, t, n)                                                                                          \
    q.n.Reset();

/** Define a macro that resets a member vector */
#define UBCC1PI_MACRO_RESET_MEMBER_VECTOR(p, q, r, t, n)                                                                                   \
    p##_##n##_vect->clear();                                                                                                               \
    p##_##n##_isSet_vect->clear();

/** Define a macro that fills the output vectors */
#define UBCC1PI_MACRO_FILL_MEMBER_VECTOR(p, q, r, t, n)                                                                                    \
    p##_##n##_vect->push_back(*q.n.m_pValue);                                                                                              \
    p##_##n##_isSet_vect->push_back(q.n.m_isSet);

/** Define a macro that reads the vectors in the tree */
/*#define UBCC1PI_MACRO_READ_MEMBER_VECTOR(p, q, r, t, n)                                                                                    \
    std::cout << "Ubcc1pi Reading " << #q "_" #p "_" #n<<":";                                                                              \
    for (unsigned int i = 0; i < p##_##n##_inputVect->size(); ++i)                                                                         \
    {                                                                                                                                      \
        if (p##_##n##_isSet_inputVect->at(i))                                                                                              \
        {                                                                                                                                  \
            q.at(i).n.Set(p##_##n##_inputVect->at(i));                                                                                     \
            std::cout <<" (" << i << ")" << p##_##n##_inputVect->at(i);                                                                    \
        }                                                                                                                                  \
    }                                                                                                                                      \
    std::cout << std::endl;*/

#define UBCC1PI_MACRO_READ_MEMBER_VECTOR(p, q, r, t, n)                                                                                    \
    for (unsigned int i = 0; i < p##_##n##_inputVect->size(); ++i)                                                                         \
    {                                                                                                                                      \
        if (p##_##n##_isSet_inputVect->at(i))                                                                                              \
        {                                                                                                                                  \
            q.at(i).n.Set(p##_##n##_inputVect->at(i));                                                                                     \
        }                                                                                                                                  \
    }                                                                                                                                      \

/** Define a macro to count the size of the member vector - here q will receive the size of the vector */
#define UBCC1PI_MACRO_GET_MEMBER_VECTOR_SIZE(p, q, r, t, n)                                                                                \
    *q = p##_##n##_inputVect->size();

/** Define a macro that prints each of the member variables */
#define UBCC1PI_MACRO_PRINT_MEMBER(p, q, r, t, n)                                                                                          \
    std::cout << std::setw(28) << "(" #t ")" << "  ";                                                                                      \
    std::cout << std::setw(44) << (#q "." #n "()") << "  ";                                                                                \
    std::cout << q.n.ToString() << std::endl;

/** Define a macro to bind a member variable to an output branch */
#define UBCC1PI_MACRO_BIND_OUTPUT_BRANCH(p, q, r, t, n)                                                                                    \
    pTree->Branch(#p "_" #n, q.n.m_pValue.get());                                                                                          \
    pTree->Branch(#p "_" #n "_isSet", &q.n.m_isSet);

/** Define a macro to bind a member variable to an output vector branch */
#define UBCC1PI_MACRO_BIND_OUTPUT_VECTOR_BRANCH(p, q, r, t, n)                                                                             \
    pTree->Branch(#p "_" #n "_vect", p##_##n##_vect.get());                                                                                \
    pTree->Branch(#p "_" #n "_isSet_vect", p##_##n##_isSet_vect.get());

/** Define a macro to bind a member variable to an input branch */
#define UBCC1PI_MACRO_BIND_INPUT_BRANCH(p, q, r, t, n)                                                                                     \
    pTree->SetBranchAddress(#p "_" #n "_isSet", &q.n.m_isSet);                                                                             \
    if (r) {pTree->SetBranchAddress(#p "_" #n, &(q.n.m_pAddress));} else {pTree->SetBranchAddress(#p "_" #n, q.n.m_pAddress);}

/** Define a macro to bind a member variable to an input vector branch */
#define UBCC1PI_MACRO_BIND_INPUT_VECTOR_BRANCH(p, q, r, t, n)                                                                              \
    pTree->SetBranchAddress(#p "_" #n "_vect", &p##_##n##_inputVect);                                                                      \
    pTree->SetBranchAddress(#p "_" #n "_isSet_vect", &p##_##n##_isSet_inputVect);





/** Define a macro that declares a member variable and an associated boolean to check if that variable has been set */
#define PELEE_MACRO_DECLARE_MEMBER(p, q, r, t, n)                                                                                        \
    Member<t> n = Member<t>(#n);

/** Define a macro that declares a member vector, used to bind particles to vectors in trees */
#define PELEE_MACRO_DECLARE_MEMBER_VECTOR(p, q, r, t, n)                                                                                   \
    std::shared_ptr< std::vector<t> > p##_##n##_vect;                                                                                      \
    std::vector<t> * p##_##n##_inputVect;

/** Define a macro that initializes the member vectors */
#define PELEE_MACRO_INIT_MEMBER_VECTOR(p, q, r, t, n)                                                                                      \
    p##_##n##_vect = std::make_shared< std::vector<t> >();                                                                                 \
    p##_##n##_inputVect = nullptr;

/** Define a macro that resets a member variable */
#define PELEE_MACRO_RESET_MEMBER(p, q, r, t, n)                                                                                            \
    q.n.Reset();                                                                                                                           \
    q.n.m_isSet = true;

/** Define a macro that resets a member vector */
#define PELEE_MACRO_RESET_MEMBER_VECTOR(p, q, r, t, n)                                                                                     \
    p##_##n##_vect->clear();

/** Define a macro that fills the output vectors */
#define PELEE_MACRO_FILL_MEMBER_VECTOR(p, q, r, t, n)                                                                                      \
    p##_##n##_vect->push_back(*q.n.m_pValue);                                                                                              \
    std::cout<<#q "_" #p "_" #n "_vect: "<<*q.n.m_pValue<<std::endl;

/** Define a macro that reads the vectors in the tree */
/*#define PELEE_MACRO_READ_MEMBER_VECTOR(p, q, r, t, n)                                                                                      \
    std::cout << "Pelee Reading " << #q "_" #p "_" #n<<":";                                                                                \
    for (unsigned int i = 0; i < p##_##n##_inputVect->size(); ++i)                                                                         \
    {                                                                                                                                      \
        q.at(i).n.Set(p##_##n##_inputVect->at(i));                                                                                         \
        std::cout <<" ("<<i<<")"<< p##_##n##_inputVect->at(i);                                                                             \
    }                                                                                                                                      \
    std::cout << std::endl;*/

/** Define a macro that reads the vectors in the tree */
#define PELEE_MACRO_READ_MEMBER_VECTOR(p, q, r, t, n)                                                                                      \
    for (unsigned int i = 0; i < p##_##n##_inputVect->size(); ++i)                                                                         \
    {                                                                                                                                      \
        q.at(i).n.Set(p##_##n##_inputVect->at(i));                                                                                         \
    }                                                                                                                                      \

/** Define a macro to count the size of the member vector - here q will receive the size of the vector */
#define PELEE_MACRO_GET_MEMBER_VECTOR_SIZE(p, q, r, t, n)                                                                                  \
    *q = p##_##n##_inputVect->size();

/** Define a macro that prints each of the member variables */
#define PELEE_MACRO_PRINT_MEMBER(p, q, r, t, n)                                                                                            \
    std::cout << std::setw(28) << "(" #t ")" << "  ";                                                                                      \
    std::cout << std::setw(44) << (#q "." #n "()") << "  ";                                                                                \
    std::cout << q.n.ToString() << std::endl;

/** Define a macro to bind a member variable to an output branch */
#define PELEE_MACRO_BIND_OUTPUT_BRANCH(p, q, r, t, n)                                                                                      \
    q.n.m_isSet = true;                                                                                                                    \
    pTree->Branch(#n, q.n.m_pValue.get());

/** Define a macro to bind a member variable to an output vector branch */
#define PELEE_MACRO_BIND_OUTPUT_VECTOR_BRANCH(p, q, r, t, n)                                                                               \
    pTree->Branch(#n, p##_##n##_vect.get());

/** Define a macro to bind a member variable to an input branch */
#define PELEE_MACRO_BIND_INPUT_BRANCH(p, q, r, t, n)                                                                                       \
    q.n.m_isSet = true;                                                                                                                    \
    if (r) {pTree->SetBranchAddress(#n, &(q.n.m_pAddress));} else {pTree->SetBranchAddress(#n, q.n.m_pAddress);}

/** Define a macro to bind a member variable to an input vector branch */
#define PELEE_MACRO_BIND_INPUT_VECTOR_BRANCH(p, q, r, t, n)                                                                                \
    pTree->SetBranchAddress(#n, &p##_##n##_inputVect);



/** Define a macro to convert between members of different classes */
#define PELEE_TO_UBCC1PI_MEMBER_CONVERSION(p, q, r)                                                                                        \
    q.Set(p);

} // namespace ubcc1pi

#endif
