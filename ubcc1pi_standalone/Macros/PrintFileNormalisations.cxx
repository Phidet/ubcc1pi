/**
 *  @file  ubcc1pi_standalone/Macros/PrintFileNormalisations.cxx
 *
 *  @brief The implementation file of the PrintFileNormalisations macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void PrintFileNormalisations(const Config &config)
{
    for (const auto &[fileRun, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        std::cout << std::setprecision(std::numeric_limits<float>::max_digits10) << normalisation << " : " << filePath << std::endl;
    }

    std::cout << "\nFinished" << std::endl;
}

} // namespace ubcc1pi_plots
