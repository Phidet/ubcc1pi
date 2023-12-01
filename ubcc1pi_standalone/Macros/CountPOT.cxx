/**
 *  @file  ubcc1pi_standalone/Macros/CountPOT.cxx
 *
 *  @brief The implementation file of the CountPOT macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include "ubcc1pi_standalone/Objects/FileReader.h"

#include "ubcc1pi_standalone/Helpers/AnalysisHelper.h"

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void CountPOT(const Config &config)
{
    for (const auto &[run, normalisation, sampleType, useThisFile, filePath] : config.input.files)
    {
        if(!useThisFile) continue;
        if(sampleType==AnalysisHelper::DataBNB || sampleType==AnalysisHelper::DataEXT) continue;

        FileReader<EventPeLEE, SubrunPeLEE> readerPeLEE(filePath, true);
        const auto nEvents = readerPeLEE.GetNumberOfEvents();
        const auto nSubruns = readerPeLEE.GetNumberOfSubruns();
        const auto pEventPeLEE = readerPeLEE.GetBoundEventAddress();
        const auto pSubrun = readerPeLEE.GetBoundSubrunAddress();

        // Count the total POT
        std::cout << "Processing file: " << filePath << std::endl;
        std::cout << "  - Getting total POT over " << nSubruns << " sub-runs." << std::endl;

        float totalPOT = 0.f;
        for (unsigned int i = 0; i < nSubruns; ++i)
        {
            readerPeLEE.LoadSubrun(i);

            if (pSubrun->pot.IsSet())
                totalPOT += pSubrun->pot();
        }

        std::cout << "  - POT = " << totalPOT << std::endl;

        // // Count the total number of events
        // std::cout << "  - Getting event counts over " << nEvents << " events." << std::endl;
        // float nEventsPassingCCInc = 0.f;
        // float nEventsTrueCCInc = 0.f;
        // float nEventsWeighted = 0.f;
        // float nEventsPassingCCIncWeighted = 0.f;
        // float nEventsTrueCCIncWeighted = 0.f;

        // for (unsigned int i = 0; i < nEvents; ++i)
        // {
        //     readerPeLEE.LoadEvent(i);

        //     const auto weight = AnalysisHelper::GetNominalEventWeight(pEvent);
        //     const auto passesCCInc = pEvent->reco.passesCCInclusive.IsSet() && pEvent->reco.passesCCInclusive();
        //     const auto trueCCInc = AnalysisHelper::IsTrueCCInclusive(pEvent, true);

        //     nEventsPassingCCInc += (passesCCInc ? 1 : 0);
        //     nEventsTrueCCInc += (trueCCInc ? 1 : 0);
        //     nEventsWeighted += weight;
        //     nEventsPassingCCIncWeighted += (passesCCInc ? weight : 0);
        //     nEventsTrueCCIncWeighted += (trueCCInc ? weight : 0);
        // }

        // std::cout << "  - nEvents (unweighted) = " << nEvents << std::endl;
        // std::cout << "  - nEvents passing CCInc. (unweighted) = " << nEventsPassingCCInc << std::endl;
        // std::cout << "  - nEvents true CCInc. (unweighted) = " << nEventsTrueCCInc << std::endl;
        // std::cout << "  - nEvents (nominal weight) = " << nEventsWeighted << std::endl;
        // std::cout << "  - nEvents passing CCInc. (nominal weight) = " << nEventsPassingCCIncWeighted << std::endl;
        // std::cout << "  - nEvents true CCInc. (nominal weight) = " << nEventsTrueCCIncWeighted << std::endl;
    }
}

} // namespace ubcc1pi_macros
