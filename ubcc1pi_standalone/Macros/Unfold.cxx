/**
 *  @file  ubcc1pi_standalone/Macros/Unfold.cxx
 *
 *  @brief The implementation file of the Unfold macro
 */

#include "ubcc1pi_standalone/Macros/Macros.h"

#include <stdexcept>
#include <fstream>

using namespace ubcc1pi;

namespace ubcc1pi_macros
{

void Unfold(const Config &config)
{
    // ifstream f("/uboone/app/users/zarko/getDataInfo.py");
    // const bool zarkosToolExists = f.good(); 
    const std::string path = "$MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/unfolding-scripts/";
    const std::vector<std::string> commands {
        "./univmake " + path + "file_list.txt " + path + "ubcc1pi_bin_config.txt /uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output.root " + path + "file_properties.txt" // [FILE_PROPERTIES_CONFIG_FILE]
    };
    // ./univmake $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/unfolding-scripts/file_list.txt $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/unfolding-scripts/ubcc1pi_bin_config.txt /uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output.root $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/unfolding-scripts/file_properties.txt

    // make clean && make && rm -f /uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output.root && ./univmake $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/xsec_analyzer/nuwro_file_properties.txt $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/xsec_analyzer/ubcc1pi_bin_config.txt /uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output.root $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/xsec_analyzer/nuwro_file_properties.txt
    // make clean && make && rm -f /uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output.root && ./univmake $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/xsec_analyzer/nuwro_file_properties_2Percent.txt $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/xsec_analyzer/ubcc1pi_bin_config.txt /uboone/data/users/jdetje/ubcc1pi_univmake/univmake_output.root $MRB_SOURCE/ubcc1pi/ubcc1pi_standalone/xsec_analyzer/nuwro_file_properties_2Percent.txt

    for (const auto command : commands)
    {
        const auto status = std::system(command.c_str());
        if (status != 0) throw std::runtime_error("Failed to execute command: " + command);
    }
}

} // namespace ubcc1pi_macros
