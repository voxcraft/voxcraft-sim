#if !defined(CTOOL_H)
#define CTOOL_H



#include <unistd.h>
#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <boost/process.hpp>
#include <boost/process/start_dir.hpp>
namespace bp = boost::process;
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;
#include <boost/foreach.hpp>

namespace ctool {


    inline std::string u_format_now(std::string format) {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        std::stringstream folderName;
        folderName << std::put_time(std::localtime(&in_time_t), format.c_str());
        return folderName.str();
    }

    inline bool u_with_ext(fs::path file, std::string ext) {

        std::string ext_file = file.filename().extension().string();
        boost::to_upper(ext);
        boost::to_upper(ext_file);

        return ext==ext_file;
    }


    std::string get_exe_path() {
        char exe_path[256];
        size_t end = readlink("/proc/self/exe", &exe_path[0], 255);
        exe_path[end] = '\0';
        std::string ret(exe_path);
        return ret;
    }

    void ptree_merge(const pt::ptree& source, pt::ptree& destination) {
        pt::ptree vxd = source.get_child("VXD", source);
        BOOST_FOREACH( pt::ptree::value_type const&v, vxd.get_child("") ) {
            std::string replace = v.second.get<std::string>("<xmlattr>.replace", "");
            if (replace.length()>0) {
                destination.put_child( replace , v.second);
            }
        }
    }
}


#endif // CTOOL_H