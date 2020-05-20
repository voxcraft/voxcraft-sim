#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <thread>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <boost/asio/ip/host_name.hpp>
#include <boost/foreach.hpp> 
// #include <boost/process.hpp> //no! boost::process cannot use in cu file, use std::system instead.

#include "VX3.cuh"

int main() {
    int nDevices;
    VcudaGetDeviceCount(&nDevices);
    if (nDevices<=0) {
        printf("This daemon should run on a node(a Server with GPU).\n");
        return 1;
    }

    //Setup a workspace folder
    fs::path hostname(boost::asio::ip::host_name());
    fs::path done("done");

    try {
        boost::filesystem::create_directory(hostname);
        boost::filesystem::create_directory(hostname/done);
    } catch (...) {}

    fs::path vxh = hostname/"heartbeats.vxh";
    std::string str_worker = "vx3_node_worker";
    fs::path worker(str_worker);

    FILE *fp;
    time_t now;
    while(1) {
        //Monitor
        BOOST_FOREACH(auto const &file, fs::directory_iterator(hostname))   
        { 
            if(fs::is_regular_file(file) && boost::algorithm::to_lower_copy(file.path().extension().string())==".vxt")
            {
                // do something with file
                printf("executing: %s\n", file.path().filename().c_str());
                std::string command = worker.string() + " -i " + file.path().string() + " -o x";

                std::cout << command << "\n";
                int ret = std::system(command.c_str());
                if (ret!=0) {
                    printf(COLORCODE_BOLD_RED "ERROR: executable returns an error.\n");
                }
                
                fs::rename(file.path(), hostname/done/file.path().filename());
            } 
        }

        //Heart beat
        time(&now);
        fp = fopen(vxh.c_str(), "w");
        fprintf(fp, "%ld", now);
        fclose(fp);
        std::this_thread::sleep_for(std::chrono::seconds(2));
    }
}
