#if !defined(VX3_SIMULATION_MANAGER)
#define VX3_SIMULATION_MANAGER
#include <boost/filesystem.hpp>
#include <iostream>
#include <thread>
#include <utility>
#include <vector>
namespace fs = boost::filesystem;
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include "VX3_VoxelyzeKernel.cuh"

class VX3_SimulationManager {
  public:
    VX3_SimulationManager(std::vector<std::vector<fs::path>> in_sub_batches, fs::path in_base, fs::path in_input_dir,
                          int in_num_of_devices);
    ~VX3_SimulationManager();

    void start();
    void readVXD(fs::path base, std::vector<fs::path> files, int device_index);
    std::vector<std::vector<fs::path>> splitIntoSubBatches();
    void startKernel(int num_tasks, int device_index);
    void collectResults(int num_simulation, int device_index);
    void sortResults();
    void enlargeGPUHeapSize();
    void ParseMathTree(VX3_MathTreeToken *field_ptr, size_t max_length, std::string node_address, pt::ptree &tree);

    /* DATA */
    int num_of_devices;                              // Total number of GPUs on one single node. One
                                                     // DeepGreen node has 8 GPUs.
    std::vector<VX3_VoxelyzeKernel *> d_voxelyze_3s; // Multiple device memory passing to different device.

    std::vector<std::vector<fs::path>> sub_batches;
    fs::path input_dir;
    fs::path base;
    VX3_VoxelyzeKernel h_d_base;

    // fs::path input_directory;
    // fs::path output_file;

    std::vector<VX3_SimulationResult> h_results;

    double HeapSize=1;
};

#endif // VX3_SIMULATION_MANAGER
