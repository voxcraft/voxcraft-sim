#include "VX3_SimulationManager.cuh"
#include "ctool.h"
#include <boost/algorithm/string/case_conv.hpp>
#include <queue>
#include <stack>
#include <utility>

#include "VX3_VoxelyzeKernel.cuh"
#include "VX_Sim.h" //readVXA

__global__ void CUDA_Simulation(VX3_VoxelyzeKernel *d_voxelyze_3, int num_simulation, int device_index) {
    int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_index < num_simulation) {
        VX3_VoxelyzeKernel *d_v3 = &d_voxelyze_3[thread_index];
        if (d_v3->num_d_links == 0 and d_v3->num_d_voxels == 0) {
            printf(COLORCODE_BOLD_RED "No links and no voxels. Simulation %d (%s) abort.\n" COLORCODE_RESET, thread_index,
                   d_v3->vxa_filename);
            return;
        }
        d_v3->syncVectors();           // Everytime we pass a class with VX3_vectors in
                                       // it, we should sync hd_vector to d_vector first.
        d_v3->saveInitialPosition();
        d_v3->isSurfaceChanged = true; // trigger surface regenerating and calculate normal thrust for the first time
        d_v3->registerTargets();
        printf(COLORCODE_GREEN "%d) Simulation %d runs: %s.\n" COLORCODE_RESET, device_index, thread_index, d_v3->vxa_filename);
        // printf("%d) Simulation %d: links %d, voxels %d.\n", device_index, i,
        // d_v3->num_d_links, d_v3->num_d_voxels); printf("%d) Simulation %d
        // enableAttach %d.\n", device_index, i, d_v3->enableAttach);
        //
        // print check regenerateSurfaceVoxels() is correct. (TODO: shouldn't
        // this be tested in seperate test code? :) printf("all voxels:"); for
        // (int j=0;j<d_v3->num_d_voxels;j++) {
        //     printf(" [%d]%p ", j, &d_v3->d_voxels[j]);
        // }
        // printf("\nsurface:");
        // for (int j=0;j<d_v3->num_d_surface_voxels;j++) {
        //     printf(" [%d]%p ", j, d_v3->d_surface_voxels[j]);
        // }
        //
        if (d_v3->RecordStepSize) { // output History file
            // rescale the whole space. so history file can contain less digits. ( e.g. not 0.000221, but 2.21 )
            printf("\n{{{setting}}}<rescale>0.001</rescale>\n");
            // materials' color
            for (int i = 0; i < d_v3->num_d_voxelMats; i++) {
                auto &mat = d_v3->d_voxelMats[i];
                printf("{{{setting}}}<matcolor><id>%d</id><r>%.2f</r><g>%.2f</g><b>%.2f</b><a>%.2f</a></matcolor>\n", mat.matid,
                       mat.r / 255., mat.g / 255., mat.b / 255., mat.a / 255.);
            }
            printf("\n{{{setting}}}<voxel_size>%f</voxel_size>\n", d_v3->voxSize);
        }

        double vs = 1 / 0.001;

        d_v3->updateCurrentCenterOfMass();
        d_v3->InitializeCenterOfMass();
        int real_stepsize = int(d_v3->RecordStepSize / (10000 * d_v3->recommendedTimeStep() * d_v3->DtFrac))+1;
        printf("real_stepsize: %d ; recommendedTimeStep %f; d_v3->DtFrac %f . \n", real_stepsize, d_v3->recommendedTimeStep(),
               d_v3->DtFrac);
        // printf("Initial CoM: %f %f %f mm\n",
        // d_v3->initialCenterOfMass.x*1000, d_v3->initialCenterOfMass.y*1000,
        // d_v3->initialCenterOfMass.z*1000);
        for (int j = 0; j < 1000000; j++) { // Maximum Steps 1000000
            if (d_v3->StopConditionMet())
                break;
            if (!d_v3->doTimeStep()) {
                printf(COLORCODE_BOLD_RED "\n%d) Simulation %d Diverged: %s.\n" COLORCODE_RESET, device_index, thread_index,
                       d_v3->vxa_filename);
                break;
            }
            if (d_v3->RecordStepSize) { // output History file
                if (j % real_stepsize == 0) {
                    if (d_v3->RecordVoxel) {
                        // Voxels
                        printf("<<<Step%d Time:%f>>>", j, d_v3->currentTime);
                        for (int i = 0; i < d_v3->num_d_surface_voxels; i++) {
                            auto v = d_v3->d_surface_voxels[i];
                            if (v->removed)
                                continue;
                            if (v->isSurface()) {
                                printf("%.1f,%.1f,%.1f,", v->pos.x * vs, v->pos.y * vs, v->pos.z * vs);
                                printf("%.1f,%.2f,%.2f,%.2f,", v->orient.AngleDegrees(), v->orient.x, v->orient.y, v->orient.z);
                                VX3_Vec3D<double> ppp, nnn;
                                nnn = v->cornerOffset(NNN);
                                ppp = v->cornerOffset(PPP);
                                printf("%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,", nnn.x * vs, nnn.y * vs, nnn.z * vs, ppp.x * vs, ppp.y * vs,
                                       ppp.z * vs);
                                printf("%d,", v->mat->matid); // for coloring
                                printf("%.1f,", v->localSignal);  // for coloring as well.
                                printf(";");
                            }
                        }
                        printf("<<<>>>");
                    }
                    if (d_v3->RecordLink) {
                        // Links
                        printf("|[[[%d]]]", j);
                        for (int i = 0; i < d_v3->d_v_links.size(); i++) {
                            auto l = d_v3->d_v_links[i];
                            if (l->removed)
                                continue;
                            // only draw links that are not detached.
                            if (!l->isDetached) {
                                auto v1 = l->pVPos;
                                printf("%.4f,%.4f,%.4f,", v1->pos.x, v1->pos.y, v1->pos.z);
                                auto v2 = l->pVNeg;
                                printf("%.4f,%.4f,%.4f,", v2->pos.x, v2->pos.y, v2->pos.z);
                                printf(";");
                            }
                        }
                        printf("[[[]]]");
                    }
                    printf("\n");
                }
            }
        }
        d_v3->updateCurrentCenterOfMass();
        d_v3->computeFitness();
        printf(COLORCODE_BLUE "%d) Simulation %d ends: %s Time: %f, angleSampleTimes: %d.\n" COLORCODE_RESET, device_index, thread_index,
               d_v3->vxa_filename, d_v3->currentTime, d_v3->angleSampleTimes);
    }
}

VX3_SimulationManager::VX3_SimulationManager(std::vector<std::vector<fs::path>> in_sub_batches, fs::path in_base, fs::path in_input_dir,
                                             int in_num_of_devices)
    : sub_batches(in_sub_batches), base(in_base), num_of_devices(in_num_of_devices), input_dir(in_input_dir) {
    d_voxelyze_3s.resize(num_of_devices);
    for (int i = 0; i < num_of_devices; i++) {
        d_voxelyze_3s[i] = NULL;
    }
}
VX3_SimulationManager::~VX3_SimulationManager() {
    for (auto d : d_voxelyze_3s) {
        if (d)
            VcudaFree(d);
    }
}

void VX3_SimulationManager::start() {
    for (int device_index = 0; device_index < num_of_devices; device_index++) { // multi GPUs
        auto files = sub_batches[device_index];
        if (files.size()) {
            VcudaSetDevice(device_index);
            printf("=== set device to %d for %ld simulations ===\n", device_index, files.size());
            // readVXA(base)
            readVXD(base, files, device_index);
            startKernel(files.size(), device_index);
        }
    }
    VcudaDeviceSynchronize();
    for (int device_index = 0; device_index < num_of_devices; device_index++) { // multi GPUs
        auto files = sub_batches[device_index];
        collectResults(files.size(), device_index);
    }
    sortResults();
}

void VX3_SimulationManager::ParseMathTree(VX3_MathTreeToken *field_ptr, size_t max_length, std::string node_address, pt::ptree &tree) {
    // Classic BFS, push all token into stack
    std::queue<pt::ptree> frontier;
    std::stack<std::pair<std::string, std::string>> tokens;
    tokens.push(make_pair((std::string) "mtEND", (std::string) ""));
    auto root = tree.get_child_optional(node_address);
    if (!root) {
        // printf(COLORCODE_BOLD_RED "ERROR: No ParseMathTree %s in VXA.\n", node_address.c_str());
        return;
    }
    frontier.push(tree.get_child(node_address));
    while (!frontier.empty()) {
        std::queue<pt::ptree> next_frontier;
        auto t = frontier.front();
        frontier.pop();
        BOOST_FOREACH (pt::ptree::value_type &v_child, t.get_child("")) {
            std::string value = v_child.second.data();
            boost::trim_right(value);
            std::string op = v_child.first.data();
            boost::trim_right(op);

            // std::cout << op << ":" << value << "\n";
            tokens.push(make_pair(op, value));
            frontier.push(v_child.second);
        }
    }
    // pop from stack to VX3_MathTreeToken* (so we get a reversed order)
    int i = 0;
    while (!tokens.empty()) {
        if (i > max_length) {
            printf(COLORCODE_BOLD_RED "ERROR: Token size overflow.\n");
            return;
        }
        std::pair<std::string, std::string> tok = tokens.top();
        VX3_MathTreeToken *p = &field_ptr[i];
        if (tok.first == "mtEND") {
            p->op = mtEND;
        } else if (tok.first == "mtVAR") {
            p->op = mtVAR;
            if (tok.second == "x") {
                p->value = 0;
            } else if (tok.second == "y") {
                p->value = 1;
            } else if (tok.second == "z") {
                p->value = 2;
            } else if (tok.second == "hit") {
                p->value = 3;
            } else if (tok.second == "t") {
                p->value = 4;
            } else if (tok.second == "angle") {
                p->value = 5;
            } else if (tok.second == "targetCloseness") {
                p->value = 6;
            } else if (tok.second == "numClosePairs") {
                p->value = 7;
            } else if (tok.second == "num_voxel") {
                p->value = 8;
            } else {
                printf(COLORCODE_BOLD_RED "ERROR: No such variable.\n");
                break;
            }
        } else if (tok.first == "mtCONST") {
            p->op = mtCONST;
            try {
                p->value = std::stod(tok.second);
            } catch(...) {
                printf(COLORCODE_BOLD_RED "ERROR: mtCONST with no number.\n");
                break;
            }
        } else if (tok.first == "mtADD") {
            p->op = mtADD;
        } else if (tok.first == "mtSUB") {
            p->op = mtSUB;
        } else if (tok.first == "mtMUL") {
            p->op = mtMUL;
        } else if (tok.first == "mtDIV") {
            p->op = mtDIV;
        } else if (tok.first == "mtPOW") {
            p->op = mtPOW;
        } else if (tok.first == "mtSQRT") {
            p->op = mtSQRT;
        } else if (tok.first == "mtE") {
            p->op = mtE;
        } else if (tok.first == "mtPI") {
            p->op = mtPI;
        } else if (tok.first == "mtSIN") {
            p->op = mtSIN;
        } else if (tok.first == "mtCOS") {
            p->op = mtCOS;
        } else if (tok.first == "mtTAN") {
            p->op = mtTAN;
        } else if (tok.first == "mtATAN") {
            p->op = mtATAN;
        } else if (tok.first == "mtLOG") {
            p->op = mtLOG;
        } else if (tok.first == "mtINT") {
            p->op = mtINT;
        } else if (tok.first == "mtABS") {
            p->op = mtABS;
        } else if (tok.first == "mtNOT") {
            p->op = mtNOT;
        } else if (tok.first == "mtGREATERTHAN") {
            p->op = mtGREATERTHAN;
        } else if (tok.first == "mtLESSTHAN") {
            p->op = mtLESSTHAN;
        } else if (tok.first == "mtAND") {
            p->op = mtAND;
        } else if (tok.first == "mtOR") {
            p->op = mtOR;
        } else if (tok.first == "mtNORMALCDF") {
            p->op = mtNORMALCDF;
        } else {
            printf(COLORCODE_BOLD_RED "ERROR: Token Operation not implemented.\n");
            break;
        }
        i++;
        tokens.pop();
    }
}

void VX3_SimulationManager::readVXD(fs::path base, std::vector<fs::path> files, int device_index) {
    pt::ptree pt_baseVXA;
    pt::read_xml(base.string(), pt_baseVXA);

    int num_simulation = files.size();

    VcudaMalloc((void **)&d_voxelyze_3s[device_index], num_simulation * sizeof(VX3_VoxelyzeKernel));

    int i = 0;
    for (auto &file : files) {
        // Read VXD file, clone base VXA, replace parts specified in VXD, send
        // to MainSim.ReadVXA to process. printf("reading %s\n",
        // (input_dir/file).c_str());
        pt::ptree pt_VXD;
        pt::read_xml((input_dir / file).string(), pt_VXD);
        pt::ptree pt_merged = pt_baseVXA;
        ctool::ptree_merge(pt_VXD, pt_merged);
        std::ostringstream stream_merged;
        std::string str_merged;
        pt::write_xml(stream_merged, pt_merged);
        str_merged = stream_merged.str();
        CXML_Rip XML;
        XML.fromXMLText(&str_merged);
        CVX_Environment MainEnv;
        CVX_Sim MainSim;
        CVX_Object MainObj;
        MainEnv.pObj = &MainObj; // connect environment to object
        MainSim.pEnv = &MainEnv; // connect Simulation to envirnment
        std::string RetMessage;
        // std::cout<<str_merged;
        MainSim.ReadVXA(&XML, &RetMessage);
        MainSim.Import(NULL, NULL, &RetMessage);
        if (!RetMessage.empty()) {
            printf(COLORCODE_BOLD_RED "%s\n" COLORCODE_RESET, RetMessage.c_str());
        }
        // for (auto m:MainSim.Vx.voxelMats) {
        //     int i=0;
        //     for (auto mm:m->dependentMaterials) {
        //         printf("m:%p %d/%ld -> mm: %p\n", m, i,
        //         m->dependentMaterials.size(), mm); i++;
        //     }
        // }
        VX3_VoxelyzeKernel h_d_tmp(&MainSim);
        // More VXA settings which is new in VX3
        strcpy(h_d_tmp.vxa_filename, file.filename().c_str());

        std::string RawPrint = pt_merged.get<std::string>("VXA.RawPrint", "");
        std::cout << RawPrint << "\n";

        ParseMathTree(h_d_tmp.StopConditionFormula, sizeof(h_d_tmp.StopConditionFormula),
                      "VXA.Simulator.StopCondition.StopConditionFormula", pt_merged);
        h_d_tmp.EnableCollision = pt_merged.get<bool>("VXA.Simulator.AttachDetach.EnableCollision", true);
        h_d_tmp.enableAttach = pt_merged.get<bool>("VXA.Simulator.AttachDetach.EnableAttach", false);
        h_d_tmp.enableDetach = pt_merged.get<bool>("VXA.Simulator.AttachDetach.EnableDetach", false);
        h_d_tmp.watchDistance = pt_merged.get<double>("VXA.Simulator.AttachDetach.watchDistance", 1.0);
        h_d_tmp.boundingRadius = pt_merged.get<double>("VXA.Simulator.AttachDetach.boundingRadius", 0.75);
        h_d_tmp.SafetyGuard = pt_merged.get<int>("VXA.Simulator.AttachDetach.SafetyGuard", 500);
        ParseMathTree(h_d_tmp.AttachCondition[0], sizeof(h_d_tmp.AttachCondition[0]),
                      "VXA.Simulator.AttachDetach.AttachCondition.Condition_0", pt_merged);
        ParseMathTree(h_d_tmp.AttachCondition[1], sizeof(h_d_tmp.AttachCondition[1]),
                      "VXA.Simulator.AttachDetach.AttachCondition.Condition_1", pt_merged);
        ParseMathTree(h_d_tmp.AttachCondition[2], sizeof(h_d_tmp.AttachCondition[2]),
                      "VXA.Simulator.AttachDetach.AttachCondition.Condition_2", pt_merged);
        ParseMathTree(h_d_tmp.AttachCondition[3], sizeof(h_d_tmp.AttachCondition[3]),
                      "VXA.Simulator.AttachDetach.AttachCondition.Condition_3", pt_merged);
        ParseMathTree(h_d_tmp.AttachCondition[4], sizeof(h_d_tmp.AttachCondition[4]),
                      "VXA.Simulator.AttachDetach.AttachCondition.Condition_4", pt_merged);
        h_d_tmp.RecordStepSize = pt_merged.get<int>("VXA.Simulator.RecordHistory.RecordStepSize", 0);
        h_d_tmp.RecordLink = pt_merged.get<int>("VXA.Simulator.RecordHistory.RecordLink", 0);
        h_d_tmp.RecordVoxel = pt_merged.get<int>("VXA.Simulator.RecordHistory.RecordVoxel", 1);

        ParseMathTree(h_d_tmp.fitness_function, sizeof(h_d_tmp.fitness_function), "VXA.Simulator.FitnessFunction", pt_merged);
        ParseMathTree(h_d_tmp.force_field.token_x_forcefield, sizeof(h_d_tmp.force_field.token_x_forcefield),
                      "VXA.Simulator.ForceField.x_forcefield", pt_merged);
        ParseMathTree(h_d_tmp.force_field.token_y_forcefield, sizeof(h_d_tmp.force_field.token_y_forcefield),
                      "VXA.Simulator.ForceField.y_forcefield", pt_merged);
        ParseMathTree(h_d_tmp.force_field.token_z_forcefield, sizeof(h_d_tmp.force_field.token_z_forcefield),
                      "VXA.Simulator.ForceField.z_forcefield", pt_merged);

        // h_d_tmp.EnableTargetCloseness = pt_merged.get<int>("VXA.Simulator.EnableTargetCloseness", 0); abandoned.
        h_d_tmp.SavePositionOfAllVoxels = pt_merged.get<int>("VXA.Simulator.SavePositionOfAllVoxels", 0);
        h_d_tmp.MaxDistInVoxelLengthsToCountAsPair = pt_merged.get<double>("VXA.Simulator.MaxDistInVoxelLengthsToCountAsPair", 0);

        h_d_tmp.EnableCilia = pt_merged.get<int>("VXA.Simulator.EnableCilia", 0);
        h_d_tmp.EnableSignals = pt_merged.get<int>("VXA.Simulator.EnableSignals", 0);
        
        // for Secondary Experiment
        h_d_tmp.SecondaryExperiment = pt_merged.get<int>("VXA.Simulator.SecondaryExperiment", 0);
        h_d_tmp.ReinitializeInitialPositionAfterThisManySeconds = pt_merged.get<double>("VXA.Simulator.ReinitializeInitialPositionAfterThisManySeconds", 0.0);

        h_d_tmp.EnableExpansion = pt_merged.get<int>("VXA.Simulator.EnableExpansion", 0);

        HeapSize = pt_merged.get<double>("VXA.GPU.HeapSize", 0.5);
        if (HeapSize > 1.0) {
            HeapSize = 0.99;
        }
        if (HeapSize < 0.01) {
            HeapSize = 0.01;
        }

        VcudaMemcpy(d_voxelyze_3s[device_index] + i, &h_d_tmp, sizeof(VX3_VoxelyzeKernel), cudaMemcpyHostToDevice);
        i++;
    }
}

// GPU Heap is for in-kernel malloc(). Refer to
// https://stackoverflow.com/a/34795830/7001199
void VX3_SimulationManager::enlargeGPUHeapSize() {
    size_t HeapSizeInBytes;
    size_t free, total;
    VcudaMemGetInfo(&free, &total);
    printf("Total GPU memory %ld bytes.\n", total);
    HeapSizeInBytes = HeapSize * total; // add some additional size
    printf("Set GPU heap size to be %ld bytes.\n", HeapSizeInBytes);
    VcudaDeviceSetLimit(cudaLimitMallocHeapSize,
                        HeapSizeInBytes); // Set Heap Memory to 1G, instead of merely 8M.

    // if "Lane User Stack Overflow" ocurs, maybe Stack Size too small, can try this:
    // VcudaDeviceSetLimit(cudaLimitStackSize, 2048);
}

void VX3_SimulationManager::startKernel(int num_simulation, int device_index) {
    int threadsPerBlock = 512;
    int numBlocks = (num_simulation + threadsPerBlock - 1) / threadsPerBlock;
    if (numBlocks == 1)
        threadsPerBlock = num_simulation;
    // printf("Starting kernel on device %d. passing d_voxelyze_3s[device_index]
    // %p.\n", device_index, d_voxelyze_3s[device_index]);
    // VX3_VoxelyzeKernel *result_voxelyze_kernel = (VX3_VoxelyzeKernel
    // *)malloc(
    //     num_simulation * sizeof(VX3_VoxelyzeKernel));
    // VcudaMemcpy(result_voxelyze_kernel, d_voxelyze_3s[device_index],
    //             num_simulation * sizeof(VX3_VoxelyzeKernel),
    //             cudaMemcpyDeviceToHost);
    enlargeGPUHeapSize();
    CUDA_Simulation<<<numBlocks, threadsPerBlock>>>(d_voxelyze_3s[device_index], num_simulation, device_index);
    CUDA_CHECK_AFTER_CALL();
}

void VX3_SimulationManager::collectResults(int num_simulation, int device_index) {
    // insert results to h_results
    VX3_VoxelyzeKernel *result_voxelyze_kernel = (VX3_VoxelyzeKernel *)malloc(num_simulation * sizeof(VX3_VoxelyzeKernel));
    VcudaMemcpy(result_voxelyze_kernel, d_voxelyze_3s[device_index], num_simulation * sizeof(VX3_VoxelyzeKernel), cudaMemcpyDeviceToHost);
    for (int i = 0; i < num_simulation; i++) {
        VX3_SimulationResult tmp;
        tmp.currentTime = result_voxelyze_kernel[i].currentTime;
        tmp.fitness_score = result_voxelyze_kernel[i].fitness_score;
        tmp.x = result_voxelyze_kernel[i].currentCenterOfMass.x;
        tmp.y = result_voxelyze_kernel[i].currentCenterOfMass.y;
        tmp.z = result_voxelyze_kernel[i].currentCenterOfMass.z;
        result_voxelyze_kernel[i].initialCenterOfMass.copyTo(tmp.initialCenterOfMass);
        result_voxelyze_kernel[i].currentCenterOfMass.copyTo(tmp.currentCenterOfMass);

        tmp.numClosePairs = result_voxelyze_kernel[i].numClosePairs;
        tmp.voxSize = result_voxelyze_kernel[i].voxSize;
        tmp.num_voxel = result_voxelyze_kernel[i].num_d_voxels;
        tmp.vxa_filename = result_voxelyze_kernel[i].vxa_filename;
        VX3_Voxel *tmp_v;
        tmp_v = (VX3_Voxel *)malloc(result_voxelyze_kernel[i].num_d_voxels * sizeof(VX3_Voxel));
        VcudaMemcpy(tmp_v, result_voxelyze_kernel[i].d_voxels, result_voxelyze_kernel[i].num_d_voxels * sizeof(VX3_Voxel),
                   cudaMemcpyDeviceToHost);
        tmp.SavePositionOfAllVoxels = result_voxelyze_kernel[i].SavePositionOfAllVoxels;
        VX3_Vec3D<>* tmp_init;
        tmp_init = (VX3_Vec3D<>*)malloc(result_voxelyze_kernel[i].num_d_voxels * sizeof(VX3_Vec3D<>));
        VcudaMemcpy(tmp_init, result_voxelyze_kernel[i].d_initialPosition, result_voxelyze_kernel[i].num_d_voxels * sizeof(VX3_Vec3D<>), cudaMemcpyDeviceToHost);
        tmp.num_measured_voxel = 0;
        tmp.total_distance_of_all_voxels = 0.0;
        for (int j = 0; j < result_voxelyze_kernel[i].num_d_voxels; j++) {
            tmp.voxel_init_pos.push_back(Vec3D<>(tmp_init[j].x, tmp_init[j].y, tmp_init[j].z));
            tmp.voxel_position.push_back(Vec3D<>(tmp_v[j].pos.x, tmp_v[j].pos.y, tmp_v[j].pos.z));
            tmp.voxel_mats.push_back(tmp_v[j].matid);
            if (tmp_v[j].isMeasured) {
                tmp.num_measured_voxel ++;
                tmp.total_distance_of_all_voxels += tmp.voxel_position.back().Dist(tmp.voxel_init_pos.back());
            }
        }
        delete tmp_v;

        // tmp.computeFitness();
        h_results.push_back(tmp);
    }
}

void VX3_SimulationManager::sortResults() { sort(h_results.begin(), h_results.end(), VX3_SimulationResult::compareFitnessScore); }
