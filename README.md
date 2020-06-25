# Voxcraft

This repo is one of the three parts of voxcraft software.

1. [voxcraft-sim](https://github.com/voxcraft/voxcraft-sim): A highly parallelized physics engine that can simulate the voxel-based soft robots. This part utilizes CUDA and GPU.

2. voxcraft-evo: The evolutionary algorithms that can automatically design voxel-based soft robots.

3. [voxcraft-viz](https://github.com/voxcraft/voxcraft-viz): The visualization tool that can playback the history file produced by voxcraft-sim, and it can used for manually design bots and run it in a CPU-based physics engine.

Learn more about the whole voxcraft project (not just software) to get a bigger picture, please refer to: https://voxcraft.github.io/

# Dependencies
* Thrust: Install as `git clone --recurse-submodules https://github.com/thrust/thrust`
    * compile with cmake -DThrust_DIR=/path/to/thrust/thrust/cmake 
    
	

# Installation
* Modify the instances of "compute_70" and "sm_70" in the CMakeLists.txt file to your compute capability represented as an integer. For example if you had a 6.1 GPU you would have a compute capability of 61. To figure out your compute capability see Wikipedia's CUDA  [article](https://en.wikipedia.org/wiki/CUDA).


## On DeepGreen

DeepGreen is UVM’s GPU cluster. We have already compiled everything on DeepGreen, so it will be quite easy to use `voxcraft-sim` on DeepGreen.

Follow a five-minute instruction here: https://github.com/liusida/gpuVoxels-dg-installation (Sorry about the old name there, since we recently renamed our project.)

## On Google Colab

[Google Colab](https://colab.research.google.com/) provides a free online GPU environment.

Create a new notebook and go to Menu->Runtime->Change runtime type, select GPU.

Then, run the script:
```python
!git clone https://github.com/voxcraft/voxcraft-sim.git; cd voxcraft-sim/;

print("Source code downloaded.")

!cd voxcraft-sim; rm build -rf; mkdir build; cd build; cmake -DCMAKE_BUILD_TYPE=Release -DCUDA_DEBUG=OFF ..; make -j 10;

print("Executables built.")

!cd voxcraft-sim/build; ./voxcraft-sim -i ../demos/basic/ -o output.xml -f > ../../a.history

print("Simulation done.")

!ls

from google.colab import files
files.download('a.history')
```

Here is a [readonly example notebook](https://colab.research.google.com/drive/1yiqw7Uq3W3CgYCinXq4t808M2l7uuLv1?usp=sharing)

## On Your Desktop/Laptop

The most difficult part of compiling this project from scratch is installing the CUDA environment. First, make sure you have NVidia Graphic Cards, then [download CUDA 10.1](https://developer.nvidia.com/cuda-10.1-download-archive-base) and install it.

Once you have the CUDA 10.1 environment, the rest is easy:

```
sudo apt-get update
sudo apt-get install -y git cmake libboost-all-dev

git clone https://github.com/voxcraft/voxcraft-sim.git
cd voxcraft-sim
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCUDA_DEBUG=OFF ..
make -j 10
```

### Notes
* Building been tested on a CMAKE versions 3.12 and newer.
* Versions of boost known to work: 1.65.1 and 1.67.0

Now you will have two executables: `voxcraft-sim` and `vx3_node_worker`. Copy them to your environment and the simulation is ready to use.

Try one of the demos:

```
./voxcraft-sim -i ../demos/basic/ > demo_basic.history
```

It will produce a `demo_basic.history` file that is for [voxcraft-viz](https://github.com/voxcraft/voxcraft-viz) to visualize.

# Citation

If you need to cite our work, here is the format:

```
@MISC{liu_voxcraft_2020,
	title = {Voxcraft-sim, a GPU-accelerated voxel-based physics engine},
	url = {https://github.com/voxcraft/voxcraft-sim},
	howpublished = {\url{https://github.com/voxcraft/voxcraft-sim}},
	author = {Sida Liu and David Matthews and Sam Kriegman and Josh Bongard},
	year = {2020}
	doi = {10.5281/zenodo.3835152},
}
```
[![DOI](https://zenodo.org/badge/265434971.svg)](https://zenodo.org/badge/latestdoi/265434971)
