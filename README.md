# voxcraft-sim

_A GPU accelerated voxel-based physics engine._

**Documentation:** https://voxcraft.readthedocs.io/en/latest/

**More information:** https://voxcraft.github.io/design


# Installation

## Docker (recommended)

#### Requirements
When building voxcraft-sim the makefile checks if a GPU is available. Thus it is necessary for docker build to be able to see your GPU. To that end install and configure the [nvidia-container-runtime](https://stackoverflow.com/questions/59691207/docker-build-with-nvidia-runtime).

#### Installing nvidia runtime

```
distribution=$(. /etc/os-release;echo $ID$VERSION_ID) \
   && curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add - \
   && curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | sudo tee /etc/apt/sources.list.d/nvidia-docker.list
   
sudo apt-get update && sudo apt-get install -y nvidia-docker2
sudo systemctl restart docker
```

#### Build
```bash
cd voxcraft-sim
docker build -t voxcraft-sim .
```

#### Run
```bash
docker run -it --gpus all voxcraft-sim
```


## Google Colab

[Google Colab](https://colab.research.google.com/) provides a free online GPU environment.

Create a new notebook and go to Menu->Runtime->Change runtime type, select GPU.

Then, run the script:
```python
!git clone https://github.com/voxcraft/voxcraft-sim.git; cd voxcraft-sim/;

print("Source code downloaded.")

!cd voxcraft-sim; rm build -rf; mkdir build; cd build; cmake ..; make -j 10;

print("Executables built.")

!cd voxcraft-sim/build; ./voxcraft-sim -i ../demos/basic/ -o output.xml -f > ../../a.history

print("Simulation done.")

!ls

from google.colab import files
files.download('a.history')
```

Here is a [readonly example notebook](https://colab.research.google.com/drive/1yiqw7Uq3W3CgYCinXq4t808M2l7uuLv1?usp=sharing)

## Local install

The most difficult part of compiling this project from scratch is installing the CUDA environment. First, make sure you have NVidia Graphic Cards, then [download CUDA 10.1](https://developer.nvidia.com/cuda-10.1-download-archive-base) and install it.

Once you have the CUDA 10.1 environment, the rest is easy:

```
sudo apt-get update
sudo apt-get install -y git cmake libboost-all-dev

git clone https://github.com/voxcraft/voxcraft-sim.git
cd voxcraft-sim
mkdir build
cd build
cmake ..
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
