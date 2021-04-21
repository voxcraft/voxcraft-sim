FROM nvcr.io/nvidia/cuda:10.2-devel-ubuntu18.04 AS cuda

ARG GID
ARG UID
ENV GID=${GID:-2000}
ENV UID=${UID:-2000}

WORKDIR "/root"

RUN apt-get update
RUN apt-get remove -y --purge cmake
RUN apt-get install -y git libboost-all-dev wget screen

# Install Miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

# Build voxcraft-sim
RUN conda install -c anaconda cmake

# Add a user that aligned to the outer system
RUN addgroup -S app --gid $GID \
    && adduser -S -G app -u $UID -h /home/app -D app
# Switch to that user
USER app
RUN mkdir /home/app/voxcraft-sim
VOLUME ["/home/app/voxcraft-sim"]
WORKDIR /home/app
# Compile the source code for the outer system
RUN cd voxcraft-sim && mkdir build && cd build && cmake .. && make -j 10
