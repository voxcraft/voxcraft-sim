FROM nvcr.io/nvidia/cuda:10.2-devel-ubuntu18.04 AS cuda

ARG GID
ARG UID
ENV GID=${GID:-2000}
ENV UID=${UID:-2000}

WORKDIR "/root"

RUN apt-get update
RUN apt-get remove -y --purge cmake
RUN apt-get install -y git libboost-all-dev wget screen

# Clean all apt list
RUN apt-get update
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

# Add a user that aligned to the outer system
RUN addgroup --gid $GID app
RUN adduser --gid $GID --uid $UID --home /home/app --shell /bin/bash --disabled-password --gecos "" app
# Switch to that user
USER app
WORKDIR /home/app
RUN mkdir /home/app/voxcraft-sim

# Install Miniconda
ENV PATH="~/miniconda3/bin:${PATH}"
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir ~/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

# Install cmake for compilation
RUN /home/app/miniconda3/bin/conda install -c anaconda cmake

# # Mount current folder into container
# VOLUME . "/home/app/voxcraft-sim"
# # Compile the source code for the outer system
# RUN cd voxcraft-sim && rm build -rf && mkdir build && cd build && /home/app/miniconda3/bin/cmake .. && make -j 10
# RUN echo "Done"