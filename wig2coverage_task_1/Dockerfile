FROM ubuntu:16.04

# Install any libraries necessary to run the task
# RUN set -ex \
#     && apt-get update \
#     && apt-get install -y --no-install-recommends \
#         libcairo2-dev \
#         libxt-dev \
#     && rm -rf /var/lib/apt/lists/* \
#     && install2.r -e \
#         Cairo \
#         data.table \
#     && rm -rf /tmp/*
RUN apt-get update && apt-get install -y python sudo dstat
RUN apt-get install -y python-pip wget zlib1g-dev
RUN pip install -U pip
RUN pip install -U setuptools
RUN pip install numpy scipy pysam cython
RUN pip install h5py pandas ArgumentParser
RUN pip install lxml
RUN pip install wiggelen


# Copy the built tool and any supporting files into the image
COPY src /opt/src

# Set the working directory
WORKDIR /opt/src
