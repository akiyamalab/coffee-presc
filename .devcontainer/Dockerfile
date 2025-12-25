FROM kyanagis/restretto:1.1

# Set build arguments
ARG PYTHON_VERSION=3.12
ARG MINICONDA_VERSION=24.7.1-0

# Install Miniconda
RUN PYTHON_VERSION_NO_DOT=$(echo ${PYTHON_VERSION} | tr -d .) && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py${PYTHON_VERSION_NO_DOT}_${MINICONDA_VERSION}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set PATH environment variable to include conda
ENV PATH=/opt/conda/bin:$PATH

# Install Python packages
RUN pip install "numpy>=2.0.0"
RUN pip install "rdkit==2024.3.5" rdkit-stubs
RUN pip install h5py pandas scipy
RUN pip install ipykernel isort black
RUN pip install mypy pytest
