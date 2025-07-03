# pedigree-simulator

[IBDsims project by Jean Morrison](https://github.com/jean997/IBDsims)


## Installation

First, click the green `Code` button at the top of this page and select a cloning option.

Dependency and reference data installation takes place using Docker. First, however, you must install and launch the Docker client on your machine. Instructions to install Docker Engine on your system can be found [here](https://docs.docker.com/engine/install/).

Navigate into the project directory cloned from GitHub:

```bash
cd pedigree-simulator
```

Build the Docker image:

```bash
docker build -f Dockerfile.github -t pedigreesim .
```

Note: The build process may take about 10-20 minutes due to the size of the reference data being downloaded (approx. 40 gigabytes after being unzipped).



## Execution

Run (in interactive mode):

```bash
# Step 1: Set entrypoint to bring you into the Docker image location
docker run -v \
    /local/path/to/pedigree_sim_repo/output:/usr/src/output \
    -it --entrypoint /bin/bash pedigreesim:latest 

# Step 2: Run 
perl main.pl 100 uniform3 20 EUR parallel

```
