
# Pedigree Simulator

Pedigree Simulator is an improved version of a tool originally developed by [Staples _et al._ (2014)](https://www.cell.com/ajhg/fulltext/S0002-9297(14)00427-3) to benchmark PRIMUS in its original publication. It also uses code from the [IBDsims program](https://github.com/jean997/IBDsims) developed by [Morrison (2013)](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21737). This version was developed to introduce runtime optimizations and containerization for use in a variety of contexts, notably including the benchmarking of [COMPADRE](https://compadre.dev), a tool that unifies PRIMUS, ERSA, and PADRE.


## Installation

First, click the green `Code` button at the top of this page and select a cloning option.

Dependency and reference data installation takes place using Docker, which must first be [installed and launched](https://docs.docker.com/engine/install/) on your machine.

Navigate into the project directory cloned from GitHub:

```bash
cd pedigree-simulator
```

Build the Docker image:

```bash
docker build -f Dockerfile.github -t pedigreesim .
```

Note: The build process will take between 10-20 minutes due to the size of the reference data being downloaded (approx. 40 gigabytes after being unzipped).



## Execution

After building the Docker image, enter the container using `docker run`. During this step, you should also set your local volume mount (for writing the output) specified by the `-v` flag. Note: Make sure to provide the absolute path to your local output folder in this flag. Even if you run the `docker run` command from inside the top level of the repository folder, you must use `./output` instead of just `output`.

```bash
docker run -v \
    /your/path/to/pedigree-simulator/output:/usr/src/output \
    -it --entrypoint /bin/bash pedigreesim:latest 
```


Once inside the container, you can run the tool from the command line:

```bash
perl main.pl 100 uniform3 20 EUR parallel
```

### Arguments

The `main.pl` script takes several positional arguments. The following descriptions use the above command as an example.

#### Required:
- `100`: The simulation "number", or the unique identifier for the output folder/files
- `uniform3`: The simulation "type." Currently, the script supports `uniform3`, `uniform2`, and `halfsib3`. The key distinction here is that `halfsib3` offers half-sibling relationships in the pedigree. The trailing number represents the average number of offspring per node in the pedigree.  
- `20`: The number of individuals in the pedigree.
- `EUR`: The 1000 Genomes superpopulation from which founder genotypes are drawn. Currently, the script supports EUR (European) and AMR (Admixed American) superpopulation seeding.

#### Optional:
- `parallel` enables parallel processing of the genotype adding step with 22 threads (one per chromosome). This is much faster but RAM intensive and not recommended outside of HPC/server environments. This step will run using a single thread if this argument is not used.


### Notes

- This tool generates a full pedigree as well as incrementally missing versions (up to 20% of all pedigree nodes). This is an artifact of the code left over from our developement done in line with the COMPADRE benchmarking, where we evaluated pedigree reconstruction success as pedigrees became more sparse. If you want to change the maximum % of samples removed in this incremental process, please update the global `$missing_denominator` variable in line 45 of `src/main.pl` _before_ building the Docker image. The default value of 5 divides the total pedigree size by 5, removing 1/5th of all nodes in the last incrementally missing version of the pedigree. If you want more missingness than 20%, consider decreasing the value to 4 or 2, and if you want more, increase it.  



## Questions?

Please email <strong><i>contact AT compadre DOT dev</strong></i> with the subject line "Pedigree Simulator Help" or [submit an issue report/pull request on GitHub](https://github.com/belowlab/pedigree-simulator/issues). 

If you use Pedigree Simulator in your research, please cite the following:
```
Evans GF, Baker JT, Petty LE, Petty AS, Polikowsky HG, Bohlender RJ, Chen HH, Chou CY, 
Viljoen KZ, Beilby JM, Kraft SJ, Zhu W, Landman JM, Morrow AR, Bian D, Scartozzi AC, 
Huff CD, Below JE. COMPADRE: Combined Pedigree-aware Distant Relatedness Estimation 
for improved pedigree reconstruction using integrated relationship estimation approaches 
[Publication details forthcoming]
```


## License

Pedigree Simulator was developed by the [Below Lab](https://thebelowlab.com) in the Division of Genetic Medicine at Vanderbilt University Medical Center, Nashville, TN, USA. 

Pedigree Simulator is distributed under the following APACHE 2.0 license: https://compadre.dev/licenses/sim_license.txt
