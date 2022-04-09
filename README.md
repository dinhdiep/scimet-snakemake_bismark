# scimet-snakemake_bismark


## View the pipeline

```bash
snakemake --forceall --dag | dot -Tpng > dag.png
```

## Preliminary steps

Make sure that [Miniconda 3] is installed. When prompted, answer yes to add conda to the PATH environment variable.
```bash
   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
```

Download scimet-snakemake_bismark and sciMETv2
```bash
   git clone https://github.com/hdinhdp/scimet-snakemake_bismark
   git clone https://github.com/adeylab/sciMETv2
```

Go to the  directory
```bash
   cd scimet-snakemake_bismark
```

Use the yaml file to install all required software into an isolated Conda environment with the name **bisread-sra**
```bash
   conda env create --name scimet-snakemake --file environment.yml
```

Activate the environment
```bash
   conda activate scimet-snakemake
```

Note: to exit the environment, just execute the following command. Don't do this until you are done working.
```bash
   conda deactivate
```

## Generating samples.json and configuring config.yml

Make the samples.json file with all the samples and the location of their respective read 1 and read 2 files. Modify the parameters in the config.yml file to the correct paths for 'SAMPLES_JSON', 'BOWTIE_DIR', 'BISMARK_REF', and set the desired 'MAP_SCORE' filter.


## Run the snakemake

It is possible to speedup the process by allowing snakemake to spawn multiple jobs in parallel whenever possible by specifying the job number.
```bash 
   snakemake --jobs 8
```

