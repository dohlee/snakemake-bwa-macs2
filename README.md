# snakemake-bwa-macs2

BWA-MACS2 ChIP-seq analysis pipeline in snakemake.

## Quickstart

1. Clone the repo.
```
$ git clone https://github.com/dohlee/snakemake-bwa-macs2.git
$ cd snakemake-bwa-macs2
```

2. Modify configurations in `config.yaml`, sample information in `manifest.csv`, and treatment-control relationship `relation.csv` as you want.

3. Run the pipeline.

If you already have snakemake installed, it might be useful to just use `--use-conda` option. Tweak `-j` parameter according to the number of available cores on your system. Also adjust `network` resource to limit the number of concurrent downloads of SRA files.
```
$ snakemake --use-conda -p -j 32 --resources network=2
```

Or you can create separate environment for this pipeline and run it.
```
$ conda env create -f environment.yaml
$ snakemake -p -j 32 --resources network=2
```

