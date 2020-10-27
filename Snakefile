import glob
import os
import pandas as pd

from pathlib import Path

wildcard_constraints:
    run='SRR[0-9]+'

onsuccess: shell('push ChIP-seq pipeline OK')
onerror: shell('push ChIP-seq pipeline ERROR')

configfile: 'config.yaml'
manifest = pd.read_csv(config['manifest'])
relation = pd.read_csv(config['relation'])

runs = manifest['run_accession'].values
run2lib = {rec.run_accession:rec.library_layout for rec in manifest.to_records()}
single_runs = [r for r in runs if run2lib[r].upper().startswith('SINGLE')]
paired_runs = [r for r in runs if run2lib[r].upper().startswith('PAIRED')]

narrow_treats = [r.treat for r in relation.to_records() if r.peak_type.upper().startswith('NARROW')]
narrow_controls = [r.control for r in relation.to_records() if r.peak_type.upper().startswith('NARROW')]
broad_treats = [r.treat for r in relation.to_records() if r.peak_type.upper().startswith('BROAD')]
broad_controls = [r.control for r in relation.to_records() if r.peak_type.upper().startswith('BROAD')]

DATA_DIR = Path(config['data_dir'])
RESULT_DIR = Path(config['result_dir'])

include: 'rules/download.smk'
include: 'rules/fastqc.smk'
include: 'rules/bwa.smk'
include: 'rules/trim-galore.smk'
include: 'rules/macs2.smk'

# Raw FASTQ files.
RAW_FASTQ_SINGLE = expand(str(DATA_DIR / '{run}.fastq.gz'), run=single_runs)
RAW_FASTQ_PAIRED = expand(str(DATA_DIR / '{run}.read1.fastq.gz'), run=paired_runs)
RAW_FASTQ_SINGLE_QC = expand(str(DATA_DIR / '{run}_fastqc.zip'), run=single_runs)
RAW_FASTQ_PAIRED_QC = expand(str(DATA_DIR / '{run}.read1_fastqc.zip'), run=paired_runs)

# Quality-controlled FASTQ files.
FASTQ_SINGLE = expand(str(RESULT_DIR / '01_trim_galore' / '{run}.trimmed.fastq.gz'), run=single_runs)
FASTQ_PAIRED = expand(str(RESULT_DIR / '01_trim_galore' / '{run}.read1.trimmed.fastq.gz'), run=paired_runs)
FASTQ_SINGLE_QC = expand(str(RESULT_DIR / '01_trim_galore' / '{run}.trimmed_fastqc.zip'), run=single_runs)
FASTQ_PAIRED_QC = expand(str(RESULT_DIR / '01_trim_galore' / '{run}.read1.trimmed_fastqc.zip'), run=paired_runs)

# Called peaks.
NARROWPEAKS = expand(str(RESULT_DIR / '04_macs2_callpeak' / '{treat}_vs_{control}_peaks.narrowPeak'), zip, treat=narrow_treats, control=narrow_controls)
BROADPEAKS = expand(str(RESULT_DIR / '04_macs2_callpeak' / '{treat}_vs_{control}_peaks.broadPeak'), zip, treat=broad_treats, control=broad_controls)

# bwa alignments.
BAM = expand(str(RESULT_DIR / '02_bwa' / '{run}.sorted.bam'), run=single_runs+paired_runs)

ALL = []
ALL.append(RAW_FASTQ_SINGLE_QC)
ALL.append(RAW_FASTQ_PAIRED_QC)
ALL.append(FASTQ_SINGLE)
ALL.append(FASTQ_PAIRED)
ALL.append(FASTQ_SINGLE_QC)
ALL.append(FASTQ_PAIRED_QC)
ALL.append(BAM)
ALL.append(NARROWPEAKS)
ALL.append(BROADPEAKS)

rule all:
    input: ALL
