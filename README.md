## Installation

##### Requirements
- LUMPY
    * g++ compiler
    * CMake
- LUMPY Express (optional)
    * Samtools (0.1.18+) ([htslib.org/](http://www.htslib.org/))
    * SAMBLASTER (0.1.19+) ([github repo](https://github.com/GregoryFaust/samblaster))
    * Python 2.7 ([python.org/](https://www.python.org/)) with pysam (0.8.3+) and NumPy (1.8.1+)
    * sambamba ([gihub repo](https://github.com/lomereiter/sambamba))
    * gawk ([GNU project](https://www.gnu.org/software/gawk/))

##### Install

Default method to install:

```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```





#### Configuration
LUMPY Express runs several external program whose paths are specified in
[scripts/lumpyexpress.config](scripts/lumpyexpress.config). This config
must reside in the same directory as lumpyexpress, or be specified explicitly
with the -K flag.

The installation Makefile auto-generates a lumpyexpress.config file
and places it in the "bin" directory.

#### Input
LUMPY Express expects BWA-MEM aligned BAM files as input.
It automatically parses sample, library, and read group information using the @RG
tags in the BAM header.
Each BAM file is expected to contain exactly one sample.

The minimum input is a coordinate-sorted BAM file (-B), from which LUMPY Express
extracts splitters and discordants using SAMBLASTER before running LUMPY.
Optionally, users may supply coordinate-sorted splitter (-S) and discordant (-D)
BAM files which will bypass SAMBLASTER extraction for faster analysis.

#### Output
LUMPY Express produces a VCF file according to [VCF spec 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

## LUMPY (traditional) usage
Flexible and customizable breakpoint detection for advanced users.

```
usage:    lumpy [options]
```

**Options**
```
-g       Genome file (defines chromosome order)
-e       Show evidence for each call
-w       File read windows size (default 1000000)
-mw      minimum weight across all samples for a call
-msw     minimum per-sample weight for a call
-tt      trim threshold
-x       exclude file bed file
-t       temp file prefix, must be to a writeable directory
-P       output probability curve for each variant
-b       output as BEDPE instead of VCF

-sr      bam_file:<file name>,
         id:<sample name>,
       	 back_distance:<distance>,
         min_mapping_threshold:<mapping quality>,
         weight:<sample weight>,
         min_clip:<minimum clip length>,
         read_group:<string>

-pe      bam_file:<file name>,
         id:<sample name>,
         histo_file:<file name>,
         mean:<value>,
         stdev:<value>,
         read_length:<length>,
         min_non_overlap:<length>,
         discordant_z:<z value>,
         back_distance:<distance>,
         min_mapping_threshold:<mapping quality>,
         weight:<sample weight>,
         read_group:<string>

-bedpe   bedpe_file:<bedpe file>,
         id:<sample name>,
         weight:<sample weight>
```

## Example workflows

#### Pre-processing
We recommend aligning data with [SpeedSeq](https://github.com/cc2qe/speedseq), which
performs BWA-MEM alignment, marks duplicates and extracts split and discordant
read-pairs.
```
speedseq align -R "@RG\tID:id\tSM:sample\tLB:lib" \
    human_g1k_v37.fasta \
    sample.1.fq \
    sample.2.fq
```

Otherwise, data may be aligned with BWA-MEM.

```
# Align the data
bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" human_g1k_v37.fasta sample.1.fq sample.2.fq \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > sample.bam

# Extract the discordant paired-end alignments.
samtools view -b -F 1294 sample.bam > sample.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h sample.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > sample.splitters.unsorted.bam

# Sort both alignments
samtools sort sample.discordants.unsorted.bam sample.discordants
samtools sort sample.splitters.unsorted.bam sample.splitters
```

#### Running LUMPY
LUMPY has two distinct execution alternatives. LUMPY Express is a simplified wrapper for standard analyses.
LUMPY (traditional) is more customizable, for advanced users and specialized experiments.

##### LUMPY Express
- Run LUMPY Express on a single sample with pre-extracted splitters and discordants
    ```
    lumpyexpress \
        -B sample.bam \
        -S sample.splitters.bam \
        -D sample.discordants.bam \
        -o sample.vcf
    ```

- Run LUMPY Express jointly on multiple samples with pre-extracted splitters and discordants
    ```
    lumpyexpress \
        -B sample1.bam,sample2.bam,sample3.bam \
        -S sample1.splitters.bam,sample2.splitters.bam,sample3.splitters.bam \
        -D sample1.discordants.bam,sample2.discordants.bam,sample3.discordants.bam \
        -o multi_sample.vcf
    ```



#### Post-processing
[SVTyper](https://github.com/hall-lab/svtyper) can call genotypes on LUMPY output VCF files
using a Bayesian maximum likelihood algorithm.
```
svtyper \      
    -B sample.bam \
    -S sample.splitters.bam \
    -i sample.vcf
    > sample.gt.vcf
```
