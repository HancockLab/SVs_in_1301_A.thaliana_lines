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


## Pipeline

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

##### LUMPY Express
- Run LUMPY Express on a single sample with pre-extracted splitters and discordants
    ```
    lumpyexpress \
        -B sample.bam \
        -S sample.splitters.bam \
        -D sample.discordants.bam \
        -o sample.vcf
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
