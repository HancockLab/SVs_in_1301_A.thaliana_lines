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
```
git clone --recursive git@github.com:arq5x/lumpy-sv.git
cd lumpy-sv
make
cp bin/* /usr/local/bin/.
```


## Pipeline

```
# Align paired-end reads to Reference
bwa mem -R "@RG\tID:id\tSM:sample\tLB:lib" reference.fasta sample.R1.fastq sample.R2.fastq \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -Sb - \
    > sample.bam

# Extract the discordant paired-end alignments from BAM file
samtools view -b -F 1294 sample.bam > sample.discordants.unsorted.bam

# Extract the split-read alignments from BAM file
samtools view -h sample.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > sample.splitters.unsorted.bam

# Sort discordant and splitters bam files
samtools sort sample.discordants.unsorted.bam -o sample.discordants.bam
samtools sort sample.splitters.unsorted.bam - o sample.splitters.bam

# Run LUMPYExpress on a single sample with pre-extracted splitters and discordants
lumpyexpress -B sample.bam -S sample.splitters.bam -D sample.discordants.bam -o sample.vcf
    
# Genotyping individual Samples with SVTyper
svtyper -B sample.bam -i sample.vcf -l sample.json > sample.gt.vcf
```
