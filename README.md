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

# Run LUMPYEXPRESS on a single sample with pre-extracted splitters and discordants
lumpyexpress -P -B sample.bam -S sample.splitters.bam -D sample.discordants.bam -o sample.vcf
    
# Genotyping individual samples with SVTyper
svtyper -B sample.bam -i sample.vcf -l sample.json > sample.gt.vcf

# Sorting,compressing and indexing VCF files
vcf-sort sample.gt.vcf > sample_sorted.gt.vcf
bgzip sample_sorted.gt.vcf
tabix sample_sorted.gt.vcf.gz

# Run SVTOOLS LSORT to combine and sort variants from multiple samples
svtools lsort sample1_sorted.gt.vcf.gz sample2_sorted.gt.vcf.sample3_sorted.gt.vcf.gz \
| bgzip -c > sorted.vcf.gz

# Run SVTOOLS LMERGE to merge variant calls likely representing same variant in the sorted VCF
zcat sorted.vcf.gz | svtools lmerge -i /dev/stdin -f 20 | bgzip -c > merged.vcf.gz

# Run SVTOOLS GENOTYPE to genotype all samples present in the merged set
mkdir -p gt
zcat merged.vcf.gz | vawk --header '{  $6="."; print }' | svtools genotype -B NA12877.bam -l NA12877.bam.json \
| sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - > gt/NA12877.vcf

```
