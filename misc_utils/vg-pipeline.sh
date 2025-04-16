# Pipeline for experiments (VG)

## NOTE: It is assumed that you have `hg38.fa` reference genome and `pggb.vcf` variant calls files. 
## NOTE: If you don't please get them from https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/pggb/vcfs/
## NOTE: It is assumed that you have `hg38.chr22.fa` reference genome and `pggb.chr22.vcf` variant call files.
## NOTE: If you don't, extract chromosome 22 from hg38.fa and run `bcftools view pggb.vcf --regions chr22 > pggb.chr22.vcf`
## NOTE: For the alignment, you have to have following files in read directory:
##          - reads/hg002.chr22.85.fa
##          - reads/hg002.chr22.90.fa
##          - reads/hg002.chr22.95.fa
##          - reads/hg002.chr22.hifi.fa
## NOTE: HiFi reads can be get from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/ 
## NOTE: HG002 reference (diploid) genome can be get from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz 
## Processing step for reads:
##  HiFi:
##      samtools view -b hg002.hifi.sorted.bam chr22 | samtools fastq - > hg002.chr22.hifi.fq
##      samtools view -u -f 4 hg002.hifi.sorted.bam | samtools fastq - > hg002.hifi.unmapped.fq
##      cat hg002.hifi.unmapped.fq >> hg002.chr22.hifi.fq
##      sed -n '1~4s/^@/>/p;2~4p' hg002.chr22.hifi.fq > hg002.chr22.hifi.fa
##      rm hg002.chr22.hifi.fq hg002.hifi.unmapped.fq
##  Sim:
##      pbsim --strategy wgs --method qshmm --qshmm data/QSHMM-RSII.model --depth 15 --genome ../hg002.chr22.fa --prefix hg002.chr22.85 --accuracy-mean 0.85 > hg002.chr22.85.out 2>&1
##      cat hg002.chr22.85_0002.fastq >> hg002.chr22.85_0001.fastq && rm hg002.chr22.85_0002.fastq
##      cat hg002.chr22.85_0002.maf >> hg002.chr22.85_0001.maf && rm hg002.chr22.85_0002.maf
##      sed -n '1~4s/^@/>/p;2~4p' hg002.chr22.85_0001.fastq > hg002.chr22.85.fa
##      pbsim --strategy wgs --method qshmm --qshmm data/QSHMM-RSII.model --depth 15 --genome ../hg002.chr22.fa --prefix hg002.chr22.90 --accuracy-mean 0.90 > hg002.chr22.90.out 2>&1
##      cat hg002.chr22.90_0002.fastq >> hg002.chr22.90_0001.fastq && rm hg002.chr22.90_0002.fastq
##      cat hg002.chr22.90_0002.maf >> hg002.chr22.90_0001.maf && rm hg002.chr22.90_0002.maf
##      sed -n '1~4s/^@/>/p;2~4p' hg002.chr22.90_0001.fastq > hg002.chr22.90.fa
##      pbsim --strategy wgs --method qshmm --qshmm data/QSHMM-RSII.model --depth 15 --genome ../hg002.chr22.fa --prefix hg002.chr22.95 --accuracy-mean 0.95 > hg002.chr22.95.out 2>&1
##      cat hg002.chr22.95_0002.fastq >> hg002.chr22.95_0001.fastq && rm hg002.chr22.95_0002.fastq
##      cat hg002.chr22.95_0002.maf >> hg002.chr22.95_0001.maf && rm hg002.chr22.95_0002.maf
##      sed -n '1~4s/^@/>/p;2~4p' hg002.chr22.95_0001.fastq > hg002.chr22.95.fa
##      rm *.ref *.fastq
## NOTE: since HG002 is diploid genome, depth 15 means depth 30 of reads in simulation.
##
## The required data and program executables structure should be as follows:
##          ├── akhal
##          ├── lcpan
##          ├── lcpan-merge.sh
##          ├── hg38.chr22.fa
##          ├── hg38.chr22.fa.fai
##          ├── hg38.fa
##          ├── hg38.fa.fai
##          ├── pggb.chr22.vcf
##          ├── pggb.vcf
##          └── reads
##              ├── hg002.chr22.85.fa
##              ├── hg002.chr22.90.fa
##              ├── hg002.chr22.95.fa
##              └── hg002.chr22.hifi.fa
##
## NOTE: samtools and bcftools are required (system wide installed) for variant calling

### Thread scaling analyses based on non-overlapping gfa graph

echo "Experiment on different thread numbers started (nov-gfa)..."

mkdir lcpan-threads
cd lcpan-threads

for t in 1 2 4 8 16; do \
    /bin/time -v ../lcpan -vg \
        -r ../hg38.fa \
        -v ../pggb.vcf \
        -p hg38.pggb.lcpan.t${t} \
        -t $t \
        --gfa \
        --verbose > hg38.pggb.lcpan.t${t}.out 2>&1; \
    /bin/time -v bash ../lcpan-merge.sh hg38.pggb.lcpan.t${t}.log >> hg38.pggb.lcpan.t${t}.out 2>&1; \
done
parallel '/bin/time -v ../akhal parse hg38.pggb.lcpan.t{}.gfa >> hg38.pggb.lcpan.t{}.out 2>&1' ::: 1 2 4 8 16
stat --format="%s %n" *.gfa 2>/dev/null | tee sizes.txt > /dev/null && rm -f *.gfa

cd ..

### LCP levels analyses based on non-overlapping rgfa graph

echo "Experiment on different LCP levels started (nov-rgfa)..."

mkdir lcpan-levels-nov-rgfa
cd lcpan-levels-nov-rgfa

for l in 4 5 6 7; do \
    /bin/time -v ../lcpan -vg \
        -r ../hg38.fa \
        -v ../pggb.vcf \
        -p hg38.pggb.lcpan.l${l} \
        -l ${l} \
        --verbose > hg38.pggb.lcpan.l${l}.out 2>&1; \
    /bin/time -v bash ../lcpan-merge.sh hg38.pggb.lcpan.l${l}.log >> hg38.pggb.lcpan.l${l}.out 2>&1; \
done
parallel '/bin/time -v ../akhal parse hg38.pggb.lcpan.l{}.rgfa >> hg38.pggb.lcpan.l{}.out 2>&1' ::: 4 5 6 7
parallel '/bin/time -v ../akhal stats hg38.pggb.lcpan.l{}.rgfa >> hg38.pggb.lcpan.l{}.out 2>&1' ::: 4 5 6 7
stat --format="%s %n" *.rgfa 2>/dev/null | tee sizes.txt > /dev/null && rm -f *.rgfa

cd ..

### LCP levels analyses based on non-overlapping gfa graph

echo "Experiment on different LCP levels started (nov-gfa)..."

mkdir lcpan-levels-nov-gfa
cd lcpan-levels-nov-gfa

for l in 4 5 6 7; do \
    /bin/time -v ../lcpan -vg \
        -r ../hg38.fa \
        -v ../pggb.vcf \
        -p hg38.pggb.lcpan.l${l} \
        -l ${l} \
        --verbose \
        --gfa > hg38.pggb.lcpan.l${l}.out 2>&1; \
    /bin/time -v bash ../lcpan-merge.sh hg38.pggb.lcpan.l${l}.log >> hg38.pggb.lcpan.l${l}.out 2>&1; \
done
parallel '/bin/time -v ../akhal parse hg38.pggb.lcpan.l{}.gfa >> hg38.pggb.lcpan.l{}.out 2>&1' ::: 4 5 6 7
parallel '/bin/time -v ../akhal stats hg38.pggb.lcpan.l{}.gfa >> hg38.pggb.lcpan.l{}.out 2>&1' ::: 4 5 6 7
stat --format="%s %n" *.gfa 2>/dev/null | tee sizes.txt > /dev/null && rm -f *.gfa

cd ..

### VG construction analyses with chuning principle

echo "Experiment on VG tool graph construction ..."

mkdir vg-construction
cd vg-construction

CHUNK_SIZE=10000000
rm -f hg38.chunks.txt

while read -r chr chr_length _; do \
	for ((i = 0; i * CHUNK_SIZE < chr_length; i++)); do \
		start=$((i * CHUNK_SIZE + 1)); \
		end=$(((i + 1) * CHUNK_SIZE)); \
		[[ $end -gt $chr_length ]] && end=$chr_length; \
		echo "${chr}:${start}-${end}" >> hg38.chunks.txt; \
	done; \
done < ../hg38.fa.fai

/bin/time -v bash -c "i=0; while read -r region; do output_file='chunk.'\"\${i}\"'.vg'; vg construct -r ../hg38.fa -v ../pggb.vcf.gz -f -R \"\$region\" > \"\$output_file\"; ((i++)); done < hg38.chunks.txt" > hg38.pggb.vg.out 2>&1
/bin/time -v vg combine -p chunk.*.vg > hg38.pggb.vg 2>> hg38.pggb.vg.out
/bin/time -v vg convert -f hg38.pggb.vg > hg38.pggb.vg.gfa 2>> hg38.pggb.vg.out
rm chunk.*.vg;
/bin/time -v ../akhal stats hg38.pggb.vg.gfa >> hg38.pggb.vg.out 2>&1
stat --format="%s %n" *.gfa 2>/dev/null | tee sizes.txt > /dev/null && rm -f *.gfa

cd ..

### Alignment experiment on chr22

echo "Experiment on alignment for hg002.chr22 ..."

mkdir alignment-chr22
cd alignment-chr22

#### Create vg graph
echo "VG chr22 graph construction"
/bin/time -v vg construct -r ../hg38.chr22.fa -v ../pggb.chr22.vcf -m 64 > hg38.chr22.pggb.vg 2> hg38.chr22.pggb.vg.out
/bin/time -v vg convert -f hg38.chr22.pggb.vg > hg38.chr22.pggb.vg.gfa 2>> hg38.chr22.pggb.vg.out

#### Create nov-gfa lcpan graph
echo "LCPan chr22 graph construction"
/bin/time -v ../lcpan -vg -r ../hg38.chr22.fa -v ../pggb.chr22.vcf -p hg38.chr22.pggb.lcpan --gfa --verbose > hg38.chr22.pggb.lcpan.out 2>&1
/bin/time -v bash ../lcpan-merge.sh hg38.chr22.pggb.lcpan.log >> hg38.chr22.pggb.lcpan.out 2>&1

#### Align HiFi reads of (hg002.chr22)
echo "HiFi"
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.lcpan.gfa -f ../reads/hg002.chr22.hifi.fa -a hg002.chr22.lcpan.hifi.gaf -x vg -t 32 > hg002.chr22.lcpan.hifi.out 2>&1
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.vg.gfa -f ../reads/hg002.chr22.hifi.fa -a hg002.chr22.vg.hifi.gaf -x vg -t 32 > hg002.chr22.vg.hifi.out 2>&1

## LCPan SNP and indel calling
echo "LCPan Variant Calling"
/bin/time -v ../akhal gaf2sam hg38.chr22.pggb.lcpan.gfa hg002.chr22.lcpan.hifi.gaf ../reads/hg002.chr22.hifi.fa hg002.chr22.lcpan.hifi.sam >> hg002.chr22.lcpan.hifi.out 2>&1
/bin/time -v ../akhal sampoke ../hg38.chr22.fa hg002.chr22.lcpan.hifi.sam hg002.chr22.lcpan.hifi.qual.sam >> hg002.chr22.lcpan.hifi.out 2>&1
rm hg002.chr22.lcpan.hifi.sam
samtools view -bS hg002.chr22.lcpan.hifi.qual.sam | samtools sort -o hg002.chr22.lcpan.hifi.qual.bam && samtools index hg002.chr22.lcpan.hifi.qual.bam >> hg002.chr22.lcpan.hifi.out 2>&1
bcftools mpileup -f ../hg38.chr22.fa hg002.chr22.lcpan.hifi.qual.bam | bcftools call -mv -Oz -o hg002.chr22.lcpan.hifi.vcf.gz >> hg002.chr22.lcpan.hifi.out 2>&1

## VG SNP and indel calling
echo "VG Variant Calling"
/bin/time -v ../akhal gaf2sam hg38.chr22.pggb.vg.gfa hg002.chr22.vg.hifi.gaf ../reads/hg002.chr22.hifi.fa hg002.chr22.vg.hifi.sam >> hg002.chr22.vg.hifi.out 2>&1
/bin/time -v ../akhal sampoke ../hg38.chr22.fa hg002.chr22.vg.hifi.sam hg002.chr22.vg.hifi.qual.sam >> hg002.chr22.vg.hifi.out 2>&1
rm hg002.chr22.vg.hifi.sam
samtools view -bS hg002.chr22.vg.hifi.qual.sam | samtools sort -o hg002.chr22.vg.hifi.qual.bam && samtools index hg002.chr22.vg.hifi.qual.bam >> hg002.chr22.vg.hifi.out 2>&1
bcftools mpileup -f ../hg38.chr22.fa hg002.chr22.vg.hifi.qual.bam | bcftools call -mv -Oz -o hg002.chr22.vg.hifi.vcf.gz >> hg002.chr22.vg.hifi.out 2>&1

#### Align PacBio-Sim (85) reads of (hg002.chr22)
echo "PacBio-sim acc:85"
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.lcpan.gfa -f ../reads/hg002.chr22.85.fa -a hg002.chr22.lcpan.85.gaf -x vg -t 32 > hg002.chr22.lcpan.85.out 2>&1
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.vg.gfa -f ../reads/hg002.chr22.85.fa -a hg002.chr22.vg.85.gaf -x vg -t 32 > hg002.chr22.vg.85.out 2>&1

#### Align PacBio-Sim (90) reads of (hg002.chr22)
echo "PacBio-sim acc:90"
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.lcpan.gfa -f ../reads/hg002.chr22.90.fa -a hg002.chr22.lcpan.90.gaf -x vg -t 32 > hg002.chr22.lcpan.90.out 2>&1
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.vg.gfa -f ../reads/hg002.chr22.90.fa -a hg002.chr22.vg.90.gaf -x vg -t 32 > hg002.chr22.vg.90.out 2>&1

#### Align PacBio-Sim (95) reads of (hg002.chr22)
echo "PacBio-sim acc:95"
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.lcpan.gfa -f ../reads/hg002.chr22.95.fa -a hg002.chr22.lcpan.95.gaf -x vg -t 32 > hg002.chr22.lcpan.95.out 2>&1
/bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g hg38.chr22.pggb.vg.gfa -f ../reads/hg002.chr22.95.fa -a hg002.chr22.vg.95.gaf -x vg -t 32 > hg002.chr22.vg.95.out 2>&1

cd ..

# mv results to seperate folder
mkdir lcpan.vg.results
mv lcpan-threads lcpan.vg.results
mv lcpan-levels-nov-rgfa lcpan.vg.results
mv lcpan-levels-nov-gfa lcpan.vg.results
mv vg-construction lcpan.vg.results
mv alignment-chr22 lcpan.vg.results

### gprof

## If you want to profile the LCPan-VG with gprof, please modify Makefiles by removing optimization (from -O3 to -O0) and adding `-pg` flag in `lcpan` and `lcptools`. Then, please run:

# mkdir lcpan-vg-gprof
# cd lcpan-vg-gprof
# 
# for t in 1 2 4 8 16; do \
#     ../lcpan -r hg38.fa -v pggb.vcf -p hg38.pggb.lcpan.t${t}.gprof -t $t --verbose > hg38.pggb.lcpan.t${t}.gprof.out 2>&1; \
#     ../lcpan-merge.sh hg38.pggb.lcpan.t${t}.gprof.log >> hg38.pggb.lcpan.t${t}.gprof.out 2>&1; \
#      mv gmon.out gmon.out.t${t}; \
# done
# rm -f *.gfa
#
# for t in 1 2 4 8 16; do \
#     gprof ../lcpan gmon.out.t${t} > hg38.pggb.lcpan.t${t}.gprof; \
# done
#
# cd ..
#
# mv lcpan-vg-gprof lcpan.vg.results

### perf

## You can run the profiling using `perf` on LCPan as follows. Note that you need a `sudo` privilidge to run this script as perf requires to access core level parameters to make calculations.

# mkdir lcpan-vg-perf
# cd lcpan-vg-perf
#
# perf_stat_metrics="duration_time,branch-instructions,branch-misses,cache-misses,cache-references,cpu-cycles,instructions,mem-loads,mem-stores,cycles,instructions,branches,faults,migrations,L1-dcache-loads,L1-dcache-stores,L1-dcache-prefetch-misses,L1-dcache-load-misses,LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,LLC-prefetch-misses";
#
# for t in 1 2 4 8 16; do \
#     sudo perf stat -o "hg38.pggb.lcpan.t${t}.perf-stat.txt" -B -e $perf_stat_metrics \
#         ../lcpan -r hg38.fa -v pggb.vcf -p hg38.pggb.lcpan.t${t}.perf -t $t --verbose > hg38.pggb.lcpan.t${t}.perf.out 2>&1; \
#     sudo chown "$USER:$USER" "hg38.pggb.lcpan.t${t}.perf-stat.txt"; \
#     for ((i=1; i<=${t}; i++)); do \
#         sudo chown "$USER:$USER" hg38.pggb.lcpan.t${t}.perf.rgfa.${i}; \
#         sudo chown "$USER:$USER" hg38.pggb.lcpan.t${t}.log; \
#     done; \
#     sudo chown "$USER:$USER" hg38.pggb.lcpan.t${t}.perf.rgfa; \
#     /bin/time -v bash ../lcpan-merge.sh hg38.pggb.lcpan.t${t}.perf.log >> hg38.pggb.lcpan.t${t}.perf.out 2>&1; \
# done
# rm -f *.gfa
#
# cd ..
#
# mv lcpan-vg-perf lcpan.vg.results