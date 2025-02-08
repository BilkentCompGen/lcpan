# Pipeline for experiments

## LCPan with different threads.

It is assumed that you have `hg38.fa` reference genome and `pggb.vcf` variant calls file. If you don't please get them from [pggb](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/pggb/vcfs/)

### HG38 experimental results for LCPan on different threads

You can run following script to make thread scaling analyses based on non-overlapping rgfa graph

```sh
for t in 1 2 4 8 16; do cores=$(seq -s, 0 $((t-1))); /bin/time -v taskset -c ${cores} bin/lcpan -r hg38.fa -v pggb.vcf -o hg38.pggb.lcpan.t${t}.rgfa -t $t -p hg38.pggb.lcpan.t${t} --verbose > hg38.pggb.lcpan.t${t}.out 2>&1; /bin/time -v bash lcpan-merge.sh hg38.pggb.lcpan.t${t}.log >> hg38.pggb.lcpan.t${t}.out 2>&1; done
```

You can make statistical analyses on non-overlapping rgfa graph on the results by running:

```sh
parallel '/bin/time -v python3 misc_utils/analyze.py hg38.pggb.lcpan.t{}.rgfa >> hg38.pggb.lcpan.t{}.out 2>&1' ::: 1 2 4 8 16
```

---

### HG38 experimental results for VG

Vg crashed with large genome size (>10million bp). Hence, you need to process the genome by chunks, and then merge them. Merging parts ensuring the consistency by traversing the sub-graphs to find optimal merging.

Before running the script that divides genome into chunks, you need to index the reference genome. You can run `samtools faidx hg38.fa`.

Now, the following command will chunk the genome into regions (10,000,000) to feed the vg:

```sh
CHUNK_SIZE=10000000; rm -f hg38.chunks.txt; while read -r chr chr_length _; do for ((i = 0; i * CHUNK_SIZE < chr_length; i++)); do start=$((i * CHUNK_SIZE + 1)); end=$(((i + 1) * CHUNK_SIZE)); [[ $end -gt $chr_length ]] && end=$chr_length; echo "${chr}:${start}-${end}" >> hg38.chunks.txt; done; done < ../hg38.fa.fai
```

VG takes compressed and indexed VCF file. So if you don't have that and its index, you can run

```sh
bgzip -c pggb.vcf > pggb.vcf.gz
tabix -p vcf pggb.vcf.gz
```

You can run the following script to process fasta file chunk by chunk in different thread settings:

```sh
for t in 1 2 4 8 16; do cores=$(seq -s, 0 $((t-1))); i=0; /bin/time -v taskset -c ${cores} bash -c "while read -r region; do output_file='chunk.t${t}.'\"\${i}\"'.vg'; vg construct -r hg38.fa -v pggb.vcf.gz -f -R \"\$region\" -t ${t} > \"\$output_file\"; ((i++)); done < hg38.chunks.txt"  > hg38.pggb.vg.t${t}.out 2>&1; done
for t in 1 2 4 8 16; do /bin/time -v vg combine -p chunk.t${t}.*.vg > hg38.pggb.t${t}.vg; rm chunk.t${t}.*.vg; done
```

---

### HG38 experimental results for LCPan on different LCP levels

```sh
for l in 4 5 6 7; do taskset -c 0 /bin/time -v bin/lcpan -r hg38.fa -v pggb.vcf -o hg38.pggb.lcpan.l${l}.rgfa -l ${l} -p hg38.pggb.lcpan.l${l} --verbose > hg38.pggb.lcpan.l${l}.out 2>&1; /bin/time -v bash lcpan-merge.sh hg38.pggb.lcpan.l${l}.log >> hg38.pggb.lcpan.l${l}.out 2>&1; done
```

Calculate stats of GFAs

```sh
parallel '/bin/time -v python3 misc_utils/analyze.py hg38.pggb.lcpan.l{}.rgfa >> hg38.pggb.lcpan.l{}.out 2>&1' ::: 4 5 6 7
```

---

#### Experimental results on overlapping graph

```sh
for l in 4 5 6 7; do taskset -c 0 /bin/time -v bin/lcpan -r hg38.fa -v pggb.vcf -o hg38.pggb.lcpan.l${l}.rgfa -l ${l} -p hg38.pggb.lcpan.l${l} -s --verbose > hg38.pggb.lcpan.l${l}.out 2>&1; /bin/time -v bash lcpan-merge.sh hg38.pggb.lcpan.l${l}.log >> hg38.pggb.lcpan.l${l}.out 2>&1; done
```

Calculate stats of GFAs

```sh
parallel '/bin/time -v python3 misc_utils/analyze.py hg38.pggb.lcpan.l{}.rgfa >> hg38.pggb.lcpan.l{}.out 2>&1' ::: 4 5 6 7
```

---

### LCPan profiling

#### perf

You can run the profiling using `perf` on LCPan as follows. Note that you need a `sudo` privilidge to run this script as perf requires to access core level parameters to make calculations.

```sh
perf_stat_metrics="duration_time,branch-instructions,branch-misses,cache-misses,cache-references,cpu-cycles,instructions,mem-loads,mem-stores,cycles,instructions,branches,faults,migrations,L1-dcache-loads,L1-dcache-stores,L1-dcache-prefetch-misses,L1-dcache-load-misses,LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,LLC-prefetch-misses"; for t in 1 2 4 8 16; do cores=$(seq -s, 0 $((t-1))); sudo perf stat -o "hg38.pggb.lcpan.t${t}.perf-stat.txt" -B -e $perf_stat_metrics  taskset -c ${cores} bin/lcpan -r hg38.fa -v pggb.vcf -o hg38.pggb.lcpan.t${t}.perf.rgfa -t $t -p hg38.pggb.lcpan.t${t}.perf --verbose > hg38.pggb.lcpan.t${t}.perf.out 2>&1; sudo chown "$USER:$USER" "hg38.pggb.lcpan.t${t}.perf-stat.txt"; for ((i=1; i<=${t}; i++)); do sudo chown "$USER:$USER" hg38.pggb.lcpan.t${t}.perf.rgfa.${i}; sudo chown "$USER:$USER" hg38.pggb.lcpan.t${t}.log; done; sudo chown "$USER:$USER" hg38.pggb.lcpan.t${t}.perf.rgfa; /bin/time -v bash lcpan-merge.sh hg38.pggb.lcpan.t${t}.perf.log >> hg38.pggb.lcpan.t${t}.perf.out 2>&1; done
```

If you want to profile the LCPan with gprof, please modify Makefiles by removing optimization (from -O3 to -O0) and adding `-pg` flag in `lcpan` and `lcptools`. Then, please run:

#### gprof

```sh
for t in 1 2 4 8 16; do cores=$(seq -s, 0 $((t-1))); /bin/time -v taskset -c ${cores} bin/lcpan -r hg38.fa -v pggb.vcf -o hg38.pggb.lcpan.t${t}.gprof.rgfa -t $t -p hg38.pggb.lcpan.t${t}.gprof --verbose > hg38.pggb.lcpan.t${t}.gprof.out 2>&1; mv gmon.out gmon.out.t${t}; /bin/time -v bash lcpan-merge.sh hg38.pggb.lcpan.t${t}.gprof.log >> hg38.pggb.lcpan.t${t}.gprof.out 2>&1; done; for t in 1 2 4 8 16; do gprof bin/lcpan gmon.out.t${t} > hg38.pggb.lcpan.t${t}.gprof; done
```

## Alignment

### Construct GFA with LCP on HG38.chr1

```sh
# construct gfa graph using LCPan
bin/lcpan -r hg38.chr1.fa -v pggb.chr1.vcf -o lcpan.chr1.gfa -p lcpan.chr1 --verbose > lcpan.chr1.out 2>&1

bash lcpan-merge.sh lcpan.chr1.log >> lcpan.chr1.out 2>&1

# construct gfa graph using VG
CHUNK_SIZE=10000000; rm -f hg38.chr1.chunks.txt; while read -r chr chr_length _; do for ((i = 0; i * CHUNK_SIZE < chr_length; i++)); do start=$((i * CHUNK_SIZE + 1)); end=$(((i + 1) * CHUNK_SIZE)); [[ $end -gt $chr_length ]] && end=$chr_length; echo "${chr}:${start}-${end}" >> hg38.chr1.chunks.txt; done; done < ../hg38.chr1.fa.fai

i=0; while read -r region; do output_file="chunk.${i}.vg"; echo "Processing ${region}"; vg construct -r hg38.chr1.fa -v pggb.chr1.vcf.gz -f -R ${region} --threads 1 > ${output_file} || echo "Warning: vg construct failed for ${region}, skipping..."; ((i++)); done < hg38.chr1.chunks.txt

vg combine -p chunk.*.vg > hg38.chr1.vg

# align with GraphAligner
# hifi reads
cores=$(seq -s, 0 39); taskset -c ${cores} /bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g ../lcpan.chr1.gfa -f hg002.chr1.subsampled.hifi.fa -a hg002.lcpan.hifi.gaf -x vg -t 40 > hg002.lcpan.hifi.out 2>&1; taskset -c ${cores} /bin/time -v  ~/tools/GraphAligner/bin/GraphAligner -g ../vg.chr1.gfa -f hg002.chr1.subsampled.hifi.fa -a hg002.vg.hifi.gaf -x vg -t 40 > vg.out 2>&1;

cores=$(seq -s, 70 109); taskset -c ${cores} /bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g ../lcpan.chr1.gfa -f hg002.chr1.subsampled.85.fa -a hg002.lcpan.85.gaf -x vg -t 40 > hg002.lcpan.85.out 2>&1; taskset -c ${cores} /bin/time -v  ~/tools/GraphAligner/bin/GraphAligner -g ../vg.chr1.gfa -f hg002.chr1.subsampled.85.fa -a hg002.vg.85.gaf -x vg -t 40 > vg.out 2>&1;

# simulated reads with %85 accuracy
cores=$(seq -s, 0 39); taskset -c ${cores} /bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g ../lcpan.chr1.gfa -f hg002.chr1.subsampled.85.fa -a hg002.lcpan.85.gaf -x vg -t 40 > hg002.lcpan.85.out 2>&1;
cores=$(seq -s, 70 109); taskset -c ${cores} /bin/time -v  ~/tools/GraphAligner/bin/GraphAligner -g ../vg.chr1.gfa -f hg002.chr1.subsampled.85.fa -a hg002.vg.85.gaf -x vg -t 40 > hg002.vg.85.out 2>&1;

# simulated reads with %90 accuracy
cores=$(seq -s, 0 39); taskset -c ${cores} /bin/time -v ~/tools/GraphAligner/bin/GraphAligner -g ../lcpan.chr1.gfa -f hg002.chr1.subsampled.90.fa -a hg002.lcpan.90.gaf -x vg -t 40 > hg002.lcpan.90.out 2>&1;
cores=$(seq -s, 70 109); taskset -c ${cores} /bin/time -v  ~/tools/GraphAligner/bin/GraphAligner -g ../vg.chr1.gfa -f hg002.chr1.subsampled.90.fa -a hg002.vg.90.gaf -x vg -t 40 > hg002.vg.90.out 2>&1;
```

---

```sh
ln -s ~/programs/lcp_vg/bin/ .
ln -s ~/programs/lcp_vg/lcpan-merge.sh .
ln -s ~/programs/lcp_vg/misc_utils .
ln -s /tmp/lcpan/hg38.fa .
ln -s /tmp/lcpan/hg38.fa.fai .
ln -s /tmp/lcpan/pggb.vcf .
```

---

You can make validation of the graphs, whether links, segments, overlaps matches (this will make sense when you construct the overlapping rgfa graph), by running:

```sh
parallel '/bin/time -v python3 misc_utils/validate.py hg38.pggb.lcpan.t{}.rgfa >> hg38.pggb.lcpan.t{}.out 2>&1' ::: 1 2 4 8 16
```

You can compare GFAs whether they are same or not. No that this script considers segment strings (not ids) and the links in between them.

```sh
parallel -j 4 '/bin/time -v python3 misc_utils/compare-ov2ov.py hg38.pggb.lcpan.t1.rgfa hg38.pggb.lcpan.t{}.rgfa >> hg38.pggb.lcpan.t${t}.out 2>&1' ::: 2 4 8 16
```