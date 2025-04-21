### Variant calling experiment on chr22
echo "Variant calling experiment on hg002.chr22 alignments ..."

mkdir sv-eval-chr22
cd sv-eval-chr22

## LCPan gaf to sam convertion
echo "LCPan gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.lcpan.gfa \
    ../alignment-chr22/hg002.chr22.lcpan.hifi.gaf \
    ../../reads/hg002.chr22.hifi.fa \
    hg002.chr22.lcpan.hifi.sam \
    --simple > hg002.chr22.lcpan.hifi.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.lcpan.hifi.sam \
    hg002.chr22.lcpan.hifi.qual.sam >> hg002.chr22.lcpan.hifi.out 2>&1
samtools view -bS hg002.chr22.lcpan.hifi.qual.sam | samtools sort -o hg002.chr22.lcpan.hifi.bam && samtools index hg002.chr22.lcpan.hifi.bam >> hg002.chr22.lcpan.hifi.out 2>&1
rm hg002.chr22.lcpan.hifi.sam
rm hg002.chr22.lcpan.hifi.qual.sam

## VG gaf to sam convertion
echo "VG gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.vg.gfa \
    ../alignment-chr22/hg002.chr22.vg.hifi.gaf \
    ../../reads/hg002.chr22.hifi.fa \
    hg002.chr22.vg.hifi.sam \
    --simple > hg002.chr22.vg.hifi.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.vg.hifi.sam \
    hg002.chr22.vg.hifi.qual.sam >> hg002.chr22.vg.hifi.out 2>&1
samtools view -bS hg002.chr22.vg.hifi.qual.sam | samtools sort -o hg002.chr22.vg.hifi.bam && samtools index hg002.chr22.vg.hifi.bam >> hg002.chr22.vg.hifi.out 2>&1
rm hg002.chr22.vg.hifi.sam
rm hg002.chr22.vg.hifi.qual.sam

## LCPan Sim (95) gaf to sam convertion
echo "LCPan Sim (95) gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.lcpan.gfa \
    ../alignment-chr22/hg002.chr22.lcpan.95.gaf \
    ../../reads/hg002.chr22.95.fa \
    hg002.chr22.lcpan.95.sam \
    --simple > hg002.chr22.lcpan.95.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.lcpan.95.sam \
    hg002.chr22.lcpan.95.qual.sam >> hg002.chr22.lcpan.95.out 2>&1
samtools view -bS hg002.chr22.lcpan.95.qual.sam | samtools sort -o hg002.chr22.lcpan.95.bam && samtools index hg002.chr22.lcpan.95.bam >> hg002.chr22.lcpan.95.out 2>&1
rm hg002.chr22.lcpan.95.sam
rm hg002.chr22.lcpan.95.qual.sam

## VG Sim (95) gaf to sam convertion
echo "VG Sim (95) gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.vg.gfa \
    ../alignment-chr22/hg002.chr22.vg.95.gaf \
    ../../reads/hg002.chr22.95.fa \
    hg002.chr22.vg.95.sam \
    --simple > hg002.chr22.vg.95.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.vg.95.sam \
    hg002.chr22.vg.95.qual.sam >> hg002.chr22.vg.95.out 2>&1
samtools view -bS hg002.chr22.vg.95.qual.sam | samtools sort -o hg002.chr22.vg.95.bam && samtools index hg002.chr22.vg.95.bam >> hg002.chr22.vg.95.out 2>&1
rm hg002.chr22.vg.95.sam
rm hg002.chr22.vg.95.qual.sam

## LCPan Sim (90) gaf to sam convertion
echo "LCPan Sim (90) gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.lcpan.gfa \
    ../alignment-chr22/hg002.chr22.lcpan.90.gaf \
    ../../reads/hg002.chr22.90.fa \
    hg002.chr22.lcpan.90.sam \
    --simple > hg002.chr22.lcpan.90.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.lcpan.90.sam \
    hg002.chr22.lcpan.90.qual.sam >> hg002.chr22.lcpan.90.out 2>&1
samtools view -bS hg002.chr22.lcpan.90.qual.sam | samtools sort -o hg002.chr22.lcpan.90.bam && samtools index hg002.chr22.lcpan.90.bam >> hg002.chr22.lcpan.90.out 2>&1
rm hg002.chr22.lcpan.90.sam
rm hg002.chr22.lcpan.90.qual.sam

## VG Sim (90) gaf to sam convertion
echo "VG Sim (90) gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.vg.gfa \
    ../alignment-chr22/hg002.chr22.vg.90.gaf \
    ../../reads/hg002.chr22.90.fa \
    hg002.chr22.vg.90.sam \
    --simple > hg002.chr22.vg.90.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.vg.90.sam \
    hg002.chr22.vg.90.qual.sam >> hg002.chr22.vg.90.out 2>&1
samtools view -bS hg002.chr22.vg.90.qual.sam | samtools sort -o hg002.chr22.vg.90.bam && samtools index hg002.chr22.vg.90.bam >> hg002.chr22.vg.90.out 2>&1
rm hg002.chr22.vg.90.sam
rm hg002.chr22.vg.90.qual.sam

## LCPan Sim (85) gaf to sam convertion
echo "LCPan Sim (85) gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.lcpan.gfa \
    ../alignment-chr22/hg002.chr22.lcpan.85.gaf \
    ../../reads/hg002.chr22.85.fa \
    hg002.chr22.lcpan.85.sam \
    --simple > hg002.chr22.lcpan.85.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.lcpan.85.sam \
    hg002.chr22.lcpan.85.qual.sam >> hg002.chr22.lcpan.85.out 2>&1
samtools view -bS hg002.chr22.lcpan.85.qual.sam | samtools sort -o hg002.chr22.lcpan.85.bam && samtools index hg002.chr22.lcpan.85.bam >> hg002.chr22.lcpan.85.out 2>&1
rm hg002.chr22.lcpan.85.sam
rm hg002.chr22.lcpan.85.qual.sam

## VG Sim (85) gaf to sam convertion
echo "VG Sim (85) gaf2bam"
/bin/time -v ../../akhal gaf2sam \
    ../alignment-chr22/hg38.chr22.pggb.vg.gfa \
    ../alignment-chr22/hg002.chr22.vg.85.gaf \
    ../../reads/hg002.chr22.85.fa \
    hg002.chr22.vg.85.sam \
    --simple > hg002.chr22.vg.85.out 2>&1
/bin/time -v ../../akhal sampoke \
    ../../hg38.chr22.fa \
    hg002.chr22.vg.85.sam \
    hg002.chr22.vg.85.qual.sam >> hg002.chr22.vg.85.out 2>&1
samtools view -bS hg002.chr22.vg.85.qual.sam | samtools sort -o hg002.chr22.vg.85.bam && samtools index hg002.chr22.vg.85.bam >> hg002.chr22.vg.85.out 2>&1
rm hg002.chr22.vg.85.sam
rm hg002.chr22.vg.85.qual.sam

## Run SV calling with pbsv

source ~/scripts/activate_conda

echo "LCPan Sim - (85)"

pbsv discover hg002.chr22.lcpan.85.bam hg002.chr22.lcpan.85.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.lcpan.85.svsig.gz hg002.chr22.lcpan.85.vcf
rm hg002.chr22.lcpan.85.svsig.gz

echo "VG Sim - (85)"

pbsv discover hg002.chr22.vg.85.bam hg002.chr22.vg.85.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.vg.85.svsig.gz hg002.chr22.vg.85.vcf
rm hg002.chr22.vg.85.svsig.gz

echo "LCPan Sim - (90)"

pbsv discover hg002.chr22.lcpan.90.bam hg002.chr22.lcpan.90.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.lcpan.90.svsig.gz hg002.chr22.lcpan.90.vcf
rm hg002.chr22.lcpan.90.svsig.gz

echo "VG Sim - (90)"

pbsv discover hg002.chr22.vg.90.bam hg002.chr22.vg.90.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.vg.90.svsig.gz hg002.chr22.vg.90.vcf
rm hg002.chr22.vg.90.svsig.gz

echo "LCPan Sim - (95)"

pbsv discover hg002.chr22.lcpan.95.bam hg002.chr22.lcpan.95.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.lcpan.95.svsig.gz hg002.chr22.lcpan.95.vcf
rm hg002.chr22.lcpan.95.svsig.gz

echo "VG Sim - (95)"

pbsv discover hg002.chr22.vg.95.bam hg002.chr22.vg.95.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.vg.95.svsig.gz hg002.chr22.vg.95.vcf
rm hg002.chr22.vg.95.svsig.gz

echo "LCPan HiFi"

pbsv discover hg002.chr22.lcpan.hifi.bam hg002.chr22.lcpan.hifi.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.lcpan.hifi.svsig.gz hg002.chr22.lcpan.hifi.vcf
rm hg002.chr22.lcpan.hifi.svsig.gz

echo "VG HiFi"

pbsv discover hg002.chr22.vg.hifi.bam hg002.chr22.vg.hifi.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.chr22.vg.hifi.svsig.gz hg002.chr22.vg.hifi.vcf
rm hg002.chr22.vg.hifi.svsig.gz

## Generate Gold Standard
minimap2 -ax map-hifi ../../hg38.chr22.fa ../../reads/hg002.chr22.hifi.fa -t 64 > hg002.gold_standard.sam
../../akhal sampoke ../../hg38.chr22.fa hg002.gold_standard.sam hg002.gold_standard.qual.sam
grep -v "[0-9]H" hg002.gold_standard.qual.sam > hg002.gold_standard.filtered.sam
samtools view -bS hg002.gold_standard.filtered.sam | samtools sort -o hg002.gold_standard.filtered.bam && samtools index hg002.gold_standard.filtered.bam
rm hg002.gold_standard.qual.sam
rm hg002.gold_standard.filtered.sam

pbsv discover hg002.gold_standard.filtered.bam hg002.gold_standard.filtered.svsig.gz
pbsv call ../../hg38.chr22.fa hg002.gold_standard.filtered.svsig.gz hg002.gold_standard.vcf
rm hg002.gold_standard.filtered.svsig.gz

## Variant Calling comparison

vcf2bed < hg002.gold_standard.vcf > hg002.gold_standard.bed
vcf2bed < hg002.chr22.lcpan.hifi.vcf > hg002.chr22.lcpan.hifi.bed
vcf2bed < hg002.chr22.vg.hifi.vcf > hg002.chr22.vg.hifi.bed
vcf2bed < hg002.chr22.lcpan.95.vcf > hg002.chr22.lcpan.95.bed
vcf2bed < hg002.chr22.vg.95.vcf > hg002.chr22.vg.95.bed
vcf2bed < hg002.chr22.lcpan.90.vcf > hg002.chr22.lcpan.90.bed
vcf2bed < hg002.chr22.vg.90.vcf > hg002.chr22.vg.90.bed
vcf2bed < hg002.chr22.lcpan.85.vcf > hg002.chr22.lcpan.85.bed
vcf2bed < hg002.chr22.vg.85.vcf > hg002.chr22.vg.85.bed

awk '{print $1"\t"$2}' ../../hg38.chr22.fa.fai > hg38.chr22.txt

../../bedtools slop -i hg002.gold_standard.bed -g hg38.chr22.txt -b 100 > hg002.gold_standard.expanded.bed

../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.lcpan.hifi.bed -v > lcpan.hifi.bed
../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.vg.hifi.bed -v > vg.hifi.bed
../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.lcpan.95.bed -v > lcpan.95.bed
../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.vg.95.bed -v > vg.95.bed
../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.lcpan.90.bed -v > lcpan.90.bed
../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.vg.90.bed -v > vg.90.bed
../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.lcpan.85.bed -v > lcpan.85.bed
../../bedtools intersect -a hg002.gold_standard.expanded.bed -b hg002.chr22.vg.85.bed -v > vg.85.bed

GOLD=hg002.gold_standard.expanded.bed
GOLD_COUNT=$(wc -l < $GOLD)

METHODS=("lcpan.hifi" "vg.hifi" "lcpan.95" "vg.95" "lcpan.90" "vg.90" "lcpan.85" "vg.85")

echo -e "Method\tTP\tFP\tFN\tPrecision\tRecall\tF1"

for METHOD in "${METHODS[@]}"; do
    PRED="hg002.chr22.${METHOD}.bed"
    FN_FILE="${METHOD}.bed"
  
    FN=$(wc -l < $FN_FILE)
    TP=$((GOLD_COUNT - FN))

    FP=$(bedtools intersect -v -a "$PRED" -b "$GOLD" | wc -l)

    PRECISION=$(echo "scale=5; $TP / ($TP + $FP)" | bc)
    RECALL=$(echo "scale=5; $TP / ($TP + $FN)" | bc)
    F1=$(echo "scale=5; 2 * $PRECISION * $RECALL / ($PRECISION + $RECALL)" | bc)

    echo -e "${METHOD}\t${TP}\t${FP}\t${FN}\t0${PRECISION}\t0${RECALL}\t0${F1}"
done

conda deactivate