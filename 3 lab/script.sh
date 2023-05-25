rm -f ali.sam

./minimap2 -a -t 6 GCF_000005845.2_ASM584v2_genomic.fna SRR24658890.fasta > ali.sam
echo minimap2

QUALITY=$(samtools flagstat ali.sam | python3 -c 'from sys import stdin; f=stdin.read(); import re; print(re.findall("\d+.\d+%", f)[0].replace("%", ""))')
echo QUALITY=$QUALITY

CMP=$(python3 -c "print(int(float($QUALITY) > 90))")
echo CMP=$CMP

rm -f ali.bam
samtools view -S -b ali.sam > ali.bam

rm -f ali.sorted.bam
samtools sort ali.bam  -o ali.sorted.bam

freebayes -f GCF_000005845.2_ASM584v2_genomic.fna ali.sorted.bam > var.vcf