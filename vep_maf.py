conda activate vep_maf

sed -i 's/\r$//' *.vcf

for vcf in *.vcf; do sample=$(basename "$vcf" .vcf); echo "Processing $sample"; vcf2maf.pl --input-vcf "$vcf" --output-maf "${sample}.maf" --tumor-id "$sample" --ref-fasta hg38.fa --vep-path /home/ssarkar/ensembl-vep --vep-data ~/.vep --ncbi-build GRCh38; done
