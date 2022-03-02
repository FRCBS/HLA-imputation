
# all VCFs into one file
./src/vcf_compressing.sh
bcftools concat -a ./results/R9_*_imputed.vcf.gz > ./tmp/tmp.vcf
bcftools sort -O z ./tmp/tmp.vcf > ./results/R9_imputed_HLAs_v2.vcf.gz
tabix -h ./results/R9_imputed_HLAs_v2.vcf.gz

# to plink dosage
plink2 --vcf ./results/R9_imputed_HLAs_v2.vcf.gz dosage=GP --make-pgen --out ./results/R9_imputed_HLAs
plink2 --pfile ./results/R9_imputed_HLAs --export A --out ./results/R9_imputed_HLAs


