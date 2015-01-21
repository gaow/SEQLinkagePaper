seqlink --vcf test.vcf.gz --fam test.tfam -j8 --freq AF -o InputAF \
    --prevalence 0.01 --moi AD -W 0 -M 1 --run-linkage
seqlink --vcf test.vcf.gz --fam test.tfam -j8 -o SampleAF \
    --prevalence 0.01 --moi AD -W 0 -M 1 --run-linkage
