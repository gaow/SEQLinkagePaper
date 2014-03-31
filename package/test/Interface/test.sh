seqlink --vcf test.vcf.gz --fam test.tfam -j8 --vanilla --freq AF -o InputAF \
    --prevalence 0.01 --moi AD --wild-pen 0 --muta-pen 1
seqlink --vcf test.vcf.gz --fam test.tfam -j8 --vanilla -o SampleAF \
    --prevalence 0.01 --moi AD --wild-pen 0 --muta-pen 1
