help:
	@echo "This is test program for SEQLinkage package. Use auto-complete to view options"
interface:
	@seqlink --vcf DefaultData/test.vcf.gz --fam DefaultData/test.tfam -j8 --freq AF -o InputAF --prevalence 0.01 --moi AD -W 0 -M 1 --run-linkage
	@seqlink --vcf DefaultData/test.vcf.gz --fam DefaultData/test.tfam -j8 -o SampleAF --prevalence 0.01 --moi AD -W 0 -M 1 --run-linkage
missing:
	@echo Linhai
regular:
	@echo Linhai
clean:
	rm -rf cache
