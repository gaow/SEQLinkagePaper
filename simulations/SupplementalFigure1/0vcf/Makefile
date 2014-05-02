all: demo.vcf.gz
	tabix -p vcf demo.vcf.gz
demo.vcf.gz: demo.vcf
	bgzip -c demo.vcf > demo.vcf.gz
demo.vcf: FigS1.tped
	./convert.pl FigS1.tped > demo.vcf
clean:
	rm -rf demo.vcf*
.PHONY: clean all
