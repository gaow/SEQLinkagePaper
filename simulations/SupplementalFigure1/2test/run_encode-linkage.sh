vcf=../0vcf/demo.vcf.gz
tfam=../1tfam/demo.tfam
seqlink --output-entries 20 \
		--bin 3 \
		--theta-inc 0.05 \
		--theta-max 0.5 \
		-K 0.001 \
		--wt-pen 0 \
		--mut-pen 1 \
		--moi AR \
		--fam $tfam \
		--vcf $vcf \
		--jobs 1 \
		--blueprint bp.txt \
		--output demo1c \
		--format LINKAGE \
		--run-linkage \
		--tempdir tmp
		#--vanilla
#slinco --wild_pen 0 --muta_pen 1 --moi AR --tfam $tfam --vcf $vcf --jobs 1 --blueprint bp.txt --size 3 --output demo2 --format mlink --runner mlink #--vanilla
#slinco --wild_pen 0 --muta_pen 1 --moi AR --tfam $tfam --vcf $vcf --jobs 1 --blueprint bp.txt --size 0 --output demo3 --format mlink --runner mlink #--vanilla
