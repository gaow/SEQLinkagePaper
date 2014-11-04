scp:
	scp RunGeneHunter.py csg:/mnt/nfs/zhangdi/data/0archive/seqlink/3revision1/2multipoint/genehunter/2run/
	scp LinkagePowerCalc.py csg:/mnt/nfs/zhangdi/data/0archive/seqlink/3revision1/3gene_with_recombinations/
run:
	qsub SubmitJobs.pbs
clean:
	rm -rf RUN* info/* cache/
test:
	PBS_O_WORKDIR=. PBS_ARRAYID=6 bash SubmitJobs.pbs
.PHONY: run clean test
