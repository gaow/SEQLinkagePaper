run:
	qsub SubmitJobs.pbs
clean:
	rm -rf RUN* FAM* info/* B* cache/
test:
	PBS_O_WORKDIR=. PBS_ARRAYID=6 bash SubmitJobs.pbs
.PHONY: run clean test
