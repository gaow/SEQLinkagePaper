#!/usr/bin/perl -w
use strict;
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\tI1\tI2\tII1\tII2\tII3\tII4\n";
while (<>) {
	chomp;
	my ($chr, $id, undef, $pos, @geno) = (split);
	next if $id=~/^[ABCD]$/;
	my $geno = &collapse(\@geno);
	printf "%s\t%s\t%s\t0\t1\t.\t.\t.\tGT\t%s\n", $chr, $pos, $id, join(',', $geno);
}
sub collapse {
	my $genoArray = shift;
	my @genoNew;
	for my $i (0..(@$genoArray/2 - 1)) {
		push @genoNew, join('|', $genoArray->[2 * $i] - 1, $genoArray->[2*$i+1] - 1);
	}
	return join("\t", @genoNew);
}
