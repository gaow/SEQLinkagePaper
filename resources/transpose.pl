#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(sum);

#parse options
my $tped;
my $tfam;
my $out = 'LINKAGE';
my $format = 'mlink';
GetOptions('p|tped=s' => \$tped,
	   'format=s' => \$format,
	   'f|tfam=s' => \$tfam, 
	   'o|out=s' => \$out)
	   or die "Error in command line arguments\n";

#open files for Input
open TPED, $tped or die "cannot open $tped: $!\n";
open TFAM, $tfam or die "cannot open $tfam: $!\n";


#plink output
if ($format eq 'plink') {
	open PED, ">$out.ped" or die "cannot open $out.ped: $!\n";
	open MAP, ">$out.map" or die "cannot open $out.ped: $!\n";
	my $geno;
	while (<TPED>) {
		chomp;
		my @s = split;
		printf MAP "%s %s %d %d\n", @s[0..3];
		my @g = @s[4..$#s];
		map s/NA/0/, @g;
		push @{$geno}, \@g;
	}
	close TPED or die "cannot close TPED: $!\n";
	close MAP or die "cannot close MAP: $!\n";
	
	while (<TFAM>) {
		chomp;
		s/\s+/ /g;
		print PED;
		map {my $a1 = shift @{$_}; my $a2 = shift @{$_}; print PED " $a1 $a2"} @{$geno};
		print PED "\n";
	}
	close TFAM or die "cannot close TFAM: $!\n";
	close PED or die "cannot close PED: $!\n";
} elsif ($format eq 'mlink') {
	mkdir $out;
	my @fam;
	my @founders;
	while (<TFAM>) {
		chomp;
		s/\s+/ /g;
		push @fam, $_;
		if (m/ 0 0 /) {
			push @founders, 2 * ($. - 1);
			push @founders, 2 * ($. - 1) + 1;
		}
	}
	close TFAM or die "cannot close TFAM: $!\n";
	while (<TPED>) {
		chomp;
		my ($chr, $name, $dis, $pos, @s) = split;
		map s/NA/0/, @s;
		open LOC, ">$out/$name.LOC" or die "cannot open $out/$name.LOC: $!\n";
		my @alleles = @s[@founders];
		my @stat = @{&stat(\@alleles)};
		my $num = shift @stat;
		next if $num == 0;
		print LOC ($chr =~ m/(X|Y|XY)/) ? " 2 0 0 5\n" : " 2 0 0 5\n";
		print LOC " 0 0.0 0.0 0\n";
		print LOC "1 2\n";
		print LOC "1 2\n";
		print LOC " 0.999 0.001\n";
		print LOC " 1\n";
		print LOC " 0.05 0.9  0.9\n";
		print LOC "3 $num\n";
		printf LOC " %f" x $num . "\n", @stat;
		print LOC "0 0\n";
		print LOC "0.0\n";
		print LOC "1 0.05 0.45\n";
		close LOC or die "cannot close LOC: $!\n";
		open PRE, ">$out/$name.PRE" or die "cannot open $out/$name.PRE: $!\n";
		for (0..$#fam) {
			printf PRE "%s %d %d\n", $fam[$_], $s[2 * $_], $s[2 * $_ + 1];
		}
		close PRE or die "cannot close PRE: $!\n";
	}
}
sub stat {
	my $alleles = shift;
	return [0] if scalar(@{$alleles}) == 0;
	my @cnt;
	for my $a (@{$alleles}) {
		$cnt[$a]++;
	}
	for (0..$#cnt) {
		unless ($cnt[$_]) {
			$cnt[$_] = 0;
		}
	}
	shift @cnt; #clear the missing count.
	#printf STDERR "%d %d\n", scalar(@cnt), sum(@cnt);
	#map {printf STDERR "\t%d\t%d\n", $_,$cnt[$_]} 0..$#cnt;
	return [scalar(@cnt), map($_/sum(@cnt), @cnt)];
}

exit 0;
