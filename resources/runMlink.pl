#!/usr/bin/perl -w
use strict;
use File::Copy;

#prepare dir
my $workdir = shift;
my $resdir = shift;
my $bindir = "$resdir/bin";

open ERR, ">$workdir/err";

if (! -e $workdir) {
	exit;
}
print STDERR "start running $workdir\n";

#prepare genemap
my %genemap;
my %pos;
open MAP, "$resdir/genemap.txt" or die "cannot open $resdir/genemap.txt: $!\n";
while (<MAP>) {
	my ($chr, $start, $end, $gene) = (split)[0..3];
	my $checkey = join('_', $chr, $start, $end);
	$chr = 23 if $chr eq 'X';
	$chr = 24 if $chr eq 'Y';
	$chr = 25 if $chr eq 'XY';
	if (exists $pos{$checkey}) {
		$genemap{$gene} = ['out', $chr, $start, $end];
	} else {
		$pos{$checkey} = 1;
		$genemap{$gene} = ['in', $chr, $start, $end];
	}
}
close MAP or die "cannot close MAP: $!\n";

#generate lod scores along chromosome
chdir $workdir or die "cannot chdir to $workdir: $!\n";
my @markers = <./*.PRE>;
map s/^\.\/(\S+?)\.PRE$/$1/, @markers;
open RES, ">all_lodscores.txt" or die "cannot open all_lodscores.txt: $!\n";
for my $g (sort compPos @markers) {
#	if ($genemap{$g}->[0] eq 'out') {
#		next;
#	}
	copy("$g.LOC", 'datafile.dat');
	copy("$g.PRE", 'pedfile.pre');
	system("$bindir/makeped pedfile.pre pedfile.ped n");
	system("$bindir/pedcheck -p pedfile.ped -d datafile.dat -c >/dev/null 2>/dev/null");
	copy('zeroout.dat', 'pedfile.dat');
	system("$bindir/unknown");
	system("$bindir/mlink");
	copy('outfile.dat', 'final.out');
	system("/usr/bin/rm -f *.dat final.lod");
	open LL, "| $bindir/linklods" or die "cannot pipe linklods: $!\n";
	print LL "\n";
	close LL or die "cannot close pipe linklods: $!\n";
	#system("echo -e \"\n\" | $bindir/linklods");
	open LOD, "final.lod" or die "cannot open final.lod: $!\n";
	my @lods;
	while (<LOD>) {
		if (m/THETAS\s+(\d+\.\d+)/) {
			push @lods, [$1, 'NA'];
		} elsif (m/TOTALS\s+(-?\d+\.\d+)/) {
			$lods[-1]->[1] = $1;
		}
	}
	close LOD or die "cannot close LOD: $!\n";
	system("/usr/bin/rm -f ped* name*");
	for my $l (@lods) {
		my ($basename) = ($g =~ m/^(\S+?)(?:\[\d+\])?$/);
		printf RES "$g\t%d\t%d\t%d\t%f\t%f\n", (@{$genemap{$basename}})[1..3], @$l;
	}
}
close RES or die "cannot close RES: $!\n";
print STDERR "end running $workdir\n";


sub compPos {
	my ($g1, $suffix1) = ($a =~ m/^(\S+?)(?:\[(\d+)\])?$/);
	my ($g2, $suffix2) = ($b =~ m/^(\S+?)(?:\[(\d+)\])?$/);
	if (($g1 eq $g2) and $suffix1 and $suffix2) {
		return ($suffix1 <=> $suffix2);
	}
	if (not exists $genemap{$g1} or not exists $genemap{$g2}) {
		print ERR "$a $b $g1 $g2 $suffix1 $suffix2\n";
	}
	my($c1,$s1,$e1) = (@{$genemap{$g1}})[1..3];
	my($c2,$s2,$e2) = (@{$genemap{$g2}})[1..3];
	if ($c1 < $c2) {
		return -1;
	} elsif ($c1 == $c2 and $s1 < $s2) {
		return -1;
	} elsif ($c1 == $c2 and $s1 == $s2 and $e1 < $e2) {
		return -1;
	} elsif ($c1 == $c2 and $s1 == $s2 and $e1 == $e2) {
		return ($a cmp $b);
	}
	return 1;
}
exit 0;

