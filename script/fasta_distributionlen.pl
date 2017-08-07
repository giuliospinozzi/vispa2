#!/usr/bin/perl

#
# distribution length of a fasta file
#

use warnings;
use strict;
use Bio::SeqIO;

my $i = $ARGV[0];
my $o = $ARGV[1];

open (OUT, '>',$o);

my $seq_in  = Bio::SeqIO->new( -format => 'fasta',-file => $i);

my %length = ();
while( my $seq = $seq_in->next_seq() ) {
	my $seq1 = $seq->seq;	
	chomp($seq1);
	#print($seq1);
	my $l = length($seq1);
	#print $l."\n";
	my $tmp =  $length {$l};
	$length {$l} = $tmp+1;  
}
my @x = keys(%length);
my $sum = 0;
foreach my $k (@x) {
   $sum = $sum + $length{$k};
}

@x = sort{$a <=> $b}@x;
print(@x);
foreach my $k (@x) {
   print OUT $k."\t".$length{$k}."\t".$length{$k}*100/$sum."\n";
}

close (OUT);
