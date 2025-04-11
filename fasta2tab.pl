#!/usr/bin/env perl

# fasta2tab: Usage: perl fasta2tab.pl infile.fa > output_file
# Script from Johan Nylander: https://github.com/nylander/fasta-tab/blob/master/scripts/fasta2tab

local $/ = '>';
while(<>) {
    chomp;
    next if($_ eq '');
    my ($h, @S) = split /\n/;
    my $s = join('', @S);
    print STDOUT "$h\t$s\n" unless (!$h);
}