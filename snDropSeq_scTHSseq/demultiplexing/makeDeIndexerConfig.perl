#!/usr/bin/perl -w

use strict;


my %combo;
my $i7_table = $ARGV[0];
my $i5_table = $ARGV[1];
my $r5_table = $ARGV[2];

my @i7;
my @i5;
my @r5;
open(IN, $i7_table);
while(my $line = <IN>){
	chomp($line);
	push(@i7, $line);
}
close(IN);

open(IN, $i5_table);
while(my $line = <IN>){
	chomp($line);
	push(@i5, $line);
}
close(IN);

open(IN, $r5_table);
while(my $line = <IN>){
	chomp($line);
	push(@r5, $line);
}
close(IN);

for(my $i = 0; $i <scalar(@i7); $i++){
	for(my $j=0; $j < scalar(@i5); $j++){
		for(my $k=0; $k < scalar(@r5); $k++){
			print $i+1, "\t", $j+1, "\t", $k+1, "\t", $i7[$i], "_", $i5[$j], "_", $r5[$k], "\n";
		}
	}
}
