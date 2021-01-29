#!/usr/bin/perl

use warnings;
use strict;

my $NGROUPS = 5;
my $M = 10000;

# Group index file
open F, ">test.gri" or die $!;
for my $i (0..($M-1)) {
    printf F ("%d %d\n", $i, rand($NGROUPS));
}
close F;

# Group mixture file (single line)
open F, ">test.grm" or die $!;
for my $i (1..$NGROUPS) {
    printf F ("%s;", "0.0001,0.001,0.01");
}
close F;

