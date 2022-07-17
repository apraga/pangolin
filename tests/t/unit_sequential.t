#!/usr/bin/perl

# Sequential unit tests
# Change the number of first and last partitions here. 

use lib "t";
use UnitTests;
use Test::Simple v1.0; # This is found in deps
use Test::More;


my $n_min = 1;
my $n_max = 1197; # Default

# Compute the number of tests
if ($n_min == 1) { plan tests => $n_max/3+1;}
else { plan tests => ($n_max-$n_min)/3 + 1;}

my $test = new UnitTests($n_min, $n_max);
$test->set_type("seq");
$test->run();
