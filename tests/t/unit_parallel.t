#!/usr/bin/perl

# Parallel unit tests
# Change the number of first and last partitions here. No arguments is the
# default.

use lib "t";
use UnitTests;
use Test::Simple v1.0; # This is found in deps
use Test::More;

my $n_min = 66; 
my $n_max = 126; 

# Compute the number of tests
plan tests => ($n_max-$n_min)/3 + 1;

my $test = new UnitTests($n_min, $n_max);
$test->set_type("paral");
$test->run();
