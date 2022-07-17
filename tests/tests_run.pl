#!/usr/bin/perl

#-------------------------------------------------------------------------------
# Script starting all the tests
# Default : no tests
# Use --unit or --func to enable unit or functional tests respectively (you can
# use both)
# Disable parallel unit tests with --no-parallel
#-------------------------------------------------------------------------------

use TAP::Harness;
use Getopt::Long; # Command line arguments

my $unit = 0;  
my $func = 0 ;
my $no_parallel = 0 ;
GetOptions ( 'unit' => \$unit, 'func' => \$func,
             'no-parallel' => \$no_parallel) or die "Error in arguments";

my @tests;
if (!$func && !$unit) { print "Use --func or --unit option\n"; exit;}

if ($unit) { 
  push(@tests, ['t/unit_sequential.t', 'Sequential unit tests']);
  if (!$no_parallel) {
    push(@tests, ['t/unit_parallel.t', 'Parallel unit tests']);
  }
}

if ($func) { push(@tests, 
    ['t/functional_hourdin.t', 'Functional test (Hourdin)'],
    ['t/functional_lauritzen.t', 'Functional test (Lauritzen)'],
    ['t/functional_io_hdf5.t', 'Functional tests (HD5 I/O)']
  );}

my $harness = TAP::Harness->new(
  #{verbosity=>1} # if you want more details
);

$harness->runtests(@tests);
