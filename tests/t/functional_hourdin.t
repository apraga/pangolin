#!/usr/bin/perl

#-------------------------------------------------------------------------------
# Functional testing with Hourdin test case.
# Used by tests_run.pl.
# Can be run in standalone with 
# [PERL5LIB=deps/Test-Simple/] perl t/functional_hourdin.t 
#-------------------------------------------------------------------------------

use lib "t";
use Functional;
use Test::Simple v1.0; # This is found in deps
use Test::More;

my $n_min = 3;
my $n_max = 126;
plan tests => (($n_max-$n_min)/3+1);

# !! For parallel I/O you need to do it on the /scratch on Neptune !!
my $folder = "/scratch/ae/praga/";
my $folder_in = $folder."input/hourdin_80lat/";
my $folder_out = $folder."tests/output_hourdin/";

$test = new Functional("hourdin");
#$test->set_io("ascii");
$test->set_nmin($n_min);
$test->set_nmax($n_max);
$test->set_input_folders($folder_in, $folder_in);
$test->set_output_folder($folder_out);
$test->set_nbtracers(1);
$test->run();
