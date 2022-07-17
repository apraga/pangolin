#-------------------------------------------------------------------------------
# Functional testing with to check HDF5 I/O. It will run tests for different
# nb_lat2. 80 lat2 is mandatory and it will try to launch tests for 160 and
# 320lat.
# Used by tests_run.pl.
# Can be run in standalone with 
# [PERL5LIB=deps/Test-Simple/] perl t/functional_hourdin.t 
#-------------------------------------------------------------------------------

use lib "t";
use Functional;
use Test::Simple v1.0; # This is found in deps
use Test::More;

my $n_min = 1;
my $n_max = 126;

# !! For parallel I/O you need to do it on the /scratch on Neptune !!
my $folder = "/scratch/ae/praga/";
#my $folder = "/wkdir/pae2/praga/";

# Define the common test
$test = new Functional("hdf5_io");
$test->set_nmin($n_min);
$test->set_nmax($n_max);
$test->set_nbtracers(1);

#my $ratio_in = $folder."input/gaussianhills_80lat/";
my $ratio_in = $folder."input/random/80lat/";
my $winds_in = $folder."input/cv_winds/80lat/";
my $folder_out = $folder."tests/output_hdf5";
subtest "80lat" => sub {
  $test->set_input_folders($ratio_in, $winds_in);
  $test->set_output_folder($folder_out);
  $test->set_nb_lat2(80);
  $test->run();
};


$ratio_in =~ s/80/160/;
$winds_in =~ s/80/160/;
$folder_out = $folder."tests/output_hdf5_160";
subtest "160lat" => sub {
  plan skip_all => "Missing directory ${ratio_in}" if !(-d $ratio_in);
  plan skip_all => "Missing directory ${winds_in}" if !(-d $winds_in);
  $test->set_nb_lat2(160);
  $test->set_input_folders($ratio_in, $winds_in);
  $test->set_output_folder($folder_out);
  $test->run();
};

$ratio_in =~ s/160/320/;
$winds_in =~ s/160/320/;
$folder_out =~ s/160/320/;

subtest "320lat" => sub {
  plan skip_all => "Missing directory ${ratio_in}" if !(-d $ratio_in);
  plan skip_all => "Missing directory ${winds_in}" if !(-d $winds_in);
  $test->set_nb_lat2(320);
  $test->set_input_folders($ratio_in, $winds_in);
  $test->set_output_folder($folder_out);
  $test->run();
};

done_testing();
