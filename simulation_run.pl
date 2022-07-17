#! /usr/bin/perl
#-------------------------------------------------------------------------------
# Script for starting the simulation. Without arguments, will ask for the first
# and last number of partitions. You can override it by passing it as arguments.
# Usage : perl simulation_run.pl [--start=NB --end=NB2] [--startend==NB]
#-------------------------------------------------------------------------------


use warnings;
use strict;
use SimulationRun; # Contains $bash and $mpirun variables
use Getopt::Long; # Command-line arguments
use File::Copy; 

my $n_start = '-1';  
my $n_end = '-1';  
my $no_run = '-1';  
my $config = 'config';  
my $logdir = 'log';


#-------------------------------------------------------------------------------
# Functions definitions
#-------------------------------------------------------------------------------

# Check the nb of partitions are sorted
sub get_correct_nb_run() {
  my ($n_start, $n_end) = get_nb_run();
  while ($n_start > $n_end) {
    print "n_start should be less than n_end\n";
    ($n_start, $n_end) = get_nb_run();
  }
  return $n_start, $n_end;
}


# Get first and last number of partitions
sub get_nb_run() {
  print "Number of partitions (start): \n";
  chomp(my $n_start = <>);
  while (!check_nb_parts($n_start)) { chomp($n_start = <>); }

  print "Number of partitions (end): \n";
  chomp(my $n_end = <>);
  while (!check_nb_parts($n_end)) { chomp($n_end = <>); }
  return ($n_start, $n_end);
}

# Check number of partitions. Return 0 if wrong
sub check_nb_parts {
  while ($_[0] < 1) {
    print "Value should be > 0\n"; 
    return 0;
  }
  while ($_[0] > 1 && $_[0] %3 != 0) {
    print "Value should be multiple of 3 or 1\n"; 
    return 0;
  }
  return 1;
}

# Start a mpi run directly
# Arguments : n_start, n_end, config, no-run, logdir
sub run_simu  {
  print "Simulation \n";
  my $conf_new = $_[2].".new";
  my $flags = "--config=$conf_new";
  if ($_[3] > 0) { $flags .= " --no-run"; }

  for (my $n = $_[0]; $n <= $_[1]; $n=$n+3) {
    # Log file : append perl PID number to have a unique log name
    my $log = $_[4]."/log_$n.$$";

    print "Nb partitions = $n\n";
    substitute($_[2], $conf_new, "nb_partitions.*","nb_partitions = $n", );

    # Allow log file to be updated 
    #system("$mpirun -np $n src/pangolin $flags &> $log ");
    system("$mpirun -np $n bin/pangolin $flags");
  }
}

# Arguments : n_start, n_end, config
sub generate_config {
  print "Generating config files\n";
  my $folder = "tmp";
  mkdir $folder;
  for (my $n = $_[0]; $n <= $_[1]; $n=$n+3) {
    my $conf_new = "$folder/config_$n";
    my @args = ("nb_partitions.*","nb_partitions = $n", 
      "ratio_paral.*","ratio_paral_$n",
      "output_u.*","output_u = u_$n",
      "output_v.*","output_v = v_$n");
    substitute($_[2], $conf_new, @args);
  }
}

# Start a run with a batch job
# Arguments : n_start, n_end, config, no-run
sub run_batch {

  # Log file is set inside the batch job
  my $batch_new = $batch.".new";
  my @args = ("nb_min=.*","nb_min=$_[0]",
    "nb_max=.*","nb_max=$_[1]", 
    "select=.*","select=$_[1]");
  substitute($batch, $batch_new, @args);

  system("$mpirun $batch_new");
}

# Substitute a list of strings in a file. Create a backup which will used for
# the simulation. 
# Arguments : file, new_file, [string_in, string_out, ...]
sub substitute {
  open my $in,  '<',  $_[0] or die "Can't read $_[0]: $!";
  my $fnew = defined($_[1])? "$_[1]" : "$_[0].new" ;
  open my $out, '>', $fnew or die "Can't write $fnew: $!";

  while( <$in> ) {
    # Replace several strings
    for (my $i = 2; $i < $#_; $i+=2) {
      s/$_[$i]/$_[$i+1]/g;
    }
    print $out $_;
  }

  close $in;
  close $out;
}


#-------------------------------------------------------------------------------
# Main 
#-------------------------------------------------------------------------------
GetOptions ('start=i' => \$n_start, 
'end=i' => \$n_end,
          'startend=i' => sub {$n_start=$_[1]; $n_end=$n_start},
          'no-run' => \$no_run,
        'config=s' => \$config,
      'logdir=s' => \$logdir);

  if ($n_start < 0 || $n_end < 0) {
    print "Enter the number of processus: \n";
    ($n_start, $n_end) = get_correct_nb_run();
  }
  else {
    die if (!check_nb_parts($n_start));
    die if (!check_nb_parts($n_end));
  }

  if ($batch eq "") { 
    run_simu($n_start, $n_end, $config, $no_run, $logdir);
  }
  else {
    generate_config($n_start, $n_end, $config);
    run_batch($n_start, $n_end, $config, $no_run);
  }
