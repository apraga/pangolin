#!/usr/bin/perl

=head1 NAME

UnitTests - module for running sequential and parallel tests

=head1 SYNOPSIS

  use UnitTests;

Usage : set the first and last number of partitions and sequential/parallel :

  my $test = new UnitTests($n_min, $n_max);
  $test->set_type("paral");
  $test->run();

Warning : the number of partititions must be 1 or a multiple of 3.

=cut

package UnitTests;  

use strict;
use warnings;

# Export main subroutine
use Exporter;
our @ISA = 'Exporter';
our @EXPORT = qw(run_all_tests generate_config);

use lib '../'; # Simulation run variables, relative to current folder !
use SimulationRun;
use Test::Simple v1.0; # This is found in deps
use Test::More;

# Constructor, Argument : n_min, m_max
sub new
{
    my $class = shift;
    my $self = {
        n_min => shift,
        n_max => shift,
        mpirun => "$mpidir/mpirun",# For clusters, must have direct node access
        nb_traces => 1,
        type => undef,
        folder => undef,
    };
    bless $self, $class;
    check_nb_parts($self->{n_min});
    check_nb_parts($self->{n_max});

    return $self;
}


# Check nb of partitions
sub check_nb_parts {
  die "$_[0] should be 1 or multiple of 3\n" if ($_[0] != 1 && $_[0] % 3 != 0);
}

# 'Seq' or 'paral'
sub set_type {
  my $self = shift;
  my $tmp = shift;
  if ($tmp ne "seq" && $tmp ne "paral") {
    die "Wrong first argument : must be seq or paral";
  }
  if ($tmp eq "paral" && $self->{n_min} < 3) {
    die "Parallel test must start at 3 partitions";
  }
  $self->{type} = $tmp;
  $self->set_exe();
}

sub set_nb_procs {
  my $self = shift;
  $self->{nb_procs} = 1;
  if ($self->{type} eq "paral") {$self->{nb_procs} = $_[0];}
}

sub set_exe {
  my $self = shift;

  my $folder = "src/";
  $self->{exe} = $folder."tests_sequential";
  if ($self->{type} eq "paral") {$self->{exe} = $folder."tests_parallel";}
}

# Write a config file 
# Arguments : config_file nb_partitions
sub generate_config {
  die "Needs 2 arguments" if ($#_ != 1) ;

  open my $file,  '+>',  $_[0]      or die "Can't read $_[0]: $!";
  print $file "nb_partitions = $_[1]\n";
  print $file "nb_lat2 = 90\n";
  print $file "test_case = zonal\n";
  print $file "output_dir = output\n";
  close $file;
}


# Run tests for a given number of partition and check the output.
# Arguments : nb_partitions
sub run_test {
  my $self = shift;
  my $conf="config.tmp";
  mkdir "output" unless -d "output";
  generate_config($conf, $_[0]);

  my $args = "-np $self->{nb_procs} $self->{exe} --config=$conf --no-run";
  # Pipe the output of the run and check for error
  # We also parse possible errors
  open PS, "$self->{mpirun} $args 2>&1 |";
  while (<PS>) {
    my $cur = $_;
    if ($cur =~ /error/){ BAIL_OUT($cur);}

    # If the line begins by "Checking", then we search for the string "failed"
    if ($cur =~ /^\s*Checking/){
      chomp($cur);
      # Name of the current check
      my @name = split(/\.\.\./, $cur);
      my $test_name = $name[0]." ($_[0] partition(s))" ;
      my $failed = $cur =~ /fails/;
      isnt($failed, 1, $test_name );
    }
  }
  close PS;
}

sub check_type {
  my $self = shift;
  if (!defined $self->{type}) {
    die "You need to specify the type : parall or seq";
  }
}

# Run all tests
sub run {
  my $self = shift;

  $self->check_type();

  print "Running tests from ".$self->{n_min}." to ".$self->{n_max}."\n";

  for (my $n = $self->{n_min}; $n <= $self->{n_max}; $n = $n+3){
    $self->set_nb_procs($n);
  
    subtest "Nb partitions = $n" => sub { 
      $self->run_test($n); 
    };
    # We increment by 3, so special case for 1 partition
    if ($n == 1) { $n = 0;}
  }
}

1;
