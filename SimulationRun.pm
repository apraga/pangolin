#!/usr/bin/perl

=pod
=head1 NAME

SimulationRun - Small module for sharing simulation variables (used by simulation_run.pl for
example)

=head1 SYNOPSIS

  use SimulationRun;

then use $batch and $mpirun. 

Note : these variables are updated by automake (during make install).
 
=cut

package SimulationRun;  
use strict;
use warnings;

# Export $batch and $mpirun
use Exporter;
our @ISA = 'Exporter';
our @EXPORT = qw($batch $mpirun $mpidir);

our ($batch, $mpirun, $mpidir);
$batch="";
$mpirun="mpirun";
$mpidir="";

1;
