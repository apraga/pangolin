package CompareAsciiFiles;

=head1 NAME

CompareAsciiFiles - Perl extension for comparing two data files (ASCII)

=head1 SYNOPSIS

  use CompareFiles;

Compare two files and store errors in a log:

  compare_ascii_files("file1", "file2", $log);

=head1 DESCRIPTION

Compare concentration by using absolute error (less than 2e-12 by default but
can be overriden).
Read two ASCII-formatted data files and compare the values lines by lines.
Skips comments (#) and empty lines.
Data must be of the format :

  value x_coord y_coord

We only compure the values, not the coordinates

=cut

use strict;
use warnings;

# Export the main function
require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(compare_ascii_files);

my $max_error = 2e-12;

#-------------------------------------------------------------------------------
# Count the number of lines in files, ignoring comments and blank lines
#-------------------------------------------------------------------------------
sub nb_lines {
  open my $read, $_[0]or die "Could not open $_[0]: $!";
  my $nb = 0;
  while (my $line = <$read>) {
    if ($line !~ /#/ and $line !~ /^\s*$/) { $nb++;}
  }
  close $read;
  return $nb;
}

#-------------------------------------------------------------------------------
# Main function : compare two files and return the number of errors.
# Compare the first column of two files and check the difference is less than
# $max error. Also, appends the error to a log.
# Arguments : file1 file2 filelog [max_error]
#-------------------------------------------------------------------------------
sub compare_ascii_files {

  print "Comparing $_[0], $_[1]\n";
  print "Checking nb lines : ";
  my $n1 = nb_lines($_[0]);
  my $n2 = nb_lines($_[1]);
  if ($n1 != $n2) { print "number of lines differ : $n1 $n2 \n"; return;}
  print "ok \n";

  if (defined $_[3]) { $max_error = $_[3]; }

  print "Checking values : ";
  open my $read1, $_[0]or die "Could not open $_[0]: $!";
  open my $read2, $_[1] or die "Could not open $_[1]: $!";

  open my $log, ">$_[2]" or die "Could not open $_[2]: $!";
  #print $log "########### $_[0], $_[1]\n";

  my $diff_max = -9999.;
  my $line_max = 0;
  my $nb_errors = 0;

  while(my $line1 = <$read1>)  {   
    my $line2 = <$read2>;

    # Skip comment
    if ($line1 =~ /#/) { $line1 = <$read1>;}
    if ($line2 =~ /#/) { $line2 = <$read2>;}
    if ($line2 =~ /\*\*/) { die "File contains undefined values";}

    my @values1 = split(' ',$line1);
    my $q1 = $values1[0];

    my @values2 = split(' ',$line2);
    my $q2 = $values2[0];

    # Compare coordinates : latitude
    my $diff = abs($values1[1] - $values2[1]);
    if ($diff > 2e-13) { 
      print $log "$. : lat $diff\n" ; 
      $nb_errors = $nb_errors + 1;
    }

    # Compare coordinates : longitude
    $diff = abs($values1[2] - $values2[2]);
    if ($diff > 2e-13) { 
      print $log "$. : lon $diff\n" ; 
      $nb_errors = $nb_errors + 1;
    }

    # Absolute error
    $diff = $q2 - $q1;

    # Compare concentration by using relative error, except when data is small
#    if (abs($q2 - $q1) < 1e-13){ next; }
#    if (abs($q2) > abs($q1)) {
#      $diff = ($q2 - $q1)/$q2;
#    }
#    else {
#      $diff = ($q2 - $q1)/$q1;
#    }
    # Scientific notation
    $diff = sprintf("%e", $diff);

    if (abs($diff) > $max_error) {
      print $log "$. : $diff $values1[1] $values1[2]\n" ;
      $nb_errors = $nb_errors + 1;

      if (abs($diff) > $diff_max) {
        $diff_max = $diff;
        $line_max = $.;
      }
    }
    #last if $. == 10;
  }

  print "$nb_errors errors\n";
  if ($nb_errors > 0) { warn "Max diff $diff_max at line $line_max\n"};
  close $read1;
  close $read2;
  
  close $log;
  return $nb_errors;
}

1;
