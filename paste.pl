# Extract ratio data (first columns) from 2 files for scatter plot
use strict;
use warnings;

use File::Slurp qw(read_file write_file);

die "Needs 2 arguments" if ($#ARGV != 1)  ;

my $cfl='0.95';
my $folder="/wkdir/pae2/praga/parallel/mixing";
my $file1="$folder/ratio_$ARGV[0]lat_$ARGV[1]_cos_CFL${cfl}_T2.dat";
my $file2=$file1;
$file2 =~ s/cos/corr/;
my $file_out = "$folder/mix-pangolin-$ARGV[0]lat-$ARGV[1]-CFL${cfl}.dat";

open my $data1, '<', $file1;
open my $data2, '<', $file2;
open my $out, '>', $file_out;
while (!eof($data1) and !eof($data2)) {
  my @q1 = split(/\s+/, <$data1>);
  my @q2 = split(/\s+/, <$data2>);
  print $out $q1[1], " ", $q2[1], "\n";
}
