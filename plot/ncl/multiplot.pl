#Replace input files in contour plot for multiple plots
use strict;
use warnings;

use File::Slurp qw(read_file write_file);

my $ncl = "/usr/local/ncl/bin/ncl";
my $in = 'plot_contour.ncl';
my $out = 'tmp.ncl';
my $data = read_file $in, {binmode => ':utf8'};

# Generate config for different nb of procs
for (my $i = 3; $i <= 126; $i+=3) {
  print "i=$i\n";
  my $tmp = $data;
  $tmp =~ s/[0-9]+(?=cores)/$i/g;
  $tmp =~ s/nb_partitions = \K.*[0-9]*/$i/g;
  write_file $out, {binmode => ':utf8'}, $tmp;
  system("$ncl $out");
}
