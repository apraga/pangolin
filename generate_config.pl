# Generate a list of configuration files in tmp/ from the config
use strict;
use warnings;

use File::Slurp qw(read_file write_file);

my $in = 'config';
my $data = read_file $in, {binmode => ':utf8'};

# Generate config for different nb of procs
sub generate_nb_procs {
  for (my $i = 3; $i < 64; $i+=3) {
    my $out = 'tmp/config_'.$i;
    my $tmp = $data;
    $tmp =~ s/[0-9]+(?=cores)/$i/g;
    $tmp =~ s/nb_partitions = \K.*[0-9]*/$i/g;
    write_file $out, {binmode => ':utf8'}, $tmp;
  }
}


# Generate config for different nb of lat and dt
sub generate_nb_lat_dt {
  my @nblat2 = (80, 160, 320);
#  my @dt = (60, 30, 15, 7.5, 4, 2);
  my @dt = (20, 10, 5);
  foreach my $i (0..@nblat2-1) {
    my $out = 'tmp/config_'.$nblat2[$i];
    my $tmp = $data;
    $tmp =~ s/dt=\K[0-9]*/$dt[$i]/g;
    $tmp =~ s/nb_lat2 = \K.*[0-9]*/$nblat2[$i]/g;
    $tmp =~ s/gaussianhills_\K[0-9]*/$nblat2[$i]/g;
    $tmp =~ s/output(?=\n)/tmp_$nblat2[$i]/g;
    write_file $out, {binmode => ':utf8'}, $tmp;
  }
}

# Generate config for different dt
sub generate_dt {
  my @dt = (10, 9, 8, 6, 5, 4);
  foreach(@dt){
    my $out = 'tmp/config_'.$_;
    my $tmp = $data;
    $tmp =~ s/dt=\K[0-9]*/$_/g;
    $tmp =~ s/dt_\K[0-9]*/$_/g;
    write_file $out, {binmode => ':utf8'}, $tmp;
  }
}

#-------------------------------------------------------------------------------

#generate_dt();
generate_nb_procs();
