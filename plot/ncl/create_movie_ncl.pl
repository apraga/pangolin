use File::Find;

# Generates pictures from data set and create a movie by joining the pictures
my @dates;

# Add to array
sub add_movie { 
  # For .dat
  if ($_ =~ /dat/) {
    # Extract date from file name
    $_ =~ s/ratio_//;
    $_ =~ s/\.dat//;
    push @dates, $_;
  }
}

# Argument : date input output
sub replace {
  open my $in, '<', $_[1] or die "Can't open file: $!";
  open my $out, '>', $_[2] or die "Can't write file: $!";
  while (<$in>) {
    $_ =~ s/date =.*$/date = \"$_[0]\"/;
    print $out $_;
  }
  close($in);
  close($out);
}

#-------------------------------------------------------------------------------

my $dir = "/wkdir/pae2/praga/parallel_results/movie/";
my $tmp = "tmp.ncl";
my $file = "plot_contour.ncl";
my $output_dir = "files";

my $ncl = "/usr/local/ncl/bin/ncl";
mkdir $output_dir unless -e $output_dir;

# Generate date array
find(\&add_movie, $dir);

# Generate all plots
foreach(@dates) {
  my $date = $_;
  print "$date\n";

  # Replace in config
  replace($date, $file, $tmp);
 
  # Run ncl
  system("$ncl $tmp >/dev/null");
  system("convert -trim $output_dir/$date.png $output_dir/$date.png");
}

#Rename files to 001.dat 002.dat for ffmpeg
my $cur;
my $k = 0;
my @pics =  glob "$output_dir/*.png";
foreach(@pics) {
  $cur = sprintf("%08d", $k);
  rename $_, $output_dir."/".$cur.".png";
  $k++;
}

my $opts =  "-r 25 -sameq -s hd1080";
system("ffmpeg -i $output_dir/%08d.png $opts out.mp4");
