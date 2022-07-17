# Create a diff between two simulations on the Pangolin grid

# Arguments : file1 file2 output
# file 1 is the referenc
sub create_diff {
  open my $read1, $_[0]or die "Could not open $_[0]: $!";
  open my $read2, $_[1] or die "Could not open $_[1]: $!";

  open my $output, ">$_[2]" or die "Could not open $_[2]: $!";

  my $log = "diff.log";
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
    print $output $diff, " ", $values1[1], " " , $values1[2], "\n";
  } 
  close($read1);
  close($read2);
  close($output);
} 

#my $date = "201301030000";
my $date = "T";
#my $folder = "/wkdir/pae2/praga/parallel/cv_rate/";
#my $file1 = $folder."ratio_160lat_un_CFL0.7_0.dat";
#my $file2 = $folder."ratio_160lat_un_CFL0.7_$date.dat";
my $folder = "../../output_normal/";
my $file1 = $folder."ratio_1_201301010000.dat";
my $file2 = $folder."ratio_1_201301130000.dat";
create_diff($file1, $file2, "diff_normal.dat");

#my @dt = (9,7.5,6,5,4,3,0.6);
#foreach(@dt) {
#print "$_\n";
#my $file2 = $folder."ratio_dt_$_\_1_201301$date.dat";
#my $file1 = $folder."ratio_anal_201301$date.dat";
#(my $file2 = $file1)  =~ s/anal/1/;
#create_diff($file1, $file2, "diff_$_.dat");
#}
