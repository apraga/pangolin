=pod

=head1 NAME

Functional - a class for Pangolin functional tests.

=head1 SYNOPSIS
 
  use Functional;
  my $test = new Functional("hourdin");
  $test->run();

=head1 DESCRIPTION

This is a main class defining functional test : it will run a defined test case
for a varying number of cores, which must be 1 or a multiple of 3.
Default is 24.
The configuration is defined by the constructor. See the examples for more
informations.

It compares parallel results to sequential.
Parallel simulation are run for 1 and then rom n_min to n_max, except when output files already
exists. So you can skip a run by creating yourself output data.

A run is blocking, so it may take a long time.
For parallel run, you must be able to start direct MPI jobs (node access) or put
the perl script in a job.
Needs initial conditions : ratio and constant winds in input/ folder
It is your responsibility to ensure these file exists : 
ratio_201301010000.h5 and for Hourdin u_201301010000.dat, v_201301010000.dat

By default, the data format is hdf5. This can be changed with set_io().

Warning: for parallel I/O, ensure your filesystem is configured for it !
On Neptune (CERFACS), this means you have to read and write from/to the
/scratch.

All the tests run with a fixed configuration (80lat) but that can be changed
with set_nb_lat2(80).

=head1 EXAMPLES

  use Functional;
  my $test = new Functional("hourdin");
  $test->set_nmax(9);
  $test->set_nbtracers(9);
  $test->run();

Change the I/O format to ascii :

  use Functional;
  my $test = new Functional("hourdin");
  $test->set_io("ascii");
 
=cut

package Functional;
use CompareAsciiFiles; # Comparing files
use Test::Simple v1.0; # This is found in deps
use Test::More;
use lib "../";
use SimulationRun; 

# Argument : "hourdin", "lauritzen"
sub new
{
    my $class = shift;
    my $self = {
        type => shift, 
        mpirun => "$mpidir/mpirun",# For clusters, must have direct node access
        n_min => 3,
        n_max => 24,
        nb_traces => 1,
        nb_lat2 => 80,
        test_case => undef,
        input_dir => undef,
        input_winds_dir => undef,
        io => "hdf5", #default but we can have ascii too
        ext => ".h5", # file extension, set by io type
    };
    bless $self, $class;

    if ($self->{type} eq "hourdin") { 
      $self->hourdin(); 
    }
    elsif ($self->{type} eq "lauritzen") { 
      $self->lauritzen(); 
    }
    elsif ($self->{type} eq "hdf5_io") { 
      $self->hdf5_io(); 
    }
    else {
      die "Wrong type: must be hourdin, lauritzen or hdf5_io";
    }
    return $self;
}

# Hourdin configuration : 2 days of simulation
# Fixed winds
sub hourdin {
  my $self = shift;
  print "Hourdin case\n";
  $self->{t_start} = "201301010000";
  $self->{t_end} = "201301030000";
  $self->{dt} = "3";
  $self->{output_folder} = "output_hourdin";
}

# Lauritzen configuration : 6 days of simulation
# Time-dependant winds
sub lauritzen {
  my $self = shift;
  print "Lauritzen case\n";
  $self->{t_start} = "201301010000";
  $self->{t_end} = "201301070000";
  $self->{dt} = "15";
  $self->{test_case} = "testsuite_cv";
  $self->{output_folder} = "output_lauritzen";
}

# Testing hdf5 I/O
sub hdf5_io {
  my $self = shift;
  print "hdf5 I/O\n";
  $self->{t_start} = "201301010000";
  $self->{t_end} = "201301010000";
  $self->{dt} = "1";
  $self->{output_folder} = "output_hdf5";
}

sub set_io {
  my $self = shift;
  $self->{io} = shift;

  # Ascii format
  $self->{ext} = ".dat" if !( $self->{io} eq "hdf5") ;
}

sub set_nmin {
  my $self = shift;
  $self->{n_min} = shift;
}

sub set_nmax {
  my $self = shift;
  $self->{n_max} = shift;
}

sub set_nbtracers {
  my $self = shift;
  $self->{nb_tracers} = shift;
  $self->check_nbtracers();
}

sub set_nb_lat2 {
  my $self = shift;
  $self->{nb_lat2} = shift;
}

# Set ratio and winds directory
sub set_input_folders {
  my $self = shift;
  $self->{input_dir} = shift;
  $self->{input_winds_dir} = shift;
}

sub set_output_folder {
  my $self = shift;
  $self->{output_folder} = shift;
}

# Ensure the nb of tracers is the same in the code
sub check_nbtracers {
  my $self = shift;
  open my $file, '<', '../src/parameters.F90'   or die "Can't open: $!";

  while (<$file>) {
    while (/NB_TRACERS = /g) {
      $_ =~ s/\D//g;
      die "Nb of tracers different from the code (should be $_)" if ($_ != $self->{nb_tracers});
    }
  }
}

# Generate a config file for a complete run 
# Arguments : config_file nb_cores
sub generate_config {
  my $self = shift;
  open my $config, ">$_[0]" or die "Could not open $_[0]: $!";

  print $config "nb_partitions = $_[1]\n";
  print $config "nb_lat2 = $self->{nb_lat2}\n";
  print $config "t_start = $self->{t_start}\n";
  print $config "t_end = $self->{t_end}\n";
  print $config "dt = $self->{dt}\n";
  print $config "T_winds = 0000\n";
  print $config "T_output = 130000\n";

  print $config "input_dir = $self->{input_dir}\n";
  print $config "input_winds_dir = $self->{input_winds_dir}\n";

  print $config "output_dir = $self->{output_folder}/$_[1]cores\n";
  print $config "output_winds_dir = $self->{output_folder}/$_[1]cores\n";
  if (defined $self->{test_case}) {
    print $config "test_case = $self->{test_case}\n";
  }

  close $config;
}

# Run sequential and parallel to compare them
sub run_parallel {
  my $self = shift;
  my $config = $_[0];

  $self->run_simulation(1, $config);
  my $n_min = $self->{n_min};
  my $n_max = $self->{n_max};

  for (my $n = $n_min; $n <= $n_max; $n=$n+3) {
    subtest "Nb partitions = $n" => sub {
      $self->run_simulation($n, $config);
      $self->compare_ratio_seq($n);
    }
  }
}

# Fake run to check the output is the same as the input
sub run_fake {
  my $self = shift;
  my $config = $_[0];

  my $n_min = $self->{n_min};
  my $n_max = $self->{n_max};

  my $n = $n_min;
  while ($n <= $n_max) {
    subtest "Nb partitions = $n" => sub {
      $self->run_simulation($n, $config);
      $self->compare_to_input($n);
    };
    if ($n == 1) { $n = 3;}
    else { $n += 3; }
  }
}

# Starts all run 
sub run {
  my $self = shift;
  my $config = "config.tmp";
  $self->check_input();
  mkdir $self->{output_folder} unless -e $self->{output_folder};

  # Run either fake or real tests
  if ($self->{type} eq "hdf5_io") { $self->run_fake($config); }
  else { $self->run_parallel($config);}

}

# Start simulation if output files does not exist
# Arguments : nb_cores config_file
sub run_simulation {
  my $self = shift;
  print "Run for $_[0] cores... ";
  $self->generate_config($_[1], $_[0]) ;

  my $folder = "$self->{output_folder}/".$_[0]."cores";
  mkdir $folder unless -e $folder;

  if (-e $folder."/ratio_1_$self->{t_end}$self->{ext}") {
    print "testing folder", $folder, "\n";
    print "Skip run, file exists\n";
    return;
  }
  print "\n";
  system("$self->{mpirun} -np $_[0] ../bin/pangolin --config=$_[1] > /dev/null");
  #system("$self->{mpirun} -np $_[0] ../src/pangolin --config=$_[1]");

  # Clean
  unlink <$folder/profiling*>;
}

sub check_file { die "Missing input data: $_[0]" unless -e $_[0]; }
sub check_folder { die "Missing folder: $_[0]" unless -d $_[0]; }

# Check input directory and files exist
sub check_input {
  my $self = shift;

  my $folder = $self->{input_dir};
  check_folder($folder);
  my $ext = $self->{ext};

  for (my $i=1;$i <= $self->{nb_tracers}; ++$i){
    check_file($folder."/ratio_$i\_201301010000".$ext);
  }

  $folder = $self->{input_winds_dir};
  check_folder($folder);
  if (!defined($self->{test_case})) {
      check_file($folder."/u_201301010000".$ext);
      check_file($folder."/v_201301010000".$ext);
    }
}

# Start our executable and check return code (0 if success)
# Arguments : file1 file2 logfile
sub compare_hdf5_files {
  my $exec = "../bin/compare_hdf5";
  die "${exec} not found or not executable" unless (-e $exec && -x $exec);

  system("${exec} ${_[0]} ${_[1]} > ${_[2]}");
  return $?;
}

# Check the simulation write is identical to input
sub compare_to_input {
  my $self = shift;
  my $nb = $_[0];
  my $id = 1; # Only compare first tracer

  subtest 'All data '.$self->{test_case} => sub {
    foreach $file ("ratio_${id}", "u", "v") {
      SKIP: {
        my $file1 = "$self->{input_dir}${file}_$self->{t_end}";
        if ($file eq "u" || $file eq "v") {
          $file1 = "$self->{input_winds_dir}${file}_$self->{t_end}";
        }
        my $file2 = "$self->{output_folder}/".$nb."cores/${file}_$self->{t_end}";
        my $log = "$self->{output_folder}/diff_${nb}_${id}.log";
        my $test = "Nb cores = $nb";

        $file1 = $file1.$self->{ext};
        $file2 = $file2.$self->{ext};
        print "Comparing ", $file1, " ", $file2, "\n";

        skip "Missing input ${file} file", 1, if (!-e $file1);
        skip "Missing output ${file} file", 1, if (!-e $file2);
        cmp_ok(compare_hdf5_files($file1, $file2, $log), "==", 0, $test) ;
      }
    }
  }
}

# Compare ratio files
# Argument : nb_cores, tracer_id
sub compare_ratio_seq_files {
  my $self = shift;
  my $nb = $_[0];
  my $id = $_[1];
  my $file1 = "$self->{output_folder}/1cores/ratio_".$id."_$self->{t_end}";
  my $file2 = "$self->{output_folder}/".$nb."cores/ratio_$id\_$self->{t_end}";
  my $log = "$self->{output_folder}/diff_".$nb."_$id.log";
  my $test = "Nb cores = $nb";

  $file1 = $file1.$self->{ext};
  $file2 = $file2.$self->{ext};

  if ($self->{io} eq "hdf5") {
    cmp_ok(compare_hdf5_files($file1, $file2, $log), "==", 0, $test) ;
  }
  else {
    cmp_ok(compare_ascii_files($file1, $file2, $log, $max_err), "==", 0, $test) ;
  }
}

# Compare parallel output to sequential output
# Arguments : nb_cores
sub compare_ratio_seq {
  my $self = shift;
  my $nb = $_[0];
  my $max_err = 5e-11;
 
  # One subtest for each tracer
  subtest 'All tracers '.$self->{test_case} => sub {
    plan tests => $self->{nb_tracers};
    foreach my $id (1..$self->{nb_tracers}) {
      $self->compare_ratio_seq_files($nb, $id);
    }
  };
}

1;
