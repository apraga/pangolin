#my $exec = "ffmpeg -r 1  -pattern_type glob -i '*.png' -r 30 -pix_fmt yuv420p out.mp4";
my $ffmpeg = "/space/praga/softwares/ffmpeg/bin/ffmpeg";
my $exec = "$ffmpeg -r 25  -pattern_type glob -i 'files/*.png' out.mp4";
#my $exec = "$ffmpeg -r 1  -pattern_type glob -i '*.png' -crf 1 out.mp4";
system($exec);
