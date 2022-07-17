use PDF::API2;

die "Needs 3 arguments\n" if ($#ARGV+1 != 2);
my $pdf = PDF::API2->open(@ARGV[0]);
my $page = $pdf->openpage(1);
$page->rotate(-90);
$pdf->saveas(@ARGV[1]);
