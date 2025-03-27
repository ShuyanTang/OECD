use strict;
use warnings;
use Parallel::ForkManager;

my $pheno= shift @ARGV;
my $thread= shift @ARGV;

my $ForkManager=Parallel::ForkManager->new($thread);
foreach my $i (1..1000){
    $ForkManager->start and next;;
    system "Rscript permutation.R  $pheno $i";
    $ForkManager->finish;
}
$ForkManager->wait_all_children;