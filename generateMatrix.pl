use strict;
use warnings;

my $input= "OE.case_control.LoF_A.snv"; 
my $output="OE.case_control.LoF_A.matrix"; 


my (%final,%matrix,%sample)=();


open FILE,"$input" or die "$input\t$!";
my $head=<FILE>;
chomp($head);
my $i=0;
my %heads;
foreach my $h(split /\t/,$head){
    $heads{$h}=$i++;
}
while(my $content = <FILE>){
    $content =~ s/[\r\n]//g;
    my @a = split /\t/,$content;
    my $id=$a[$heads{'Gene ID'}];
    foreach my $cc ('Case_depth','Control_depth'){
        foreach  my $geno( split ";",$a[$heads{$cc}]){
            next if ($geno =~ /^[\s\.]$/);
            my ($sam,$gt,$ref_depth,$alt_depth)=split ",",$geno;
            my $allele=0;
            if ($gt eq 'HET'){
                $allele=1;
            }elsif($gt eq 'HOMA'){
                $allele=2;
            }
            $sample{$sam}=1;
                $final{$id}{$sam}+=$allele;
        }
    }
    print SAVE "$content\n";
}
close FILE;
close SAVE;


## write matrix
my @samsort;
foreach my $sam (sort keys %sample){
    push @samsort,$sam;
}

my $allsamples=join "\t",@samsort;

open MATRIX,">$output";
print MATRIX "ID\t$allsamples\n";
foreach my $id (sort keys %final){
    print MATRIX "$id";
    foreach my $case (@samsort){
        my $n=$final{$id}{$case}; 
        print MATRIX "\t$n";
    }
    print MATRIX "\n";
}
close MATRIX;






