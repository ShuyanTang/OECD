use strict;
use warnings;
use feature qw(say);


my $input_file= 'OE.case_control.all.snv';
my $outputfile = "OE.case_control.LoF_A.snv";


open FILE,"$input_file" or die "$input_file\t$!";
open SAVE,">$outputfile";
my $head=<FILE>;chomp($head);
print SAVE "$head\n";
my $i=0;
my %heads;
foreach my $h(split /\t/,$head){
    $heads{$h}=$i++;
}
while(my $content = <FILE>){
    $content =~ s/[\r\n]//g;
    my @a = split /\t/,$content;
    my $function = $a[$heads{'Function'}];
    my $length = $a[$heads{'Protein Length'}];
    my $aa = $a[$heads{'AA_Change_Pos'}];

    my $chr=$a[$heads{'Chrs'}];
    my $start=$a[$heads{'Start'}];
    my $end=$a[$heads{'End'}];
    my $ref=$a[$heads{'Ref'}];
    my $alt=$a[$heads{'Alt'}];
    my $snv = join "\t",($chr,$start,$end,$ref,$alt);

    ## remove freq>0.001
    my $freq=0;
    foreach my $k(split /,/,"gnomAD_t_EAS_AF,gnomAD_total_AF"){ 
        $freq++ if ($a[$heads{$k}] !~ /^[\s\.]*$/ and $a[$heads{$k}] >= 0.001 ); 
    }
    next if $freq>0;

    ## select varaints
    my $mut_type=0;
    $mut_type++ if( $function =~  /stop_gain|frameshift|splice site/ and $aa/$length <0.98 ); 
    $mut_type++ if ($function =~ /missense/ and $a[$heads{'AlphaMissense_Class'}] =~ /pathogenic/);
    next if $mut_type==0;

    print SAVE "$content\n";

}
close FILE;
close SAVE;