#!perl -w

use strict;

my $input = shift(@ARGV) or die;#Filtering/RALL_screened_by_distance
##RALL    TRINITY_DN45844_c0_g1_i1        201     3L      1       11262052        11262253        11262253.11262052
##I want to sort these by chromosomal location to make it easier to look at - also find duplicates / multiple transcripts

my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %ch = ();
my %data = ();

open(A, "<$input");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);

    my $pos = $a[4] . "." . $a[5];

    if(!(exists($data{$a[2] . "\t" . $a[4] . "\t" . $a[5]}))){
	$data{$a[2] . "\t" . $a[4] . "\t" . $a[5]} = $line;
	push(@{$ch{$a[2]}}, $pos);
    }elsif(exists($data{$a[2] . "\t" . $a[4] . "\t" . $a[5]})){

	print $line, "\n";
    }


}
close A;

open(B, ">>$output");
while((my $k, my $v) = each(%ch)){

    my @sort = sort {$a <=> $b} @{$ch{$k}};

    foreach my $s (@sort){
	my @m = split(/\./, $s);
	
	print B $data{$k . "\t" . $m[0] . "\t" . $m[1]}, "\n";

    }
}
close B;
