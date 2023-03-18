#!perl -w

use strict;

my $input = shift(@ARGV) or die; #Sorted_AG_candidates
my $dist = shift(@ARGV) or die; #5kb using
my $fasta = shift(@ARGV) or die;#the concat fasta file for dmel
my $output = shift(@ARGV) or die;
unlink(qq{$output});
my $output2 = shift(@ARGV) or die;
unlink(qq{$output2}); #the list of expanded positions 

my %id = ();
open(A, "<$input");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    if($line =~ m/TRINITY/){
	if($a[4] < $a[3]){
	    print "switched\n";
	}
	if(!(exists($id{$a[1]}))){
	    $id{$a[1]} = $a[2] . "\t" . $a[3] . "\t" . $a[4];
	}
    }
}
close A;

my %ch = ();
my $tmp = 0;
open(B, "<$fasta");
while(my $line2 = <B>){
    chomp $line2;
    if($line2 =~ m/^>/){
	$line2 =~ s/>//;
	$tmp = $line2;
    }else{
	$ch{$tmp} = $line2;
    }
}
close B;
open(C, ">>$output");
open(D, ">>$output2");
while((my $k, my $v) = each(%id)){
    my @x = split(/\t/, $v);
    print D $k, "\t", $v, "\n";
    my $start = ($x[1] - $dist);
    my $len = (($x[2] + $dist) - $start);
    
    my $fasta = substr($ch{$x[0]}, $start, $len); 

    print C ">", $k, "\n", $fasta, "\n";

}
close C;
close D;
