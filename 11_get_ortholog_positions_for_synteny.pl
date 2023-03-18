#!perl -w

use strict;

my $sp = shift(@ARGV) or die; #sim or yak
my $input = shift(@ARGV) or die;
my $ranges = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
unlink(qq{$output});
my $output2 = shift(@ARGV) or die;
unlink(qq{$output2});

my %up = ();
my %down = ();
my %revup = ();
my %revdown = ();
my $col1 = 0;
my $col2 = 0;
if($sp eq "sim"){
    $col1 = 9;
    $col2 = 10;
}elsif($sp eq "yak"){
    $col1 = 11;
    $col2 = 12;
}else{
    die;
}
open(A, "<$input"); ##pull out the right genes
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    if($line !~ m/Type/){
	$up{$a[1]} = $a[$col1];
	push(@{$revup{$a[$col1]}}, $a[1]); #reverse lookup more than one candidate can have the same nearest gene
	$down{$a[1]} = $a[$col2];
	push(@{$revdown{$a[$col2]}}, $a[1]);
    }
}
close A;
my %upch = ();
my %uppos = ();
my %downch = ();
my %downpos = ();
open(B, "<$ranges");
while(my $line2 = <B>){
    chomp $line2;
    my @b = split(/\t/, $line2);
    $b[1] =~ s/Scf_//;
    my $start = $b[2];
    my $stop = $b[3];
    if($b[2] > $b[3]){
	$start = $b[3];
	$stop = $b[2];
    }
    if(exists($revup{$b[0]})){	#need this one?
	foreach my $x (@{$revup{$b[0]}}){
	    $upch{$x} = $b[1];
	    $uppos{$x} = $start;
	}
    }
    if(exists($revdown{$b[0]})){
	foreach my $y (@{$revdown{$b[0]}}){
	    $downch{$y} = $b[1];
	    $downpos{$y} = $stop;
	}
    }
}
close B;
open(C, ">>$output");
open(D, ">>$output2");
while((my $k, my $v) = each(%up)){
    
    if(($up{$k} ne "NA") and ($down{$k} ne "NA")){
	if($upch{$k} eq $downch{$k}){ #sanity check here
	    print C $k, "\t", $v, "\t", $down{$k}, "\t", $upch{$k}, "\t", $uppos{$k}, "\t", $downpos{$k}, "\n";
	    
	}else{ ##this is the set I have to check by hand
	    print D $k, "\t", $v, "\t", $down{$k}, "\t", $upch{$k}, "\t", $downch{$k},  "\n";
	}
    }else{
	if($up{$k} eq "NA"){
	    $upch{$k} = "NA";
	    $uppos{$k} = "NA";
	}
	if($down{$k} eq "NA"){
	    $downpos{$k} = "NA";
	}
	##Print the positions I have for ones where there is not an ortholog next
	print D $k, "\t", $v, "\t", $down{$k}, "\t", $upch{$k}, "\t", $uppos{$k}, "\t", $downpos{$k}, "\n";
	    
    }
}
close C;
close D;
