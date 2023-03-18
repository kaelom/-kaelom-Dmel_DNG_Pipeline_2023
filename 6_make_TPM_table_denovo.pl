#!perl -w

use strict;

my $indir = shift(@ARGV) or die;
my $min = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
unlink(qq{$output});

opendir DIR, "$indir";
my @abund = grep {/abund\.tab/} readdir DIR;
closedir DIR;

my %gene = ();
my %data = ();
my @lines = ();
foreach my $abund (@abund){
    print $abund, "\n";
    my @ab = split(/_/, $abund);
    my $id = $ab[0] . "_" . $ab[1] . "_" . $ab[2];
    push (@lines, $id);

    $abund = $indir . $abund;
    open(A, "<$abund");
    while(my $line = <A>){
	chomp $line;
	my @a = split(/\t/, $line);
	if($line =~ m/TRINITY/){
	    $gene{$a[0]} = 1;
	    $data{$id . "\t" . $a[0]} = $a[8];
	}
    }
    close A;
}
my %min = ();
while((my $k, my $v) = each(%gene)){
    $min{$k} = 0;
    foreach my $ral (@lines){
	if($data{$ral . "\t" . $k} >= $min){
	    $min{$k}++;
	}
    }
}

my @sorted = sort (@lines);

open(B, ">>$output");
print B "ID";
foreach my $ral (@sorted){    
    print B "\t", $ral;
}
print B "\n";
while((my $k, my $v) = each(%gene)){
    if($min{$k} >= 1){
	print B $k;
	foreach my $ral (@sorted){
	    print B "\t", $data{$ral . "\t" . $k};
	}
	print B "\n";
    }
}
close B;
