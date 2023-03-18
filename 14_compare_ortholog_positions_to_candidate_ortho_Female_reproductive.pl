#!perl -w

use strict;

my $input = shift(@ARGV) or die; #D*_synteny_positions - this is the larger range - it has the furthest positions of the upstream and downstream genes
###TRINITY_DN33386_c0_g1   FBgn0181924     FBgn0182359     2R      5328670 5381181
my $blast = shift(@ARGV) or die; #Synteny_ortholog_positions
#Dsim    TRINITY_DN42156_c0_g1   3R      8518107 8527615 9508
#Dyak    TRINITY_DN42156_c0_g1   3R      11650582        11660466        9884
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %pos = ();
my $count = 0;
open(A, "<$input");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    if($line !~ m/NA/){
	my $start = $a[4];
	my $stop = $a[5];
	#  print $line, "\n";
	if($a[5] < $a[4]){
	    $start = $a[5];
	    $stop = $a[4];
	}
	$pos{$a[0]} = $a[3] . "\t" . ($start - 10000) . "\t" . ($stop + 10000);
	$count++;
    }
}
close A;
print "Total Found = ", $count, "\n";

my @w = split(/_/, $input);
print $w[2], "\n";

my %found = ();

open(C, ">>$output");
open(B, "<$blast"); #this is the summarized positions of the blast output of the candidates +/- some to the ortholog genomes
while(my $line = <B>){
    chomp $line;
    my @b = split(/\t/, $line);
    if($b[0] eq $w[2]){ #looking at the right species	
	if(exists($pos{$b[1]})){	    
	    my @p = split(/\t/, $pos{$b[1]});##Ch    Start    Stop
	    if(($p[0] eq $b[2]) and ($p[1] < $b[3]) and ($p[2] > $b[4])){
		$found{$b[1]} = 1;
		print C $pos{$b[1]}, "\t", $line, "\n";
		
	    }else{
		print "no overlap\t",  $b[1], "\n";
	    }
	}else{
#	    print "Missing\t", $b[1], "\n";
	}
    }
}
close B;
close C;

while((my $k, my $v) = each(%pos)){
    if(!(exists($found{$k}))){
	print $k, "\n";
    }
}
