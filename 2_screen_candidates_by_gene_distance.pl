#!perl -w

use strict;

####Edit this to find the maximum size of a gene based on all transcripts and do filtering based on those positions!!!!!

##I need to take the positions of genes in the R??.positions file and screen them against being within 500 bp of a gene
my $screen = shift(@ARGV) or die; #Dmel_gene_500extend
my $input = shift(@ARGV) or die; #/data/julie/AG_Denovo/RALL/*.positions
my $output = shift(@ARGV) or die;
unlink(qq{$output});

##Dmel_gene_500extend
#FBgn0031081     X       19960797        19969823        +

##positions
#TRINITY_DN3610_c0_g1_i1 441     2L      1       13152277        13152717        13152717.13152277

my %ch = ();
my %start = ();
my %stop = ();

open(A, "<$screen"); ##this stores the information for regions 500bp up and downstream from genes (taken from the flybase file)
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    $ch{$a[0]} = $a[1];
    $start{$a[0]} = $a[2];
    $stop{$a[0]} = $a[3];
}
close A;

my %remove = (); #set to screen out due to proximity
my %all = (); #the set of all candidates

###First go through and screen and record ones to screen out - then go through and print

open(P, "<$input");
while(my $line = <P>){
    chomp $line;
    my @b = split(/\t/, $line);
    if($line !~ m/ID/){
	$all{$b[0]} = $line;
	my $id = $b[0];
	$id =~ s/_i\d+//;
	
	while((my $k, my $v) = each(%ch)){
	    if($v eq $b[2]){ ##same chrom
		if((($start{$k} <= $b[4]) and ($stop{$k} >= $b[4])) or
		   (($start{$k} <= $b[5]) and ($stop{$k} >= $b[5])) or
		   (($start{$k} <= $b[4]) and ($stop{$k} >= $b[5])) or
		   (($start{$k} >= $b[4]) and ($stop{$k} <= $b[5]))){
		    
		    if(!(exists($remove{$id}))){
			#	print $pos[0], "\t", $b[0], "\n";
			$remove{$id} = 1;
			#	print  $line, "\t", $k, "\t", $v, "\t", $start{$k}, "\t", $stop{$k}, "\n";
		    }
		}
	    }
	}
    }
}
close P;


open(D, ">>$output");
while((my $k, my $v) = each(%all)){
    #   print $k, "\n";
    $k =~ s/_i\d+$//;
    # print $k, "\n";
    if(!(exists($remove{$k}))){
	print D $v, "\n";
    }
}
close D;

