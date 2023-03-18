#!perl -w

use strict;

my $RALL = shift(@ARGV) or die; #RALL.90_500.500screened (or whatever cutoffs I go with
my $output = shift(@ARGV) or die; #I want a list of candidate genes by RALL with transcripts
unlink(qq{$output});

##RALL
##TRINITY_DN41900_c0_g1_i1        214     2L      2       21649996        21650155        21649996.21650085,21650089.21650155

my %RALLsets = (); #I need the base gene info plus all transcripts for each RALL

my %data = (); #RALL information
my %direction = (); #exon strand information
my %matches = (); #individual lines matching RALL
my %mdata = (); #info for individual data
my %remove = ();
my %longest = ();
my %transcripts = ();

open(A, "<$RALL");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    
    my $gene = $a[0];
    $gene =~ s/_i\d+//;
    push(@{$longest{$gene}}, $a[1]); #saving the lengths
    my @g = split(/_/, $a[0]);
    my $trans = $g[4];
    if(!(exists($transcripts{$gene . "_" . $trans}))){ #store each transcript under its gene id
	push(@{$RALLsets{$gene}}, $trans); 
	$transcripts{$gene . "_" . $trans} = 1;
	#print "Found\n";
    }
    my @b = split(/,/, $a[6]);
    
    if(scalar(@b) > 1){
	#	print $line, "\n";
	my @direction = ();
	my %tmp = ();
	foreach my $b (@b){
	    my @c = split(/\./, $b);
	    
	    if($c[0] < $c[1]){
		if(!(exists($tmp{"+"}))){
		    $tmp{"+"} = 1;
		    push (@{$direction{$a[0]}}, "+");
		}
	    }elsif($c[0] > $c[1]){
		if(!(exists($tmp{"-"}))){
		    $tmp{"-"} = 1;
		    push (@{$direction{$a[0]}}, "-");
		}
	    }
	}
    }
    $data{$a[0]} = $line;	
}
close A;

###I want to drop all genes where the longest transcript is < 300bp
while((my $k, my $v) = each(%longest)){
    my $long = 0;
    foreach my $l (@{$longest{$k}}){
	if($l > $long){
	    $long = $l
	}
    }
    if($long < 300){
	$remove{$k} = 1;
    }
}

open(C, ">>$output");
###I need to get rid of duplicate entries
while((my $k, my $v) = each(%RALLsets)){
    #while((my $k, my $v) = each(%data)){
    ###first print the alternate transcripts from RAL
    if(!(exists($remove{$k}))){
	print C "NEW RECORD\n";
	foreach my $t (@{$RALLsets{$k}}){
	    
	    my $trans = $k . "_" . $t;
	    if(!(exists($direction{$trans}))){
		print C $data{$trans}, "\t", "N", "\n";
	    }elsif(exists($direction{$trans})){
		print C $data{$trans}, "\t", join(",",@{$direction{$trans}}), "\n";
	    }
	}
    }
}
close C;
    
