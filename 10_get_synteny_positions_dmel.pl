#!perl -w

use strict;

my $input = shift(@ARGV) or die; #Sorted_AG_candidates
my $ortho = shift(@ARGV) or die; #/data/FlyRef/dmel_orthologs_in_drosophila_species_fb_2021_02.tsv
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %ch = ();
my %start = ();
my %stop = ();
my %type = ();

open(A, "<$input");
#Type    ID      Ch      Start   Stop    Dist
while(my $line = <A>){
    chomp $line;
    if($line !~ m/Ch/){
	my @a = split(/\t/, $line);
	$ch{$a[1]} = $a[2];
	if($a[4] < $a[3]){
	    print "problem\n";
	}
	$start{$a[1]} = $a[3];
	$stop{$a[1]} = $a[4];
	$type{$a[1]} = $a[0];
    }
}
close A;

my %up = (); #nearest upstream
my %down = ();#nearest downstream
my %updist = ();#how far?
my %downdist = (); #how far?
my %sim = ();
my %yak = ();

my %used = ();

open(B, "<$ortho");
while(my $line = <B>){
    chomp $line;
    if($line !~ m/^#/){
	my @b = split(/\t/, $line);
	if($b[6] =~ m/Dsim/){
	    $sim{$b[0]} = $b[5];
	}
	if($b[6] =~ m/Dyak/){
	    $yak{$b[0]} = $b[5];
	}
	
	if(!(exists($used{$b[0]}))){ #only check each Dmel gene's position once
	    #print $b[0], "\n";
	    $used{$b[0]} = 1;
	    #FBgn0000008     a       2R      22136968..22172834      1       FBgn0141564     Dmoj\GI18825                            EOG091902HN
	    while((my $k, my $v) = each(%ch)){
		if($v eq $b[2]){ #right chrom - check if nearest
		    
		    my @c = split(/\.\./, $b[3]);		    
		    if($c[0] > $c[1]){
			print "Swap\n";
		    }
		    if($c[1] < $start{$k}){			
			my $dist = abs($start{$k} - $c[1]);#how far
			if(!(exists($up{$k}))){
			    $up{$k} = $b[0]; #which gene
			    $updist{$k} = $dist;
			}else{
			    if($dist < $updist{$k}){#is that nearer?
				$up{$k} = $b[0];
				$updist{$k} = $dist;
			    }
			}
		    }
		    
		    if($c[0] > $stop{$k}){
			my $dist = abs($stop{$k} - $c[0]); #how far
			if(!(exists($down{$k}))){
			    $down{$k} = $b[0];
			    $downdist{$k} = $dist;
			}else{
			    if($dist < $downdist{$k}){#is it nearer?
				$down{$k} = $b[0];
				$downdist{$k} = $dist;
			    }
			}
		    }
		}
	    }
	}
    }
}
close B;

open(C, ">>$output");
print C "Type\tID\tCh\tStart\tStop\tMelUp\tDistUp\tMelDown\tDistDown\tSimUp\tSimDown\tYakUp\tYakDown\n";
while((my $k, my $v) = each(%ch)){
    print $k, "\t", $v, "\n";
    if(!(exists($sim{$up{$k}}))){
	$sim{$up{$k}} = "NA";
    }
    if(!(exists($sim{$down{$k}}))){
	$sim{$down{$k}} = "NA";
    }
    if(!(exists($yak{$up{$k}}))){
	$yak{$up{$k}} = "NA";
    }
    if(!(exists($yak{$down{$k}}))){
	$yak{$down{$k}} = "NA";
    }
    print C $type{$k}, "\t", $k, "\t", $v, "\t", $start{$k}, "\t", $stop{$k}, "\t", $up{$k}, "\t", $updist{$k}, "\t", $down{$k}, "\t", $downdist{$k}, "\t", $sim{$up{$k}}, "\t", $sim{$down{$k}}, "\t", $yak{$up{$k}}, "\t", $yak{$down{$k}}, "\n";
}
close C;
