#!perl -w

use strict;

#I want to check that regions are syntenic - look for nearby genes in some window as an initial pass
#For intergenic ones I need nearest gene - for intronic ones I should be able to use the intron info to find orthologs

my $set = shift(@ARGV) or die; #Intron_synteny_regions
my $indir = shift(@ARGV) or die; #the directory with the alignments to the outgroups & dmel - I need the positions so I can look at genes in that region
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %g = (); #the list of genes
open(A, "<$set");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    $g{$a[0]} = $a[1] . "\t" . $a[2] . "\t" . $a[3];
}
close A;

my @sp = (); #hold species specific info for chromosome matches
opendir DIR, "$indir";
my @files = grep {/align/} readdir DIR;
closedir DIR;
my %ch = ();
my %pos = ();
my %found = ();
foreach my $file (@files){  ##AG_plus5kb.to.Dsim.align
    my @f = split(/\./, $file);

#    print $f[2], "\n";
    push(@sp, $f[2]); ##this is species info
    $file = $indir . $file;
    open(F, "<$file");
    while(my $line = <F>){
	chomp $line;      
	my @b = split(/\t/, $line);
	$b[0] =~ s/_i\d+//;
	if(exists($g{$b[0]})){
	    my @c = split(/\t/, $g{$b[0]});
	    $b[1] =~ s/Scf_//;	   
	    if($b[1] =~ m/$c[0]/){#I want to only look at things on matching chromosomes - since I'm looking at much bigger regions than just the denovo transcript
#		print $b[1], "\t", $c[0], "\n";
		if(($b[1] !~ m/random/) and ($b[1] !~ m/NODE/) and ($b[1] !~ m/chrUn/) and ($b[1] !~ m/chrUh/)){ #not considering matches to small unplaced chromosomes - mostly b/c many of these seem to be not real, but parts of other chromosomes
		    if(($b[2] >= 80) and ($b[3] >= 1000)){
			
			if(!(exists($found{$f[2] . "\t" . $b[0] . "\t" . $b[1]}))){
			    push(@{$ch{$f[2] . "\t" . $b[0]}}, $b[1]);
			    $found{$f[2] . "\t" . $b[0] . "\t" . $b[1]} = 1;
			}	       		    
			my $pos = $b[8] . "." . $b[9];		    
			if($b[9] < $b[8]){
			    $pos = $b[9] . "." . $b[8];
			}		
			###species \t geneID
			push(@{$pos{$f[2] . "\t" . $b[0]}}, $pos);
		    }
		}
	    }
	}
    }
    close F;
}

open(B, ">>$output");
my %tstart = ();
my %tstop = ();
while((my $k, my $v) = each(%g)){
    my %t = ();
    foreach my $s (@sp){		
	if(exists($ch{$s . "\t" . $k})){ ##there is an outgroup match
	    
	    if(scalar(@{$ch{$s . "\t" . $k}}) == 1){ #make sure there is only one chromosome
		
		my @sorted = sort {$a <=> $b} @{$pos{$s . "\t" . $k}};		
		
		foreach my $sort (@sorted){
		    my @x = split(/\./, $sort); #these position pairs are already sorted
		    if(!(exists($tstart{$s . "\t" . $k}))){
			$tstart{$s . "\t" . $k} = $x[0];	
			$tstop{$s . "\t" . $k} = $x[1];
		    }elsif(exists($tstart{$s . "\t" . $k})){
			foreach my $x (@x){
			    if($x < $tstart{$s . "\t" . $k}){
				$tstart{$s . "\t" . $k} = $x;
			    }
			    if($x > $tstop{$s . "\t" . $k}){
				$tstop{$s . "\t" . $k} = $x;
			    }
			}
		    }
		}
		print B $s, "\t", $k, "\t", @{$ch{$s . "\t" . $k}}[0], "\t", $tstart{$s . "\t" . $k}, "\t", $tstop{$s . "\t" . $k}, "\t", abs($tstart{$s . "\t" . $k} - $tstop{$s . "\t" . $k}), "\n"; #"\t", join(",", @{$pos{$s . "\t" . $k}}), "\n";
	    }else{ #multiple chromosome matches - let's look at these
	#	print $s, "\t", $k, "\t", join(",", @{$ch{$s . "\t" . $k}}), "\n";
	    }
	}else{
	    print $s, "\t", $k, "\t", "missing\n";
	}
    }
    
    %t = ();
    
}
close B;
