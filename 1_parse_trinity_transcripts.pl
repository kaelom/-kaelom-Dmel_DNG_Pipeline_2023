#!perl -w

use strict;

##What I want to do is screen out all transcripts that match coding sequence or TEs then get a set of candidates to screen against other species as well as check for positions and compare to stringtie results

##Make sure that the level of screening is at the gene not the transcript

my $dir = shift(@ARGV) or die; #/data/julie/AG_Denovo/blast/
my $percent = shift(@ARGV) or die; #the percent alignment match
my $length = shift(@ARGV) or die; #the length of the  match
my $r = shift(@ARGV) or die; #the base name
my $out = shift(@ARGV) or die; #output directory
my $fasta = shift(@ARGV) or die; #fullpath to fasta
my $tag = $percent . "_" . $length;

#TRINITY_DN15687_c0_g1_i1        X       100.000 525     0       0       1       525     7013435 7012911 0.0     970

#foreach my $r (@ral){ #storing the transcripts from the fasta file
print $r, "\n";
#my $fasta = "/home/kdlombardo/hayleys/assemblies_and_blastout/" . $r . "/" . $r . ".Trinity.fasta";

my %len = ();
my %seq = ();

my $tlen = ();
open(F, "<$fasta");
while(my $line = <F>){
    chomp $line;
    if($line =~ m/^>/){
        my @x = split(/\s+/, $line);
        $x[0] =~ s/^>//;
        $tlen = $x[0];
    }else{
        $len{$tlen} = length($line);
        $seq{$tlen} = $line;
    }
}
close F;

    
my %transcript = ();#all transripts -dmel
my %transch = ();##set of chromosomes for good transcript hits
my %transdata = (); #store info to print if not screened out
my %drop = (); #this is the set of transcripts to drop due to alignment problems

my %intergenic = ();#intergenic - should be  good

my %used = (); #all other
my %screen = ();

#my $ch = $r . ".chromosome_outfmt6.txt";
#my $intergenic = $r . ".intergenic_outfmt6.txt";

opendir DIR, "$dir";
my @align = grep {/outfmt6.txt/} readdir DIR;
closedir DIR;

my $output =  $out . $r . "." . $tag . ".transcripts"; #the set of things to keep
my $screenout = $out . $r . "." . $tag . ".screened"; #the reason things that were screened out were screened out
my $position = $out . $r . "." . $tag . ".positions"; #positions of events
unlink(qq{$output});
unlink(qq{$screenout});
unlink(qq{$position});

foreach my $align (@align){  #for each alignment file Flybase and outgroup trinity transcripts
    my @l = split(/\./, $align);
    $l[1] =~ s/\-r\d+//;
    print $l[0], "\t", $l[1], "\n";
    my $id = $l[0];
    $align = $dir . $align;
    if($id eq $r){
	
	if($align =~ m/chromosome/){ #####update here to do what previous positon file did !!!!!
	    ####I want to keep track of alignments with multiple hits
	    if($align =~ m/dmel/){  #to look at the Dmel candidates
		print $align, "\n";
		my %ch = (); #temporary storage for ch
		my %pos = (); #temporary storage for pos
		my %cset = ();
		open(A, "<$align"); #store alignment info
		while(my $line = <A>){
		    chomp $line;
		    my @a = split(/\t/, $line);
		    if(($a[2] > 80) and ($a[3] > 30)){##Keep this constant as this is matching to the Dmel chromosome!!!
			$transcript{$a[0]} = 1; #just storing for now
			$drop{$a[0]} = 0;#initialize
			if(!(exists($cset{$a[0] . "\t" . $a[1]}))){##this checks for multiple chromosome matches
			    push(@{$ch{$a[0]}}, $a[1]);
			    $cset{$a[0] . "\t" . $a[1]} = 1;
			}			    
			my $pos = $a[8] . "." . $a[9];
			push(@{$pos{$a[0]}}, $pos);			    
		    }
		}
		
		###Now go through all transcripts and drop
		
		while((my $k, my $v) = each(%transcript)){
		    if(scalar(@{$ch{$k}}) > 1){
			$drop{$k} = 1;
			#	    print $r, " Drop ", $k, "\n";
		    }else{
			$transch{$k} = $ch{$k}[0];
		    }
		    
		    my $first = "X";
		    my $last = "X";
		    
		    foreach my $s (@{$pos{$k}}){
			my @k = split(/\./, $s);
			if($first eq "X"){				
			    $first = $k[0];
			    $last = $k[1];
			    if($k[0] > $k[1]){
				$first = $k[1];
				$last = $k[0];
			    }
			}else{
			    foreach my $k (@k){
				if($k < $first){
				    $first = $k;
				}
				if($k > $last){
				    $last = $k;
				}
			    }
			}
		    }
		    
		    if(($last - $first) > 30000){
			$drop{$k} = 1;
		    }
		    if($drop{$k} == 0){
			my @sort = sort {$a <=> $b} @{$pos{$k}};
			
			$transdata{$k} = $k . "\t" . $len{$k} . "\t" . $transch{$k} . "\t" . scalar(@sort) . "\t" . $first . "\t" . $last . "\t" . join(",", @sort);
			
		    }else{
			#print $k, "\t", $first, "\t", $last, "\t", join(",", @{$ch{$k}}), "\t", join(",", @{$pos{$k}}), "\n";
		    }
		}
		%ch = ();
		%pos = ();
	    }
	    
	}elsif($align =~ m/intergenic/){
	    print $align, "\n";
	    open(B, "<$align");
	    while(my $line = <B>){
		chomp $line;
		my @b = split(/\t/, $line);
		$intergenic{$b[0]} = "intergenic"; #want to check these - may be useful
	    }
	    
	}else{###now we are doing the screening out process - I want to do this at the gene level rather than the transcript level
	    #	print $align, "\n";
	    open(C, "<$align");
	    while(my $line = <C>){
		chomp $line;
		my @c = split(/\t/, $line);
		if($align !~ m/dmel/){
		    if(($c[2] >= $percent) and ($c[3] >= $length)){
			my $gene = $c[0];
			$gene =~ s/_i\d+//; #convert to gene ID here
			if(!(exists($used{$gene . "." . $l[1]}))){
			    $used{$gene . "." . $l[1]} = 1;			    
			    push(@{$screen{$gene}}, $l[1]);
			}
		    }
		}else{
		    if(($c[2] >= 80) and ($c[3] >= 30)){
			my $gene = $c[0];
			$gene =~ s/_i\d+//; #convert to gene ID here
			if(!(exists($used{$gene . "." . $l[1]}))){
			    $used{$gene . "." . $l[1]} = 1;			    
			    push(@{$screen{$gene}}, $l[1]);
			}
		    }
		}
	    }
	}
    }
}	
open(D, ">>$output");
open(E, ">>$screenout");
open(X, ">>$position");
my %printed = ();
print X "ID\tTransLen\tCh\tExonNum\tStart\tStop\tExons\n";
while((my $k, my $v) = each(%transcript)){
    my $gene = $k;
    $gene =~ s/_i\d+//;
    if(exists($transdata{$k})){
	my @t = ();
	if(!(exists($screen{$gene}))){
	    
	    print D ">", $k, ":";
	    if(exists($intergenic{$k})){
		print D $intergenic{$k}, ":";
	    }else{
		print D "NA:";
	    }
	    print D $len{$k}, "\n";
	    print D $seq{$k}, "\n";
	    
	    print X $transdata{$k}, "\n";
	}else{
	    if(!(exists($printed{$gene}))){
		print E $k, "\t", join(",", @{$screen{$gene}}), "\n";
		$printed{$gene} = 1;
	    }
	}
    }else{
	#   print $r, "\t", $k, "\n";
    }
}
close D;
close E;
close X;
%used = ();
%transcript = ();
%transch = ();
%transdata = ();
%drop = ();
%screen = ();
%intergenic = ();
%len = ();
%seq = ();

