#!perl -w

use strict;

my $input = shift(@ARGV) or die; ##Table_of_candidates_by_RALL
my $tag = shift(@ARGV) or die; ##RALL, test, etc
my $output = shift(@ARGV) or die; ##GTF formatted AG denovo genes to concatenate to the dmel gtf
unlink(qq{$output});

#RALL    TRINITY_DN35597_c0_g1_i1        596     2R      1       9491047 9491640 9491047.9491640 N

my %gene = ();
my %data = ();

open(A, "<$input");
open(B, ">>$output");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
   # print $a[0], "\n";
    if($a[0] =~ m/$tag/){  #this is a RALL combined transcript
	###base alternate transcripts on different exon models
#	print $a[0], "\n";
       	my $start = $a[4];
	my $stop = $a[5];
	my $strand = ".";
	if($a[4] > $a[5]){
	    $start = $a[5];
	    $stop = $a[4];
	}

	my $gid = $a[0];
	$gid =~ s/_i\d+$//;

	if($a[7] eq "+"){
	    $strand = "+";
	}elsif($a[7] eq "-"){
	    $strand = "-";
	}elsif($a[7] eq "N"){
	    $strand = ".";
	}else{
	    $strand = ".";
	}
	
	if(!(exists($gene{$gid}))){ #print gene info first time seeing new id

	    $gene{$gid} = 1;
	    
	    print B $a[2], "\t", "RALL", "\t", "gene", "\t", $start, "\t", $stop, "\t", ".", "\t", $strand, "\t", ".", "\t", "gene_id ", "\"", $gid, "\"; ", "gene_symbol ", "\"", $gid, "\";", "\n";

	
	    print B $a[2], "\t", "RALL", "\t", "mRNA", "\t", $start, "\t", $stop, "\t", ".", "\t", $strand, "\t", ".", "\t", "gene_id ", "\"", $gid, "\"; ", "gene_symbol ", "\"", $gid, "\"; ", "transcript_id ", "\"TR_", $a[0], "\"; ", "transcript_symbol \"Dmel\\TR_", $a[0], "\";", "\n";

	    my @ex = split(",", $a[6]);

	    foreach my $ex (@ex){
		my @x = split(/\./, $ex);
		
		my $exstart = $x[0];
		my $exstop = $x[1];
		
		if($x[1] < $x[0]){
		    $exstart = $x[1];
		    $exstop = $x[0];
		}

		print B $a[2], "\t", "RALL", "\t", "exon", "\t", $exstart, "\t", $exstop, "\t", ".", "\t", $strand, "\t", ".", "\t", "gene_id ", "\"", $gid, "\"; ", "gene_symbol ", "\"", $gid, "\"; ", "transcript_id ", "\"TR_", $a[0], "\"; ", "transcript_symbol \"Dmel\\TR_", $a[0], "\";", "\n";

	    }
	}elsif(exists($gene{$gid})){
	    
	    print B $a[2], "\t", "RALL", "\t", "mRNA", "\t", $start, "\t", $stop, "\t", ".", "\t", $strand, "\t", ".", "\t", "gene_id ", "\"", $gid, "\"; ", "gene_symbol ", "\"", $gid, "\"; ", "transcript_id ", "\"TR_", $a[0], "\"; ", "transcript_symbol \"Dmel\\TR_", $a[0], "\";", "\n";
	    
	    my @ex = split(",", $a[6]);
	    
	    foreach my $ex (@ex){
		my @x = split(/\./, $ex);
		
		my $exstart = $x[0];
		my $exstop = $x[1];
		
		if($x[1] < $x[0]){
		    $exstart = $x[1];
		    $exstop = $x[0];
		}
		
		print B $a[2], "\t", "RALL", "\t", "exon", "\t", $exstart, "\t", $exstop, "\t", ".", "\t", $strand, "\t", ".", "\t", "gene_id ", "\"", $gid, "\"; ", "gene_symbol ", "\"", $gid, "\"; ", "transcript_id ", "\"TR_", $a[0], "\"; ", "transcript_symbol \"Dmel\\TR_", $a[0], "\";", "\n";
	    }
	}
    }
}
close A;
close B;
