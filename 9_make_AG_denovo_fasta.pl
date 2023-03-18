#!perl -w

use strict;

my $input = shift(@ARGV) or die;#Sorted_AG_candidates
my $fasta = shift(@ARGV) or die; #Screened_fastas/RALL_filtered.fasta
my $output = shift(@ARGV) or die; #AG_denovo_candidate_TPM*.fasta
unlink(qq{$output});

my %genes = ();

open(A, "<$input");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    if($line !~ m/ID/){
	$genes{$a[1]} = 1; #this is all genes that are kept - now print out all transcripts.	
    }
}
close A;

open(B, "<$fasta");
open(C, ">>$output");
my $found = 0;
while(my $line2 = <B>){
    chomp $line2;
    if($line2 =~ m/^>/){
	$line2 =~ s/^>//;
	my @x = split(/\s+/, $line2);
	my @b = split(/:/, $x[0]);
	$b[0] =~ s/_i\d+$//;
	print $b[0], "\n";
	if(exists($genes{$b[0]})){
	    
            $line2 =~ s/:/-/go;
	    print C ">", $line2, "\n";
	    $found = 1;
	}
    }else{
	if($found == 1){
	    print C $line2, "\n";
	    $found = 0;
	}
    }
}
close B;
close C;
