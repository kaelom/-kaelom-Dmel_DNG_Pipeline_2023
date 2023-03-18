#!perl -w

use strict;

my $input = shift(@ARGV) or die; #AG_denovo_candidates_TPM1
my $table = shift(@ARGV) or die; #Filtering/Table_of_Candidates_by_RALL
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %list = ();
my $count = 0;
open(A, "<$input");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    if($line !~ m/ID/){
	$list{$a[0]} = 1;
	$count++;
    }
}
close A;

print "Total = ", $count, "\n";
my $print = 0;

my %double = ();

#RALL    TRINITY_DN10646_c0_g1_i1        353     2L      1       20047979        20048331        20048331.20047979       -


open(B, "<$table");
open(C, ">>$output");
LOOP:while(my $line2 = <B>){
    chomp $line2;
    my @b = split(/\t/, $line2);
    
    if($line2 =~ m/NEW/){
	$print = 0;		
    }else{
	my $name = $b[0];
	$name =~ s/_i\d+//;
	if($print == 0){
	    if(exists($list{$name})){
		$print = 1;
		print C "NEW RECORD\n";
		if($b[3] > 1){			       
		    print C $line2, "\n";
		}elsif($b[3] == 1){
		    
		    print C $b[0], "\t", $b[1], "\t", $b[2], "\t", $b[3], "\t", $b[4], "\t", $b[5], "\t", $b[6], "\t", "N", "\n";
		}
	    }
	}elsif($print == 1){
	    if($b[3] > 1){			       
		print C $line2, "\n";
	    }elsif($b[3] == 1){
		
		print C $b[0], "\t", $b[1], "\t", $b[2], "\t", $b[3], "\t", $b[4], "\t", $b[5], "\t", $b[6], "\t", "N", "\n";
	    }
	}
    }
}
close B;
close C;
    
