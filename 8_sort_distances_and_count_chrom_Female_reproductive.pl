#!perl -w

use strict;

my $in1 = shift(@ARGV) or die;
my $drop = shift(@ARGV) or die; #the distance to drop candidates if they are this close or closer
my $t = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
unlink(qq{$output});



my %ch = ();
my %pos = ();
my %data = ();
my %type = ();

my @input = ($in1);
foreach my $input (@input){
    my @f = split(/\./, $input);
    open(A, "<$input");
    while(my $line = <A>){
	chomp $line;
	my @a = split(/\t/, $line);
	if($line =~ m/TRINITY/){
	    my $id = $a[0];
	    $id =~ s/_i\d+//;
	    if(!(exists($pos{$id}))){
		$type{$id} = $f[0];
		$pos{$id} = $a[2] . "." . $a[4] . "." . $a[5];
	    }elsif(exists($pos{$id})){
		my @t = split(/\./, $pos{$id});
		if($a[4] > $a[5]){
		    print "problem\n";
		}			
		#start with previous values
		my $tstart = $t[1];
		my $tstop = $t[2];		
		if($a[4] < $t[1]){
		    $tstart = $a[4];
		}
		if($a[5] > $t[2]){
		    $tstop = $a[5];
		}
		$pos{$id} = $t[0] . "." . $tstart . "." . $tstop;
	    }
	}
    }
}
close A;

while((my $k, my $v) = each(%pos)){
    my @p = split(/\./, $v);
    my $pos = $p[1] . "." . $p[2];
    push(@{$ch{$p[0]}}, $pos);
    $data{$pos} = $k;
}

open(B, ">>$t");
#print B "Type\tGeneID\tCh\tStart\tStop\tDist\n";
while((my $k, my $v) = each(%ch)){
    my @sort = sort {$a <=> $b} @{$ch{$k}};
    my $dist = -99;
    my $count = 0;
    foreach my $s (@sort){
	$count++;
        my @m = split(/\./, $s);
	if(exists($sort[$count])){
	    my @n = split(/\./, $sort[$count]);
	    
	    $dist = $n[0] - $m[1];
	    print B $type{$data{$s}}, "\t", $data{$s}, "\t", $k, "\t", $m[0], "\t", $m[1], "\t", $dist, "\n";
	}else{
	    print B $type{$data{$s}}, "\t", $data{$s}, "\t", $k, "\t", $m[0], "\t", $m[1], "\t", "-99", "\n";
	}
    }
}
close B;

##I need to drop events withing 500bp of each other and also less than 300bp in length
##Intronic        TRINITY_DN34701_c0_g1   2R      7167911 7168732 302
my %list = ();##the set
my %drop = ();

my $previous = 0;

open(A, "<$t");
 LOOP:while(my $line = <A>){
     chomp $line;
     my @a = split(/\t/, $line);
     #     if($a[2] !~ m/Ch/){
     $list{$a[1]} = 1;
     if($previous == 1){ ##the previous one was too close
	 $drop{$a[1]} = 1; #drop the current one too
	 if($a[5] >= $drop){
	     $previous = 0; #next one is good
	     next LOOP;
	 }elsif(($a[5] < $drop) and ($a[5] ne -99)){
	     $previous = 1; #next one is also too close
	     next LOOP;
	 }elsif($a[5] eq -99){
	     $previous = 0; #reset and move to the next chromosome
	     next LOOP;
	 }
     }elsif($previous == 0){
	 if($a[5] >= $drop){ #next one is fine
	     $previous = 0;
	     next LOOP;
	 }elsif(($a[5] < $drop) and ($a[5] ne -99)){
	     $drop{$a[1]} = 1;
	     $previous = 1;
	     next LOOP;
	 }elsif($a[5] eq -99){
	     $previous = 0;
	     next LOOP;
	 }
     }
}
close A;

open(A, "<$t");
open(B, ">>$output");
while(my $line = <A>){
    chomp $line;
    my @b = split(/\t/, $line);
    if(!(exists($drop{$b[1]}))){
	print B $line, "\n";
    }
}
close B;
close A;

unlink(qq{$t});
