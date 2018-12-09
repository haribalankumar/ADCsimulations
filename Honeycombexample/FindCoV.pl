#!/usr/bin/perl

# uses data point locations to place and scale an initial mesh

use strict;

my ($node, $CoV);
my (@concen);


################################################
############## INPUTS #################

my $nodefile1  = 'Time0p4/HoneycombAdvDiffVel5p0_0p4.part0.exnode';

my $time1 = 0.4;

#######################################################


print "open initial node file $nodefile1 \n";

open EXNODE, "<$nodefile1" or die "\033[31mError: Can't open node file\033[0m ";
my $NumberofNodes=0;

my $sumconc =0.0;
my $sdconc =0.0;

my $line = <EXNODE>;
while ($line) {

    if ($line =~ /Node:\s+(\d+)/) {
	$node = $1;
	$NumberofNodes=$NumberofNodes+1;
    
	for (my $i=1;$i<=3;$i++) { ## for x,y,z
	    $line = <EXNODE>;
        }
	    $line = <EXNODE>;  

            $concen[$node] = $line;
            ## print "t1.  $concen  \n";
            $sumconc = $sumconc + $concen[$node];
	}
	
    $line = <EXNODE>;    
}

close EXNODE;

my $avg = $sumconc/$NumberofNodes;

print "number of nodes =   $NumberofNodes \n";
print "\n";
print "------printing average-------\n";
print "$time1  $avg  \n";

for (my $i=1;$i<=$NumberofNodes;$i++) {
    $sdconc = $sdconc + ($concen[$i]-$avg)*($concen[$i]-$avg);
}

$sdconc = sqrt($sdconc/$NumberofNodes);
$CoV = $sdconc / $avg;


print "------printing average, SD, CoV -------\n";
print "$time1 $avg  $sdconc  $CoV \n";



