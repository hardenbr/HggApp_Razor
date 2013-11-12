#!/usr/bin/perl

use Time::Local;
use Getopt::Std;

getopts('c:');
$castorList = $opt_c;

$mergedList = "$castorList"."_NoDup.list";

system("head -1 $castorList > $mergedList");

open (CASTORLIST,"<$castorList");
@filecastor=<CASTORLIST>;
close CASTORLIST;

open (LOCALLIST,"$mergedList");
@filelocal=<LOCALLIST>;
close LOCALLIST;

open (MERGEDLIST,">>$mergedList");

foreach $linecastor (@filecastor) { 
   chomp($linecastor);
    if($linecastor =~ /\S+\_(\S+)\_\S+\_\S+\.root/) {
        $jobcastor = $1;
        $found = 0;
        foreach $linelocal (@filelocal) {
            chomp($linelocal);
            if($linelocal =~ /\S+\_(\S+)\_\S+\_\S+\.root/) {
                $joblocal = $1;
                if($jobcastor == $joblocal) { $found = 1; break; }
            }
        }
        if($found == 0) { 
            push(@filelocal,"$linecastor");
            print MERGEDLIST "$linecastor\n"; 
        }
    }
}
