#!/usr/bin/perl

# Stand-alone prototype version of PlatinumCRISPr
# Free for academic use. Please get in contact for commercial licencing options 
# Disclaimer: we decline any responsibility regarding functionality of sgRNAs.

# The programs, library and source code of PlatinumCRISPr are free
# software. They are distributed in the hope that they will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

# Permission is granted for research, educational, and commercial use
# and modification so long as 1) the package and any derived works are not
# redistributed for any fee, other than media costs, 2) proper credit is
# given to the authors and the University of Birmingham

# If you want to include this software in a commercial product, please contact 
# the authors. 



# Author: Roland Arnold r.arnold.2@bham.ac.uk
# for latest versions of the software, please visit https://github.com/rolandA1234/PlatinumCrispR


use strict; use Data::Dumper;
use Time::HiRes qw( time );
my $start_time; my $end_time;
use Cwd 'abs_path';
use File::Find;
#system("cp -r /opt/platinum-crispr/www/local_perl/lib/perl5/lib /opt/platinum-crispr/www/local_perl/lib/perl5/lib2")
my $template="";
while (<IN>){$template.=$_;} close IN;

my $basis_standard="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU";

my $result="";
my $guide="";
my $basis="";
my $origguide="";

my $runmode="standard";
my $rand=rand(1000000);

$guide=$ARGV[0]; chomp($guide);
$guide=~s/[^ATGCatgcUu]//g;
$guide=~s/\s//g;
$origguide=$guide;
$basis=$basis_standard; 



my @accu_result;

my $execpath="";

my %expl;
$expl{'S'}="S: stem";
$expl{'I'}="I: interior loop";
$expl{'B'}="B: bulge";
$expl{'E'}="E: dangling end";
$expl{'H'}="H: hairpin";
$expl{'M'}="M: multi-loop";
$expl{'K'}="K: pseudoKnot";
$expl{'X'}="X: external loop";

my %pairings;
$pairings{'A'}->{'U'}=1;
$pairings{'U'}->{'A'}=1;
$pairings{'G'}->{'C'}=1;
$pairings{'C'}->{'G'}=1;
$pairings{'G'}->{'U'}=1;
$pairings{'U'}->{'G'}=1;

my %expl;


push @accu_result, "guide\t$guide";

$guide=uc($guide);
$guide=~s/ //g;
$guide=~s/T/U/g;
my $valid=0;


if (length($guide)<18) {
	push @accu_result,("Error: too short guide sequence (" .length($guide) . ") - minimum is 18 nucleotides");
	$valid++;
}
if ($guide=~m/[^ATGCU]/) {
	push @accu_result,("Error: unknown nucleotide in sequence");
	$valid++;
}
if (!($guide=~/^G/)){
	push @accu_result,("<b>Warning: no G at start of guide-sequence. This might be correct, depending on the CRISPR system used.</b>");
	#$guide="G$guide";
	push @accu_result, "hadGStart\tno";
	#$valid++;
} else {
	push @accu_result, "hadGStart\tyes";

	#$guide="G$guide"; # remove in later versions, this is for Renn et al benchmark
}

my $seq=$guide.$basis;

push @accu_result, "seq\t$seq";
my $randtmp="$ARGV[1]/$rand/";
mkdir $ARGV[1];
		



if ($valid==0){
	

	(mkdir $randtmp) || die ("cannot generate working directory" . " for guide $guide\n $!");

	$seq=~s/[^ATCGU]//g;
	my $tmpfile1="$rand" . "input";
	chdir ($randtmp); # let's make sure not to get into conflicts with the 'static' filenames 
	(open OUT, ">$tmpfile1") || die ("could not create tmp-file $tmpfile1"); print OUT "$seq\n"; close OUT;
	$start_time = time();

	(open IN, "RNAfold -p $tmpfile1 |") || die ("Error executing RNAfold");  my $seq=<IN>; my $res=<IN>; my $extra=<IN>; my $centroid=<IN>; close IN;
	
	#system("rm $tmpfile1");	
	my @accu_init;
	for (@accu_result){push @accu_init, $_;}

	analyze($res, "");

	$end_time = time();
	


	} else {
		addline("Error: invalid input. Please make sure that the guide RNA length is between 18 and 23, starts with G and does not contain non nucleic-acid elements!")
}

my $elapsed_time=("%.2f\n", $end_time - $start_time);
$result.="Time needed for scan: $elapsed_time seconds<br>\n";

for my $qq (@accu_result){
	if ($qq=~/^overall_decision/){
		print "$guide\t$qq\n";
	}
}


sub analyze(){
	

	my $decision="high";

	my $res=shift;
	

	# ANALYSIS START.
	my ($br, $free)=split(" ", $res);

	$free=~s/[\)\()]//g;
	$free=~s/[\}\{]//g;
	push @accu_result, "freeEnergy\t$free";

	print abs_path($rand . "_output.dbn") . "\n";
	
	(my $tmpfile2=abs_path( $rand . "_output.dbn")) || die ("cannot execute abs_path $!");


	(open OUT, ">$tmpfile2") || die ("could not create file at $tmpfile2"); print OUT ">test\n$seq\n$br\n"; close OUT;
	$seq=~s/[\r\n]//g;
	die ("no dbn file created") if (!(-e $tmpfile2));

	(open IN2, "perl /opt/platinum-crispr/www/scripts/bpRNA.pl $tmpfile2 |") || die ("WARNING cannot execute bpRNA");
   
     my @struc=<IN2>; close IN2;
	

	for (my $i=0;$i<scalar(@struc);$i++){chomp($struc[$i]);}

	my @l=split("",$seq);
	my @l2=split("",$struc[5]);
	
	my $struc_string=join("", @l2);


	my $rev=reverse($struc_string);
	my $revs=reverse($seq);
	my @coloured;
	my %expl;
	
	push @accu_result, "length\t" . length($guide);
	push @accu_result, "bracket\t$br\n";
	push @accu_result, "structure\t" . $struc_string;
	push @accu_result, "sequence\t" . reverse($revs);
	my $structural_2_3="high";
	
	my $tetra_result=test_tetraloop($br, $free);
	$decision="low" if ($tetra_result eq "low"); ### 1
	
	
	my $third_loop_result=test_3dloop($rev);
	
	$structural_2_3="low" if ($third_loop_result eq "low"); 

	my $second_loop_result=test_secondloop($rev);


	
	$structural_2_3="low" if ($second_loop_result eq "low");
	if ($structural_2_3 eq "low"){
		$decision="low"; ### 2
	}
	
	my $mimic=test_bulge_mimic($rev, $revs);
	$decision="low" if ($mimic eq "low"); ###3

	my $hairpin4=test_hairpin_4nt_or_more_in_G19($rev, $revs);
	my $rule6_2=rule6_self_complementarity2($rev, $revs);
	my $rule6_3=rule6_self_complementarity3($rev, $revs);
	my $rule6_4=rule6_self_complementarity4($rev, $revs);
	my $rule6_5=rule6_self_complementarity5($rev, $revs);

	my $accu_rule6="high";
	if ($hairpin4 eq "low"){$accu_rule6="low";}
	if ($rule6_2 eq "low"){$accu_rule6="low";}
	if ($rule6_3 eq "low"){$accu_rule6="low";}
	if ($rule6_4 eq "low"){$accu_rule6="low";}
	if ($rule6_5 eq "low"){$accu_rule6="low";}

	$decision="low" if ($accu_rule6 eq "low"); ###4

	if ($accu_rule6 eq "high"){
			push @accu_result,"rule_combined_6\tcombined_rule6\thigh";
		} else {
			push @accu_result,"rule_combined_6\tcombined_rule6\tlow";
	}


	$decision="low" if (test_gc_seed_rule7_4_or_5($revs) eq "low"); ### 5
	
	

	$decision="low" if (test_gc_seed_boundary($revs, 0.40, 0.80) eq "low"); # ### 6
	

	my $rule10=test_4ntbase($revs);
	$decision="low" if ($rule10 eq "low"); ### 7

	
	my $graf=test_graf_motif_low($revs, $rev, $br);
	$decision="low" if ($graf eq "low"); ### 8
	
	my $rule13=test_13($br);
	$decision="low" if ($rule13 eq "low"); ### 9


	
	push @accu_result, "overall_decision\t$decision";
	return ($decision);
		
}



sub test_gc_seed_rule7_4_or_5(){

	my $revs=shift;	
	my $number=shift;
	my $res="high";
	my $seed=substr($revs,82,6);
	my $gc=gc_contentCount($seed);
	if (($gc>5)||($gc<4)){
		$res="low";
		push @accu_result, "rule7_4_or_5\trule7 $gc $seed\tno\tlow";
	} else {
		push @accu_result, "rule7_4_or_5\trule7 $gc $seed\tyes\thigh";
	}

	return $res;
}


sub test_4ntbase(){
	my $revs=shift;

	my $r=substr($revs,82,4);
	my $res="high";
	

	if ($r=~/UU[UC][UC]/){
		$res="low";
	}
	push @accu_result, "rule10\t4ntbase $r\t$res";
	return $res;
}


sub test_gc_seed_boundary(){
		my $revs=shift;
		my $min=shift;
		my $max=shift;
		my $seed_region=substr($revs,82,102-83);
		my $gc=gc_content($seed_region);
		my $res="high";
		
		if ($gc<=$min){
			$res="low";
			
		}
		if ($gc>=$max){
			$res="low";
			
		}
		
		push @accu_result, "rule8_$min" . "-$max\tgc_seed\t$res";
		return $res;
}

sub test_graf_motif_low{
	my $revs=shift;
	my $rev=shift;
	my $klammer=shift;
	my $klammer_rev=reverse($klammer);
	my $res="high";
	my @log;
	
	my $hp_length=0;
	
	my $totests=substr($revs,95,length($revs)-95);
	my $totest=substr($rev,95,length($rev)-95);
	my $hp_length=0;
	
	for (my $i=0;$i<length($totest);$i++){
			$hp_length++ if (substr($totest,$i,1) eq 'S'); # CHECK
	}
	
	my $totest=substr($klammer_rev,39,3);
	my $sss=$totest;
	
	
	my $totests=substr($revs,82,4);
	
	
	if ($totests=~/UC[CU]G/){
			
		

		if ($sss eq "(((" || $sss eq ")))"){
			$res="low";
			push @log, "39 to 3 do basepairing, UC[CU]G matches at 82";
		}
	}

	
	my $totests=substr($revs,82,4);
	if ($totests=~/C[CU]G[GA]/){
			
		if ($sss eq "(((" || $sss eq ")))"){
			$res="low";
			push @log, "39 to 3 do basepairing, C[CU]G[GA] matches at 82"

		}
		 
	}
	
	my $log=join("-", @log);
	
	push @accu_result, "rule12\tGrafMotif\t$log\t$res";
	return $res;
}


sub test_13(){ 
	my $res="high";
	my $klammer=shift;
	my $klammer_reverse=reverse($klammer);
	my $totest1=substr($klammer_reverse,50,2);
	my $totest2=substr($klammer_reverse,39,2);
	my $pairings=0;
	if ($totest1 eq "(("){
		$pairings++;
	}
	if ($totest1 eq "))"){
		$pairings++;
	}
	if ($totest2 eq "(("){
		$pairings++;
	}
	if ($totest2 eq "))"){
		$pairings++;
	}
	my $res="high";
	$res="low" if ($pairings==2);
	push @accu_result, "rule13\t51_52_40_41 $pairings\tpaired\t$res";

	return $res;
}



sub test_hairpin_4nt_or_more_in_G19(){
	my $rev=shift;
	my $revs=shift;
	
	my $totest=substr($rev, 82, 19);
	my $totests=substr($revs, 82, 19);

		my $res="high";
		
			if ($totest=~m/SSSSHH?H?SSSS/){
				
				$res="low";
		} else {
				
		}
		return $res;

}


sub rule6_self_complementarity2(){
	my $rev=shift;
	my $revs=shift;
	
	my $totest=substr($rev, 82, 19);
	my $totests=substr($revs, 82, 19);

		my $res="high";
		
			if ($totest=~m/SSSSSHHHHSSSSS/){
				
				$res="low";
		} else {
				
		}
		return $res;

}


sub rule6_self_complementarity3(){
	my $rev=shift;
	my $revs=shift;
	
	my $totest=substr($rev, 82, 19);
	my $totests=substr($revs, 82, 19);

		my $res="high";
		
			if ($totest=~m/SSSSSSHHHHHSSSSSS/){
				
				$res="low";
		} else {
				
		}
		return $res;

}

sub rule6_self_complementarity4(){
	my $rev=shift;
	my $revs=shift;
	
	my $totest=substr($rev, 82, 19);
	my $totests=substr($revs, 82, 19);

		my $res="high";
	
			if ($totest=~m/SSSSSSSHHHHHHSSSSSSS/){
				
				$res="low";
		} else {
			
		}
		return $res;

}
sub rule6_self_complementarity5(){
	my $rev=shift;
	my $revs=shift;
	
	my $totest=substr($rev, 82, 19);
	my $totests=substr($revs, 82, 19);

		my $res="high";
		
			if ($totest=~m/SSSSSSSSHHHHHHHSSSSSSSS/){
				
				$res="low";
		} else {
				
		}
		return $res;

}



sub test_tetraloop() {

	my $br=shift;
	my $energy=shift;
	my $br_rev=reverse($br);
	push @accu_result, "bracket_rev\t$br_rev\n";

	my $tetra=reverse( substr($br_rev,52,30) ); 
	push @accu_result, "tetraloop\t$tetra\twith energy $energy\n";
	my $res="high";
		
	if ($tetra eq "(((((((.((((....))))...)))))))"){
		   
	} else {
			
			
			$res="low";
		
	}

	push @accu_result, "rule1\ttetraloopOK\t$tetra\t$res";
	return $res;
}



sub test_secondloop() {
	my $rev=shift;
	my $res="high";
	my $totest=substr($rev,23,10);
	push @accu_result, "secondloop\t$totest\t$rev\n";
	if ($totest eq "SSSHHHHSSS"){
				
			} else {
				$res="low";
		}
	return $res;
}

sub test_3dloop() {
		my $rev=shift;
		
		my $totest=substr($rev,6,15);
		my $res="high";
		push @accu_result,"loop3\t$totest";

		my $goodloop=   "SSSSSSHHHSSSSSS";
		if ($totest eq $goodloop){

				push @accu_result, "rule3\tloop3_OK\tloop3_intact\thigh";
			} else {
				$res="low";
				push @accu_result, "rule3\tloop3_OK -> $rev $totest\tloop3_not_intact\tlow\n";
			
		}

		
	return $res;
}


sub test_bulge_mimic(){
		my $rev=shift;
		my $revs=shift;
		my $res="high";
		my $totest=substr($rev,82,4);
		my $totests=substr($revs,82,4);
	
	
		if ($totests=~/[AG][CU]G[AG]/){ 
			
			
			$res="low";
			push @accu_result, "rule4\tbulgeMimic\tis_mimic\tlow";
		} else {
			
			push @accu_result, "rule4\tbulgeMimic\tis_no_mimic\thigh";

		}
		return $res;
}


sub gc_content($){
	my $seq=shift;
	my $c=0;
	for (my $i=0;$i<length($seq);$i++){
		$c++ if (substr($seq, $i, 1) eq "G");
		$c++ if (substr($seq, $i, 1) eq "C");
	}
	#push @accu_result, "$c $seq <----";
	return $c / length($seq);
}
sub u_contentCount(){
	my $seq=shift;
	my $c=0;
	for (my $i=0;$i<length($seq);$i++){
		
		$c++ if (substr($seq, $i, 1) eq "U");
		$c++ if (substr($seq, $i, 1) eq "T");
	}
	return $c;
}
sub gc_contentCount(){
	my $seq=shift;
	my $c=0;
	for (my $i=0;$i<length($seq);$i++){
		$c++ if (substr($seq, $i, 1) eq "G");
		$c++ if (substr($seq, $i, 1) eq "C");
	}
	return $c;
}
sub g_content(){
	my $seq=shift;
	my $c=0;
	for (my $i=0;$i<length($seq);$i++){
		$c++ if (substr($seq, $i, 1) eq "G");
	}
	return $c / (length($seq));
}

sub addline_nop(){
	my $r=shift;
	$result.="$r";
}
sub addline(){
	my $r=shift;
	$result.="<p>$r</p>";
}
