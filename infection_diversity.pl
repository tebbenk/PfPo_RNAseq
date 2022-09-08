#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

##open mpileup files 
open (SAMPLES, "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/sample_lists/ovale_IDs.txt");
my @sample_names;

while(<SAMPLES>) {
	chomp;
	push @sample_names, $_;
}

print join "\n", @sample_names;

my %major_allele;
my $sample;
my $location;
my $n = 0;


##Add each position in each sample to the hash as keys
for my $sample (@sample_names){
	print "sample: $sample\n";
	my $mpileup_path = "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/${sample}.pf.subset.rmdup_subset.mpileup";
	
#	print "$mpileup_path\n";
 	open(IN, $mpileup_path);
# 	print "opened\n";
 	$n++;
 	##Assign major allele for each position in each mpileup file and store as hash values
 	
 	while(my $line = <IN>){
 	chomp $line;
#	print "$line\n";
#	sleep(1);
	$c++;
	if ($c == 10000) {last}

	my @position = split(/\t+/, $line);
	my $ref_allele = $position[2];
	my $chr = $position[0];
	my $pos = $position[1];
	$location = $chr."$pos";
	my $allele;
		if($position[3] >= 20){
			my $reference = 0;
			my $RAF;
			for(my $call = 0; $call < length $position[4]; $call++){
				$allele = substr($position[4], $call, 1);
				if($allele eq "\." or $allele eq "\," or $allele eq "\>" or $allele eq "\<"){
					$reference++;
#					print "$reference\n";
				}
			}
			$RAF = $reference/$position[3];
			if($RAF >= 0.5){
				$major_allele{$location}[$n] = $ref_allele;
#				print "$sample\n";
#				print "reference: $major_allele{$location}[$n]\n";
#				print "key: $location\n";
#				sleep(1);
			}
			else{
				if($allele =~ /[ACTGactg]/){
					$major_allele{$location}[$n] = $allele;
#					print "alt allele: $major_allele{$location}[$n]\n";
#					sleep(1);
				}
			}		
		}	
	}
	close(IN);
}

my %tot;
my %diff;
my $matrix_pos;

##Iterate through samples in an array and compare them 
foreach my $bp (keys %major_allele){
#	print "position: $bp\n";
	for(my $s = 1; $s <= $n; $s++){
		for(my $y = ($s+1); $y <=$n; $y++){
#				print "Sample 1: $major_allele{$bp}[$s]\n";
#				print "Sample 2 $major_allele{$bp}[$y]\n";
			if((defined $major_allele{$bp}[$s]) && (defined $major_allele{$bp}[$y])){
#				print "Sample 1: $major_allele{$bp}[$s]\n";
#				print "Sample 2 $major_allele{$bp}[$y]\n";
				$matrix_pos = $s."_$y";
				$tot{$matrix_pos}++;
#				print "Total: $tot{$matrix_pos}\n";
				if($major_allele{$bp}[$s] ne $major_allele{$bp}[$y]){
					$diff{$matrix_pos}++;
#					print "Sample 1: $major_allele{$bp}[$s]\n";
#					print "Sample 2: $major_allele{$bp}[$y]\n";
#					print "Differences: $bp\t $matrix_pos $diff{$matrix_pos}\n";
#					print OUT "Differences: $bp\t $matrix_pos $diff{$matrix_pos}\n";
				}
			}
			
		}			
	}
}

open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/number_differences_symptomatic.txt");

for(my $s = 1; $s <= $n; $s++){
	for(my $y = ($s+1); $y <=$n; $y++){
		$matrix_pos = $s."_$y";
		print "Differences: $matrix_pos $diff{$matrix_pos}\n";
		print "Total: $matrix_pos $tot{$matrix_pos}\n";
		print OUT "Differences: $matrix_pos $diff{$matrix_pos}\n";
		print OUT "Total: $matrix_pos $tot{$matrix_pos}\n";
	}
}

close(OUT);
exit; 




