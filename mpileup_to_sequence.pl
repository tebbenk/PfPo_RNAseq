#!/usr/bin/perl -w                                                                                                                                             \
                                                                                                                                                                            
use strict;
use warnings;

my @mpileup_files = </local/projects-t3/SerreDLab-3/kieran.tebben/mali/ovale/po_snp_files/*sorted_subset.mpileup>;
my $mpileup;
my @read;
my $position;
my $allele;
my $RAF; 
my $reference;


for my $mpileup (@mpileup_files) {
#	print "$mpileup\n";
	my @files = split /\//, $mpileup;
	my $sample_name = $files[8];	
	(my $sample = $sample_name) =~ s/\.[^.]+$//;
	$sample =~ s/\.mpileup*//;
	print "$sample\n";
	
	my $out = $sample."_cytb_sequence.txt";
	
	open(OUT, ">", "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/ovale/po_snp_files/$out");
	
	chomp $sample;
	my $mpileup_path = "/local/projects-t3/SerreDLab-3/kieran.tebben/mali/ovale/po_snp_files/$sample_name";
	open(IN, $mpileup_path);

	my $last = 144;	
	while(my $line = <IN>){
		my $reference = 0;
		my $A = 0;
		my $T = 0;
		my $C = 0;
		my $G = 0;
#		print "$last\n";
		chomp $line;
		@read = split(/\t+/, $line);  ## splitting line into array
		if($read[0] eq "PocGH01_MIT_v2"){
#			print "$last\n";
#			print "$read[0]\n";
			if(($read[1] > 145) && ($read[1] <= 1275)){
				$last++;
				if($read[1] != ($last+1)){
#					print "$last\t";
					print OUT "N";
				}
				else{
					for($position = 0; $position < length $read[4]; $position++){
						$allele = substr($read[4], $position, 1);
						if($allele eq "\." or $allele eq "\," or $allele eq "\>" or $allele eq "\<"){
							$reference++;
						}
						elsif($allele =~ "[Aa]"){
							$A++;
						}
						elsif($allele =~ "[Tt]"){
							$T++;
						}
						elsif($allele =~ "[Cc]"){
							$C++;
						}
						elsif($allele =~ "[Gg]"){
							$G++;
						}
					}
				}
				#print "$A, $C, $G, $T\n";
				$RAF = $reference/$read[3];
				if($RAF >= 0.5){
					print OUT "$read[2]";
#					print "Ref: $read[1]\t$read[2]\t$read[2]\n";
				}
				elsif(($A > $T) && ($A > $C) && ($A > $G)){
					print OUT "A";
#					print "Alt: $read[1]\t$read[2]\tA\n";
				}
				elsif(($T > $A) && ($T > $C) && ($T > $G)){
					print OUT "T";
#					print "Alt: $read[1]\t$read[2]\tT\n";
				}
				elsif(($G > $A) && ($G > $C) && ($G > $T)){
					print OUT "G";
#					print "Alt: $read[1]\t$read[2]\tG\n";
				}
				elsif(($C > $A) && ($C > $G) && ($C > $T)){
					print OUT "C";
#					print "Alt: $read[1]\t$read[2]\tC\n";
				}
			}		
		}
	}
	close(OUT);
}
exit;