#!/bin/perl -w

use Bio::SearchIO;

$evalue=0.00001;
my $arquivo = 'annotation_v2.gff';
open(my $fh, '>', $arquivo) or die "It was not possible open the file '$arquivo' $!";
print"\ncreating gff file...\n";
$fasign="";
foreach $argnum (0 .. $#ARGV) {
	print "\nprocessing ",$ARGV[$argnum]," file...\n";

	$file = $ARGV[$argnum];	
	my $searchio = Bio::SearchIO->new( -format => 'blast', -file   => $file );
	while ( my $result = $searchio->next_result() ) {
		#print $result->query_name."\t".$result->query_description."\n";
		$aux2=$result->query_description;
		while( my $hit = $result->next_hit ) {
			#print $hit->description."\n";    
	    		$aux= $hit->description;
			#print $hit->frac_aligned_query()."\n";
			if($aux=~ m/hypothetical/ig){
				next;
			}else{
				if($hit->frac_aligned_query() >= 0.95 && $hit->significance <= $evalue ) {
					#copy description to /product
					$fasign="assigned function";
						#last;
				}else{
					if($hit->frac_aligned_query() <=0.95 && $hit->significance <= $evalue ){
						#putative
						$fasign="putative";						
						#last;
					}else{
						if($hit->frac_aligned_query() <=0.95 &&$hit->significance >= $evalue ) {
							#hypothetical
							$fasign="hypothetical";
							#last;
						}
					}
				}

				#writing on gff

				print"\nwriting ",$result->query_name," on gff \n";	
				$aux2=$result->query_description;
				$hebra="-";
				if($aux2 =~ m/plus/ig){
					$hebra="+";
				}	
				$locations=$result->query_description;
				#print $locations,"\n";				
				#get name of hit before [mycoplasma genitalium]
				#$name = $hit->name;
				#process string to get locations of start and end
				($start)= ($locations =~ /([0-9]+)/);
				($end)=($locations =~ m/(\d+)\s*---/);
				($fr)=(substr($locations, -1));
        			print $fh $result->query_name,"\tbiomagic\t.\t", $start,"\t",$end,"\t.\t",$hebra,"\t",$fr,"\t",$fasign,"\n";
				last;	
			}


		}
			#print $hit->description."\n";
	}
}
print "\n---finished!\n";
close $fh;
