#!/bin/perl -w

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;

my $file = $ARGV[0];
my $archive = new Bio::Seq();
my $searchio = Bio::SearchIO->new( -format => 'blast', -file => $file);
my $IO =new Bio::SeqIO(-format =>'GenBank',-file => '>salida.gbff');
my $feature;
my $location;

while ( my $result = $searchio->next_result() ) {
	#print $result->query_name."\t".$result->query_description."\n";
	$location = $result->query_description;	
	#($start) = ($location =~ /([0-9]+)/);
#	($end) = ($location =~ m/(\d+)s*---/);
	while( my $hit = $result->next_hit ) {
    		my $aux= $hit->description;
		my $hsp = $hit->next_hsp;
		#print $hit->frac_aligned_query()."\n";
		if($aux=~ m/hypothetical/ig){
			next;
		}else{
			if($hit->frac_aligned_query() >= 0.95 && $hit->significance <= 10**(-5) ) {
				#Description en /product
				$feature = new Bio::SeqFeature::Generic(
					-start => 1,#($location =~ /([0-9]+)/),
					-end =>3, # ($location =~ m/(\d+)s*---/),
					-primary_tag => 'gene');
				
				$archive->add_SeqFeature($feature);

				$feature = new Bio::SeqFeature::Generic(
					-start => 1,
					-end => 4,
					-primary_tag => 'CDS',
					-tag => {product => $hit->description,
						note => $hit->significance,
						db_xref => $hit->accession,
						translation => $hsp->hit_string});
				$archive->add_SeqFeature($feature);

			
			}else{
				if($hit->frac_aligned_query() <=0.95 && $hit->significance <= 10**(-5)){
					#Description a /product + putative
					#/note descripcion de posiciones
					$feature = new Bio::SeqFeature::Generic(
						-start => 2,
						-end => 2,
						-primary_tag => 'gene');
					$archive->add_SeqFeature($feature);

					$feature = new Bio::SeqFeature::Generic(
						-start =>2,
						-end => 4,
						-primary_tag => 'CDS',
						-tag =>{product => $hit->description.", putative",
							note =>	' ',
							db_xref => $hit->accession,
							translation => $hsp->hit_string});
					$archive->add_SeqFeature($feature);
					
				}else{
					if($hit->frac_aligned_query() <=0.95 && $hit->significance >= 10**(-5)) {
						#/product = "hypotetical protein"
						#/note = no significative hits found"
						$feature = new Bio::SeqFeature::Generic(
							-start =>3,
							-end => 5,
							-primary_tag => 'gene');
						$archive->add_SeqFeature($feature);

						$feature =new Bio::SeqFeature::Generic(
							-start => 3,
							-end => 6,
							-primary_tag => 'CDS',
							-tag => {product => 'hypotetical protein',
								note => 'no significative hits found',
								db_xref => $hit->accession,
								translation => $hsp->hit_string});
						$archive->add_SeqFeature($feature);
					}
				}
			}
			$IO->write_seq($archive);

		}


	}
#		print $hit->description."\n";
#	while( my $hsp = $hit->next_hsp ) {
        	 # process the Bio::Search::HSP::HSPI object
#        }
}
