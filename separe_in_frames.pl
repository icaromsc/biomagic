#!/bin/perl -w
use Bio::SeqIO;

$file = $ARGV[0];
$name = $file;
$name =~ s/\..*//; 

$seqio_obj = Bio::SeqIO->new(-file => $file, 
                             -format => "fasta" );
$seq_objP = $seqio_obj->next_seq;
$seq_objN = $seq_objP->revcom; 
print"\ncreating frames...\n";
for ($i = 0; $i < 3; $i++) {
	print "\n finding ORFs from frame ".$i."...";
	$prot_obj = $seq_objP->translate(-frame => $i);
	$temp_seq_obj = Bio::Seq->new(-seq=>$prot_obj->seq,
                         -display_id => $seq_objP->id,
                         -desc =>$seq_objP->desc." 5'3 frame ".$i." translated",
                         -alphabet => "protein");
	
	$aux =  $temp_seq_obj->seq;
	$largo = length ($seq_objP->seq);
	$seqio_obj2 = Bio::SeqIO->new(-file => '>'.$name.'_p_strand_fr_'.$i.'.fasta',-format => 'fasta' );
	$v=0;
	$r=0;	
		while ($aux =~ m/(\w{80,})/g ){
			$v++;
			if($i == 0){
				$inicio = ($-[0])*3;
			   	$fin = ($+[0])*3;
			   	#print $aux;
		   	}else{
				if($i == 1){ 
			   		$inicio = (($-[0])*3)-1;
					$fin = (($+[0])*3)-1;
			   	}else{
					if($i == 2){
						$inicio = (($-[0])*3)-2;
					   	$fin = (($+[0])*3)-2;
				   	}
			   	}
		   	}	   
	  		$orf_obj = Bio::Seq->new(-seq=>$1,
			-display_id => "seq".$v,
			-desc=>$inicio." ".$fin." --- plus strand frame ".$i);
			$seqio_obj2->write_seq($orf_obj);
	  	}
	print "\n  --finded ".$v." ORFs on plus strand\n"; 
   		
	$prot_obj = $seq_objN->translate(-frame => $i);
 	$temp_seq_obj = Bio::Seq->new(-seq=>$prot_obj->seq,
                         -display_id => $seq_objP->id,
                         -desc =>$seq_objP->desc." 3'5 frame ".$i." translated",
                         -alphabet => "protein");
	$aux2 = $temp_seq_obj->seq;
	$seqio_obj2 = Bio::SeqIO->new(-file => '>'.$name.'_n_strand_fr_'.$i.'.fasta',-format => 'fasta' );  
	
	while ($aux2 =~ m/(\w{80,})/g ){
		$r++;
		if($i == 0){
			$inicio = $largo - (($-[0])*3);
			$fin = $largo - (($+[0])*3);
		}else{
			if($i == 1){
				$inicio = $largo - ((($-[0])*3)+1);
			        $fin = $largo - ((($+[0])*3)+1);		
			}else{
				if($i == 2){
					$inicio = $largo - ((($-[0])*3)+2);
					$fin = $largo - ((($+[0])*3)+2);
				}
			}
		}
		$orf_obj = Bio::Seq->new(-seq=>$1,
		-display_id => "seq".$r,
		-desc=>$inicio." ".$fin." --- minus strand frame ".$i);
		$seqio_obj2->write_seq($orf_obj);
	}
	print "  --finded ".$r." ORFs on minus strand\n";
	print "\nsaving .fasta file...\n"; 
}
print "\nfinished!!!\n";
