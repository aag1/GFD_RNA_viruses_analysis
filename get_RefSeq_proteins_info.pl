use strict;
use warnings;
use Bio::DB::GenPept;
use Bio::Annotation::Collection;
use Bio::SeqIO;


my @IDS = @ARGV;


my $db = Bio::DB::GenPept->new();


print "RefSeq_prot_id\tRefSeq_prot_len\tRefSeq_prot_desc\tRefSeq_prot_taxo\n";

while (@IDS) {

    my @ids = splice @IDS, 0, 50;

	my $seqio = $db->get_Stream_by_id( [@ids] );

	while( my $entry  =  $seqio->next_seq ) {

		my $id = $entry->id . '.' . $entry->version;
		my $len = $entry->length;
		my $desc = $entry->desc;
		my $taxo = join ";", reverse($entry->species->classification);
		print $id . "\t" . $len . "\t" . $desc . "\t" . $taxo . "\n";

	}

}
