#!/usr/bin/env perl

package MyApp::AnnotateAssemblyBlast;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SeqFeature::Generic;
use Bio::Tools::Run::StandAloneBlastPlus;
use Data::Printer;
use File::Temp;
use namespace::autoclean;


    command_long_description q[Align a fasta sequence to a genbank reference using HMMER or CLUSTALW and transfer the feature annotations.];

    has_file 'input_file' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 1,
        must_exist    => 1,
        documentation => q[FASTA file to be annotated],
    );

    has_file 'reference_file' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 1,
        must_exist    => 1,
        documentation => q[Reference file],
    );

    has_file 'output_genbank' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 1,
        documentation => 'Output Genbank',
    );

    has_file 'alignment_output_file' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 0,
        documentation => 'Output HMMER alignment file (optional)',
    );


    has 'reference' => (
        is         => 'ro',
        isa        => 'Object',
        lazy_build => 1
    );

    has 'input_fasta' => (
        is         => 'ro',
        isa        => 'Object',
        lazy_build => 1
    );


    sub get_seqobj {
        my ( $self, $file ) = @_;
        my $in = Bio::SeqIO->new( -file => $file );
        my $seqobj = $in->next_seq;
        my $seq = $seqobj->seq;
        #remove non IUPAC bases
        $seq =~ s/[^ACGTURYSWKMBDHVN]+//gi;
        $seqobj->seq($seq);
        return $seqobj;
    }

    sub _build_reference {
        my ($self) = @_;
        return $self->get_seqobj( $self->reference_file->stringify );
    }

    sub _build_input_fasta {
        my ($self) = @_;
        return $self->get_seqobj( $self->input_file->stringify );
    }

    sub get_bl2seq_alignment {
        my ($self) = @_;
        my $fh = File::Temp->new(UNLINK => 1);
        my $ref_fh = File::Temp->new(UNLINK => 1);
        my $fname = $fh->filename;
        my $ref_name = $ref_fh->filename;
        my $factory = Bio::Tools::Run::StandAloneBlastPlus->new();

        my $out_ref = Bio::SeqIO->new(-file =>  ">" . $ref_name, -format => 'fasta' );
        $out_ref->write_seq($self->reference);

        #my $string_gen = String::Random->new();
        #my $tmp_name = $string_gen->randregex('\w\w\w\w\w');
        #my $tmp_fas = $tmp_name + '.fas';

        my $bl2seq_report = $factory->bl2seq(
            -method      => 'blastn',
            -query       => $ref_name,
            -subject     => $self->input_file->stringify,
            -outfile     => $fname,
            -method_args => [
                '-task'            => 'dc-megablast',
                '-xdrop_gap_final' => '500'
            ],
        );

        #$aln->map_chars('\.','-');
        # Note that report is a Bio::SearchIO object

        # Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
        my $str =
          Bio::AlignIO->new( -file => $fname, '-format' => 'bl2seq' );
    
        my @aln_array;
        while (my $aln = $str->next_aln){
            push @aln_array, $aln;
        }
        return \@aln_array;
    }

    sub transfer_simple_location {
        my ( $self, $aln, $location_object ) = @_;
        my $seq_object = $self->reference;
        my $ref_seq = $aln->get_seq_by_pos(1);
        my $aligned_seq = $aln->get_seq_by_pos(2);
        
        # check if the feature location is present in alignment
        my $broken_start;
        my $broken_end;
        if ($location_object->end <= $ref_seq->start || $location_object->start  >= $ref_seq->end ){
            return 0;
        }
        
        if ($location_object->start < $ref_seq->start){
            $broken_start = $ref_seq->start;
        }
        if ($location_object->end > $ref_seq->end){
            $broken_end = $ref_seq->end;
        }
        
        my $start_col;
        if ($broken_start) {
            $start_col = $aln->column_from_residue_number( $seq_object->id,
                $broken_start );

        }
        else {
            $start_col = $aln->column_from_residue_number( $seq_object->id,
                $location_object->start );

        }

        my $end_col;
        if ($broken_end) {
            $end_col = $aln->column_from_residue_number( $seq_object->id,
                $broken_end );
        }
        else {
            $end_col = $aln->column_from_residue_number( $seq_object->id,
                $location_object->end );
        }

        my $loc_start = $aligned_seq->location_from_column($start_col);
        my $loc_end   = $aligned_seq->location_from_column($end_col);

        #my $aln_slice = $aln->slice($start_col, $end_col);
        #say "\tidentity: ".$aln_slice->percentage_identity;
        if ( $loc_start && $loc_end ) {
            my $new_loc = $loc_start->union($loc_end);
            $new_loc->strand($aligned_seq->strand);

            #skip features small broken features
            if ( $new_loc->end - $new_loc->start  <= 2  && $loc_start->location_type =~ /IN\-BETWEEN/ && $loc_start->location_type =~ /IN\-BETWEEN/ ){
                return 0;
            }
            if (   $loc_start->location_type =~ /IN\-BETWEEN/
                || $loc_end->location_type =~ /IN\-BETWEEN/ || $broken_start || $broken_end )
            {
                # Create a fuzzy location object
                my $fuzzy_loc = Bio::Location::Fuzzy->new(-start => $new_loc->start, -end => $new_loc->end, -strand => $new_loc->strand);
                if ( $loc_start->location_type =~ /IN\-BETWEEN/ || $broken_start) {
                    $fuzzy_loc->start( '<' . $new_loc->start );
                }
                if ( $loc_end->location_type =~ /IN\-BETWEEN/ || $broken_end ) {
                    $fuzzy_loc->end( '>' . $new_loc->end );
                }
                $new_loc = $fuzzy_loc;
            }
            if ($new_loc->end > $self->input_fasta->length){
                $new_loc->end($self->input_fasta->length);
            }
            return $new_loc;
        }
        else{
            return 0;
        }
    }

    sub transfer_split_location {
        my ( $self, $aln, $location_object ) = @_;
        my $count;
        my $splitlocation = Bio::Location::Split->new();
        for my $sublocation ($location_object->sub_Location){
            my $mod_loc = $self->transfer_simple_location($aln, $sublocation);
            if ($mod_loc) {
                $count++;
                $splitlocation->add_sub_Location($mod_loc);
               
            }
        }
        
        if ($count){
            return($splitlocation);
        }else{
            return 0;
        }
    }

    sub run {
        my ($self)     = @_;
        my $aligner = 'BL2SEQ';
        my $aln_array = $self->get_bl2seq_alignment;
        #p $aln;
        # Check if HMMER could not properly align the sequence
        my $seq_object =  $self->reference;
        my $new_seq = $self->input_fasta;
        for my $feat_object ( $seq_object->get_SeqFeatures ) {
            for my $aln ( @{$aln_array} ) {
                my $mod_location;
                if ( $feat_object->location->isa('Bio::Location::Simple') ) {
                    $mod_location = $self->transfer_simple_location( $aln,
                        $feat_object->location );
                }
                elsif ( $feat_object->location->isa('Bio::Location::Split') ) {
                    $mod_location = $self->transfer_split_location( $aln,
                        $feat_object->location );
                }

                # If there  is a location for the feature in the FASTA sequence
                if ($mod_location) {
                    my $new_feat = Bio::SeqFeature::Generic->new();
                    $new_feat->location($mod_location);
                    $new_feat->primary_tag( $feat_object->primary_tag );
                    for my $tag ( $feat_object->get_all_tags ) {
                        for my $value ( $feat_object->get_tag_values($tag) ) {
                            $new_feat->add_tag_value( $tag, $value );
                        }
                    }

                    $new_seq->add_SeqFeature($new_feat);
                }

            }
        }

        # Add tranlation to the new sequences
        for my $feat ( $new_seq->get_SeqFeatures ) {
            if ($feat->primary_tag =~ /CDS/ig){
                $feat->add_tag_value('translation',$feat->spliced_seq->translate->seq)
            }
        }
 
        my $comment = Bio::Annotation::Comment->new;
        $comment->text("Aligned with ".$aligner);
        my $coll = Bio::Annotation::Collection->new;
        $coll->add_Annotation( 'comment', $comment );
        $new_seq->annotation($coll);
        my $out = Bio::SeqIO->new(
            -file   => '>' . $self->output_genbank->stringify,
            -format => 'Genbank'
        );
        $out->write_seq($new_seq);
    }

 __PACKAGE__->meta->make_immutable;
