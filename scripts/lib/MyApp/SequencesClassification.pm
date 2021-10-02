#!/usr/bin/env perl

use warnings;
use strict;
package MyApp::SequencesClassification;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
with 'MyApp::Role::JSON';
use MooseX::FileAttribute;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use namespace::autoclean;
use Data::Printer;
use Bio::SeqIO;
use File::Basename;
use File::Path 'make_path';

    ###################################################

    has_directory 'output_dir' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[Very important option!],
    );

    has_file 'genbank_file' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[Very important option!],
    );

    has 'genbank' => (
        is         => 'rw',
        isa        => 'Bio::Seq',
        lazy_build => 1,
    );

    has 'location' => (
        is         => 'rw',
        isa        => 'HashRef',
        lazy_build => 1
    );

    has 'sample_id' => (
        is         => 'rw',
        isa        => 'Str',
        lazy_build => 1,
    );

    has 'file_id' => (
        is         => 'rw',
        isa        => 'Str',
        lazy_build => 1,
    );

    has 'genes_correct_order' => (
        is      => 'ro',
        isa     => 'ArrayRef',
        default => sub { return [ 'psi', 'gag', 'pol', 'vif', 'vpr', 'tat', 'rev', 'vpu', 'env', 'nef' ] }
    );

    has 'gene_range_allowed' => (
        is      => 'ro',
        isa     => 'HashRef[ArrayRef]',
        default => sub {
            return {
                'gag' => [ 473 * 3, 532 * 3 ],
                'pol' => [ 885 * 3, 923 * 3 ],
                'env' => [ 810 * 3, 913 * 3 ],
                'nef' => [ 170 * 3, 232 * 3 ],
                'vif' => [ 188 * 3, 199 * 3 ],
                'rev' => [ 85 * 3,  130 * 3 ],
                'vpu' => [ 65 * 3,  94 * 3 ],
                'tat' => [ 84 * 3,  105 * 3 ],
                'vpr' => [ 84 * 3,  103 * 3 ]
            };
        }
    );


#--------------------------------------
# Moose attribute builders
#--------------------------------------
    sub _build_genbank {
        my ( $self ) = @_;
        my $in  = Bio::SeqIO->new( -file => $self->genbank_file->stringify, -format => 'genbank' );
        my $seq = $in->next_seq;
        return $seq;
    }


    sub _build_location {
        my ( $self ) = @_;
        my %hash_location;
        for my $feat_obj ( $self->genbank->get_SeqFeatures ) {
            if ( $feat_obj->primary_tag eq "gene" || $feat_obj->primary_tag eq "regulatory" ) {
                my @feat_name = $feat_obj->get_tag_values( 'standard_name' );
                my $gene_name;
                foreach my $normalized_gene_name ( @{ $self->genes_correct_order } ) {
                    if ( $feat_name[ 0 ] =~ m/\Q$normalized_gene_name\E/i ) {
                        $gene_name = $normalized_gene_name;
                    }
                }
                push @{ $hash_location{ $self->genbank->id }{ $gene_name } }, $feat_obj->location;
            }
        }

        $self->_merge_close_location( \%hash_location );
        return \%hash_location;
    }


    sub _build_sample_id {
        my ($self) = @_;

        my $modified_genbank_id = $self->genbank->id;
        $modified_genbank_id =~ s/_S.*//g;

        return($modified_genbank_id)

    }


    sub _build_file_id {
        my ($self) = @_;

        my $file_id = $self->genbank->id;

        return($file_id)

    }


#--------------------------------------
# Moose Private Methods
#--------------------------------------
    sub _check_location {
        my ( $self, $current_location, $next_location ) = @_;

        if ( ( $current_location->start < $next_location->start ) ) {

            return ( 1 );
        }
        elsif ( ( $current_location->start == $next_location->start ) ) {

            return ( 1 );
        }
        else {
            return ( 0 );
        }

    }


    sub _can_be_merged {
        my ( $self, $current_location, $next_location ) = @_;
        if ( $next_location->start - $current_location->end < 10 && $current_location->strand == $next_location->strand ) {
            return 1;
        }
        else {
            return 0;
        }
    }


    sub _merge_close_location {
        my ( $self, $hash_location ) = @_;
        for my $seq ( keys %{ $hash_location } ) {
            for my $gene ( sort { $a cmp $b } keys %{ $hash_location->{ $seq } } ) {

                # Get list of locations
                my @locations = @{ $hash_location->{ $seq }{ $gene } };

                # Let bioperl merge all overlaping rangesqq
                my @maximum_merged_locations = Bio::Location::Simple->disconnected_ranges( @locations );
                @maximum_merged_locations = sort { $a->start <=> $b->start } @maximum_merged_locations;

                my @new_locations;
                my $merged_location;
                while ( my $current_location = shift @maximum_merged_locations ) {
                    if ( @maximum_merged_locations ) {
                        unless ( $merged_location ) {
                            my $next_location = $maximum_merged_locations[ 0 ];
                            if ( $self->_can_be_merged( $current_location, $next_location ) ) {
                                $merged_location = Bio::Location::Simple->union( $current_location, $next_location );
                                shift @maximum_merged_locations;
                            }
                            else {
                                push @new_locations, $current_location;
                            }
                        }
                        else {
                            if ( $self->_can_be_merged( $merged_location, $current_location ) ) {
                                $merged_location->union( $current_location );
                                shift @maximum_merged_locations;
                            }
                            else {
                                push @new_locations, $merged_location;
                                $merged_location = undef;
                            }
                        }
                    }
                    else {
                        if ( $merged_location ) {
                            if ( $self->_can_be_merged( $merged_location, $current_location ) ) {
                                $merged_location->union( $current_location );
                            }
                            else {
                                push @new_locations, $current_location;
                            }
                            push @new_locations, $merged_location;
                            $merged_location = undef;
                        }
                        else {
                            push @new_locations, $current_location;
                        }
                    }
                }
                push @new_locations, $merged_location if $merged_location;
                $hash_location->{ $seq }{ $gene } = \@new_locations;
            }
        }


    }


    sub _get_features {
        my ( $self, $primary_tag, $tag, $value ) = @_;
        my @features;
        foreach my $feat_obj ( $self->genbank->get_SeqFeatures ) {
            if ( $feat_obj->primary_tag eq $primary_tag ) {
                my @values = $feat_obj->get_tag_values( $tag );

                if ( $values[ 0 ] =~ /\Q$value\E/i ) {
                    push @features, $feat_obj;
                }
            }
        }

        # sort feature by start location
        @features = sort { $a->location->start <=> $b->location->start } @features;
        return @features;
    }


    sub _create_feature_from_location_object {
        my ( $self, $location ) = @_;
        my $final_location;
        if ( ref $location eq 'ARRAY' ) {
            $final_location = Bio::Location::Split->new();
            $final_location->add_sub_Location( @{ $location } );
        }
        elsif ( $location->isa( 'Bio::Location::Simple' ) ) {
            $final_location = $location;
        }
        my $feature = Bio::SeqFeature::Generic->new( -location => $final_location );
        my $seq     = Bio::PrimarySeq->new( -id => $self->genbank->id, -seq => $self->genbank->seq );
        $feature->attach_seq( $seq );
        return $feature;
    }


    sub _check_genes {
        my ($self, $productive_genes, $unproductive_genes, $MSD_seq, $psi_status) = @_;

        for my $feat_obj ( $self->genbank->get_SeqFeatures ) {

            if ( $feat_obj->has_tag( 'standard_name' ) ) {
                my @cur_feat_name = $feat_obj->get_tag_values( 'standard_name' );

                if ( $feat_obj->primary_tag eq "misc_feature" and $cur_feat_name[ 0 ] eq "MSD site" ) {

                    ${$MSD_seq} = $feat_obj->seq->seq;

                    if ( ${$MSD_seq} eq "GT" ) {
                        ${$psi_status} = "intact";
                        my $text = '>>> psi_status:  intact';
                        p $text;
                    }
                }

            }
        }

        for my $seq_id ( sort { $a cmp $b } keys %{ $self->location } ) {
            for my $gene ( sort { $a cmp $b } keys %{ $self->location->{ $seq_id } } ) {
                my $bio_loc;

                if ( $gene eq "tat" or $gene eq "rev" ) {
                    $bio_loc = $self->get_split_location( 'CDS', 'standard_name', $gene );
                }
                elsif ( $gene eq "psi" ) {
                    next;
                }
                else {
                    $bio_loc = $self->location->{ $seq_id }{ $gene }->[ 0 ];
                }

                my $cds_intact;
                if ( $bio_loc ) {
                    $cds_intact = $self->check_intact( $bio_loc, $gene );
                }

                if ( $cds_intact ) {
                    push( @{$productive_genes}, $gene );

                    #p $translation;
                    say ">>> " . $gene . ": intact";
                }
                else {
                    #p $translation;
                    say ">>> " . $gene . ": not intact";
                    push( @{$unproductive_genes}, $gene );
                }
            }
        }

    }

#--------------------------------------
# Moose Public Methods
#--------------------------------------
    sub check_psi {
        my ($self, $json_report_hash) = @_;
        my $found_psi = 0;

        for my $feat_obj ($self->genbank->get_SeqFeatures) {

            if ( $feat_obj->primary_tag eq "regulatory" ) {
               $found_psi++; 
            }
        }

        $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"input"} = $self->genbank_file->stringify;
        if ($found_psi == 0) {
            $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_psi"} = "No";
        } else {
            $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_psi"} = "Yes";
        }

    }


    sub check_inversions {
        my ($self, $json_report_hash) = @_;
        my @inverted_loci;
        my %inverted_loci_type;
        my %correct_loci_type;

        for my $feat_obj ($self->genbank->get_SeqFeatures) {
            if ($feat_obj->primary_tag eq "CDS" or $feat_obj->primary_tag eq "LTR") {
                my @gene_id     = $feat_obj->get_tag_values('standard_name');
                my $gene_strand = $feat_obj->location->strand;
                if ($gene_strand == -1) {
                    push @inverted_loci, $gene_id[0];
                    push @{$inverted_loci_type{$feat_obj->primary_tag}}, $feat_obj;
                }else{
                    push @{$correct_loci_type{$feat_obj->primary_tag}}, $feat_obj;
                }
            }
        }

        if (scalar @inverted_loci != 0) {
            if ($inverted_loci_type{'LTR'} && !$inverted_loci_type{'CDS'} && $correct_loci_type{'LTR'} ){
                my @reverse_sorted_ltr = sort {$b->start <=> $a->end}  @{$inverted_loci_type{'LTR'}};
                if ($self->genbank->length - $reverse_sorted_ltr[0]->start <= 150 ){
                    $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_invertions"} = "No";
                }
                else {
                    $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"inverted_genes"} = \@inverted_loci;
                    $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_invertions"} = "Yes";
                }
            }
            else {
                $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"inverted_genes"} = \@inverted_loci;
                $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_invertions"} = "Yes";
            }      
        } else {
            $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_invertions"} = "No";
        }
    }


    sub check_gene_order {
        my ( $self, $json_report_hash ) = @_;

        for my $seq_id ( keys %{ $self->location } ) {
            #say $seq_id;
            my @aux_genes = @{ $self->genes_correct_order };

            my @clean_genes;
            foreach my $gene ( @aux_genes ) {
                if ( $self->location->{ $seq_id }{ $gene } ) {
                    push @clean_genes, $gene;
                }
            }
            
            my $not_correct_order_count = 0;
            @aux_genes = @clean_genes;
            while ( my $gene = shift @aux_genes ) {
                if ( $self->location->{ $seq_id }{ $gene } ) {

                    if ( @aux_genes ) {
                        my $next_gene = $aux_genes[ 0 ];
                        #say "\tCurrent $gene => next $next_gene";
                        for my $current_gene_location ( @{ $self->location->{ $seq_id }{ $gene } } ) {

                            for my $next_gene_location ( @{ $self->location->{ $seq_id }{ $next_gene } } ) {
                                #say "\t\t" . $gene . " " . $current_gene_location->end . " " . $next_gene . " " . $next_gene_location->start;
                                my $verif_location = $self->_check_location( $current_gene_location, $next_gene_location );
                                if ($verif_location == 1) {
                                    $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"genes_order"}{$gene.'_'.$current_gene_location->start.'->'.$next_gene.'_'.$next_gene_location->start} = "Correct";
                                } else {
                                    $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"genes_order"}{$gene.'_'.$current_gene_location->start.'->'.$next_gene.'_'.$next_gene_location->start} = "Not_correct";
                                    $not_correct_order_count++;
                                }
                            }
                        }
                    }
                }
            }

            if ($not_correct_order_count > 0) {
                $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_order"} = "No";
            } else {
                $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"has_order"} = "Yes";
            }

        }
    }  


    sub check_duplication {
        my ($self, $json_report_hash) = @_;
        my $duplication = 0;

        for my $seq_id (keys %{$self->location}){
            for my $gene_id (keys %{$self->location->{$seq_id}}){

                my $gene_count = scalar @{$self->location->{$seq_id}{$gene_id}};
                $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{"gene_count"}{$gene_id} = $gene_count;

                if ($gene_count > 1 and $duplication == 0) {
                    $duplication = 1;
                }
            }
        }

        if ($duplication == 1) {
            $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{'has_duplication'} = 'Yes';
        } else {
            $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{'has_duplication'} = 'No';
        }

    }


    sub check_ltr {
        my ($self, $json_report_hash) = @_;
        my $has_ltr = 0;

        for my $feat_obj ( $self->genbank->get_SeqFeatures ) {

            if ( $feat_obj->primary_tag eq "LTR" ) {
                $has_ltr++;
                my @feat_name = $feat_obj->get_tag_values('standard_name');

                if ($feat_name[0] =~ /^5/) {
                    push @{ $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{'ltr'} }, '5_prime_start_'.$feat_obj->location->start.'_end_'.$feat_obj->location->end.'_size_'.(($feat_obj->location->end-$feat_obj->location->start) + 1);
                } elsif ($feat_name[0] =~ /^3/) {
                    push @{ $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{'ltr'} }, '3_prime_start_'.$feat_obj->location->start.'_end_'.$feat_obj->location->end.'_size_'.(($feat_obj->location->end-$feat_obj->location->start) + 1);
                }
            }
        }

        if ($has_ltr == 0) {
            $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{'has_ltr'} = 'No';
        } else {
            $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step}{'has_ltr'} = 'Yes';
        }
    }


    sub count_stop_codons {
        my ( $self, $translation ) = @_;
        my $stop_codons_count = 0;

        $stop_codons_count = grep { /\*/ } split( "", $translation );

        return $stop_codons_count;
    }


    sub feature_reannotation {
        my ( $self, $new_translation, $gene ) = @_;
        my ( $cds, $cds_error )  = $self->_get_features( 'CDS',  'standard_name', $gene );
        my ( $gene_feature, $gene_feature_error ) = $self->_get_features( 'gene', 'standard_name', $gene );
        
        if ($cds_error or $gene_feature_error) {
            return 0;
        }

        if ( $cds->location->isa( 'Bio::Location::Split' ) ) {
            my @sublocations = $cds->location->sub_Location;
            @sublocations = sort { $a->end <=> $b->end } @sublocations;
            my $first_exon_nt = ( ( $sublocations[ 0 ]->end - $sublocations[ 0 ]->start ) + 1 );
            $sublocations[ 1 ]->end( ( $sublocations[ 1 ]->start + ( ( ( length( $new_translation ) * 3 ) - $first_exon_nt ) ) - 1 ) );
            $gene_feature->location->end($cds->end);
        }
        elsif ( $cds->location->isa( 'Bio::Location::Simple' ) ) {
            $cds->location->end( ( $cds->location->start + ( length( $new_translation ) * 3 ) ) - 1 );
            $gene_feature->location->end( ( $gene_feature->location->start + ( length( $new_translation ) * 3 ) ) - 1 );
        }
    }


    sub check_intact {

        my ( $self, $bio_loc, $gene ) = @_;
        my $match_unknown_aa;

        my $seq_feature = Bio::SeqFeature::Generic->new( -location => $bio_loc );
        my $seq         = Bio::PrimarySeq->new( -id => $self->genbank->id, -seq => $self->genbank->seq );
        $seq_feature->attach_seq( $seq );

        my $translation = $seq_feature->spliced_seq->translate->seq;
        #p $translation;
        say $translation;

        #Checking existence of unknown amino acid
        $match_unknown_aa = index( uc $translation, "X" );

        #If there is one unkown aminoacid or more than two stop codons on the gene
        if ( $match_unknown_aa != -1) {
            return 0;
        }

        my $sub_translation;
        my $sub_count;
        if ( $translation =~ /^(\S+)\S{20}$/ ) {
            $sub_translation = $1;
            $sub_count       = $self->count_stop_codons( $sub_translation );
        }
        else {
            return 0;
        }

        if ( $sub_count >= 1 ) {
            return 0;
        }

        my @split_translation = split( "", $translation );
        if ( $split_translation[ -1 ] eq "*" ) {


            my @features                   = $self->_get_features( 'CDS', 'standard_name', $gene );
            my $last_feature               = $features[ $#features ];


            if ( $last_feature->location->isa( 'Bio::Location::Fuzzy' ) ) {
                return 0 if $last_feature->location->start_pos_type eq 'BEFORE';

            }
            elsif ( $last_feature->location->isa( 'Bio::Location::Split' ) ) {
                my @sublocations = $last_feature->location->sub_Location;
                @sublocations = sort { $a->end <=> $b->end } @sublocations;

                foreach my $sublocation ( @sublocations ) {
                    if ( $sublocation->isa( 'Bio::Location::Fuzzy' ) ) {
                        return 0 if $sublocation->start_pos_type eq 'BEFORE';

                    }

                }

            }


            my $new_translation;
            if ($translation =~ /(\S+?\*)/) {
                $new_translation = $1;
            }

            my $has_passed = $self->has_expected_genes_size( ( length( $new_translation ) * 3 ), $gene );
            if ( $has_passed ) {
                # Reannotation
                $self->feature_reannotation($new_translation, $gene);
                return 1;
            }
            else {
                return 0;
            }

        }
        elsif ( $split_translation[ -1 ] ne "*" ) {
            #p $bio_loc;
            my $extended_translation;
            my @current_gene_allowed_range = @{ $self->gene_range_allowed->{ $gene } };
            my @features                   = $self->_get_features( 'CDS', 'standard_name', $gene );
            my $last_feature               = $features[ $#features ];

            # check if feature has broken end;
            if ( $last_feature->location->isa( 'Bio::Location::Split' ) ) {
                my @sublocations = $last_feature->location->sub_Location;
                @sublocations = sort { $b->end <=> $a->end } @sublocations;
                my $last_sublocation = $sublocations[ 0 ];
                if ( $last_sublocation->isa( 'Bio::Location::Fuzzy' ) ) {
                    return 0 if $last_sublocation->end_pos_type eq 'AFTER';
                }
            }
            elsif ( $last_feature->location->isa( 'Bio::Location::Fuzzy' ) ) {
                return 0 if $last_feature->location->end_pos_type eq 'AFTER';
            }

            # Search if there is a stop codon within the maximum expected range
            if ( $bio_loc->isa( "Bio::Location::Simple" ) ) {
                my $max_extension_end = (($bio_loc->start + $current_gene_allowed_range[ 1 ]) - 1) ;
                if ( $max_extension_end > $self->genbank->length ) {
                    $max_extension_end = $self->genbank->length;
                }

                $extended_translation = $self->genbank->trunc( $bio_loc->start, $max_extension_end )->translate->seq;
            }
            elsif ( $bio_loc->isa( "Bio::Location::Split" ) ) { 
                my $total_residues;
                my @sublocations = $bio_loc->sub_Location;
                my @extend_sublocations;
                while ( my $sub_location = shift @sublocations ) {
                    my $residues_count = (( $sub_location->end - $sub_location->start ) + 1);
                    $total_residues += $residues_count;
                    if ( @sublocations ) {
                        push @extend_sublocations, $sub_location;
                    }
                    else {
                        if ( $total_residues <= $current_gene_allowed_range[ 1 ] ) {

                            my $new_end = $sub_location->end + ( $current_gene_allowed_range[ 1 ] - $total_residues );
                            if ($new_end > $self->genbank->length) {
                                $new_end = $self->genbank->length;
                            }

                            $sub_location->end( $new_end );
                            push @extend_sublocations, $sub_location;
                        }
                        else {
                            return 0;
                        }
                    }
                }
                my $split_feature = $self->_create_feature_from_location_object( \@extend_sublocations );
                $extended_translation = $split_feature->spliced_seq->translate->seq;
            }
            my $debug = "Extended Translation:\n$extended_translation";
            #p $debug;
            if ( $extended_translation =~ /(\S+?\*)/ ) {
                my $extended_translation_up_to_first_stop_codon = $1;
                my $has_passed = $self->has_expected_genes_size( ( length( $extended_translation_up_to_first_stop_codon ) * 3 ), $gene );
                if ( $has_passed ) {
                    # Reannotation
                    $self->feature_reannotation($extended_translation_up_to_first_stop_codon, $gene);

                    return $extended_translation;
                }
                else {
                    return 0;
                }
            }
            else {
                return 0;
            }
        }
    }


    sub get_split_location {
        my ($self, $primary_tag, $tag, $value) = @_;
        my @parts;
        my $split_location;
        my @features = $self->_get_features($primary_tag, $tag, $value);
        foreach my $feat_obj ( @features ) {
            my $location_obj = $feat_obj->location;
            if ( $location_obj->isa( 'Bio::Location::Simple' ) ) {
                push @parts, $location_obj;
            }
            elsif ( $location_obj->isa( 'Bio::Location::Split' ) ) {
                $split_location = $location_obj;
            }
            elsif ( $location_obj->isa( 'Bio::Location::Fuzzy' ) ) {
                push @parts, $location_obj;
            }
            else {
                die "We have a problem!";
            }
        }
        if (@parts){
            $split_location = Bio::Location::Split->new();
            $split_location->add_sub_Location(@parts);
        }
        return $split_location;
    }


    sub sequence_classification {
        my ( $self, $json_report_hash ) = @_;

        my $json_report_current = $json_report_hash->{'passed'}{ $self->sample_id }{$self->file_id}{'STEP'.$self->current_step};
        $json_report_current->{ 'output' } = $self->json_file->stringify;

        my @productive_genes;
        my @unproductive_genes;
        my $MSD_seq;
        my $psi_status = "mutation";
        my $classification;

        $self->_check_genes( \@productive_genes, \@unproductive_genes, \$MSD_seq, \$psi_status );

        $json_report_current->{ 'productive_genes' }   = join( ", ", @productive_genes );
        $json_report_current->{ 'unproductive_genes' } = join( ", ", @unproductive_genes );
        $json_report_current->{ 'psi_status' }         = $psi_status;

        # description
        my $non_functional_genes = 'non-functional genes:' . join ",", @unproductive_genes;
        $self->genbank->description( $non_functional_genes );

        my $genes_count = 0;
        foreach my $gene ( keys %{ $json_report_current->{ 'gene_count' } } ) {
            $genes_count += $json_report_current->{ 'gene_count' }{ $gene };
        }

        if ( $json_report_current->{ 'has_order' } eq 'No' ) {
            $classification = 'duplications_or_inversion';
        }
        elsif ( $json_report_current->{ 'has_invertions' } eq 'Yes' and $json_report_current->{ 'has_ltr' } eq 'Yes') {
            $classification = 'duplications_or_inversion';
        }
        elsif ( $json_report_current->{ 'has_duplication' } eq 'Yes' ) {
            $classification = 'duplications_or_inversion';
        }
        # missing_ltr_or_psi classification
        elsif ( ( $json_report_current->{ 'has_ltr' } eq "No" ) ) {
            $classification = 'missing_ltr_or_psi'
        }
        else {
            my $has_five_prime_ltr  = 0;
            my $has_three_prime_ltr = 0;

            my @found_ltrs      = @{ $json_report_current->{ 'ltr' } };
            my @ltr_five_prime  = grep { /^5/ } @found_ltrs;
            my @ltr_three_prime = grep { /^3/ } @found_ltrs;

            if ( scalar @ltr_five_prime > 0 ) {
                foreach my $mapped_ltr ( @ltr_five_prime ) {
                    my $ltr_start = $mapped_ltr;
                    my $ltr_end   = $mapped_ltr;
                    my $ltr_size  = $mapped_ltr;

                    if ( $ltr_start =~ /.*start_(\d+)_.*/ ) {
                        $ltr_start = $1;
                    }

                    if ( $ltr_end =~ /.*_end_(\d+)_.*/ ) {
                        $ltr_end = $1;
                    }

                    if ( $ltr_size =~ /.*size_(\d+)$/ ) {
                        $ltr_size = $1;
                    }

                    if ( $ltr_start <= 100 ) {
                        $has_five_prime_ltr = 1;
                    }

                }
            }

            if ( scalar @ltr_three_prime > 0 ) {
                foreach my $mapped_ltr ( @ltr_three_prime ) {
                    my $ltr_start = $mapped_ltr;
                    my $ltr_end   = $mapped_ltr;
                    my $ltr_size  = $mapped_ltr;

                    if ( $ltr_start =~ /.*start_(\d+)_.*/ ) {
                        $ltr_start = $1;
                    }

                    if ( $ltr_end =~ /.*_end_(\d+)_.*/ ) {
                        $ltr_end = $1;
                    }

                    if ( $ltr_size =~ /.*size_(\d+)$/ ) {
                        $ltr_size = $1;
                    }

                    if ( $ltr_end >= ( $self->genbank->length - 100 ) ) {
                        $has_three_prime_ltr = 1;
                    }

                }
            }

            if ( $has_three_prime_ltr == 0 ) {
                $classification = 'missing_ltr_or_psi';
            }
            elsif ( $has_five_prime_ltr == 0 and $json_report_current->{ 'has_psi' } eq "No" ) {
                $classification = 'missing_ltr_or_psi';
            }
            else {
                if ( $genes_count < 10 ) {
                    $classification = 'missing_internal_genes';
                }
                elsif ( $genes_count == 10 ) {

                    if ( scalar @productive_genes == 9 ) {
                        if ( $psi_status eq "intact" ) {
                            $classification = 'intact';
                        }
                        elsif ( $psi_status eq "mutation" ) {
                            $classification = 'msd_mutation';
                        }
                    }
                    else {
                        $classification = 'non_functional';
                    }
                }
                else {
                    die "MORE THAN 10 FEATURES " . $self->genbank->id;
                }
            }
        }

        # Add classification
        $json_report_current->{ 'output_gb_classification' } = $self->output_dir->stringify . '/'.$classification.'/' . $self->genbank->id . '.gb';
        my $classif_extraction = $json_report_current->{ 'output_gb_classification' };
        if ( $classif_extraction =~ /.*\/(.*)\/(.*\.gb)/ ) {
            $json_report_current->{ 'final_classification' } = $1;
        }
    }


    sub has_expected_genes_size {
        my ($self,$gene_length,  $gene) = @_;
        my $expected_minimum_length = $self->gene_range_allowed->{$gene}->[0];
        my $expected_maximum_length = $self->gene_range_allowed->{$gene}->[1];

        if ($expected_minimum_length <= $gene_length and $gene_length <= $expected_maximum_length) {
            return 1;
        } else {
            return 0;
        }

    }


    sub export_genbank {

        my ($self, $json_report_hash) = @_;

        foreach my $seq (keys %{$json_report_hash->{'passed'}}) {
            my $folder = dirname($json_report_hash->{'passed'}{$seq}{$self->file_id}{'STEP'.$self->current_step}{'output_gb_classification'});
            make_path($folder);


            my $classif_text;
            if ($folder =~ /.*\/(.*)$/) {
                $classif_text = $1;
            }

            my $classif = Bio::SeqFeature::Generic->new(-start => 1, -end => $self->genbank->length, -primary_tag => 'region');
            $classif->add_tag_value('standard_name', $classif_text);
            $self->genbank->add_SeqFeature($classif);
            my $out = Bio::SeqIO->new(-file => '>'. $json_report_hash->{'passed'}{$seq}{$self->file_id}{'STEP'.$self->current_step}{'output_gb_classification'}, -format => 'genbank');
            $out->write_seq($self->genbank);
        }


    }


    sub run { 
        my ($self) = @_;
        my %json_report;

        # Checking psi existance
        $self->check_psi(\%json_report);

        # Checking 5' and 3' LTR
        $self->check_ltr(\%json_report);

        # Checking inversions
        $self->check_inversions(\%json_report);

        # Check gene order
        $self->check_gene_order(\%json_report);

        #Check duplications
        $self->check_duplication(\%json_report);

        # Sequence classification
        $self->sequence_classification(\%json_report);

        # Exporting genbank
        $self->export_genbank(\%json_report);

        # Exporting JSON
        $self->export_json(\%json_report);

        #p %json_report;
    }

 __PACKAGE__->meta->make_immutable;


