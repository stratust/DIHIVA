#!/usr/bin/env perl

package MyApp::CreateJSONReference;
use MooseX::App::Command;
extends 'MyApp';
with 'MyApp::Role::JSON';
use namespace::autoclean;
use Bio::SeqIO;

command_short_description
  q[Create a JSON file from all fasta files generated by Spades.];


sub get_reference_info {
    my ($self, $foldername) = @_;
    my $fasta_ref = "results/STEP3/annot/". $foldername . ".fasta";

    my $seqio_object = Bio::SeqIO->new(-file => $fasta_ref, -format => "fasta");
    my $seq_record = $seqio_object->next_seq;
    my $seq_ref_length = $seq_record->length();
    my $ref_id = $seq_record->id;

    return($ref_id, $seq_ref_length);

}

sub filter_samples_fasta_step_three {

    my ( $self, @fasta_files, ) = @_;
    my %json_hash;
    my $filename;
    my $foldername;
    my $verification;

    foreach my $file (@fasta_files) {

        my $sample_id;
        if ( $file =~ /^.*\/(\S+)\_S\d+\_.*/ ) {
            $sample_id = $1;
        }
        else {
            die "Cannot parse filename and retrive sample_id and pair information: $file";
        }

        my $reads = 0;
        my $cmd   = 'grep -P "^>\S+" ' . $file . '| wc -l';
        $reads = qx/$cmd/;
        chomp $reads;

        if ( $file =~ /^.*spades\/(.*)\/clean_scaffolds.*/ ) {
            $foldername = $1;
        }
        if ( $file =~ /^.*annot\/(.*)\.fasta$/ ) {
            $foldername = $1;
        }

        if ( $reads >= $self->min_reads ) {
            my @seq_info = $self->get_reference_info($foldername);
            $json_hash{'passed'}{$sample_id}{$foldername}{ 'STEP' . $self->current_step }{'input'}{ 'results/STEP2/spades/' . $foldername . '/scaffolds.fasta' } = $reads;
            $json_hash{'passed'}{$sample_id}{$foldername}{ 'STEP' . $self->current_step }{'output'}{$file} = $reads;
            $json_hash{'passed'}{$sample_id}{$foldername}{ 'STEP' . $self->current_step }{'ref_id'} = $seq_info[0];
            $json_hash{'passed'}{$sample_id}{$foldername}{ 'STEP' . $self->current_step }{'ref_length'} = $seq_info[1];
        } else {
            $json_hash{'failed'}{$sample_id}{$foldername}{ 'STEP' . $self->current_step }{'input'}{ 'results/STEP2/spades/' . $foldername . '/scaffolds.fasta' } = $reads;
            $json_hash{'failed'}{$sample_id}{$foldername}{ 'STEP' . $self->current_step }{'output'}{$file} = $reads;
        }

    }

    $self->export_json( \%json_hash );

}

sub run {
    my ($self) = @_;
    my @fasta_files;

    # Get subdirs
    my @dirs = glob( $self->input_directory . '/*' );
    foreach my $dir (@dirs) {
        if ( -d $dir ) {
            my @aux_fasta_files = glob( $dir . '/*.' . $self->file_extension );
            push @fasta_files, @aux_fasta_files if @aux_fasta_files;
        }
    }

    $self->filter_samples_fasta_step_three(@fasta_files);

}

__PACKAGE__->meta->make_immutable;
1;