#!/usr/bin/env perl

package MyApp::CreateJSONContigRef;
use MooseX::App::Command;
extends 'MyApp';
with 'MyApp::Role::JSON';
use namespace::autoclean;

command_short_description
  q[Create a JSON file from all fasta files generated by Spades.];

sub run {
    my ($self) = @_;
    my @fasta_files;

    @fasta_files =
      glob( $self->input_directory . '/*.' . $self->file_extension );

    $self->filter_samples_fasta_modified(@fasta_files);
}

__PACKAGE__->meta->make_immutable;
1;