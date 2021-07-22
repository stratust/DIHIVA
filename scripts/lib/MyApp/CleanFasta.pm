#!/usr/bin/env perl

package MyApp::CleanFasta;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use Bio::SeqIO;
use namespace::autoclean;


    command_short_description q[This command is awesome];
    command_long_description q[This command is so awesome, yadda yadda yadda];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        documentation => q[Fasta Sequence!],
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(o)],
        documentation => q[Output Fasta Sequence!],
    );

    sub run {
        my ($self) = @_;
        my $output;
        if (! $self->output_file){
            if ( $self->input_file->stringify =~ /(.*)\/(.*)/){
                $output = "$1/clean_$2";
            }
        }
        else {
            $output = $self->output_file->stringify;
        }

        my $in = Bio::SeqIO->new(-file => $self->input_file->stringify, -format => 'fasta');

        my $out = Bio::SeqIO->new(-file => ">$output", -format => 'fasta');
        while (my  $seq = $in->next_seq){
            my $str_seq = $seq->seq;
            $str_seq =~ s/n//ig;
            $seq->seq($str_seq);
            $out->write_seq($seq);
        }
    }

__PACKAGE__->meta->make_immutable;
