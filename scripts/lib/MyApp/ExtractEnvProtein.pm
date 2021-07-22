#!/usr/bin/env perl

package MyApp::ExtractEnvProtein;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use Bio::SeqIO;
use Data::Printer;
use namespace::autoclean;

command_short_description q[This command is awesome];
command_long_description q[This command is so awesome, yadda yadda yadda];

has_directory 'genbank_directory' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    cmd_aliases   => [qw(d)],
    required      => 1,
    must_exist    => 1,
    default       => "results/STEP9/defective",
    documentation => q[Very important option!],
);

has_directory 'fasta_directory' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    cmd_aliases   => [qw(o)],
    required      => 1,
    must_exist    => 0,
    default       => "results/STEP9/env_fasta",
    documentation => q[Very important option!],
);

sub run {
    my ($self) = @_;
    my @files = glob( $self->genbank_directory . '/*.gb' );
    foreach my $file (@files) {
        my $in  = Bio::SeqIO->new( -file => $file, -format => 'genbank' );
        my $seq = $in->next_seq;
        my @cds_features =
          grep { $_->primary_tag eq 'gene' } $seq->get_SeqFeatures;
        foreach my $feat (@cds_features) {
            my @values = $feat->get_tag_values('gene');
            if ( $values[0] eq 'env' ) {
                my $prot     = $feat->seq->translate;
                my $prot_str = $prot->seq;
                $prot_str =~ s/\*$//g;
                if ( $prot_str !~ /\*/ ) {
                    mkdir $self->fasta_directory
                      unless ( -e $self->fasta_directory );
                    my $fname;
                    $fname = $1 . ".fasta" if $file =~ /.*\/(.*)\.gb/;
                    my $out = Bio::SeqIO->new(
                        -file   => ">" . $self->fasta_directory . '/' . $fname,
                        -format => 'fasta'
                    );
                    $prot->desc('');
                    $out->write_seq($prot);
                }
            }
        }
    }
}

__PACKAGE__->meta->make_immutable;
1;
