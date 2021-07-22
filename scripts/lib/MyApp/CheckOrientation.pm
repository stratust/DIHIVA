#!/usr/bin/env perl

package MyApp::CheckOrientation;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use namespace::autoclean;
use Bio::Tools::Run::StandAloneBlastPlus;
use Data::Printer;

    command_short_description q[ Get HIV sequence strand ];
    command_long_description q[ Get HIV sequence strand ];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        documentation => q[Contig file!]
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(o)],
        required      => 1,
        documentation => q[Output filename!]
    );

    has_file 'reference_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[HXB2 reference]
    );

    has_file 'report_file' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        cmd_aliases   => [ qw(r) ],
        required      => 0,
        documentation => q[ Output report filename ],
    );

    option 'min_hsp_length' => (
        is            => 'rw',
        isa           => 'Int',
        required      => '1',
        default       => '250',
        documentation => q[Minimum HSP length ],
    );

    option 'min_hsp_identity' => (
        is            => 'rw',
        isa           => 'Int',
        required      => '1',
        default       => '50',
        documentation => q[Minimum HSP percent_identity ],
    );

    has contigs => (
        is         => 'rw',
        isa        => 'HashRef',
        lazy_build => 1,
        documentation => q[HashRef of Bio::Seq contigs],
    );


    sub _build_contigs {
        my ($self) = @_;
        my %hash;
        my $in = Bio::SeqIO->new(
            -file   => $self->input_file->stringify,
            -format => 'fasta'
        );
        while ( my $seq  = $in->next_seq ) {
            $hash{ $seq->id } = $seq;
        }
        return \%hash;
    }


    sub run_blast {
        my ( $self ) = @_;

        my $fac;

        # $fac = Bio::Tools::Run::StandAloneBlastPlus->new( -DB_NAME => $self->blast_db_directory->stringify . '/full' );
        $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
            -db_data => $self->reference_file->stringify,
            -create  => 1
        );

        my $blast_out_file = $self->report_file->stringify;

        $self->log->warn( "Starting BLAST ..." );
        $fac->blastn(
            -query       => $self->input_file->stringify,
            -outfile     => $blast_out_file,
            -method_args => [
                '-task' => 'blastn',
                '-num_threads'   => 1,
                '-perc_identity' => 50,
            ]
        );

        $fac->rewind_results;
        $self->log->warn( "BLAST done!" );

        return $fac;
    }


    sub check_strand {

        my ( $self, $fac ) = @_;

        RESULT: while ( my $result = $fac->next_result ) {

            INTERATION: while ( my $iteration = $result->next_iteration ) {

                HIT: while ( my $hit = $iteration->next_hit ) {

                    HSP: while ( my $hsp = $hit->next_hsp ) {

                            if ( $hsp->length( 'total' ) > $self->min_hsp_length and $hsp->percent_identity >= $self->min_hsp_identity ) {

                                my $contig_ref = $self->contigs;
                                my ($seq_obj)  = sort { $b->length <=> $a->length } values %{ $contig_ref };

                                if ( $hsp->strand( 'query' ) ne $hsp->strand( 'hit' ) ) {
                                    $seq_obj = $seq_obj->revcom;
                                    $seq_obj->id( $seq_obj->id . '_rv' );
                                }

                                my $outseq = Bio::SeqIO->new(
                                    -file   => ">" . $self->output_file->stringify,
                                    -format => 'fasta',
                                );
                                $outseq->write_seq( $seq_obj );

                                last RESULT;

                            }
                        }
                    }
                }
            }
    }


    sub run {
        my ($self) = @_;

        my $fac = $self->run_blast();
        $fac->rewind_results;

        $self->check_strand($fac); 

        $fac->cleanup;
    }

    __PACKAGE__->meta->make_immutable;
