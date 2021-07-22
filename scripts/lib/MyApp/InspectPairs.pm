#!/usr/bin/env perl

package MyApp::InspectPairs;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use namespace::autoclean;
use File::Basename;
use Data::Printer;
use File::Find::Rule;

command_short_description q[This command is awesome];
command_long_description q[This command is so awesome, yadda yadda yadda];

    has_directory 'data_dir' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[Very important option!],
    );


    has_file 'output' => (
        traits        => [ 'AppOption' ],
        cmd_type      => 'option',
        required      => 1,
        documentation => q[Very important option!],
    );


    sub check_pairs {

        my ( $self, $array_ref ) = @_;

        my %hash_files;
        my @files = File::Find::Rule->file()->name( "*.fastq.gz" )->in( $self->data_dir );

        foreach my $file ( @files ) {

            my $filename = basename( $file );

            if ( $filename =~ m/(\S+)_R[12]_(00.*)\.fastq\.gz/ ) {

                push @{ $hash_files{ $1.'_'.$2 } }, $file;
            }
        }


        foreach my $pair ( sort { $a cmp $b } keys %hash_files ) {

            my $sample_count = scalar @{ $hash_files{ $pair } };

            if ( $sample_count != 2 ) {

                push @{ $array_ref }, $pair;

            }


        }
    }


    sub run {
        my ( $self ) = @_;
        my @unpaired_samples;
        my $scalar_unpaired_samples;

        $self->check_pairs( \@unpaired_samples );

        $scalar_unpaired_samples = scalar @unpaired_samples;

        if ( $scalar_unpaired_samples == 0 ) {

            my $cmd = "touch " . $self->output;
            system( $cmd);

        }
        else {

            say "SAMPLES WITH PROBLEMS";
            say join( "\n", @unpaired_samples );

        }


    }

__PACKAGE__->meta->make_immutable;
