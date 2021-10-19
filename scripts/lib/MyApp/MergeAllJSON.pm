#!/usr/bin/env perl

package MyApp::MergeAllJSON;
use MooseX::App::Command;
use MooseX::FileAttribute;
use feature "say";
extends 'MyApp';
use namespace::autoclean;
with 'MyApp::Role::JSON';
use File::JSON::Slurper qw/ read_json write_json  /;
use Hash::Merge qw / merge /;
use Data::Printer;
use List::MoreUtils qw(uniq);
use File::Copy;
use File::Basename;

has_directory 'results_dir' => (
    traits   => ['AppOption'],
    cmd_type => 'option',
    required => 1
);


sub merge_json {
    my ($self, $json_file) = @_;
    my $merged_json = read_json($self->results_dir->stringify . '/JSON_STEP0.json' );
    #my $merged_json = read_json(shift @{$json_file});

    foreach my $json ( @{$json_file} ) {
        $merged_json = merge($merged_json, read_json($json) );
    }

    return($merged_json);
}


sub format_json {
    my ( $self, $json_file ) = @_;
    my @passed_samples = keys( %{ $json_file->{ 'passed' } } );
    my @failed_samples = keys( %{ $json_file->{ 'failed' } } );
    my @all_samples    = uniq( ( @passed_samples, @failed_samples ) );
    my %formatted_json;

    foreach my $sample ( @all_samples ) {

        my ( $formatted_merged_json );

        if ( $json_file->{ 'passed' }{ $sample } and $json_file->{ 'failed' }{ $sample } ) {
            $formatted_merged_json = merge( $json_file->{ 'passed' }{ $sample }, $json_file->{ 'failed' }{ $sample } );
        }
        elsif ( $json_file->{ 'passed' }{ $sample } ) {
            $formatted_merged_json = $json_file->{ 'passed' }{ $sample };
        }
        elsif ( $json_file->{ 'failed' }{ $sample } ) {
            $formatted_merged_json = $json_file->{ 'failed' }{ $sample };
        }
        $formatted_json{ $sample } = $formatted_merged_json;
    }

    return \%formatted_json;
}


sub choose_representative_file {
    my ( $self, $json_by_sample ) = @_;

    # Go over json and index/choose file to represent sample
    my %aux_hash;
    foreach my $sample ( sort { $a cmp $b } keys %{ $json_by_sample } ) {
        my @files = sort { $a cmp $b } keys %{ $json_by_sample->{ $sample } };
        # if there is one file perl sample
        if ( scalar @files == 1 ) {
            $json_by_sample->{ $sample }{ 'chosen' } = $files[ 0 ];
        }
        else {

            #create a hash indexing sample->classification->file
            foreach my $file ( @files ) {
                my $cur_file = $json_by_sample->{ $sample }{ $file };

                my $final_classification;
                if ($cur_file->{ 'STEP6' }){
                    $final_classification = $cur_file->{ 'STEP6' }{ 'final_classification' };
                }
                if ( $final_classification ) {
                    push @{ $aux_hash{ $sample }{ $final_classification } }, $file;
                }
                else {
                    my $classification;
                    if ($cur_file->{ 'STEP5' }){
                        $classification = $cur_file->{ 'STEP5' }{ 'classification' };
                    }
                    if ( $classification and $classification =~ /bad/i ) {
                        push @{ $aux_hash{ $sample }{ $classification } }, $file;
                    }
                    else {
                        my $classification;
                        if ($cur_file->{ 'STEP2' }){
                            $classification = $cur_file->{ 'STEP2' }{ 'classification' };
                        }
                        if ( $classification ) {
                            push @{ $aux_hash{ $sample }{ $classification } }, $file;
                        }
                        else {
                            my $classification;
                            if ($cur_file->{ 'STEP1' }){
                                $classification = $cur_file->{ 'STEP1' }{ 'classification' };
                            }
                            if ( $classification ) {
                                push @{ $aux_hash{ $sample }{ $classification . ' STEP1' } }, $file;
                            }
                            else {
                                if ($cur_file->{ 'STEP0' }){
                                    $classification = $cur_file->{ 'STEP0' }{ 'classification' };
                                }
                                if ( $classification ) {
                                    push @{ $aux_hash{ $sample }{ $classification . ' STEP0' } }, $file;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #chose if there is more than one file per sample
    my @classif_priorities = ( 'intact', 'msd_mutation', 'non_functional', 'structural_variations', 'missing_internal_genes', 'missing_ltr_or_psi', 'BAD', 'problems_in_assembly' , 'EMPTY WELL STEP1', 'EMPTY WELL STEP0' );
    foreach my $sample_hash ( sort { $a cmp $b } keys %aux_hash ) {

        foreach my $classif ( @classif_priorities ) {

            if ( $aux_hash{ $sample_hash }{ $classif } ) {

                my @aux_hash_samples = @{$aux_hash{ $sample_hash }{ $classif }};
                my $chosen_sample = shift @aux_hash_samples;

                $json_by_sample->{ $sample_hash }{ 'chosen' } = $chosen_sample;
                $json_by_sample->{ $sample_hash }{ 'not_chosen' } = join(",", @aux_hash_samples);
                last;

            }
        }

    }
}


sub create_table_from_json {
    my ($self, $json_by_sample) = @_;
    my @rows;
    my @header = qw/sample filename well final_classification classification psi_status reason steps total_reads filtered_reads productive_genes unproductive_genes env_double_peaks n_double_peaks_sites n_low_coverage_sites percent_sequence_coverage percent_sequence_with_high_coverage ref_id ref_length has_duplication has_invertions has_ltr has_order has_psi output_step7 chosen/;
    foreach my $sample (sort {$a cmp $b} keys %{$json_by_sample}){
        my $chosen = $json_by_sample->{$sample}->{chosen};
        foreach my $filename (sort {$a cmp $b} keys %{$json_by_sample->{$sample}}){
            next if $filename =~ /^chosen|not_chosen$/;

            my $current_file = $json_by_sample->{$sample}->{$filename};

            my @passed_steps;
            # create hash with all header keys and empty values
            my %columns = map { $_ => '' } @header;

            $columns{sample} = $sample;
            $columns{filename} = $filename;
            $columns{chosen} = 1 if $filename eq $chosen;

            # STEP6
            if ($current_file->{STEP6}){
                if ( $current_file->{ STEP6 }->{ productive_genes } or $current_file->{ STEP6 }->{ unproductive_genes } ) {
                    my $productive_genes   = $current_file->{ STEP6 }->{ productive_genes };
                    my $unproductive_genes = $current_file->{ STEP6 }->{ unproductive_genes };
                    $productive_genes   =~ s/, /|/g if $current_file->{ STEP6 }->{ productive_genes };
                    $unproductive_genes =~ s/, /|/g if $current_file->{ STEP6 }->{ unproductive_genes };
                    $columns{ 'productive_genes' }   = $productive_genes;
                    $columns{ 'unproductive_genes' } = $unproductive_genes;

                }
                $columns{'has_duplication'} = $current_file->{STEP6}->{has_duplication};
                $columns{'has_invertions'} = $current_file->{STEP6}->{has_invertions};
                $columns{'has_ltr'} = $current_file->{STEP6}->{has_ltr};
                $columns{'has_order'} = $current_file->{STEP6}->{has_order};
                $columns{'has_psi'} = $current_file->{STEP6}->{has_psi};
                $columns{'final_classification'} = $current_file->{STEP6}->{final_classification};
                $columns{'psi_status'} = $current_file->{STEP6}->{psi_status} if $current_file->{STEP6}->{psi_status};


                my $out_gb_classification =  $current_file->{STEP6}->{output_gb_classification};
                if ($out_gb_classification =~ m/^.*\/(results\/.*)/) {
                    $out_gb_classification = $1;
                    $out_gb_classification =~ s/STEP6/STEP7/g;

                }
                $columns{'output_step7'} = $out_gb_classification;

                push @passed_steps, 6;
            }

            # STEP5
            if ($current_file->{STEP5}){
                $columns{'classification'} = $current_file->{STEP5}->{classification};
                $columns{'ref_id'} = $current_file->{STEP5}->{ref_id};
                $columns{'ref_length'} = $current_file->{STEP5}->{ref_length};
                $columns{'n_double_peaks_sites'} = $current_file->{STEP5}->{n_double_peaks_sites};

                my $env_double_peaks = scalar keys %{$current_file->{STEP5}->{double_peaks_sites_by_feature}{'env gene'}} if $current_file->{STEP5}->{double_peaks_sites_by_feature}{'env gene'};
                if ($env_double_peaks) {
                    $columns{'env_double_peaks'} = $env_double_peaks;
                } else {
                    $columns{'env_double_peaks'} = 0;
                }
                $columns{'n_low_coverage_sites'} = $current_file->{STEP5}->{n_low_coverage_sites};
                $columns{'percent_sequence_coverage'} = $current_file->{STEP5}->{percent_sequence_coverage};
                $columns{'percent_sequence_with_high_coverage'} = $current_file->{STEP5}->{percent_sequence_with_high_coverage};
                $columns{'reason'} = $current_file->{STEP5}->{reason} if $current_file->{STEP5}->{reason};
                push @passed_steps, 5;
            }

            # STEP4
            if ($current_file->{STEP4}){
                push @passed_steps, 4;
            }

            # STEP3
            if ($current_file->{STEP3}){
                push @passed_steps, 3;
            }

            # STEP2
            if ($current_file->{STEP2}){
                $columns{'reason'} = $current_file->{STEP2}->{reason} if $current_file->{STEP2}->{reason};
                $columns{'classification'} = $current_file->{STEP2}->{classification} unless $columns{'classification'};
                push @passed_steps, 2;
            }

            # STEP1
            if ($current_file->{STEP1}){
                $columns{'filtered_reads'} = $current_file->{STEP1}->{output}->{R1}->{reads};
                $columns{'reason'} = $current_file->{STEP1}->{reason} if ($current_file->{STEP1}->{reason});
                $columns{'classification'} = $current_file->{STEP1}->{classification} unless $columns{'classification'};
                push @passed_steps, 1;
            }

            # STEP0
            if ($current_file->{STEP0}){
                my $well = $current_file->{STEP0}->{well};
                $well =~ s/\D+//g;
                $columns{'well'} = $well;

                $columns{'total_reads'} = $current_file->{STEP0}->{input}->{R1}->{reads};
                $columns{'reason'} = $current_file->{STEP0}->{reason} if ($current_file->{STEP0}->{reason});
                $columns{'classification'} = $current_file->{STEP0}->{classification} unless $columns{'classification'};
                push @passed_steps, 0;
            }

            $columns{steps} = shift @passed_steps;
            push @rows, \%columns;
        }
    }
    #print table
    say join ",", @header;
    foreach my $row (@rows){
        #foreach my $col (@header){
            say  join ",", @{$row}{@header};
#             p $aux;
            # p @header;
            #p $row;
        #}
    }
}


sub copy_chosen_gb_to_folder {

    my ( $self, $json_by_sample ) = @_;

    foreach my $sample ( sort { $a cmp $b } keys %{ $json_by_sample } ) {

        my $chosen_file = $json_by_sample->{ $sample }->{ chosen };

        if ( $json_by_sample->{ $sample }->{ $chosen_file }->{ STEP6 } ) {
            my $out_gb_classification = $json_by_sample->{ $sample }->{ $chosen_file }->{ STEP6 }->{ output_gb_classification };
            my $output_step7          = $out_gb_classification;

            $output_step7 =~ s/STEP6/STEP7/g;

            my $dirname = dirname( $output_step7 );

            if ( !-e $dirname ) {
                mkdir $dirname;
            }

            copy( $out_gb_classification, $output_step7 );

        }

    }
}


sub run {

    my ($self) = @_;
    my @json_file = glob($self->results_dir->stringify . '/STEP[123456]/*.json');
    my $ref = $self->merge_json(\@json_file);

    my $json_by_sample = $self->format_json($ref);

    $self->choose_representative_file($json_by_sample);

    $self->create_table_from_json($json_by_sample);

    $self->copy_chosen_gb_to_folder($json_by_sample);

    $self->export_json($json_by_sample);


}
