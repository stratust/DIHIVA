#!/usr/bin/env perl

package MyApp::SearchReference;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use Bio::Tools::Run::StandAloneBlastPlus;
use List::MoreUtils qw(uniq);
use Data::Printer;
use namespace::autoclean;

    command_short_description q[Search the best HIV reference using blast.];
    command_long_description q[Search the best HIV reference using blast.
        The default behavior is to search for ENV references.
        For full genomes use: --full_genome option;
        For vpu-env-nef use: --vpu_to_nef option;];

    has_file 'input_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(i)],
        required      => 1,
        must_exist    => 1,
        default       => 'nt.fa',
        documentation => q[Contigs File!],
    );

    has_file 'output_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(o)],
        required      => 1,
        documentation => q[Output filename  ],
    );

    has_file 'report_file' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(r)],
        required      => 0,
        documentation => q[Output report filename  ],
    );

    has_directory 'blast_db_directory' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        cmd_aliases   => [qw(d)],
        required      => 1,
        must_exist    => 1,
        default       => '/work/blast/db2',
        documentation => q[Blast DB directory!],
    );

    has_directory 'hiv_genbank_directory' => (
        traits        => ['AppOption'],
        cmd_type      => 'option',
        required      => 1,
        must_exist    => 1,
        default       => '',
        documentation => q[Directory contained all genbank files used to create the blastdb!],
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
        default       => '75',
        documentation => q[Minimum HSP percent_identity ],
    );

    option 'max_contigs_to_analyze' => (
        is            => 'rw',
        isa           => 'Int',
        required      => '1',
        default       => '10',
        documentation => q[Maximum number of large contigs to analyze ],
    );

     option 'full_genome' => (
        is            => 'rw',
        isa           => 'Bool',
        required      => '0',
        documentation => q[Search only for full sequenced genomes],
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
        my ($self) = @_;

        # existing blastdb:
        #my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
        #    -DB_NAME => $self->blast_db_directory->stringify . '/nt' );

        my $fac;

        $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
            -DB_NAME => $self->blast_db_directory->stringify . '/full'
        );

        #my $fac = Bio::Tools::Run::StandAloneBlastPlus->new(
        #    -DB_NAME => 'nt', remote => 1 );

        my $blast_out_file = $self->report_file->stringify;
        $blast_out_file =~ s/\..*//g;
        $blast_out_file .= "_blast_report.txt";


        $self->log->warn("Starting BLAST ...");
        $fac->blastn(
            -query       => $self->input_file->stringify,
            -outfile     => $blast_out_file,
            -method_args => [
                '-task' => 'megablast',
                #'-num_alignments' => 8,
                #'-num_descriptions' => 8, #hits
                '-num_threads'    => 1,
                '-evalue'         => 1e-10,
                '-perc_identity'   => 80,
            ]
        );

        $fac->rewind_results;
        $self->log->warn("BLAST done!");
        return $fac;
    }


    sub get_cds_sequence_for_accession {
        my ( $self, $seq_obj ) = @_;

        my %parts;
        say "\tChecking: " .$seq_obj->accession;

        for my $feat_object ( $seq_obj->get_SeqFeatures ) {
            if ( $feat_object->primary_tag eq "CDS" ) {
                if ( $feat_object->has_tag('gene') ) {
                    for my $val ( $feat_object->get_tag_values('gene') ) {
                        if ( $val =~ /^env|gag|nef|vpu$/i ) {
                            say "\t\tFor $val:";
                            if ($val =~/gag/i ) {
                                # check accession reference contains psi packaging region
                                my $gag_start = $feat_object->location->start;
                                if ($gag_start > 130 ) {
                                    say "\t\t\tHas Psi";
                                }
                                else {
                                    say "\t\t\tNot Psi found";
                                    next;
                                }

                            }
                            # check if the feature is complete
                            if ( $feat_object->location->start_pos_type eq
                                'EXACT'
                                && $feat_object->location->end_pos_type eq
                                'EXACT' )
                            {

                                my $protein = $feat_object->seq->translate->seq;

                                next unless $protein =~ /\*$/;

                                say "\t\t\tHas stopcodon at the end!";
                                # remove last stop codon if exists;
                                $protein =~ s/^(\S+)\*$/$1/g;

                                unless  ( $protein =~ /\*/  ) {
                                    $parts{lc($val)} = $feat_object;
                                    say "\t\t\tNo internal stopcodons!";
                                }
                            }
                            else {
                                say "\t\t\tIncomplete $val!";
                            }
                        }
                    }
                }
            }
        }

        if ( %parts &&  $parts{gag} && $parts{env} && $parts{nef} ){
                return $seq_obj;
        }
        else {
            return 0;
        }
    }


    sub run {
        my ($self) = @_;
        my $cmd;
        $cmd = $1 if __PACKAGE__ =~ /\:\:(.*)$/;
        $self->log->warn("==> Starting $cmd <==");

        # Code Here
        my $fac = $self->run_blast();
        $fac->rewind_results;

        my @accessions;

        my $found_ref = 0;

        my $n_results = 0;
        RESULT: while ( my $result = $fac->next_result ) {

            ITERATION: while ( my $iteration = $result->next_iteration ) {

                HIT: while ( my $hit = $iteration->next_hit ) {

                    HSP: while ( my $hsp = $hit->next_hsp ) {

                        if ( $hsp->length('total') > $self->min_hsp_length ) {
                            if ( $hsp->percent_identity >=
                                $self->min_hsp_identity )
                            {
                                push @accessions, $hit->accession;
                                say $hit->accession;
                            }

                        }

                    }

                }

            }

            $n_results++;
            last RESULT if $n_results >= $self->max_contigs_to_analyze;
        }

        $fac->rewind_results;

        #Get genbank file
        #        my $db_obj  = Bio::DB::GenBank->new;
        @accessions = uniq @accessions;
        @accessions = @accessions[0..49] if scalar @accessions > 50;

        $self->log->warn("Retrieving accession info by Entrez batch ...");
        
        my %seq_objs;
        foreach my $acc (@accessions){
        
            # check if genbank file exist
            my $gb_file = $self->hiv_genbank_directory->stringify."/".$acc.".gb";
            if ( -e $gb_file ){
                my $in = Bio::SeqIO->new(-file => $gb_file, -format => 'Genbank');
                my $seq_obj = $in->next_seq;
                say $seq_obj->accession ," => " ,$seq_obj->length;
                $seq_objs{$seq_obj->accession} = $seq_obj;
            }
            else {
                die "Cannot find file $gb_file!!";
            }
        }
        #my $seq_obj_stream = $db_obj->get_Stream_by_id(\@accessions);

        $self->log->warn("All available info retrieved!");


        my %contamination;
        my %selected;

        my %contig_hits;
        my %contig_hits_info;

        $n_results = 0;

        RESULT: while ( my $result = $fac->next_result ) {

            ITERATION: while ( my $iteration = $result->next_iteration ) {

                # process all search results within the input stream
                HIT: while ( my $hit = $iteration->next_hit ) {
                    # insert code here for hit processing
                    ## $hit is a Bio::Search::Hit::HitI compliant object
                    HSP: while ( my $hsp = $hit->next_hsp ) {
                        ## $hsp is a Bio::Search::HSP::HSPI compliant object
                        if ( $hsp->length('total') > $self->min_hsp_length ) {
                            if ( $hsp->percent_identity >= $self->min_hsp_identity ) {
                                                  # Check contamination
                                my $seq_obj = $seq_objs{$hit->accession};
                                next HSP unless $seq_obj;

                                if ( $seq_obj->species->scientific_name =~ /Human immunodeficiency virus/i ){
                                    push @{ $contig_hits{$result->query_name} }, $seq_obj;
                                    push @{ $contig_hits_info{$result->query_name}{$seq_obj->accession}}, $hsp;
 
                                    push @{ $selected{$result->query_name} }, $hit->accession . ' (' . $seq_obj->species->scientific_name  .')'. ' - '. $hit->description;
                                    last HSP;
                                }
                                else{
                                    push @{ $contamination{$result->query_name} }, $hit->accession . ' (' . $seq_obj->species->scientific_name  .')'. ' - '. $hit->description;
                                    last HSP;
                                }
                           }

                        }

                    }

                }

            }

            $n_results++;
            last RESULT if $n_results >= $self->max_contigs_to_analyze;
        }
        $fac->cleanup;


        my $selected_acc = "None";
        CONTIG: foreach my $contig_name (
            sort { ( $a =~ /NODE_(\d+)_length/i )[0] <=> ( $b =~ /NODE_(\d+)_length/i )[0] }
            keys %contig_hits
          )
        {
            ACC: foreach my $seq_obj ( @{ $contig_hits{$contig_name} } ) {
                my $sequence;
            
                $sequence = $self->get_cds_sequence_for_accession($seq_obj);
                if ($sequence) {
                    $selected_acc = $seq_obj->accession;
                    $sequence->id($selected_acc);
                    # Get the longest contig from spades
                    my $contig_seqobj = $self->contigs->{$contig_name};
                    my $hsp = $contig_hits_info{$contig_name}{$selected_acc}[0];

                    if ($hsp->strand('query') ne $hsp->strand('hit')){
                        $contig_seqobj = $contig_seqobj->revcom;
                        $contig_seqobj->id($contig_seqobj->id.'_rv')
                    }

                    $contig_seqobj->id(
                        $contig_seqobj->id.
                        '_'.
                        $selected_acc
                    );

                    $contig_seqobj->desc(
                        'HSP:'.$hsp->length('hit')."/".$hsp->query->seqlength);


                    my $outseq = Bio::SeqIO->new(
                        -file   => ">" . $self->output_file->stringify,
                        -format => 'fasta',
                    );
                    #$outseq->write_seq($sequence);
                    $outseq->write_seq($contig_seqobj);
                    last CONTIG;
                }

            }

        }

        my $report_file = $self->report_file;
        $report_file = $self->input_file->stringify.'.report.txt' unless $report_file;


        open( my $out, '>', $report_file )  || die "Cannot open/write file " . $report_file . "!";

        say $out "--------------------------------------------";
        say $out $self->input_file->stringify;
        say $out "--------------------------------------------";
        say $out "SELECTED: ".$selected_acc;


        say $out "List of selected:";
        foreach my $k (sort { ( $a =~ /NODE_(\d+)_length/i )[0] <=> ( $b =~ /NODE_(\d+)_length/i )[0]} keys %selected) {
            say $out "$k:";
            foreach my $row (@{$selected{$k}}) {
                say $out "\t$row";
            }
        }

        say $out "--------------------------------------------";
        say $out "Contamination: ";
        say $out "--------------------------------------------";

        foreach my $k (sort { ( $a =~ /NODE_(\d+)_length/ )[0] <=> ( $b =~ /NODE_(\d+)_length/ )[0]} keys %contamination) {
            say $out "$k:";
            foreach my $row (@{$contamination{$k}}) {
                say $out "\t$row";
            }
        }

        close( $out );

        $self->log->warn("==> END $cmd <==");
    }

__PACKAGE__->meta->make_immutable;

