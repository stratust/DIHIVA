package HIV::Assembly::Bam {
    use Moose;
    use Bio::DB::Sam;
    use Bio::SeqIO;
    use Data::Printer;
    use File::Basename;
    use namespace::autoclean;
    use feature 'say';

    # Attributes
    # --------------------------------------------------------------------------------------------------
    has 'bam' => (
        is            => 'ro',
        required      => 1,
        documentation => 'BAM file path',
        trigger => \&_build_coverage_and_pile,
    );

    has 'ref_genbank_file' => (
        is            => 'ro',
        required      => 0,
        documentation => 'Reference file path in genbank format',
    );

    has '_bam_obj' => (
        is            => 'ro',
        isa           => 'Bio::DB::Sam',
        documentation => 'Read BAM file and returns a Bio::DB::Sam',
        lazy_build    => 1
    );

    has 'ref_length' => (is => 'ro', isa => 'Int', documentation => 'BAM reference length', lazy_build => 1);

    has 'ref_id' => (is => 'ro', isa => 'Str', documentation => 'BAM reference name', lazy_build => 1);

    has 'n_reads' => (is => 'ro', isa => 'Int', documentation => 'Number of read in the BAM file', lazy_build => 1);

    has 'coverage' => (
        is            => 'ro',
        isa           => 'HashRef',
        default       => sub { return {} },
        documentation => 'Coverage HashRef'
    );

    has 'coverage_by_feature' => (
        is            => 'ro',
        isa           => 'HashRef',
        default       => sub { return {} },
        documentation => 'CDS coverage HashRef'
    );

    has 'pile' => (
        is            => 'ro',
        isa           => 'HashRef',
        default       => sub { return {} },
        documentation => 'Pile HashRef'
    );

    has 'pile_by_feature' => (
        is            => 'ro',
        isa           => 'HashRef',
        default       => sub { return {} },
        documentation => 'CDS pile HashRef'
    );

    has 'dp_threshold' => (
        is       => 'rw',
        isa      => 'Num',
        default  => 0.7,
        required => 1,
        documentation =>
            'double peak threshold.(e.g. 0.7 means at least 70% reads agree)'
    );

    has 'mininum_coverage' => (
        is            => 'rw',
        isa           => 'Int',
        default       => 10,
        required      => 1,
        documentation => 'Mininimum coverage'
    );

    has 'min_reads' => (
        is => 'rw', 
        isa => 'Int', 
        default => 500,
        required => 1,
        documentation => 'Mininimum number of reads'
    );

    has 'primary_tag' => (
        is => 'rw', 
        isa => 'Str', 
        default => 'CDS',
        required => 1,
        documentation => 'Primary tag to search for features'
    );

    has 'tag' => (
        is => 'rw', 
        isa => 'Str', 
        default => 'standard_name',
        required => 1,
        documentation => 'Tag name to search for features'
    );

    has 'no_ambiguities' => (
        is => 'rw', 
        isa => 'Bool', 
        default => 0,
        required => 1,
        documentation => 'Avoid ambiguities'
    );



    # Builders and auxiliar subrotines
    # --------------------------------------------------------------------------------------------------
 
    sub callback {
        my ($self, $ref_id, $pos, $pileups) = @_;
        $self->coverage->{$pos} = 0 unless $self->coverage->{$pos};
        my $refbase = $self->_bam_obj->segment($ref_id,$pos,$pos)->dna;
        foreach my $pileup (@{$pileups}) {
            my $aln   = $pileup->alignment;
            my $qBase;
            if ($pileup->indel > 0 ){
                $qBase = substr($aln->qseq, $pileup->qpos, 1 + $pileup->indel);
            }
            elsif ($pileup->indel < 0){
                $qBase = '-' x abs($pileup->indel);
            }
            else {
                $qBase = substr($aln->qseq, $pileup->qpos, 1);
            }
            die "Cannot find base for $pos!" unless $qBase; 
            #say $aln->query->name. ': '. $pileup->indel.' - ' .$qBase;
            $self->pile->{$pos}->{$qBase}++;
            $self->coverage->{$pos}++;
        }
    }

    sub _build_cds_coverage_and_pile {
        my ($self ) = @_;
        my $in = Bio::SeqIO->new(-file => $self->ref_genbank_file, -format => 'genbank');
        # check if genbank has the same id as the bam file reference
        my $seq = $in->next_seq;
        if ($seq->id ne  $self->ref_id){
            die "Genbank sequence ID (".$seq->id.") doesn't match BAM reference name (".$self->ref_id.")!";
        }
        foreach my $feat ($seq->get_SeqFeatures) {
            if ($feat->primary_tag eq $self->primary_tag){
                my %pile_aux;
                my %cov_aux;
                my @values = $feat->get_tag_values($self->tag);
                my $name = shift @values;
                my @positions = $feat->start .. $feat->end; 
                @pile_aux{@positions} = @{$self->pile}{@positions};
                $self->pile_by_feature->{$name} = \%pile_aux;
                @cov_aux{@positions} = @{$self->coverage}{@positions};
                $self->coverage_by_feature->{$name} = \%cov_aux;
            }
        }
    }

    sub _build__bam_obj {
        my ($self)     = @_;
        my $bam        = Bio::DB::Sam->new(-bam => $self->bam, -autoindex => 1);
        return $bam;
    }

    sub _build_ref_length {
        my ($self)     = @_;
        my $header     = $self->_bam_obj->header;
        return $header->target_len->[0];
    }

    sub _build_ref_id {
        my ($self)     = @_;
        my $header     = $self->_bam_obj->header;
        my @targets = @{$header->target_name};
        if (scalar @targets > 1 ){
            die "There are more than one reference in the BAM file header. Please verify your BAM file."
        } else {
            return $targets[0]; 
        }
    }

    sub _build_n_reads {
        my ($self)     = @_;
        my $n_reads =0;
        my $low_level_bam = Bio::DB::Bam->open($self->bam);
        my $header = $low_level_bam->header; # necessary to use read1 function;
        while (my $aln = $low_level_bam->read1){
            $n_reads++;
        }
        return $n_reads;
    }


    sub _build_coverage_and_pile {
        my ($self) = @_;
        my $region = $self->ref_id . ":1-" . $self->ref_length;
        $self->_bam_obj->pileup($region, sub { $self->callback(@_) });
        
        # populate cds pile and coverage
        if ($self->ref_genbank_file){
            $self->_build_cds_coverage_and_pile;
        }
    }

    # Methods
    # --------------------------------------------------------------------------------------------------

    sub find_low_coverage_sites {
        my ($self) = @_;
        my %low_cov;
        foreach my $pos (keys %{$self->coverage}) {
            if ($self->coverage->{$pos} < $self->mininum_coverage){
                $low_cov{$pos} = $self->coverage->{$pos};
            }
        }
        return \%low_cov;
    }


    sub find_double_peaks_sites {
        my ($self) = @_;
        my %double_peaks;
        foreach my $pos (keys %{$self->pile}) {
            my @bases                 = keys %{$self->pile->{$pos}};
            my $bases_above_threshold = 0;
            foreach my $base (@bases) {
                my $fraction = $self->pile->{$pos}->{$base} / $self->coverage->{$pos};
                if ($fraction >= $self->dp_threshold) {
                    $bases_above_threshold++;
                }
            }
            if ($bases_above_threshold == 0) {
                $double_peaks{$pos} = $self->pile->{$pos};
            }
        }
        return \%double_peaks;
    }


    sub find_double_peaks_sites_by_feature {
        my ($self, $simplify) = @_;
        my %double_peaks_by_feature;
        my %double_peaks_by_feature_final;
        my $double_peaks = $self->find_double_peaks_sites;
        foreach my $region (keys %{$self->pile_by_feature}) {
            #$double_peaks_by_feature{$region} = {};
            foreach my $pos (keys %{$self->pile_by_feature->{$region}}) {
                if ($double_peaks->{$pos}){
                    $double_peaks_by_feature{$region}{$pos} = $double_peaks->{$pos};
                }
            }
        }

        if ($simplify){
            foreach my $region (keys %double_peaks_by_feature) {
                $double_peaks_by_feature_final{$region} = scalar keys %{$double_peaks_by_feature{$region}};
            }
        }
        else {
            %double_peaks_by_feature_final = %double_peaks_by_feature;
        }

        return \%double_peaks_by_feature_final;
    }

    sub find_low_coverage_sites_by_feature {
        my ($self, $simplify) = @_;
        my %low_coverage_by_feature;
        my %low_coverage_by_feature_final;
        my $low_coverage = $self->find_low_coverage_sites;
        foreach my $region (keys %{$self->coverage_by_feature}) {
            #$double_peaks_by_feature{$region} = {};
            foreach my $pos (keys %{$self->coverage_by_feature->{$region}}) {
                if ($low_coverage->{$pos}){
                    $low_coverage_by_feature{$region}{$pos} = $low_coverage->{$pos};
                }
            }
        }

        if ($simplify){
            foreach my $region (keys %low_coverage_by_feature) {
                $low_coverage_by_feature_final{$region} = scalar keys %{$low_coverage_by_feature{$region}};
            }
        }
        else {
            %low_coverage_by_feature_final = %low_coverage_by_feature;
        }

        return \%low_coverage_by_feature_final;
    }

    sub generate_consensus {
        my ($self ) = @_;
        my $consensus;
        my %IUPAC = (
            "A" => "A",
            "T" => "T",
            "C" => "C",
            "G" => "G",
            "AC" => "M",
            "AG" => "R",
            "AT" => "W",
            "CG" => "S",
            "CT" => "Y",
            "GT" => "K",
            "ACG" => "V",
            "ACT" => "H",
            "AGT" => "D",
            "CGT" => "B",
            "ACGT" => "N",
        );
        foreach my $pos (sort {$a <=> $b} keys %{$self->pile}) {
            # sort bases by frequency (desc)
            my @bases = sort { $self->pile->{$pos}->{$b} <=> $self->pile->{$pos}->{$a} or $a cmp $b } keys %{$self->pile->{$pos}};
            my @bases_above_threshold;
            my @bases_below_threshold;
            foreach my $base (@bases) {
                my $fraction = $self->pile->{$pos}->{$base} / $self->coverage->{$pos};
                if ($fraction >= $self->dp_threshold) {
                    push @bases_above_threshold, $base;
                }
                else {
                    push @bases_below_threshold, $base;
                }
            }
            if (scalar @bases_above_threshold == 1) {
                $consensus .= $bases_above_threshold[0];
            }
            elsif (scalar @bases_above_threshold > 1){
                if ($self->no_ambiguities){
                    $consensus .= $bases_above_threshold[0];
                }
                else {
                    my $sorted_bases = join "", sort {$a cmp $b} @bases_above_threshold;
                    # remove non ACGT from key
                    $sorted_bases =~ s/[^acgt]//gi;
                    if ($IUPAC{$sorted_bases}){
                        $consensus .= $IUPAC{$sorted_bases};
                    }
                    else {
                        die "Cannot find a IUPAC ambiguity code for: $sorted_bases !";
                    }
                }
            } # This should never be executed
            elsif (scalar @bases_below_threshold == 1){
                $consensus .= $bases_below_threshold[0];
            }
            elsif (scalar @bases_below_threshold > 1){
                if ($self->no_ambiguities){
                    $consensus .= $bases_below_threshold[0];
                }
                else {
                    my $sorted_bases = join "", sort {$a cmp $b} @bases_below_threshold;
                    $sorted_bases =~ s/[^acgt]//gi;
                    if ($IUPAC{$sorted_bases}){
                        $consensus .= $IUPAC{$sorted_bases};
                    }
                    else {
                        die "Cannot find a IUPAC ambiguity code for: $sorted_bases !";
                    }
                }
            }

            else{
                p $self->pile->{$pos};
                die "Cannot find base for consensus at position: ". $pos;
            }
              
        }

        # Remove gaps from consensus
        $consensus =~ s/\-//g;
        return $consensus;
    }

    sub consensus_as_fasta_file {
        my ($self, $fasta_file ) = @_;
        # seq id  will be the as bam file without extension
        my $filename = fileparse($self->bam,(".bam"));
        my $consensus = Bio::PrimarySeq->new(-id => $filename, -seq => $self->generate_consensus);
        my $out_consensus = Bio::SeqIO->new(-file => '>'.$fasta_file, -format => 'fasta');
        $out_consensus->write_seq($consensus);
    }

    sub get_bam_info {
        my ($self) =  @_;
        my %bam_info;

        my $n_sites_covered = scalar keys %{$self->coverage};
        my $n_low_coverage_sites = scalar keys %{$self->find_low_coverage_sites};
        my $percent_sequence_coverage = $n_sites_covered/$self->ref_length * 100;
        my $percent_sequence_high_coverage = ($n_sites_covered - $n_low_coverage_sites)/$self->ref_length * 100;
        
        $bam_info{'ref_id'}             = $self->ref_id;
        $bam_info{'ref_length'}         = $self->ref_length;
        $bam_info{'n_reads'}            = $self->n_reads;
        $bam_info{'n_low_coverage_sites'} = $n_low_coverage_sites;
        $bam_info{'low_coverage_sites'} = $self->find_low_coverage_sites;
        $bam_info{'n_double_peaks_sites'} = scalar keys %{$self->find_double_peaks_sites};
        $bam_info{'double_peaks_sites'} = $self->find_double_peaks_sites;
        $bam_info{'percent_sequence_coverage'} = sprintf("%.2f",$percent_sequence_coverage);
        $bam_info{'percent_sequence_with_high_coverage'} = sprintf("%.2f",$percent_sequence_high_coverage);
        # by feature
        if ($self->ref_genbank_file){
            $bam_info{'double_peaks_sites_by_feature'} = $self->find_double_peaks_sites_by_feature;
            my $n_double_peaks_sites_by_feature = $self->find_double_peaks_sites_by_feature(1);
            $bam_info{'n_double_peaks_sites_by_feature'} = $n_double_peaks_sites_by_feature;
            $bam_info{'low_coverage_sites_by_feature'} = $self->find_low_coverage_sites_by_feature;
            $bam_info{'n_low_coverage_sites_by_feature'} = $self->find_low_coverage_sites_by_feature(1);

            foreach my $region (keys %{$n_double_peaks_sites_by_feature}) {
                if ($region =~ /env/i){
                    my $dbp = $n_double_peaks_sites_by_feature->{$region};
                    if ($dbp > 0){
                        $bam_info{'classification'} = 'BAD';
                        $bam_info{'reason'} = 'Env double peaks';
                    }
                }
            }
        }

        unless ($bam_info{'classification'}){
            if ($bam_info{'n_double_peaks_sites'} > 10){
                $bam_info{'classification'} = 'BAD';
                $bam_info{'reason'} = 'double peaks';
            }
            elsif ($self->n_reads < $self->min_reads){
                $bam_info{'classification'} = 'BAD';
                $bam_info{'reason'} = "Less than ". $self->min_reads." reads";
            }
            else {
                $bam_info{'classification'} = 'GOOD';
                $bam_info{'reason'} = '';
            }
        }

        return \%bam_info;
    }

    __PACKAGE__->meta->make_immutable;
}

