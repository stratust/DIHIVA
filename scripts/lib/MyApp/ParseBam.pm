package MyApp::ParseBam;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use HIV::Assembly::Bam;
use namespace::autoclean;

command_short_description q[This command parses a HIV Assembly in BAM format and generate a consensus fasta file.];


has_file 'bam_file' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    required      => 1,
    documentation => q[Input BAM file.],
);


has_file 'output_consensus_file' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    required      => 1,
    documentation => q[Consensus fasta file.]
);


sub run {
    my ($self) = @_;
    my %bam_info;
    my $filename = $self->bam_file->basename;

    # remove extension
    $filename =~ s/^(.*)\.\S+$/$1/g;

    my $bam = HIV::Assembly::Bam->new(
        bam              => $self->bam_file,
        dp_threshold     => 0.7,
        minimum_coverage => 10,
        min_reads        => 500,
        no_ambiguities   => 1,
    );

    $bam->consensus_as_fasta_file($self->output_consensus_file->stringify);
}

__PACKAGE__->meta->make_immutable;
1;
