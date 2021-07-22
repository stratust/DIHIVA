package MyApp::ClassifyBam;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use HIV::Assembly::Bam;
use JSON::PP;
use Bio::SeqIO;
use namespace::autoclean;

command_short_description q[This command parses a HIV Assembly in BAM format.];
command_long_description q[This command parses a HIV Assembly in BAM format.];


has_file 'bam_file' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    required      => 1,
    documentation => q[Input BAM file.],
);


has_file 'reference_genbank_file' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    required      => 1,
    documentation => q[Annotated reference sequence in genbank format.],
);


has_file 'output_consensus_file' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    required      => 1,
    documentation => q[Consensus fasta file.]
);

has 'json' => (
    is            => 'ro',
    isa           => 'HashRef',
    documentation => 'JSON output file',
    default       => sub { return {} }
);

with 'MyApp::Role::JSON';

sub run {
    my ($self) = @_;
    my %bam_info;
    my $filename = $self->bam_file->basename;

    # remove extension
    $filename =~ s/^(.*)\.\S+$/$1/g;
    my $sample= $filename;
    $sample =~ s/(.*)\_S\d+\_.*/$1/g;

    my $bam = HIV::Assembly::Bam->new(
        bam              => $self->bam_file,
        dp_threshold     => 0.7,
        minimum_coverage => 10,
        min_reads        => 500,
        ref_genbank_file => $self->reference_genbank_file,
        primary_tag      => 'gene',
        tag              => 'standard_name',
        no_ambiguities   => 1,
    );

    my $bam_info = $bam->get_bam_info;
    my $status = 'passed';
    $status = 'failed' if $bam_info->{classification} eq 'BAD';
    $self->json->{$status}{$sample}{$filename}{'STEP' . $self->current_step}= $bam_info;

    my $json        = JSON::PP->new->ascii->pretty->allow_nonref;
    my $json_pretty = $json->encode($self->json);
    open(my $out, '>', $self->json_file)
        || die "Cannot open/write file " . $self->json_file . "!";
    say $out $json_pretty;
    close($out);
    
    $bam->consensus_as_fasta_file($self->output_consensus_file->stringify);
}

__PACKAGE__->meta->make_immutable;
1;
