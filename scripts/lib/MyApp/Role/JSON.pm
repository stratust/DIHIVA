package MyApp::Role::JSON;
use feature qw(say);
use MooseX::App::Role;
use MooseX::FileAttribute;
use JSON::MaybeXS;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);

requires 'log';

has_file 'json_file' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    cmd_flag      => 'jsonfile',
    required      => 1,
    documentation => q[Output JSON file for samples containing LESS than the minimum number of reads],
);

has_directory 'input_directory' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    cmd_aliases   => [qw(i)],
    required      => 0,
    documentation => q[FASTQ directory]
);

has 'file_extension' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    is            => 'rw',
    isa           => 'Str',
    required      => 0,
    default       => 'fastq.gz',
    documentation => q[File extension to search]
);

has 'min_reads' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    is            => 'rw',
    isa           => 'Int',
    required      => 0,
    documentation => q[Mininum number of reads required.]
);

has 'current_step' => (
    traits   => ['AppOption'],
    cmd_type => 'option',
    cmd_flag => 'step',
    is       => 'rw',
    isa      => 'Int',
    required => 0
);


sub export_json {
    my ($self, $json_hash) = @_;
    my $json = JSON::MaybeXS->new(utf8 => 1, pretty => 1, sort_by => 1);
    my $json_hash_export = $json->encode($json_hash);

    open(my $out, '>', $self->json_file)
        || die "Cannot open/write file " . $self->json_file . "!";

    say $out $json_hash_export;

    close($out);
}

1;
