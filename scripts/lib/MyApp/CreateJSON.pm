package MyApp::CreateJSON;
use MooseX::App::Command;
extends 'MyApp';
with 'MyApp::Role::JSON';
use namespace::autoclean;

command_short_description
    q[Create a JSON file from all fastq files in a directory.];
command_long_description
    q[Create a JSON file from all fastq files in a directory.];

sub filter_samples_pairs {
    my ($self, @fastq_files) = @_;
    my %json_hash;

    foreach my $file (@fastq_files) {

        my ($sample_id, $pair, $filename, $well);
        if ($file =~ /^.*\/((\S+)\_(S\d+)\_.*(R[12])\_.*)/) {
            $filename  = $1;
            $sample_id = $2;
            $well      = $3;
            $pair      = $4;
        } else {
            die
                "Cannot parse filename and retrive sample_id and pair information: $file";
        }

        my $reads = 0;
        my $cmd   = 'gunzip -c ' . $file . '| wc -l';
        my $lines = qx/$cmd/;
        chomp $lines;

        $reads = $lines/4;

        # check if $reads contain a integer number
        unless ( $reads =~ /^-?\d+\z/ ) { 
            die "Error, FASTQ file might be corrupted! Number of reads ($reads) is not a integer number of reads!";
        }

        my $shortname;
        if ($filename =~ /^(.*)_(R[12]_)(.*)\.fastq\.gz/) {
            $shortname = $1 . '_' . $3;
        }

        if ($reads >= $self->min_reads) {
            $json_hash{'passed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'input'}{$pair}{file} = $file;
            $json_hash{'passed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'input'}{$pair}{reads} = $reads;
            $json_hash{'passed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'output'}{$pair}{file} = $file;
            $json_hash{'passed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'output'}{$pair}{reads} = $reads;
            $json_hash{'passed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'well'} = $well;
            $json_hash{'passed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'reason'} = '';
            $json_hash{'passed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'classification'} = '';
        } else {
            $json_hash{'failed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'well'} = $well;
            $json_hash{'failed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'input'}{$pair}{file} = $file;
            $json_hash{'failed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'input'}{$pair}{reads} = $reads;
            $json_hash{'failed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'output'}{$pair}{file} = $file;
            $json_hash{'failed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'output'}{$pair}{reads} = $reads;
            $json_hash{'failed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'reason'} = 'not enough reads';
            $json_hash{'failed'}{$sample_id}{$shortname}{'STEP' . $self->current_step}{'classification'} = 'EMPTY WELL';
        }

    }
    $self->export_json(\%json_hash);
}

sub run {
    my ($self) = @_;
    my @fastq_files = glob($self->input_directory . '/*.' . $self->file_extension);
    $self->filter_samples_pairs(@fastq_files);
}

__PACKAGE__->meta->make_immutable;
1;
