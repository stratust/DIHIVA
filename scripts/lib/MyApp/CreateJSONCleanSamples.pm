package MyApp::CreateJSONCleanSamples;
use MooseX::App::Command;
extends 'MyApp';
with 'MyApp::Role::JSON';
use namespace::autoclean;

command_short_description
    q[Create a JSON file from all fastq files after removing contaminats.];


sub filter_samples_pairs {

    my ($self, @fastq_files) = @_;
    my %json_hash;
    my $step = 'STEP' . $self->current_step;
    foreach my $file (@fastq_files) {

        # get sample ID from file
        my ($sample_id, $pair, $filename, $well);
        if ($file =~ /^.*\/((\S+)\_(S\d+)\_.*(R[12])\_.*)/) {
            $filename  = $1;
            $sample_id = $2;
            $pair      = $4;
        } else {
            die
                "Cannot parse filename and retrive sample_id and pair information: $file";
        }

        my $reads = 0;
        my $cmd   = 'gunzip -c ' . $file . '| wc -l';
        my $lines = qx/$cmd/;
        chomp $lines;

        $reads = $lines / 4;

        # check if $reads contain a integer number
        unless ($reads =~ /^-?\d+\z/) {
            die
                "Error, FASTQ file might be corrupted! Number of reads ($reads) is not a integer number of reads!";
        }

        my $shortname;
        if ($filename =~ /^(.*)_(R[12]_)(.*)\.fastq\.gz/) {
            $shortname = $1 . '_' . $3;
        }

        my @k =  qw/step_failed classification reason/;

        my %aux = map {$_ => '' } @k;

        my $status = 'passed';

        if ( $reads < $self->min_reads ) {
            %aux = (
                'step_failed'    => 1,
                'classification' => "EMPTY WELL",
                'reason'         => "not enough reads after filtering"
                );

            $status = 'failed';
        }
        $json_hash{ $status }{ $sample_id }{ $shortname }{ $step }{ 'input' }{ $pair }{ file }   = '/data/' . $filename;
        $json_hash{ $status }{ $sample_id }{ $shortname }{ $step }{ 'input' }{ $pair }{ reads }  = '';
        $json_hash{ $status }{ $sample_id }{ $shortname }{ $step }{ 'output' }{ $pair }{ file }  = $file;
        $json_hash{ $status }{ $sample_id }{ $shortname }{ $step }{ 'output' }{ $pair }{ reads } = $reads;
        my $ref = $json_hash{ $status }{ $sample_id }{ $shortname }{ $step };
        map {$ref->{$_} = $aux{$_}} @k;

    }
    $self->export_json(\%json_hash);
}

sub run {
    my ($self) = @_;
    my @fastq_files;

    # Get subdirs
    my @dirs = glob($self->input_directory . '/*');
    foreach my $dir (@dirs) {
        if (-d $dir) {
            my @aux_fastq_files = glob($dir . '/*_R[12]_*' . $self->file_extension);
            push @fastq_files, @aux_fastq_files if @aux_fastq_files;
        }
    }
    $self->filter_samples_pairs(@fastq_files);
}

__PACKAGE__->meta->make_immutable;
1;
