package MyApp::CheckDuplicatedFiles;
use feature qw(say);
use MooseX::App::Command;
extends 'MyApp';    # inherit log
use MooseX::FileAttribute;
use Cwd 'abs_path';
use namespace::autoclean;

command_short_description q[Checks md5sum file to duplicated files.];


has_file 'md5sum_file' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    required      => 1,
    documentation => q[Input MD5SUM file.],
);


has_file 'duplicated_files_list' => (
    traits        => ['AppOption'],
    cmd_type      => 'option',
    required      => 1,
    documentation => q[Output file where duplicated files should be written.]
);


sub run {
    my ($self) = @_;
    my %h;
    open( my $in, '<', $self->md5sum_file) 
        || die "Cannot open/read file " . $self->md5sum_file . "!";
    while ( my $row = <$in> ){
        chomp $row;
        my ($md5sum, $file) = split /\s+/, $row;
        push @{$h{$md5sum}}, $file;
    }
    close( $in );

    my @duplicated_files;
    foreach my $md5sum (sort keys %h) {
        if (scalar @{$h{$md5sum}} > 1){
            foreach my $file (@{$h{$md5sum}}) {
                # check if files are empty
                my $cmd = 'gunzip -c ' . abs_path($file) . ' | wc -l';
                my $n_lines = qx($cmd);
                push @duplicated_files, "$md5sum  $file  $n_lines" if $n_lines > 0;
            }
        }
    }
    if (@duplicated_files){
        open( my $out, '>', $self->duplicated_files_list )
            || die "Cannot open/write file " . $self->duplicated_files_list . "!";
        foreach my $file (@duplicated_files) {
            say $out $file;
        }
        close( $out );
    }

}

__PACKAGE__->meta->make_immutable;
1;
