package MyApp::MergeJSONStep6;
use MooseX::App::Command;
use MooseX::FileAttribute;
extends 'MyApp';
with 'MyApp::Role::JSON';
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);
use JSON::MaybeXS;
use File::JSON::Slurper qw/ read_json write_json  /;
use Data::Printer;
use Hash::Merge qw/ merge /;
use namespace::autoclean;

has_directory 'json_directory' => (
    traits      => ['AppOption'],
    cmd_type    => 'option',
    cmd_flag    => 'json_directory',
    required    => 1,
    documentation => 'q[]'
);


sub merge_json {
    my ($self, $json_files_ref) = @_;

    my $first_json = shift(@{$json_files_ref});
    my $return_merged_json = read_json($first_json);


    foreach my $json_file (@{$json_files_ref}){
        my $ref = read_json($json_file);
        $return_merged_json =  merge($return_merged_json, $ref);
    }

    return($return_merged_json);
}


sub run {
    my ($self) = @_;
    my @json_files = glob($self->json_directory->stringify."/*.json");
    my $merged_json;

    $merged_json = $self->merge_json(\@json_files);
    $self->export_json($merged_json);

}
