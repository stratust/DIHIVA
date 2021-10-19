package MyApp::ChooseGBFile;
use MooseX::App::Command;
use MooseX::FileAttribute;
extends 'MyApp';
with 'MyApp::Role::JSON';
use JSON::MaybeXS;
use File::JSON::Slurper qw/ read_json write_json /;
use Data::Printer;
use File::Basename;
use File::Copy;
use namespace::autoclean;

has_file 'json_step6' => (
    traits   => ['AppOption'],
    cmd_type => 'option',
    required => 1
);


sub create_aux_hash {

    my ( $self, $aux ) = @_;
    my %classif_hash;
    my @samples;

    @samples = sort {$a cmp $b} keys %{$aux};
    foreach my $sample (@samples) {

        my $classif = $aux->{$sample}{'final_classification'};
        push @{ $classif_hash{$classif} }, $sample;

    }

    return(\%classif_hash);

} 

sub select_genbank { 

    my ($self, $json_step6_ref) = @_;
    my @samples;
    my $ref_classif;
    my %final_json;
    my @classif_priorities = ('intact', 'msd_mutation', 'non_functional', 'structural_variations', 'missing_internal_genes', 'missing_ltr_or_psi');
    @samples = keys %{$json_step6_ref->{'passed'}};

    foreach my $sample (@samples) {
        my @sub_hash_keys = keys %{$json_step6_ref->{'passed'}{$sample}{'STEP6'}} ;

        $ref_classif = $self->create_aux_hash($json_step6_ref->{'passed'}{$sample}{'STEP6'});
        my $not_chosen_yet = 1;
        foreach my $classification (@classif_priorities){

            if ($ref_classif->{$classification} and $not_chosen_yet){
                # remove the chosen sample from the hash
                my $chosen_sample = shift @{$ref_classif->{$classification}};

                my $output_step7 = $json_step6_ref->{'passed'}{$sample}{'STEP6'}{$chosen_sample}{'output_gb_classification'};
                $output_step7 =~ s/STEP6/STEP7/g;

                $final_json{'passed'}{$sample}{'STEP7'}{'chosen'}{$chosen_sample} = $json_step6_ref->{'passed'}{$sample}{'STEP6'}{$chosen_sample};
                $final_json{'passed'}{$sample}{'STEP7'}{'chosen'}{$chosen_sample}{'output_step7'} = $output_step7;

                # unset not_chosen_yet flag
                $not_chosen_yet = 0;
            }

            if ($ref_classif->{$classification}){
                foreach my $not_chosen_sample ( @{$ref_classif->{$classification}} ) {
                    $final_json{'passed'}{$sample}{'STEP7'}{'not_chosen'}{$not_chosen_sample} = $json_step6_ref->{'passed'}{$sample}{'STEP6'}{$not_chosen_sample};
                }
            }
        }
    }

    return(\%final_json);

}

sub move_gb_to_folders {

    my ($self, $json_step6_ref) = @_;
    my ($dirname, $step6_output, $step7_output, $sample_classif, @samples, @chosen_file, %seen );

    @samples = sort {$a cmp $b} keys %{$json_step6_ref->{'passed'}} ;
    foreach my $sample (@samples) {

        @chosen_file = keys %{$json_step6_ref->{'passed'}{$sample}{'STEP7'}{'chosen'}};
        foreach my $file (@chosen_file) {

            $step6_output = $json_step6_ref->{'passed'}{$sample}{'STEP7'}{'chosen'}{$file}{'output_gb_classification'};
            $step7_output = $json_step6_ref->{'passed'}{$sample}{'STEP7'}{'chosen'}{$file}{'output_step7'};
            $sample_classif =  $json_step6_ref->{'passed'}{$sample}{'STEP7'}{'chosen'}{$file}{'final_classification'};

            if ($seen{$sample_classif}) {
                copy($step6_output, $step7_output);

            } else {
                $dirname = dirname($step7_output);
                mkdir($dirname);
                copy($step6_output, $step7_output);
                $seen{$sample_classif} = 1;
            }
        }
    }
}

sub run {
    my ($self) = @_;
    my $final_json;
    my $json_step6_ref = read_json($self->json_step6->stringify);

    $final_json = $self->select_genbank($json_step6_ref);
    $self->move_gb_to_folders($final_json);
    $self->export_json($final_json);
}
