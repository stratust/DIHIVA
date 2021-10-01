#===============================================================================
# STEP 6: Check productivity for Near Full Lenght sequences
#===============================================================================
import glob, os, re, json

# Globals ---------------------------------------------------------------------
configfile: 'config.yml'

CWD = os.getcwd()

OUTPUTDIR = config['OUTPUTDIR']

DATA_DIR = config['DATADIR']

SAMPLES_JSON = OUTPUTDIR + '/STEP5/JSON_STEP5.json'
JSON_FILE = OUTPUTDIR + '/STEP6/JSON_STEP6.json'

GENBANK_FILES = glob.glob(OUTPUTDIR + '/STEP5/classified/good/genbank_consensus/*.gb')

if len(GENBANK_FILES) == 0:
    output_html = CWD + '/results/STEP5/summary.html'
    output_csv = CWD + '/results/STEP5/all_seq_table.csv'
    merged_json = CWD + '/results/STEP5/MergedJSONException.json'
    script_path = "'" + CWD + "/scripts/generate-summary.Rmd" + "'"

    merge_script_cmd = "PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl merge_all_json --results_dir " + CWD + "/results --jsonfile " + merged_json + " > " + output_csv
    print(merge_script_cmd)
    os.system(merge_script_cmd)

    report_cmd = 'Rscript -e " thetitle='+ "'HIVA Report'; rmarkdown::render(" + script_path + " , params=list(all_seq_table_csv = '" + output_csv + "' ), output_file = '" + output_html + "'" +')" '
    print(report_cmd)
    os.system(report_cmd)
    sys.exit(1)


PROD_FILES = list(
                set(
                    [
                        CWD + "/" + OUTPUTDIR + '/STEP6/productivity_check/' +
                        os.path.basename(file).split('.')[0] +
                        '.txt' for file in GENBANK_FILES
                    ]
                )
            )


# Rules -----------------------------------------------------------------------
rule run_all:
    input: [ PROD_FILES ]
    output:
        JSON_FILE
    shell:
        """ 
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl merge_json_step6 \
            --json_directory results/STEP6/json_samples \
            --jsonfile {output}
        """


rule check_productivity:
    input:
        dp_call=CWD+"/"+OUTPUTDIR+'/STEP5/quality_calls/{sample}.json',
        classfied_file=CWD+"/"+OUTPUTDIR+'/STEP5/classified/log/{sample}.txt'
    params:
        genbank=CWD+"/"+OUTPUTDIR+"/STEP5/classified/good/genbank_consensus/{sample}.gb",
        hg_genbank=CWD+"/"+OUTPUTDIR+"/STEP5/classified/good_hg/genbank_consensus/{sample}.gb",
        script_path=config['BLAST_SCRIPT_PATH'],
        intact_dir=CWD+"/"+OUTPUTDIR+"/STEP6/intact",
        nonfunctional_dir=CWD+"/"+OUTPUTDIR+"/STEP6/non_functional",
        msd_mut_dir=CWD+"/"+OUTPUTDIR+"/STEP6/msd_mutation",
        missing_internal_genes_dir=CWD+"/"+OUTPUTDIR+"/STEP6/missing_internal_genes",
        output_dir=CWD+"/"+OUTPUTDIR+"/STEP6",
        json_dir=CWD+"/"+OUTPUTDIR+"/STEP6/json_samples",
        json_file=CWD+"/"+OUTPUTDIR+"/STEP6/json_samples/{sample}.json"
    output:
        prod_file=CWD+"/"+OUTPUTDIR+"/STEP6/productivity_check/{sample}.txt"
    shell:
        """
        if [[ ! -e {params.intact_dir} ]]; then
            mkdir {params.json_dir}
            mkdir {params.intact_dir}
            mkdir {params.msd_mut_dir}
            mkdir {params.nonfunctional_dir}
            mkdir {params.missing_internal_genes_dir}
        fi
        if [[ -e {params.genbank} ]]; then

            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl sequences_classification \
               --genbank_file {params.genbank} \
               --output_dir {params.output_dir} \
               --step 6 \
               --jsonfile {params.json_file} \
               &> {output.prod_file} \

        else
            echo "bad assembly, annotations not checked" > {output.prod_file}
        fi
        """

