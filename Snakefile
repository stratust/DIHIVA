import glob, os, re, json

configfile: 'config.yml'
CWD = os.getcwd()
OUTPUTDIR = config['OUTPUTDIR']

subworkflow step_zero:
    snakefile: 'subworkflows/step0.Snakefile'

subworkflow step_one:
    snakefile: 'subworkflows/step1.Snakefile'

subworkflow step_two:
    snakefile: 'subworkflows/step2.Snakefile'

subworkflow step_three:
    snakefile: 'subworkflows/step3.Snakefile'

subworkflow step_four:
    snakefile: 'subworkflows/step4.Snakefile'

subworkflow step_five:
    snakefile: 'subworkflows/step5.Snakefile'

subworkflow step_six:
    snakefile: 'subworkflows/step6.Snakefile'


rule all_subworkflows:
    input:
        step_zero(OUTPUTDIR + '/check_duplicated.out'),
        step_one(OUTPUTDIR + '/STEP1/JSON_STEP1.json'),
        step_two(OUTPUTDIR + '/STEP2/JSON_STEP2.json'),
        step_three(OUTPUTDIR + '/STEP3/JSON_STEP3.json'),
        step_four(OUTPUTDIR + '/STEP4/JSON_STEP4.json'),
        step_five(OUTPUTDIR + '/STEP5/JSON_STEP5.json'),
        CWD + '/' + OUTPUTDIR + '/STEP7/summary.html'


rule move_files_and_generate_summary:
    input:
        step_six(OUTPUTDIR + '/STEP6/JSON_STEP6.json')
    output:
        output_one = CWD + '/' + OUTPUTDIR + '/STEP7/summary.html',
        output_two = CWD + '/' + OUTPUTDIR + '/STEP7/all_seq_table.csv',
        tmpdir = temp(directory(CWD + '/' + OUTPUTDIR + '/STEP7/.tmp'))
    params:
        json_merged = CWD + '/' + OUTPUTDIR + '/STEP7/MergedJSON.json'
    shell:
        """
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl merge_all_json --results_dir results --jsonfile {params.json_merged} > {output.output_two}

            Rscript -e " thetitle='DIHIVA Report'; rmarkdown::render('scripts/generate-summary.Rmd', params=list(all_seq_table_csv = '{output.output_two}' ), output_file = '{output.output_one}', intermediates_dir = '{output.tmpdir}')  "
        """

