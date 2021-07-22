#===============================================================================
# STEP 2: Denovo Assembly using De Brujn Graph (SPades)
#===============================================================================
import os, re, json

# Globals ---------------------------------------------------------------------
configfile: 'config.yml'
CWD = os.getcwd()

OUTPUTDIR = config['OUTPUTDIR']

DATA_DIR = config['DATADIR']

JSON_INPUT = OUTPUTDIR + '/STEP1/JSON_STEP1.json'

JSON_FILE = OUTPUTDIR + '/STEP2/JSON_STEP2.json'


FASTQ_FILES = list()
if os.path.isfile(JSON_INPUT):
    with open(JSON_INPUT) as json_file:
        data = json.load(json_file)
        if('passed' not in data.keys()):
            output_html = CWD + '/results/STEP1/summary.html'
            output_csv = CWD + '/results/STEP1/all_seq_table.csv'
            merged_json = CWD + '/results/STEP1/MergedJSONException.json'
            script_path = "'" + CWD + "/scripts/generate-summary.Rmd" + "'"

            merge_script_cmd = "PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl merge_all_json --results_dir " + CWD + "/results --jsonfile " + merged_json + " > " + output_csv
            print(merge_script_cmd)
            os.system(merge_script_cmd)

            report_cmd = 'Rscript -e " thetitle='+ "'HIVA Report'; rmarkdown::render(" + script_path + " , params=list(all_seq_table_csv = '" + output_csv + "' ), output_file = '" + output_html + "'" +')" '
            print(report_cmd)
            os.system(report_cmd)
            sys.exit(1)

        else:
#             for sample in data['passed']:
                # for pair in data['passed'][sample]['STEP1']['output']:
                    # for file in data['passed'][sample]['STEP1']['output'][pair].keys():
                        # FASTQ_FILES.append( file )


            for sample in data['passed']:
                for f in data['passed'][sample]:
                    for pair in data['passed'][sample][f]['STEP1']['output']:
                        file = data['passed'][sample][f]['STEP1']['output'][pair]['file']
                        FASTQ_FILES.append(file)


prog=re.compile('_R\d+')

SPADES_FILES = list(
                set(
                    [
                        CWD + "/" + OUTPUTDIR + '/STEP2/spades/' +
                        os.path.basename(prog.sub('', file).split('.')[0]) +
                        '/scaffolds.fasta' for file in FASTQ_FILES
                    ]
                )
            )

# Rules -----------------------------------------------------------------------
rule run_all:
    input: [ SPADES_FILES ]
    output:
        #samples = OUTPUTDIR + "/STEP2/samples.json",
        #low_reads = OUTPUTDIR + "/STEP2/low_reads_samples.json"
        json_out = OUTPUTDIR + '/STEP2/JSON_STEP2.json'
    params:
        dir = OUTPUTDIR + '/STEP2/spades'
    shell:
        """
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib  scripts/HIVAtoolkit.pl create_json_spades \
            -i {params.dir} \
            --min_scaffolds {config[SPADES_MIN_SCAFFOLDS]} \
            --min_length 1000 \
            --file_extension fasta \
            --jsonfile {output.json_out} \
            --step 2
        """

rule run_spades:
    input:
        #fastq = expand(CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{{prefix}}_{{suffix}}/{{prefix}}_R{pair}_{{suffix}}.fastq.gz", pair=['1','2']),
        fastq = expand(CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{{prefix}}_{{suffix}}/{{prefix}}_R{pair}_{{suffix}}.fastq.gz", pair=['1','2']),
        merged = CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_{suffix}_merged.fastq.gz"
    output: 
        file = CWD + "/" + OUTPUTDIR + "/STEP2/spades/{prefix}_{suffix}/scaffolds.fasta",
        clean_file = CWD + "/" + OUTPUTDIR + "/STEP2/spades/{prefix}_{suffix}/clean_scaffolds.fasta"
    log: 
        CWD+"/"+OUTPUTDIR+"/STEP2/spades/{prefix}_{suffix}_spades.log"
    params: 
        dir = CWD+"/"+OUTPUTDIR+"/STEP2/spades/{prefix}_{suffix}",
        dir_error = CWD+"/"+OUTPUTDIR+"/STEP2/spades_error/{prefix}_{suffix}"
    threads: 4
    run:
        shell("""
            set +e
            spades.py --tmp-dir /tmp --phred-offset 33 --careful -k 21,33,55,77 --only-assembler -o {params.dir} --pe1-1 {input.fastq[0]} --pe1-2 {input.fastq[1]} --merged {input.merged} -t {threads} &> {log}
            if [[ ! -e {output.file} ]];then
                touch {output.file}
                touch {output.clean_file}
            else
                PERL5LIB=$PERL5LIB:$PWD/scripts/lib  scripts/HIVAtoolkit.pl clean_fasta -i {output.file} -o {output.clean_file}
            fi
            find {params.dir} ! -wholename '{output.file}' -and ! -wholename '{output.clean_file}'  -type f -print0 | xargs -0 rm -f
            find {params.dir} ! -wholename '{params.dir}' -type d -print0 | xargs -0 rm -rf

            exitcode=$?
            if [ $exitcode -eq 1 ];then
                exit 0
                if [[ ! -e {output.file} ]];then
                    touch {output.file}
                fi
            else
                exit 0
            fi
        """)

