#===============================================================================
# STEP 3: Search for the closest references using denovo contigs
#===============================================================================
import os, re, json

# Globals ---------------------------------------------------------------------
configfile: 'config.yml'
CWD = os.getcwd()

OUTPUTDIR = config['OUTPUTDIR']

DATA_DIR = config['DATADIR']

SAMPLES_JSON = OUTPUTDIR + '/STEP2/JSON_STEP2.json'

# BLAST database
BLAST_NT_DB = config['NT_DATABASE_PATH']
BLAST_SCRIPT_PATH = config['BLAST_SCRIPT_PATH']


JSON_FILE = OUTPUTDIR + '/STEP3/JSON_STEP3.json'

SPADES_FILES = list()
if os.path.isfile(SAMPLES_JSON):
    with open(SAMPLES_JSON) as json_file:
        data = json.load(json_file)
        if('passed' not in data.keys()):
            output_html = CWD + '/results/STEP2/summary.html'
            output_csv = CWD + '/results/STEP2/all_seq_table.csv'
            merged_json = CWD + '/results/STEP2/MergedJSONException.json'
            script_path = "'" + CWD + "/scripts/generate-summary.Rmd" + "'"

            merge_script_cmd = "PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl merge_all_json --results_dir " + CWD + "/results --jsonfile " + merged_json + " > " + output_csv
            print(merge_script_cmd)
            os.system(merge_script_cmd)

            report_cmd = 'Rscript -e " thetitle='+ "'HIVA Report'; rmarkdown::render(" + script_path + " , params=list(all_seq_table_csv = '" + output_csv + "' ), output_file = '" + output_html + "'" +')" '
            print(report_cmd)
            os.system(report_cmd)
            sys.exit(1)
        else:
            for sample in data['passed']:
                for file in data['passed'][sample]:
                    f = data['passed'][sample][file]['STEP2']['output']['file']
                    if re.match('.*clean_scaffolds.*', f):
                        SPADES_FILES.append(f)

prog=re.compile('.*\/(\S+)\/clean_scaffolds.fasta')

REFERENCE_FILES = list(
                set(
                    [
                        CWD + "/" + OUTPUTDIR + '/STEP3/annot/' +
                        prog.sub(r'\1', file) +
                        '.fasta' for file in SPADES_FILES
                    ]
                )
            )

# Rules -----------------------------------------------------------------------
rule run_all:
    input: [ REFERENCE_FILES ]
    output:
        #samples = OUTPUTDIR + "/STEP3/samples.json", # reference found
        #low_reads = OUTPUTDIR + "/STEP3/low_reads_samples.json" # no reference found
        json_out = JSON_FILE
    params:
        dir = OUTPUTDIR + '/STEP3',
        json_out = JSON_FILE
    shell:
        """
                PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl create_json_reference \
                -i {params.dir} \
                --min_reads {config[MIN_HIV_REFERENCES]} \
                --file_extension fasta \
                --jsonfile {params.json_out} \
                --step 3

        """

rule check_orientation:
    input: CWD+"/"+OUTPUTDIR+"/STEP2/spades/{prefix}_{suffix}/clean_scaffolds.fasta"
    output:
        fasta = CWD+"/"+OUTPUTDIR+"/STEP3/annot/{prefix}_{suffix}.fasta",
        report = CWD+"/"+OUTPUTDIR+"/STEP3/annot/{prefix}_{suffix}_blast_report.txt"
    threads: 2
    log: CWD+"/"+OUTPUTDIR+"/STEP3/annot/{prefix}_{suffix}.log" 
    params:
        blast_script = BLAST_SCRIPT_PATH,
        blast_db = BLAST_NT_DB,
    shell:
        """
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl check_orientation --input_file {input} --report_file {output.report} --output_file {output.fasta} --reference_file {params.blast_db}  &> {log}

        if [[ ! -e {output.fasta} ]];then
            touch {output.fasta}
        fi
        """
