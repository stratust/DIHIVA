import os, json, re

# Globals ---------------------------------------------------------------------
configfile: 'config.yml'
CWD = os.getcwd()

OUTPUTDIR = config['OUTPUTDIR']

DATA_DIR =  config['DATADIR']

# Intitial list of samples
SAMPLES_JSON = OUTPUTDIR + '/JSON_STEP0.json'

if not os.path.isdir(OUTPUTDIR):
    os.mkdir(OUTPUTDIR)

if not os.path.isfile(SAMPLES_JSON):
    os.system(
       ' PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl create_json -i ' + DATA_DIR +
       ' --min_reads ' + str(config['MIN_INITIAL_READ_PAIRS']) +
       ' --jsonfile ' + SAMPLES_JSON  +
       ' --step 0'
       )

FASTQ_FILES = list()
if os.path.isfile(SAMPLES_JSON):
    with open(SAMPLES_JSON) as json_file:
        data = json.load(json_file)
        if('passed' not in data.keys()):
            output_html = CWD + '/results/summary.html'
            output_csv = CWD + '/results/all_seq_table.csv'
            merged_json = CWD + '/results/MergedJSONException.json'
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
                for f in data['passed'][sample]:
                    for pair in data['passed'][sample][f]['STEP0']['input']:
                        file = data['passed'][sample][f]['STEP0']['input'][pair]['file']
                        FASTQ_FILES.append(file)

else:
    print("** Please genenerate samples.json\n")
    sys.exit(1)

prog=re.compile('_R\d+')


MERGED_FILES = list(
                set(
                    [
                        CWD + "/" + OUTPUTDIR + '/STEP1/bbmerge/' +
                        os.path.basename(prog.sub('', file).split('.')[0]) +'/'+os.path.basename(file) for file in FASTQ_FILES
                    ]
                )
            )


subworkflow step_zero:
    snakefile: 'step0.Snakefile'

# Rules -----------------------------------------------------------------------
rule run_all:
    input: MERGED_FILES
    output:
        #samples = OUTPUTDIR + '/STEP1/samples.json',
        #low_reads = OUTPUTDIR + '/STEP1/low_reads_samples.json'
        json_out = OUTPUTDIR + '/STEP1/JSON_STEP1.json'
    params:
        dir =  OUTPUTDIR + '/STEP1/remove_contaminants_with_bbduk' if config['BBDUK_HIV_REF'] else OUTPUTDIR + '/STEP1/remove_contaminants'
    shell:
        """
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl create_json_clean_samples -i {params.dir} \
                --min_reads {config[MIN_CLEAN_READ_PAIRS]} \
                --jsonfile {output.json_out} \
                --step 1
        """

#===============================================================================
# STEP 1: Quality Check
#===============================================================================

# Remove PCR amplification and dedup optical duplicates
rule clumpify_dedup:
    input: 
        expand(DATA_DIR+"/{{prefix}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2'])
    output:
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/{{prefix}}_R{pair}_{{suffix}}.dedupe.fastq.gz",pair=['1','2']),
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/{{prefix}}_R{pair}_{{suffix}}.corrected.fastq.gz",pair=['1','2']),
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/{{prefix}}_R{pair}_{{suffix}}.unsorted.fastq.gz",pair=['1','2']),
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/{{prefix}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2'])
    params:
        fq1=CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/{prefix}_R1_{suffix}.fastq",
        fq2=CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/{prefix}_R2_{suffix}.fastq"
    log: 
        CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/logs/{prefix}_{suffix}.log"
    threads: 3,
    resources:
        mem_mb=8000
    shell: 
        """
            _JAVA_OPTIONS="-XX:ParallelGCThreads=1 -XX:+UseParallelGC" clumpify.sh in1={input[0]} \
                in2={input[1]}  \
                out1={output[0]} \
                out2={output[1]} \
                seed=1 \
                dedupe=t \
                containment=t \
                k=21 \
                threads={threads} \
                -Xmx{resources.mem_mb}M \
                -da \
                > {log} 2>&1
            _JAVA_OPTIONS="-XX:ParallelGCThreads=1 -XX:+UseParallelGC" clumpify.sh in1={output[0]} \
                in2={output[1]}  \
                out1={output[2]} \
                out2={output[3]} \
                seed=1 \
                k=21 \
                passes=6 \
                subs=5 \
                ecc=t \
                ecco=t \
                unpair=t \
                repair=t \
                threads={threads} \
                -Xmx{resources.mem_mb}M \
                -da \
                >> {log} 2>&1
            _JAVA_OPTIONS="-XX:ParallelGCThreads=1 -XX:+UseParallelGC" clumpify.sh in1={output[2]} \
                in2={output[3]}  \
                out1={output[4]} \
                out2={output[5]} \
                seed=1 \
                dedupe=t \
                containment=t \
                k=21 \
                threads={threads} \
                -Xmx{resources.mem_mb}M \
                -da \
                >> {log} 2>&1
            gunzip -c {output[4]} | fastq-sort --id  > {params.fq1}
            gunzip -c {output[5]} | fastq-sort --id  > {params.fq2}
            gzip -n "{params.fq1}" "{params.fq2}"

        """


rule trim_galore:
    input: 
        #expand(DATA_DIR+"/{{prefix}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2'])
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/clumpify/{{prefix}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2'])
    output:
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/trim_galore/{{prefix}}_R{pair}_{{suffix}}_val_{pair}.fq.gz",pair=['1','2'])
    params:
        output_dir=CWD+"/"+OUTPUTDIR+'/STEP1/trim_galore'
    log: 
        CWD+"/"+OUTPUTDIR+"/STEP1/trim_galore/{prefix}_{suffix}.log"
    threads: 2
    resources:
        mem_mb=8000
    shell: 
        """
            _JAVA_OPTIONS='' trim_galore -e 0.1 --paired --quality 20 --length 100 --stringency 1 --fastqc --dont_gzip --output_dir={params.output_dir} {input} > {log} 2>&1
            F1=$(echo {output[0]} | perl -plane 's/\.gz//g')
            F2=$(echo {output[1]} | perl -plane 's/\.gz//g')
            gzip -n "$F1" "$F2"
        """


rule remove_contaminants_with_bbduk:
    input:
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/adapter_trimmed/{{prefix}}_R{pair}_{{suffix}}_val_{pair}.fq.gz",pair=['1','2']) if os.path.isfile(CWD+'/adapter.fasta') else expand(CWD+"/"+OUTPUTDIR+"/STEP1/trim_galore/{{prefix}}_R{pair}_{{suffix}}_val_{pair}.fq.gz",pair=['1','2'])
    output:
        final_files = expand(CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{{prefix}}_{{suffix}}/{{prefix}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2']),
        pre_final_files = temp( expand(CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{{prefix}}_{{suffix}}/{{prefix}}_R{pair}_{{suffix}}_pre_final.fastq.gz",pair=['1','2']) )
    params:
        ref = config['BBDUK_HIV_REF'],
        univec = config['UNIVEC_DB'],
        adapters = config['ADAPTERS_DB'],
        primers = config['PRIMERS_DB'],
        primers_rev_comp = config['PRIMERS_REV_COMP_DB'],
        correct = config['CORRECTION'],
        common_reads_list = CWD + "/" + OUTPUTDIR + "/STEP1/remove_contaminants_with_bbduk/{prefix}_{suffix}/common_reads.list"
    log:
        hiv=CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{prefix}_{suffix}/logs/{prefix}_{suffix}_hiv.log",
        univec=CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{prefix}_{suffix}/logs/{prefix}_{suffix}_univec.log",
        adapters=CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{prefix}_{suffix}/logs/{prefix}_{suffix}_adapters.log",
        primers=CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{prefix}_{suffix}/logs/{prefix}_{suffix}_primers.log"
    threads: 7
    resources:
        mem_mb=3000
    shell:
        """
            F1=$(echo {output.final_files[0]} | perl -plane 's/\.gz//g')
            F2=$(echo {output.final_files[1]} | perl -plane 's/\.gz//g')

            F1_AUX=$(echo {output.final_files[0]} | perl -plane 's/\.fastq\.gz/\_aux\.fastq/g')
            F2_AUX=$(echo {output.final_files[1]} | perl -plane 's/\.fastq\.gz/\_aux\.fastq/g')
            F1_AUX2=$(echo {output.final_files[0]} | perl -plane 's/\.fastq\.gz/\_aux2\.fastq/g')
            F2_AUX2=$(echo {output.final_files[1]} | perl -plane 's/\.fastq\.gz/\_aux2\.fastq/g')
            F1_AUX3=$(echo {output.final_files[0]} | perl -plane 's/\.fastq\.gz/\_aux3\.fastq/g')
            F2_AUX3=$(echo {output.final_files[1]} | perl -plane 's/\.fastq\.gz/\_aux3\.fastq/g')
            F1_AUX4=$(echo {output.final_files[0]} | perl -plane 's/\.fastq\.gz/\_aux4\.fastq/g')
            F2_AUX4=$(echo {output.final_files[1]} | perl -plane 's/\.fastq\.gz/\_aux4\.fastq/g')
            F1_AUX5=$(echo {output.final_files[0]} | perl -plane 's/\.fastq\.gz/\_aux5\.fastq/g')
            F2_AUX5=$(echo {output.final_files[1]} | perl -plane 's/\.fastq\.gz/\_aux5\.fastq/g')


            bbduk.sh -Xmx{resources.mem_mb}M in1={input[0]} in2={input[1]}  outm1=$F1_AUX  outm2=$F2_AUX ref={params.ref} k=31 hdist=1 threads={threads} > {log.hiv} 2>&1
            bbduk.sh -Xmx{resources.mem_mb}M in1=$F1_AUX in2=$F2_AUX  out1=$F1_AUX2  out2=$F2_AUX2 ref={params.adapters} k=16 hdist=0 mink=15 ktrim='N' threads={threads} > {log.adapters} 2>&1
            bbduk.sh -Xmx{resources.mem_mb}M in1=$F1_AUX2 in2=$F2_AUX2  out1=$F1_AUX3  out2=$F2_AUX3 ref={params.univec} k=27 hdist=0 mink=15 ktrim='N' threads={threads} > {log.univec} 2>&1

            if [[ {params.correct} == 'True' ]];then

                cutadapt --action=trim --no-indels -g file:{params.primers} -O 10 -e 0.1  -o $F1_AUX4 $F1_AUX3 >> {log.primers} 2>&1
                cutadapt --action=trim --no-indels -g file:{params.primers} -O 10 -e 0.1  -o $F2_AUX4 $F2_AUX3 >> {log.primers} 2>&1

                cutadapt --action=trim --no-indels -a file:{params.primers_rev_comp} -O 10 -e 0.1  -o $F1_AUX5 $F1_AUX4 >> {log.primers} 2>&1
                cutadapt --action=trim --no-indels -a file:{params.primers_rev_comp} -O 10 -e 0.1  -o $F2_AUX5 $F2_AUX4 >> {log.primers} 2>&1

                seqtk sample -s {config[SAMPLING_SEED]} <( fastq-sort --id $F1_AUX5 ) {config[N_READS_SAMPLED]} > {output.pre_final_files[0]}
                seqtk sample -s {config[SAMPLING_SEED]} <( fastq-sort --id $F2_AUX5 ) {config[N_READS_SAMPLED]} > {output.pre_final_files[1]}

            else

                seqtk sample -s {config[SAMPLING_SEED]} <( fastq-sort --id $F1_AUX3 ) {config[N_READS_SAMPLED]}  > {output.pre_final_files[0]}
                seqtk sample -s {config[SAMPLING_SEED]} <( fastq-sort --id $F2_AUX3 ) {config[N_READS_SAMPLED]}  > {output.pre_final_files[1]}

            fi


            cat <( seqkit seq --name --only-id {output.pre_final_files[0]}  ) <( seqkit seq --name --only-id {output.pre_final_files[1]} ) | sort | uniq -c | perl -lane 'print $F[1] if $F[0] > 1 ' > {params.common_reads_list}

            seqtk subseq {output.pre_final_files[0]} {params.common_reads_list} > $F1
            seqtk subseq {output.pre_final_files[1]} {params.common_reads_list} > $F2


            gzip -n "$F1" "$F2"

            rm -rf $F1_AUX
            rm -rf $F1_AUX2
            rm -rf $F1_AUX3
            rm -rf $F2_AUX
            rm -rf $F2_AUX2
            rm -rf $F2_AUX3
            rm -rf $F1_AUX4
            rm -rf $F2_AUX4
            rm -rf $F1_AUX5
            rm -rf $F2_AUX5


        """


rule merge_paired_reads:
    input:
        expand(CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{{prefix}}_{{suffix}}/{{prefix}}_R{pair}_{{suffix}}.fastq.gz",pair=['1','2']),
    output:
        r1=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_R1_{suffix}.fastq.gz",
        r2=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_R2_{suffix}.fastq.gz",
        merged=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_{suffix}_merged.fastq.gz"
    params:
        sr1=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_R1_{suffix}.fastq",
        sr2=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_R2_{suffix}.fastq",
        smerged=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_{suffix}_merged.fastq",
        ur1=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_R1_{suffix}_unsorted.fastq.gz",
        ur2=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_R2_{suffix}_unsorted.fastq.gz",
        umerged=CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/{prefix}_{suffix}_merged_unsorted.fastq.gz"
    log:
        CWD+"/"+OUTPUTDIR+"/STEP1/bbmerge/{prefix}_{suffix}/logs/{prefix}_{suffix}.log"
    threads: 2
    resources:
        mem_mb=8000
    shell:
        """
            bbmerge-auto.sh -Xmx{resources.mem_mb}M in1={input[0]} in2={input[1]} out={params.umerged} outu1={params.ur1} outu2={params.ur2} threads={threads} pfilter=1 rem k=62 vstrict > {log} 2>&1
            gunzip -c {params.umerged} | fastq-sort --id  > {params.smerged}
            gunzip -c {params.ur1} | fastq-sort --id  > {params.sr1}
            gunzip -c {params.ur2} | fastq-sort --id  > {params.sr2}
            gzip -n "{params.smerged}" "{params.sr1}" "{params.sr2}"
            rm -rf {params.umerged} {params.ur1} {params.ur2}
        """
