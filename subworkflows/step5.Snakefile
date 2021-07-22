#===============================================================================
# STEP 5: Assembly Quality Evalutation
#===============================================================================
import glob, os, re, json

# Globals ---------------------------------------------------------------------
configfile: 'config.yml'

CWD = os.getcwd()

OUTPUTDIR = config['OUTPUTDIR']

DATA_DIR = config['DATADIR']

SAMPLES_JSON = OUTPUTDIR + '/STEP4/JSON_STEP4.json'

JSON_FILE = OUTPUTDIR + '/STEP5/JSON_STEP5.json'

GENBANK_FILES = list()
if os.path.isfile(SAMPLES_JSON):
    with open(SAMPLES_JSON) as json_file:
        data = json.load(json_file)
        for sample in data['passed']:
            for file in data['passed'][sample]:
                for f in data['passed'][sample][file]['STEP4']['output']:
                    GENBANK_FILES.append(f)

prog=re.compile('.*/(\S+).gb')

CLASSIFICATION_FILES = list(
                set(
                    [
                        CWD + "/" + OUTPUTDIR + '/STEP5/classified/log/' +
                        prog.sub(r'\1', file) + '.txt' for file in GENBANK_FILES
                    ]
                )
            )


# Rules -----------------------------------------------------------------------
rule run_all:
    input: [ CLASSIFICATION_FILES ]
    output:
        #samples = OUTPUTDIR + "/STEP5/samples.json", # reference found
        #low_reads = OUTPUTDIR + "/STEP5/low_reads_samples.json" # no reference found
        json_out = JSON_FILE
    params:
        dir = OUTPUTDIR + '/STEP5'
    shell:
        """
                PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl  create_json_good_assemblies \
                -i {params.dir} \
                --file_extension gb \
                --step 5 \
                --jsonfile {output.json_out}
        """


rule call_quality_bam:
    input:
        bam=CWD+"/"+OUTPUTDIR+'/STEP4/minimap2_final_assembly/{sample}/{sample}.bam',
        unpadded_fasta=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{sample}/{sample}.fasta",
        unpadded_gb=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{sample}/{sample}.gb"
    output:
        cons_fasta=CWD+"/"+OUTPUTDIR+"/STEP5/consensus_fasta/{sample}_cons.fasta",
        json=CWD+"/"+OUTPUTDIR+"/STEP5/quality_calls/{sample}.json"
    threads: 1
    #log: 
        #analyze_log=CWD+"/"+OUTPUTDIR+"/STEP5/consensus_fasta/{sample}_analyze.log"
        #draft_log=CWD+"/"+OUTPUTDIR+"/STEP5/consensus_fasta/{sample}_draft.log"
    shell:
        """
        samtools index {input.bam}
 
        PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl classify_bam \
            --bam_file {input.bam} \
            --reference_genbank_file {input.unpadded_gb} \
            --output_consensus_file {output.cons_fasta} \
            --step 5 \
            --jsonfile {output.json}

       """


rule annotate_consensus:
    input:
        cons_fasta=CWD+"/"+OUTPUTDIR+"/STEP5/consensus_fasta/{sample}_cons.fasta",
    output:
        cons_gb=CWD+"/"+OUTPUTDIR+"/STEP5/consensus_fasta/{sample}_cons.gb"
    params:
        hxb2_genbank_file=config['HXB2_ANNOTATED_REF']
    log: 
        analyze_log=CWD+"/"+OUTPUTDIR+"/STEP5/consensus_fasta/{sample}_analyze.log"
    shell:
        """
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl annotate_assembly_blast \
            --reference_file {params.hxb2_genbank_file} \
            --input_file {input.cons_fasta} \
            --output_genbank {output} &> {log}
        """


rule classify_quality_assembly:
    input:
        bam=CWD+"/"+OUTPUTDIR+'/STEP4/minimap2_final_assembly/{sample}/{sample}.bam',
        unpadded_gb=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{sample}/{sample}.gb",
        json=CWD+"/"+OUTPUTDIR+'/STEP5/quality_calls/{sample}.json',
        cons_gb=CWD+"/"+OUTPUTDIR+"/STEP5/consensus_fasta/{sample}_cons.gb"
    params:
        projectname='{sample}',
        good_dir=CWD+"/"+OUTPUTDIR+"/STEP5/classified/good",
        bad_dir=CWD+"/"+OUTPUTDIR+"/STEP5/classified/bad"
    output:
        classfied_file=CWD+"/"+OUTPUTDIR+'/STEP5/classified/log/{sample}.txt'
    shell:
        """
        export quality=$(perl -ane 'print $1 if /.*classification.*?(\w+).*/' {input.json})

        if [[ ! -e {params.good_dir} ]]; then
            mkdir -p {params.good_dir}/assembly; mkdir -p {params.bad_dir}/assembly;
            mkdir -p {params.good_dir}/genbank_consensus; mkdir -p {params.bad_dir}/genbank_consensus;
        fi

        if [[ $quality == "GOOD" ]]; then
            cp {input.bam} {params.good_dir}/assembly
            cp {input.unpadded_gb} {params.good_dir}/assembly/{params.projectname}.gb
            cp {input.cons_gb} {params.good_dir}/genbank_consensus/{params.projectname}.gb
        else
            cp {input.bam} {params.bad_dir}/assembly
            cp {input.unpadded_gb} {params.bad_dir}/assembly/{params.projectname}.gb
            cp {input.cons_gb} {params.bad_dir}/genbank_consensus/{params.projectname}.gb
        fi
        touch {output.classfied_file}
        """
