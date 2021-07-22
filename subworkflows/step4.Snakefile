#===============================================================================
# STEP 4: Create Modified reference
#===============================================================================
import glob, os, re, json

# Globals ---------------------------------------------------------------------
configfile: 'config.yml'
CWD = os.getcwd()

OUTPUTDIR = config['OUTPUTDIR']

DATA_DIR = config['DATADIR']

SAMPLES_JSON = OUTPUTDIR + '/STEP3/JSON_STEP3.json'
JSON_FILE = OUTPUTDIR + '/STEP4/JSON_STEP4.json'

DENOVO = True if 'DENOVO' in config and config['DENOVO'] else False


REF_FILES = list()
if os.path.isfile(SAMPLES_JSON):
    with open(SAMPLES_JSON) as json_file:
        data = json.load(json_file)
        for sample in data['passed']:
            for file in data['passed'][sample]:
                for f in data['passed'][sample][file]['STEP3']['output']:
                    REF_FILES.append(f)


GENBANK_FILES = list(
                set(
                    [
                        CWD + "/" + OUTPUTDIR + '/STEP4/minimap2_final_assembly/' +
                        os.path.basename(file).split('.')[0] + '/' +
                        os.path.basename(file).split('.')[0] + '.gb' for file in REF_FILES
                    ]
                )
            )


# Rules -----------------------------------------------------------------------
rule run_all:
    input: [ GENBANK_FILES ]
    output:
        #samples = OUTPUTDIR + "/STEP4/samples.json", # reference found
        #low_reads = OUTPUTDIR + "/STEP4/low_reads_samples.json" # no reference found
        json_out = JSON_FILE
    params:
        dir = OUTPUTDIR + '/STEP4/minimap2_final_assembly'
    shell:
        """
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl create_json_annotated_reference \
            -i {params.dir} \
            --min_reads 1 \
            --file_extension gb \
            --jsonfile {output.json_out} \
            --step 4
        """


rule minimap2_final_assembly:
    input: 
        p1 = CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{prefix}_{suffix}/{prefix}_R1_{suffix}.fastq.gz",
        p2 = CWD+"/"+OUTPUTDIR+"/STEP1/remove_contaminants_with_bbduk/{prefix}_{suffix}/{prefix}_R2_{suffix}.fastq.gz",
        ref = CWD+"/"+OUTPUTDIR+"/STEP3/annot/{prefix}_{suffix}.fasta"
    output: 
        sam_file_draft=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}/{prefix}_{suffix}_draft.sam",
        bam_file_draft=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}/{prefix}_{suffix}_draft.bam",
        sam_file=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}/{prefix}_{suffix}.sam",
        bam_file=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}/{prefix}_{suffix}.bam",
        consensus_fasta=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}/{prefix}_{suffix}.fasta"
    params:
        outputdir=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}",
        prefix="{prefix}_{suffix}"
    threads: 1
    log: 
        minimap2_draft=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}_minimap2_draft.log",
        minimap2=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}_minimap2.log",
        sam2consensus=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{prefix}_{suffix}_sam2consensus.log",
    shell:
        """
        minimap2 -ax sr {input.ref} {input.p1} {input.p2} -o {output.sam_file_draft} > {log.minimap2_draft} 2>&1
        samtools view -h -F4 -Sb {output.sam_file_draft} | samtools sort - > {output.bam_file_draft} 2>> {log.minimap2_draft}
        sleep 1;
        samtools index {output.bam_file_draft} >> {log.minimap2_draft} 2>&1
        PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl parse_bam \
            --bam_file {output.bam_file_draft} \
            --output_consensus_file {output.consensus_fasta} 2> {log.sam2consensus}
        minimap2 -ax sr {output.consensus_fasta} {input.p1} {input.p2} -o {output.sam_file}.aux > {log.minimap2} 2>&1
        samtools view -h -F 4  {output.sam_file}.aux -o {output.sam_file} >> {log.minimap2} 2>&1
        rm {output.sam_file}.*
        samtools sort {output.sam_file} -o {output.bam_file} >> {log.minimap2} 2>&1
        """


#            --input_model {input.hmmer_model} \
rule annotate_minimap2_assembly:
    input:
        #unpadded_fasta=CWD+"/"+OUTPUTDIR+"/STEP4/mira_final_assembly/{sample}/{sample}_assembly/{sample}_d_results/{sample}_out_AllStrains.unpadded.fasta",
        unpadded_fasta=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{sample}/{sample}.fasta",
    output:
        #unpadded_gb=CWD+"/"+OUTPUTDIR+"/STEP4/mira_final_assembly/{sample}/{sample}_assembly/{sample}_d_results/{sample}.gb",
        unpadded_gb=CWD+"/"+OUTPUTDIR+"/STEP4/minimap2_final_assembly/{sample}/{sample}.gb",
    params:
        hxb2_genbank_file=config['HXB2_ANNOTATED_REF'],
    log: 
        draft_log=CWD+"/"+OUTPUTDIR+"/STEP4/consensus_fasta/{sample}_draft.log"
    shell:
        """
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl annotate_assembly_blast \
            --reference_file {params.hxb2_genbank_file} \
            --input_file {input.unpadded_fasta} \
            --output_genbank {output} &> {log.draft_log}
        """

