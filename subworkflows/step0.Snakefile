import glob, os

configfile: 'config.yml'
OUTPUTDIR = config['OUTPUTDIR']


rule all:
    input:
        [OUTPUTDIR + "/check_pairs.out",
         OUTPUTDIR + "/check_duplicated.out"]


rule check_pairs:
    input:
        data_dir = config['DATADIR']
    output:
        OUTPUTDIR + "/check_pairs.out"
    log:
        OUTPUTDIR + "/check_pairs.log"
    shell:
        """

            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl inspect_pairs --data_dir {input} --output {output}

        """


rule check_duplicated_files:
    input:
        data_dir = config['DATADIR'],
        out_first = OUTPUTDIR + "/check_pairs.out"
    output:
        outfile = OUTPUTDIR + "/check_duplicated.out"
    params:
        md5_file = OUTPUTDIR + "/md5_parallel.txt"
    log:
        OUTPUTDIR + "/check_duplicated.log"
    shell:
        """
            find {input.data_dir}/*.gz | parallel  -j 60 md5sum {{}} > {params.md5_file}
            PERL5LIB=$PERL5LIB:$PWD/scripts/lib scripts/HIVAtoolkit.pl check_duplicated_files \
                --md5sum_file {params.md5_file} \
                --duplicated_files_list {log}
            if [[ ! -e {log}  ]];then
                echo 'NO DUPLICATED FILES!!'
                touch {output.outfile}
            else
                echo 'DUPLICATED FILES !!'
                exit 1;
            fi
        """
