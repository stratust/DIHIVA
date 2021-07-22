#!/bin/bash

INI_DIR=$PWD

DIR_NAME=$(basename $PWD);

echo "Currently in: $PWD"

check_execution()
{
    STEP=$1
    shift

    cd $INI_DIR
    mkdir $PWD/results_dihiva

    if [ $STEP -eq 7 ]; then
        for i in {5..7}; do
            if [ $i -eq 6 ]; then
                continue

            else

                if [ $i -eq 5 ];then
                    cp -r $INI_DIR/results/STEP$i $INI_DIR/results_dihiva/STEP$i-Assembly_Classification
                    rm $INI_DIR/results_dihiva/STEP$i-Assembly_Classification/JSON_STEP$i.json
                    rm -rf $INI_DIR/results_dihiva/STEP$i
                fi

                if [ $i -eq 7 ];then
                    cp -r $INI_DIR/results/STEP$i $INI_DIR/results_dihiva/STEP$i-Sequence_Classification
                    rm  $INI_DIR/results_dihiva/STEP$i-Sequence_Classification/MergedJSON.json
                    rm -rf $INI_DIR/results_dihiva/STEP$i
                fi

            fi

        done


    elif [ $STEP -eq 0 ]; then
        cp -r $INI_DIR/results/summary.html $INI_DIR/results_dihiva

    elif [ $STEP -eq 1 ]; then
        cp -r $INI_DIR/results/STEP1/summary.html $INI_DIR/results_dihiva

    elif [ $STEP -eq 2 ]; then
        cp -r $INI_DIR/results/STEP2/summary.html $INI_DIR/results_dihiva

    elif [ $STEP -eq 5 ]; then
        cp -r $INI_DIR/results/STEP5/summary.html $INI_DIR/results_dihiva
    fi

    zip -r $INI_DIR/$analysisname.zip results_dihiva/

    echo "DIHIVAFINISHED"
    touch DIHIVAFINISHED.txt

}


cores=$1
analysisname=$2
mem=$3

if snakemake -s local.Snakefile -j $cores --greediness 0.8 --resources mem_mb=$mem ;then
    check_execution 7

elif [ -f "results/summary.html" ]; then
    check_execution 0

elif [ -f "results/STEP1/summary.html" ]; then
    check_execution 1

elif [ -f "results/STEP2/summary.html" ]; then
    check_execution 2

elif [ -f "results/STEP5/summary.html" ]; then
    check_execution 5

else
    echo "We have a problem"
fi


