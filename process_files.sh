#!/bin/bash

cur_path=$(pwd)'/data/experiment'
pwd=$(pwd)
echo $cur_path
mkdir -p ./data/experiment

# Change next two lines
cd "C:\Users\Andrej\Desktop\Magisterij\2_letnik\NPMP\NPMP-inference\reviewAndAssessment"
conda activate dl

myArray=('results')
for option in ${myArray[@]}; do
    cd $option
    mkdir -p $cur_path'/'$option
    for DIR in *; do 
        mkdir -p $cur_path'/'$option'/'$DIR; 
        cd $DIR
        for DIR_n in *; do
            if [ -d "$DIR_n" ]; then 
                cd $DIR_n
                #data/Ecoli/16
                mkdir -p $cur_path'/'$option'/'$DIR'/'$DIR_n
                for FILE in *structure.tsv; do
                    if [ -f "$FILE" ]; then
                        echo $FILE;
                        python $pwd/to_pajek.py --path_in $(pwd)'/' --filename $FILE --path_out $cur_path'/'$option'/'$DIR'/'$DIR_n
                    fi
                done
                for METHOD_DIR in *; do
                    if [ -d "$METHOD_DIR" ]; then
                        mkdir -p $cur_path'/'$option'/'$DIR'/'$DIR_n'/'$METHOD_DIR
                        cd $METHOD_DIR
                        for FILE in *structure.tsv; do
                            if [ -f "$FILE" ]; then
                                echo $FILE;
                                python $pwd/to_pajek.py --path_in $(pwd)'/' --filename $FILE --path_out $cur_path'/'$option'/'$DIR'/'$DIR_n'/'$METHOD_DIR
                            fi
                        done
                        cd ..
                    fi
                done
                cd ..
            fi
        done
        cd ..
    done
    cd ..
done