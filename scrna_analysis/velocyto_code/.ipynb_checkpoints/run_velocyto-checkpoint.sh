#!/bin/bash

set -euo pipefail

#eval "$(conda shell.bash hook)"
#conda activate velocyto

start=`date +%s`

# Set the directory path where your .zip files are located
aligned_path="/ictstr01/home/icb/daniel.garger/Atherosclerosis/aligned_data"
output_dir="/ictstr01/home/icb/daniel.garger/Atherosclerosis/velocyto/output"
raw_data_path="/lustre/groups/cbm01/datasets/daniel.garger/plaque_raw_data"


## Set batches and reference genomes to loop over
batches=("Batch6")
ref_genomes=("GRCh38-p13-Gencode_v33" "GRCh38-p14-Gencode_v44")
#ref_genomes=("GRCh38-p14-Gencode_v44")
sample_name_list=("CAR17_D") #  CAR16_D CAR17_D
	

echo "$start"

for batch_name  in "${batches[@]}"; do
   
   
    for ref_genome in "${ref_genomes[@]}"; do
       
       # Change to the specified directory
       batch_path="$aligned_path/$ref_genome/$batch_name"
       cd "$batch_path"

       for sample_directory in "$batch_path"/*; do

            sample_name=$(basename "$sample_directory")
            
            ## Drop the "Sample_" prefix from the Sample name, and check if the sample is in the list of Samples to run cellranger on
            #sample_name_short="${sample_name#Sample_}"
	    sample_name_short="${sample_name%_count}"
            echo $sample_name_short
	    
            # Check if the sample name is in the list of valid names
            if [[ " ${sample_name_list[@]} " =~ "${sample_name_short}" ]]; then
            
                echo $sample_directory

                # cd into the outs directory of sample
                cd "$sample_directory/outs"

                echo $sample_name

                # Check for .BAM files, if they are there, start velocyto
                if ls *.bam 1> /dev/null 2>&1; then

                  echo "BAM file exists, starting velocyto for sample $sample_name with genome ${ref_genome}"

                  ## Create path to extract zip files of sample into   
                  output_path="$output_dir/$ref_genome/$batch_name/$sample_name"
                  echo "Creating output directory: $output_path"
                  mkdir -p "$output_path"
                  echo "Output directory created: $output_path"



                    ## IF NO .BAM FILES PRESENT IN THE FOLDER, SKIP TO THE NEXT SAMPLE
                    else
                    echo "There are no .bam files in the directory of sample $sample_name."
                    continue  # Skip to the next iteration of the loop
                fi        

                ## UNZIP BARCODES.TSV IF NECESSARY
                filt_barcodes_path="$sample_directory/outs/filtered_feature_bc_matrix/barcodes.tsv"

                if [ ! -e $filt_barcodes_path ]; then
                    echo "Extracting barcodes.tsv"
                    filt_barcodes_gz_path="$sample_directory/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
                    gzip -dk "$filt_barcodes_gz_path"

                else
                    echo "Barcodes.tsv already exists"
                fi



                ## CREATE RUN_ID
                run_id=${sample_name}

                ## Create repeat_mask_path
                repeat_mask_path="/ictstr01/home/icb/daniel.garger/Atherosclerosis/velocyto/hg38_repeat_mask.gtf"

                ## CREATE .BAM PATH
                bam_path="$sample_directory/outs/possorted_genome_bam.bam"

                 ## CREATE ANNOTATION .GTF PATH
                 if [[ $ref_genome == *v33 ]]; then
                    annot_path="/ictstr01/home/icb/daniel.garger/cellranger-7.1.0/reference_sources/gencode.v33.primary_assembly.annotation.gtf"

                elif [[ $ref_genome == *v44 ]]; then
                    annot_path="/ictstr01/home/icb/daniel.garger/cellranger-7.1.0/reference_sources/gencode.v44.primary_assembly.annotation.gtf"

                fi

                ## RUN VELOCYTO
                velocyto run -b $filt_barcodes_path -o $output_path -m $repeat_mask_path $bam_path $annot_path -vvv

                rm $filt_barcodes_path


                echo "Velocyto run finished for for sample $run_id with genome ${ref_genome}" 

            fi
            
        done

    done

done

end=`date +%s`
runtime=$((end-start))

echo "Runtime: $runtime"
