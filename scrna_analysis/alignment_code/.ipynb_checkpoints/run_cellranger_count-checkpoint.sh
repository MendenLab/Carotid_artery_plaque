#!/bin/bash

set -eo pipefail


start=`date +%s`

# Set the directory path where your .zip files are located
output_path="/home/icb/daniel.garger/Atherosclerosis/aligned_data"
raw_data_path="/lustre/groups/cbm01/datasets/daniel.garger/plaque_raw_data"


## Set batches and reference genomes to loop over
batches=("Batch2" "Batch3" "Batch4") #"Batch1"
batches=( "Batch6")
ref_genomes=("GRCh38-p13-Gencode_v33" "GRCh38-p14-Gencode_v44")
#ref_genomes=("GRCh38-p14-Gencode_v44")
sample_name_list=("CAR17_D") #  CAR15_H CAR16_D CAR17_D already done: CAR11_D



echo "$start"

for batch_name  in "${batches[@]}"; do

   # Change to the specified directory
   batch_path="$raw_data_path/$batch_name"
   cd "$batch_path"

   for sample_directory in "$batch_path"/*; do

        sample_name=$(basename "$sample_directory")
        #echo "Sample name is: $sample_name"
        
        ## Drop the "Sample_" prefix from the Sample name, and check if the sample is in the list of Samples to run cellranger on
        sample_name_short="${sample_name#Sample_}"
        
         
       # Check if the sample name is in the list of valid names
        if [[ " ${sample_name_list[@]} " =~ "${sample_name_short}" ]]; then
        
        
            # Loop over each .zip file and unzip them
            cd "$sample_directory"

            # Check for .zip files -> if there are loop over them and extract them
            if ls *.zip 1> /dev/null 2>&1; then

              echo "Unzipping: processing sample subdirectory: $sample_name"

              ## Create path to extract zip files of sample into   
              extract_path="$output_path/temp/$sample_name"
              echo "Creating directory: $extract_path"
              mkdir -p "$extract_path"
              echo "Directory created: $extract_path"

              for zip_file in *.zip; do
                  ## Extract ZIP file
                  echo $zip_file 	     
                  cd $extract_path
                  jar xvf "$sample_directory/$zip_file"
                  cd $sample_directory

                  ## If no errors (jar command returns 0), print it to terminal
                  if [ $? -eq 0 ]; then
                    echo "$(basename "$zip_file") successfully extracted"

                  else
                     echo "Error extracting: $(basename "$zip_file")"
                  fi

               done

            ## IF NO .ZIP FILES ON THE FOLDER ONLY .FASTQ.GZ FILES, COPY THEM THE THE TEMP/SAMPLE_NAME DIRECTORY, WHICH WILL BE UASED AS INPUT FILES FOR CELLRANGER
            elif ls *.fastq.gz 1> /dev/null 2>&1; then

               echo "No.zip files found,only fastq.gz files.Copying fast.gz files to temp: $sample_name"
               ## Create path to copy fastq.gz files of sample into
               extract_path="$output_path/temp/$sample_name"
               echo "Directory created: $extract_path"
               mkdir -p "$extract_path"

               for fastq_file in *.fastq.gz; do
                scp -r "$sample_directory/$fastq_file" "$extract_path"
               done

            ## IF NO .ZIP OR FASTQ.GZ FILES PRESENT IN THE FOLDER, SKIP TO THE NEXT SAMPLE
            else
                echo "No .zip or .fastq.gz files found. Jumping to next sample"
                continue  # Skip to the next iteration of the loop
            fi        


            ## REMOVE _I1_001.FAST.GZ FILES => THEY ARE USED FOR DEMULTIPLEXING (ALREADY DONE) AND GIVE AN ERROR WITH CELLRANGER COUNT
            find "$extract_path" -type f -name '*_I1_00*' -exec rm {} +


            ## RUN CELLRANGER FOR EACH REFERNCE GENOME
            for ref_genome in "${ref_genomes[@]}"; do

                ## DELETE "Sample_" FROM SAMPLE NAME -> THIS WAY CELLRANGER RECOGNIZES SAMPLE NAME AUTOMATICALLY
                sample_name_short="${sample_name//Sample_/}"

                ## CREATE RUN_ID
                run_id=${sample_name_short}_count
                echo "Running cellranger for sample $sample_name_short with genome ${ref_genome}"

                ## Create REF_GENOME_PATH
                ref_genome_path="/home/icb/daniel.garger/Atherosclerosis/star_ref_genome/$ref_genome"

                ## CREATE FOLDER TO COLLECT CELLRANGER OUTPUT FOR SAMPLE => cd into it
                sample_cellr_output_path="$output_path/$ref_genome/$batch_name"
                mkdir -p "$sample_cellr_output_path"
                cd "$sample_cellr_output_path"

                # echo "run_id $run_id"
                # echo "extract_path $extract_path"
                # echo "sample_name_short $sample_name_short"
                # echo "sample_cellr_output_path $sample_cellr_output_path"

                ## RUN CELLRANGE COUNT
                cellranger count --id=$run_id \
                    --fastqs=$extract_path \
                    --sample=$sample_name_short \
                    --include-introns=false \
                    --transcriptome=$ref_genome_path #\
                #	    --no-bam

                echo "Cellranger run finished for for sample $run_id with genome ${ref_genome}" 

            done	

            ## DELETE EXTRACTED FASTQ FILES of SAMPLE  
            echo "Removing Fastq files of sample $sample_name from temp directory"
            rm -r $extract_path
        
    
    
        fi
    done

done

## Remove temp directory
echo "Removing $output_path/temp directory"
rm -rf "$output_path/temp"

end=`date +%s`
runtime=$((end-start))

echo "Runtime: $runtime"
