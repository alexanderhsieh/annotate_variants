## Copyright Broad Institute, 2020
## This script should take as input a set of raw de novo SNVs and annotates them using
## Ensembl Variant Effect Predictor (VEP).  This script also parses the annotation files and
## appends user-selected columns to the original raw de novo SNV VCF File.
## Output includes (1) original raw de novo SNVs annotated and (2) raw VEP output 
## NOTE: scatter-gather VEP step requires ~8GB memory (docker on mac allocates 2GB by default; need to adjust memory limit)
## 
## TESTED: 
## Versions of other tools on this image at the time of testing:
##
## LICENSING : This script is released under the WDL source code license (BSD-3) (see LICENSE in https://github.com/broadinstitute/wdl). 
## Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. 
## Please see the docker for detailed licensing information pertaining to the included programs.
##

###########################################################################
#Workflow Definition
###########################################################################
workflow annotate_variants {
  
  File variants
  String ref_ver
  String cache_dir
  String cache_version

  File convert_script

  File parser_script
  String parser_cols  

  String final_prefix = basename(variants, ".raw.txt")


  parameter_meta {
    variants: "input raw de novo SNVs in tab-separated text format (output of gvcf_to_denovo)"
    ref_ver: "reference genome version; e.g. GRCh37, GRCh38"
    cache_dir: "path to VEP cache download location"
    cache_version: "cache version being used, important if cache version doesn't match latest release"
    parser_script: "path to script used to parse and append VEP columns to original input file"
    parser_cols: "comma-separated string listing VEP columns to parse; default: SYMBOL,Gene,BIOTYPE,Consequence,Existing_variation,MAX_AF,MAX_AF_POPS"
  }


  meta {
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }



  # Step 0: convert txt to vcf
  call txt_to_vcf {
    input:
    variants = variants,
    script = convert_script
  }
  
  # Step 1: generate VEP annotations
  call run_vep {
    input:
    ref = ref_ver,
    vcf = txt_to_vcf.out,
    cache_dir = cache_dir,
    cache_version = cache_version
  }



  # Step 2: Parse and append VEP columns to original vcf file
  call add_vep_cols {
    input:
    original_variants = variants,
    vep_vcf = run_vep.vep_out,
    script = parser_script,
    cols = parser_cols
  }



  #Outputs (1) original VCF with updated annotation columns and (2) raw VEP output
  output {
    
    File annotated_variants = add_vep_cols.out

  }

}


###########################################################################
#Task Definitions
###########################################################################
# converts tab-separated variants file into vcf format for VEP
task txt_to_vcf {
  File variants
  String outprefix = basename(variants, '.raw.txt')  
  String outfname = "${outprefix}.vcf"

  File script

  command <<<


    python ${script} -i ${variants} -o ${outfname}

  >>>

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File out = "${outfname}"
  }

}

#Runs VEP in the same task to preserve environment
## NOTE: REQUIRES ~8GB memory (docker on mac allocates 2GB by default) - need to increase memory limit if running locally
## see: https://gatkforums.broadinstitute.org/wdl/discussion/11522/the-job-was-aborted-from-outside-cromwell-sporadic-failures-when-running-cromwell-locally
task run_vep {
  
  String ref # GRCh37 or GRCh38
  String cache_dir # path to location of cache files
  String cache_version
  File vcf

  String outprefix = basename(vcf, '.vcf')  
  String outfname = "VEP_raw.${outprefix}.vcf"
  
  command <<<
    
    vcf_name=`basename ${vcf}`
    input_path=`dirname ${vcf}` 
    output_path=$(pwd)


    # run VEP
    ## note: may need to handle potential disagreements in VEP and cache versions via arguments
    ## used workaround described in https://gatkforums.broadinstitute.org/gatk/discussion/comment/50056#Comment_50056
    
    docker run -v $input_path:/home/var/ -v $output_path:/home/out/ -v ${cache_dir}:/home/vep/.vep ensemblorg/ensembl-vep:latest ./vep \
      --cache \
      --dir_cache "/home/vep/.vep" \
      --cache_version ${cache_version} \
      --offline \
      --format vcf \
      --vcf \
      --force_overwrite \
      --assembly ${ref} \
      --input_file "/home/var/$vcf_name" \
      --output_file "/home/out/${outfname}" \
      --no_stats \
      --pick \
      --gencode_basic \
      --hgvs \
      --symbol \
      --transcript_version \
      --canonical \
      --biotype \
      --numbers \
      --max_af \
      --af_gnomad

  >>>

  output {
    File vep_out = "${outfname}"
  }

}


#Parses VEP output and appends user-defined columns to the original VCF file
task add_vep_cols {
  File original_variants
  File vep_vcf
  File script
  String cols
  
  String outprefix = basename(original_variants, '.raw.txt')


  command <<<

    python ${script} -i ${original_variants} -v ${vep_vcf} -c ${cols} -o "${outprefix}.VEP.txt"
  >>>

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File out = "${outprefix}.VEP.txt"
  }

}


