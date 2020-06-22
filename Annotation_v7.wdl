## Copyright Broad Institute, 2020
## This script should take as input a set of raw de novo SNVs and adds annotations as new columns that are used downstream for filtering/QC.
## Outputs the original raw de novo SNVs file w/ additional columns annotations
## VEP annotations: gene symbol, gnomAD max population frequency, biotype, consequence
## Other annotations: GATK RankSum Test values (i.e. PV4; mapQ bias, baseQ bias, readpos bias), 
## strand bias p-value and odds ratio, FDR-based minimum Nalt, repeat regions (mappability, segdup, 
## low-complexity-regions), variant clustering
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
  File cache_dir
  String cache_version

  File convert_script

  File parser_script
  String parser_cols  


  File script_pv4 

  File script_sb 
  Float sb_or 
  Float sb_p 

  File script_fdr 
  Float fdr_min_vaf
  Int fdr_size 
  Float fdr_e_fp 
  Float fdr_seq_err 

  File script_rr_bash 
  File script_rr_parse 
  
  File rr_map 
  File rr_seg 
  File rr_lcr 

  File script_vc 
  Int vc_dist 


  parameter_meta {
    variants: "input raw de novo SNVs in tab-separated text format (output of gvcf_to_denovo)"
    ref_ver: "reference genome version; e.g. GRCh37, GRCh38"
    cache_dir: "gzipped folder containing VEP cache download location"
    cache_version: "cache version being used, important if cache version doesn't match latest release"
    parser_script: "path to script used to parse and append VEP columns to original input file"
    parser_cols: "comma-separated string listing VEP columns to parse; default: SYMBOL,Gene,BIOTYPE,Consequence,Existing_variation,MAX_AF,MAX_AF_POPS"
    script_pv4: "full path to script filter_GATK_RankSum.py"
    script_sb: "full path to script filter_strandbias.R"
    sb_or: "strand bias test Fisher's Exact Test OR threshold value"
    sb_p: "strand bias test Fisher's Exact Test p-value threshold"
    script_fdr: "full path to script filter_fdrmin.R"
    fdr_min_vaf: "floor value for variant allele faction in FDR-based minimum Nalt threshold calculation"
    fdr_size: "region size for FDR-based min Nalt threshold calculation"
    fdr_e_fp: "threshold for expected number of false positives across entire region size"
    fdr_seq_err: "sequencing error rate for FDR-based minimum Nalt threshold calculation"
    script_rr_bash: "full path to bash script run_bedtools.sh"
    script_rr_parse: "full path to python script parse_bedtools_isec.py"
    rr_map: "full path to mappability BEDfile"
    rr_seg: "full path to segmental duplication BEDfile"
    rr_lcr: "full path to LCR BEDfile"
    script_vc: "full path to script flag_vclust.py"
    vc_dist: "distance in bp to define variant clusters"
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

  #run PV4 filter
  call filter_PV4 {
    input:
    infile = add_vep_cols.out,
    script = script_pv4
  }

  #run SB (strand bias) filter
  call filter_SB {
    input:
    infile = filter_PV4.outfile,
    script = script_sb,
    cutoff_or = sb_or,
    cutoff_p = sb_p
  }

  #run FDR (FDR-based min altdp) filter
  call filter_FDR {
    input:
    infile = filter_SB.outfile,
    script = script_fdr,
    min_vaf = fdr_min_vaf,
    size = fdr_size,
    e_fp = fdr_e_fp,
    seq_err = fdr_seq_err
  }

  #run RR (repeat region) filter
  call filter_RR {
    input:
    infile = filter_FDR.outfile,
    script_bash = script_rr_bash,
    script_parse = script_rr_parse,
    map = rr_map,
    seg = rr_seg,
    lcr = rr_lcr
  }

  #run VC (variant cluster) filter
  call filter_VC {
    input:
    infile = filter_RR.outfile,
    script = script_vc,
    dist = vc_dist
  }

  #Outputs original variants file with updated annotation columns 
  output {
    
    File annotated_variants = filter_VC.outfile

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
  File cache_dir # path to location of cache files
  String cache_version
  File vcf

  String cache_dirname = basename(cache_dir, '.tar.gz')
  String outprefix = basename(vcf, '.vcf')  
  String outfname = "VEP_raw.${outprefix}.vcf"

  Int disk_size = 200 # test 100G for VEP?
  
  command <<<

    echo "## extracting and localizing cache directory"
    time tar -xf ${cache_dir} -C ~/
    CACHE_PATH=$(cd ~/${cache_dirname}; pwd)
    echo "## success; see cache directory -- "$CACHE_PATH



    ls -thlR ~/ > home_file_listing.txt

    /opt/vep/src/ensembl-vep/vep \
      --cache \
      --dir_cache "$CACHE_PATH" \
      --cache_version ${cache_version} \
      --offline \
      --format vcf \
      --vcf \
      --force_overwrite \
      --assembly ${ref} \
      --input_file ${vcf} \
      --output_file ${outfname} \
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


    rm -rf ${cache_dir}
    rm -rf $CACHE_PATH

  >>>

  runtime {
    docker: "ensemblorg/ensembl-vep:latest"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: disk_size
  }

  output {
    File flist = "home_file_listing.txt"

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

#Adds PV4 columns to variants file for downstream filtering
task filter_PV4 {
  File infile
  File script
  String outprefix = basename(infile, '.txt')

  command {
    python ${script} ${infile} "${outprefix}.PV4.txt" 
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File outfile = "${outprefix}.PV4.txt"
  }
}

#Adds strand_bias columns to variants file for downstream filtering
task filter_SB {
  File infile
  String outprefix = basename(infile, '.txt')
  File script
  Float cutoff_or
  Float cutoff_p

  command {
    Rscript ${script} ${infile} ${outprefix} ${cutoff_or} ${cutoff_p}
  }
  
  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }
  
  output {
    File outfile = "${outprefix}.SB.txt"
  }
}

#Adds FDR-based minimum Nalt filter-related columns to variants file for downstream filtering
task filter_FDR {
  File infile
  String outprefix = basename(infile, '.txt')
  File script
  Float min_vaf
  Int size 
  Float e_fp
  Float seq_err 

  command {
    Rscript ${script} ${infile} ${outprefix} ${min_vaf} ${size} ${e_fp} ${seq_err}
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File outfile = "${outprefix}.FDR.txt"
  }
}

#Adds Repeat Region (Mappability, SegDup, LCR) filter-related columns to variants file for downstream filtering
task filter_RR {
  File infile
  File script_bash
  File script_parse

  File map
  File seg
  File lcr
  String outprefix = basename(infile, '.txt')

  command {
    git clone https://github.com/arq5x/bedtools2.git

    BEDPATH=$(cd bedtools2/src; pwd)

    sh ${script_bash} -v ${infile} -b "$BEDPATH" -m ${map} -s ${seg} -l ${lcr} -p ${script_parse}
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File tmpbed = "tmp.bed" # temporary BEDfile created from variants file
    File isec_file = "bed.isec.out.txt" # BEDtools intersect output
    File outfile = "${outprefix}.RR.txt"
  }
}

#Adds Variant Cluster filter-related columns to variants file for downstream filtering
task filter_VC {
  File infile
  File script 
  Int dist # distance (in bp) used to define a "cluster"
  String outprefix = basename(infile, '.txt')

  command {
    python ${script} ${infile} ${dist} "${outprefix}.VC.txt"
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File outfile = "${outprefix}.VC.txt"
  }
}

