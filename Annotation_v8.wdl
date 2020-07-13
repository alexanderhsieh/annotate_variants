version 1.0
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
	input {
		File variants
		String ref_ver
		File cache_dir
		#String cache_version

		File convert_script

		File parser_script
		#String parser_cols


		File script_pv4 

		File script_sb 
		#Float sb_or 
		#Float sb_p 

		File script_fdr 
		#Float fdr_min_vaf
		#Int fdr_size 
		#Float fdr_e_fp 
		#Float fdr_seq_err 

		File script_rr_bash 
		File script_rr_parse 

		File rr_map 
		File rr_seg 
		File rr_lcr 

		File script_vc 
		#Int vc_dist 
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
			vcf_idx = txt_to_vcf.idx,
			cache_dir = cache_dir
			#cache_version = cache_version
	}

	# Step 2: Parse and append VEP columns to original vcf file
	call add_vep_cols {
		input:
			original_variants = variants,
			vep_vcf = run_vep.vep_out,
			script = parser_script
			#cols = parser_cols
	}

	#run PV4 filter
	call flag_PV4 {
		input:
			infile = add_vep_cols.out,
			script = script_pv4
	}

	#run SB (strand bias) filter
	call flag_SB {
		input:
			infile = flag_PV4.out,
			script = script_sb
			#cutoff_or = sb_or,
			#cutoff_p = sb_p
	}

	#run FDR (FDR-based min altdp) filter
	call flag_FDR {
		input:
			infile = flag_SB.out,
			script = script_fdr
			#min_vaf = fdr_min_vaf,
			#size = fdr_size,
			#e_fp = fdr_e_fp,
			#seq_err = fdr_seq_err
	}

	#run RR (repeat region) filter
	call flag_RR {
		input:
			infile = flag_FDR.out,
			script_bash = script_rr_bash,
			script_parse = script_rr_parse,
			map = rr_map,
			seg = rr_seg,
			lcr = rr_lcr
	}

	#run VC (variant cluster) filter
	call flag_VC {
		input:
			infile = flag_RR.out,
			script = script_vc
			#dist = vc_dist
	}

	#Outputs original variants file with updated annotation columns 
	output {

		File annotated_variants = flag_VC.out

	}

}


###########################################################################
#Task Definitions
###########################################################################
# converts tab-separated variants file into vcf format for VEP
task txt_to_vcf {
	input {
		File variants
		File script
	}

	String outprefix = basename(variants, '.raw.txt')  
	String outfname = "~{outprefix}.vcf"
	String outfname_gz = "~{outprefix}.vcf.gz"

	command {
		python ~{script} -i ~{variants} -o ~{outfname}

		bgzip -c ~{outfname} > ~{outfname_gz}

		tabix -p vcf ~{outfname_gz}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outfname_gz}"
		File idx = "~{outfname_gz}.tbi"
	}

}


#Runs VEP in the same task to preserve environment
## NOTE: REQUIRES ~8GB memory (docker on mac allocates 2GB by default) - need to increase memory limit if running locally
## see: https://gatkforums.broadinstitute.org/wdl/discussion/11522/the-job-was-aborted-from-outside-cromwell-sporadic-failures-when-running-cromwell-locally
task run_vep {
	input {
		String ref # GRCh37 or GRCh38
		File cache_dir # path to location of cache files
		String cache_version = "100"
		File vcf
		File vcf_idx
		Int disk_size = 125 # test 100G for VEP?
	}


	String cache_dirname = basename(cache_dir, '.tar.gz')
	#String outprefix = basename(vcf, '.vcf')  
	String outprefix = basename(vcf, '.vcf.gz')  
	String outfname = "VEP_raw.~{outprefix}.vcf"

	command {

		echo "## extracting and localizing cache directory"
		time tar -xf ~{cache_dir} -C ~/
		CACHE_PATH=$(cd ~/~{cache_dirname}; pwd)
		echo "## success; see cache directory -- "$CACHE_PATH

		/opt/vep/src/ensembl-vep/vep \
		--cache \
		--dir_cache "$CACHE_PATH" \
		--cache_version ~{cache_version} \
		--offline \
		--format vcf \
		--vcf \
		--assembly ~{ref} \
		--input_file ~{vcf} \
		--output_file ~{outfname} \
		--no_stats \
		--pick \
		--symbol \
		--canonical \
		--biotype \
		--max_af 


		rm -rf ~{cache_dir}
		rm -rf $CACHE_PATH

	}

	runtime {
		docker: "ensemblorg/ensembl-vep:latest"
		disks: "local-disk " + disk_size + " HDD"
		bootDiskSizeGb: disk_size
		memory: "16G"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File vep_out = "~{outfname}"
	}

}

#Parses VEP output and appends user-defined columns to the original VCF file
task add_vep_cols {
	input {
		File original_variants
		File vep_vcf
		File script
		String cols = "SYMBOL,Gene,BIOTYPE,Consequence,Existing_variation,MAX_AF,MAX_AF_POPS"
	}

	String outprefix = basename(original_variants, '.raw.txt')


	command {
		python ~{script} -i ~{original_variants} -v ~{vep_vcf} -c ~{cols} -o "${outprefix}.VEP.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.VEP.txt"
	}

}

#Adds PV4 columns to variants file for downstream filtering
task flag_PV4 {
	input {
		File infile
		File script
	}

	String outprefix = basename(infile, '.txt')

	command {
		python ~{script} ~{infile} "~{outprefix}.PV4.txt" 
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.PV4.txt"
	}
}

#Adds strand_bias columns to variants file for downstream filtering
task flag_SB {
	input {
		File infile
		File script
		Float cutoff_or = 3.0
		Float cutoff_p = 0.001
	}
	String outprefix = basename(infile, '.txt')

	command {
		Rscript ~{script} ~{infile} ~{outprefix} ~{cutoff_or} ~{cutoff_p}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.SB.txt"
	}
}

#Adds FDR-based minimum Nalt filter-related columns to variants file for downstream filtering
task flag_FDR {
	input {
		File infile
		File script
		Float min_vaf = 0.01
		Int size = 30000000
		Float e_fp = 0.01
		Float seq_err = 0.005
	}
	String outprefix = basename(infile, '.txt')

	command {
		Rscript ~{script} ~{infile} ~{outprefix} ~{min_vaf} ~{size} ~{e_fp} ~{seq_err}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.FDR.txt"
	}
}

#Adds Repeat Region (Mappability, SegDup, LCR) filter-related columns to variants file for downstream filtering
task flag_RR {
	input {
		File infile
		File script_bash
		File script_parse

		File map
		File seg
		File lcr
	}

	String outprefix = basename(infile, '.txt')

	command {
		git clone https://github.com/arq5x/bedtools2.git

		BEDPATH=$(cd bedtools2/src; pwd)

		sh ~{script_bash} -v ~{infile} -b "$BEDPATH" -m ~{map} -s ~{seg} -l ~{lcr} -p ~{script_parse}
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File tmpbed = "tmp.bed" # temporary BEDfile created from variants file
		File isec_file = "bed.isec.out.txt" # BEDtools intersect output
		File out = "~{outprefix}.RR.txt"
	}
}


#Adds Variant Cluster filter-related columns to variants file for downstream filtering
task flag_VC {
	input {
		File infile
		File script 
		Int dist = 10# distance (in bp) used to define a "cluster"
	}

	String outprefix = basename(infile, '.txt')

	command {
		python ~{script} ~{infile} ~{dist} "~{outprefix}.VC.txt"
	}

	runtime {
		docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
		preemptible: 3
		maxRetries: 3
	}

	output {
		File out = "~{outprefix}.VC.txt"
	}
}

