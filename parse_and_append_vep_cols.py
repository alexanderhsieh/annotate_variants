## Purpose: parse VEP CSQ output and append designated columns to original denovos file
'''
Usage: parse_and_append_vep_cols.py [options] > output file

Options:
  -h, --help            show this help message and exit
  -i input_file, --input=input_file
                        input tab-separated variants file
  -v vep_file, --vep=vep_file
                        VEP-annotated VCF
  -c cols, --col=cols
                    	comma-separated list of columns to parse and append
  -o output_file, --output=output_file
                    	output tab-separated variants file with additional columns
'''
import sys
import yaml
from optparse import OptionParser
import subprocess
import os


## handle arguments
parser = OptionParser()
parser.add_option('-i', '--input', dest='input_file',help='input tab-separated variants file')
parser.add_option('-v', '--vep', dest='vep_file',help='VEP-annotated VCF')
parser.add_option('-c', '--col', dest='cols', help='comma-separated list of columns to parse and append')
parser.add_option('-y', '--yml', dest='yml_file',help='config YML file')
parser.add_option('-o', '--output', dest='output_file',help='output tab-separated variants file')
(options, args) = parser.parse_args()

## check all arguments present
if (options.input_file == None or options.vep_file == None or options.cols == None or options.yml_file == None or options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()


input_file = options.input_file
vep_file = options.vep_file
cols = options.cols
yml = options.yml_file
output_file = options.output_file



## read config file for parsing header: id, chr, pos, ref, alt, etc
with open(yml, 'r') as ymlf:
	cfg = yaml.load(ymlf, Loader=yaml.FullLoader)

## parse VEP file and store as dict

vepd = {} # { id:chr:pos : csq } ; csq= { colname : colval }
with open(vep_file, 'r') as vf:
	for line in vf:
		## parse CSQ format
		if line.startswith('##INFO'):
			csqfmt = line.split('Format: ')[1].split('">')[0].split('|')
		## iterate over non-INFO lines
		elif not line.startswith('##'):
			tmp = line.strip().split('\t')
			if tmp[0] == '#CHROM':
				idx = {col:index for index, col in enumerate(tmp)}
			else:
				id = tmp[idx['ID']]
				chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]

				csqval = tmp[idx['INFO']].split('CSQ=')[1].split('|')

				csq = dict(zip(csqfmt, csqval))

				#print '\n'.join(sorted(csq.keys()))
				#print '\t'.join([csq[c] for c in cols])

				key = ':'.join([id, chr, pos])

				vepd[key] = csq


## read input file and append designated columns

## list of output columns
outcols = cols.split(',')

## read header line to get first column value of header
head = subprocess.Popen("head -n 1 %s"%(input_file), shell=True, stdout=subprocess.PIPE).stdout.read()
head_token = head.split('\t')[0]

## added 03/26/2020; handles alternative input file format (VCF instead of txt)
if input_file.endswith('.vcf'):
	head_token = '#CHROM'


## added 04/16/2020: write to output_file instead of print to stdout()
outf = open(output_file, 'w')

with open(input_file, 'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if line.startswith('##'): ## handle VCF header lines
			#print '\t'.join(tmp)
			outf.write('\t'.join(tmp) + '\n')
		else:
			if tmp[0] == head_token or tmp[0].lower() in cfg['head']['id']:
				idx = {col.lower():index for index, col in enumerate(tmp)} # original header dict
				idxmap = {} # re-mapped header dict
				## re-map key fields (id, chr, pos, ref, alt)
				## print error and exit if there is a problem with column header parsing
				try:
					idxmap['id'] = idx[list(set(idx.keys()).intersection(cfg['head']['id']))[0]]
					idxmap['chr'] = idx[list(set(idx.keys()).intersection(cfg['head']['chr']))[0]]
					idxmap['pos'] = idx[list(set(idx.keys()).intersection(cfg['head']['pos']))[0]]
					idxmap['ref'] = idx[list(set(idx.keys()).intersection(cfg['head']['ref']))[0]]
					idxmap['alt'] = idx[list(set(idx.keys()).intersection(cfg['head']['alt']))[0]]
				except:
					print('\n' + '## ERROR: problem with header column mapping')
					print('##        check that all relevant columns are mapped in config.yml' + '\n')
					sys.exit()

				## print header of output file
				#print '\t'.join(tmp) + '\t' + '\t'.join(outcols)
				outf.write('\t'.join(tmp) + '\t' + '\t'.join(outcols) + '\n')
			
			else:
				
				#id = tmp[idx['ID']]
				#chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]

				#id = tmp[idx[head_token]]
				#chr, pos, ref, alt = tmp[idx['CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]
				#chr, pos, ref, alt = tmp[idx['chrom']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]
				#chr, pos, ref, alt = tmp[idx['chr']], tmp[idx['pos']], tmp[idx['ref']], tmp[idx['alt']]
				id = tmp[idxmap['id']]
				chr = tmp[idxmap['chr']]
				pos = tmp[idxmap['pos']]
				ref = tmp[idxmap['ref']]
				alt = tmp[idxmap['alt']]

				key = ':'.join([id, chr, pos])

				## fetch CSQ dictionary for this variant
				tmpd = vepd[key]

				## get user-specificed column values
				to_append = [tmpd[c] for c in outcols]

				## replace missing annotations with '.'
				to_append_fixed = ['.' if len(w) == 0 else w for w in to_append]  

				#print '\t'.join(tmp) + '\t' + '\t'.join(to_append_fixed)
				outf.write('\t'.join(tmp) + '\t' + '\t'.join(to_append_fixed) + '\n')


outf.close()
