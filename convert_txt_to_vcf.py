## Purpose: take raw variants file and format as ANNOVAR input file (.avinput)
## 			chr, start, end, ref, alt, otherinfo
'''
Usage: prep_avinput.py [options]

Options:
  -h, --help            show this help message and exit
  -i input_file, --input=input_file
                        input tab-separated variants file
  -c config_file, --config=config_file
                        configuration file
  -o output_file, --output=output_file
'''
import sys
import yaml
from optparse import OptionParser
import os

## handle arguments
parser = OptionParser()
parser.add_option('-i', '--input', dest='input_file',help='input tab-separated variants file')
parser.add_option('-c', '--config', dest='config_file',help='configuration yml file')
parser.add_option('-o', '--output', dest='output_file',help='output filename')
(options, args) = parser.parse_args()

## check all arguments present
if (options.input_file == None or options.config_file == None or options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	sys.exit()


input_file = options.input_file
config_file = options.config_file
output_file = options.output_file


with open(config_file, 'r') as ymlf:
	cfg = yaml.load(ymlf, Loader=yaml.FullLoader)

'''
for section in cfg:
    print(section)
    for key, value in cfg[section].items():
        print '     %s: %s' % (key, value)
'''

vcfout = open(output_file, 'w')

with open(input_file,'r') as f1:
	for line in f1:
		tmp = line.strip().split('\t')
		if not tmp[0].startswith('##'):
			## first col of header is either (1) VCF-style starting with #CHROM or (2) some id field
			if tmp[0].startswith('#') or tmp[0].lower() in cfg['head']['id']:
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
					print '\n' + '## ERROR: problem with header column mapping' 
					print '##        check that all relevant columns are mapped in config.yml' + '\n'
					sys.exit()
				## write out VCF header
				vcfout.write("##fileformat=VCFv4.1"+'\n')
				head = '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
				#head = '\t'.join(['#CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
				vcfout.write(head+'\n')
			else:
				id = tmp[idxmap['id']]
				chr = tmp[idxmap['chr']]
				pos = tmp[idxmap['pos']]
				ref = tmp[idxmap['ref']]
				alt = tmp[idxmap['alt']]

				qual, filt, info, format = '.', '.','.','GT:AD:DP:GQ:PL'

				outstring = '\t'.join([chr, pos,id,  ref, alt, qual, filt, info, format])
				#outstring = '\t'.join([chr, pos, ref, alt, qual, filt, info, format])
				vcfout.write(outstring+'\n')

vcfout.close()


