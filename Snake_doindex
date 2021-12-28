# Download reference genome files and run STAR to generate INDEX
# johaGL, 2021


species = config['species']
urlfa = config['urls']['fasta']
urlgtf = config['urls']['gtf']
urlgff3 = config['urls']['gff3']
refdir = config['refdir']

# extracted file names
ofa = urlfa.split("/")[-1].replace(".gz","") 
ogt = urlgtf.split("/")[-1].replace(".gz","") 
ogff = urlgff3.split("/")[-1].replace(".gz","")

# Rules
rule all:
	input : f'{refdir}{species}/{ofa}',
		f'{refdir}{species}/STARIndex_{species}/README.txt'

rule getfiles:
	output : fa = f'{refdir}{species}/{ofa}',
		gtf = f'{refdir}{species}/{ogt}',
		gff3 = f'{refdir}{species}/{ogff}'
	shell : "wget -O {output.fa}.gz {urlfa}; gunzip -f {output.fa}.gz > {output.fa} ; \
		wget -O {output.gtf}.gz {urlgtf}; gunzip -f {output.gtf}.gz > {output.gtf} ; \
		wget -O {output.gff3}.gz {urlgff3}; gunzip -f {output.gff3}.gz > {output.gff3} "	
	
rule BuildStarIndex:
	input : fa = f'{refdir}{species}/{ofa}',
		gtf = f'{refdir}{species}/{ogt}'
	output : infoindex = f'{refdir}{species}/STARIndex_{species}/README.txt'
	shell : "printf 'Files used to bluild index {species} : \n1. {input.fa} \n2. {input.gtf}\n' > {output.infoindex}; \
		STAR --runThreadN 6 --runMode genomeGenerate \
		  --genomeDir {refdir}{species}/STARIndex_{species} \
		  --genomeFastaFiles {input.fa} \
		  --sjdbGTFfile {input.gtf} \
		  --sjdbOverhang 100 "
