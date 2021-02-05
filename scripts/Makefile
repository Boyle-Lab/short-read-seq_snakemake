SAMPLEINFO = /data/projects/protamine/1619-SS/runinfo.tsv
ROOTDIR = /data/projects/protamine/1619-SS
GENES = /data/projects/protamine/data/mm10_refGene.bed
BAMDIR = ${ROOTDIR}/prune
MACSDIR = ${ROOTDIR}/fseq2
BWDIR = ${ROOTDIR}/bigwig
RESDIR = ${ROOTDIR}/results

CPUS=12

NAME = 1619-SS

SAMPLES = $$(for file in $$(ls ${BAMDIR}/*.rep*.pruned.bam | grep -v 'Input'); do basename $$file .pruned.bam; done)
BAMFILES = $$(ls ${BAMDIR}/*.pruned.bam)
SORTED_BAMS = $$(for file in $$(ls ${BAMDIR}/*.rep*.pruned.bam | grep -v 'Input'); do printf "%s.pruned.sorted.bam\n" $$(basename $$file .pruned.bam); done)
SORTED_INPUT = $$(echo $$(for SAMPLE in $(SAMPLES); do ls $(BAMDIR)/$$(echo $$SAMPLE | awk -F '.' '{print $$1}').Input.*.pruned.bam | sed 's/pruned\.bam/pruned.sorted.bam/'; done) | sed 's/ /\n/g' | sort | uniq)
BW = $$(for sample in ${SAMPLES}; do printf "%s/%s.bw " ${BWDIR} $$sample; done)

OTHERFILES = /data/projects/protamine/18342R/bigwig/P2V5_7-14-20.rep1.bw /data/projects/protamine/18342R/bigwig/P2V5_7-14-20.rep2.bw /data/projects/protamine/18342R/bigwig/P1V5_7-24-20.rep1.bw
LABELS = $(SAMPLES) P2V5_7-14-20.rep1 P2V5_7-14-20.rep2 P1V5_7-24-20.rep1

all: link_bam_files sort foldenrich summary correlation heatmap

peaks: link_bam_files sort foldenrich

matrices: summary correlation heatmap

test: $($SAMPLES)
	for SAMPLE in $(SAMPLES); do \
		ls $(BAMDIR)/$$SAMPLE.pruned.bam;\
		for ctrl in $$(ls $(BAMDIR)/$$(echo $$SAMPLE | awk -F '.' '{print $$1}').Input.*.pruned.bam) ; do\
			if [ ! -f $$ctrl ] ; then \
				echo $$ctrl ;\
			fi;\
		done;\
	done ;\
	ls $(OTHERFILES)

test_inputs:
	for input in $(SORTED_INPUT); do \
		echo $$input ;\
	done

test_labels:
	for input in $(LABELS); do \
		echo $$input ;\
	done

link_bam_fles:
	ls $(BAMFILES) > /dev/null 2> /dev/null
	if [ $? != 0 ]; then \
		bash scripts/link_bams.sh $(SAMPLEINFO) $(BAMDIR) ;\
	else \
		echo "Links to bam files are already present. Nothing to do." ;\
	fi


sort: $($BAMFILES)
	for SAMPLE in $(SAMPLES); do \
		echo $$SAMPLE ;\
		if [ ! -f $(BAMDIR)/$$SAMPLE.pruned.sorted.bam ] ; then \
			printf "\tSorting %s into %s.\n" $(BAMDIR)/$$SAMPLE.pruned.bam $(BAMDIR)/$$SAMPLE.pruned.sorted.bam ; \
			samtools sort -@ $(CPUS) -n -o $(BAMDIR)/$$SAMPLE.pruned.sorted.bam $(BAMDIR)/$$SAMPLE.pruned.bam ;\
		else \
			printf "\t%s already exists. Nothing to be done.\n" $(BAMDIR)/$$SAMPLE.pruned.sorted.bam ;\
		fi ;\
	done ;\
	for SI in $(SORTED_INPUT); do\
		if [ ! -f $$SI ] ; then \
			printf "Sorting %s into %s.\n\t" $$(echo $$SI | sed 's/\.sorted//') $$SI ;\
			samtools sort -@ $(CPUS) -n -o $$SI $$(echo $$SI | sed 's/\.sorted//') ;\
		else \
			printf "%s already exists. Nothing to be done.\n" $$SI ;\
		fi ;\
	done


foldenrich: $($SORTED_BAMS) $($SORTED_INPUT)
	for SAMPLE in $(SAMPLES); do \
		echo $$SAMPLE ; \
		if [ ! -f $(MACSDIR)/$$SAMPLE.bw ] ; then \
			fseq2 callpeak \
				-control_file $$(ls $(BAMDIR)/$$(echo $$SAMPLE | awk -F '.' '{print $$1}').Input.*.pruned.sorted.bam) \
				-pe -o $(MACSDIR) -cpus $(CPUS)\
				-chrom_size_file /data/genomes/mm10/mm10.chrom.sizes \
				-sig_format bigwig \
				-name $$SAMPLE \
				$(BAMDIR)/$$SAMPLE.pruned.sorted.bam ;\
		fi ;\
	done ;\
	mv $(MACSDIR)/*.bw $(BWDIR)

summary: $($BW) $($OTHERFILES)
	multiBigwigSummary bins -b $(BW) $(OTHERFILES) -o $(RESDIR)/$(NAME).summary.npz

correlation: $(RESDIR)/$(NAME).summary.npz
	plotCorrelation --corData $(RESDIR)/$(NAME).summary.npz \
		--labels $(LABELS) --corMethod spearman \
		--skipZeros --removeOutliers --whatToPlot heatmap \
		--plotFile $(RESDIR)/$(NAME).summary.pdf

heatmap: $($BW)
	computeMatrix scale-regions -S $(BW) \
		-R $(GENES) \
		-b 2000 -a 2000 -m 1000 \
		-o $(RESDIR)/$(NAME).matrix.refgene.gz
	plotHeatmap --matrixFile $(RESDIR)/$(NAME).matrix.refgene.gz \
		--outFileName $(RESDIR)/$(NAME).matrix.refgene.pdf \
		--outFileSortedRegions $(RESDIR)/$(NAME).matrix.refgene.sr

cleanall: cleansort cleanfe cleanbigwigs cleanmatrices cleanplots

cleansort:
	-rm -f $(BAMDIR)/*.sorted.bam

cleanfe:
	-rm -f $(MACSDIR)/*

cleanbigwigs:
	-rm -f $(BWDIR)/*.bw

cleanmatrices:
	-rm -f $(RESDIR)/$(NAME).summary.npz *.matrix.refgene.gz

cleanplots:
	-rm -f $(RESDIR)/*.matrix.refgene.pdf $(RESDIR)/*.matrix.refgene.sr
	-rm -f $(RESDIR)/$(NAME).summary.pdf
