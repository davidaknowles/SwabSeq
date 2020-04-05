SHELL := /bin/bash
.DELETE_ON_ERROR:

#===============================================================================
# VARIABLES

# sequencing run ids
RUNS := $(shell find data -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)

#===============================================================================
# RECIPIES

all: conds star
conds: $(addprefix pipeline/, $(addsuffix /conditions.csv, $(RUNS)))
star: $(addprefix pipeline/, $(addsuffix /starcode.csv, $(RUNS)))

# cleanup
clean:
	rm -f pipeline/*
deep_clean:
	rm -f pipeline/* results/*
.PRECIOUS: $(addprefix pipeline/, %/conditions.csv %/starcode.csv)
.SECONDARY:

#===============================================================================
# Run starcode on each fastq, using gnu parallel to parse the SampleSheet.csv
pipeline/%/starcode.csv: data/% pipeline/%/conditions.csv
	@echo "Counting BCs for all fastq's in $<"
	@parallel --header : --colsep "," \
		zcat $</"{Sample_ID}"_S*_R1_001.fastq.gz \
		\| awk -v bc_len="{bc_len}" -f src/count-bcs.awk \
		\| starcode -d2 -t1 --sphere --print-clusters 2> /dev/null \
		\| python src/tidy-star.py \
		\| awk -v name="{Sample_ID}" \''{print name, $$1, $$2, $$3}'\' OFS="," \
	:::: $(word 2, $^) 2> $(@:.csv=.err) > $(@:.csv=.tmp) \
	&& echo "Sample_ID,Centroid,Count,barcode" \
	| cat - $(@:.csv=.tmp) > $@ \
	&& rm $(@:.csv=.tmp)

# grab relevant section of samplesheet (make sure to catch windows return)
pipeline/%/conditions.csv: data/%/SampleSheet.csv pipeline/%/bc-map.csv
	@echo "Parsing $<"
	@python src/strip-windows.py $< \
		| awk '/Sample_ID/{seen=1} seen{print}' \
		| Rscript src/bc-lengths.R $(word 2, $^) $@ \
		2> $(@:.csv=.err)

# copy barcode map from data/ folder.
pipeline/%/bc-map.csv:
	@echo "Grabbing $@"
	@mkdir -p $(dir $@)
	@cp data/barcode-map.csv $@

