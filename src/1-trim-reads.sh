###.............................................................................
# (c) Miguel Camacho Sánchez
# miguelcamachosanchez@gmail.com // miguelcamachosanchez.weebly.com
# https://scholar.google.co.uk/citations?user=1M02-S4AAAAJ&hl=en
# June 2020
###.............................................................................
#GOAL: demultiplex amplicon data based on loci primers
#PROJECT: amplicon-genotyping
###.............................................................................
#batch renaming after patterns in data/raw/id-match
mmv < data/raw/id-match

#write list of samples
ls data/raw/*fastq* | sed -e 's|.[12].fastq.gz||;s|data/raw/||' | sort | uniq -c \
| grep '^.*2 ' | awk '{print $2}' > data/intermediate/samples-list

#dual demultiplexing of reads.
#dual demultiplexing is only enabled from cutadapt v2.1
cat data/intermediate/samples-list | while read sample
do
  cutadapt \
      -e 0.15 --no-indels \
      --discard-untrimmed \
      --pair-adapters \
      -q 10 \
      -m 150 \
      -g file:data/intermediate/forward \
      -G file:data/intermediate/reverse \
      -o data/intermediate/$sample.{name}.1.fastq.gz \
      -p data/intermediate/$sample.{name}.2.fastq.gz \
      data/raw/$sample.1.fastq.gz \
      data/raw/$sample.2.fastq.gz
done > output/cutadapt.txt
