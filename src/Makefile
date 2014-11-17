name=extract
sam_src=~/tools-src/samtools/samtools-0.1.19
.PHONY: all

all: $(name)

$(name): $(name).cpp
	g++ -O3 -Wall -I $(sam_src) -L $(sam_src) $(name).cpp -lm -lbam -lpthread -lz -o $(name)
	chmod 755 $(name) 
