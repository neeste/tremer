# Makefile for ydna files 
CC=gcc
CFLAGS= -Wall
GFLAGS= 
FTDNA=~/FTDNA
PROJ=../Neely_project
YDNA=../Neely_project
YSTR=$(FTDNA)/YSTR
PGMS= strmerge strdata snptree
STRS= strdata1.txt strdata2.txt strdata3.txt \
      strdata4.txt strdata5.txt strdata6.txt strdata7.txt
DATS= Neely_s1.txt Neely_s2.txt Neely_s12.txt Neely_s3.txt \
      Neely_s4.txt Neely_s5.txt Neely_s6.txt Neely_s7.txt
.PHONY: build snp update dist

build : $(PGMS) $(DATS) snp

snp : snptree
	./snptree $(DATS)

update : 
	python3 update_tremer_ages.py ../Neely_project/Neely_s*.txt
	cp Neely_s*.txt ../Neely_project

dist :
	rm -rf build dist update_tremer_ages.spec
	pyinstaller --onefile update_tremer_ages.py
	
Neely_s1.txt : strdata1.txt 
	./strdata -o $@ $(YDNA)/Neely_s1.txt 
	cp $@ $(PROJ) 

Neely_s2.txt : strdata2.txt 
	./strdata -o $@ $(YDNA)/Neely_s2.txt 
	cp $@ $(PROJ) 

Neely_s12.txt : strdata12.txt 
	./strdata -o $@ $(YDNA)/Neely_s12.txt 
	cp $@ $(PROJ) 

Neely_s3.txt : strdata3.txt 
	./strdata -o $@ $(YDNA)/Neely_s3.txt 
	cp $@ $(PROJ) 

Neely_s4.txt : strdata4.txt 
	./strdata -o $@ $(YDNA)/Neely_s4.txt 
	cp $@ $(PROJ) 

Neely_s5.txt : strdata5.txt 
	./strdata -o $@ $(YDNA)/Neely_s5.txt 
	cp $@ $(PROJ) 

Neely_s6.txt : strdata6.txt 
	./strdata -o $@ $(YDNA)/Neely_s6.txt 
	cp $@ $(PROJ) 

Neely_s7.txt : strdata7.txt 
	./strdata -o $@ $(YDNA)/Neely_s7.txt 
	cp $@ $(PROJ) 

strdata1.txt : $(YSTR)/group_1/*.csv
	./strmerge -o $@ $^

strdata2.txt : $(YSTR)/group_2/*.csv
	./strmerge -o $@ $^

strdata3.txt : $(YSTR)/group_3/*.csv 
	./strmerge -o $@ $^

strdata4.txt : $(YSTR)/group_4/*.csv 
	./strmerge -o $@ $^

strdata5.txt : $(YSTR)/group_5/*.csv 
	./strmerge -o $@ $^

strdata6.txt : $(YSTR)/group_6/*.csv 
	./strmerge -o $@ $^

strdata7.txt : $(YSTR)/group_7/*.csv 
	./strmerge -o $@ $^

strdata12.txt : strdata1.txt strdata2.txt
	cat strdata1.txt strdata2.txt > strdata12.txt

clean:
	rm -f *.o $(PGMS) *.exe *~ 
	rm -f Neely*.{txt,svg,png} {str,cal,snpy}*.txt 
	rm -f tree[1-6].txt 
	rm -rf build dist update_tremer_ages.spec

