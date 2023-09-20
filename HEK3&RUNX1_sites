# seqkit seq --threads 30 --min-len 7000 20230713-NPL230912-P7-PAQ57405-sup.pass.fastq.gz > Passed_7Kmore.fastq
# ~/minibar.py -P index4minibar6_HEK3_trim -p 1 -E 5 -F –T index4minibar6_HEK3.txt Passed_7Kmore.fastq 
# ~/minibar.py -P index4minibar6_RUNX1_trim -p 1 -E 5 -F -T index4minibar6_RUNX1.txt Passed_7Kmore.fastq
setwd(dir = "/data/pan/Zhou_He/AE_nanopore/raw_data_2023-07-31/Nanopore/HZ-0630/20230713-NPL230912-P7-PAQ57405-sup/5e093197_qc_report/Passed_7Kmore_demultiplexed/Blastn/")
library(tidyverse)
library(data.table)
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-1K-1.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-1K -out index4minibar6_HEK3_trimHEK3-1K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-1K-2.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-1K -out index4minibar6_HEK3_trimHEK3-1K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-1K-3.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-1K -out index4minibar6_HEK3_trimHEK3-1K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-200-1.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-200 -out index4minibar6_HEK3_trimHEK3-200-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-200-2.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-200 -out index4minibar6_HEK3_trimHEK3-200-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-200-3.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-200 -out index4minibar6_HEK3_trimHEK3-200-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-8K-1.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-8K -out index4minibar6_HEK3_trimHEK3-8K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-8K-2.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-8K -out index4minibar6_HEK3_trimHEK3-8K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_HEK3_trimHEK3-8K-3.fasta -task 'dc-megablast' -db ../../duplicationArea.HEK3-8K -out index4minibar6_HEK3_trimHEK3-8K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-1K-1.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-1K -out index4minibar6_RUNX1_trimRUNX1-1K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-1K-2.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-1K -out index4minibar6_RUNX1_trimRUNX1-1K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-1K-3.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-1K -out index4minibar6_RUNX1_trimRUNX1-1K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-200-1.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-200 -out index4minibar6_RUNX1_trimRUNX1-200-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-200-2.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-200 -out index4minibar6_RUNX1_trimRUNX1-200-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-200-3.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-200 -out index4minibar6_RUNX1_trimRUNX1-200-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-8K-1.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-8K -out index4minibar6_RUNX1_trimRUNX1-8K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-8K-2.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-8K -out index4minibar6_RUNX1_trimRUNX1-8K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
# /usr/bin/blastn -query index4minibar6_RUNX1_trimRUNX1-8K-3.fasta -task 'dc-megablast' -db ../../duplicationArea.RUNX1-8K -out index4minibar6_RUNX1_trimRUNX1-8K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 #
RUNX1_200.length=269
RUNX1_1K.length=1132
RUNX1_8K.length=7677
HEK3_200.length=234
HEK3_1K.length=1098
HEK3_8K.length=8207
RUNX1_200.insert.length=30
RUNX1_1K.insert.length=30
RUNX1_8K.insert.length=50
HEK3_200.insert.length=20
HEK3_1K.insert.length=50
HEK3_8K.insert.length=50
filtered.F<-function(blastnTab.rawFileName=NULL,Amp.length=NULL,
                     threshold.length=0.90,threshold.pident=90,
                     distance.threshold=NULL) {
  blastnTab.raw<-fread(file = blastnTab.rawFileName,skip = 1)
  colnames(blastnTab.raw)<-c("qaccver","saccver","pident","length","mismatch",
                             "gapopen","qstart","qend","sstart","send",
                             "evalue","bitscore","qlen","qcovhsp","qcovus")
  blastn.Tab<-blastnTab.raw%>%filter(length>threshold.length*Amp.length,
                                     pident>threshold.pident)
  blastn.Tab.dupNames<-as.data.frame(table(blastn.Tab$qaccver))%>%
    filter(Freq>1)%>%.[,1]%>%as.character()
  blastn.Tab<-blastn.Tab[blastn.Tab$qaccver%in%blastn.Tab.dupNames,]
  #nrow
  read_nrow<-blastn.Tab%>%group_by(qaccver)%>%
    mutate(strand=ifelse(sstart<send,1,-1))%>%
    arrange(qstart,.by_group = T)%>%summarise(read.nrow=n())
  #strand.sum
  strand_sum<-blastn.Tab%>%group_by(qaccver)%>%
    mutate(strand=ifelse(sstart<send,1,-1))%>%
    arrange(qstart,.by_group = T)%>%summarise(strand.sum=sum(strand))
  remove1<-cbind(read_nrow,strand_sum[,2])%>%
    mutate(m=ifelse(read.nrow==abs(strand.sum),1,0))%>%
    filter(m==0)%>%select(qaccver)%>%unique()
  #distances
  remove2<-blastn.Tab%>%group_by(qaccver)%>%
    mutate(strand=ifelse(sstart<send,1,-1))%>%arrange(qstart,.by_group = T)%>%
    mutate(distance=qstart-lag(qend))%>%filter(distance>distance.threshold)%>%
    select(qaccver)%>%unique()
  removeRead<-union(remove1,remove2)
  blastn.Tab.f<-blastn.Tab%>%filter(!qaccver%in%removeRead$qaccver)
  #short alignment
  # blastn.Tab.SA<-blastnTab.raw%>%filter(pident>threshold.pident,
  #                                       length<threshold.length*Amp.length)%>%
  #   select(qaccver)%>%unique()
  # blastn.Tab.SAcheck<-blastnTab.raw%>%
  #   filter(qaccver%in%intersect(blastn.Tab.f$qaccver,blastn.Tab.SA$qaccver))
  #short alignment #nrow
  # read_nrow.SAcheck<-blastn.Tab.SAcheck%>%group_by(qaccver)%>%
  #   mutate(strand=ifelse(sstart<send,1,-1))%>%
  #   arrange(qstart,.by_group = T)%>%summarise(read.nrow=n())
  #short alignment #strand.sum
  # strand_sum.SAcheck<-blastn.Tab.SAcheck%>%group_by(qaccver)%>%
  #   mutate(strand=ifelse(sstart<send,1,-1))%>%
  #   arrange(qstart,.by_group = T)%>%summarise(strand.sum=sum(strand))
  # remove1.SAcheck<-cbind(read_nrow.SAcheck,strand_sum.SAcheck[,2])%>%
  #   mutate(m=ifelse(read.nrow==abs(strand.sum),1,0))%>%
  #   filter(m==0)%>%select(qaccver)%>%unique()
  #short alignment #distances
  # remove2.SAcheck<-blastn.Tab.SAcheck%>%group_by(qaccver)%>%
  #   mutate(strand=ifelse(sstart<send,1,-1))%>%arrange(qstart,.by_group = T)%>%
  #   mutate(distance=qstart-lag(qend))%>%filter(distance>distance.threshold)%>%
  #   select(qaccver)%>%unique()
  # removeRead.SAcheck<-union(remove1.SAcheck,remove2.SAcheck)
  # blastn.Tab.SA.removed<-blastn.Tab.SAcheck%>%filter(!qaccver%in%removeRead.SAcheck$qaccver)
  # blastn.Tab.f2<-blastn.Tab.f%>%filter(!qaccver%in%blastn.Tab.SA.removed$qaccver)
  print(blastnTab.rawFileName)
  blastn.Tab.f2<-blastn.Tab.f
  blastn.Tab.f2$qaccver%>%unique()%>%length()%>%print()
  as.data.frame(table(blastn.Tab.f2$qaccver))%>%filter(Freq>1)%>%.[,2]%>%summary()%>%print()
  blastn.Tab.f2$qaccver%>%table()%>%sort(decreasing = T)%>%as.data.frame()%>%
    write.table(file=paste(blastnTab.rawFileName,".filteredDupFreq.txt",sep = ""),
                quote = F,sep = "\t",row.names = F,col.names = F)
  return(blastn.Tab.f2)
}
blastn_HEK3_200_1.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-200-1.blastnRes",
                                Amp.length=HEK3_200.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_200.insert.length+10))
blastn_HEK3_200_2.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-200-2.blastnRes",
                                Amp.length=HEK3_200.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_200.insert.length+10))
blastn_HEK3_200_3.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-200-3.blastnRes",
                                Amp.length=HEK3_200.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_200.insert.length+10))
blastn_HEK3_1K_1.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-1K-1.blastnRes",
                                Amp.length=HEK3_1K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_1K.insert.length+10))
blastn_HEK3_1K_2.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-1K-2.blastnRes",
                                Amp.length=HEK3_1K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_1K.insert.length+10))
blastn_HEK3_1K_3.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-1K-3.blastnRes",
                                Amp.length=HEK3_1K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_1K.insert.length+10))
blastn_HEK3_8K_1.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-8K-1.blastnRes",
                               Amp.length=HEK3_8K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_8K.insert.length+10))
blastn_HEK3_8K_2.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-8K-2.blastnRes",
                               Amp.length=HEK3_8K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_8K.insert.length+10))
blastn_HEK3_8K_3.f<-filtered.F(blastnTab.rawFileName="index4minibar6_HEK3_trimHEK3-8K-3.blastnRes",
                               Amp.length=HEK3_8K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(HEK3_8K.insert.length+10))

blastn_RUNX1_200_1.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-200-1.blastnRes",
                                 Amp.length=RUNX1_200.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_200.insert.length+10))
blastn_RUNX1_200_2.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-200-2.blastnRes",
                                 Amp.length=RUNX1_200.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_200.insert.length+10))
blastn_RUNX1_200_3.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-200-3.blastnRes",
                                 Amp.length=RUNX1_200.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_200.insert.length+10))
blastn_RUNX1_1K_1.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-1K-1.blastnRes",
                                Amp.length=RUNX1_1K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_1K.insert.length+10))
blastn_RUNX1_1K_2.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-1K-2.blastnRes",
                                Amp.length=RUNX1_1K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_1K.insert.length+10))
blastn_RUNX1_1K_3.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-1K-3.blastnRes",
                                Amp.length=RUNX1_1K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_1K.insert.length+10))
blastn_RUNX1_8K_1.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-8K-1.blastnRes",
                                Amp.length=RUNX1_8K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_8K.insert.length+10))
blastn_RUNX1_8K_2.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-8K-2.blastnRes",
                                Amp.length=RUNX1_8K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_8K.insert.length+10))
blastn_RUNX1_8K_3.f<-filtered.F(blastnTab.rawFileName="index4minibar6_RUNX1_trimRUNX1-8K-3.blastnRes",
                                Amp.length=RUNX1_8K.length,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_8K.insert.length+10))

for (i in 1:18) {
  tt<-fread(file = list.files(pattern = "blastnRes.filteredDupFreq.txt")[i])
  print(list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[i])
  tt%>%pull(V2)%>%table()%>%print()
}
# for i in *.blastnRes.filteredDupFreq.txt; do awk -F"\t" '{print$1}' $i > ${i%%Freq.txt}Seq.txt; seqkit grep -j 20 -f ${i%%Freq.txt}Seq.txt ${i%%blastnRes.filteredDupFreq.txt}fasta > ${i%%blastnRes.filteredDupFreq.txt}filteredDupSeq.fasta; rm ${i%%Freq.txt}Seq.txt; done

library(Biostrings)
HEK3_200_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion.fasta",recursive=T)[2])),
                   as.character(readDNAStringSet("../../duplicationArea.HEK3-200.fasta")),sep="")
for (i in c(2:21,23:29,34)) {
  fullRef.HEK3_200<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[1])),
    as.character(readDNAStringSet("../../duplicationArea.HEK3-200.fasta")),
    paste(rep(HEK3_200_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[2])),
    sep=""))
  writeXStringSet(x=fullRef.HEK3_200,
                  filepath=paste("FullReference_HEK3_200","_",i,"A.fasta",sep=""))
}
RUNX1_200_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion.fasta",recursive=T)[5])),
                    as.character(readDNAStringSet("../../duplicationArea.RUNX1-200.fasta")),sep="")
for (i in c(2:21,23)) {
  fullRef.RUNX1_200<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[2])),
    as.character(readDNAStringSet("../../duplicationArea.RUNX1-200.fasta")),
    paste(rep(RUNX1_200_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[5])),
    sep=""))
  writeXStringSet(x=fullRef.RUNX1_200,
                  filepath=paste("FullReference_RUNX1_200","_",i,"A.fasta",sep=""))
}
HEK3_1K_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion.fasta",recursive=T)[1])),
                  as.character(readDNAStringSet("../../duplicationArea.HEK3-1K.fasta")),sep="")
for (i in c(2:8,10,11)) {
  fullRef.HEK3_1K<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[1])),
    as.character(readDNAStringSet("../../duplicationArea.HEK3-1K.fasta")),
    paste(rep(HEK3_1K_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[1])),
    sep=""))
  writeXStringSet(x=fullRef.HEK3_1K,
                  filepath=paste("FullReference_HEK3_1K","_",i,"A.fasta",sep=""))
}
RUNX1_1K_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion.fasta",recursive=T)[4])),
                   as.character(readDNAStringSet("../../duplicationArea.RUNX1-1K.fasta")),sep="")
for (i in 2:8) {
  fullRef.RUNX1_1K<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[2])),
    as.character(readDNAStringSet("../../duplicationArea.RUNX1-1K.fasta")),
    paste(rep(RUNX1_1K_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[4])),
    sep=""))
  writeXStringSet(x=fullRef.RUNX1_1K,
                  filepath=paste("FullReference_RUNX1_1K","_",i,"A.fasta",sep=""))
}
HEK3_8K_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion.fasta",recursive=T)[3])),
                  as.character(readDNAStringSet("../../duplicationArea.HEK3-8K.fasta")),sep="")
for (i in c(2)) {
  fullRef.HEK3_8K<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[1])),
    as.character(readDNAStringSet("../../duplicationArea.HEK3-8K.fasta")),
    paste(rep(HEK3_8K_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[3])),
    sep=""))
  writeXStringSet(x=fullRef.HEK3_8K,
                  filepath=paste("FullReference_HEK3_8K","_",i,"A.fasta",sep=""))
}
RUNX1_8K_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion.fasta",recursive=T)[6])),
                   as.character(readDNAStringSet("../../duplicationArea.RUNX1-8K.fasta")),sep="")
for (i in c(2)) {
  fullRef.RUNX1_8K<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[2])),
    as.character(readDNAStringSet("../../duplicationArea.RUNX1-8K.fasta")),
    paste(rep(RUNX1_8K_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[6])),
    sep=""))
  writeXStringSet(x=fullRef.RUNX1_8K,
                  filepath=paste("FullReference_RUNX1_8K","_",i,"A.fasta",sep=""))
}
#########################################
HEK3_200_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[4])
for (i in c(2:15,17:19,21,24,25,34)) {
  HEK3_200_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[4]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_200_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[5])
for (i in c(2:17,19:21,24,25,27,29)) {
  HEK3_200_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[5]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_200_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[6])
for (i in c(2:17,20,23:24,26:29)) {
  HEK3_200_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[6]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_1K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[1])
for (i in c(2:8)) {
  HEK3_1K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[1]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_1K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[2])
for (i in c(2:8,10,11)) {
  HEK3_1K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[2]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_1K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[3])
for (i in c(2:6)) {
  HEK3_1K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[3]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_8K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[7])
for (i in c(2)) {
  HEK3_8K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[7]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_8K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[8])
for (i in c(2)) {
  HEK3_8K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[8]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_8K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[9])
for (i in c(2)) {
  HEK3_8K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[9]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_200_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[13])
for (i in c(2:21,23)) {
  RUNX1_200_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[13]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_200_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[14])
for (i in c(2:14,16:19)) {
  RUNX1_200_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[14]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_200_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[15])
for (i in c(2:15,20)) {
  RUNX1_200_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[15]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_1K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[10])
for (i in c(2:8)) {
  RUNX1_1K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[10]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_1K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[11])
for (i in c(2:7)) {
  RUNX1_1K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[11]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_1K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[12])
for (i in c(2:5)) {
  RUNX1_1K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[12]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_8K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[16])
for (i in c(2)) {
  RUNX1_8K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[16]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_8K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[17])
for (i in c(2)) {
  RUNX1_8K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[17]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_8K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[18])
for (i in c(2)) {
  RUNX1_8K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[18]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
# for i in *.blastnRes.filtered_[2-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_[2-9]A.txt}fasta > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_1[0-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_1[0-9]A.txt}fasta > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_2[0-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_2[0-9]A.txt}fasta > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_3[0-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_3[0-9]A.txt}fasta > ${i%txt}fasta; done
# for i in FullReference_*.fasta; do echo $i; /usr/bin/makeblastdb -in $i -dbtype 'nucl' -out ${i%%.fasta}; done
###blastn#####
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_9A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_9A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_10A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_10A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_11A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_11A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_12A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_12A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_13A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_13A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_14A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_14A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_15A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_15A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[2,3].blastnRes.filtered_16A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_16A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_17A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_17A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-1.blastnRes.filtered_18A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_18A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1,2].blastnRes.filtered_19A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_19A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[2,3].blastnRes.filtered_20A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_20A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1,2].blastnRes.filtered_21A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_21A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-3.blastnRes.filtered_23A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_23A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1-3].blastnRes.filtered_24A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_24A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[1,2].blastnRes.filtered_25A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_25A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-3.blastnRes.filtered_26A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_26A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[2,3].blastnRes.filtered_27A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_27A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-3.blastnRes.filtered_28A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_28A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-[2,3].blastnRes.filtered_29A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_29A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-200-1.blastnRes.filtered_34A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_200_34A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

# for i in index4minibar6_HEK3_trimHEK3-1K-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-[1,2].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-[1,2].blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-2.blastnRes.filtered_10A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_10A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_HEK3_trimHEK3-1K-2.blastnRes.filtered_11A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_1K_11A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

# for i in index4minibar6_HEK3_trimHEK3-8K-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK3_8K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_9A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_9A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_10A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_10A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_11A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_11A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_12A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_12A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_13A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_13A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1-3].blastnRes.filtered_14A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_14A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1,3].blastnRes.filtered_15A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_15A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1,2].blastnRes.filtered_16A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_16A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1,2].blastnRes.filtered_17A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_17A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1,2].blastnRes.filtered_18A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_18A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1,2].blastnRes.filtered_19A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_19A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-[1,3].blastnRes.filtered_20A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_20A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-1.blastnRes.filtered_21A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_21A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-1.blastnRes.filtered_23A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_200_23A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

# for i in index4minibar6_RUNX1_trimRUNX1-1K-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_1K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_1K_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_1K_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_1K_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-[1,2].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_1K_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-[1,2].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_1K_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-1.blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_1K_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# 
# for i in index4minibar6_RUNX1_trimRUNX1-8K-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_RUNX1_8K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

for (i in c(1:164)) {
  rr<-fread(file = list.files(pattern = "A.blastnRes$",recursive = F)[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  print(list.files(pattern = "A.blastnRes",recursive = T)[i])
  # rr%>%mutate(D=slen-qlen)%>%mutate(D2=slen-length)%>%filter(abs(D2)<100)%>%
  #   filter(mismatch<100)%>%arrange(desc(abs(D2)))%>%head()%>%print()
  rr.f<-rr%>%mutate(D=slen-qlen)%>%mutate(D2=slen-length)%>%
    mutate(D3=qlen-length)%>%filter(abs(D2)<0.02*slen)%>%
    filter(abs(D3)<0.02*slen)%>%filter(mismatch<0.02*slen)
  rr.f%>%pull(qaccver)%>%unique()%>%length()%>%print()
  rr.f%>%pull(qaccver)%>%unique()%>%as.data.frame()%>%
    write.table(file=paste(list.files(pattern="A.blastnRes$",recursive=F)[i],".filteredRead.txt",sep=""),
                quote = F,sep = "\t",row.names = F,col.names = F)
}

# for i in index4minibar6_HEK3_trimHEK3-200-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-200-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-200-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-200-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-200-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-200-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-1K-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-1K-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-1K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-1K-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-1K-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-1K-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-8K-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-8K-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-8K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-8K-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_HEK3_trimHEK3-8K-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_HEK3_trimHEK3-8K-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done

# for i in index4minibar6_RUNX1_trimRUNX1-200-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-200-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-200-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-200-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-200-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-1K-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-1K-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-1K-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-1K-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-8K-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-8K-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-8K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-8K-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in index4minibar6_RUNX1_trimRUNX1-8K-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i index4minibar6_RUNX1_trimRUNX1-8K-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done

