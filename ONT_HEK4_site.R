# seqkit seq --threads 30 --min-len 7000 --max-len 20000 20230904-NPL2300089-P4-PAQ83671-sup.pass.fastq.gz > Passed_7_20K.fastq
# ~/minibar.py -P demultiplexed_ -p 1 -E 5 -F –T index_HEK4.txt Passed_7_20K.fastq
setwd(dir="/data/pan/Zhou_He/AE_nanopore/raw_data_2023-09-15/Nanopore/HZ-0727/20230904-NPL2300089-P4-PAQ83671-sup/cb17ded0_qc_report/")
library(tidyverse)
library(data.table)
# /usr/bin/blastn -query demultiplexed_HEK4-1641-1.fasta -task 'dc-megablast' -db duplicationArea.HEK4-1641 -out HEK4-1641-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-1641-2.fasta -task 'dc-megablast' -db duplicationArea.HEK4-1641 -out HEK4-1641-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-1641-3.fasta -task 'dc-megablast' -db duplicationArea.HEK4-1641 -out HEK4-1641-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-262-1.fasta -task 'dc-megablast' -db duplicationArea.HEK4-262 -out HEK4-262-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-262-2.fasta -task 'dc-megablast' -db duplicationArea.HEK4-262 -out HEK4-262-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-262-3.fasta -task 'dc-megablast' -db duplicationArea.HEK4-262 -out HEK4-262-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-8224-1.fasta -task 'dc-megablast' -db duplicationArea.HEK4-8224 -out HEK4-8224-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-8224-2.fasta -task 'dc-megablast' -db duplicationArea.HEK4-8224 -out HEK4-8224-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query demultiplexed_HEK4-8224-3.fasta -task 'dc-megablast' -db duplicationArea.HEK4-8224 -out HEK4-8224-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
  
HEK4_200.A.len=262
HEK4_1K.A.len=1641
HEK4_8K.A.len=8224
HEK4_200.I.len=30
HEK4_1K.I.len=30
HEK4_8K.I.len=50
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
blastn_HEK4_200_1.f<-filtered.F(blastnTab.rawFileName="HEK4-262-1.blastnRes",
                                Amp.length=HEK4_200.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 40)
blastn_HEK4_200_2.f<-filtered.F(blastnTab.rawFileName="HEK4-262-2.blastnRes",
                                Amp.length=HEK4_200.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 40)
blastn_HEK4_200_3.f<-filtered.F(blastnTab.rawFileName="HEK4-262-3.blastnRes",
                                Amp.length=HEK4_200.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 40)
blastn_HEK4_1K_1.f<-filtered.F(blastnTab.rawFileName="HEK4-1641-1.blastnRes",
                               Amp.length=HEK4_1K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 40)
blastn_HEK4_1K_2.f<-filtered.F(blastnTab.rawFileName="HEK4-1641-2.blastnRes",
                               Amp.length=HEK4_1K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 40)
blastn_HEK4_1K_3.f<-filtered.F(blastnTab.rawFileName="HEK4-1641-3.blastnRes",
                               Amp.length=HEK4_1K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 40)
blastn_HEK4_8K_1.f<-filtered.F(blastnTab.rawFileName="HEK4-8224-1.blastnRes",
                               Amp.length=HEK4_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 60)
blastn_HEK4_8K_2.f<-filtered.F(blastnTab.rawFileName="HEK4-8224-2.blastnRes",
                               Amp.length=HEK4_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 60)
blastn_HEK4_8K_3.f<-filtered.F(blastnTab.rawFileName="HEK4-8224-3.blastnRes",
                               Amp.length=HEK4_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 60)
for (i in 1:9) {
  tt<-fread(file = list.files(pattern = "blastnRes.filteredDupFreq.txt")[i])
  print(list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[i])
  tt%>%pull(V2)%>%table()%>%print()
}
# for i in *.blastnRes.filteredDupFreq.txt; do awk -F"\t" '{print$1}' $i > ${i%%Freq.txt}Seq.txt; seqkit grep -j 20 -f ${i%%Freq.txt}Seq.txt demultiplexed_${i%%blastnRes.filteredDupFreq.txt}fasta > ${i%%blastnRes.filteredDupFreq.txt}filteredDupSeq.fasta; rm ${i%%Freq.txt}Seq.txt; done
library(Biostrings)
HEK4_200_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion_HEK4",recursive=T)[2])),
                   as.character(readDNAStringSet(list.files(pattern="duplicationArea.HEK4.*.fasta",recursive=T)[2])),sep="")
for (i in c(2:28,31,32)) {
  fullRef.HEK4_200<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[2])),
    as.character(readDNAStringSet(list.files(pattern="duplicationArea.HEK4.*.fasta",recursive=T)[2])),
    paste(rep(HEK4_200_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[1])),
    sep=""))
  writeXStringSet(x=fullRef.HEK4_200,
                  filepath=paste("FullReference_HEK4_262","_",i,"A.fasta",sep=""))
}
HEK4_1K_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion_HEK4",recursive=T)[1])),
                  as.character(readDNAStringSet(list.files(pattern="duplicationArea.HEK4.*.fasta",recursive=T)[1])),sep="")
for (i in c(2:7)) {
  fullRef.HEK4_1K<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[1])),
    as.character(readDNAStringSet(list.files(pattern="duplicationArea.HEK4.*.fasta",recursive=T)[1])),
    paste(rep(HEK4_1K_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[1])),
    sep=""))
  writeXStringSet(x=fullRef.HEK4_1K,
                  filepath=paste("FullReference_HEK4_1641","_",i,"A.fasta",sep=""))
}
HEK4_8K_IA<-paste(as.character(readDNAStringSet(list.files(pattern="insertion_HEK4",recursive=T)[3])),
                  as.character(readDNAStringSet(list.files(pattern="duplicationArea.HEK4.*.fasta",recursive=T)[3])),sep="")
for (i in c(2)) {
  fullRef.HEK4_8K<-DNAStringSet(paste(
    as.character(readDNAStringSet(list.files(pattern="L.fasta",recursive=T)[3])),
    as.character(readDNAStringSet(list.files(pattern="duplicationArea.HEK4.*.fasta",recursive=T)[3])),
    paste(rep(HEK4_8K_IA,times=c(i-1)),collapse = ""),
    as.character(readDNAStringSet(list.files(pattern="R.fasta",recursive=T)[1])),
    sep=""))
  writeXStringSet(x=fullRef.HEK4_8K,
                  filepath=paste("FullReference_HEK4_8224","_",i,"A.fasta",sep=""))
}
###############################################################################
HEK4_200_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[4])
for (i in c(2:23,26,27,31,32)) {
  HEK4_200_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[4]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_200_2.DupFreq<-fread(file=list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[5])
for (i in c(2:19,23:25,31)) {
  HEK4_200_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[5]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_200_3.DupFreq<-fread(file=list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[6])
for (i in c(2:20,24,25,28,32)) {
  HEK4_200_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[6]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_1K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[1])
for (i in c(2:5,7)) {
  HEK4_1K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[1]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_1K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[2])
for (i in c(2:7)) {
  HEK4_1K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[2]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_1K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[3])
for (i in c(2:5)) {
  HEK4_1K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[3]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_8K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[7])
for (i in c(2)) {
  HEK4_8K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[7]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_8K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[8])
for (i in c(2)) {
  HEK4_8K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[8]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_8K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[9])
for (i in c(2)) {
  HEK4_8K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[9]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
# for i in *.blastnRes.filtered_[2-9]A.txt; do echo $i; seqkit grep -f $i demultiplexed_${i%%blastnRes.filtered_[2-9]A.txt}fasta > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_1[0-9]A.txt; do echo $i; seqkit grep -f $i demultiplexed_${i%%blastnRes.filtered_1[0-9]A.txt}fasta > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_2[0-9]A.txt; do echo $i; seqkit grep -f $i demultiplexed_${i%%blastnRes.filtered_2[0-9]A.txt}fasta > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_3[0-9]A.txt; do echo $i; seqkit grep -f $i demultiplexed_${i%%blastnRes.filtered_3[0-9]A.txt}fasta > ${i%txt}fasta; done
# for i in FullReference_*.fasta; do echo $i; /usr/bin/makeblastdb -in $i -dbtype 'nucl' -out ${i%%.fasta}; done
###blastn#####
# for i in HEK4-262-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_9A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_9A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_10A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_10A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_11A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_11A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_12A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_12A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_13A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_13A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_14A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_14A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_15A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_15A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_16A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_16A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_17A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_17A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_18A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_18A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_19A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_19A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1,3].blastnRes.filtered_20A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_20A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-1.blastnRes.filtered_21A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_21A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-1.blastnRes.filtered_22A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_22A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1,2].blastnRes.filtered_23A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_23A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[2,3].blastnRes.filtered_24A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_24A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[2,3].blastnRes.filtered_25A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_25A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-1.blastnRes.filtered_26A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_26A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-1.blastnRes.filtered_27A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_27A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-3.blastnRes.filtered_28A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_28A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1,2].blastnRes.filtered_31A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_31A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1,3].blastnRes.filtered_32A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_262_32A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_1641_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_1641_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_1641_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_1641_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-2.blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_1641_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1,2].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_1641_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-8224-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db FullReference_HEK4_8224_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

for (i in c(1:10,12:26,28:72,74:89)) {
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

# for i in HEK4-262-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-262-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-262-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-262-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-262-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-262-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-1641-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-1641-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-1641-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-1641-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-1641-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-1641-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-8224-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-8224-1.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-8224-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-8224-2.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-8224-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i demultiplexed_HEK4-8224-3.fasta > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done



