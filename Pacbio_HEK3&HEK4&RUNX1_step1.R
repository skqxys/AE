setwd(dir = "/data/pan/Zhou_He/HZ-0129/")
# for i in duplicationArea.*.fasta; do /usr/bin/makeblastdb -dbtype nucl -in $i; echo $i "done"; done
# cd demultiplexed_V2/
# /usr/bin/blastn -query HEK4-1641-1.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-1641.fasta -out HEK4-1641-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-1641-2.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-1641.fasta -out HEK4-1641-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-1641-3.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-1641.fasta -out HEK4-1641-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-262-1.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-262.fasta -out HEK4-262-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-262-2.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-262.fasta -out HEK4-262-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-262-3.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-262.fasta -out HEK4-262-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-8224-1.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-8224.fasta -out HEK4-8224-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-8224-2.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-8224.fasta -out HEK4-8224-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK4-8224-3.fa -task 'dc-megablast' -db ../duplicationArea.HEK4-8224.fasta -out HEK4-8224-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-1K-1.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-1K.fasta -out HEK3-1K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-1K-2.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-1K.fasta -out HEK3-1K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-1K-3.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-1K.fasta -out HEK3-1K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-200-1.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-200.fasta -out HEK3-200-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-200-2.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-200.fasta -out HEK3-200-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-200-3.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-200.fasta -out HEK3-200-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-8K-1.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-8K.fasta -out HEK3-8K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-8K-2.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-8K.fasta -out HEK3-8K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query HEK3-8K-3.fa -task 'dc-megablast' -db ../duplicationArea.HEK3-8K.fasta -out HEK3-8K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-1K-1.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-1K.fasta -out RUNX1-1K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-1K-2.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-1K.fasta -out RUNX1-1K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-1K-3.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-1K.fasta -out RUNX1-1K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-200-1.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-200.fasta -out RUNX1-200-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-200-2.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-200.fasta -out RUNX1-200-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-200-3.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-200.fasta -out RUNX1-200-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-8K-1.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-8K.fasta -out RUNX1-8K-1.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-8K-2.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-8K.fasta -out RUNX1-8K-2.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
# /usr/bin/blastn -query RUNX1-8K-3.fa -task 'dc-megablast' -db ../duplicationArea.RUNX1-8K.fasta -out RUNX1-8K-3.blastnRes -outfmt '6 std qlen qcovhsp qcovus' -perc_identity 80 -num_threads 20 &
setwd(dir = "/data/pan/Zhou_He/HZ-0129/demultiplexed_V2/")
library(tidyverse)
library(data.table)
RUNX1_200.A.len=269
RUNX1_1K.A.len=1132
RUNX1_8K.A.len=7677
HEK3_200.A.len=234
HEK3_1K.A.len=1098
HEK3_8K.A.len=8207
HEK4_200.A.len=262
HEK4_1K.A.len=1641
HEK4_8K.A.len=8224
RUNX1_200.I.len=30
RUNX1_1K.I.len=30
RUNX1_8K.I.len=50
HEK3_200.I.len=20
HEK3_1K.I.len=50
HEK3_8K.I.len=50
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
  print(blastnTab.rawFileName)
  blastn.Tab.f2<-blastn.Tab.f
  blastn.Tab.f2$qaccver%>%unique()%>%length()%>%print()
  as.data.frame(table(blastn.Tab.f2$qaccver))%>%filter(Freq>1)%>%.[,2]%>%summary()%>%print()
  blastn.Tab.f2$qaccver%>%table()%>%sort(decreasing = T)%>%as.data.frame()%>%
    write.table(file=paste(blastnTab.rawFileName,".filteredDupFreq.txt",sep = ""),
                quote = F,sep = "\t",row.names = F,col.names = F)
  return(blastn.Tab.f2)
}
blastn_HEK3_200_1.f<-filtered.F(blastnTab.rawFileName="HEK3-200-1.blastnRes",
                                Amp.length=HEK3_200.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_200.I.len+10))
blastn_HEK3_200_2.f<-filtered.F(blastnTab.rawFileName="HEK3-200-2.blastnRes",
                                Amp.length=HEK3_200.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_200.I.len+10))
blastn_HEK3_200_3.f<-filtered.F(blastnTab.rawFileName="HEK3-200-3.blastnRes",
                                Amp.length=HEK3_200.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_200.I.len+10))
blastn_HEK3_1K_1.f<-filtered.F(blastnTab.rawFileName="HEK3-1K-1.blastnRes",
                               Amp.length=HEK3_1K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_1K.I.len+10))
blastn_HEK3_1K_2.f<-filtered.F(blastnTab.rawFileName="HEK3-1K-2.blastnRes",
                               Amp.length=HEK3_1K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_1K.I.len+10))
blastn_HEK3_1K_3.f<-filtered.F(blastnTab.rawFileName="HEK3-1K-3.blastnRes",
                               Amp.length=HEK3_1K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_1K.I.len+10))
# blastn_HEK3_8K_1.f<-filtered.F(blastnTab.rawFileName="HEK3-8K-1.blastnRes",
#                                Amp.length=HEK3_8K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_8K.I.len+10))
blastn_HEK3_8K_2.f<-filtered.F(blastnTab.rawFileName="HEK3-8K-2.blastnRes",
                               Amp.length=HEK3_8K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_8K.I.len+10))
blastn_HEK3_8K_2.f$qaccver%>%table()%>%sort(decreasing = T)%>%as.data.frame()%>%
  write.table(file="HEK3-8K-2.blastnRes.filteredDupFreq.txt",
              quote = F,sep = "\t",row.names = T,col.names = F)
# blastn_HEK3_8K_3.f<-filtered.F(blastnTab.rawFileName="HEK3-8K-3.blastnRes",
#                                Amp.length=HEK3_8K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_8K.I.len+10))

blastn_RUNX1_200_1.f<-filtered.F(blastnTab.rawFileName="RUNX1-200-1.blastnRes",
                                 Amp.length=RUNX1_200.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_200.I.len+10))
blastn_RUNX1_200_2.f<-filtered.F(blastnTab.rawFileName="RUNX1-200-2.blastnRes",
                                 Amp.length=RUNX1_200.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_200.I.len+10))
blastn_RUNX1_200_3.f<-filtered.F(blastnTab.rawFileName="RUNX1-200-3.blastnRes",
                                 Amp.length=RUNX1_200.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_200.I.len+10))
blastn_RUNX1_1K_1.f<-filtered.F(blastnTab.rawFileName="RUNX1-1K-1.blastnRes",
                                Amp.length=RUNX1_1K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_1K.I.len+10))
blastn_RUNX1_1K_2.f<-filtered.F(blastnTab.rawFileName="RUNX1-1K-2.blastnRes",
                                Amp.length=RUNX1_1K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_1K.I.len+10))
blastn_RUNX1_1K_3.f<-filtered.F(blastnTab.rawFileName="RUNX1-1K-3.blastnRes",
                                Amp.length=RUNX1_1K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_1K.I.len+10))
blastn_RUNX1_8K_1.f<-filtered.F(blastnTab.rawFileName="RUNX1-8K-1.blastnRes",
                                Amp.length=RUNX1_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_8K.I.len+10))
blastn_RUNX1_8K_2.f<-filtered.F(blastnTab.rawFileName="RUNX1-8K-2.blastnRes",
                                Amp.length=RUNX1_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_8K.I.len+10))
# blastn_RUNX1_8K_3.f<-filtered.F(blastnTab.rawFileName="RUNX1-8K-3.blastnRes", 
#                                 Amp.length=RUNX1_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_8K.I.len+10))

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
blastn_HEK4_8K_1.f$qaccver%>%table()%>%sort(decreasing = T)%>%as.data.frame()%>%
  write.table(file="HEK4-8224-1.blastnRes.filteredDupFreq.txt",
              quote = F,sep = "\t",row.names = T,col.names = F)
blastn_HEK4_8K_2.f<-filtered.F(blastnTab.rawFileName="HEK4-8224-2.blastnRes",
                               Amp.length=HEK4_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 60)
blastn_HEK4_8K_3.f<-filtered.F(blastnTab.rawFileName="HEK4-8224-3.blastnRes",
                               Amp.length=HEK4_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 60)

for (i in 1:23) {
  tt<-fread(file = list.files(pattern = "blastnRes.filteredDupFreq.txt")[i])
  print(list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[i])
  tt%>%pull(V2)%>%table()%>%print()
}
# for i in *.blastnRes.filteredDupFreq.txt; do awk -F"\t" '{print$1}' $i > ${i%%Freq.txt}Seq.txt; seqkit grep -j 20 -f ${i%%Freq.txt}Seq.txt ${i%%blastnRes.filteredDupFreq.txt}fa > ${i%%blastnRes.filteredDupFreq.txt}filteredDupSeq.fasta ; rm ${i%%Freq.txt}Seq.txt; done
library(Biostrings)
# HEK3_200 2-28A 33A 34A
# RUNX1_200 2-18A
# HEK4_200 2-24A
# HEK3_1K 2-8A
# RUNX1_1K 2-6A 9A
# HEK4_1K 2-5A
###############################################################################
HEK3_200_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[4])
for (i in c(2:26,28)) {
  HEK3_200_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[4]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_200_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[5])
for (i in c(2:26,28,33,43)) {
  HEK3_200_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[5]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_200_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[6])
for (i in c(2:23,26,27,35)) {
  HEK3_200_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[6]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_1K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[1])
for (i in c(2:6,8)) {
  HEK3_1K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[1]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_1K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[2])
for (i in c(2:8)) {
  HEK3_1K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[2]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_1K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[3])
for (i in c(2:7)) {
  HEK3_1K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[3]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK3_8K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[7])
for (i in c(2)) {
  HEK3_8K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[7]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_200_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[19])
for (i in c(2:13,15)) {
  RUNX1_200_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[19]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_200_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[20])
for (i in c(2:14,16,18)) {
  RUNX1_200_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[20]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_200_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[21])
for (i in c(2:12,14,15,17)) {
  RUNX1_200_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[21]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_1K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[16])
for (i in c(2:6)) {
  RUNX1_1K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[16]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_1K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[17])
for (i in c(2:6,9)) {
  RUNX1_1K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[17]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_1K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[18])
for (i in c(2:6)) {
  RUNX1_1K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[18]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_8K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[22])
for (i in c(2)) {
  RUNX1_8K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[22]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
RUNX1_8K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[23])
for (i in c(2)) {
  RUNX1_8K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=F)[23]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_200_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[11])
for (i in c(2:14,16,19,24)) {
  HEK4_200_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[11]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_200_2.DupFreq<-fread(file=list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[12])
for (i in c(2:15,17,19)) {
  HEK4_200_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[12]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_200_3.DupFreq<-fread(file=list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[13])
for (i in c(2:18,21)) {
  HEK4_200_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[13]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_1K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[8])
for (i in c(2:5)) {
  HEK4_1K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[8]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_1K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[9])
for (i in c(2:5)) {
  HEK4_1K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[9]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_1K_3.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[10])
for (i in c(2:4)) {
  HEK4_1K_3.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[10]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_8K_1.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[14])
for (i in c(2)) {
  HEK4_8K_1.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[14]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
HEK4_8K_2.DupFreq<-fread(file=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[15])
for (i in c(2)) {
  HEK4_8K_2.DupFreq%>%filter(V2==i)%>%select(V1)%>%
    write.table(file=sub(pattern="DupFreq.txt",
                         replacement=paste("_",i,"A.txt",sep=""),
                         x=list.files(pattern="blastnRes.filteredDupFreq.txt",recursive=T)[15]),
                quote = F,sep = "\t",row.names = F,col.names = F)
}

# for i in *.blastnRes.filtered_[2-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_[2-9]A.txt}fa > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_1[0-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_1[0-9]A.txt}fa > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_2[0-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_2[0-9]A.txt}fa > ${i%txt}fasta; done
# for i in *.blastnRes.filtered_3[0-9]A.txt; do echo $i; seqkit grep -f $i ${i%%blastnRes.filtered_3[0-9]A.txt}fa > ${i%txt}fasta; done
# cd ../
# for i in FullReference_*.fasta; do echo $i; /usr/bin/makeblastdb -in $i -dbtype 'nucl' -out ${i%%.fasta}; done
###blastn#####
# cd demultiplexed_V2/
# for i in HEK3-200-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_9A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_9A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_10A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_10A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_11A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_11A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_12A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_12A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_13A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_13A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-200-[1-3].blastnRes.filtered_14A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_200_14A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_1K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_1K_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_1K_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_1K_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_1K_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK3-8K-2.blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK3_8K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

# for i in RUNX1-200-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_200_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_200_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_200_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_200_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_200_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_200_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_200_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-1K-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_1K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-1K-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_1K_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-1K-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_1K_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-1K-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_1K_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in RUNX1-8K-[1-2].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_RUNX1_8K_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

# for i in HEK4-262-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_6A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_6A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_7A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_7A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_8A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_8A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_9A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_9A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_10A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_10A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_11A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_11A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_12A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_12A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_13A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_13A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-262-[1-3].blastnRes.filtered_14A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_262_14A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_1641_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_3A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_1641_3A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_4A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_1641_4A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_5A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_1641_5A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done
# for i in HEK4-8224-[1-2].blastnRes.filtered_2A.fasta; do echo $i; /usr/bin/blastn -query $i -task 'dc-megablast' -db ../FullReference_HEK4_8224_2A -out ${i%%fasta}blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80; done

for (i in c(1:142)) {
  rr<-fread(file = list.files(pattern = "A.blastnRes$",recursive = F)[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  print(list.files(pattern = "A.blastnRes$",recursive = T)[i])
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

# for i in HEK3-200-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK3-200-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK3-200-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK3-200-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK3-200-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK3-200-3.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK3-1K-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK3-1K-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK3-1K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK3-1K-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK3-1K-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK3-1K-3.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK3-8K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK3-8K-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done

# for i in RUNX1-200-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-200-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in RUNX1-200-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-200-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in RUNX1-200-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-200-3.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in RUNX1-1K-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-1K-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in RUNX1-1K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-1K-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in RUNX1-1K-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-1K-3.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in RUNX1-8K-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-8K-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in RUNX1-8K-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i RUNX1-8K-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done

# for i in HEK4-262-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-262-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-262-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-262-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-262-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-262-3.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-1641-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-1641-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-1641-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-1641-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-1641-3.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-1641-3.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-8224-1.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-8224-1.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# for i in HEK4-8224-2.*.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i HEK4-8224-2.fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done

