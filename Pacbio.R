setwd(dir = "/data/pan/Zhou_He/HZ-0129/")
# lima --ccs --min-score 80 --min-end-score 50 --min-ref-span 0.5 --different --min-scoring-regions 2 m84179_240228_131041_s4.hifi_reads.bc1001.bam barcode.fasta demultiplexed_v2.bam
# bam2fastq --split-barcodes -o AE demultiplexed_v2.bam
# rename demultiplexed fastq.gz files and move them into demultiplexed_v2 folder
# cd demultiplexed_V2/
# for i in *.fastq.gz;do seqkit fq2fa $i >${i%%stq.gz};done
# cd ../
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
#                                Amp.length=HEK3_8K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_8K.I.len+10)) # 0 passed read
blastn_HEK3_8K_2.f<-filtered.F(blastnTab.rawFileName="HEK3-8K-2.blastnRes",
                               Amp.length=HEK3_8K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_8K.I.len+10))
blastn_HEK3_8K_2.f$qaccver%>%table()%>%sort(decreasing = T)%>%as.data.frame()%>%
  write.table(file="HEK3-8K-2.blastnRes.filteredDupFreq.txt",
              quote = F,sep = "\t",row.names = T,col.names = F)
# blastn_HEK3_8K_3.f<-filtered.F(blastnTab.rawFileName="HEK3-8K-3.blastnRes",
#                                Amp.length=HEK3_8K.A.len,threshold.length=0.80,threshold.pident=80,distance.threshold=c(HEK3_8K.I.len+10)) # 0 passed read

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
#                                 Amp.length=RUNX1_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = c(RUNX1_8K.I.len+10)) # 0 passed read

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
# blastn_HEK4_8K_3.f<-filtered.F(blastnTab.rawFileName="HEK4-8224-3.blastnRes",
#                                Amp.length=HEK4_8K.A.len,threshold.length = 0.80,threshold.pident = 80,distance.threshold = 60) # 0 passed read

for (i in 1:23) {
  tt<-fread(file = list.files(pattern = "blastnRes.filteredDupFreq.txt")[i])
  print(list.files(pattern = "blastnRes.filteredDupFreq.txt",recursive = T)[i])
  tt%>%pull(V2)%>%table()%>%print()
}
# for i in *.blastnRes.filteredDupFreq.txt; do awk -F"\t" '{print$1}' $i > ${i%%Freq.txt}Seq.txt; seqkit grep -j 20 -f ${i%%Freq.txt}Seq.txt ${i%%blastnRes.filteredDupFreq.txt}fa > ${i%%blastnRes.filteredDupFreq.txt}filteredDupSeq.fasta ; rm ${i%%Freq.txt}Seq.txt; done
# construct different length full-length references (or copy previous full-length reference files)
# HEK3_200 2-28A 33A 34A
# RUNX1_200 2-18A
# HEK4_200 2-24A
# HEK3_1K 2-8A
# RUNX1_1K 2-6A 9A
# HEK4_1K 2-5A
#####################separate different "A" events################################
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
##### blastn against full-length reference #####
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
# filter reads based blastn results(full-length reference)
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
# filter WT
# /usr/bin/blastn -query HEK3-WT.fa -task 'dc-megablast' -db ../FullReference_HEK3_WT -out HEK3-WT.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80
# /usr/bin/blastn -query HEK4-WT.fa -task 'dc-megablast' -db ../FullReference_HEK4_WT -out HEK4-WT.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80
# /usr/bin/blastn -query RUNX1-WT.fa -task 'dc-megablast' -db ../FullReference_RUNX1_WT -out RUNX1-WT.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 80 -qcov_hsp_perc 80
for (i in c(1:3)) {
  rr<-fread(file = list.files(pattern = "WT.blastnRes$",recursive = F)[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  print(list.files(pattern = "WT.blastnRes$",recursive = T)[i])
  rr.f<-rr%>%mutate(D=slen-qlen)%>%mutate(D2=slen-length)%>%
    mutate(D3=qlen-length)%>%filter(abs(D2)<0.02*slen)%>%
    filter(abs(D3)<0.02*slen)%>%filter(mismatch<0.02*slen)
  rr.f%>%pull(qaccver)%>%unique()%>%length()%>%print()
  rr.f%>%pull(qaccver)%>%unique()%>%as.data.frame()%>%
    write.table(file=paste(list.files(pattern="WT.blastnRes$",recursive=F)[i],".filteredRead.txt",sep=""),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
# for i in *WT.blastnRes.filteredRead.txt; do seqkit grep -j 20 -f $i ${i%%blastnRes.filteredRead.txt}fa > ${i%%blastnRes.filteredRead.txt}filteredRead.fasta; done
# Last round filtration
# mkdir LastRoundFilteration; cd LastRoundFilteration
# create referenceFile: RUNX1_left1.fa RUNX1_left2.fa RUNX1_right1.fa HEK3_left.fa HEK3_right.fa HEK4_left.fa HEK4_right.fa
# for i in *.fa; do /usr/bin/makeblastdb -dbtype nucl -in $i; echo $i "done"; done
# for i in ../HEK3*.filteredRead.fasta; do m=${i##../}; echo $m; /usr/bin/blastn -query $i -task 'megablast' -db HEK3_left.fa -out ${m%%filteredRead.fasta}L.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in ../HEK3*.filteredRead.fasta; do m=${i##../}; echo $m; /usr/bin/blastn -query $i -task 'megablast' -db HEK3_right.fa -out ${m%%filteredRead.fasta}R.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query ../HEK3-WT.filteredRead.fasta -task 'megablast' -db HEK3_left.fa -out HEK3-WT.L.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query ../HEK3-WT.filteredRead.fasta -task 'megablast' -db HEK3_right.fa -out HEK3-WT.R.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# for i in ../HEK4*.filteredRead.fasta; do m=${i##../}; echo $m; /usr/bin/blastn -query $i -task 'megablast' -db HEK4_left.fa -out ${m%%filteredRead.fasta}L.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in ../HEK4*.filteredRead.fasta; do m=${i##../}; echo $m; /usr/bin/blastn -query $i -task 'megablast' -db HEK4_right.fa -out ${m%%filteredRead.fasta}R.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query ../HEK4-WT.filteredRead.fasta -task 'megablast' -db HEK4_left.fa -out HEK4-WT.L.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query ../HEK4-WT.filteredRead.fasta -task 'megablast' -db HEK4_right.fa -out HEK4-WT.R.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# for i in ../RUNX1*.filteredRead.fasta; do m=${i##../}; echo $m; /usr/bin/blastn -query $i -task 'megablast' -db RUNX1_left1.fa -out ${m%%filteredRead.fasta}L1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in ../RUNX1*.filteredRead.fasta; do m=${i##../}; echo $m; /usr/bin/blastn -query $i -task 'megablast' -db RUNX1_left2.fa -out ${m%%filteredRead.fasta}L2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in ../RUNX1*.filteredRead.fasta; do m=${i##../}; echo $m; /usr/bin/blastn -query $i -task 'megablast' -db RUNX1_right1.fa -out ${m%%filteredRead.fasta}R.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query ../RUNX1-WT.filteredRead.fasta -task 'megablast' -db RUNX1_left1.fa -out RUNX1-WT.L1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query ../RUNX1-WT.filteredRead.fasta -task 'megablast' -db RUNX1_left2.fa -out RUNX1-WT.L2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query ../RUNX1-WT.filteredRead.fasta -task 'megablast' -db RUNX1_right1.fa -out RUNX1-WT.R.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# rm HEK4-1641-1.blastnRes.filtered_5A.*
setwd("LastRoundFilteration/")
for (i in c(1:56)) {
  rrL<-fread(file = list.files(pattern = "HEK3.*.L.blastnRes$",recursive = F)[i])
  rrR<-fread(file = list.files(pattern = "HEK3.*.R.blastnRes$",recursive = F)[i])
  colnames(rrL)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  colnames(rrR)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rrL.f<-rrL%>%filter(length==slen)
  rrR.f<-rrR%>%filter(length==slen)
  print(i)
  print(list.files(pattern = "HEK3.*.L.blastnRes$",recursive = F)[i])
  rr.f<-intersect(rrL.f%>%pull(qaccver)%>%unique(),rrR.f%>%pull(qaccver)%>%unique())
  rr.f%>%length()%>%print()
  rr.f%>%as.data.frame()%>%
    write.table(file=paste(sub(pattern="L.blastnRes",replacement="",x=list.files(pattern = "HEK3.*.L.blastnRes$",recursive=F)[i]),
                           "filteredRead.txt",sep=""),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
for (i in c(1:6,8:52)) { # i=7 ;File 'HEK4-1641-2.blastnRes.filtered_5A.R.blastnRes' has size 0
  rrL<-fread(file = list.files(pattern = "HEK4.*.L.blastnRes$",recursive = F)[i])
  rrR<-fread(file = list.files(pattern = "HEK4.*.R.blastnRes$",recursive = F)[i])
  colnames(rrL)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  colnames(rrR)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rrL.f<-rrL%>%filter(length==slen)
  rrR.f<-rrR%>%filter(length==slen)
  print(i)
  print(list.files(pattern = "HEK4.*.L.blastnRes$",recursive = F)[i])
  rr.f<-intersect(rrL.f%>%pull(qaccver)%>%unique(),rrR.f%>%pull(qaccver)%>%unique())
  rr.f%>%length()%>%print()
  rr.f%>%as.data.frame()%>%
    write.table(file=paste(sub(pattern="L.blastnRes",replacement="",x=list.files(pattern = "HEK4.*.L.blastnRes$",recursive=F)[i]),
                           "filteredRead.txt",sep=""),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
for (i in c(1:36)) {
  rrL1<-fread(file = list.files(pattern = "RUNX1.*.L1.blastnRes$",recursive = F)[i])
  rrL2<-fread(file = list.files(pattern = "RUNX1.*.L2.blastnRes$",recursive = F)[i])
  rrR<-fread(file = list.files(pattern = "RUNX1.*.R.blastnRes$",recursive = F)[i])
  colnames(rrL1)<-c("qaccver","saccver","pident","length","mismatch",
                    "gapopen","qstart","qend","sstart","send",
                    "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  colnames(rrL2)<-c("qaccver","saccver","pident","length","mismatch",
                    "gapopen","qstart","qend","sstart","send",
                    "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  colnames(rrR)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rrL1.f<-rrL1%>%filter(length==slen)
  rrL2.f<-rrL2%>%filter(length==slen)
  rrR.f<-rrR%>%filter(length==slen)
  print(i)
  print(list.files(pattern = "RUNX1.*.R.blastnRes$",recursive = F)[i])
  rr.f<-intersect(intersect(rrL1.f%>%pull(qaccver)%>%unique(),
                            rrL2.f%>%pull(qaccver)%>%unique()),
                  rrR.f%>%pull(qaccver)%>%unique())
  rr.f%>%length()%>%print()
  rr.f%>%as.data.frame()%>%
    write.table(file=paste(sub(pattern="R.blastnRes",replacement="",x=list.files(pattern="RUNX1.*.R.blastnRes$",recursive=F)[i]),
                           "filteredRead.txt",sep=""),
                quote = F,sep = "\t",row.names = F,col.names = F)
}
# for i in *.blastnRes.filtered_[1-9]A.filteredRead.txt; do echo $i; seqkit grep -n -f $i ../${i%%blastnRes.filtered_[1-9]A.filteredRead.txt}fastq.gz > ${i%%txt}fq.gz; done
# for i in *.blastnRes.filtered_1[0-4]A.filteredRead.txt; do echo $i; seqkit grep -n -f $i ../${i%%blastnRes.filtered_1[0-4]A.filteredRead.txt}fastq.gz > ${i%%txt}fq.gz; done
# for i in *WT.filteredRead.txt;do echo $i; seqkit grep -n -f $i ../${i%%filteredRead.txt}fastq.gz > ${i%%txt}fq.gz; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_2A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_200_2A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_2A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_2A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_2A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_2A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_3A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_200_3A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_3A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_3A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_3A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_3A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_4A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_200_4A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_4A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_4A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_4A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_4A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_5A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_200_5A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_5A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_5A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_5A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_5A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_6A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_200_6A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_6A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_6A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_6A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_6A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_7A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_200_7A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_7A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_7A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_7A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_7A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-200-[1-3].blastnRes.filtered_8A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_200_8A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_8A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_8A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_8A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_8A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_9A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_9A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_9A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_9A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_10A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_10A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_10A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_10A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_11A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_11A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_11A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_11A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_12A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_12A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_12A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_12A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_13A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_13A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_13A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_13A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-200-[1-3].blastnRes.filtered_14A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_200_14A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-262-[1-3].blastnRes.filtered_14A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_262_14A.fasta $i ${i%%fq.gz}bam; done

# for i in RUNX1-1K-[1-3].blastnRes.filtered_2A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_1K_2A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_2A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_1K_2A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_2A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_1641_2A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-1K-[1-3].blastnRes.filtered_3A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_1K_3A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_3A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_1K_3A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_3A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_1641_3A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-1K-[1-3].blastnRes.filtered_4A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_1K_4A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_4A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_1K_4A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK4-1641-[1-3].blastnRes.filtered_4A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_1641_4A.fasta $i ${i%%fq.gz}bam; done
# for i in RUNX1-1K-[1-3].blastnRes.filtered_5A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_1K_5A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_5A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_1K_5A.fasta $i ${i%%fq.gz}bam; done
# for i in HEK3-1K-[1-3].blastnRes.filtered_6A.filteredRead.fq.gz; do pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_1K_6A.fasta $i ${i%%fq.gz}bam; done
# pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK3_WT.fasta HEK3-WT.filteredRead.fq.gz HEK3-WT.filteredRead.bam
# pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_HEK4_WT.fasta HEK4-WT.filteredRead.fq.gz HEK4-WT.filteredRead.bam
# pbmm2 align --preset CCS --sort --log-level INFO --bam-index BAI -N 1 ../../FullReference_RUNX1_WT.fasta RUNX1-WT.filteredRead.fq.gz RUNX1-WT.filteredRead.bam
## calculate base frequencies
## library(Rsamtools)
## library(tidyverse)
## for (i in 1:138) { #19 20 32 33 45 46 
##   print(c(i,list.files(pattern="bam$")[i]))
##   pTab<-pileup(file=list.files(pattern="bam$")[i],
##                pileupParam=PileupParam(max_depth=200000,include_deletions=TRUE, include_insertions=TRUE,distinguish_strands=F,min_base_quality=0))
##   pTab2<-data.frame(pos=rep(1:max(pTab$pos),each=8),nucleotide=rep(levels(pTab$nucleotide),times=max(pTab$pos)))
##   pTab3<-full_join(x=pTab2,y = pTab[-1])
##   pTab3[is.na(pTab3)]<-0
##   pTab3<-pTab3%>%group_by(pos)%>%mutate(Freq=100*count/sum(count))
##   write.table(x=pTab3,file=sub(pattern="filteredRead.bam",replacement="nucleotideFreq.txt",
##                                x=sub(pattern=".blastnRes.filtered",replacement="",
##                                      x=list.files(pattern="bam$")[i])),
##               quote=F,sep="\t",row.names=F,col.names=T)
## }
## for (i in 1:138) {
##   ff<-fread(file = list.files(pattern = "Freq.txt")[i])
##   ff.p<-ff%>%group_by(pos)%>%slice_max(order_by = Freq)%>%
##     filter(Freq<99.9)%>%pull(pos)%>%unique()
##   ff.f<-ff%>%filter(pos%in%ff.p)%>%filter(count!=0)
##   write.table(x=ff.f,file=sub(pattern="txt",replacement="f1.txt",x=list.files(pattern = "Freq.txt")[i]),
##               quote=F,sep="\t",row.names=F,col.names=T)
## }

# for i in *.filteredRead.fq.gz; do seqkit fq2fa -j 20 $i -o ${i%%fq.gz}fa; done
#####focus on amplification area#####
# for i in RUNX1-200-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db ../../duplicationArea.RUNX1-200.fasta -out ${i%%fa}dup.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query RUNX1-WT.filteredRead.fa -task 'megablast' -db ../../duplicationArea.RUNX1-200.fasta -out RUNX1-WT.filteredRead.dup200.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
for (i in 1:7) {
  print(list.files(pattern="RUNX1-200-1.*filteredRead.dup.blastnRes")[i])
  m=i+1
  rr<-fread(file=list.files(pattern="RUNX1-200-1.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="RUNX1-200-1.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:7) {
  print(list.files(pattern="RUNX1-200-2.*filteredRead.dup.blastnRes")[i])
  m=i+1
  rr<-fread(file=list.files(pattern="RUNX1-200-2.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="RUNX1-200-2.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:7) {
  print(list.files(pattern="RUNX1-200-3.*filteredRead.dup.blastnRes")[i])
  m=i+1
  rr<-fread(file=list.files(pattern="RUNX1-200-3.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="RUNX1-200-3.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
rr<-fread(file=list.files(pattern="RUNX1-WT.*filteredRead.dup200.blastnRes"))
colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                "gapopen","qstart","qend","sstart","send",
                "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.f<-rr%>%filter(length==slen)%>%select(qaccver)
print(c(nrow(rr.f)))
write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                            x=list.files(pattern="RUNX1-WT.*filteredRead.dup200.blastnRes")),
            quote=F,sep="\t",row.names=F,col.names=F)
# vi duplicationArea.RUNX1-1K.v1.fasta
# vi duplicationArea.RUNX1-1K.v2.fasta
# for i in duplicationArea.RUNX1-1K.v*; do /usr/bin/makeblastdb -in $i -dbtype 'nucl' -out ${i%%.fasta}; done
# for i in RUNX1-1K-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.RUNX1-1K.v1 -out ${i%%fa}dup.v1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in RUNX1-1K-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.RUNX1-1K.v2 -out ${i%%fa}dup.v2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query RUNX1-WT.filteredRead.fa -task 'megablast' -db duplicationArea.RUNX1-1K.v1 -out RUNX1-WT.filteredRead.dup1K.v1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query RUNX1-WT.filteredRead.fa -task 'megablast' -db duplicationArea.RUNX1-1K.v2 -out RUNX1-WT.filteredRead.dup1K.v2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
for (i in 1:4) {
  print(list.files(pattern="RUNX1-1K-1.*filteredRead.dup.v1.blastnRes")[i])
  m=i+1
  rr.v1<-fread(file=list.files(pattern="RUNX1-1K-1.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="RUNX1-1K-1.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                              x=list.files(pattern="RUNX1-1K-1.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:4) {
  print(list.files(pattern="RUNX1-1K-2.*filteredRead.dup.v1.blastnRes")[i])
  m=i+1
  rr.v1<-fread(file=list.files(pattern="RUNX1-1K-2.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="RUNX1-1K-2.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="RUNX1-1K-2.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:4) {
  print(list.files(pattern="RUNX1-1K-3.*filteredRead.dup.v1.blastnRes")[i])
  m=i+1
  rr.v1<-fread(file=list.files(pattern="RUNX1-1K-3.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="RUNX1-1K-3.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="RUNX1-1K-3.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
rr.v1<-fread(file=list.files(pattern="RUNX1-WT.*filteredRead.dup1K.v1.blastnRes"))
colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                "gapopen","qstart","qend","sstart","send",
                "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)
rr.v2<-fread(file=list.files(pattern="RUNX1-WT.*filteredRead.dup1K.v2.blastnRes"))
colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)
print(length(union(rr.v1.f,rr.v2.f)))
write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                            x=list.files(pattern="RUNX1-WT.*filteredRead.dup1K.v1.blastnRes")),
            quote=F,sep="\t",row.names=F,col.names=F)

# vi duplicationArea.HEK3-200.v1.fasta
# vi duplicationArea.HEK3-200.v2.fasta
# for i in duplicationArea.HEK3-200.v*; do /usr/bin/makeblastdb -in $i -dbtype 'nucl' -out ${i%%.fasta}; done
# for i in HEK3-200-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.HEK3-200.v1 -out ${i%%fa}dup.v1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in HEK3-200-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.HEK3-200.v2 -out ${i%%fa}dup.v2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query HEK3-WT.filteredRead.fa -task 'megablast' -db duplicationArea.HEK3-200.v1 -out HEK3-WT.filteredRead.dup200.v1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query HEK3-WT.filteredRead.fa -task 'megablast' -db duplicationArea.HEK3-200.v2 -out HEK3-WT.filteredRead.dup200.v2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
for (i in 6:13) {
  print(list.files(pattern="HEK3-200-1.*filteredRead.dup.v1.blastnRes")[i])
  m=(i-5)+1
  rr.v1<-fread(file=list.files(pattern="HEK3-200-1.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-200-1.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-200-1.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:5) {
  print(list.files(pattern="HEK3-200-1.*filteredRead.dup.v1.blastnRes")[i])
  m=(i+8)+1
  rr.v1<-fread(file=list.files(pattern="HEK3-200-1.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-200-1.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-200-1.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 6:13) {
  print(list.files(pattern="HEK3-200-2.*filteredRead.dup.v1.blastnRes")[i])
  m=(i-5)+1
  rr.v1<-fread(file=list.files(pattern="HEK3-200-2.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-200-2.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-200-2.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:5) {
  print(list.files(pattern="HEK3-200-2.*filteredRead.dup.v1.blastnRes")[i])
  m=(i+8)+1
  rr.v1<-fread(file=list.files(pattern="HEK3-200-2.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-200-2.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-200-2.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 6:13) {
  print(list.files(pattern="HEK3-200-3.*filteredRead.dup.v1.blastnRes")[i])
  m=(i-5)+1
  rr.v1<-fread(file=list.files(pattern="HEK3-200-3.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-200-3.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-200-3.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:5) {
  print(list.files(pattern="HEK3-200-3.*filteredRead.dup.v1.blastnRes")[i])
  m=(i+8)+1
  rr.v1<-fread(file=list.files(pattern="HEK3-200-3.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-200-3.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-200-3.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
rr.v1<-fread(file=list.files(pattern="HEK3-WT.*filteredRead.dup200.v1.blastnRes"))
colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)
rr.v2<-fread(file=list.files(pattern="HEK3-WT.*filteredRead.dup200.v2.blastnRes"))
colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)
print(length(union(rr.v1.f,rr.v2.f)))
write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                             x=list.files(pattern="HEK3-WT.*filteredRead.dup200.v1.blastnRes")),
            quote=F,sep="\t",row.names=F,col.names=F)
# vi duplicationArea.HEK3-1K.v1.fasta
# vi duplicationArea.HEK3-1K.v2.fasta
# for i in duplicationArea.HEK3-1K.v*; do /usr/bin/makeblastdb -in $i -dbtype 'nucl' -out ${i%%.fasta}; done
# for i in HEK3-1K-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.HEK3-1K.v1 -out ${i%%fa}dup.v1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in HEK3-1K-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.HEK3-1K.v2 -out ${i%%fa}dup.v2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query HEK3-WT.filteredRead.fa -task 'megablast' -db duplicationArea.HEK3-1K.v1 -out HEK3-WT.filteredRead.dup1K.v1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query HEK3-WT.filteredRead.fa -task 'megablast' -db duplicationArea.HEK3-1K.v2 -out HEK3-WT.filteredRead.dup1K.v2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
for (i in 1:5) { 
  print(list.files(pattern="HEK3-1K-1.*filteredRead.dup.v1.blastnRes")[i])
  m=i+1
  rr.v1<-fread(file=list.files(pattern="HEK3-1K-1.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-1K-1.*filteredRead.dup.v2.blastnRes")[i]) #i=5 6A no passed read
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-1K-1.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:5) { 
  print(list.files(pattern="HEK3-1K-2.*filteredRead.dup.v1.blastnRes")[i])
  m=i+1
  rr.v1<-fread(file=list.files(pattern="HEK3-1K-2.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-1K-2.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-1K-2.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:5) { 
  print(list.files(pattern="HEK3-1K-3.*filteredRead.dup.v1.blastnRes")[i])
  m=i+1
  rr.v1<-fread(file=list.files(pattern="HEK3-1K-3.*filteredRead.dup.v1.blastnRes")[i])
  colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.v2<-fread(file=list.files(pattern="HEK3-1K-3.*filteredRead.dup.v2.blastnRes")[i])
  colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(union(rr.v1.f,rr.v2.f))))
  write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                               x=list.files(pattern="HEK3-1K-3.*filteredRead.dup.v1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
rr.v1<-fread(file=list.files(pattern="HEK3-WT.*filteredRead.dup1K.v1.blastnRes"))
colnames(rr.v1)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.v1.f<-rr.v1%>%filter(length==slen)%>%pull(qaccver)
rr.v2<-fread(file=list.files(pattern="HEK3-WT.*filteredRead.dup1K.v2.blastnRes"))
colnames(rr.v2)<-c("qaccver","saccver","pident","length","mismatch",
                   "gapopen","qstart","qend","sstart","send",
                   "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.v2.f<-rr.v2%>%filter(length==slen)%>%pull(qaccver)
print(length(union(rr.v1.f,rr.v2.f)))
write.table(x=as.data.frame(union(rr.v1.f,rr.v2.f)),file=sub(pattern=".v1.blastnRes$",replacement="F.txt",
                                                             x=list.files(pattern="HEK3-WT.*filteredRead.dup1K.v1.blastnRes")),
            quote=F,sep="\t",row.names=F,col.names=F)

# for i in HEK4-262-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db ../../duplicationArea.HEK4-262.fasta -out ${i%%fa}dup.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query HEK4-WT.filteredRead.fa -task 'megablast' -db ../../duplicationArea.HEK4-262.fasta -out HEK4-WT.filteredRead.dup262.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
list.files(pattern="HEK4-262-1.*filteredRead.dup.blastnRes")[1:5]
for (i in 6:13) {
  print(list.files(pattern="HEK4-262-1.*filteredRead.dup.blastnRes")[i])
  m=(i-5)+1
  rr<-fread(file=list.files(pattern="HEK4-262-1.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="HEK4-262-1.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=T)
}
for (i in 1:5) {
  print(list.files(pattern="HEK4-262-1.*filteredRead.dup.blastnRes")[i])
  m=(i+8)+1
  rr<-fread(file=list.files(pattern="HEK4-262-1.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="HEK4-262-1.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=T)
}
for (i in 6:13) {
  print(list.files(pattern="HEK4-262-2.*filteredRead.dup.blastnRes")[i])
  m=(i-5)+1
  rr<-fread(file=list.files(pattern="HEK4-262-2.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="HEK4-262-2.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=T)
}
for (i in 1:5) {
  print(list.files(pattern="HEK4-262-2.*filteredRead.dup.blastnRes")[i])
  m=(i+8)+1
  rr<-fread(file=list.files(pattern="HEK4-262-2.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="HEK4-262-2.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=T)
}
for (i in 6:13) {
  print(list.files(pattern="HEK4-262-3.*filteredRead.dup.blastnRes")[i])
  m=(i-5)+1
  rr<-fread(file=list.files(pattern="HEK4-262-3.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="HEK4-262-3.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=T)
}
for (i in 1:5) {
  print(list.files(pattern="HEK4-262-3.*filteredRead.dup.blastnRes")[i])
  m=(i+8)+1
  rr<-fread(file=list.files(pattern="HEK4-262-3.*filteredRead.dup.blastnRes")[i])
  colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                  "gapopen","qstart","qend","sstart","send",
                  "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.f<-rr%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%select(1)
  print(c(m,"A",nrow(rr.f)))
  write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                              x=list.files(pattern="HEK4-262-3.*filteredRead.dup.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=T)
}
rr<-fread(file=list.files(pattern="HEK4-WT.*filteredRead.dup262.blastnRes"))
colnames(rr)<-c("qaccver","saccver","pident","length","mismatch",
                "gapopen","qstart","qend","sstart","send",
                "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.f<-rr%>%filter(length==slen)%>%select(qaccver)
print(c(nrow(rr.f)))
write.table(x=rr.f,file=sub(pattern=".blastnRes$",replacement="F.txt",
                            x=list.files(pattern="HEK4-WT.*filteredRead.dup262.blastnRes")),
            quote=F,sep="\t",row.names=F,col.names=T)
# vi duplicationArea.HEK4-1K.fg1.fasta
# vi duplicationArea.HEK4-1K.fg2.fasta
# vi duplicationArea.HEK4-1K.fg3.fasta
# for i in duplicationArea.HEK4-1K.fg*; do /usr/bin/makeblastdb -in $i -dbtype 'nucl' -out ${i%%.fasta}; done
# for i in HEK4-1641-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.HEK4-1K.fg1 -out ${i%%fa}dup.fg1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in HEK4-1641-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.HEK4-1K.fg2 -out ${i%%fa}dup.fg2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# for i in HEK4-1641-*.fa; do echo $i; /usr/bin/blastn -query $i -task 'megablast' -db duplicationArea.HEK4-1K.fg3 -out ${i%%fa}dup.fg3.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100; done
# /usr/bin/blastn -query HEK4-WT.filteredRead.fa -task 'megablast' -db duplicationArea.HEK4-1K.fg1 -out HEK4-WT.filteredRead.dup1641.fg1.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query HEK4-WT.filteredRead.fa -task 'megablast' -db duplicationArea.HEK4-1K.fg2 -out HEK4-WT.filteredRead.dup1641.fg2.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
# /usr/bin/blastn -query HEK4-WT.filteredRead.fa -task 'megablast' -db duplicationArea.HEK4-1K.fg3 -out HEK4-WT.filteredRead.dup1641.fg3.blastnRes -outfmt '6 std qlen slen qcovhsp qcovus' -perc_identity 100
for (i in 1:3) { 
  print(list.files(pattern="HEK4-1641-1.*filteredRead.dup.fg1.blastnRes")[i])
  m=i+1
  rr.fg1<-fread(file=list.files(pattern="HEK4-1641-1.*filteredRead.dup.fg1.blastnRes")[i])
  colnames(rr.fg1)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg1.f<-rr.fg1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.fg2<-fread(file=list.files(pattern="HEK4-1641-1.*filteredRead.dup.fg2.blastnRes")[i])
  colnames(rr.fg2)<-c("qaccver","saccver","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send",
                     "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg2.f<-rr.fg2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.fg3<-fread(file=list.files(pattern="HEK4-1641-1.*filteredRead.dup.fg3.blastnRes")[i])
  colnames(rr.fg3)<-c("qaccver","saccver","pident","length","mismatch",
                      "gapopen","qstart","qend","sstart","send",
                      "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg3.f<-rr.fg3%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f))))
  write.table(x=as.data.frame(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f)),
              file=sub(pattern=".fg1.blastnRes$",replacement="F.txt",
                       x=list.files(pattern="HEK4-1641-1.*filteredRead.dup.fg1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:3) { 
  print(list.files(pattern="HEK4-1641-2.*filteredRead.dup.fg1.blastnRes")[i])
  m=i+1
  rr.fg1<-fread(file=list.files(pattern="HEK4-1641-2.*filteredRead.dup.fg1.blastnRes")[i])
  colnames(rr.fg1)<-c("qaccver","saccver","pident","length","mismatch",
                      "gapopen","qstart","qend","sstart","send",
                      "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg1.f<-rr.fg1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.fg2<-fread(file=list.files(pattern="HEK4-1641-2.*filteredRead.dup.fg2.blastnRes")[i])
  colnames(rr.fg2)<-c("qaccver","saccver","pident","length","mismatch",
                      "gapopen","qstart","qend","sstart","send",
                      "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg2.f<-rr.fg2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.fg3<-fread(file=list.files(pattern="HEK4-1641-2.*filteredRead.dup.fg3.blastnRes")[i])
  colnames(rr.fg3)<-c("qaccver","saccver","pident","length","mismatch",
                      "gapopen","qstart","qend","sstart","send",
                      "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg3.f<-rr.fg3%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f))))
  write.table(x=as.data.frame(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f)),
              file=sub(pattern=".fg1.blastnRes$",replacement="F.txt",
                       x=list.files(pattern="HEK4-1641-2.*filteredRead.dup.fg1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
for (i in 1:3) { 
  print(list.files(pattern="HEK4-1641-3.*filteredRead.dup.fg1.blastnRes")[i])
  m=i+1
  rr.fg1<-fread(file=list.files(pattern="HEK4-1641-3.*filteredRead.dup.fg1.blastnRes")[i])
  colnames(rr.fg1)<-c("qaccver","saccver","pident","length","mismatch",
                      "gapopen","qstart","qend","sstart","send",
                      "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg1.f<-rr.fg1%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.fg2<-fread(file=list.files(pattern="HEK4-1641-3.*filteredRead.dup.fg2.blastnRes")[i])
  colnames(rr.fg2)<-c("qaccver","saccver","pident","length","mismatch",
                      "gapopen","qstart","qend","sstart","send",
                      "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg2.f<-rr.fg2%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  rr.fg3<-fread(file=list.files(pattern="HEK4-1641-3.*filteredRead.dup.fg3.blastnRes")[i])
  colnames(rr.fg3)<-c("qaccver","saccver","pident","length","mismatch",
                      "gapopen","qstart","qend","sstart","send",
                      "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
  rr.fg3.f<-rr.fg3%>%filter(length==slen)%>%pull(qaccver)%>%table()%>%as.data.frame()%>%
    filter(Freq==m)%>%pull(1)
  print(c(m,"A",length(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f))))
  write.table(x=as.data.frame(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f)),
              file=sub(pattern=".fg1.blastnRes$",replacement="F.txt",
                       x=list.files(pattern="HEK4-1641-3.*filteredRead.dup.fg1.blastnRes")[i]),
              quote=F,sep="\t",row.names=F,col.names=F)
}
rr.fg1<-fread(file=list.files(pattern="HEK4-WT.*filteredRead.dup1641.fg1.blastnRes"))
colnames(rr.fg1)<-c("qaccver","saccver","pident","length","mismatch",
                    "gapopen","qstart","qend","sstart","send",
                    "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.fg1.f<-rr.fg1%>%filter(length==slen)%>%pull(qaccver)
rr.fg2<-fread(file=list.files(pattern="HEK4-WT.*filteredRead.dup1641.fg2.blastnRes"))
colnames(rr.fg2)<-c("qaccver","saccver","pident","length","mismatch",
                    "gapopen","qstart","qend","sstart","send",
                    "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.fg2.f<-rr.fg2%>%filter(length==slen)%>%pull(qaccver)
rr.fg3<-fread(file=list.files(pattern="HEK4-WT.*filteredRead.dup1641.fg3.blastnRes"))
colnames(rr.fg3)<-c("qaccver","saccver","pident","length","mismatch",
                    "gapopen","qstart","qend","sstart","send",
                    "evalue","bitscore","qlen","slen","qcovhsp","qcovus")
rr.fg3.f<-rr.fg3%>%filter(length==slen)%>%pull(qaccver)
print(c(length(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f))))
write.table(x=as.data.frame(intersect(intersect(rr.fg1.f,rr.fg2.f),rr.fg3.f)),
            file=sub(pattern=".fg1.blastnRes$",replacement="F.txt",
                     x=list.files(pattern="HEK4-WT.*filteredRead.dup1641.fg1.blastnRes")),
            quote=F,sep="\t",row.names=F,col.names=F)
