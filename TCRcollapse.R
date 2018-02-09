# d160316<-read.table("/Users/yingy_adm/Documents/imgtout/tcr160316.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
# d160316$batch<-rep("1",dim(d160316)[1])
# d220616<-read.table("/Users/yingy_adm/Documents/imgtout/tcr220616.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
# d220616$batch<-rep("2",dim(d220616)[1])
# d1116 <- read.table("/Users/yingy_adm/Documents/imgtout/tcr2016nov.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
# d1116$batch<-rep("3",dim(d1116)[1])
# d170120<-read.table("/Users/yingy_adm/Documents/TCRdata/170120/tcr2017.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
# d170120$batch<-rep("4",dim(d170120)[1])
#tcrdata <- rbind(t1,t2,t3)
#tcrdata <- tcrdata[,-2]
#rm("t1","t2","t3")
d160316<-read.table("/Users/yingy_adm/Documents/imgtout/r160316/pd160316.txt",sep="\t",stringsAsFactors=FALSE)
d160316$batch<-rep("1",dim(d160316)[1])
d220616<-read.table("/Users/yingy_adm/Documents/imgtout/r220616/pd220616.txt",sep="\t",stringsAsFactors=FALSE)
d220616$batch<-rep("2",dim(d220616)[1])
d1116 <- read.table("/Users/yingy_adm/Documents/imgtout/r161027/pd161027.txt",sep="\t",stringsAsFactors=FALSE)
d1116$batch<-rep("3",dim(d1116)[1])
d170120<-read.table("/Users/yingy_adm/Documents/imgtout/r231216/pd231216.txt",sep="\t",stringsAsFactors=FALSE)
d170120$batch<-rep("4",dim(d170120)[1])
tcrdata <- rbind(d160316,d220616,d1116,d170120)
#substitute "_" with "|" in Sequence.ID:
#tcrdata$Sequence.ID=gsub("_","|",tcrdata$Sequence.ID)

#extract dupcount from seq ID################################
library(stringr)
seqID=str_split_fixed(tcrdata$Sequence.ID, "[|]", 4)#num,barcode,umi,dupcount
#apply(seqID,2,function(x) sum(is.na(x)))#check NA
barcode=seqID[,2]
dupc=str_split_fixed(seqID[,4],"=",2)[,2]#convert dupcount to numric and add it to dataset
tcrdata$dupcount<- as.numeric(dupc)
tcrdata$V.GENE=str_split_fixed(tcrdata$V.GENE.and.allele,'[*]',2)[,1]#new v-gene column
tcrdata$J.GENE=str_split_fixed(tcrdata$J.GENE.and.allele,'[*]',2)[,1]#new J-gene column
#tcrdata$D.GENE=str_split_fixed(tcrdata$D.GENE.and.allele,'[*]',2)[,1]#new D-gene column
tcrdata$R1=str_split_fixed(str_split_fixed(barcode, ",", 3)[,1],"=",2)[,2]
tcrdata$R2=str_split_fixed(barcode, ",", 3)[,2]
tcrdata$rb=str_split_fixed(barcode, ",", 4)[,3]
tcrdata$umi=str_split_fixed(seqID[,3],"=",2)[,2]
rm(seqID,barcode,dupc)
tcrdata=tcrdata[,c(-1,-2,-3,-4)]#remove the old Seq.ID,Functionality,V,J column
colnames(tcrdata)<-c("CDR3nt","CDR3a","batch","dupcount","V.GENE","J.GENE","R1barcode","R2barcode","RowBarcode","UMI")

#######################################
#correct the sequences such that allow one nt mismatch if records with the same flag
tcrdata$flag<-paste(tcrdata$batch,tcrdata$R1barcode,tcrdata$R2barcode,tcrdata$RowBarcode,tcrdata$UMI,tcrdata$V.GENE,tcrdata$J.GENE,sep="|")
length(table(tcrdata$flag))#299,715
flagtab<-as.data.frame(table(tcrdata$flag))
flagtab<-subset(flagtab,Freq>1)# select 105,468 of 299,715
n=0 #i="3|B2|2|G|AAGTGG|Homsap TRBV15|Homsap TRBJ2-3"(91459); j=1273931
for (i in flagtab$Var1[30932:105468])
{ index<-which(tcrdata$flag==i)
maxdup<-max(tcrdata[index,]$dupcount)
if (maxdup>1)
{ #function to correct sequences with consensus
  ind_maxdup<-which.max(tcrdata[index,]$dupcount)
  ref<-tcrdata[index,]$CDR3nt[ind_maxdup]
  index <- index[-ind_maxdup]
  for (j in index)
  {
    if (adist(tcrdata[j,1],ref)==1)
    {
      tcrdata[j,1]<-ref #change the nt seq
      tcrdata[j,2]<-tcrdata[index,]$CDR3a[ind_maxdup] #change the amino acid seq
      n=n+1}
  }
}
}
write.table(tcrdata, file = "/Users/yingy_adm/Documents/TCRdata/corrected_tcrdata.csv", sep="\t", row.names = T, qmethod = "double")
tcrdata<-read.table(file = "/Users/yingy_adm/Documents/TCRdata/corrected_tcrdata.csv", sep="\t")
tcrdata$flag<-paste(tcrdata$flag,tcrdata$CDR3nt,tcrdata$CDR3a,sep = "|")
corr_tcrdata<-aggregate(dupcount~flag,data=tcrdata,FUN=sum)
# a5_data$R1barcode=str_split_fixed(a5_data$flag5, "[|]", 5)[,1]
# a5_data$R2barcode=str_split_fixed(a5_data$flag5, "[|]", 5)[,2]

#write.table(a4_data, file = "/Users/yingy_adm/Documents/TCRdata/newtcr4.csv", sep="\t", row.names = T, qmethod = "double")
write.table(a5_data, file = "/Users/yingy_adm/Documents/TCRdata/newtcr5.csv", sep="\t", row.names = T, qmethod = "double")

#tcrdata$number=as.character(seqID[,1])#read's idx
#tcrdata=tcrdata[,c(-1)]#remove the old Sequence.number
#tcrdata$seq.ID=paste(seqID[,2],seqID[,3], sep="|")#new seqID with barcode,umi

#remove duplicated reads with identical "seq.ID","batch","V.GENE","J.GENE","CDR3"(collapse UMI)
library(plyr)
#remove UMI in flag
flag=str_split_fixed(corr_tcrdata$flag,"[|]", 9)
corr_tcrdata$flag<-paste(flag[,1],flag[,2],flag[,3],flag[,4],flag[,6],flag[,7],flag[,8],flag[,9],sep='|') 
sdata<-subset(corr_tcrdata,dupcount>1)#sah:selected dataset removing reads of dupcount=1 ???
sdata$molcount<-rep(1,dim(sdata)[1])#initiate molecular count as 1
#tcrdata$flag5<-paste(tcrdata$R1barcode,tcrdata$R2barcode,tcrdata$V.GENE,tcrdata$J.GENE,tcrdata$CDR3,sep='|') # 5 conditions
#new=ddply(tcrdata,.(J.GENE,CDR3,seq.ID,V.GENE),nrow,.drop = TRUE)#collapsed data with count
#tcrdata %>% distinct(flag, .keep_all = TRUE)#collapsed data
#tcrdata %>% group_by(flag) %>% summarise(dupcount = sum(dupcount))#flag with num.count
#ddply(tcrdata,"flag",numcolwise(sum))
#a4_data<-aggregate(dupcount~flag4,data=tcrdata,FUN=sum)

#reduce to original condition and get the count of umi types (vivo proliferation)
ah_data<-aggregate(molcount~flag,data=sdata,FUN=sum)
ah_data$batch=str_split_fixed(ah_data$flag, "[|]", 8)[,1]
ah_data$R1barcode=str_split_fixed(ah_data$flag, "[|]", 8)[,2]
ah_data$R2barcode=str_split_fixed(ah_data$flag, "[|]", 8)[,3]
ah_data$RowBarcode=str_split_fixed(ah_data$flag, "[|]", 8)[,4]
ah_data$V.GENE=str_split_fixed(ah_data$flag, "[|]", 8)[,5]
ah_data$J.GENE=str_split_fixed(ah_data$flag, "[|]", 8)[,6]
ah_data$CDR3nt=str_split_fixed(ah_data$flag, "[|]", 8)[,7]
ah_data$CDR3a=str_split_fixed(ah_data$flag, "[|]", 8)[,8]
ah_data <- ah_data[,-1]
#patient id batch 1
ii1357=which(ah_data$batch=="1" & (ah_data$R2barcode %in% c("1","2","3") & ah_data$R1barcode %in% c("A1","B1"))|
               (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A2","B2")))
ii1358=which(ah_data$batch=="1" & (ah_data$R2barcode %in% c("2","3") & ah_data$R1barcode %in% c("A2","B2"))|
               (ah_data$R2barcode %in% c("1","2") & ah_data$R1barcode %in% c("A3","B3"))|
               (ah_data$R2barcode %in% c("2","3") & ah_data$R1barcode %in% c("A4","B4")))
ii1364=which(ah_data$batch=="1" & (ah_data$R2barcode %in% c("1","2","3"))&ah_data$R1barcode %in% c("A5","B5")|
               (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A6","B6")))
ii1368=which(ah_data$batch=="1" & (ah_data$R2barcode %in% c("2","3") & ah_data$R1barcode %in% c("A6","B6"))|
               (ah_data$R2barcode %in% c("1","2") & ah_data$R1barcode %in% c("A7","B7")))
ii1370=which(ah_data$batch=="1" & (ah_data$R2barcode =="3" & ah_data$R1barcode %in% c("A7","B7"))|
               (ah_data$R2barcode %in% c("1","2","3") & ah_data$R1barcode %in% c("A8","B8")))
ii1386=which(ah_data$batch=="1" & (ah_data$R2barcode %in% c("1","2","3","4") & ah_data$R1barcode %in% c("A10","B10"))|
               (ah_data$R2barcode %in% c("4","5","6") & ah_data$R1barcode %in% c("A9","B9")))
ii1390=which(ah_data$batch=="1" & (ah_data$R2barcode %in% c("5","6","7","8") & ah_data$R1barcode %in% c("A10","B10"))|
               (ah_data$R2barcode %in% c("7","8","9") & ah_data$R1barcode %in% c("A9","B9")))
ii1393=which(ah_data$batch=="1" & (ah_data$R2barcode %in% c("9","10","11","12") & ah_data$R1barcode %in% c("A10","B10"))|
               (ah_data$R2barcode %in% c("10","11","12") & ah_data$R1barcode %in% c("A9","B9")))
ii1341=which(ah_data$batch=="1" & ah_data$R2barcode =="1" & ah_data$R1barcode %in% c("A9","B9"))
ii1347=which(ah_data$batch=="1" & ah_data$R2barcode =="2" & ah_data$R1barcode %in% c("A9","B9"))
ii1331=which(ah_data$batch=="1" & ah_data$R2barcode =="3" & ah_data$R1barcode %in% c("A9","B9"))
ii1335_7=which(ah_data$batch=="1" & ah_data$R2barcode =="3" & ah_data$R1barcode %in% c("A3","B3"))#TCC1335.7 
ii1005_2_58=which(ah_data$batch=="1" & ah_data$R2barcode =="1" & ah_data$R1barcode %in% c("A4","B4"))#TCC1005.2.58
#patients id for batch 2
ahi1357=which(ah_data$batch=="2" & (ah_data$R2barcode %in% c("5","6","7","8"))&ah_data$R1barcode %in% c("A1","B1"))
ahi1393=which(ah_data$batch=="2" & (ah_data$R2barcode %in% c("5","6","7","8"))&ah_data$R1barcode %in% c("A2","B2"))
ahi1409=which(ah_data$batch=="2" & (ah_data$R2barcode %in% c("5","6","7","8"))&ah_data$R1barcode %in% c("A3","B3"))
ahi1386=which(ah_data$batch=="2" & (ah_data$R2barcode %in% c("9","10","11","12"))&ah_data$R1barcode %in% c("A1","B1"))
ahi1408=which(ah_data$batch=="2" & (ah_data$R2barcode %in% c("9","10","11","12"))&ah_data$R1barcode %in% c("A2","B2"))
ahi1422=which(ah_data$batch=="2" & (ah_data$R2barcode %in% c("9","10","11","12"))&ah_data$R1barcode %in% c("A3","B3"))
#patients id for batch 3
i1426=which(ah_data$batch=="3" & ah_data$R2barcode=="1" &
            (ah_data$R1barcode %in% c("A1","B1","A2","B2","A3","B3","A4","B4")))
i1431=which(ah_data$batch=="3" & ah_data$R2barcode=="2" &
            (ah_data$R1barcode %in% c("A1","B1","A2","B2","A3","B3","A4","B4")))
i1433=which(ah_data$batch=="3" & ah_data$R2barcode=="3" &
            (ah_data$R1barcode %in% c("A1","B1","A2","B2","A3","B3","A4","B4")))
i1440=which(ah_data$batch=="3" & ah_data$R2barcode %in% c("8","10") &
            (ah_data$R1barcode %in% c("A1","B1","A2","B2","A3","B3","A4","B4")))
i1443=which(ah_data$batch=="3" & ah_data$R2barcode=="9" &
              (ah_data$R1barcode %in% c("A1","B1","A2","B2","A3","B3","A4","B4")))
i1486=which(ah_data$batch=="3" & ah_data$R2barcode=="5" &
            (ah_data$R1barcode %in% c("A5","B5","A6","B6","A7","B7","A8","B8")))
i1409=which(ah_data$batch=="3" & ah_data$R2barcode=="6" &
            (ah_data$R1barcode %in% c("A5","B5","A6","B6","A7","B7","A8","B8")))
i1422=which(ah_data$batch=="3" & ah_data$R2barcode=="7" &
            (ah_data$R1barcode %in% c("A5","B5","A6","B6","A7","B7","A8","B8")))
i1450=which(ah_data$batch=="3" & ah_data$R2barcode=="10" &
            (ah_data$R1barcode %in% c("A6","B6","A7","B7","A8","B8","A9","B9")))
i1451=which(ah_data$batch=="3" & ah_data$R2barcode=="11" &
            (ah_data$R1barcode %in% c("A6","B6","A7","B7","A8","B8","A9","B9")))
i1453=which(ah_data$batch=="3" & ah_data$R2barcode=="12" &
            (ah_data$R1barcode %in% c("A6","B6","A7","B7","A8","B8","A9","B9")))
ip1440=which(ah_data$batch=="4" & 
            ((ah_data$R2barcode %in% c("8","10") & ah_data$R1barcode %in% c("A1","B1","A2","B2","A3","B3","A4","B4"))|
            (ah_data$R2barcode == "5" & ah_data$R1barcode %in% c("A5","B5","A6","B6","A7","B7","A8","B8"))))
ip1443=which(ah_data$batch=="4" & 
            ((ah_data$R2barcode =="9" & ah_data$R1barcode %in% c("A1","B1","A2","B2","A3","B3","A4","B4"))|
            (ah_data$R2barcode == "6" & ah_data$R1barcode %in% c("A5","B5","A6","B6","A7","B7","A8","B8"))))
ip1465=which(ah_data$batch=="4" & ah_data$R2barcode =="7" & 
               ah_data$R1barcode %in% c("A5","B5","A6","B6","A7","B7","A8","B8"))

#assign patients id to ah_data
ah_data$patient=NA
ah_data$patient[ii1357]="CD1357"
ah_data$patient[ii1358]="CD1358"
ah_data$patient[ii1364]="CD1364"
ah_data$patient[ii1368]="CD1368"
ah_data$patient[ii1370]="CD1370"
ah_data$patient[ii1386]="CD1386"
ah_data$patient[ii1390]="CD1390"
ah_data$patient[ii1393]="CD1393"
ah_data$patient[ii1341]="CD1341"
ah_data$patient[ii1347]="CD1347"
ah_data$patient[ii1331]="CD1331"
ah_data$patient[ii1335_7]="TCC1335.7"
ah_data$patient[ii1005_2_58]="TCC1005.2.58"
ah_data$patient[ahi1357]="CD1357"
ah_data$patient[ahi1393]="CD1393"
ah_data$patient[ahi1409]="CD1409"
ah_data$patient[ahi1386]="CD1386"
ah_data$patient[ahi1408]="CD1408"
ah_data$patient[ahi1422]="CD1422"
ah_data$patient[i1426]="CD1426"
ah_data$patient[i1431]="CD1431"
ah_data$patient[i1433]="CD1433"
ah_data$patient[i1440]="CD1440"
ah_data$patient[i1486]="CD1486"
ah_data$patient[i1409]="CD1409"
ah_data$patient[i1422]="CD1422"
ah_data$patient[i1450]="CD1450"
ah_data$patient[i1451]="CD1451"
ah_data$patient[i1453]="CD1453"
ah_data$patient[ip1440]="CD1440"
ah_data$patient[ip1443]="CD1443"
ah_data$patient[ip1465]="CD1465"
ah_data<-ah_data[-which(ah_data$patient %in% NA),] #(95,593obs > 95,153 obs)
#dilution rate for batch 1
id108=which(ah_data$batch=="1" & (
              (ah_data$R2barcode =="1" & ah_data$R1barcode %in% c("A1","B1"))|
              (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A2","B2"))|
              (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A5","B5"))|
              (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A6","B6"))|
              (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A7","B7"))|
              (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A10","B10"))|
              (ah_data$R2barcode=="4" & ah_data$R1barcode %in% c("A9","B9"))|
              (ah_data$R2barcode=="5" & ah_data$R1barcode %in% c("A10","B10"))|
              (ah_data$R2barcode=="7" & ah_data$R1barcode %in% c("A9","B9"))|
              (ah_data$R2barcode=="9" & ah_data$R1barcode %in% c("A10","B10"))|
              (ah_data$R2barcode=="10" & ah_data$R1barcode %in% c("A9","B9"))|
              (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A9","B9"))|
              (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A9","B9"))|
              (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A9","B9"))|
              (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A3","B3"))|
              (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A4","B4"))))
id36=which(ah_data$batch=="1" & (
             (ah_data$R2barcode =="2" & ah_data$R1barcode %in% c("A1","B1"))|
             (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A2","B2"))|
             (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A4","B4"))|
             (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A5","B5"))|
             (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A6","B6"))|
             (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A8","B8"))|
             (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A10","B10"))|
             (ah_data$R2barcode=="5" & ah_data$R1barcode %in% c("A9","B9"))|
             (ah_data$R2barcode=="6" & ah_data$R1barcode %in% c("A10","B10"))|
             (ah_data$R2barcode=="8" & ah_data$R1barcode %in% c("A9","B9"))|
             (ah_data$R2barcode=="10" & ah_data$R1barcode %in% c("A10","B10"))|
             (ah_data$R2barcode=="11" & ah_data$R1barcode %in% c("A9","B9"))))
id18=which(ah_data$batch=="1" & (
             (ah_data$R2barcode == "3" & ah_data$R1barcode %in% c("A1","B1"))|
             (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A3","B3"))|
             (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A4","B4"))|
             (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A5","B5"))|
             (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A7","B7"))|
             (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A8","B8"))|
             (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A10","B10"))|
             (ah_data$R2barcode=="6" & ah_data$R1barcode %in% c("A9","B9"))|
             (ah_data$R2barcode=="7" & ah_data$R1barcode %in% c("A10","B10"))|
             (ah_data$R2barcode=="9" & ah_data$R1barcode %in% c("A9","B9"))|
             (ah_data$R2barcode=="11" & ah_data$R1barcode %in% c("A10","B10"))|
             (ah_data$R2barcode=="12" & ah_data$R1barcode %in% c("A9","B9"))))
id9=which(ah_data$batch=="1" & (
            (ah_data$R2barcode == "1" & ah_data$R1barcode %in% c("A2","B2"))|
            (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A3","B3"))|
            (ah_data$R2barcode=="1" & ah_data$R1barcode %in% c("A6","B6"))|
            (ah_data$R2barcode=="2" & ah_data$R1barcode %in% c("A7","B7"))|
            (ah_data$R2barcode=="3" & ah_data$R1barcode %in% c("A8","B8"))|
            (ah_data$R2barcode=="4" & ah_data$R1barcode %in% c("A10","B10"))|
            (ah_data$R2barcode=="8" & ah_data$R1barcode %in% c("A10","B10"))|
            (ah_data$R2barcode=="12" & ah_data$R1barcode %in% c("A10","B10"))))
#dilution rate for batch 3
idd108=which(ah_data$batch=="3" & 
               ((ah_data$R2barcode %in% c("1","2","3","8","9","10") & ah_data$R1barcode %in% c("A1","B1"))|
               (ah_data$R2barcode %in% c("5","6","7") & ah_data$R1barcode %in% c("A5","B5"))|
               (ah_data$R2barcode %in% c("10","11","12") & ah_data$R1barcode %in% c("A6","B6"))))
idd36=which(ah_data$batch=="3" & 
              ((ah_data$R2barcode %in% c("1","2","3","8","9","10") & ah_data$R1barcode %in% c("A2","B2"))|
                 (ah_data$R2barcode %in% c("5","6","7") & ah_data$R1barcode %in% c("A6","B6"))|
                 (ah_data$R2barcode %in% c("10","11","12") & ah_data$R1barcode %in% c("A7","B7"))))
idd18=which(ah_data$batch=="3" & 
              ((ah_data$R2barcode %in% c("1","2","3") & ah_data$R1barcode %in% c("A3","B4"))|
                 (ah_data$R2barcode %in% c("8","9","10") & ah_data$R1barcode %in% c("A3","B3"))|
                 (ah_data$R2barcode %in% c("5","6","7") & ah_data$R1barcode %in% c("A7","B7"))|
                 (ah_data$R2barcode %in% c("10","11","12") & ah_data$R1barcode %in% c("A8","B8"))))
idd9=which(ah_data$batch=="3" & 
              ((ah_data$R2barcode %in% c("1","2","3") & ah_data$R1barcode %in% c("A4","B3"))|
                 (ah_data$R2barcode %in% c("8","9","10") & ah_data$R1barcode %in% c("A4","B4"))|
                 (ah_data$R2barcode %in% c("5","6","7") & ah_data$R1barcode %in% c("A8","B8"))|
                 (ah_data$R2barcode %in% c("10","11","12") & ah_data$R1barcode %in% c("A9","B9"))))
#assign dilution to wells batch 1
ah_data$dilution=NA
ah_data$dilution[id108]="d1080"
ah_data$dilution[id36]="d360"
ah_data$dilution[id18]="d180"
ah_data$dilution[id9]="d90"
#assign dilution to wells batch 2
ah_data$dilution[which(ah_data$batch=="2" & ah_data$R2barcode %in% c("5","9"))]="d1080"
ah_data$dilution[which(ah_data$batch=="2" & ah_data$R2barcode %in% c("6","10"))]="d360"
ah_data$dilution[which(ah_data$batch=="2" & ah_data$R2barcode %in% c("7","11"))]="d180"
ah_data$dilution[which(ah_data$batch=="2" & ah_data$R2barcode %in% c("8","12"))]="d90"
#assign dilution to wells batch 3
ah_data$dilution[idd108]="d1080"
ah_data$dilution[idd36]="d360"
ah_data$dilution[idd18]="d180"
ah_data$dilution[idd9]="d90"
#assign dilution to wells batch 4
ah_data$dilution[which(ah_data$batch=="4" & ah_data$R1barcode %in% c("A1","B1","A5","B5"))]="d1080"
ah_data$dilution[which(ah_data$batch=="4" & ah_data$R1barcode %in% c("A2","B2","A6","B6"))]="d360"
ah_data$dilution[which(ah_data$batch=="4" & ah_data$R1barcode %in% c("A3","B3","A7","B7"))]="d180"
ah_data$dilution[which(ah_data$batch=="4" & ah_data$R1barcode %in% c("A4","B4","A8","B8"))]="d90"
ah_data<-ah_data[-which(ah_data$dilution %in% NA),] #(70,814 obs > 70,787 obs)
#write.table(ah_data, file = "/Users/yingy_adm/Documents/TCRdata/rmdupc_4batch.csv")
write.table(ah_data, file = "/Users/yingy_adm/Documents/TCRdata/rmdupc.csv")

control1408<-subset(rmdupc_3batch,patient=="CD1408")[,8]
p1357<-subset(rmdupc_3batch,patient=="CD1357")[,8]
write.table(control1408, file = "/Users/yingy_adm/Documents/TCRdata/control1408.txt")
write.table(p1357, file = "/Users/yingy_adm/Documents/TCRdata/p1357.txt")
##########################

#assign patient id to new column and remove the NA rows
# i1357=which((a5_data$R2barcode %in% c("5","6","7","8"))&a5_data$R1barcode %in% c("A1","B1"))
# i1393=which((a5_data$R2barcode %in% c("5","6","7","8"))&a5_data$R1barcode %in% c("A2","B2"))
# i1409=which((a5_data$R2barcode %in% c("5","6","7","8"))&a5_data$R1barcode %in% c("A3","B3"))
# i1386=which((a5_data$R2barcode %in% c("9","10","11","12"))&a5_data$R1barcode %in% c("A1","B1"))
# i1408=which((a5_data$R2barcode %in% c("9","10","11","12"))&a5_data$R1barcode %in% c("A2","B2"))
# i1422=which((a5_data$R2barcode %in% c("9","10","11","12"))&a5_data$R1barcode %in% c("A3","B3"))

# a5_data$patient=NA
# a5_data$patient[i1357]="CD1357"
# a5_data$patient[i1393]="CD1393"
# a5_data$patient[i1409]="CD1409"
# a5_data$patient[i1386]="CD1386"
# a5_data$patient[i1408]="CD1408"
# a5_data$patient[i1422]="CD1422"
# a5<-a5_data[-which(a5_data$patient %in% NA),]

ah<-ah_data[-which(ah_data$patient %in% NA),]

ah$dilution[which(ah$R2barcode %in% c("5","9"))]="d1080"
ah$dilution[which(ah$R2barcode %in% c("6","10"))]="d360"
ah$dilution[which(ah$R2barcode %in% c("7","11"))]="d180"
ah$dilution[which(ah$R2barcode %in% c("8","12"))]="d90"

#ah$Vgene=str_split_fixed(ah$ah, "[|]", 6)[,4]
#ah$Jgene=str_split_fixed(ah$ah, "[|]", 6)[,5]
#ah$CDR3=str_split_fixed(ah$ah, "[|]", 6)[,6]
#ah<-ah[,-1]

#count number of unique TCR for each patients(220616 dupcount>0)
a5$flag5=str_split_fixed(a5$flag5,"[|]",3)[,3]
s220616_0<-paste(a5$patient,a5$flag5,sep="_")
su220616_0<-unique(s220616_0)
table(str_split_fixed(su220616_0,"_",2)[,1])
#table(str_split_fixed(su220616_0,"V",2)[,1])
ah$tcr<-paste(ah$V.GENE,ah$J.GENE,ah$CDR3nt,sep="_")
s220616_0<-paste(ah$patient,ah$tcr,sep="_")
su220616_0<-unique(s220616_0)
table(str_split_fixed(su220616_0,"_",2)[,1])

#search for obs.with flag=="PRIMER=B1,5,G,|BARCODE=CGCTTA|Homsap TRAV20|Homsap TRBJ2-1*01 F|ASSLGVALSSYNEQF"
#example<-tcrdata[which(tcrdata$flag=="PRIMER=B1,5,G,|BARCODE=CGCTTA|Homsap TRAV20|Homsap TRBJ2-1*01 F|ASSLGVALSSYNEQF"),]
#ex1<-tcrdata[which(tcrdata$flag4=="PRIMER=B1,5,G,|BARCODE=CGCTTA|Homsap TRAV20|Homsap TRBJ2-1*01 F|ASSLGVALSSYNEQF"),]
#ex2<-tcrdata[which(tcrdata$flag5=="PRIMER=B1,5,G,|BARCODE=AAAACT|Homsap TRAV20|Homsap TRBJ1-5*01 F|ASAGDSQTQH|Homsap TRBD2*01 F"),]
###################************************
# library(raster)#coefficient of variation (cv)&&
#count number of unique TCR for each well(patient and dilution) in Batch 2 (lib=d220616)
 ah_data$tcr<-paste(ah_data$V.GENE,ah_data$J.GENE,ah_data$CDR3nt,sep="_")
 mat<-matrix(0,6,8)
 abmat<-matrix(0,6,4)
 R1=list("A1","B1","A2","B2","A3","B3")
 R2=list("5","6","7","8","9","10","11","12")
 cd=list("CD1357","CD1393","CD1409","CD1386","CD1408","CD1422")
 dilution=list("d1080","d360","d180","d90")
# cvmat<-matrix(0,6,8)
# barcode12<-rep(NA,48)
# cv<-rep(0,48)
i=j=k=0
for (r1 in R1)
# for (r1 in cd)
   {i=i+1
    j=0
  for (r2 in R2)
# for (r2 in dilution)
   {j=j+1
#   k=k+1
#   ss=subset(ah,patient==r1 & dilution==r2)
    ss=subset(ah_data,R1barcode==r1 & R2barcode==r2 & batch=="2")
    mat[i,j]=length(unique(ss$tcr))
#   abmat[i,j]=length(unique(ss$tcr))
#   print r1
#   print r2
#   print table(ss$RowBarcode)
#   a<-floor(cv(table(ss$RowBarcode)))
#   cvmat[i,j]<-a
#   barcode12[k]<-paste(r1,r2,sep = '_')
#   cv[k]<-a
        }
 }
# cvdf<-cbind(barcode12,cv)
colnames(mat)=c("5","6","7","8","9","10","11","12")
rownames(mat)=c("A1","B1","A2","B2","A3","B3")
colnames(abmat)=c("d1080","d360","d180","d90")
rownames(abmat)=c("CD1357","CD1393","CD1409","CD1386","CD1408","CD1422")
write.table(a_bmat2, file = "/Users/yingy_adm/Documents/TCRdata/no_uniqTCR_for_p+d.csv", sep="\t", row.names = T, qmethod = "double")
# write.table(cvmat, file = "/Users/yingy_adm/Documents/TCRdata/cv.csv", sep="\t", row.names = T, qmethod = "double")
# write.table(scvmat, file = "/Users/yingy_adm/Documents/TCRdata/cv_dupcount2.csv", sep="\t", row.names = T, qmethod = "double")

#count number of unique TCR for each patients+dilution in Batch 1 (lib 160316 dupcount>0)
ah_data$dilutionab<-paste(ah_data$dilution,substr(ah_data$R1barcode,1,1),sep="")
abmat1<-matrix(0,11,4)
#a_bmat1<-matrix(0,11,8)
cd_batch1=list("CD1357","CD1358","CD1364","CD1368","CD1370","CD1386","CD1390","CD1393","CD1341","CD1347","CD1331")
dilutionAB=list("d1080A","d1080B","d360A","d360B","d180A","d180B","d90A","d90B")
dilution=list("d1080","d360","d180","d90")

i=j=0
for (r1 in cd_batch1)
{i=i+1
j=0
for (r2 in dilutionAB)
  #  for (r2 in dilution)
{j=j+1
      ss=subset(ah_data, patient==r1 & dilutionab==r2 & batch=="1")
      #ss=subset(aah,patient==r1 & dilution==r2)
      a_bmat1[i,j]=length(unique(ss$tcr))
#abmat1[i,j]=length(unique(ss$flag5))
  }
}
rownames(a_bmat1)<- cd_batch1
#rownames(abmat1)<- cd_batch1
colnames(a_bmat1)<-dilutionAB
#colnames(abmat1)<-dilution

##############*************************
# plot of vsdf
#library("ggplot2")
#qplot(cvdf[,1], cvdf[,2],col="red")
#ggplot(diamonds, aes(clarity, fill=cut)) + geom_bar(position="dodge")
#ggplot(cvdf, aes(cvdf[,1], fill=cut)) + geom_bar(position="dodge")

sa5<-subset(a5,dupcount>1)
#sa5$flag5=str_split_fixed(sa5$flag5,"[|]",3)[,3]
s220616_1<-paste(sa5$patient,sa5$flag5,sep="_")
su220616_1<-unique(s220616_1)
table(str_split_fixed(su220616_1,"_",2)[,1])
table(str_split_fixed(su220616_1,"V",2)[,1])
write.table(sa5[,c("patient","tcrflag","dilution")], file = "/Users/yingy_adm/Documents/TCRdata/s220616.csv", sep="\t", row.names = T, qmethod = "double")
#library(raster)#coefficient of variation (cv)
smat<-matrix(0,6,8)
#scvmat<-matrix(0,6,8)
#sbarcode12<-rep(NA,48)
#scv<-rep(0,48)
i=j=k=0
for (r1 in R1) {
  i=i+1
  j=0
  for (r2 in R2) {j=j+1
#  k=k+1
  ss=subset(sa5,R1barcode==r1 & R2barcode==r2)
#  a<-floor(cv(table(ss$RowBarcode)))
#  scvmat[i,j]<-a
#  sbarcode12[k]<-paste(r1,r2,sep = '_')
#  scv[k]<-a
  #ss$flag=paste(ss$,ss$UMI)
  #t=table(ss$flag)
  smat[i,j]=dim(ss)[1]
  }
}
#cvsdf<-cbind(rep("dupcount>1",48),sbarcode12,scv)
#cvdataframe<-rbind(cvdf,cvsdf)
colnames(smat)=colnames(cvsmat)=c("5","6","7","8","9","10","11","12")
rownames(smat)=rownames(cvsmat)=c("A1","B1","A2","B2","A3","B3")
write.table(smat, file = "/Users/yingy_adm/Documents/TCRdata/no_uniqTCR_1.csv", sep="\t", row.names = T, qmethod = "double")
# plot of cvsdf
#library("ggplot2")
#qplot(cvsdf[,1], cvsdf[,2],main="CV(Dupcount>1)")

# matrix for ACTG count from UMI
 uniqumi<-dim(table(a5$UMI))
 umi<-str_split_fixed(a5$UMI, "",6)
 umimat<-cbind(table(umi[,1]),table(umi[,2]),table(umi[,3]),table(umi[,4]),table(umi[,5]),table(umi[,6]))
 umimat<-umimat/159020
# plot for ACTG count from UMI
barplot(umimat, main="Frequency of ATCG in UMI",
        xlab="Position in UMI", col=c("blue","red","green","yellow"),)
par(xpd=TRUE)
legend(7,1,rownames(umimat), col=c("blue","red","green","yellow"),pch = 15)

#uniq tcr counting matrix for all
tcr<-str_split_fixed(a5$flag5,"[|]",8)
a5$tcrflag<-paste(tcr[,5],tcr[,6],tcr[,7],tcr[,8],sep = '|')
no_tcr_mat<-matrix(0,6,8)
i=j=0
for (r1 in R1) {
  i=i+1
  j=0
  for (r2 in R2) {j=j+1
  ss=subset(a5,R1barcode==r1 & R2barcode==r2)
  t=table(ss$tcrflag5)
  no_tcr_mat[i,j]=dim(t)
  }
}
colnames(no_tcr_mat)=c("5","6","7","8","9","10","11","12")
rownames(no_tcr_mat)=c("A1","A2","A3","B1","B2","B3")
write.table(no_tcr_mat, file = "/Users/yingy_adm/Documents/TCRdata/no_uniq_TCR.csv", sep="\t", row.names = T, qmethod = "double")

#uniq tcr counting matrix for those dupcount>1
stcr<-str_split_fixed(sa5$flag5,"[|]",8)
sa5$tcrflag<-paste(stcr[,5],stcr[,6],stcr[,7],stcr[,8],sep = '|')
no_tcr_smat<-matrix(0,6,8)
i=j=0
for (r1 in R1) {
  i=i+1
  j=0
  for (r2 in R2) {j=j+1
  ss=subset(sa5,R1barcode==r1 & R2barcode==r2)
  t=table(ss$tcrflag)
  no_tcr_smat[i,j]=dim(t)
  }
}
colnames(no_tcr_smat)=c("5","6","7","8","9","10","11","12")
rownames(no_tcr_smat)=c("A1","A2","A3","B1","B2","B3")
write.table(no_tcr_smat, file = "/Users/yingy_adm/Documents/TCRdata/no_uniq_sTCR.csv", sep="\t", row.names = T, qmethod = "double")

#find all sequence that were found in more than one individuals
cd1357<-unique(subset(a5,patient=="CD1357")$tcrflag)
cd1393<-unique(subset(a5,patient=="CD1393")$tcrflag)
cd1409<-unique(subset(a5,patient=="CD1409")$tcrflag)
cd1386<-unique(subset(a5,patient=="CD1386")$tcrflag)
cd1408<-unique(subset(a5,patient=="CD1408")$tcrflag)
cd1422<-unique(subset(a5,patient=="CD1422")$tcrflag)
fre_tcr<-table(c(cd1357,cd1393,cd1409,cd1386,cd1408,cd1422))
write.table(fre_tcr, file = "/Users/yingy1/Documents/TCRdata/TCRfrequency_patient.csv", sep="\t", row.names = T, qmethod = "double")

a5$UMIBC<-paste(tcr[,1],tcr[,2],tcr[,3],tcr[,4],sep = '|')
m<-data.frame(NA,2575,3)#initiate list for all the 343 wells
j=0
well=unique(a5$wellflag)
for (i in 1:length(well))#343 of 384 possible combination
{ ss<-subset(a5, a5$wellflag==well[i])
umifortcr<-aggregate(UMIBC ~ tcrflag, data = ss, c)
umifortcr$no_umi<-lengths(umifortcr$UMIBC, use.names = TRUE)
j=j+dim(umifortcr)[1]
m[1:j,1]<-"well"
m[1:j,2]<-umifortcr$no_umi
m[1:j,3]<-umifortcr$UMIBC
}
#countRowbarcode2 for QC
countRowbarcode2<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(countRowbarcode2) <- paste0("Well","Rowbarcode","Freq")
for (r1 in R1)
{
  for (r2 in R2) 
  {
    ss=subset(sa5,R1barcode==r1 & R2barcode==r2)
    t = as.data.frame(table(ss$RowBarcode))
    well = rep(paste(r1,r2,sep="_"),dim(t)[1])
    t = cbind(well,t)
    countRowbarcode2=rbind(countRowbarcode2,t)
  }
}
write.table(countRowbarcode2, file = "/Users/yingy_adm/Documents/TCRdata/countRowbarcode2.csv")
