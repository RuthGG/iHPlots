#!/usr/bin/R
# Carla Giner
# 25.10.16

args = c("tmp/2020-11-19_02_ihplot/00_Complete/HsInv1057/positions.txt" ,
         "tmp/2020-11-19_02_ihplot/00_Complete/HsInv1057/haplotypes.txt" ,
         "tmp/2020-11-19_02_ihplot/00_Complete/HsInv1057/samples.txt",
         2, "tmp/2020-11-19_02_ihplot/01_Region/HsInv0389/", 0)


# READ OPTIONS ############################################################################################################
    args = commandArgs(trailingOnly=TRUE)
    # test if there is at least one argument: if not, return an error
    Usage<-function(){
        write("\nUsage:\nRscript SelectRegion.R <PositionsFile> <HaplotypesFile> <SamplesFile> <MinAlleleCount> <Prefix> <FlankingRegion>", stderr())
        q()
    }

    if (length(args)==0 | length(args)<5 | args[1]=="help") {
    
       # ***************************************************
       Usage()
       # ***************************************************

    } else {
        PositionsFile<-args[1]
        HaplotypesFile<-args[2]
        SamplesFile<-args[3]
        MinCount<-as.numeric(args[4])
        Prefix<-args[5]
        FlankingRegion<-if(args[6]>0){as.integer(args[6])}else{0}
    }

    # ***************************************************
    # write(paste0("\n",
    #     "\n\tPositions file: ",PositionsFile,
    #     "\n\tHaplotypes file: ",HaplotypesFile,
    #     "\n\tSamples file: ",SamplesFile,
    #     "\n\tMinimum Allele count to be included: ",MinCount,
    #     "\n\tOutput prefix: ",Prefix,
    #     "\n\tFlanking region included (bps): ",FlankingRegion,"\n"),stdout())
    # ***************************************************

# READ DATA ###############################################################################################################

    # Read data, create 3 data structures: Haplotypes matrix, Positions info, Samples info
    Positions<-read.table(PositionsFile,header=TRUE,stringsAsFactors=FALSE)
    rownames(Positions)<-Positions$ID
    Haplotypes<-read.table(HaplotypesFile,stringsAsFactors=FALSE)
    Samples<-read.table(SamplesFile,header=TRUE,stringsAsFactors=FALSE)

# SUBSET INSIDE AND NO-SINGLETON ##########################################################################################
    # Exclude less than X counts (and positions with NA in the AltCount.)
    # Update AltCount in positions file, just in case
    Positions$AltCount<-apply(Haplotypes,2,function(v){sum(as.integer(v),na.rm=TRUE)})

    CountThreshold<-sapply(Positions$AltCount,function(x){if(is.na(x)){FALSE}else if(x<MinCount | (nrow(Haplotypes)-x)<MinCount){FALSE}else{TRUE}})

    SelectedCols<-Positions[with(Positions, !(VariantType %in% "Breakpoint") & abs(as.integer(RelativePosition))<=FlankingRegion & CountThreshold),"ID"]

    if(length(SelectedCols)==0){
        write("\nThere are no SNPs inside the inversion!\n",stdout())
        q()
    }else if(length(SelectedCols)==1){
        write("\nThere is only one SNP inside the inversion!\n",stdout())
        q()
    }
    # Make a table with haplotype strings + ID.
    StringHaps<-apply(Haplotypes[,SelectedCols],1,paste,collapse="")
    tempTable<-table(StringHaps)
    HaplotypesTable<-data.frame(HapID=factor(paste0("Hap",seq(1:nrow(tempTable))),levels=paste0("Hap",seq(1:nrow(tempTable)))),Hap01=names(tempTable),Count=as.vector(tempTable),stringsAsFactors=FALSE)
    rownames(HaplotypesTable)<-HaplotypesTable$Hap01

    Samples$HapID<-HaplotypesTable[StringHaps,"HapID"]
    
    library(reshape2)
    library(plyr)
    GenoCounts<-dcast(Samples,formula=list(.(HapID),.(Genotype)),value.var="HapID",fun.aggregate=length)
    if( is.null(GenoCounts$I )){GenoCounts$I<-0}
    if( is.null(GenoCounts$S )){GenoCounts$I<-0}
    HaplotypesTable<-cbind(HaplotypesTable,GenoCounts[, c("S","I")])

    rownames(HaplotypesTable)<-NULL
    HaplotypesTable$Shared<-apply(HaplotypesTable[,c("S","I")],1,function(v){v[1]!=0 & v[2]!=0})

     
    HaplotypesTable$HapATGC<-sapply(HaplotypesTable$Hap01,function(s){
        v<-strsplit(s,split="")[[1]]
        newv<-v
        for (i in 1:length(SelectedCols)) {
            AA<-Positions[SelectedCols[i],"AncestralAllele"]
            Allele<-v[i]
            newv[i]<-if(is.na(AA)){
                if(Allele == "0"){
                    "G"
                }else if(Allele == "1"){
                    "C"
                }else{
                    "N"
                }
            }else if(AA == 0){
                if(Allele == "0"){
                    "A"
                }else if(Allele == "1"){
                    "T"
                }else{
                    "N"
                }
            }else if(AA == 1){
                if(Allele == "1"){
                    "A"
                }else if(Allele == "0"){
                    "T"
                }else{
                    "N"
                }
            }else{
                "N"
            }
        }
        paste(newv,collapse="")
    })

    ReferenceHaplotype<-paste0(ifelse(is.na(Positions[SelectedCols,"AncestralAllele"]),"G","A"),collapse="")

# SAVE DATA ###############################################################################################################
    #  1. HaplotypesInside.fasta ">Ref\nAAAAAACAAAAAAACAAA"">Hap1-SI\nATACAG"
    #  2. SamplesInfo.txt "Indiv_1\t[I/S]\tHap1-SI"
    #  3. PositionsInfo.txt Positions included with BP etc...
    #  4. HaplotypeInfo.txt Shared, count...

    # HaplotypesInside.fasta
    write(">Reference",paste0(Prefix,"Haplotypes.fasta"))
    write(ReferenceHaplotype,paste0(Prefix,"Haplotypes.fasta"),append=TRUE)
    
    for (j in 1:nrow(HaplotypesTable)) {
        write(paste0(">",HaplotypesTable[j,"HapID"],"_S",HaplotypesTable[j,"S"],"_I",HaplotypesTable[j,"I"]),paste0(Prefix,"Haplotypes.fasta"),append=TRUE)
        write(HaplotypesTable[j,"HapATGC"],paste0(Prefix,"Haplotypes.fasta"),append=TRUE)
    }

    # SamplesInfo.txt
    write.table(Samples,file=paste0(Prefix,"SampleInfo.txt"),row.names=FALSE,sep="\t")
    # PositionsInfo.txt
    write.table(Positions[SelectedCols,c("ID","Chr","Position","rsID","Ref","Alt","VariantType","RelativePosition","AncestralAllele","AltCount")],file=paste0(Prefix,"PositionInfo.txt"),row.names=FALSE,sep="\t",quote=FALSE)
    # HaplotypeInfo.txt
    write.table(HaplotypesTable[,c("HapID","S","I","Shared","Hap01","HapATGC")],file=paste0(Prefix,"HaplotypeInfo.txt"),row.names=FALSE,sep="\t")
    
q()
