#!/usr/bin/R
# Ruth Gómez, based on another script by Carla Giner
# 19.11.2020

# args = c("6:167610402" , "6:167771360" ,
#          "analysis/2020-11-19_01_phase/HsInv1057/HsInv1057_phased.vcf",
#          "analysis/2020-11-19_01_phase/HsInv1057/1KGP_haplotypes.txt",
#          NA, "tmp/2020-11-19_02_ihplot/HsInv1057/")

# READ OPTIONS ############################################################################################################
    args = commandArgs(trailingOnly=TRUE)

      
# test if there is at least one argument: if not, return an error
Usage<-function(){
    write("\nUsage:\nRscript Tidy_Ph3_Haplotypes.R <BP1 Chr:Position> <BP2 Chr:Position> <PhasedInvFile> <HaplotypesFile> <AncOrientation> <Prefix>", stderr())
    q()
}

if (length(args)==0 | length(args)<6 | args[1]=="help") {

   # ***************************************************
   Usage()
   # ***************************************************

} else {
        BP1Pos<-args[1]
        BP2Pos<-args[2]
        PhasedFile<-args[3]
        HaplotypesFile<-args[4]
        ancOrient<-args[5]
        Prefix<-args[6]
    }


# READ DATA AND CREATE TIDY DATA STRUCTURES ################################################################################
    
    # Inversion coordinates info
    Chromosome<-strsplit(BP1Pos,split=":")[[1]][1]
    BP1<-as.integer(strsplit(BP1Pos,split=":")[[1]][2])
    BP2<-as.integer(strsplit(BP2Pos,split=":")[[1]][2])
    
    # Data from 1000 genomes 
    Haplotypes<-read.table(HaplotypesFile,header=TRUE,stringsAsFactors=FALSE)
    
    # PHASED GENOTYPES
    # ----------------------------------------------------------------------
    
    # Phased Inversion data -> parse to give the same format as Haplotypes
    lenPhased <- length(readLines(PhasedFile))
    Phased<-read.table(PhasedFile, stringsAsFactors=FALSE, comment.char = "", skip = lenPhased-2, header = TRUE)
    colnames(Phased)[c(1:7)]<-c("CHROM","POS","ID","REF","ALT","VT","AA")
    Phased[1,]<-gsub(":.*$" , "", Phased)
  
    Males<-read.table("data/raw/1KGP/integrated_call_male_samples_v3.20130502.ALL.panel", header = TRUE)
    
    if (Chromosome == "X"){
      for(m in Males$sample){
        if(m %in% colnames(Phased)){
          Phased[,m] <- ifelse(Phased[,m] == "0|0", "0", "1")
        }
      }
    }
    
    # Fill in some gaps in the phased info
    Phased[1,c("CHROM","POS","ID","REF","ALT","VT","AA")]<-c(Chromosome,BP1,NA,"Std","Inv","Breakpoint",ancOrient)
    Phased[2,]<-Phased[1,]
    Phased[2,c("POS")]<-c(BP2)

    # Add phase info to general info
    Haplotypes<-rbind(Haplotypes, Phased[, colnames(Haplotypes)])
    
    # Check if there are duplicated lines. If so, remove them
    Haplotypes<-Haplotypes[!duplicated(Haplotypes),]
    
    # HAPLOTYPES 
    # ----------------------------------------------------------------------
    DiploidHaps<-Haplotypes[,setdiff(names(Haplotypes),c("CHROM","POS","ID","REF","ALT","VT","AA"))]
    Ploidy<-sapply(DiploidHaps[1,],function(x){if(x %in% c(0,1)){"Haploid"}else if(x %in% c("0|0","0|1","1|0","1|1")){"Diploid"}else{NA}})
    
    HaploidHaps<-NULL
    HaploidNames<-NULL
    for (i in 1:ncol(DiploidHaps)) {
        if(Ploidy[i] %in% "Diploid"){
            IndivName<-names(Ploidy[i])
            NewHaps<-t(sapply(DiploidHaps[,i],function(x)strsplit(x,split="|",fixed=TRUE)[[1]]))
            rownames(NewHaps)<-seq(1,nrow(NewHaps),1)
            HaploidHaps<-cbind(HaploidHaps,NewHaps)
            HaploidNames<-c(HaploidNames,paste0(IndivName,"_1"),paste0(IndivName,"_2"))
        }else{
            IndivName<-names(Ploidy[i])
            HaploidHaps<-cbind(HaploidHaps,DiploidHaps[,i])
            HaploidNames<-c(HaploidNames,paste0(IndivName,"_1"))
        }
    }
    colnames(HaploidHaps)<-HaploidNames


    # POSITIONS
    # ----------------------------------------------------------------------
    Positions<-Haplotypes[,c("CHROM","POS","ID","REF","ALT","VT","AA")]
    names(Positions)<-c("Chr","Position","rsID","Ref","Alt","VariantType","AA")

    # Add Relative position, AltCount and AncestralAllele[0|1]
    Positions$Position<-as.numeric(Positions$Position)
    Positions$RelativePosition<-sapply(Positions$Position,function(x){
        if(x<BP1){x-BP1}else if(x>=BP1 & x<=BP2){0}else{x-BP2}
        })

    Positions$AncestralAllele<-apply(Positions[,c("Ref","Alt","AA")],1,function(v){
            v[3]<-toupper(strsplit(v[3],split="|")[[1]][1])
            v<-toupper(v)
            if(v[3] %in% v[1]){
                0
            }else if(v[3] %in% v[2]){
                1
            }else{
                NA
            }
            })
    
    Positions$AncestralAllele[Positions$VariantType %in% "Breakpoint"]<-as.integer(as.character(factor(ancOrient,levels=c("Std","Inv"),labels=c("0","1"))))

    # Order hap matrix and Positions According to Position
    HaploidHaps<-HaploidHaps[order(Positions$Position),]

    Positions<-Positions[order(Positions$Position),]
    Positions$ID<-paste0("Pos",1:nrow(Positions))

    rownames(HaploidHaps)<-Positions$ID
    Positions$AltCount<-apply(HaploidHaps,1,function(v)sum(as.numeric(v),na.rm=TRUE))

    # SAMPLES 
    # ----------------------------------------------------------------------

    # Make samples from haplotypes
    Samples<-data.frame(Sample=HaploidNames,stringsAsFactors=FALSE)
    
    # Take individuals ID
    Samples$ID<-sapply(Samples$Sample,function(x){strsplit(x,split="_")[[1]][1]})

    # Take populations from 1KGP
    IndivKGP<-read.table("data/raw/1KGP/integrated_call_samples_v3.20130502.ALL.panel",header=TRUE)[,1:3]
    names(IndivKGP)<-c("sample","population","superpopulation")
    rownames(IndivKGP)<-IndivKGP$sample

    Samples$Population<-IndivKGP[Samples$ID,"population"]
    Samples$Superpopulation<-IndivKGP[Samples$ID,"superpopulation"]

    # Take genotypes info from Phased info : 0 = STD, 1 = INV
    Samples$FullGenotype<-NA
    Samples$GenotypeScore<-NA
    Samples$Genotype<-NA
    
    Samples$FullGenotype<-sapply(Samples$ID,  function(x){ 
      if( Phased[1,x] %in% c("0|0", "0")){"STD"}else if(Phased[1,x] %in% c("1|1", "1")){"INV"}else{"HET"}
      })
    
    Samples$Genotype<-substr(Samples$FullGenotype, 1, 1)
    Samples$Genotype[Samples$Genotype == "H"]<-"SI"
    
    
# WRITE DATA ################################################################################
    

    write.table(Positions[,c("ID","Chr","Position","rsID","Ref","Alt","VariantType","RelativePosition","AncestralAllele","AltCount")],file=paste0(Prefix,"positions.txt"),row.names=FALSE,sep="\t", quote=FALSE)
    write.table(t(HaploidHaps),file=paste0(Prefix,"haplotypes.txt"),sep="\t",quote=FALSE)
    write.table(Samples,file=paste0(Prefix,"samples.txt"),row.names=FALSE,sep="\t",quote=FALSE)


q()
