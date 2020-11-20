#!/usr/bin/R
# Carla Giner
# 25.10.16


# READ OPTIONS ############################################################################################################
    args = commandArgs(trailingOnly=TRUE)
    # test if there is at least one argument: if not, return an error
    Usage<-function(){
        write("\nUsage:\nRscript AnalyseClusters.R <PositionsFile> <HaplotypesFile> <SamplesFile> <MainOrientationThreshold[0-1]> <Prefix> <Mutations[Derived|Minor]> <TreeFraction>", stderr())
        q()
    }

    if (length(args)==0 | length(args)<7 | args[1]=="help") {
    
       # ***************************************************
       Usage()
       # ***************************************************

    } else {
        PositionsFile<-args[1]
        HaplotypesFile<-args[2]
        SamplesFile<-args[3]
        MainThreshold<-as.numeric(args[4])
        Prefix<-args[5]
        Mutations<-args[6]
        MinMutations<-ceiling(eval(parse(text=args[7])))
    }
    if(MainThreshold<=0 | MainThreshold >1){
        write("MainThreshold has to be in the range (0,1]! Setting it to 0.75\n")
        MainThreshold <- 0.75
    }
    # ***************************************************
    write(paste0("\n",
        "\n\tPositions file: ",PositionsFile,
        "\n\tHaplotypes file: ",HaplotypesFile,
        "\n\tSamples file: ",SamplesFile,
        "\n\tFraction of chromosomes to define main orientation in a cluster: ",MainThreshold,
        "\n\tOutput prefix: ",Prefix,
        "\n\tMutations shown in haplotypes: ",Mutations,
        "\n\tMinimum number of mutations clusters need to have: ",MinMutations,"\n"),stdout())
    # ***************************************************

Inv<-regmatches(m=gregexpr(text=PositionsFile,pattern="HsInv[0-9]{4}",perl=TRUE),x=PositionsFile)[[1]]

library(ggplot2)
library(ggdendro)
library(RColorBrewer)
library(reshape2)
library(cowplot)

# # 0. Options
#   inv<-"HsInv0389"
#   MinCount<-"MinCount7"
#   Prefix<-"."
#   FlankingRegion<-0
#   MinMutations<-57

#   # filenames
#   HaplotypesFile<-paste0("results/02_Ph1_HaplotypesInside/",MinCount,"/",inv,"/HaplotypeInfo.txt")
#   PositionsFile<-paste0("results/02_Ph1_HaplotypesInside/",MinCount,"/",inv,"/PositionInfo.txt")
#   SamplesFile<-paste0("results/02_Ph1_HaplotypesInside/",MinCount,"/",inv,"/SampleInfo.txt")

# 1. Load data
    # read files
    HaplotypeInfo<-read.table(HaplotypesFile,header = TRUE,stringsAsFactors = FALSE,colClasses = c("character","integer","integer","logical","character","character"))
    PositionInfo<-read.table(PositionsFile,header = TRUE,stringsAsFactors = FALSE)
    SampleInfo<-read.table(SamplesFile,header = TRUE,stringsAsFactors = FALSE)
    
    # Create Type variable for haplotypes
    HaplotypeInfo$Type<-NA
    HaplotypeInfo$Type<-ifelse(HaplotypeInfo$S>0,"S","I")
    HaplotypeInfo$Type[HaplotypeInfo$Shared]<-"SI"
    
    # Create population variables to know the continent of each haplotype
    PopLevels<-sort(unique(SampleInfo$Population))
    if(identical(PopLevels,c("CEU","CHB","GIH","JPT","LWK","TSI","YRI"))){
      # Desired order:
      PopLevels<-c("LWK","YRI","CHB","JPT","CEU","TSI","GIH")
    }else if(identical(PopLevels,c("CEU","CHB","JPT","LWK","TSI","YRI"))){
      PopLevels<-c("LWK","YRI","CHB","JPT","CEU","TSI")
    }


    HapPopCount<-t(sapply(split(SampleInfo$Population,f=factor(SampleInfo$HapID,levels=HaplotypeInfo$HapID)),function(x){table(factor(x,levels=PopLevels))}))
    HaplotypeInfo<-cbind(HaplotypeInfo,HapPopCount)

    # Create Variables with info about the other chromosome in SampleInfo
    SampleInfo$OtherGenotype<-sapply(SampleInfo$Sample,function(x){
      indivName<-strsplit(x,split="_")[[1]]
      if(indivName[2] == "1"){
        othergeno<-SampleInfo$Genotype[grep(x = SampleInfo$Sample,pattern = paste0(indivName[1],"_","2"))]
      }else{
        othergeno<-SampleInfo$Genotype[grep(x = SampleInfo$Sample,pattern = paste0(indivName[1],"_","1"))]
      }
      if(length(othergeno)<1){
        othergeno<-NA
      }
      othergeno
    })
    
    SampleInfo$OtherHapID<-sapply(SampleInfo$Sample,function(x){
      indivName<-strsplit(x,split="_")[[1]]
      if(indivName[2] == "1"){
        otherhap<-SampleInfo$HapID[grep(x = SampleInfo$Sample,pattern = paste0(indivName[1],"_","2"))]
      }else{
        otherhap<-SampleInfo$HapID[grep(x = SampleInfo$Sample,pattern = paste0(indivName[1],"_","1"))]
      }
      if(length(otherhap)<1){
        otherhap<-NA
      }
      otherhap
    })
    
    SampleInfo$IndivType<-mapply(x=SampleInfo$Genotype,y=SampleInfo$OtherGenotype,function(x,y){if(is.na(y)){"Hemizygous"}else if(x==y & x!="SI"){"Homozygous"}else{"Heterozygous"}})

    
# 2. Estimate distances between sequences
    HapMatrix<-t(sapply(HaplotypeInfo$HapATGC,function(x){as.character(strsplit(x,split="")[[1]])}))
    rownames(HapMatrix)<-NULL
    
    distances<-as.dist(apply(HapMatrix,1,function(v){
      apply(HapMatrix,1,function(w){
        sum(v!=w)
      })
    }))
  
    # Change Haplotype Coding by MinorAllele 1 and MajorAllele 0
    HapMatrix01<-t(sapply(HaplotypeInfo$Hap01,function(x){as.character(strsplit(x,split="")[[1]])}))
    rownames(HapMatrix01)<-NULL

    HapMatrixMinMaj<-apply(HapMatrix01,2,function(v){
      Ref<-sum(v %in% "0")
      Alt<-sum(v %in% "1")
      if(Alt>Ref){
        vout<-gsub(x = v, pattern= "0", replacement = "Min")
        vout<-gsub(x = vout, pattern= "1", replacement = "Maj")
      }else{
        vout<-gsub(x = v, pattern= "0", replacement = "Maj")
        vout<-gsub(x = vout, pattern= "1", replacement = "Min")
      }
      vout
    })
    
# 3. Cluster the data based on the distances, divide them in subclusters
    Cluster<-hclust(distances)
    
    # SUBCLUSTER
    maxHeight<-max(Cluster$height)
    if(MinMutations<=0 | MinMutations>=maxHeight){
        write("Not a valid number of mutations! Setting it to 0.75*maxHeight\n",stdout())
        MinMutations <- 0.75*maxHeight
    }
    SubClusters<-cutree(Cluster,h = MinMutations)
    # SubCluster returns the original order.
    
    HaplotypeInfo$OriginalOrder<-as.integer(seq(1:nrow(HaplotypeInfo)))
    HaplotypeInfo$NewOrder[Cluster$order]<-1:nrow(HaplotypeInfo)
    HaplotypeInfo$SubCluster<-as.character(SubClusters)

    # Add subcluster info to the SampleInfo table
    SampleInfo$Subcluster<-sapply(SampleInfo$HapID,function(x){HaplotypeInfo$SubCluster[HaplotypeInfo$HapID %in% x]})
    
# 4. Find Subcluster outliers (S in I clusters, I in S clusters)
    # MainThreshold<-0.75
    SubClusterClassif<-data.frame()
    DiscordantSamples<-data.frame()
    # Cluster by cluster analyse main orientation and 
    for(i in c(unique(SubClusters))){
      
      m<-HaplotypeInfo[ HaplotypeInfo$SubCluster %in% i ,]
      
      # Main Orientation
      NumChrStd<-sum(m[,"S"])
      NumChrInv<-sum(m[,"I"])
      NumHapStd<-sum(m[,"S"]>0)
      NumHapInv<-sum(m[,"I"]>0)
      TotalNumChrs<-sum(c(NumChrStd,NumChrInv))
      MainOrientation<-if(TotalNumChrs==0){
        NA
      }else if(NumChrStd/TotalNumChrs>=MainThreshold){
        "S"
      }else if(NumChrInv/TotalNumChrs>=MainThreshold){
        "I"
      }else{
        "SI"
      }
      newdf<-data.frame(SubCluster=i,NumChrStd=NumChrStd,NumChrInv=NumChrInv,NumHapStd=NumHapStd,NumHapInv=NumHapInv,MainOrientation=MainOrientation)
      SubClusterClassif<-rbind(SubClusterClassif,newdf)
      
      # Discordant Haplotypes
      if(MainOrientation %in% "SI"){
        # No discordant individuals, everything is discordant!
      }else{
        discordantHaps<-m$HapID[!m$Type %in% MainOrientation]
        for (hap in discordantHaps){
          discdf<-data.frame(SampleInfo[ SampleInfo$HapID %in% hap & !SampleInfo$Genotype %in% MainOrientation, c("Sample","Population","Genotype","HapID","OtherGenotype","OtherHapID","IndivType")])
          discdf$SubCluster<-i
          discdf$SubClusterMainOrientation<-MainOrientation
          DiscordantSamples<-rbind(DiscordantSamples,discdf)
        }
      }
    }
    
    # Find the subcluster of the other chromosome in discordant samples
    DiscordantSamples$OtherSubCluster<-sapply(DiscordantSamples$OtherHapID,function(x){ 
      if(is.na(x)){
        othersubcl<-NA
      }else{
        othersubcl<-HaplotypeInfo$SubCluster[HaplotypeInfo$HapID %in% x]
      }
      othersubcl
      })
    DiscordantSamples$OtherSubClusterMainOrientation<-sapply(DiscordantSamples$OtherSubCluster,function(x){
      if(is.na(x)){
        othermain<-NA
      }else{
        othermain<-as.character(SubClusterClassif$MainOrientation[SubClusterClassif$SubCluster %in% x])
      }
      othermain
      })
    
    # Check if discordant samples can be phasing error
    # I don't know how to check that automatically...

    # Add Population information on the Subcluster data frame
    SubClustPopCount<-t(sapply(split(SampleInfo,f=factor(SampleInfo$Subcluster,levels=SubClusterClassif$SubCluster)),function(x){table(factor(x$Population,levels=PopLevels),factor(x$Genotype,levels=c("S","I")))}))
    colnames(SubClustPopCount)<-paste(c(rep("S",length(PopLevels)),rep("I",length(PopLevels))),rep(PopLevels,2),sep="_")
    SubClusterClassif<-cbind(SubClusterClassif,SubClustPopCount)

    # Add Inversion event column of the least diverse orientation
    LeastDiverse<-if(sum(SubClusterClassif$NumHapStd) > sum(SubClusterClassif$NumHapInv)){"Inv"}else{"Std"}
    write(paste0("In inversion ",Inv," the least diverse orientation is ",LeastDiverse,".\nInversion events are going to be counted as the number of clusters where there is at least one ",LeastDiverse," chromosome"), stdout())

    SubClusterClassif$InversionEvent<-if(LeastDiverse == "Inv"){ifelse(SubClusterClassif$NumChrInv>0,"Std>Inv","None")}else{ifelse(SubClusterClassif$NumChrStd>0,"Inv>Std","None")}

    # Save them
    write.table(SubClusterClassif,file = paste0(Prefix,"SubClusterInfo.txt"),sep="\t",row.names=FALSE)
    write.table(DiscordantSamples,file = paste0(Prefix,"DiscordantSamples.txt"),sep="\t",row.names=FALSE)
    write.table(SampleInfo,file = paste0(Prefix,"SampleInfo.txt"),sep="\t",row.names=FALSE)
    write.table(HaplotypeInfo,file = paste0(Prefix,"HaplotypeInfo.txt"),sep="\t",row.names=FALSE)
    write.table(PositionInfo,file = paste0(Prefix,"PositionInfo.txt"),sep="\t",row.names=FALSE)
# 4. Represent the data
    # Dendrogram
    SubClusterRect<-as.data.frame(t(sapply(split(x=HaplotypeInfo[,c("HapID","SubCluster","NewOrder")],f=HaplotypeInfo$SubCluster),function(m)c(min(m[,"NewOrder"]),max(m[,"NewOrder"]),0,MinMutations,m[1,"SubCluster"]))),stringsAsFactors = FALSE)
    colnames(SubClusterRect)<-c("xmin","xmax","ymin","ymax","SubCluster")
    SubClusterRect$xmin<-as.integer(SubClusterRect$xmin)
    SubClusterRect$xmax<-as.integer(SubClusterRect$xmax)
    SubClusterRect$ymin<-as.integer(SubClusterRect$ymin)
    SubClusterRect$ymax<-as.integer(SubClusterRect$ymax)
    SubClusterRect$SubCluster<-factor(SubClusterRect$SubCluster,levels=as.character(1:nrow(SubClusterRect)))
    
    ClusterDendro<-dendro_data(Cluster)
    labelsdf<-label(ClusterDendro)
    labelsdf$S<-NA
    labelsdf$I<-NA
    labelsdf$S[HaplotypeInfo$NewOrder]<-HaplotypeInfo$S
    labelsdf$I[HaplotypeInfo$NewOrder]<-HaplotypeInfo$I
    if(length(PopLevels)>5){
      Populations<-brewer.pal(n=length(PopLevels),"Paired")[1:length(PopLevels)] 
    }else{
      Populations<-brewer.pal(n=10,"Paired")[seq(2,length(PopLevels)*2,2)]
    }
    names(Populations)<-PopLevels
    PopulationsPresence<-melt(HaplotypeInfo,id.vars=c("HapID","NewOrder"),measure.vars=PopLevels,variable.name="Population",value.name="Count")
    PopulationsPresence$PopOrder<-NA
    for (i in 1:length(PopLevels)) {
      PopulationsPresence$PopOrder[PopulationsPresence$Population %in% PopLevels[i]]<-i
    }

    gdendro<-ggplot() + 
      geom_rect(data=SubClusterRect, aes(xmin=xmin-0.5, ymin=ymin, xmax=xmax+0.5, ymax=ymax,fill=SubCluster))+
      geom_segment(data=segment(ClusterDendro), aes(x=x, y=y, xend=xend, yend=yend)) + 
      # geom_text(data=label(ClusterDendro), aes(x=x, y=y, label=paste0("Hap",label), hjust=0), size=3,show.legend = FALSE) +
      geom_point(data=PopulationsPresence[PopulationsPresence$Count>0,],aes(x=NewOrder,y=-2*maxHeight/8-(PopOrder*0.02-0.01)*maxHeight,color=Population))+
      scale_alpha(range=c(0.5,1))+
      scale_fill_brewer(type="qual",palette=1, name=NULL)+
      scale_color_manual(values=Populations,name=NULL)+
      coord_flip() + scale_y_reverse(expand = c(0.05, 0)) +
      theme_classic()+
      theme(axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            panel.background=element_rect(fill="white"),
            legend.background = element_rect(fill=NA),
            panel.grid=element_blank(),legend.position = c(0.1,0.8))+
      scale_x_continuous(expand = c(0,0))+
      ggtitle("Haplotype tree")
   
    if(sum(labelsdf$S>0)>0){
      gdendro<-gdendro+geom_segment(data=labelsdf[labelsdf$S>0,], aes(x=x, y=0, xend=x, yend=-maxHeight/8,alpha=S), size=0.5,show.legend = FALSE,arrow=arrow(type="closed",angle="10",length=unit(0.02,"npc")))

    }
    if(sum(labelsdf$I>0)>0){
      gdendro<-gdendro+geom_segment(data=labelsdf[labelsdf$I>0,], aes(x=x, y=-2*maxHeight/8, xend=x, yend=-maxHeight/8,alpha=I), size=0.5,show.legend = FALSE,arrow=arrow(type="closed",angle="10",length=unit(0.02,"npc")))
    }
          
    # Distances plot
    m.distances<-as.matrix(distances)
    df.distances<-melt(m.distances,varnames = c("HapA","HapB"),value.name = "Differences")
    df.distances$Comparison<-mapply(x=df.distances$HapA, y=df.distances$HapB, FUN = function(x,y){paste(c("Comparison",sort(c(x,y))),collapse="_")})
    df.distances$HapA<-factor(df.distances$HapA,levels=HaplotypeInfo$OriginalOrder[order(HaplotypeInfo$NewOrder)])
    df.distances$HapB<-factor(df.distances$HapB,levels=HaplotypeInfo$OriginalOrder[order(HaplotypeInfo$NewOrder)])
    df.distances<-df.distances[order(df.distances$HapA),]
    # df.distances<-df.distances[!duplicated(df.distances$Comparison),]
    
    gdist<-ggplot(df.distances)+
      geom_tile(aes(HapA,HapB,fill=Differences),show.legend=FALSE)+
      scale_fill_distiller(type="seq",palette = "Greys",guide=FALSE)+
      geom_rect(data=SubClusterRect,aes(xmin=xmin-0.5,xmax=xmax+0.5,ymin=xmin-0.5,ymax=xmax+0.5,color=SubCluster),alpha=0,size=1,show.legend = FALSE)+
      scale_color_brewer(type="qual",palette=1, name=NULL)+
      theme_classic()+
      theme(
        # axis.line.y = element_blank(),
            # axis.ticks.y = element_blank(),
            # axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
      )+
      ggtitle("Similarity matrix")
    
    # Position matrix
    if(Mutations %in% "Minor"){
        HapdataframeMinMaj<-melt(HapMatrixMinMaj,varnames=c("Haplotype","Position"))
        HapdataframeMinMaj$Haplotype<-factor(HapdataframeMinMaj$Haplotype,levels=Cluster$order,labels=HaplotypeInfo$HapID[Cluster$order])
        HapdataframeMinMaj$Position<-factor(HapdataframeMinMaj$Position,levels=as.character(seq(1:nrow(PositionInfo))),labels=PositionInfo$ID)
        
        ghaplo<-ggplot(HapdataframeMinMaj)+geom_tile(aes(Position,Haplotype,fill=value),color="darkgray",show.legend = FALSE)+
          scale_fill_manual(values=c("white","black"),labels=c("Maj","Min"),name="Alleles")+ylab(label = NULL)+
          geom_rect(data=SubClusterRect,aes(xmin=0.5,xmax=nrow(PositionInfo)+0.5,ymin=xmin-0.45,ymax=xmax+0.45,color=SubCluster),alpha=0,size=1,show.legend = FALSE)+
          scale_color_brewer(type="qual",palette=1, name=NULL)+
          theme_classic()+
          theme(axis.line.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.y=element_blank(),
                axis.title.y=element_blank(),
                axis.text.x = element_text(angle = 90,vjust = 0.5),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank())+
          ggtitle("Haplotypes (minor allele)")
    }else{
        Hapdataframe<-melt(HapMatrix,varnames=c("Haplotype","Position"))
        Hapdataframe$Haplotype<-factor(Hapdataframe$Haplotype,levels=Cluster$order,labels=HaplotypeInfo$HapID[Cluster$order])
        Hapdataframe$Position<-factor(Hapdataframe$Position,levels=as.character(seq(1:nrow(PositionInfo))),labels=PositionInfo$ID)
        Hapdataframe$value<-factor(Hapdataframe$value,levels=c("A","C","G","T"))
        ghaplo<-ggplot(Hapdataframe)+geom_tile(aes(Position,Haplotype,fill=value))+
          scale_fill_manual(values=c("white","gray20","gray80","black"),labels=c("Anc","Alt","Ref","Der"),name="Alleles",drop=FALSE)+ylab(label = NULL)+
          geom_rect(data=SubClusterRect,aes(xmin=0.5,xmax=nrow(PositionInfo)+0.5,ymin=xmin-0.45,ymax=xmax+0.45,color=SubCluster),alpha=0,size=1,show.legend = FALSE)+
          scale_color_brewer(type="qual",palette=1, name=NULL)+
          theme_classic()+
          theme(axis.line.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.y=element_blank(),
                axis.title.y=element_blank(),
                axis.text.x = element_text(angle = 90,vjust = 0.5),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank())+
          ggtitle("Haplotypes")
    }

    
    # Mixed plot
    
    png(filename = paste0(Prefix,"ClusterPlot.png"),width=1800,height=1000)
    # cowplot
    # plot_grid(gdendro,gSI,ghaplo,rel_widths = c(0.3,0.05,0.65),align = "h",nrow = 1)
    mainplot<-plot_grid(gdendro,gdist,ghaplo,rel_widths = c(0.2,0.3,0.5),align = "h",nrow = 1)
    
    title <- ggdraw() + draw_label(paste0(Inv,"\n distance between clusers is at least ",MinMutations," differences on average"), fontface='bold')
    plot_grid(title, mainplot, ncol=1, rel_heights=c(0.05, 1)) # rel_heights values control title margins

    dev.off()
    


# 5. Save fasta file in the order of the cluster
cmd<-paste0("rm ",paste0(Prefix,"Haplotypes.fasta"))
system(cmd)
for (i in nrow(HaplotypeInfo):1) {
  j<-HaplotypeInfo$OriginalOrder[HaplotypeInfo$NewOrder == i]
  write(paste0(">",HaplotypeInfo[j,"HapID"],"_S",HaplotypeInfo[j,"S"],"_I",HaplotypeInfo[j,"I"]),paste0(Prefix,"Haplotypes.fasta"),append=TRUE)
  write(HaplotypeInfo[j,"HapATGC"],paste0(Prefix,"Haplotypes.fasta"),append=TRUE)
}


      
q()

