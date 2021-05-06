#!/bin/bash
# Ruth Gómez Graciani
# 01 04 2020

###############################################################################
# Description:                                                                
# Run all the analyses in this project                                        
###############################################################################

# SAVE HISTORY
# Save date and command 
# =========================================================================== #
HISTORY="$@"

# USAGE 
# Help message. Don't forget to mention new options in README.md!!  
# =========================================================================== #

usage()
{
  echo "Usage: 

  $(basename $0) <COMMAND> [OPTIONS]
  "
  echo "Commands:

  setup                   Download publicly available data and required programs.
  phase                   Phase inversions with MVNcall. 
  ihplot                  Make iHPlots from phased inversions. 
  "
  echo "Options:

    -h                    Show help.
  "

}

# PARSE OPTIONS
# Parse options and get help 
# =========================================================================== #

# Parse help message
while getopts "h" OPT; do
   case "" in
      h)  usage; exit 1 ;;
   esac
done
shift $((OPTIND -1))

# Save command, if any
COMMAND=$1; shift

# Set default optional variables
#SCREENS=1;
STEP=00           # Steps

# Set empty mandatory variables

# Parse command optons
case "$COMMAND" in
  #Prepare raw data
  preprocess ) 
    while getopts "s:" OPTIONS ; do
      case "$OPTIONS" in
        s)  SCREENS=${OPTARG} ;;
      esac
    done
    shift $((OPTIND -1))
    ;;
esac

# # Check that empty mandatory variables are full
# if [ "$COMMAND" == "all"]; then
#   if [-z "${VAR}"] ||[] ; then ; ;fi
# elif []; then
# fi

# SAVE HISTORY 
# Save date and command 
# =========================================================================== #
DATE=$(date +%F)

echo "${DATE}: $0 ${HISTORY}" >> project/history.txt

# SAVE LOG 
# Save date and command 
# =========================================================================== #
  exec 3>&1 4>&2
  trap 'exec 2>&4 1>&3' 0 1 2 3
  exec 1>project/logfiles/${DATE} 2>&1

# SETUP PROGRAM
# All publicly available data is now downloaded to data/raw and programs to software/.
# =========================================================================== #

if [ "$COMMAND" == "setup" ]; then
  
  # SETUP PROGRAM - DOWNLOAD PHASING SOFTWARE: MVNCall
  # --------------------------------------------------------------------------- #

  # Download
  mkdir -p code/software/mvncall/
  wget https://mathgen.stats.ox.ac.uk/genetics_software/mvncall/mvncall_v1.0_x86_64_dynamic.tgz -P code/software/mvncall/

  # Install
  cd code/software/mvncall/ || return
  tar -zxvf mvncall_v1.0_x86_64_dynamic.tgz 
  rm mvncall_v1.0_x86_64_dynamic.tgz 
  # Program is in mvncall_v1.0_x86_64_dynamic/mvncall

  # SETUP PROGRAM - DOWNLOAD 1000 GENOMES PROJECT DATA
  # --------------------------------------------------------------------------- #

  # MASK  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/StrictMask/20140520.allAutosome.strict_mask.fasta.gz
  # MAKE USABLE VERSION OF MASK
  # GENOMES
  cd ../../../data/raw/|| return
  mkdir 1KGP
  cd 1KGP
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*
  cd ../../../

  # MASK
  cd ../../../data/raw/|| return
  mkdir Accessibility
  cd Accessibility
  # Ph3 - pilot
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.pilot_mask.whole_genome.bed
  # Ph3 - strict
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20141020.strict_mask.whole_genome.bed
  cd ../../../


fi

# PHASE INVERSION
# Phased data is required to make iHPlots
# =========================================================================== #
STEP=$(printf "%02d" $((${STEP}+1)))

if [ "$COMMAND" == "phase" ]; then

  # PHASE INVERSION - Set directories and variables
  # --------------------------------------------------------------------------- #
  # INDIR=$(echo $INDIR| grep -o '[^/]\+$')
  TMPDIR="tmp/${DATE}_${STEP}_phase"
  OUTDIR="analysis/${DATE}_${STEP}_phase"
  mkdir -p "$TMPDIR" "$OUTDIR"

# OPTION 1
 GENOTYPES="data/use/inversions_info/Genotypes_20062018.txt"
#  OPTION 2 
  # GENOTYPES="data/use/inversions_info/LongGenotypes.txt"
  INVCOORD="data/use/inversions_info/Coordinates_hg19.txt"
  GENOMES="data/raw/1KGP"
  MASK="data/use/Accessibility/20141020.strict_mask.whole_genome.bed"

  # LIST SAMPLES

    SAMPLES_REF=$(cut -f1 ${GENOMES}/integrated_call_samples_v3.20130502.ALL.panel | tail -n+2)
    SAMPLES_REF_MALE=$(cut -f1 ${GENOMES}/integrated_call_male_samples_v3.20130502.ALL.panel | tail -n+2)

  # LIST INVERSIONS - I should think how I need coordinates
  # IF THERE IS A HEADER:  | tail -n+2)
    INVERSIONS=$(cut -f1 $INVCOORD ) 

  # FOR EACH INVERSION
    for INV in $INVERSIONS; do
        # OPTION 1
    SAMPLES_INV=$(cut -f2 $GENOTYPES | tail -n+2)
    # OPTION 2  
    # SAMPLES_INV=$(grep $INV $GENOTYPES | cut -f1  | tail -n+2)
    # OPTION 3
    SAMPLES_FILTER=$(cut -f2  data/raw/CarlaExamples/Samples  | tail -n+2)
    SAMPLES_INV=$(echo "$SAMPLES_INV" "$SAMPLES_FILTER"  | tr " " "\n" | sort | uniq -d)
   

      # PHASE INVERSION - Make MVNcall input
      # --------------------------------------------------------------------------- #
      # MAKE FOLDER
        MVN_INPUT="${TMPDIR}/${INV}"
        MVN_OUTPUT="${OUTDIR}/${INV}"
        mkdir -p "$MVN_INPUT" "$MVN_OUTPUT"

      # MAKE OUTPUT FILE NAMES
        SAMPLE="${MVN_INPUT}/sample-file.txt"
        SCAFFOLD="${MVN_INPUT}/scaffold-file.txt"
        GLFS="${MVN_INPUT}/gffs-file.txt"
      
      # TAKE COORDINATES (general info)
        # (SD1 end +1; SD2 start -1)
        # I assume this format: HsInv0389	X	153564286	153575614	153613228	153624563
        CHR=$(grep -e "$INV" $INVCOORD | cut -f2)
        START=$(($(grep -e "${INV}" ${INVCOORD} | cut -f4) + 1))
        END=$(($(grep -e "${INV}" ${INVCOORD} | cut -f5) - 1))
    
      # ALSO PREPARE INPUT FOR A POSSIBLE FUTURE IHPLOT
        FLANKING=2000
        OUT_START=$(($(grep -e "${INV}" ${INVCOORD} | cut -f3) - ${FLANKING}))
        OUT_END=$(($(grep -e "${INV}" ${INVCOORD} | cut -f6) + ${FLANKING}))

      # 1000 GENOMES DATA: Genotypes and shared samples
        #  Depending on chromosome, take bcf file name and shared samples list
        if [[ "$CHR" == "Y" ]]; then
          BCF_FILE="${GENOMES}/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz"
          SAMPLES_SHARED=$(echo "$SAMPLES_INV" "$SAMPLES_REF_MALE" | tr " " "\n" | sort | uniq -d)
        elif  [[ "$CHR" == "X" ]]; then 
          BCF_FILE="${GENOMES}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
          SAMPLES_SHARED=$(echo "$SAMPLES_INV" "$SAMPLES_REF"  | tr " " "\n" | sort | uniq -d)
        else
          BCF_FILE="${GENOMES}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
          SAMPLES_SHARED=$(echo "$SAMPLES_INV" "$SAMPLES_REF"  | tr " " "\n" | sort | uniq -d)
        fi
   
      # TAKE 1000 GENOMES DATA FROM INDIVIDUALS AND REGION OF INTEREST
        # --min-ac 1:minor (remove monomorphic sites) -m2 -M2 (keep only biallelic) -v snps (only SNP variants)
        # OJO -- en el codigo de Carla filtraba aquellos que solo tuviesen 1 alelo en el AC, es decir, 1 heterocigoto entre todos los invs. Por eso yo he puesto --min-ac 2:minor
        bcftools view -r${CHR}:${START}-${END} -Ov -s $(echo $SAMPLES_SHARED | tr " " ",") --min-ac 2:minor -M2 -m2 -v snps "${BCF_FILE}" > ${MVN_INPUT}/1kgp_filter.txt 

        # Apply STRICT accessibility mask (or any mask, actually, changing the file)
        bedtools intersect -wa -header -a ${MVN_INPUT}/1kgp_filter.txt -b $MASK |\
          # Also make cusom format
        bcftools query -H -f'%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%GT]\n'   > ${MVN_INPUT}/1kgp_filter_mask_format.txt 

       # ---------------------CODE IN TESTING PHASE ---------------------------
        # TO LIMIT VARIANTS ACCORDING TO EXAMPLE:
        mv ${MVN_INPUT}/1kgp_filter_mask_format.txt  ${MVN_INPUT}/1kgp_filter_mask_format_unfiltered.txt 
        cut -f1 data/raw/CarlaExamples/Positions > "${MVN_INPUT}/posfilter.txt"
        grep -f "${MVN_INPUT}/posfilter.txt" ${MVN_INPUT}/1kgp_filter_mask_format_unfiltered.txt  > ${MVN_INPUT}/1kgp_filter_mask_format.txt

      # -----------------------------------------------------------

        # ALSO PREPARE INPUT FOR A POSSIBLE FUTURE IHPLOT - outside region, unfiltered variants
        TMP=$(echo $SAMPLES_SHARED | tr " " ",")
        echo "CHROM	POS	ID	REF	ALT	VT	AA	${TMP}" | tr "," "	" > ${MVN_OUTPUT}/1KGP_haplotypes.txt
        bcftools view -r${CHR}:${OUT_START}-${OUT_END} -Ov -v snps -s $(echo $SAMPLES_SHARED | tr " " ",")  "${BCF_FILE}" |\
          bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%VT\t%AA[\t%GT]\n' >> ${MVN_OUTPUT}/1KGP_haplotypes.txt 
        


      # -----------------------------------------------------------
      #  FUTURE OPTIONS 
      #  This is a basic version of the program. At this stage, there can be a checkpoint in which I check that there are between 50 and 200 SNPs to study
      #  If it is too small, I can apply Pilot accessibility mask or no mask
      #  If it is too large, I can use only OMNI SNPs
      # 
      #  Also, if there are inversions with INVERTED REPEATS inside the region, I should exclude those SNPs inside the IR. 
      #  This can be easily filtered with bedtools intersect and grep -v -f (file with POSITIONS! not IDs.)
      # -------------------------------------------

      # MAKE SAMPLE FILE - This is a file with the list of individuals
        # head -n 1 "${MVN_INPUT}/1kgp_filter_mask_format.txt" | cut -f6- | sed 's/:GT//g' | sed  's/\[[0-9]*\]//g'|sed 's/\t/\n/g' > $SAMPLE
        echo $SAMPLES_SHARED | tr " " "\n"> ${SAMPLE}

      # MAKE SCAFFOLD FILE - This is a file with reference haplotypes
        tail -n+2 "${MVN_INPUT}/1kgp_filter_mask_format.txt"|\
          # Transformation for X and Y chromosomes 
          awk -v OFS=" " '{for(i=1; i<=NF; i++) {if($i=="0") $i = "0|0"}}{for(i=1; i<=NF; i++) {if($i=="1") $i = "1|1"}}{print$0}' |\
          # Final touches for .haps format
          tr "|" " " | tr "\t" " "> $SCAFFOLD
  
      #  MAKE GLFS FILE - This is a VCF with the inversion. Position of inversion must be in the middle. 
        
        # Take header from random file - THIS STARTS GLFS FILE
        bcftools view -h -s $(echo $SAMPLES_SHARED | tr " " ",") "${GENOMES}/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" > $GLFS
        
        # Take inversion genotypes from genotypes file only for common individuals
        # OPTION 1
          #  Which column is the inversion in
          INV_COLUMN=$(head -1 $GENOTYPES | tr "\t" "\n" | cat -n | grep $INV | cut -f1 | tr -d " ")
          # Grep only common individuals, stored in $SAMPLE, format as in VCF
          cut -f2,$INV_COLUMN $GENOTYPES | grep -f $SAMPLE |\
            #  Standard, het and inverted are normal
            sed 's/Std\/Std/0\/0:0.00,-5.00,-5.00/g' | sed 's/Std\/Inv/0\/1:-5.00,0.00,-5.00/g' | sed 's/Inv\/Inv/1\/1:-5.00,-5.00,0.00/g' | \
            # Missing genotypes
            sed 's/ND/.\/.:-0.48,-0.48,-0.48/g' |\
            #  DEL genotypes
            sed 's/Inv\/Del/.\/1:-2.70,-0.30,-0.30/g' | \
            # Male genotypes in X format
            sed 's/Std/0\/0:0.00,-5.00,-5.00/g' | sed 's/Inv/1\/1:-5.00,-5.00,0.00/g' > ${MVN_INPUT}/genotypes_selected_vcformatted.txt
        # OPTION 2
          #  grep -f $SAMPLE $GENOTYPES | grep $INV | cut -f 1,7 | \
          #  # Male genotypes in X format
          #   sed 's/\t0$/\t0\/0/g' | sed 's/\t1$/\t1\/1/g'| \
          #  # Add genotype likelihoods 
          #   sed 's/0\/0/0\/0:0.00,-5.00,-5.00/g' | sed 's/1\/1/1\/1:-5.00,-5.00,0.00/g' >  ${MVN_INPUT}/genotypes_selected_vcformatted.txt

        # Select inversion position: inv pos should be a position between the two central SNPs, because MVNCall will use the same number of SNPs at either side
          # Set number of snps (useful for later)
          NUMSNPS=$(( $(cat $SCAFFOLD  | wc -l  )/2))         
          # Center position +1
          POSITION=$(( $(head -$NUMSNPS $SCAFFOLD | tail -1 | cut -d' ' -f3) + 1))
          # if pCenter position +1 exists Center position +2
          if [[ "$(grep $POSITION $SCAFFOLD)" != "" ]]; then
              POSITION=$(($POSITION+1))
          fi
        
        # Get genotypes sorted the same as in other files
         GENOTYPES_ORDERED=$(for SAMPLE in $SAMPLES_SHARED; do grep "$SAMPLE" ${MVN_INPUT}/genotypes_selected_vcformatted.txt | cut -f2; done)

        # Fill VCF         

        echo -e "${CHR}\t${POSITION}\t${INV}\tA\tT\t.\tPASS\tSTART=${START};END=${END}\tGT:GL\t$(echo ${GENOTYPES_ORDERED} | tr ' ' '\t')" >>${GLFS}

      # PHASE INVERSION - Make MVNcall output
      # --------------------------------------------------------------------------- #
      

      # EXECUTE 
        ./code/software/mvncall/mvncall_v1.0_x86_64_dynamic/mvncall \
          --sample-file $SAMPLE \
          --glfs $GLFS \
          --scaffold-file $SCAFFOLD \
          --o "${MVN_OUTPUT}/${INV}_phased.vcf" \
          --int $POSITION $POSITION \
          --numsnps $NUMSNPS
    done
            

fi


# MAKE NETWORK ANALYSES
# This code will finally create the iHPlots
# =========================================================================== #
STEP=$(printf "%02d" $((${STEP}+1)))

if [ "$COMMAND" == "ihplot" ]; then

  # MAKE NETWORK ANALYSES - Set directories and variables
  # --------------------------------------------------------------------------- #

  # SET DIRECTORIES
    TMPDIR="tmp/${DATE}_${STEP}_ihplot"
    OUTDIR="analysis/${DATE}_${STEP}_ihplot"
    mkdir -p "$TMPDIR" "$OUTDIR"
    INVCOORD="data/use/inversions_info/Coordinates_hg19.txt"
    
  # SET ANALYSIS VARIABLES
    MINCOUNT="2" # At least two counts
    FLANKING="0" # Only inside, no flanking region
    MAINTHRESHOLD="0.6" # Decide the cluster type if one orientation is at 60% or more
    MUTATIONS="Derived" # Show derived mutations in the plot when possible
    MINMUTKB="1.5" # 1.5/kbp is the minimum number of differences (with >MinCount) between sequences to be considered different clusters

  # START LOOP
  INVERSIONS=$(cut -f1 $INVCOORD ) 

  # FOR EACH INVERSION
  for INV in $INVERSIONS; do
  
    # TAKE GENERAL INFO FROM COORDINATES FILE 
    # (SD1 end +1; SD2 start -1)
    # I assume this format: HsInv0389	X	153564286	153575614	153613228	153624563
      CHR=$(grep -e "$INV" $INVCOORD | cut -f2)
      START=$(($(grep -e "${INV}" ${INVCOORD} | cut -f3) + 1))
      END=$(($(grep -e "${INV}" ${INVCOORD} | cut -f6) - 1))

    # MAKE NETWORK ANALYSES - Include phased breakpoints into phased 1KGP info
    # --------------------------------------------------------------------------- #
    

    # SET FOLDERS
      INDIR="analysis/2020-11-29_01_phase/${INV}"
      STEPDIR="${TMPDIR}/00_Complete/${INV}"
      mkdir -p ${STEPDIR}

    # ANCESTRAL ORIENTATION (Optional)
      ANCORIENT=$(cat data/use/inversions_info/BP_SD_positions_45inv_RevisedHG18HG19_v3.1.txt | grep $INV | cut -f19) 
      if [ "$ANCORIENT" == "" ]; then
        ANCORIENT="NA"
      fi

    # SCRIPT TO INCLUDE PHASING AND BEAUTIFY
      Rscript code/rscript/TidyHaplotypes.R "${CHR}:${START}" "${CHR}:${END}" ${INDIR}/${INV}_phased.vcf ${INDIR}/1KGP_haplotypes.txt ${ANCORIENT} "${STEPDIR}/"

    # ---------------------CODE IN TESTING PHASE ---------------------------
    cut -f1 data/raw/CarlaExamples/Positions > "${STEPDIR}/posfilter.txt"
    cut -f1 data/raw/CarlaExamples/Samples > "${STEPDIR}/samfilter.txt"

    mv ${STEPDIR}/samples.txt ${STEPDIR}/samples_unfiltered.txt 
    mv ${STEPDIR}/haplotypes.txt ${STEPDIR}/haplotypes_unfiltered.txt 
    mv ${STEPDIR}/positions.txt ${STEPDIR}/positions_unfiltered.txt 

    head -1 ${STEPDIR}/samples_unfiltered.txt >  ${STEPDIR}/samples.txt
    head -1 ${STEPDIR}/haplotypes_unfiltered.txt >  ${STEPDIR}/haplotypes.txt
    head -1 ${STEPDIR}/positions_unfiltered.txt >  ${STEPDIR}/positions.txt

    grep -f ${STEPDIR}/posfilter.txt ${STEPDIR}/positions_unfiltered.txt > ${STEPDIR}/positions.txt
    grep -f ${STEPDIR}/samfilter.txt ${STEPDIR}/haplotypes_unfiltered.txt >> ${STEPDIR}/haplotypes.txt
    grep -f ${STEPDIR}/samfilter.txt ${STEPDIR}/samples_unfiltered.txt > ${STEPDIR}/samples.txt

    # ------------------------------------------------


    # MAKE NETWORK ANALYSES - Remove singletons and trim to include only inside region. Make other input files
    # --------------------------------------------------------------------------- #

    # SET FOLDERS
      INDIR=${STEPDIR}
      STEPDIR="${TMPDIR}/01_Region/${INV}"
      mkdir -p "${STEPDIR}"

    # SCRIPT TO TRIM 
      Rscript code/rscript/SelectRegion.R ${INDIR}/positions.txt ${INDIR}/haplotypes.txt ${INDIR}/samples.txt ${MINCOUNT} "${STEPDIR}/" ${FLANKING} 

    # MAKE NETWORK ANALYSES - Make cluster analysis plots!
    # --------------------------------------------------------------------------- #
    # SET FOLDERS
      INDIR=${STEPDIR}
      STEPDIR="${TMPDIR}/02_iHPlots/${INV}"
      mkdir -p ${STEPDIR}

    # CALCULATE MINIMUN  MUTATIONS
      SIZE=$(($END-$START+1+$FLANKING+$FLANKING))
      MINMUTATIONS="${MINMUTKB}*${SIZE}/1000"

    # MAKE PLOTS
      Rscript code/rscript/AnalyseClusters.R ${INDIR}/PositionInfo.txt ${INDIR}/HaplotypeInfo.txt ${INDIR}/SampleInfo.txt ${MAINTHRESHOLD} ${STEPDIR}/ ${MUTATIONS} ${MINMUTATIONS} 

  done 

fi