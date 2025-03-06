# Initialisierung
#install.packages("stringr")
#install.packages("dplyr")
#install.packages("openxlsx", dependencies = TRUE)
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")


library(stringr)
library(dplyr)
library(openxlsx)
library(Biostrings)


# Configuration
importFolder <- ""
exportFolder <- ""
importFileName <- ""
importFileType <- ".txt"
exportExcelFileName <- ""

# Importing Library Configuration
importLibTable <- read.delim(paste(importFolder,"/",importFileName,importFileType,sep=""))


# Function for generating AminoAcid Mutation Info
getMutationStrings <- function(mutationCigar, refSeq, offsetCigar2Gen, offsetCigar2RefSeq) {
  
  #validate cigarstring
  if(!grepl("^\\d+(C|A|T|G)$",mutationCigar )) return("ERROR_INVALIDCIGARSTRING")
  
  #split cigarstring in position and nucleotide info
  NT_Mut <- str_extract(mutationCigar,"(C|A|T|G)$")
  mutationNtPosCigar <- str_extract(mutationCigar,"^\\d+")
  
  #calculate position of mutated nucleotide in gene seq
  mutationNtPosGen <- as.integer(offsetCigar2Gen) + as.integer(mutationNtPosCigar)
  mutationNtPosRefSeq <- as.integer(offsetCigar2RefSeq) + as.integer(mutationNtPosCigar)
  
  #extract reference nucleotide
  NT_Ref <- str_sub(refSeq,mutationNtPosRefSeq, mutationNtPosRefSeq)
  
  #calculate relative mutation position in AA seq (0 / 1 / 2)
  mutPosInAA <- (mutationNtPosGen+2)%%3 
  
  #extract AA nt triple out of refseq
  AA_NTSeqRef <- str_sub(refSeq, mutationNtPosRefSeq-mutPosInAA, mutationNtPosRefSeq-mutPosInAA+2)
  
  #create mutated AA nt triple
  AA_NTSeqMut <- AA_NTSeqRef
  str_sub(AA_NTSeqMut,mutPosInAA+1,mutPosInAA+1) <- NT_Mut
  
  #translate AA nt seq in AA
  AA_Ref <-  translate(DNAString(AA_NTSeqRef), no.init.codon = TRUE)
  AA_Mut <-  translate(DNAString(AA_NTSeqMut), no.init.codon = TRUE)
  
  #calculate number of amino acid in gene seq
  AANumber <- ceiling(mutationNtPosGen/3)
  
  #create result strings
  codonString <- paste(NT_Ref,mutationNtPosGen,NT_Mut,sep="")
  AAString <- paste(AA_Ref,AANumber,AA_Mut,sep="")
  
  return(c(codonString, AAString, AA_NTSeqRef, AA_NTSeqMut))
} 
# getMutationStrings("106A",refSeq1,geneOffset1,-10)

# Generation of Libraries
for (libNo in c(1:nrow(importLibTable)) ) {

  cdssequence <- importLibTable[libNo,"Sequence"]
  adapter_5 <- importLibTable[libNo,"Adapter_5"]
  adapter_3 <- importLibTable[libNo,"Adapter_3"]
  headcrop_nt_count <- importLibTable[libNo,"crop_5"]
  tailcrop_nt_count <- importLibTable[libNo,"crop_3"]
  NtOffset <- importLibTable[libNo,"NtOffset"]
  
  # creating new result table
    result <- data.frame()
  
  # writing WT Line in result table
    result[1,1] <- cdssequence
    result[1,2] <- paste(adapter_5,cdssequence,adapter_3,sep="")
    result[1,3] <- paste(str_sub(adapter_5,-10),cdssequence,str_sub(adapter_3,1,10),sep="")
    result[1,4] <- "WT"
    
  for (i in c(1:str_length(cdssequence)) ) {
    for (NT in c("A","T","C","G")) {
      if(substr(result[1,1],i,i) != NT) { 
        tempseq <- cdssequence
        substr(tempseq,i,i) <- NT
        result[nrow(result) + 1,1] <- tempseq
        result[nrow(result),2] <- paste(adapter_5,tempseq,adapter_3,sep="")
        result[nrow(result),3] <- paste(str_sub(adapter_5,-10),tempseq,str_sub(adapter_3,1,10),sep="")
        result[nrow(result),4] <- paste(i+str_length(adapter_5)-headcrop_nt_count,NT,sep="")
        
        mutationInfo <- getMutationStrings(result[nrow(result),4],cdssequence,NtOffset,headcrop_nt_count-str_length(adapter_5))
        result[nrow(result),5] <- mutationInfo[1]
        result[nrow(result),6] <- mutationInfo[2]
        result[nrow(result),7] <- mutationInfo[3]
        result[nrow(result),8] <- mutationInfo[4]
        
      }  
    }
  }
    
  assign(paste("result_",importLibTable[libNo,"Name"],sep=""),result)

}




# Prepare Exporting
print("Prepare Exporting")
if(!dir.exists(exportFolder)) {
  dir.create(exportFolder)  
  print('Creating Export Folder')
}

wb <- createWorkbook()

for (libNo in c(1:nrow(importLibTable)) ) {
  addWorksheet(wb, importLibTable[libNo,"Name"])
  writeData(wb, importLibTable[libNo,"Name"], get(paste("result_",importLibTable[libNo,"Name"],sep="")),colNames = FALSE)
}
saveWorkbook(wb, file.path(exportFolder, exportExcelFileName), overwrite = TRUE)



