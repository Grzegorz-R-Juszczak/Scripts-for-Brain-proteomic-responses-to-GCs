####################################################################
####### Ensembl data import and preprocessing ##########################
####################################################################

Ensembl_data = read.delim('Final_search_results.csv', sep=',', header = TRUE) # data import
head(Ensembl_data, n = 8) # showing first n rows

Ensembl_data$Ensembl_NCBI_ID_assigment <- ifelse(grepl(",", Ensembl_data$NCBI_ID),'Ambiguous','Unique') # testing for the presence of “,” in items in the column NCBI_ID. Presence of the “,” indicates that more than one gene is assigned to the ID.
head(Ensembl_data, n=10) # showing first n rows

####################################################################
####### import Ensembl / UniProt ID list ###################################
####################################################################

Ensembl_UniProt_final = read.delim('Ensembl_UniProt_final.csv', sep=';', header = TRUE) # data import
head(Ensembl_UniProt_final, n = 8) # showing first n rows

Ensembl_UniProt_final_2 = tidyr:: drop_na(Ensembl_UniProt_final, Ensembl_ID) # usuwanie wierszy z brakującymi wartościami w kolumnie Gene.name
head(Ensembl_UniProt_final_2, n = 10) # showing first n rows

####################################################################
####### Joining the Ensembl and Ensembl_UniProt_final data and editing #########
####################################################################

Ens_EnsUni  = dplyr:: full_join(Ensembl_data, Ensembl_UniProt_final_2, by = c("Ensembl_ID" = "Ensembl_ID"), keep = FALSE, na_matches = "never") # joining the data
head(Ens_EnsUni, n=10) # showing first n rows

library(tidyr)
Ens_EnsUni_2 = Ens_EnsUni %>% separate_longer_delim(NCBI_ID, delim = ', ') # splits by a delimiter into separate rows
head(Ens_EnsUni_2, n=10) # showing first n rows

colnames(Ens_EnsUni_2) <- c("Anchored_input_gene", "Input_gene_symbols", "Anchored_output_gene_Ensembl", "Ensembl_official_symbols", "Ensembl_ambiguous_unique", "Ensembl_gene_description", "Ensembl_ID", "NCBI_ID_in_Ensembl", "MGI_ID_in_Ensembl", "Ensembl_gene_type", "UniProtSwissProt_in Ensembl", "Ensembl_NCBI_ID_assigment", "UniProtFinal_in Ensembl") # setting column names
head(Ens_EnsUni_2, n=8) # showing first n rows


####################################################################
####### UniProt data import and pre-processing ############################
####################################################################
UniProt_data <-readr::read_tsv("UniProt_data.tsv", col_names = TRUE) # data import
head(UniProt_data, n = 8) # showing first n rows

UniProt_data <- as.data.frame(UniProt_data) # changing data format to data frame. 
head(UniProt_data, n = 8) # showing first n rows

UniProt_data <- UniProt_data [, c("From", "Entry", "Entry Name", "Protein names", "Gene Names (primary)", "GeneID", "MGI", "Organism", "Organism (ID)")] # setting column order
head(UniProt_data, n = 8) # showing first n rows

colnames(UniProt_data) <- c("From", "Entry_UniProt", "Entry_Name_UniProt", "Protein_names_UniProt", "Gene_Names_UniProt",  "NCBI_Gene_ID_UniProt", "MGI_ID_UniProt", "Organism_UniProt",  "Organism_ID_UniProt") # setting column names
head(UniProt_data, n=8) # showing first n rows

UniProt_data2 = dplyr::mutate(UniProt_data, NCBI_Gene_ID_UniProt_clean  = substr(NCBI_Gene_ID_UniProt, 1, nchar(NCBI_Gene_ID_UniProt) - 1)) # removing “;” located at the end of NCBI IDs
head(UniProt_data2, n = 5) # showing first n rows

UniProt_data3 = dplyr::mutate(UniProt_data2, MGI_ID_UniProt_clean  = substr(MGI_ID_UniProt, 1, nchar(MGI_ID_UniProt) - 1)) # removing “;” located at the end of the IDs
head(UniProt_data3, n = 5) # showing first n rows

UniProt_data3$NCBI_GeneID_UniProt_assigment <- ifelse(grepl(";", UniProt_data3$ NCBI_Gene_ID_UniProt_clean),'Ambiguous','Unique') # testing for the presence of “;” in items in the column GeneID_clean. Presence of the “;” indicates that more than one gene is assigned to the ID.
head(UniProt_data3, n=5) # showing first n rows

UniProt_data3$GeneName_UniProt_assigment <- ifelse(grepl(";", UniProt_data3$Gene_Names_UniProt),'Ambiguous','Unique') # testing for the presence of “;” in items in the column Gene_Names_UniProt. Presence of the “;” indicates that more than one gene is assigned to the ID.
head(UniProt_data3, n=10) # showing first n rows

Ambiguous <- dplyr::filter(UniProt_data3, NCBI_GeneID_UniProt_assigment == "Ambiguous" | GeneName_UniProt_assigment == "Ambiguous") # showing rows with ambiguous IDs
head(Ambiguous, n=10) # showing first n rows

#### Testing duplicated UniPror IDs ####

duplicates_from_first = duplicated(UniProt_data3$From, fromLast = FALSE) # identification of duplicated gene names (forward)
head(duplicates_from_first, n=4) # showing first n rows

duplicates_from_last = duplicated(UniProt_data3$From, fromLast = TRUE) # identification of duplicated gene names (reverse)
head(duplicates_from_first, n=4) # showing first n rows

UniProt_data4 <- cbind(UniProt_data3, duplicates_from_first, duplicates_from_last) # joining the dataset with results of duplication testing
head(UniProt_data4, n=4) # showing first n rows

UniProt_data4_not_duplicatted = dplyr::filter(UniProt_data4, duplicates_from_first == "FALSE" & duplicates_from_last == "FALSE") # selection of duplicated gene symbols
head(UniProt_data4_not_duplicatted, n=10) # showing first n rows

UniProt_data4_duplicatted = dplyr::filter(UniProt_data4, duplicates_from_first == "TRUE" | duplicates_from_last == "TRUE") # selection of duplicated gene symbols
head(UniProt_data4_duplicatted, n=10) # showing first n rows

#* Check data frame UniProt_data4_duplicatted for duplicated UniProt IDs in the “From” column.  If necessary, filter required species with the next command (silenced with #). Otherwise, go to the “UniProt_data5 = dplyr:: bind_rows” action. 

#UniProt_data4_duplicatted = dplyr::filter(UniProt_data4_duplicatted, Organism_UniProt == "Mus musculus (Mouse)") # selection of mouse identifiers
#head(UniProt_data4_duplicatted, n=10) # showing first n rows

UniProt_data5 = dplyr:: bind_rows(UniProt_data4_not_duplicatted, UniProt_data4_duplicatted) # binding rows
head(UniProt_data5, n=100) # showing first n rows

UniProt_data6 = dplyr::mutate(UniProt_data5,  duplicates_from_first = NULL, duplicates_from_last = NULL) # removing unnecessary columns
head(UniProt_data6, n=5) # showing first n rows

library(tidyr)
UniProt_data7 = UniProt_data6 %>% separate_longer_delim(NCBI_Gene_ID_UniProt_clean, delim = ';') # splits by a delimiter into separate rows
head(UniProt_data7, n=10) # showing first n rows

library(tidyr)
UniProt_data8 = UniProt_data7 %>% separate_wider_delim(Entry_Name_UniProt, delim = '_', names = c("Entry_Name2_UniProt", NA), cols_remove = FALSE) # splits by a delimiter into separate columns
head(UniProt_data8, n=10) # showing first n rows

UniProt_data8 <- as.data.frame(UniProt_data8) # changing data format to data frame. 
head(UniProt_data8, n = 8) # showing first n rows

####################################################################
####### Joining the Ensembl and UniProt data #############################
####################################################################

Ensembl_Uni  = dplyr:: full_join(Ens_EnsUni_2, UniProt_data8, by = c("UniProtFinal_in Ensembl" = "From"), keep = FALSE, na_matches = "never") # joining the data
head(Ensembl_Uni, n=10) # showing first n rows

# readr::write_csv(x = Ensembl_Uni, file = "Ensembl_Uni.csv") # data export

####################################################################
####### NCBI data import and pre-processing ############################
####################################################################

ncbi_dataset <-readr::read_tsv("ncbi_dataset.tsv", col_names = TRUE) # data import
head(ncbi_dataset, n = 8) # showing first n rows

ncbi_dataset <- as.data.frame(ncbi_dataset) # changing data format to data frame. 
head(ncbi_dataset, n = 8) # showing first n rows
ncbi_dataset <- ncbi_dataset [, c("NCBI GeneID", "Symbol", "Description", "Taxonomic Name", "Common Name", "Gene Type", "Transcripts", "Gene Group Identifier", "Gene Group Method", "Chromosomes", "Nomenclature ID", "Ensembl GeneIDs", "Annotation Genomic Range Accession", "Annotation Genomic Range Start", "Annotation Genomic Range Stop", "OMIM IDs", "Orientation", "Proteins", "Synonyms", "Taxonomic ID", "SwissProt Accessions")] # setting column order
head(ncbi_dataset, n = 8) # showing first n rows

colnames(ncbi_dataset) <- c("NCBI_GeneID", "Symbol_NCBI", "Description_NCBI", "Taxonomic_Name", "Common_Name", "NCBI_Gene_Type", "Transcripts", "Gene_Group_Identifier", "Gene_Group_Method", "Chromosomes", "Nomenclature_ID_NCBI", "Ensembl_GeneIDs", "Annotation_Genomic_Range_Accession", "Annotation_Genomic_Range_Start", "Annotation_Genomic_Range_Stop", "OMIM_IDs", "Orientation", "Proteins", "Synonyms", "Taxonomic_ID_NCBI", "SwissProt_Accessions_NCBI") # setting column names
head(ncbi_dataset, n = 8) # showing first n rows

ncbi_dataset2  = dplyr::mutate(ncbi_dataset, 
                               Taxonomic_Name = NULL, 
                               Common_Name = NULL,
                               Transcripts = NULL,
                               Gene_Group_Identifier = NULL,
                               Gene_Group_Method = NULL,
                               Chromosomes = NULL,
                               Ensembl_GeneIDs = NULL,
                               Annotation_Genomic_Range_Accession = NULL,
                               Annotation_Genomic_Range_Start = NULL,
                               Annotation_Genomic_Range_Stop = NULL,
                               OMIM_IDs = NULL,
                               Orientation = NULL,
                               Proteins = NULL,
                               Synonyms = NULL) # removing unnecessary columns
head(ncbi_dataset2, n = 8) # showing first n rows

ncbi_dataset2$NCBI_GeneID <-as.character(ncbi_dataset2$NCBI_GeneID) # changing data type in selected column
ncbi_dataset2$Nomenclature_ID_NCBI <-as.character(ncbi_dataset2$Nomenclature_ID_NCBI) # changing data type in selected column
head(ncbi_dataset2, n = 8) # showing first n rows

####################################################################
####### joining Ensembl, UniProt and NCBI data based on NCBI IDs ############
####################################################################

Ensembl_Uni_NCBI = dplyr:: left_join(Ensembl_Uni, ncbi_dataset2, by = c("NCBI_ID_in_Ensembl" = "NCBI_GeneID"), keep = TRUE, na_matches = "never") # joining the data
head(Ensembl_Uni_NCBI, n=10) # showing first n rows

# readr::write_csv(x = Ensembl_Uni_NCBI, file = "Ensembl_Uni_NCBI.csv") # data export

Ensembl_Uni_NCBI_unmatched  = na_rows <- Ensembl_Uni_NCBI [is.na(Ensembl_Uni_NCBI$NCBI_GeneID), ]  # selecting rows with unmatched data
head(Ensembl_Uni_NCBI_unmatched, n=10) # showing first n rows

Ensembl_Uni_NCBI$Gene_Symbol_congruence_Ensembl_NCBI <- ifelse(
  Ensembl_Uni_NCBI$Ensembl_official_symbols == Ensembl_Uni_NCBI$Symbol_NCBI,
  'TheSame',
  'Different') # testing gene symbol congruence
head(Ensembl_Uni_NCBI, n=10) # showing first n rows


Ensembl_Uni_NCBI$Gene_Symbol_congruence_Ensembl_UniProt <- ifelse(
  Ensembl_Uni_NCBI$Ensembl_official_symbols == Ensembl_Uni_NCBI$Gene_Names_UniProt,
  'TheSame',
  'Different') # testing gene symbol congruence
head(Ensembl_Uni_NCBI, n=10) # showing first n rows

####################################################################
####### paper data import #############################################
####################################################################

Paper_data = read.delim('Paper_data.csv', sep=';', header = TRUE) # data import
head(Paper_data, n = 8) # showing first n rows

####################################################################
####### joining the dataset with paper data ###############################
####################################################################
Paper_Ensembl_Uni_NCBI = dplyr:: left_join(Paper_data, Ensembl_Uni_NCBI, by = c("Gene_symbol_paper" = "Input_gene_symbols"), keep = TRUE, na_matches = "never") # joining the data
head(Paper_Ensembl_Uni_NCBI, n=10) # showing first n rows

######################################################################
######## Final editing ###################################################
######################################################################

Paper_Ensembl_Uni_NCBI_2 = dplyr::mutate(Paper_Ensembl_Uni_NCBI, TemporaryCol =  'Gene') # Adding temporary column with the word “gene”
head(Paper_Ensembl_Uni_NCBI_2, n=4) # showing first n rows

Paper_Ensembl_Uni_NCBI_3  = tidyr::unite(Paper_Ensembl_Uni_NCBI_2, 'Anchored_NCBI_gene_symbol', c("TemporaryCol","Symbol_NCBI"), sep = "_", remove = FALSE, na.rm = FALSE) # adding new column with joint items from two other columns
head(Paper_Ensembl_Uni_NCBI_3, n=4) # showing first n rows

Paper_Ensembl_Uni_NCBI_4 = dplyr::mutate(Paper_Ensembl_Uni_NCBI_3,  TemporaryCol = NULL) # removing temporary column
head(Paper_Ensembl_Uni_NCBI_4, n=4) # showing first n rows

readr::write_csv(x = Paper_Ensembl_Uni_NCBI_4, file = "Paper_Ensembl_Uni_NCBI_final.csv") # data export

######## The End #######

