####################################################################
####### paper data import #############################################
####################################################################

Paper_data = read.delim('Paper_data.csv', sep=';', header = TRUE) # data import
head(Paper_data, n = 8) # showing first n rows

UniProt_ID_paper_final = unique(Paper_data$UniProt_ID_paper_final) # list of unique IDs
head(UniProt_ID_paper_final, n = 8) # showing first n rows

UniProt_ID_paper_final2 = data.frame(UniProt_ID_paper_final) # creating data frame
head(UniProt_ID_paper_final2, n = 8) # showing first n rows

####################################################################
####### UniProt data import and pre-processing ############################
####################################################################

UniProt_data <-readr::read_tsv("UniProt_data.tsv", col_names = TRUE) # data import
head(UniProt_data, n = 8) # showing first n rows

UniProt_data <- as.data.frame(UniProt_data) # changing data format to data frame. 
head(UniProt_data, n = 8) # showing first n rows

UniProt_data <- UniProt_data [, c("From", "Entry", "Entry Name", "Protein names", "Gene Names (primary)", "GeneID", "HGNC", "Organism", "Organism (ID)")] # setting column order
head(UniProt_data, n = 8) # showing first n rows

colnames(UniProt_data) <- c("From", "Entry_UniProt", "Entry_Name_UniProt", "Protein_names_UniProt", "Gene_Names_UniProt",  "NCBI_Gene_ID_UniProt", "HGNC_ID_UniProt", "Organism_UniProt",  "Organism_ID_UniProt") # setting column names
head(UniProt_data, n=8) # showing first n rows

UniProt_data2 = dplyr::mutate(UniProt_data, NCBI_Gene_ID_UniProt_clean  = substr(NCBI_Gene_ID_UniProt, 1, nchar(NCBI_Gene_ID_UniProt) - 1)) # removing “;” located at the end of NCBI IDs
head(UniProt_data2, n = 5) # showing first n rows

UniProt_data3 = dplyr::mutate(UniProt_data2, HGNC_ID_UniProt_clean  = substr(HGNC_ID_UniProt, 1, nchar(HGNC_ID_UniProt) - 1)) # removing “;” located at the end of the IDs
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

*Check data frame UniProt_data4_duplicatted for duplicated UniProt IDs in the “From” column.  If necessary, filter required species with the next command (silenced with #). Otherwise, go to the “UniProt_data5 = dplyr:: bind_rows” action. 
#UniProt_data4_duplicatted = dplyr::filter(UniProt_data4_duplicatted, Organism_UniProt == "Homo sapiens (Human)") # selection of human identifiers
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
####### Joining the paper and UniProt data #################################
####################################################################

Pap_Uni  = dplyr:: full_join(UniProt_ID_paper_final2, UniProt_data8, by = c("UniProt_ID_paper_final" = "From"), keep = FALSE, na_matches = "never") # joining the data
head(Pap_Uni, n=10) # showing first n rows

####################################################################
####### NCBI data import and pre-processing ############################
####################################################################

ncbi_dataset <-readr::read_tsv("ncbi_dataset.tsv", col_names = TRUE) # data import
head(ncbi_dataset, n = 8) # showing first n rows

ncbi_dataset <- as.data.frame(ncbi_dataset) # changing data format to data frame. 
head(ncbi_dataset, n = 8) # showing first n rows
ncbi_dataset <- ncbi_dataset [, c("NCBI GeneID", "Symbol", "Description", "Gene Type", "Nomenclature ID", "Taxonomic ID", "SwissProt Accessions")] # setting column order
head(ncbi_dataset, n = 8) # showing first n rows

colnames(ncbi_dataset) <- c("NCBI_GeneID", "Symbol_NCBI", "Description_NCBI", "NCBI_Gene_Type", "Nomenclature_ID_NCBI", "Taxonomic_ID_NCBI", "SwissProt_Accessions_NCBI") # setting column names
head(ncbi_dataset, n = 8) # showing first n rows

ncbi_dataset$NCBI_GeneID <-as.character(ncbi_dataset$NCBI_GeneID) # changing data type in selected column
ncbi_dataset$Nomenclature_ID_NCBI <-as.character(ncbi_dataset$Nomenclature_ID_NCBI) # changing data type in selected column
head(ncbi_dataset, n = 8) # showing first n rows


####################################################################
####### joining paper, UniProt and NCBI data based on NCBI IDs################
####################################################################

Pap_Uni_NCBI = dplyr:: left_join(Pap_Uni, ncbi_dataset, by = c("NCBI_Gene_ID_UniProt_clean" = "NCBI_GeneID"), keep = TRUE, na_matches = "never") # joining the data
head(Pap_Uni_NCBI, n=10) # showing first n rows

Pap_Uni_NCBI_matched <- dplyr::filter(Pap_Uni_NCBI, NCBI_GeneID != "NA") # selecting rows with matched data
head(Pap_Uni_NCBI_matched, n=10) # showing first n rows

####################################################################
####### joining unmatched paper / UniProt and NCBI data based on Uniprot IDs####
####################################################################

Pap_Uni_NCBI_unmatched  = na_rows <- Pap_Uni_NCBI[is.na(Pap_Uni_NCBI$NCBI_GeneID), ]  # selecting rows with unmatched data
head(Pap_Uni_NCBI_unmatched, n=10) # showing first n rows

Pap_Uni_NCBI_unmatched2  = dplyr::mutate(Pap_Uni_NCBI_unmatched, 
                                         NCBI_GeneID = NULL, 
                                         Symbol_NCBI = NULL,
                                         Description_NCBI = NULL,
                                         NCBI_Gene_Type = NULL,
                                         Nomenclature_ID_NCBI = NULL,
                                         Taxonomic_ID_NCBI = NULL,
                                         SwissProt_Accessions_NCBI = NULL) # removing NCBI columns from previous joining
head(Pap_Uni_NCBI_unmatched2, n = 8) # showing first n rows

Pap_Uni_NCBI_V2 = dplyr:: left_join(Pap_Uni_NCBI_unmatched2, ncbi_dataset, by = c("Entry_UniProt" = "SwissProt_Accessions_NCBI"), keep = TRUE, na_matches = "never") # joining the data
head(Pap_Uni_NCBI_V2, n=10) # showing first n rows

Pap_Uni_NCBI_matched$Nomenclature_ID_NCBI <-as.character(Pap_Uni_NCBI_matched$Nomenclature_ID_NCBI) # changing data type in selected column 
head(Pap_Uni_NCBI_matched, n=4) # showing first n rows

Pap_Uni_NCBI_V3 = dplyr:: bind_rows(Pap_Uni_NCBI_matched, Pap_Uni_NCBI_V2) # joining data
head(Pap_Uni_NCBI_V3, n=10) # showing first n rows

####################################################################
####### joining unmatched paper / UniProt and NCBI data based on HGNC IDs####
####################################################################
Pap_Uni_NCBI_V3_matched <- dplyr::filter(Pap_Uni_NCBI_V3, NCBI_GeneID != "NA") # selecting rows with matched data
head(Pap_Uni_NCBI_V3_matched, n=10) # showing first n rows

Pap_Uni_NCBI_V3_matched$Nomenclature_ID_NCBI <-as.character(Pap_Uni_NCBI_V3_matched$Nomenclature_ID_NCBI) # changing data type in selected column 
head(Pap_Uni_NCBI_V3_matched, n=4) # showing first n rows

Pap_Uni_NCBI_V3_unmatched  = na_rows <- Pap_Uni_NCBI_V3 [is.na(Pap_Uni_NCBI_V3$NCBI_GeneID), ]  # selecting rows with unmatched data
head(Pap_Uni_NCBI_V3_unmatched, n=10) # showing first n rows

Pap_Uni_NCBI_V3_unmatched2  = dplyr::mutate(Pap_Uni_NCBI_V3_unmatched, 
                                            NCBI_GeneID = NULL, 
                                            Symbol_NCBI = NULL,
                                            Description_NCBI = NULL,
                                            NCBI_Gene_Type = NULL,
                                            Nomenclature_ID_NCBI = NULL,
                                            Taxonomic_ID_NCBI = NULL,
                                            SwissProt_Accessions_NCBI = NULL) # removing NCBI columns from previous joining
head(Pap_Uni_NCBI_V3_unmatched2, n = 8) # showing first n rows

Pap_Uni_NCBI_V3 = dplyr:: left_join(Pap_Uni_NCBI_V3_unmatched2, ncbi_dataset, by = c("HGNC_ID_UniProt_clean" = "Nomenclature_ID_NCBI"), keep = TRUE, na_matches = "never") # joining the data
head(Pap_Uni_NCBI_V3, n=10) # showing first n rows

Pap_Uni_NCBI_V4 = dplyr:: bind_rows(Pap_Uni_NCBI_V3_matched, Pap_Uni_NCBI_V3) # joining data
head(Pap_Uni_NCBI_V4, n=10) # showing first n rows

Pap_Uni_NCBI_V4$Gene_Symbol_congruence_UniProt_NCBI <- ifelse(
  Pap_Uni_NCBI_V4$Gene_Names_UniProt == Pap_Uni_NCBI_V4$Symbol_NCBI,
  'TheSame',
  'Different') # testing gene symbol congruence
head(Pap_Uni_NCBI_V4, n=10) # showing first n rows


##################################################################
####### Downloading data from Ensembl ################################
##################################################################
library(biomaRt) # package selection 

ensembl = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl") # selection of Ensembl biomart database and dataset

Attributes = listAttributes(mart = ensembl) # list of available attribute codes
head(Attributes, n = 4) # showing first n rows
readr::write_csv(x = Attributes, file = "EnsemblAttributes.csv") # data export

Filters = listFilters(mart = ensembl) # list of available filter codes
head(Filters, n = 4) # showing first n rows
readr::write_csv(x = Filters, file = "EnsemblFilters.csv") # data export

EnsemblDatabases = listEnsembl() # list of Ensembl databases
head(EnsemblDatabases, n = 20) # showing first n rows
readr::write_csv(x = EnsemblDatabases, file = "EnsemblDatabases.csv") # data export

EnsemblDatasets = listDatasets(mart = ensembl) # listing available dataset codes
head(EnsemblDatasets, n = 20) # showing first n rows
readr::write_csv(x = EnsemblDatasets, file = "EnsemblDatasets.csv") # data export

Ensembl = getBM(attributes=c('external_gene_name', 'description', 'ensembl_gene_id', 'entrezgene_id', 'hgnc_id', 'gene_biotype'), mart = ensembl) # import of data from Ensembl 
head(Ensembl, n = 20) # showing first n rows
readr::write_csv(x = Ensembl, file = "Ensembl.csv") # data export

colnames(Ensembl) <- c("Ensembl_gene_symbol", "Ensembl_description", "Ensembl_gene_id", "NCBI_ID_Ensembl", "HGNC_ID_Ensembl", "Ensembl_gene_type") # setting column names
head(Ensembl, n = 8) # showing first n rows

Ensembl$NCBI_ID_Ensembl <-as.character(Ensembl$NCBI_ID_Ensembl) # changing data type in selected column 

##################################################################
####### Joining data #################################################
##################################################################

Pap_Uni_NCBI_Ensembl = dplyr:: left_join(Pap_Uni_NCBI_V4, Ensembl, by = c("NCBI_GeneID" = "NCBI_ID_Ensembl"), keep = TRUE, na_matches = "never") # joining the data
head(Pap_Uni_NCBI_Ensembl, n=10) # showing first n rows


Pap_Uni_NCBI_Ensembl_matched <- dplyr::filter(Pap_Uni_NCBI_Ensembl, Ensembl_gene_id != "NA") # selecting rows with matched data
head(Pap_Uni_NCBI_Ensembl_matched, n=10) # showing first n rows

####################################################################
####### joining unmatched data based on HGNC ID####
####################################################################

Pap_Uni_NCBI_Ensembl_unmatched  = na_rows <- Pap_Uni_NCBI_Ensembl [is.na(Pap_Uni_NCBI_Ensembl$Ensembl_gene_id), ]  # selecting rows with unmatched data
head(Pap_Uni_NCBI_Ensembl_unmatched, n=10) # showing first n rows

Pap_Uni_NCBI_Ensembl_unmatched2  = dplyr::mutate(Pap_Uni_NCBI_Ensembl_unmatched, 
                                                 Ensembl_gene_symbol = NULL, 
                                                 Ensembl_description = NULL,
                                                 Ensembl_gene_id = NULL,
                                                 NCBI_ID_Ensembl = NULL,
                                                 HGNC_ID_Ensembl = NULL,
                                                 Ensembl_gene_type = NULL) # removing NCBI columns from previous joining
head(Pap_Uni_NCBI_Ensembl_unmatched2, n = 8) # showing first n rows

Pap_Uni_NCBI_Ensembl_V2 = dplyr:: left_join(Pap_Uni_NCBI_Ensembl_unmatched2, Ensembl, by = c("Nomenclature_ID_NCBI" = "HGNC_ID_Ensembl"), keep = TRUE, na_matches = "never") # joining the data
head(Pap_Uni_NCBI_Ensembl_V2, n=10) # showing first n rows

Pap_Uni_NCBI_Ensembl_V2_matched <- dplyr::filter(Pap_Uni_NCBI_Ensembl_V2, Ensembl_gene_id != "NA") # selecting rows with matched data
head(Pap_Uni_NCBI_Ensembl_V2_matched, n=10) # showing first n rows

Pap_Uni_NCBI_Ensembl_V2_unmatched  = na_rows <- Pap_Uni_NCBI_Ensembl_V2 [is.na(Pap_Uni_NCBI_Ensembl_V2$Ensembl_gene_id), ]  # selecting rows with unmatched data
head(Pap_Uni_NCBI_Ensembl_V2_unmatched, n=10) # showing first n rows

Pap_Uni_NCBI_Ensembl_V3 = dplyr:: bind_rows(Pap_Uni_NCBI_Ensembl_matched, Pap_Uni_NCBI_Ensembl_V2_matched) # joining data
head(Pap_Uni_NCBI_Ensembl_V3, n=10) # showing first n rows

Pap_Uni_NCBI_Ensembl_V3$Gene_Symbol_congruence_NCBI_Ensembl <- ifelse(
  Pap_Uni_NCBI_Ensembl_V3$Symbol_NCBI == Pap_Uni_NCBI_Ensembl_V3$Ensembl_gene_symbol,
  'TheSame',
  'Different') # testing gene symbol congruence
head(Pap_Uni_NCBI_Ensembl_V3, n=10) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Pap_Uni_NCBI_Ensembl_V3, `Ensembl_gene_symbol`) ) # grouping data according to the items in Ensembl_gene_symbol column
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = ", ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL

head(grouped, n = 10) # showing first n rows

Pap_Uni_NCBI_Ensembl_V4 <- as.data.frame(grouped) # changing data format to data frame

Pap_Uni_NCBI_Ensembl_V2_unmatched = dplyr::mutate(Pap_Uni_NCBI_Ensembl_V2_unmatched, Gene_Symbol_congruence_NCBI_Ensembl =  'NA') # adding new column to enable joining by row

Pap_Uni_NCBI_Ensembl_V2_unmatched$Organism_ID_UniProt <-as.character(Pap_Uni_NCBI_Ensembl_V2_unmatched$Organism_ID_UniProt) # changing data type in selected column 

Pap_Uni_NCBI_Ensembl_V2_unmatched$Taxonomic_ID_NCBI  <-as.character(Pap_Uni_NCBI_Ensembl_V2_unmatched$Taxonomic_ID_NCBI) # changing data type in selected column 

Pap_Uni_NCBI_Ensembl_V5 = dplyr:: bind_rows(Pap_Uni_NCBI_Ensembl_V4, Pap_Uni_NCBI_Ensembl_V2_unmatched) # joining data
head(Pap_Uni_NCBI_Ensembl_V5, n=10) # showing first n rows

######################################################################
######## Ambiguity_test #################################################
######################################################################

Ambiguity_test = dplyr::select (Pap_Uni_NCBI_Ensembl_V5, c(Ensembl_gene_symbol, UniProt_ID_paper_final)) # selecting columns
head(Ambiguity_test, n=10) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Ambiguity_test, `UniProt_ID_paper_final`) ) # grouping data according to the items in UniProt_ID_paper_final column
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = ", ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL
head(grouped, n = 10) # showing first n rows

Ambiguity_test2 <- as.data.frame(grouped) # changing data format to data frame
head(Ambiguity_test2, n = 10) # showing first n rows

Ambiguity_test2$Ambiguity_UniProtID_EnsemblSymbol <- ifelse(grepl(",", Ambiguity_test2$Ensembl_gene_symbol),'Ambiguous','Unique') # testing for the presence of “,” in items in the column Ensembl_gene_symbol. Presence of the “,” indicates that more than one gene is assigned to the ID.
head(Ambiguity_test2, n=20) # showing first n rows

Ambiguity_test3 = dplyr::select (Ambiguity_test2, c(UniProt_ID_paper_final, Ambiguity_UniProtID_EnsemblSymbol)) # selecting columns
head(Ambiguity_test3, n=10) # showing first n rows

######################################################################
######## Joining the dataset with results of ambiguity_test ######################
######################################################################

Pap_Uni_NCBI_Ensembl_V6 = dplyr:: full_join(Ambiguity_test3, Pap_Uni_NCBI_Ensembl_V5, by = c("UniProt_ID_paper_final" = "UniProt_ID_paper_final"), keep = FALSE, na_matches = "never") # joining the data
head(Pap_Uni_NCBI_Ensembl_V6, n=10) # showing first n rows

######################################################################
######## Joining the dataset with paper data #################################
######################################################################

Pap_Uni_NCBI_Ensembl_V7 = dplyr:: full_join(Paper_data, Pap_Uni_NCBI_Ensembl_V6, by = c("UniProt_ID_paper_final" = "UniProt_ID_paper_final"), keep = FALSE, na_matches = "never") # joining the data
head(Pap_Uni_NCBI_Ensembl_V7, n=10) # showing first n rows

######################################################################
######## Final editing ###################################################
######################################################################

Pap_Uni_NCBI_Ensembl_V8 = dplyr::mutate(Pap_Uni_NCBI_Ensembl_V7, TemporaryCol =  'Gene') # Adding temporary column with the word “gene”
head(Pap_Uni_NCBI_Ensembl_V8, n=4) # showing first n rows

Pap_Uni_NCBI_Ensembl_V9 = tidyr::unite(Pap_Uni_NCBI_Ensembl_V8, 'Anchored_NCBI_gene_symbol', c("TemporaryCol","Symbol_NCBI"), sep = "_", remove = FALSE, na.rm = FALSE) # adding new column with joint items from two other columns
head(Pap_Uni_NCBI_Ensembl_V9, n=4) # showing first n rows

Pap_Uni_NCBI_Ensembl_V10 = dplyr::mutate(Pap_Uni_NCBI_Ensembl_V9,  TemporaryCol = NULL) # removing temporary column
head(Pap_Uni_NCBI_Ensembl_V10, n=4) # showing first n rows

readr::write_csv(x = Pap_Uni_NCBI_Ensembl_V10, file = "Pap_Uni_NCBI_Ensembl_final.csv") # data export

######## The End #######
