InputData <- read.delim('InputData.csv', sep=',', header = TRUE) # data import
head(InputData, n = 20) # showing first n rows

library(biomaRt) # package selection 

######################################################################
###Downloading information about Ensembl database version and description of biomaRt codes for selection of all Ensembl biomart databases, datasets, filter parameters and output data (attributes)
######################################################################

#### exemplary datasets ####
#dataset,description,version
#hsapiens_gene_ensembl,Human genes (GRCh38.p14),GRCh38.p14
#mmusculus_gene_ensembl,Mouse genes (GRCm39),GRCm39
#rnorvegicus_gene_ensembl,Rat genes (mRatBN7.2),mRatBN7.2
#sscrofa_gene_ensembl,Pig genes (Sscrofa11.1),Sscrofa11.1

ensembl = useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl") # selection of Ensembl biomart database and dataset

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

################################################################
#### Filtering procedure 1 – input gene symbols treated as official symbols ####
################################################################

Filter_1_Output = getBM(attributes= 'external_gene_name', filters ='external_gene_name', values = InputData, mart = ensembl) # data import from Ensembl using list of genes from InputData for data filtering
head(Filter_1_Output, n = 20) # showing first n rows

Filter_1_Output_2 = dplyr::mutate(Filter_1_Output,  Tolower_symbol = tolower(external_gene_name)) # adding new column with gene symbols written only in lowercase letters for subsequent merging of datasets
head(Filter_1_Output_2, n = 20) # showing first n rows

InputData_2 = dplyr::mutate(InputData,  Tolower_symbol = tolower(Input_gene_symbols)) # adding new column with gene symbols written only in lowercase letters for subsequent merging of datasets
head(InputData_2, n = 20) # showing first n rows

Filter_1_Output_3  = dplyr:: left_join(InputData_2, Filter_1_Output_2, by = c("Tolower_symbol" = "Tolower_symbol"), keep = FALSE, na_matches = "never") # joining the input and output data
head(Filter_1_Output_3, n=20) # showing first n rows

Filter_1_Output_4  = dplyr::mutate(Filter_1_Output_3, Tolower_symbol = NULL) # removing unnecessary columns
head(Filter_1_Output_4, n=20) # showing first n rows

Filter_1_Output_4[nrow(Filter_1_Output_4) + 1,] = NA # adding one row with “NA” to enable proper data processing in case of data frame with 0 rows. 
head(Filter_1_Output_4, n=20) # showing first n rows


grouped <- tidyr::nest( dplyr::group_by(Filter_1_Output_4, `Input_gene_symbols`) ) # grouping data according to the items in column “Input_gene_symbols”
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

head(grouped, n = 20) # showing first n rows

Filter_1_Output_grouped <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(Filter_1_Output_grouped, n = 20) # showing first n rows

colnames(Filter_1_Output_grouped) <- c("Input_gene_symbols",  "Ensembl_official_symbol") # setting column names
head(Filter_1_Output_grouped, n = 20) # showing first n rows

##################################################################
#### Filtering procedure 2 – input gene symbols treated as synonym symbols ####
##################################################################

Filter_2_Output = getBM(attributes=c('external_gene_name', 'external_synonym'), filters ='external_synonym', values = InputData, mart = ensembl) # data import from Ensembl using list of genes from InputData for data filtering

head(Filter_2_Output, n = 20) # showing first n rows

Filter_2_Output_2 = dplyr::mutate(Filter_2_Output,  Tolower_symbol = tolower(external_synonym)) # adding new column with gene symbols written only in lowercase letters for subsequent merging of datasets
head(Filter_2_Output_2, n = 20) # showing first n rows

Filter_2_Output_3  = dplyr:: left_join(InputData_2, Filter_2_Output_2, by = c("Tolower_symbol" = "Tolower_symbol"), keep = FALSE, na_matches = "never") # joining the input and output data
head(Filter_2_Output_3, n=20) # showing first n rows

Filter_2_Output_4  = dplyr::mutate(Filter_2_Output_3, Tolower_symbol = NULL) # removal of unnecessary columns 
head(Filter_2_Output_4, n=20) # showing first n rows

Filter_2_Output_4$the_same_blind_to_capitalization <- ifelse(
  test = tolower(Filter_2_Output_4[["external_gene_name"]]) == tolower(Filter_2_Output_4[["external_synonym"]]),
  yes = T,
  no = F) # testing whether items in column external_gene_name are the same as in column external_synonym regardless of lowercase and uppercase letters
head(Filter_2_Output_4, n=20) # showing first n rows

Filter_2_Output_5 <- dplyr::filter(Filter_2_Output_4, the_same_blind_to_capitalization == "FALSE") # selecting rows with different items in columns external_gene_name and external_synonym
head(Filter_2_Output_5, n=20) # showing first n rows

Filter_2_Output_5[nrow(Filter_2_Output_5) + 1,] = NA # adding one row with “NA” to enable proper data processing in case of data frame with 0 rows. 
head(Filter_2_Output_5, n=20) # showing first n rows


grouped <- tidyr::nest( dplyr::group_by(Filter_2_Output_5, `Input_gene_symbols`)) # grouping data according to the items in column Input_gene_symbols
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

head(grouped, n = 20) # showing first n rows

Filter_2_Output_grouped <- as.data.frame(grouped) # changing data format from tibble to data frame. 

head(Filter_2_Output_grouped, n = 20) # showing first n rows

Filter_2_Output_grouped_2  = dplyr::mutate(Filter_2_Output_grouped, external_synonym = NULL, the_same_blind_to_capitalization = NULL) # removing unnecessary columns
head(Filter_2_Output_grouped_2, n = 20) # showing first n rows

colnames(Filter_2_Output_grouped_2) <- c("Input_gene_symbols", "Ensembl_official_symbol") # setting column names
head(Filter_2_Output_grouped_2, n = 20) # showing first n rows

#####################################################################
#### Joining data from filtering procedure 1 and 2 ####
#####################################################################

Joined_data = dplyr:: bind_rows(Filter_1_Output_grouped, Filter_2_Output_grouped_2) # joining datasets by rows

head(Joined_data, n=20) # showing first n rows

Joined_data <- dplyr::arrange(Joined_data, Input_gene_symbols) # data sorting according to the data Input_gene_symbols column
head(Joined_data, n=20) # showing first n rows

Joined_data_2  = dplyr::filter(Joined_data, Ensembl_official_symbol != "NA") # removing rows with missing symbols in column Ensembl_official_symbol.

head(Joined_data_2, n=20) # showing first n rows

Joined_data_2[nrow(Joined_data_2) + 1,] = NA # adding one row with “NA” to enable proper data processing in case of data frame with 0 rows. 
head(Joined_data_2, n=20) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Joined_data_2, `Input_gene_symbols`) ) # grouping data according to the items in Input_gene_symbols column
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

head(grouped, n = 20) # showing first n rows 

Joined_data_grouped <- as.data.frame(grouped) # changing data format from tibble to data frame

head(Joined_data_grouped, n = 20) # showing first n rows

Joined_data_grouped$Ambiguous_Unique <- ifelse(grepl(", ", Joined_data_grouped$ Ensembl_official_symbol),'Ambiguous','Unique') # testing for the presence of “,” in items in the column Ensembl_official_symbol. Presence of the “,” indicates that more than one official gene symbol is assigned to the synonym.
head(Joined_data_grouped, n=20) # showing first n rows

library(tidyr) # package selection
Joined_data_grouped_2 = Joined_data_grouped %>% separate_longer_delim(Ensembl_official_symbol, delim = ', ') # splitting official symbols into separate rows
head(Joined_data_grouped_2, n=20) # showing first n rows

Ensembl_official_symbol = unique(Joined_data_grouped_2$Ensembl_official_symbol) # selecting unique gene symbols for downloading additional data from Ensembl biomart
head(Ensembl_official_symbol, n=20) # showing first n rows

Ensembl_official_symbol <- as.data.frame(Ensembl_official_symbol) # changing data format to data frame
head(Ensembl_official_symbol, n=20) # showing first n rows

##################################################################
#### Downloading additional data from Ensembl ####
##################################################################

Ensembl_official_symbol_2 = getBM(attributes=c('external_gene_name', 'description', 'ensembl_gene_id', 'entrezgene_id', 'mgi_id', 'gene_biotype'), filters ='external_gene_name', values = Ensembl_official_symbol, mart = ensembl) # import of additional data from Ensembl using list of genes from Ensembl_official_symbol  for data filtering

head(Ensembl_official_symbol_2, n = 20) # showing first n rows

Ensembl_official_symbol_2[nrow(Ensembl_official_symbol_2) + 1,] = NA # adding one row with “NA” to enable proper data processing in case of data frame with 0 rows. 
head(Ensembl_official_symbol_2, n=20) # showing first n rows


Ensembl_official_symbol_2b = getBM(attributes=c('external_gene_name', 'uniprotswissprot'), filters ='external_gene_name', values = Ensembl_official_symbol, mart = ensembl) # import of additional data from Ensembl using list of genes from Ensembl_official_symbol  for data filtering

head(Ensembl_official_symbol_2b, n = 20) # showing first n rows

library(naniar)
Ensembl_official_symbol_2b2 = Ensembl_official_symbol_2b %>%
  naniar::replace_with_na(replace = list(uniprotswissprot = c(""))) # marking missing values with NA
head(Ensembl_official_symbol_2b2, n=15) # showing first n rows

Ensembl_official_symbol_2b3  = dplyr::filter(Ensembl_official_symbol_2b2, uniprotswissprot != "NA") # removing rows with missing symbols in column uniprotswissprot.
head(Ensembl_official_symbol_2b3, n=15) # showing first n rows

Ensembl_official_symbol_2b4  <- Ensembl_official_symbol_2b3  [!duplicated(Ensembl_official_symbol_2b3), ] # removing duplicated rows

head(Ensembl_official_symbol_2b4, n=15) # showing first n rows


grouped <- tidyr::nest( dplyr::group_by(Ensembl_official_symbol_2, `external_gene_name`) ) # grouping data according to the items in external_gene_name column
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

head(grouped, n = 20) # showing first n rows

Ensembl_official_symbol_3 <- as.data.frame(grouped) # changing data format to data frame

head(Ensembl_official_symbol_3, n = 20) # showing first n rows

Ensembl_official_symbol_4  = dplyr:: left_join(Ensembl_official_symbol_3, Ensembl_official_symbol_2b4, by = c("external_gene_name" = "external_gene_name"), keep = FALSE, na_matches = "never") # adding UniProt IDs that otherwise cause problems. 

head(Ensembl_official_symbol_4, n = 20) # showing first n rows


Joined_data_grouped_3  = dplyr:: left_join(Joined_data_grouped_2, Ensembl_official_symbol_4, by = c("Ensembl_official_symbol" = "external_gene_name"), keep = FALSE, na_matches = "never") # joining the input and output data
head(Joined_data_grouped_3, n=20) # showing first n rows

Joined_data_grouped_4 = dplyr::mutate(Joined_data_grouped_3, TemporaryColB =  'Output_Gene') # adding temporary column with the word “Output_Gene”
head(Joined_data_grouped_4, n=20) # showing first n rows

Joined_data_grouped_5 = tidyr::unite(Joined_data_grouped_4, 'Anchored_output_gene', c("TemporaryColB", "Ensembl_official_symbol"), sep = "_", remove = FALSE, na.rm = FALSE) # adding new column with joint items from two other columns
head(Joined_data_grouped_5, n=20) # showing first n rows

Unmatched  = dplyr:: anti_join(InputData, Joined_data_grouped_5, by = c("Input_gene_symbols" = "Input_gene_symbols")) # selecting unmatched data
head(Unmatched, n=20) # showing first n rows

Search_results = dplyr:: bind_rows(Joined_data_grouped_5, Unmatched) # joining datasets by rows
head(Search_results, n=20) # showing first n rows

Search_results_1  = dplyr::filter(Search_results, Input_gene_symbols != "NA") # removing rows with missing symbols in column Input_gene_symbols.

head(Search_results_1, n=20) # showing first n rows


Search_results_2 = dplyr::mutate(Search_results_1, TemporaryColA =  'Input_Gene') # adding temporary column with the word “Input_Gene”
head(Search_results_2, n=20) # showing first n rows

Search_results_3  = tidyr::unite(Search_results_2, 'Anchored_input_gene', c("TemporaryColA", "Input_gene_symbols"), sep = "_", remove = FALSE, na.rm = FALSE) # adding new column with joint items from two other columns
head(Search_results_3, n=20) # showing first n rows

Search_results_4  = dplyr::mutate(Search_results_3, TemporaryColA = NULL, TemporaryColB = NULL) # removing temporary columns
head(Search_results_4, n=20) # showing first n rows

colnames(Search_results_4) <- c("Anchored_input_gene", "Input_gene_symbols", "Anchored_output_gene", "Ensembl_official_symbols", "Ambiguous_Unique",  "Gene_description", "Ensembl_ID", "NCBI_ID", "MGI_ID", "Gene_type", "UniProtSwissProt") # setting column names
head(Search_results_4, n=20) # showing first n rows

readr::write_csv(x = Search_results_4, file = "Final_search_results.csv") # data export

#### THE END ####
