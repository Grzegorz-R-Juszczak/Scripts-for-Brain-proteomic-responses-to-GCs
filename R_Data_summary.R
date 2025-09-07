#####################################################################
############################ Aims ####################################
#####################################################################
#### Processing of joined data assembled from datasets  retrieved from individual #papers  and standardized in terms of protein / gene nomenclature.
#### Calculation of number of reporting papers, up fraction and expression score           # # separately for NCBI and Ensembl gene symbols
#### NCBI grouped #(Final_grouped_NCBI_gene.csv) and Ensembl grouped #(Final_grouped_Ensembl_gene.csv)
#####################################################################

############################################################
########## Data import #######################################
############################################################

library(readxl)
Data <- read_excel("Data.xlsx") # data import
head(Data, n = 5) # showing first n rows

##################################################################
##################################################################
########## calculations based on Ensembl gene symbols ##################################################################
##################################################################

##################################################################
########## Counting number of papers reporting the same gene symbol ##################################################################

Data2 = dplyr::mutate(Data,  tolower_Ensembl_gene_symbol
                      = tolower(Ensembl_gene_symbol))# creating new column with gene symbols written only in lowercase letters to avoid confusion caused by a diverse usage of the uppercase and lowercase letters in gene symbols.
head(Data2, n = 5) # showing first n rows

Data2 = tidyr::unite(Data2, 'Ensembl_gene_symbol.Author', c("Ensembl_gene_symbol", "Author"), sep = "_", remove = FALSE, na.rm = FALSE) # combining gene symbol with author names to group together repeated symbols in data from single paper
head(Data2, n = 5) # showing first n rows

EnsemblIDs= dplyr::select (Data2, c(tolower_Ensembl_gene_symbol, Ensembl_gene_symbol.Author)) # selecting two columns
head(EnsemblIDs, n = 5) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(EnsemblIDs, `tolower_Ensembl_gene_symbol`) ) # grouping data according to the items in column “tolower_Ensembl_gene_symbol”
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = "/ ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL

head(grouped, n = 5) # showing first n rows

groupedA <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedA, n = 6) # showing first n rows

groupedA$PapersReplicating <- lengths(regmatches(groupedA$Ensembl_gene_symbol.Author, gregexpr('/ ', groupedA$Ensembl_gene_symbol.Author))) # counting papers that replicated gene findings based on the number of separator “/”
head(groupedA, n = 6) # showing first n rows

groupedA2 = dplyr::mutate(groupedA, Total_nb_papers = PapersReplicating + 1) # calculating total number of papers that detected the same gene. 
head(groupedA2, n = 6) # showing first n rows

groupedA3 = dplyr::mutate(groupedA2, Ensembl_gene_symbol.Author = NULL, PapersReplicating = NULL) # removing temporary columns
head(groupedA3, n = 6) # showing first n rows

Data_3  = dplyr:: left_join(Data2, groupedA3, by = c("tolower_Ensembl_gene_symbol" = "tolower_Ensembl_gene_symbol"), keep = FALSE, na_matches = "never") # joining the entire dataset with number of papers reporting the same gene
head(Data_3, n=6) # showing first n rows

###################################################################
##### calculating up and down fraction in individual studies ##################################################################

Expression = dplyr::select (Data2, c(tolower_Ensembl_gene_symbol, Expression_paper, Ensembl_gene_symbol.Author)) # selecting columns
head(Expression, n = 5) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Expression, `Ensembl_gene_symbol.Author`) ) # grouping data according to the items in column `Ensembl_gene_symbol.Author`
me_integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(column, collapse = ", ")  
      })
  } )
grouped <- cbind(grouped, me_integrated)
grouped$data <- NULL
head(grouped, n = 5) # showing first n rows

groupedB <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedB, n = 6) # showing first n rows

groupedB$Up_count <- lengths(regmatches(groupedB$Expression_paper, gregexpr('up', groupedB$Expression_paper))) # counting gene up-regulation in individual studies
head(groupedB, n = 6) # showing first n rows

groupedB$Down_count <- lengths(regmatches(groupedB$Expression_paper, gregexpr('down', groupedB$Expression_paper))) # counting gene down-regulation in individual studies
head(groupedB, n = 6) # showing first n rows

groupedB2 = dplyr::mutate(groupedB, Up_Down_count = Up_count + Down_count) # total number of detected changes in gene expression in individual studies
head(groupedB2, n = 6) # showing first n rows


groupedB3 = dplyr::mutate(groupedB2, Paper_up_fraction = Up_count / Up_Down_count) # calculating the up fraction for each gene in individual studies
head(groupedB3, n = 6) # showing first n rows

groupedB4 = dplyr::mutate(groupedB3, Paper_down_fraction = Down_count / Up_Down_count) # calculating the down fraction for each gene in individual studies
head(groupedB4, n = 6) # showing first n rows


groupedB5  = dplyr::mutate(groupedB4, tolower_Ensembl_gene_symbol = NULL, Expression_paper = NULL) # removing unnecessary columns
head(groupedB5, n=6) # showing first n rows

Data_4  = dplyr:: left_join(Data_3, groupedB5, by = c("Ensembl_gene_symbol.Author" = "Ensembl_gene_symbol.Author"), keep = FALSE, na_matches = "never") # joining the complete dataset with results of calculations
head(Data_4, n=6) # showing first n rows

######################################################################
##### calculating mean up and down fraction for all studies ######################################################################


groupedC1 = dplyr::select(Data_4, c(Ensembl_gene_symbol.Author, tolower_Ensembl_gene_symbol, Paper_up_fraction)) # selecting columns
head(groupedC1, n = 5) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(groupedC1, `Ensembl_gene_symbol.Author`) ) # grouping data according to the items in column “Ensembl_gene_symbol.Author”
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = "/ ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL
head(grouped, n = 5) # showing first n rows

groupedC2 <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedC2, n = 6) # showing first n rows

groupedC3  = dplyr::mutate(groupedC2, Ensembl_gene_symbol.Author = NULL) # removing unnecessary columns
head(groupedC3, n=20) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(groupedC3, `tolower_Ensembl_gene_symbol`) ) # grouping data according to the items in column “tolower_Ensembl_gene_symbol” after earlier grouping by unique gene-paper ID (Ensembl_gene_symbol.Author)
me_integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(column, collapse = ", ")  
      })
  } )
grouped <- cbind(grouped, me_integrated)
grouped$data <- NULL
head(grouped, n = 5) # showing first n rows

groupedC4 <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedC4, n = 6) # showing first n rows

library(tidyr)
groupedC5 = groupedC4 %>% separate_wider_delim(Paper_up_fraction, delim = ', ', names = c("X1", "X2", "X3", "X4", "X5"), too_few = "align_start", cols_remove = TRUE) # splits by a delimiter up fractions from individual studies into separate columns
head(groupedC5, n = 6) # showing first n rows

groupedC5 <- as.data.frame(groupedC5) # changing data format from tibble to data frame. 
head(groupedC5, n = 6) # showing first n rows

str(groupedC5) # checking data type in columns

groupedC5[c('X1', 'X2', 'X3', 'X4', 'X5')] <- sapply(groupedC5[c('X1', 'X2', 'X3', 'X4', 'X5')], as.numeric) # changing data format to numeric in selected columns

str(groupedC5) # checking data type in columns


groupedC6 = dplyr::mutate(groupedC5, Mean_up_fraction = rowMeans(groupedC5[,2:6], na.rm = TRUE)) # calculating average up fraction for all papers that detected the gene
head(groupedC6, n = 6) # showing first n rows

groupedC7 = dplyr::mutate(groupedC6, Mean_down_fraction = 1 - Mean_up_fraction) # calculating down fraction
head(groupedC7, n = 6) # showing first n rows

Data_5  = dplyr:: left_join(Data_4, groupedC7, by = c("tolower_Ensembl_gene_symbol" = "tolower_Ensembl_gene_symbol"), keep = FALSE, na_matches = "never") # joining the entire dataset with results of the calculation
head(Data_5, n=6) # showing first n rows

################################################################################ calculating expression score
######################################################################

Data_6 = dplyr::mutate(Data_5, Expression_up_score = Total_nb_papers * Mean_up_fraction) # calculating expression up score
head(Data_6, n = 6) # showing first n rows

Data_7 = dplyr::mutate(Data_6, Expression_down_score = Total_nb_papers * Mean_down_fraction) # calculating expression down score
head(Data_7, n = 6) # showing first n rows

Data_7$Expression_score = pmax(Data_7$Expression_up_score, Data_7$Expression_down_score) # selecting larger value as a final expression score
head(Data_7, n = 6) # showing first n rows

Data_8 = dplyr::mutate(Data_7, Up_count = NULL, Down_count = NULL, Up_Down_count = NULL, Paper_up_fraction = NULL, Paper_down_fraction = NULL,  X1 = NULL,  X2 = NULL,  X3 = NULL,  X4 = NULL,  X5 = NULL, Mean_down_fraction = NULL, Expression_up_score = NULL, Expression_down_score = NULL, Ensembl_gene_symbol.Author = NULL) # removing temporary columns
head(Data_8, n = 6) # showing first n rows

Data_8  = dplyr:: rename(Data_8, Reporting_papers_Ensembl = Total_nb_papers, Up_fraction_Ensembl = Mean_up_fraction, Expression_score_Ensembl = Expression_score) # changing column names 
head(Data_8, n = 6) # showing first n rows

##################################################################
##################################################################
########## calculations based on NCBI gene symbols ##################################################################
##################################################################


############################################################
########## Data according to NCBI symbols #######################
############################################################

library(tidyr)
Data_9 = Data_8 %>% separate_longer_delim(Symbol_NCBI, delim = ', ') # splits by a delimiter into separate rows – used to separate NCBI symbols that could be grouped during previous stages of data processing
head(Data_9, n = 5) # showing first n rows

##################################################################
########## Counting number of papers reporting the same gene symbol ##################################################################

Data_10 = dplyr::mutate(Data_9,  tolower_NCBI_gene_symbol
                        = tolower(Symbol_NCBI)) # creating new column with gene symbols written only in lowercase letters to avoid confusion caused by a diverse usage of the uppercase and lowercase letters in gene symbols.

head(Data_10, n = 5) # showing first n rows


Data_11 = tidyr::unite(Data_10, 'NCBI_gene_symbol.Author', c("Symbol_NCBI", "Author"), sep = "_", remove = FALSE, na.rm = FALSE) # combining gene symbol with author names to group together repeated symbols in data from single paper

head(Data_11, n = 5) # showing first n rows


NCBI_IDs= dplyr::select (Data_11, c(tolower_NCBI_gene_symbol, NCBI_gene_symbol.Author)) # selecting columns
head(NCBI_IDs, n = 5) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(NCBI_IDs, `tolower_NCBI_gene_symbol`) ) # grouping data according to the items in column “tolower_NCBI_gene_symbol”
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = "/ ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL

head(grouped, n = 5) # showing first n rows

groupedA <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedA, n = 6) # showing first n rows

groupedA$PapersReplicating <- lengths(regmatches(groupedA$NCBI_gene_symbol.Author, gregexpr('/ ', groupedA$NCBI_gene_symbol.Author))) # counting papers that replicated gene findings based on the number of separator “/”

groupedA2 = dplyr::mutate(groupedA, Total_nb_papers = PapersReplicating + 1) # calculating total number of papers that detected the same gene. 

head(groupedA2, n = 6) # showing first n rows

groupedA3 = dplyr::mutate(groupedA2, NCBI_gene_symbol.Author = NULL, PapersReplicating = NULL) # removing temporary columns
head(groupedA3, n = 6) # showing first n rows

Data_12  = dplyr:: left_join(Data_11, groupedA3, by = c("tolower_NCBI_gene_symbol" = "tolower_NCBI_gene_symbol"), keep = FALSE, na_matches = "never") # joining the entire dataset with number of papers reporting the same gene

head(Data_12, n=6) # showing first n rows

###################################################################
##### calculating up and down fraction in individual studies ##################################################################

Expression = dplyr::select (Data_12, c(tolower_NCBI_gene_symbol, Expression_paper, NCBI_gene_symbol.Author)) # selecting columns
head(Expression, n = 5) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Expression, `NCBI_gene_symbol.Author`) ) #grouping data according to the items in column `NCBI_gene_symbol.Author `

me_integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(column, collapse = ", ")  
      })
  } )
grouped <- cbind(grouped, me_integrated)
grouped$data <- NULL
head(grouped, n = 5) # showing first n rows

groupedB <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedB, n = 6) # showing first n rows

groupedB$Up_count <- lengths(regmatches(groupedB$Expression_paper, gregexpr('up', groupedB$Expression_paper))) # counting gene up-regulation in individual studies

head(groupedB, n = 6) # showing first n rows

groupedB$Down_count <- lengths(regmatches(groupedB$Expression_paper, gregexpr('down', groupedB$Expression_paper))) # counting gene down-regulation in individual studies

head(groupedB, n = 6) # showing first n rows

groupedB2 = dplyr::mutate(groupedB, Up_Down_count = Up_count + Down_count) #total number of detected changes in gene expression in individual studies

head(groupedB2, n = 6) # showing first n rows

groupedB3 = dplyr::mutate(groupedB2, Paper_up_fraction = Up_count / Up_Down_count) #calculating the up fraction for each gene in individual studies
head(groupedB3, n = 6) # showing first n rows

groupedB4 = dplyr::mutate(groupedB3, Paper_down_fraction = Down_count / Up_Down_count) # calculating the down fraction for each gene in individual studies
head(groupedB4, n = 6) # showing first n rows

groupedB5  = dplyr::mutate(groupedB4, tolower_NCBI_gene_symbol = NULL, Expression_paper = NULL) # removing unnecessary columns
head(groupedB5, n=6) # showing first n rows

Data_13  = dplyr:: left_join(Data_12, groupedB5, by = c("NCBI_gene_symbol.Author" = "NCBI_gene_symbol.Author"), keep = FALSE, na_matches = "never") # joining the complete dataset with results of calculations

head(Data_13, n=6) # showing first n rows

######################################################################
##### calculating mean up and down fraction for all studies ######################################################################

groupedC1 = dplyr::select(Data_13, c(NCBI_gene_symbol.Author, tolower_NCBI_gene_symbol, Paper_up_fraction)) # selecting columns
head(groupedC1, n = 5) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(groupedC1, `NCBI_gene_symbol.Author`) ) #grouping data according to the items in column “NCBI_gene_symbol.Author”
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = "/ ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL
head(grouped, n = 5) # showing first n rows

groupedC2 <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedC2, n = 6) # showing first n rows


groupedC3  = dplyr::mutate(groupedC2, NCBI_gene_symbol.Author = NULL) # removing unnecessary columns
head(groupedC3, n=20) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(groupedC3, `tolower_NCBI_gene_symbol`) ) #grouping data according to the items in column “tolower_NCBI_gene_symbol” after earlier grouping by unique gene-paper ID (NCBI_gene_symbol.Author)

me_integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(column, collapse = ", ")  
      })
  } )
grouped <- cbind(grouped, me_integrated)
grouped$data <- NULL
head(grouped, n = 5) # showing first n rows

groupedC4 <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedC4, n = 6) # showing first n rows

library(tidyr)
groupedC5 = groupedC4 %>% separate_wider_delim(Paper_up_fraction, delim = ', ', names = c("X1", "X2", "X3", "X4", "X5"), too_few = "align_start", cols_remove = TRUE) #splits by a delimiter the up fractions from individual studies into separate columns 

head(groupedC5, n = 6) # showing first n rows

groupedC5 <- as.data.frame(groupedC5) # changing data format from tibble to data frame. 
head(groupedC5, n = 6) # showing first n rows

str(groupedC5) # checking data type in columns

groupedC5[c('X1', 'X2', 'X3', 'X4', 'X5')] <- sapply(groupedC5[c('X1', 'X2', 'X3', 'X4', 'X5')], as.numeric) # changing data format to numeric in selected columns

str(groupedC5) # checking data type in columns

groupedC6 = dplyr::mutate(groupedC5, Mean_up_fraction = rowMeans(groupedC5[,2:6], na.rm = TRUE)) # calculating average up fraction for all papers that detected the gene

head(groupedC6, n = 6) # showing first n rows

groupedC7 = dplyr::mutate(groupedC6, Mean_down_fraction = 1 - Mean_up_fraction) #calculating down fraction

head(groupedC7, n = 6) # showing first n rows

Data_14  = dplyr:: left_join(Data_13, groupedC7, by = c("tolower_NCBI_gene_symbol" = "tolower_NCBI_gene_symbol"), keep = FALSE, na_matches = "never") # joining the entire dataset with results of the calculation

head(Data_14, n=6) # showing first n rows

################################################################################ calculating expression score
######################################################################

Data_15 = dplyr::mutate(Data_14, Expression_up_score = Total_nb_papers * Mean_up_fraction) # calculating expression up score
head(Data_15, n = 6) # showing first n rows

Data_16 = dplyr::mutate(Data_15, Expression_down_score = Total_nb_papers * Mean_down_fraction) # calculating expression down score
head(Data_16, n = 6) # showing first n rows

Data_16$Expression_score = pmax(Data_16$Expression_up_score, Data_16$Expression_down_score) # selecting larger value as a final expression score
head(Data_16, n = 6) # showing first n rows

Data_16 = dplyr::mutate(Data_16, Up_count = NULL, Down_count = NULL, Up_Down_count = NULL, Paper_up_fraction = NULL, Paper_down_fraction = NULL,  X1 = NULL,  X2 = NULL,  X3 = NULL,  X4 = NULL,  X5 = NULL, Mean_down_fraction = NULL, Expression_up_score = NULL, Expression_down_score = NULL, NCBI_gene_symbol.Author = NULL) # removing temporary columns
head(Data_16, n = 6) # showing first n rows

Data_16  = dplyr:: rename(Data_16, Reporting_papers_NCBI = Total_nb_papers, Up_fraction_NCBI = Mean_up_fraction, Expression_score_NCBI = Expression_score) # changing column names 

head(Data_16, n = 6) # showing first n rows

Data_16 <- as.data.frame(Data_16) # changing data format from tibble to data frame. 
head(Data_16, n = 6) # showing first n rows

grouped <- tidyr::nest( dplyr::group_by(Data_16, `tolower_Ensembl_gene_symbol`) ) # grouping data according to the items in column “tolower_Ensembl_gene_symbol”
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = "/ ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL

head(grouped, n = 5) # showing first n rows

groupedA <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedA, n = 6) # showing first n rows

groupedA <- groupedA [order(groupedA$Expression_score_Ensembl, decreasing=TRUE, na.last = TRUE),] # data sorting
head(groupedA, n = 6) # showing first n rows

Ensembl_NA  = dplyr::filter(groupedA, Ensembl_gene_symbol == "NA") # selection of row with NA in symbol column
head(Ensembl_NA, n = 20) # showing first n rows


groupedA  = dplyr::filter(groupedA, Ensembl_gene_symbol != "NA") # removing row with NA indicating missing item in symbol column
head(groupedA, n = 6) # showing first n rows


groupedA$k1 <- ifelse(groupedA$Up_fraction_Ensembl >= 0.75, 'Up', '') # adding new column with “Up” for preferentially up-regulated proteins
head(groupedA, n = 6) # showing first n rows
groupedA$k2 <- ifelse(groupedA$Up_fraction_Ensembl <= 0.25, 'Down', '') # adding new column with “Down” for preferentially down-regulated proteins
head(groupedA, n = 6) # showing first n rows
groupedA$k3 <- ifelse(groupedA$Up_fraction_Ensembl > 0.25 & groupedA$Up_fraction_Ensembl < 0.75,
                      'Up&Down', '') # adding new column with 'Up&Down' for mixed responses
head(groupedA, n = 6) # showing first n rows

groupedA = tidyr::unite(groupedA, Response_direction, k1 : k3, sep = "", remove = TRUE, na.rm = FALSE) # Joining columns

head(groupedA, n = 6) # showing first n rows


readr::write_csv(x = groupedA, file = "Final_grouped_Ensembl_gene.csv") # data export

grouped <- tidyr::nest( dplyr::group_by(Data_16, `tolower_NCBI_gene_symbol`) ) # grouping data according to the items in column “tolower_NCBI_gene_symbol”
integrated <- purrr::map_dfr(
  .x = grouped$data,
  .f = function(data_for_given_probe) {     
    purrr::map_chr(
      .x = data_for_given_probe,
      .f  = function(column)  {
        paste(unique(column), collapse = "/ ")
      })
  } )
grouped  <- cbind(grouped, integrated)
grouped$data <- NULL

head(grouped, n = 5) # showing first n rows

groupedA <- as.data.frame(grouped) # changing data format from tibble to data frame. 
head(groupedA, n = 6) # showing first n rows

groupedA <- groupedA [order(groupedA$Expression_score_NCBI, decreasing=TRUE, na.last = TRUE),] # data sorting
head(groupedA, n = 6) # showing first n rows

NCBI_NA  = dplyr::filter(groupedA,  Symbol_NCBI == "NA") # selection of row with NA in symbol column
head(NCBI_NA, n = 20) # showing first n rows


groupedA  = dplyr::filter(groupedA,  Symbol_NCBI != "NA") # removing row with NA indicating missing item in symbol column
head(groupedA, n = 6) # showing first n rows

groupedA$k1 <- ifelse(groupedA$Up_fraction_NCBI >= 0.75, 'Up', '') # adding new column with “Up” for preferentially up-regulated proteins
head(groupedA, n = 6) # showing first n rows
groupedA$k2 <- ifelse(groupedA$Up_fraction_NCBI <= 0.25, 'Down', '') # adding new column with “Down” for preferentially down-regulated proteins

head(groupedA, n = 6) # showing first n rows
groupedA$k3 <- ifelse(groupedA$Up_fraction_NCBI > 0.25 & groupedA$Up_fraction_Ensembl < 0.75,
                      'Up&Down', '') # adding new column with 'Up&Down' for mixed responses

head(groupedA, n = 6) # showing first n rows

groupedA = tidyr::unite(groupedA, Response_direction, k1 : k3, sep = "", remove = TRUE, na.rm = FALSE) # Joining columns

head(groupedA, n = 6) # showing first n rows


readr::write_csv(x = groupedA, file = "Final_grouped_NCBI_gene.csv") # data export

#### THE END ####
