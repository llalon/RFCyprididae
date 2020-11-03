# Assignment 2 - R code only. For written sections see RMD file.

# Check which packages are missing from the system and install them.
list.of.packages <- c("tidyverse", "rentrez", "randomForest", "ggplot2", "rmarkdown", "knitr", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

# The Biostrings package must be installed from BiocManager. Check if it is missing from the system and install it. 
if (!require('Biostrings')) {
  BiocManager::install("Biostrings")
}

# Import the required libraries
library(tidyverse)
library(Biostrings)
library(rentrez)
library(randomForest)
library(ggplot2)
library(rmarkdown)

# The taxanomic group of interest.
taxa <- "Cyprinidae"

# A vector containing all genes of interest.
genes <- c("COI", "CytB")

# A vector containing the k-mer lengths. A k-mer length of 1 will determine single nucleotide proportions in the sequence. Likewise, a k-mer of 4 will determine oligionucleotide proportions of length 4. K-mer sizes are only limited by compute power and time. 
kmer.lengths <- c(1:4)

# The database that will be used for sequence queries.
search.db <- "nuccore"

# The ntree size used for the random tree machine learning algorithm. 
ntree.size = 1000


A data frame was built using sequence data of each gene from NCBIs nuccore for species in the family Cyprinidae.

# Creation of the empty data frame that is used to store sequence data for each gene.
df.sequences <- data.frame()

# Iterate through all the genes of interest, query NCBI for sequence data, and assemble the dataframe
for (i in genes) {
    # Check if fasta data already exists on the file system. Only fetch new data if it is not found.
    if (!file.exists(paste(taxa, i, "seq.fasta", sep = "_"))) {
      # Fasta data was not found. Fetch the required data.

      # Create the search term to be used to query the database.
      search.term <- paste(taxa, "[ORGN] AND ", i, "[Gene] NOT (genome[TITL])", sep = "")
      
      # Query the database using that search term.
      gene.search <- entrez_search(db = search.db, term = search.term)
      
      # Get the count of the maximum hits for the query.
      search.count <- gene.search$count
      
      # Query again using the max hits as retmax.
      gene.search <- entrez_search(db = search.db, term = search.term, retmax = search.count, use_history = TRUE)
      
      # Fetch the data in fasta format.
      gene.fetch <- entrez_fetch(db = search.db, web_history = gene.search$web_history, rettype = "fasta")
      
      # Write the fasta data to disk. The location to cache fasta data is the current R working directory.
      write(gene.fetch, file = paste(taxa, i, "seq.fasta", sep = "_"), sep = "\n")
    }
    
    # Read back in the fasta data into a DNA StringSet object.
    string.set <- readDNAStringSet(paste(taxa, i, "seq.fasta", sep = "_"))
    
    # Build a dataframe to store the DNA sequences.
    df.gene <- data.frame(gene_title = names(string.set), gene_seq = paste(string.set))
    
    # Create a new column for the gene name in the data frame. 
    df.gene$gene_name <- i
    
    # Append the sequence dataframe with the queried gene data.
    df.sequences <- rbind(df.sequences, df.gene)
}

# Remove all objects that are no longer needed.
rm(search.term, gene.search, search.count, gene.fetch, string.set, df.gene)


# Create an empty data frame to store meta data about the sequences.
df.seq.metadata <- data.frame(gene = character(), n = numeric(), min = numeric(), max = numeric(), median = numeric(), missing = numeric())

# Iterate through each gene, determine sequence meta data, and store it within the meta data data frame.
for (i in genes) {
  iN <- nrow(df.sequences[df.sequences$gene_name == i, ]) # Sample size
  iMin <- min(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Lowest sequence length
  iMax <- max(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Highest sequence length
  iMed <- median(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Median sequence length
  iMiss <- max(str_count(df.sequences[df.sequences$gene_name == i, ]$gene_seq, "N")) / min(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Highest proportion of missing values. Calculated as the largest number of N's divided by the shortest sequence. These are compared between all observations and not within the lowest, highest value ones.

  # Append all this information into the meta data data frame.
  df.seq.metadata <- rbind(df.seq.metadata, data.frame(gene = i, n = iN, min = iMin, max = iMax, median = iMed, missing = iMiss))
}

# Remove all objects that are no longer needed.
rm(iN, iMin, iMax, iMed, iMiss)



# Generate a histogram of sequence lengths for each gene. All genes are displayed on a single stacked histogram.
ggplot(df.sequences, aes(nchar(gene_seq), fill = gene_name)) +
  geom_histogram(binwidth = 50) +
  labs(title = "Histogram of DNA Sequence Lengths", x = "Sequence Length", y = "Frequency", fill = "Gene") +
  theme_minimal()


# Create a new data frame with human readable names for the columns.
df.seq.metadata2 <- df.seq.metadata
names(df.seq.metadata2) <- c("Gene", "Sample Size", "Min Sequence Length", "Max Sequence Length", "Median Sequence Length", "Highest Proportion of Missing Nucleotides")



# Determine the 1st and 3rd quantile of sequence lengths for all genes.
q1 <- quantile(nchar(df.sequences$gene_seq), probs = 0.25, na.rm = TRUE)
q3 <- quantile(nchar(df.sequences$gene_seq), probs = 0.75, na.rm = TRUE)

# Remove observations that do not fall within the 1st and 3rd quantile.
df.sequences <- df.sequences %>%
  filter((str_count(gene_seq) >= q1 & str_count(gene_seq) <= q3))

# Remove all objects that are no longer needed.
rm(q1, q3)


# Remove existing data of untrimmed sequences from the meta data data frame.
df.seq.metadata <- data.frame(gene = character(), n = numeric(), min = numeric(), max = numeric(), median = numeric(), missing = numeric())

# Iterate through each gene, determine sequence meta data, and store within the meta data data frame.
for (i in genes) {
  iN <- nrow(df.sequences[df.sequences$gene_name == i, ]) # Sample size
  iMin <- min(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Lowest sequence length
  iMax <- max(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Highest sequence length
  iMed <- median(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Median sequence length
  iMiss <- max(str_count(df.sequences[df.sequences$gene_name == i, ]$gene_seq, "N")) / min(nchar(df.sequences[df.sequences$gene_name == i, ]$gene_seq)) # Highest proportion of missing values. Calculated as the largest number of N's divided by the shortest sequence.

  # Append all this information into the meta data data frame.
  df.seq.metadata <- rbind(df.seq.metadata, data.frame(gene = i, n = iN, min = iMin, max = iMax, median = iMed, missing = iMiss))
}

# Remove all objects that are no longer needed.
rm(iN, iMin, iMax, iMed, iMiss)


# Generate a histogram of sequence lengths for each gene. All genes are displayed on a single stacked histogram.
ggplot(df.sequences, aes(nchar(gene_seq), fill = gene_name)) +
  geom_histogram(binwidth = 50) +
  labs(title = "Histogram of DNA Sequence Lengths", x = "Sequence Length", y = "Frequency", fill = "Gene") +
  theme_minimal()


# Create a new data frame with human readable names for the columns.
df.seq.metadata2 <- df.seq.metadata
names(df.seq.metadata2) <- c("Gene", "Sample Size", "Min Sequence Length", "Max Sequence Length", "Median Sequence Length", "Highest Proportion of Missing Nucleotides")


# Convert the sequences from character to a DNAStringSet object in order to use them with Biostrings functions.
df.sequences$gene_seq <- DNAStringSet(df.sequences$gene_seq)

# Iterate through K-mer lengths, calculate the K-mer for that length and store its value in a new column of the dataframe. The count of each individual nucleotide is not needed as kmers of length 1 indicate proportions of each individual nucleotide.
for (i in kmer.lengths) {
  df.sequences <- cbind(df.sequences, as.data.frame(oligonucleotideFrequency(x = df.sequences$gene_seq, width = i, as.prob = TRUE)))
}

# Convert the DNAStringSet data back into character data. Biostring functions are no longer needed. Character data is easier to work with.
df.sequences$gene_seq <- as.character(df.sequences$gene_seq)


# Set the seed to enable reproducible results.
set.seed(4732)

# Generate the validation data. The sample rate will the 20% of the lowest number of observations per gene. I.e If there are 1000 observations for COI, and 2000 for CytB it will use 20% of 1000.
df.validation <- df.sequences %>%
  group_by(gene_name) %>%
  sample_n(0.2 * min(c(nrow(df.sequences[df.sequences$gene_name == "CytB", ]), nrow(df.sequences[df.sequences$gene_name == "COI", ]))))

# change the seed before creating the training data set.
set.seed(6543)

# Generate the training data. The sample rate will be 80% of the lowest number of observations per gene.
df.training <- df.sequences %>%
  # Ensure no data overlaps with the validation data set by removing everything sampled for that set.
  filter(!gene_title %in% df.validation$gene_title) %>%
  group_by(gene_name) %>%
  sample_n(0.8 * min(c(nrow(df.sequences[df.sequences$gene_name == "CytB", ]), nrow(df.sequences[df.sequences$gene_name == "COI", ]))))

# Create a function to determine index positions of K-mer proportions. This function takes the K-mer length as a argument, and returns the index positions in the data frame for that length of K-mer. I.e for a K-mer of length 2 it will return 8:23 corresponding to columns 8-23 in the data frame. This function is needed in order to build dynamic random forest models.
getKmerIndex <- function(kmer) {
  # The index positions for a kmer of length 'kmer'. Follows exponential curve with base 4 (4 nucleotides)
  start <- (1 / 3) * ((4 ^ kmer) + 8)
  end <- (1 / 3) * ((4 ^ (kmer + 1)) + 5)
  
  # Return those value as a vector as this is what the random forest algorithm requires.
  return(c(start:end))
}

# Create a data frame to store mean error rates, prediction accuracy, and oob values for random forest models generated for each K-mer length.
df.machine.data <- data.frame(kmer_length = numeric(), error_rate = numeric(), compute_time = numeric(), predict_acc = numeric(), oob = numeric())

# Iterate through each K-mer length and build a random forest model for that K-mer length using the previously made training and validating data sets.
for (i in kmer.lengths) {
  
  # In order to determine the compute time for the random forest model generation, the system time before and after will need to be compared.
  t1 <- as.numeric(Sys.time())
  
  # Build the random forest model using the previously created training data set, the K-mer length i, and the ntree size.
  rf.classifier <- randomForest(x = df.training[, getKmerIndex(i)], y = as.factor(df.training$gene_name), ntree = ntree.size, importance = TRUE)
  
  # Determine the total execution time based on how longs it has been since before the random forest began.
  t2 <- as.numeric(Sys.time()) - t1
  
  # Test the model using our validation data.
  predict.validation <- predict(rf.classifier, df.validation[, c(i, getKmerIndex(i))])
  
  # Append data on the generated model to the dataframe.  
  df.machine.data <- rbind(df.machine.data, data.frame(kmer_length = i, error_rate = mean(rf.classifier$err.rate), compute_time = t2, predict_acc = mean(predict.validation == df.validation$gene_name), oob = mean(rf.classifier$oob.times)))
}

# Removed objects that are no longer needed.
rm(t1, t2, i, predict.validation, rf.classifier)

# Create a line plot using ggplot2 to show mean error rates and compute times with respect to K-mer lengths for each gene classifier random forest model. A Vertical line on the plot will indicate the first point that has 100% accuracy with the validation data set.

# As the magnitudes of error rates and compute times differ greatly, create a variable coefficient to normalize the scales. This number is based on the maximum values of each and the max k-mer size.
scale.coef <- max(df.machine.data$compute_time) * (10 ^ max(kmer.lengths)) * max(df.machine.data$error_rate)

# Find the lowest K-mer length that yields 100% accuracy for the validation data.
lowest.kmer <- min(df.machine.data[df.machine.data$predict_acc == 1, ]$kmer_length)

# Declare colors to be used for the plot. These colors are consistent with previous plots.
color1 <- "#F8766D"
color2 <- "#00BFC4"

# Generate the plot.
ggplot(df.machine.data, aes(x=kmer_length)) +
  
  # Plot each line on the graph. Normalizing compute time based on the scale coefficient. 
  geom_line(aes(y = compute_time / scale.coef), size=2, color=color1) + 
  geom_line(aes(y = error_rate), size=2, color=color2) +
  
  scale_y_continuous(
    
    # Label for the left axis
    name = "Error Rate",
    
    # Add a second axis, normalize based on the scale coefficient, and set the label.
    sec.axis = sec_axis(~.*scale.coef, name="Compute Time (s)")
  ) +
  
  # Plot a vertical line to indicate the first K-mer value that has 100% accuracy for the validation data set.
  geom_vline(xintercept = lowest.kmer) +
  
  # Label the plot.
  labs(title = "Error rate and Compute Time with respect to K-Mer Length for a Random Forest Based Gene Classifier", x = "K-Mer Length") +
  
  # Apply the theme to the plot.
  theme(
    axis.title.y = element_text(color = color2, size=13),
    axis.title.y.right = element_text(color = color1, size=13)
  ) +
  theme_minimal()


# Remove objects that are no longer needed.
rm(scale.coef, color1, color2, lowest.kmer)

# Display the table showing information on our random forest models.

# Create a new data frame with human readable names for the columns.
df.machine.data2 <- df.machine.data
names(df.machine.data2) <- c("K-mer Length", "Mean Error Rate", "Compute Time (s)", "Predictive Model Accuracy", "Out Of Bag Error Estimate")
