library(plyr)
library(data.table)
library(tidyverse)
library(chicane)

attach(loadNamespace("chicane"), name = "chicane_all")
#' function to normalize counts based on per RNA reads
#' @description 
#'      takes input of interaction matrix file file and gives output of TPM normalised counts for each RNA
#'
#' #make rna only bed file



tpm3 <- function(interaction.data){
  interaction.data$bait.length = interaction.data$bait.end - interaction.data$bait.start
  interaction.data$bait.length <- as.numeric(interaction.data$bait.length)
  interaction.data$count <- as.numeric(interaction.data$count)
  rna.split <- split(interaction.data, interaction.data$bait.id)
  for(i in seq_along(rna.split)){
    rna.split[[i]]$rpkm.count <- rna.split[[i]]$count*1e9/(sum(rna.split[[i]]$count) * rna.split[[i]]$bait.length)
  }
  interaction.data <- do.call(rbind, rna.split)
  interaction.data$rpkm.count <- as.numeric(interaction.data$rpkm.count)
  tpm.count <- NULL
  tpm.count <- interaction.data$rpkm.count * 1e6/sum(interaction.data$rpkm.count)
  interaction.data$rpkm.count <- NULL
  #interaction.data$count <- interaction.data$tpm.count
  return(tpm.count) 
}



####################
#' chicane
#'
#' @description
#'	Run full method for detecting significant interactions in capture Hi-C experiments, starting 
#'  either from a BAM file or preprocessed data from \code{prepare.data}
#'
#' @param bam 
#' 	Path to a BAM file
#' @param baits 
#'	Path to a BED file containing the baits
#' @param fragments 
#'	Path to a BED file containing all restriction fragments in the genome
#' @param interactions 
#'	Data table or path to a text file detailing fragment interactions, typically from \code{prepare.data}. 
#'	Can be used instead of bam/baits/fragments specification if the text files have already been prepared.
#' @param replicate.merging.method
#' 	Method that should be used for merging replicates, if applicable
#' @param bait.filters 
#'	Vector of length two, where the first element corresponds to the lower-end filter and the second to the upper-end filter.
#' 	When global multiple testing correction is performed, altering the bait filtering settings may affect the number of significant results.
#' @param target.filters
#' 	Vector of length two, giving lower and higher filter, respectively. 
#'	Changing this filtering setting may affect multiple testing correction by altering the number of tests performed.
#' @param distance.bins
#' 	Number of bins to split distance into. Models are fit separately in each bin.
#' @param multiple.testing.correction
#'	String specifying how multiple testing correction should be performed, by bait or globally.
#' @param verbose
#' 	Logical indicating whether to print progress reports.
#' @param interim.data.dir
#'  Path to directory to store intermediate QC data and plots. NULL indicate skip intermediate results.
#' @inheritParams prepare.data 
#' @inheritParams combine.replicates
#' @inheritParams fit.model
#' @inheritParams fit.glm
#' 
#' @return Data table with columns
#' 	\item{target.id}{String in chrN:start-end format identifying target fragment}
#' 	\item{bait.id}{String in chrN:start-end format identifying bait fragment}
#'	\item{target.chr}{Chromosome of target fragment}
#' 	\item{target.start}{Start coordinate of target fragment (zero-based)}
#'	\item{target.end}{End coordinate of target fragment}
#'	\item{bait.chr}{Chromosome of bait fragment}
#' 	\item{bait.start}{Start coordinate of bait fragment (zero-based)}
#'	\item{bait.end}{End coordinate of bait fragment}
#'	\item{bait.to.bait}{Boolean indicating if the interaction is bait-to-bait (i.e. the fragment listed as target is also a bait)}
#' 	\item{bait.trans.count}{The number of reads linking the bait to fragments in trans (a measure of "interactibility")}
#' 	\item{target.trans.count}{The number of reads linking the target to fragments in trans (a measure of "interactibility")}
#' 	\item{distance}{Distance between the midpoints of the bait and target fragments (basepairs). NA for trans interactions}
#' 	\item{count}{The number of reads linking the two fragments}
#' 	\item{expected}{The expected number of reads linking the two fragments under the fitted model}
#'	\item{p.value}{P-value for test of the observed number of reads significantly exceeding the expected count}
#'	\item{q.value}{FDR-corrected p-value}
#'
#' @import data.table
#'
#' @examples
#' \donttest{
#' if( bedtools.installed() ) {
#'   # start from BAM file
#'   bam <- system.file('extdata', 'Bre80_2q35.bam', package = 'chicane');
#'   baits <- system.file('extdata', '2q35.bed', package = 'chicane');
#'   fragments <- system.file('extdata', 'GRCh38_HindIII_chr2.bed.gz', package = 'chicane');
#'   results <- chicane(
#'		bam = bam, 
#'		baits = baits, 
#'		fragments = fragments
#'		);
#' }
#'
#' # start from pre-processed data
#' data(bre80); 
#' results <- chicane(interactions = bre80);
#' }
#'
#' @author Erle Holgersen <Erle.Holgersen@icr.ac.uk>
#'
#' @export chicane2
chicane2 <- function(
    bam = NULL,
    baits = NULL,
    fragments = NULL,
    interactions = NULL, 
    replicate.merging.method = 'sum',
    distribution = 'negative-binomial',
    include.zeros = 'none',
    bait.filters = c(0, 1),
    target.filters = c(0, 1),
    distance.bins = NULL,
    multiple.testing.correction = c('bait-level', 'global'),
    adjustment.terms = NULL,
    remove.adjacent = FALSE,
    temp.directory = NULL,
    keep.files = FALSE,
    maxit = 100,
    epsilon = 1e-8,
    cores = 1,
    trace = FALSE,
    verbose = FALSE,
    interim.data.dir = NULL
) {
  
  # TO DO:
  #	- check format of interactions data if passed directly
  
  ### INPUT TESTS ######################################################
  
  if( is.null(interactions) && (is.null(bam) || is.null(baits) || is.null(fragments)) ) {
    stop('Must provide either interactions or bam/baits/fragments.');
  }
  
  if( !is.null(interactions) &&  !(is.null(bam) || is.null(baits) || is.null(fragments)) ) {
    stop('Cannot deal with both interactions and bam/baits/fragments. Please provide one or the other.');
  }
  
  if( !is.null(interactions) && !is.character(interactions) && !is.data.table(interactions) ) {
    stop('interactions should be either a data.table object or the path to a text file generated by prepare.data');
  }
  
  if( is.null(interactions) ) {
    
    input.files <- list(bam, baits, fragments);
    
    # all should be character strings
    object.is.character <- vapply(
      input.files,
      FUN = function(x) all( is.character(x) ),
      FUN.VALUE = FALSE 
    );
    
    if( !all(object.is.character) ) {
      stop('bam, baits, and fragments should be character strings');
    }
    
    # baits and fragments should have length 1
    # bam could have longer length if replicates should be combined
    if( !all( 1 == vapply(input.files[-1], length, FUN.VALUE = 0) ) ) {
      stop('bam, baits, and fragments should have length 1');
    }
    
    # all files should exist
    for( input.file in unlist(input.files) ) {
      if( !file.exists(input.file) ) {
        error.message <- paste('File', input.file, 'does not exist');
        stop(error.message);
      }
    }
    
  }
  
  # adding in zeros happens at the data preparation step
  # throw an error if user wants it added to pre-prepared data
  if( !is.null(interactions) && 'none' != include.zeros ) {
    stop('Zeros are added in the prepare.data step. Please provide bam/baits/fragments');
  }
  
  ### MAIN #############################################################
  
  replicate.merging.method <- match.arg(replicate.merging.method);
  multiple.testing.correction <- match.arg(multiple.testing.correction);
  
  
  # prepare data
  if( is.null(interactions) ) {
    
    if( verbose ) {
      cat('PREPARING DATA\n');
    }
    
    interaction.data <- prepare.data2(
      bam, 
      baits, 
      fragments, 
      replicate.merging.method = replicate.merging.method,
      include.zeros = include.zeros,
      remove.adjacent = remove.adjacent,
      temp.directory = temp.directory,
      keep.files = keep.files,
      verbose = verbose
    );
    
  } else if( is.character(interactions) ) { 
    
    if( !file.exists(interactions) ) {
      error.message <- paste('File', interactions, 'does not exist');
      stop(error.message);
    }
    
    if( verbose ) {
      cat('READING DATA FROM FILE\n');
    }
    
    interaction.data <- data.table::fread(interactions);
    
  } else if( is.data.table(interactions) ) {
    interaction.data <- interactions;
    
    # free up memory
    rm(interactions);
  } 
  
  if( verbose ) {
    cat('FITTING MODEL\n');
  }
  
  chicane.results <- fit.model(
    interaction.data, 
    distance.bins = distance.bins, 
    distribution = distribution,
    bait.filters = bait.filters,
    target.filters = target.filters,
    adjustment.terms = adjustment.terms,
    verbose = verbose,
    cores = cores,
    maxit = maxit,
    epsilon = epsilon,
    trace = trace,
    interim.data.dir = interim.data.dir
  );
  
  chicane.results <- multiple.testing.correct(
    chicane.results,
    bait.level = 'bait-level' == multiple.testing.correction
  );
  
  # sort by q-value
  chicane.results <- chicane.results[ order(q.value, p.value) ];
  
  return(chicane.results);
}


#----------------------------------------------------------------------#


#' prepare.data
#' 
#' @description
#'	Prepare data for running interaction calling. Takes a BAM file and baits and restriction fragments as input, and returns a data table with data ready for analysis.
#'
#' @inheritParams chicane
#' @inheritParams combine.replicates
#' @param include.zeros
#' 	String specifying what zero counts to include. Options are none (default), cis, and all.
#' @param remove.adjacent
#' 	Logical indicating whether to remove all reads mapping to adjacent restriction fragments. 
#' @param temp.directory
#'	Directory where temporary files should be stored. Defaults to current directory. 
#' @param keep.files 
#' 	Logical indicating whether to keep temporary files
#'
#' @return Data table object with columns
#' 	\item{target.id}{String in chrN:start-end format identifying target fragment}
#' 	\item{bait.id}{String in chrN:start-end format identifying bait fragment}
#'	\item{target.chr}{Chromosome of target fragment}
#' 	\item{target.start}{Start coordinate of target fragment (zero-based)}
#'	\item{target.end}{End coordinate of target fragment}
#'	\item{bait.chr}{Chromosome of bait fragment}
#' 	\item{bait.start}{Start coordinate of bait fragment (zero-based)}
#'	\item{bait.end}{End coordinate of bait fragment}
#'	\item{bait.to.bait}{Boolean indicating if the interaction is bait-to-bait (i.e. the fragment listed as target is also a bait)}
#' 	\item{count}{The number of reads linking the two fragments}
#' 	\item{bait.trans.count}{The number of reads linking the bait to fragments in trans (a measure of "interactibility")}
#' 	\item{target.trans.count}{The number of reads linking the target to fragments in trans (a measure of "interactibility")}
#' 	\item{distance}{Distance between the midpoints of the bait and target fragments (basepairs). NA for trans interactions}
#' 
#' @examples
#' \donttest{
#' if( bedtools.installed() ) {
#'   bam <- system.file('extdata', 'Bre80_2q35.bam', package = 'chicane');
#'   baits <- system.file('extdata', '2q35.bed', package = 'chicane');
#'   fragments <- system.file('extdata', 'GRCh38_HindIII_chr2.bed.gz', package = 'chicane');
#'   input.data <- prepare.data(
#'		bam = bam, 
#'		baits = baits, 
#'		fragments = fragments
#'		);
#'  }
#'	}
#'
#' @export prepare.data2
prepare.data2 <- function(
    bam,
    baits,
    fragments,
    replicate.merging.method = 'sum',
    include.zeros = c('none', 'cis', 'all'),
    remove.adjacent = FALSE,
    temp.directory = NULL,
    keep.files = FALSE,
    verbose = FALSE
) {
  
  include.zeros <- match.arg(include.zeros);
  
  ### INPUT TESTS ###########################################################
  
  # check for too many combinations 
  # better to do this now to avoid having to run through all the counting steps first
  if( include.zeros %in% c('cis', 'all') ) {
    bait.ids <- read.bed(baits);
    fragment.ids <- read.bed(fragments);
    
    # check that there aren't too many combinations
    n.combinations <- get.combination.count(
      baits = bait.ids, 
      fragments = fragment.ids,
      cis.only = 'cis' == include.zeros
    );
    
    if( n.combinations > 10^9 ) {
      stop('Too many combinations, cannot include zeros');
    }
  }
  
  ### MAIN ##################################################################
  
  # convert each replicate separately
  replicate.data <- lapply(
    bam,
    convert.bam,
    baits = baits,
    fragments = fragments,
    temp.directory = temp.directory,
    keep.files = keep.files
  );
  
  # merge replicates if required
  if( 1 == length(replicate.data) ) {
    interaction.data <- replicate.data[[ 1 ]];
  } else {
    interaction.data <- combine.replicates(
      replicate.data, 
      method = replicate.merging.method
    );
  }
  
  # add cis-interaction zeroes if requested
  if( 'cis' == include.zeros ) {
    
    # store for sanity check
    nonzero.row.count <- nrow(interaction.data);
    
    # get chromosome of each fragment
    fragment.chr <- gsub('(.*):(.*)', '\\1', fragment.ids);
    bait.chr <- gsub('(.*):(.*)', '\\1', bait.ids);
    
    filled.in.data <- list();
    
    # loop over each unique chromosome
    # okay to loop over , as the processing pipeline ensures all fragments are in fragment list 
    for( chr in unique( fragment.chr ) ) {
      
      if( verbose ) {
        cat('Filling in zero counts on', chr, '\n');
      }
      
      # subset out all interactions with bait on chromosome of interest
      chr.interaction.data <- interaction.data[ bait.chr == chr ];
      
      filled.in.data[[ chr ]] <- fill.in.zeros(
        chr.interaction.data, 
        baits = bait.ids[ chr == bait.chr ],
        fragments = fragment.ids[ chr == fragment.chr ]
      );
    }
    
    # coerce one data table object
    interaction.data <- do.call(rbind, filled.in.data);
    
    # sanity check
    if( nonzero.row.count != nrow( interaction.data[ count != 0 ]) ) {
      
      error.message <- paste(
        'Internal bug - mismatched non-zero row numbers when filling in zeros\n',
        'Before:', nonzero.row.count, 'non-zero rows\n',
        'After:', nrow(interaction.data[ count != 0 ]), 'non-zero rows\n'
      );
      
      stop(error.message);
    }
    
  } else if( 'all' == include.zeros ) {
    
    # no need to throw error for size, will be done in zero-filling function
    interaction.data <- fill.in.zeros(
      interaction.data, 
      baits = bait.ids, 
      fragments = fragment.ids
    );
  }
  
    # add TPM normalised counts  
    interaction.data$count <- tpm3(interaction.data)
    # calculate trans counts and distance between fragments
    interaction.data <- add.covariates(interaction.data);
  
  # remove self-ligated fragments (which essentially also have disatance == 0)
  interaction.data <- interaction.data[ !(bait.id == target.id & distance == 0) ];
  
  if( remove.adjacent ) {
    # remove counts between adjacent fragments
    # if data has been processed with HiCUP, these are typically one end mapped 
    # to the opposite end of the adjacent restriction fragment,
    # as re-ligations are removed.
    interaction.data <- interaction.data[ !(bait.chr == target.chr & (bait.start == target.end | target.start == bait.end) ) ];
  }
  
  # want count to appear as the last column so it shows up next to expected
  new.column.order <- c( names(interaction.data)['count' != names(interaction.data)], 'count' );
  setcolorder(interaction.data, new.column.order);
  
  # sanity check that interaction data fits expected format
  # TO DO: figure out why this isn't working!
  # verify.interaction.data(interaction.data);
  
  return(interaction.data);
}


interaction.data <- fread("~/Radicl_Seq/Radicl_seq_functions/rep1_int.txt")
interaction.data$norm.count <- tpm3(interaction.data)
interaction.data$count.1 <- interaction.data$count
interaction.data$count <- interaction.data$norm.count
interaction.data$bait.trans.count.1 <- interaction.data$bait.trans.count
interaction.data$target.trans.count.1 <- interaction.data$target.trans.count
interaction.data$target.trans.count <- NULL
interaction.data$bait.trans.count <- NULL

interaction.data <- add.covariates(interaction.data)
interaction.data <- interaction.data %>% mutate_if(is.double, as.numeric)
wo_norm <- interaction.data[,c(1,2,9,10,13,14,15)]
with_norm <- interaction.data[,c(1,2,9,10,11,16,17)]
wo_norm <- wo_norm %>%
  rename(count = count.1, target.trans.count = target.trans.count.1, bait.trans.count = bait.trans.count.1)
wo_norm$type <- "Raw Counts"
with_norm$type <- "Norm Counts"
merged_runs <- rbind(wo_norm, with_norm)
bait.length <- NULL
for(i in seq_along(merged_runs$target.id)){
  bait.length[i] <- 
    as.numeric(strsplit(strsplit(merged_runs$bait.id[i],":")[[1]][2],"-")[[1]][2]) - as.numeric(strsplit(strsplit(merged_runs$bait.id[i],":")[[1]][2],"-")[[1]][1])
}
merged_runs$bait.length <-bait.length
#merged_runs$signif <- ifelse(merged_runs$q.value<0.05, TRUE,FALSE)
plot1 <- ggplot(merged_runs
       ,aes(x=bait.length,y=log10(count), color=type))+
  geom_point()+
  theme_cowplot()+
  facet_wrap(~type)+
  scale_fill_viridis(discrete = T,alpha = 0.5)+
  theme(legend.position = "none")+
  ylab("log10(Counts)")+
  xlab("Bait Length")
png("dotplot.png",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()


plot2 <- ggplot(merged_runs
       ,aes(x=bait.length,y=log10(count), color=type))+
  geom_violin(trim = T)+
  geom_boxplot(outlier.shape = NA)+
  theme_cowplot()+
  facet_wrap(~type)+
  scale_fill_viridis(discrete = T,alpha = 0.5)+
  theme(legend.position = "none")+
  ylab("log10(Counts)")+
  xlab("Bait Length")
png("violinplot.png",width = 35,height = 21,units = 'cm',res = 300)
plot(plot2)
dev.off()


