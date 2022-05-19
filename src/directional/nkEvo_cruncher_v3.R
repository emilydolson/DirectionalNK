################################
# Draft statistics: NK model
# Data processing and statistics for Directional_NK
# 
# Sonia Singhal
# May 2022
###############################

library(ggplot2)
library(reshape)

# Read in the files
# Main folder where files are located
# THIS SHOULD BE CHANGED TO MATCH FILE LOCATIONS
base <- "C:/Users/Sonia/Documents/NK model/"
# Read in the file that contains the names of the output documents
files <- read.table(paste(base, "files_may_18.txt", sep = ""), stringsAsFactors = F)
nval <- 15
ntimes = nrow(files)
# Compile the summary statistics, populations, and (reproductive) landscapes into individual data frames
#row <- 1
filename <- files[1,]
subfolder <- strsplit(filename, "/")[[1]][1]
subfile <- strsplit(filename, "/")[[1]][2]
# read in the statistics
trace <- read.table(paste(base, subfolder, "/stats-", subfile, sep = ""), header = T, stringsAsFactors = F)
# Read in the population output
pops <- read.table(paste(base, subfolder, "/popt-", subfile, sep = ""), header = T, stringsAsFactors = F, col.names = c("Gen", "Genome", "Count", "Fitness", "ShockFitness"))
# get N, K, mu, replicate (seed), and regime arguments
subpieces <- strsplit(subfile, "-")
kval <- as.numeric(substring(subpieces[[1]][2], 2))
# seed <- as.numeric(substring(subpieces[[1]][3], 2)) # if seed is followed by another argument
seed <- as.numeric(substring(strsplit(subpieces[[1]][3], "[.]")[[1]][1], 2)) # if seed ends argument list
muval <- ifelse(length(subpieces[[1]]) == 5, as.numeric(substring(subpieces[[1]][4], 2)), as.numeric(paste(substring(subpieces[[1]][4], 2), subpieces[[1]][5], sep = "-"))) 
regparam <- subpieces[[1]][length(subpieces[[1]])]
regval <- ifelse(length(regmatches(regparam, gregexpr("S", regparam))[[1]]) > 0, "Sudden", "Gradual")
# read in the reproductive landscape (requires the N value and the K value).
# For the shock landscape, use "/sland-" instead of "/land-".
peaks <- read.table(paste(base, subfolder, "/land-", subfile, sep = ""), skip = (6+nval+2^(kval+1)), col.names = c("Genome", "Fitness"))
# Assign treatment parameters to each data frame
trace$k <- kval
pops$k <- kval
peaks$k <- kval
trace$n <- nval
pops$n <- nval
peaks$n <- nval
trace$seed <- seed
pops$seed <- seed
peaks$seed <- seed
trace$regime <- regval
pops$regime <- regval
peaks$regime <- regval
trace$mu <- muval
pops$mu <- muval
peaks$mu <- muval
for (f in 2:nrow(files)){
  filename <- files[f,]
  subfolder <- strsplit(filename, "/")[[1]][1]
  subfile <- strsplit(filename, "/")[[1]][2]
  thistrace <- read.table(paste(base, subfolder, "/stats-", subfile, sep = ""), header = T, stringsAsFactors = F)
  thispop <- read.table(paste(base, subfolder, "/popt-", subfile, sep = ""), header = T, stringsAsFactors = F, col.names = c("Gen", "Genome", "Count", "Fitness", "ShockFitness"))
  subpieces <- strsplit(subfile, "-")
  kval <- as.numeric(substring(subpieces[[1]][2], 2))
  #seed <- as.numeric(substring(subpieces[[1]][3], 2))
  seed <- as.numeric(substring(strsplit(subpieces[[1]][3], "[.]")[[1]][1], 2))
  regparam <- subpieces[[1]][length(subpieces[[1]])]
  regval <- ifelse(length(regmatches(regparam, gregexpr("S", regparam))[[1]]) > 0, "Sudden", "Gradual")
  print(paste(seed, regval))
  muval <- ifelse(length(subpieces[[1]]) == 5, as.numeric(substring(subpieces[[1]][4], 2)), as.numeric(paste(substring(subpieces[[1]][4], 2), subpieces[[1]][5], sep = "-")))
  thispeak <- read.table(paste(base, subfolder, "/land-", subfile, sep = ""), skip = (6+nval+2^(kval+1)), col.names = c("Genome", "Fitness"))
  thistrace$k <- kval
  thispop$k <- kval
  thispeak$k <- kval
  thistrace$n <- nval
  thispop$n <- nval
  thispeak$n <- nval
  thistrace$seed <- seed
  thispop$seed <- seed
  thispeak$seed <- seed
  thistrace$regime <- regval
  thispop$regime <- regval
  thispeak$regime <- regval
  thistrace$mu <- muval
  thispop$mu <- muval
  thispeak$mu <- muval
  trace <- rbind(trace, thistrace)
  pops <- rbind(pops, thispop)
  peaks <- rbind(peaks, thispeak)
}

# Extract K=0 values
trace_0 <- subset(trace, kval == 0)
#basic relationships (replication fitness):
plot(trace_0$gen, trace_0$average, col = as.factor(trace_0$regime), pch = as.numeric(as.factor(trace_0$mu)))

#paneled by treatment (Gradual/Sudden)
#for basic relationships, substitute "trace" with "trace_0"
melted_data <- melt(trace, id = c("k", "n", "seed", "regime", "mu", "gen"))
# organized with regime paneled
melted_subset <- subset(melted_data, variable %in% c("population", "uniques", "average", "max", "shockAvg", "shockMax", "diversity", "evenness"))
ggplot(data = trace, aes(x = gen, y = average, color = k, fill = k, shape = as.factor(mu), group = k)) +
  facet_grid(regime ~ ., scales = "free_y") +
  stat_summary() #default is mean_se(). Also run a plot with all points around the mean/se. Violin won't work for a time factor!

# organized with all variables paneled
# note that all k values are lumped together!
## Potentially choose a mu to depict and put the other mus in the supplement?
ggplot(data = melted_subset, aes(x = gen, y = value, color = regime, fill = regime, shape = as.factor(k), group = regime)) +
    facet_grid(variable ~ mu, scales = "free_y") +
    stat_summary()

# organized with k and variable paneled. All mu values are combined.
ggplot(data = melted_subset, aes(x = gen, y = value, color = regime, fill = regime, shape = regime, group = regime)) +
  facet_grid(variable ~ k, scales = "free_y") +
  stat_summary()

# Check number of persisting generations
ntimes <- length(unique(trace$seed))
survivals <- data.frame(n = rep(NA, ntimes), k = rep(NA, ntimes), seed = rep(NA, ntimes), regime = rep(NA, ntimes), cutoff = rep(NA, ntimes), gen = rep(NA, ntimes), uniques = rep(NA, ntimes), average = rep(NA, ntimes), stdev = rep(NA, ntimes), max = rep(NA, ntimes), shockcut = rep(NA, ntimes), shockAvg = rep(NA, ntimes), shockStd = rep(NA, ntimes), shockMax = rep(NA, ntimes), mu = rep(NA, ntimes), diversity = rep(NA, ntimes), evenness = rep(NA, ntimes), population = rep(NA, ntimes))
for(s in 1:ntimes){
  sd <- unique(trace$seed)[s]
  sub <- subset(trace, seed == sd)
  survivals[s, ] <- subset(sub, gen == max(sub$gen), select = c(n, k, seed, regime, cutoff, gen, uniques, average, stdev, max, shockCut, shockAvg, shockStd, shockMax, mu, diversity, evenness, population))
}
# Basic boxplots of persistence times
boxplot(gen ~ k + regime, data = survivals, ylab = "Persistence (generations)", las = 2)
# working graphs of endpoint values
boxplot(uniques ~ regime + mu, data = subset(survivals, kval == 0), las = 2)
boxplot(diversity ~ regime + mu, data = subset(survivals, kval == 0), las = 2)
boxplot(evenness ~ regime + mu, data = subset(survivals, kval == 0), las = 2)
boxplot(log(population) ~ regime + mu, data = subset(survivals, kval == 0), las = 2)

# Mean + stddev of persistence times, comparing gradual and sudden, over different k
pd <- position_dodge(0.5)
ggplot(data = survivals, aes(x = k, y = gen, color = regime)) +
  #facet_grid(host ~ ., switch = "y", labeller = as_labeller(labelnames)) +
  stat_summary(fun = "mean", fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "pointrange", position = pd) + 
  theme_bw(base_size = 14) + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "grey91"),
        strip.text.y = element_text(face = "italic")) +
  scale_color_hue(name = "Treatment", labels = c("Gradual", "Sudden")) +
  ylab("Population persistence (generations)") +
  xlab("K")


# Plot fitness at the end of the run as a function of K. 
# Current code uses replicatio fitness. For shock fitness, set y = shockAvg
ggplot(data = survivals, aes(x = k, y = average, color = regime)) +
  #facet_grid(host ~ ., switch = "y", labeller = as_labeller(labelnames)) +
  stat_summary(fun = "mean", fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), geom = "pointrange", position = pd) + 
  #geom_point(stat = "mean", position = pd) +
  theme_bw(base_size = 14) + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(color = "grey91"),
        strip.text.y = element_text(face = "italic")) +
  scale_color_hue(name = "Treatment", labels = unique(survivals$regime)) +
  ylab("Final population fitness") +
  xlab("K")

##########################################
# Code for hamming distances
# NOTE THAT THIS CODE REQUIRES A CLONAL ANCESTOR

# Calculate hamming distances of two integers
hammingdist <- function(x, y){
  return (sum(as.logical(xor(intToBits(x),intToBits(y)))))
}

for(i in 1:nrow(pops)){
  # get the ancestral genome
  ancgen <- unique(pops$Genome[pops$Gen == 0 & pops$seed == pops$seed[i]])
  # calculate the hamming distance between the ancestor and the current genome
  pops$ancham[i] <- hammingdist(ancgen, pops$genome[i])
  # get fitness values for the current genome
  pops$fitdiff[i] <- pops$fit[i] - trace$average[trace$seed == pops$seed[i] & trace$gen == 0]
  pops$shockfitdiff[i] <- pops$shockfit[i] - trace$shockAvg[trace$seed == pops$seed[i] & trace$gen == 0]
  reppeaks <- subset(peaks, peaks$seed == peaks$seed[i])
  peakham <- 100
  for(j in 1:nrow(reppeaks)){
    if(hammingdist(reppeaks$Genome[j], pops$genome[i]) < peakham){
      peakham <- hammingdist(reppeaks$Genome[j], pops$genome[i])
      peakgen <- reppeaks$Genome[j]
      peakfit <- reppeaks$Fitness[j]
    } 
  }
  pops$closestPeak[i] <- peakgen
  pops$peakdist[i] <- peakham
  pops$peakfit[i] <- peakfit
}
# Distance from ancestor as a function of the number of surviving generations
plot(pops$gen, pops$ancham, col = pops$k + 1, xlab = "Number of surviving generations", ylab = "Hamming distance from ancestor", pch = as.numeric(as.factor(pops$regime)))
legend("topleft", legend = unique(pops$k), col = unique(pops$k) + 1, lty =1)

# UNDER CONSTRUCTION / REVISION
peakcorr <- data.frame(N = rep(nval, ntimes), K = rep(NA, ntimes), genos = rep(NA, ntimes), peaks = rep(NA, ntimes))
# correlate number of peaks with number of genotypes
row <- 1
for(s in 1:length(unique(survivals$seed))){
  val <- unique(survivals$seed)[s]
  k <- survivals$K[survivals$seed == val]
  peakcorr$K[row] <- k
  peakcorr$genos[row] <- survivals$Gen11genos[survivals$seed == val]
  peakcorr$peaks[row] <- stats$over71[stats$seed == val]
  row <- row + 1
}
plot(peakcorr$peaks, peakcorr$genos, xlab = "Number of peaks with fitness > 0.71", ylab = "Number of genotypes post cutoff", ylim = c(0, 2000))

# Distance of genotypes from each other
# Distance of genotypes from peaks > 0.71
## How to visualize the results?? 
## For heatmaps, want the scale to stay the same rather than changing with the maximum of that matrix (i.e., always 0-20)
# heatmap - give original matrices and use hamming as the distance measure
# lattice (levelplot) - normalize by N, set scale from 0-1
# As a result of gen distance, what is the delta fitness? --> weighted distance based on loci (multiply by some number for locus i -- get fitness contribution of each bit -- possible to calculate these differences based on 1-distant genomes)
# genotype x genotype, color by fitness difference, using landscape peaks (fitness difference & hamming distance)
# expect small W diff for small hamming dist & vv
# for 2 landscapes: Can correlate in this way
dev.set(dev.next())
pdf("NK model/hamm_hists_3Jan17.pdf", onefile = T)
for(f in 1:nrow(files)){
  filename <- files[f,]
  print(filename)
  subfolder <- strsplit(filename, "/")[[1]][1]
  subfile <- strsplit(filename, "/")[[1]][2]
  kval <- as.numeric(strsplit(subfile, "-")[[1]][2])
  sval <- as.numeric(strsplit(strsplit(subfile, "-")[[1]][3], "[.]")[[1]][1])
  popl <- read.table(paste(base, subfolder, "/pop-", subfile, sep = ""), skip = 3, col.names = c("Genome", "Count"))
  if(nrow(popl) == 0){
    next
  }
  # matrix of hamming distances within the final population
  # may need to rerun as a lower diagonal or symmetric matrix instead...?
  popham <- matrix(nrow = nrow(popl), ncol = nrow(popl))
  for(i in 1:(nrow(popl)-1)){
    for(j in (i+1):nrow(popl)){
      popham[i,j] <- hammingdist(popl$Genome[i], popl$Genome[j])
    }
  }
  hist(popham, main = paste("Between-genotype hamming distances,\nN = 20, K =", kval, "& seed =", sval), xlab = "Hamming distance")
  # matrix of hamming distances to the peaks
  peaksub <- subset(peaks, seed == sval)
  peakham <- matrix(nrow = nrow(popl), ncol = nrow(peaksub))
  for(i in 1:min(nrow(popl), nrow(peaksub))){
    for(j in 1:max(nrow(peaksub), nrow(popl))){
      ifelse(nrow(popl) < nrow(peaksub), peakham[i,j] <- hammingdist(popl$Genome[i], peaksub$Genome[j]), peakham[j, i] <- hammingdist(popl$Genome[j], peaksub$Genome[i]))
    }
  }
  hist(peakham, main = paste("Genotype-peak hamming distances,\nN = 20, K =", kval, "& seed =", sval))
  browser() 
}
dev.off()

####################################################
# This code under revision

# Check peak fitness structures
peaks <- data.frame(Genome = NA, Fitness = NA, N = NA, K = NA, seed = NA)
for(f in 1:nrow(files)){
  filename <- files[f,]
  print(filename)
  subfolder <- strsplit(filename, "/")[[1]][1]
  subfile <- strsplit(filename, "/")[[1]][2]
  #kval <- as.numeric(strsplit(subfile, "-")[[1]][2])
  kval <- as.numeric(substring(strsplit(subfile, "-")[[1]][2], 2))
  nval <- as.numeric(substring(strsplit(subfile, "-")[[1]][1], 2))
  #seed <- as.numeric(strsplit(strsplit(subfile, "-")[[1]][3], "[.]")[[1]][1])
  seed <- as.numeric(substring(strsplit(strsplit(subfile, "-")[[1]][3], "[.]")[[1]][1], 2))
  # Skip the first xx lines (which contain the epistasis patterns and the random numbers)
  thesepeaks <- read.table(paste(base, subfolder, "/land-", subfile, sep = ""), skip = (6+nval+2^(kval+1)), col.names = c("Genome", "Fitness"))
  thesepeaks$N <- nval
  thesepeaks$K <- kval
  thesepeaks$seed <- seed
  peaks <- rbind(peaks, thesepeaks)
  #browser() 
}
peaks <- peaks[-1,]

stats <- data.frame(N = rep(nval, times = ntimes), K = rep(NA, times = ntimes), seed = rep(NA, times = ntimes), total = rep(NA, times = ntimes), mean = rep(NA, times = ntimes), max = rep(NA, times = ntimes), min = rep(NA, times = ntimes))
row <- 1
#for(kval in 0:(nval - 1)){
for(kval in unique(peaks$K)){
  cat("K: ", kval)
  sub1 <- subset(peaks, K == kval)
  seeds <- unique(sub1$seed)
  for(s in 1:length(seeds)){
    val <- seeds[s]
    cat("seed: ", val, "\n")
    stats$K[row] <- kval
    stats$seed[row] <- val
    sub <- subset(sub1, seed == val)
    stats$total[row] <- nrow(sub)
    stats$max[row] <- max(sub$Fitness)
    stats$min[row] <- min(sub$Fitness)
    stats$mean[row] <- mean(sub$Fitness)
    row <- row + 1
  }
}


#######################################
# Former plots
# read in population output
row <- 1
for (f in 1:nrow(files)){
  filename <- files[f,]
  subfolder <- strsplit(filename, "/")[[1]][1]
  subfile <- strsplit(filename, "/")[[1]][2]
  popfile <- paste(base, subfolder, "/popt-", subfile, sep = "")
  if(file.exists(popfile)){
    thispop <- read.table(popfile, header = T, stringsAsFactors = F)
  }
  else {
    print (paste("No file", popfile))
    next()
  }
  subpieces <- strsplit(subfile, "-")
  kval <- as.numeric(substring(subpieces[[1]][2], 2))
  #seed <- as.numeric(substring(subpieces[[1]][3], 2))
  seed <- as.numeric(substring(strsplit(subpieces[[1]][3], "[.]")[[1]][1], 2))
  regparam <- subpieces[[1]][length(subpieces[[1]])]
  regval <- ifelse(length(regmatches(regparam, gregexpr("S", regparam))[[1]]) > 0, "Sudden", "Gradual")
  print(paste(seed, regval))
  muval <- ifelse(length(subpieces[[1]]) == 5, as.numeric(substring(subpieces[[1]][4], 2)), as.numeric(paste(substring(subpieces[[1]][4], 2), subpieces[[1]][5], sep = "-")))
  
  thispop$k <- kval
  thispop$n <- nval
  thispop$seed <- seed
  thispop$mu <- muval
  #thistrace$alpha <- aval
  thispop$regime <- regval
  pops <- rbind(pops, thispop)
}
#original code for pops
pops <- data.frame(gen = NA, genome = NA, count = NA, fit = NA, shockfit = NA,n = NA, k = NA, seed = NA, regime =NA)
for (f in 1:nrow(files)){
  filename <- files[f,]
  subfolder <- strsplit(filename, "/")[[1]][1]
  subfile <- strsplit(filename, "/")[[1]][2]
  thispop <- read.table(paste(base, subfolder, "/pop-", subfile, sep = ""), stringsAsFactors = F, skip = F, col.names = c("gen", "genome", "count", "fit", "shockfit"))
  subpieces <- strsplit(subfile, "-")
  kval <- as.numeric(substring(subpieces[[1]][2], 2))
  seed <- as.numeric(substring(strsplit(subpieces[[1]][3], "[.]")[[1]][1], 2))
  regval <- ifelse(length(regmatches(subfolder, gregexpr("S", subfolder))[[1]]) > 0, "Sudden", "Gradual")
  thispop$k <- kval
  thispop$n <- nval
  thispop$seed <- seed
  thispop$regime <- regval
  pops <- rbind(pops, thispop)
}
pops <- pops[-1,]

# Plot population size over time. Colors represent K values.
plot(trace$gen, trace$population, col = trace$k+1, xlab = "Generation", ylab = "Population size", pch = ifelse(trace$regime == "Sudden", 16, 10),ylim = c(0, 100000))
legend("topleft", legend = unique(trace$k), col = unique(trace$k)+1, lwd = 1, cex = 0.8)
# ggplot version
ggplot(data = trace, aes(x = gen, y = log(population), color = k, shape = interaction(regime, mu))) + geom_point()

# histogram of persistence under Gradual or Sudden
sudhist <- hist(survivals$gen[survivals$regime == "Sudden"], plot = F)
gradhist <- hist(survivals$gen[survivals$regime == "Gradual"], plot = F)
plot(gradhist, col = "grey", xlab = "Generations", main = "")
plot(sudhist, col = "black", add = T)
legend ("topright", legend = c("Gradual", "Sudden"), col = c("grey", "black"), lwd = 1)

# Scatterplot of persistence by K
plot(survivals$k, survivals$gen, col = ifelse(survivals$regime == "Sudden", "red", "blue"), xlab = "K", ylab = "Persistence (generations)")
legend("topright", legend = unique(survivals$regime), col = c("red", "blue"), lwd = 1)

