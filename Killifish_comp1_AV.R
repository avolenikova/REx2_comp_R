library(ggplot2)

#set the working directory to Re results folder
setwd("C:/Users/volena00/Desktop/")

#read first x rows from csv file with comparative analysis (nrows = x is number of clusters you want to analyze)
results <- read.delim("./COMPARATIVE_ANALYSIS_COUNTS.csv", skip = 2, nrows = 500, header = TRUE)

#sorting columns by their name
results<- results[,order(names(results))]

#name of columns is categories
catA <- c("Nfur_1", "Nfur_2", "Nfur_3")
catB <- c("Nkad_1", "Nkad_2", "Nkad_3")
catC <- c("Nort_1", "Nort_2", "Nort_3")

#count mean read counts for each category
results$NFmean <- rowMeans(subset(results, select = catA))
results$NKmean <- rowMeans(subset(results, select = catB))
results$NOmean <- rowMeans(subset(results, select = catC))


#One way ANOVA test for each row between Females and MAles (3 biological samples per group)
results$anovaFKO <- NA
for (i in 1:nrow(results))
{
Nf <- c(results[i, 2:4])
Nk <- c(results[i, 5:7])
No <- c(results[i, 8:10])
Nfko <- stack(data.frame(cbind(Nf, Nk, No)))
ANOVAresults <- aov(values ~ ind, data = Nfko)
results[i, "anovaFKO"] <- summary(ANOVAresults)[[1]][[1, "Pr(>F)"]]
}

#if ANOVA suggested significant differences assign value 1 to the column
results$significance <- NA
for (i in 1:nrow(results)) 
{
if (is.na(results[i, "anovaFKO"])) next()
else if (results[i, "anovaFKO"] < 0.005) results[i, "significance"]<- 1
else results[i, "significance"]<- NA
}

#print clusters with reads only in Nfur samples
for (i in 1:nrow(results))
{
if (is.na(results[i, "significance"])) next()
else if (results[i, "NKmean"] == 0 && results[i, "NOmean"] == 0) print(results[i, "cluster"])
}

#print clusters with reads only in Nkad samples
for (i in 1:nrow(results))
{
if (is.na(results[i, "significance"])) next()
else if (results[i, "NFmean"] == 0 && results[i, "NOmean"] == 0) print(results[i, "cluster"])
}

#print clusters with reads only in Nort samples
for (i in 1:nrow(results))
{
if (is.na(results[i, "significance"])) next()
else if (results[i, "NKmean"] == 0 && results[i, "NFmean"] == 0) print(results[i, "cluster"])
}

#print clusters with no reads in Nfur
for (i in 1:nrow(results))
{
  if (is.na(results[i, "significance"])) next()
  else if (results[i, "NFmean"] == 0) print(results[i, "cluster"])
}

#print clusters with no reads in Nkad
for (i in 1:nrow(results))
{
  if (is.na(results[i, "significance"])) next()
  else if (results[i, "NKmean"] == 0) print(results[i, "cluster"])
}

#print clusters with no reads in Nort
for (i in 1:nrow(results))
{
  if (is.na(results[i, "significance"])) next()
  else if (results[i, "NOmean"] == 0) print(results[i, "cluster"])
}

#count ratios per species
results$NFratio <- NA
results$NKratio <- NA
results$NOratio <- NA

for (i in 1:nrow(results))
{
  if (is.na(results[i, "significance"])) {next()}
  else if (results[i, "significance"] == 1 && results[i, "NFmean"] != 0 && results[i, "NKmean"] != 0 && results[i, "NOmean"] != 0)
    {meansAll <- as.numeric(c(results[i, 12:14]))
    lowest <- min(meansAll)
    results[i, "NFratio"] <- log((results[i, "NFmean"]/lowest), base = 2)
    results[i, "NKratio"] <- log((results[i, "NKmean"]/lowest), base = 2)
    results[i, "NOratio"] <- log((results[i, "NOmean"]/lowest), base = 2)}
}
  

#Nfur enriched clusters (can be adjusted to fit other species as well)
for (i in 1:nrow(results))
{
  if (is.na(results[i, "significance"])) next()
  else if (is.na(results[i, "NFratio"]) == FALSE && results[i,"NFratio"] > 3)
    print(results[i, "cluster"])
}

#relevantClusters
relevantClusters <- c() 
for (i in 1:nrow(results))
{
  if (is.na(results[i, "significance"])) next()
  else if (results[i, "significance"] == 1 && is.na(results[i, "NFratio"]) == FALSE)
    relevantClusters <- c(relevantClusters, i)
  }

#graph
colors <- c("Nfur" = "blue", "Nkad" = "green4", "Nort" = "magenta")
p4 <- ggplot(results, aes(x = cluster)) +
  geom_vline(xintercept = c(seq(from=5, to=500, by=5)), color = "white", size = 0.5) +
  geom_point(aes(y = NFratio, color = "Nfur")) +
  geom_point(aes(y = NKratio, color = "Nkad")) +
  geom_point(aes(y = NOratio, color = "Nort")) +
  labs(x = "RepeatExplorer2 cluster", y = "log2 sample/lowest mean read count", color = "Species") + scale_colour_manual(values = colors)
p4 <- p4 + theme_grey () + theme(axis.title = element_text(size = 13)) + theme(axis.text = element_text(size = 13))
p4 <- p4 + scale_x_continuous(breaks = seq(0, 500, by = 20))
p4
  

#output
write.csv(results, file = "Results_killifish_comparative1.csv")
dev.off( )
pdf("Plot_killifish_comparative1.pdf", width = 18)
p4
dev.off( )

#removing all objects and releasing memory
rm(list = ls())
gc()