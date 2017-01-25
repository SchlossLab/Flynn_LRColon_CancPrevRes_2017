#LR Colon analysis
#Kaitlin 1 24

#bring in files
meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
#simpson <- read.table(file='kws.an.0.03.subsample.groups.summary', header = T)
nmds <- read.table(file='data/mothur/kws_final.an.thetayc.0.03.lt.ave.nmds.axes', header = T)
simps <- read.table(file='data/mothur/kws_final.an.groups.summary', header = T)
#fullshan <- merge(shannon, metadata)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)

#load niel's script for properly , reading in thetayc distances
#source(file = 'read.dist.R')
#tyc <-read.dist(file='kws.an.thetayc.0.03.lt.dist', input = "lt", make.square=F, diag=NA)

#nmds analysis
#first separate by organ
left_mucosa <- nmds[grep('LB', nmds$group), c(2,3)]
right_mucosa <- nmds[grep('RB', nmds$group), c(2,3)]
left_stool <- nmds[grep('LS', nmds$group), c(2,3)]
right_stool <- nmds[grep('RS', nmds$group), c(2,3)]
spon_stool <- nmds[grep('SS', nmds$group), c(2,3)]

#make plot, this is probably not the best way to do it
plot(nmds$axis1, nmds$axis2)

nmds <- merge(nmds, metadata)

sharedmet <- merge(shared_file, meta_file, by.x="row.names", by.y='group')


#ok this is looking better. if can make for all sides/comparisons, then can determine if donors are closer to each other or locations. actually probably need to subset shared first, then do nmds each time 
leftside <- subset(nmds, side=='left')
rightsidenmds <- subset(nmds, side=='right')

ggplot(leftside, aes(x=axis1, y=axis2)) + geom_point(aes(color=as.factor(patient), shape=site, size = 2)) +theme_bw() + ggtitle("Left colon comparison")

ggplot(rightsidenmds, aes(x=axis1, y=axis2)) + geom_point(aes(color=as.factor(patient), shape=site, size = 2)) +theme_bw() + ggtitle("Right colon comparison")

#subset shared and then do distance, nmds to get comparisons of just those samples 

left <- meta_file$group[meta_file$side=='left']
leftshared <- shared_file[row.names(shared_file) %in% left,]
leftshared <- add_rownames(leftshared, var = "group")
leftshared <- leftshared[-2]
leftshared <- cbind(label=0.03, leftshared[1:ncol(leftshared)])

right <- meta_file$group[meta_file$side=='right']
rightshared <- shared_file[row.names(shared_file) %in% right,]
rightshared <- add_rownames(rightshared, var = "group")
rightshared <- rightshared[-2]
rightshared <- cbind(label=0.03, rightshared[1:ncol(rightshared)])

mucosa <- meta_file$group[meta_file$site=='mucosa']
mucosashared <- shared_file[row.names(shared_file) %in% mucosa,]
mucosashared <- add_rownames(mucosashared, var = "group")
mucosashared <- mucosashared[-2]
mucosashared <- cbind(label=0.03, mucosashared[1:ncol(mucosashared)])

stool <- meta_file$group[meta_file$site=='stool']
stoolshared <- shared_file[row.names(shared_file) %in% stool,]
stoolshared <- add_rownames(stoolshared, var = "group")
stoolshared <- stoolshared[-2]
stoolshared <- cbind(label=0.03, stoolshared[1:ncol(stoolshared)])


write.table(leftshared, file='data/process/leftshared.shared', sep='\t', quote=F, row.names=F)
write.table(rightshared, file='data/process/rightshared.shared', sep='\t', quote=F, row.names=F)
write.table(mucosashared, file='data/process/mucosashared.shared', sep='\t', quote=F, row.names=F)
write.table(stoolshared, file='data/process/stoolshared.shared', sep='\t', quote=F, row.names=F)

#nmds plot of distances faceted by location
ggplot(nmds, aes(axis1, axis2)) +geom_point(aes(color=patient, shape=side)) + facet_wrap(~location)

#color points by organ 
points(left_mucosa, pch=nmds$patient, col = "blue")
points(right_mucosa, pch=nmds$patient,col = "red")
points(left_stool, pch=nmds$patient, col = "green")
points(right_stool, pch=nmds$patient, col = "orange")
points(spon_stool, pch=nmds$patient, col = "purple")

legend <- c("left_mucosa", "right_mucosa", "left_stool", "right_stool", "spon_stool")
legend(x="topright", legend, col = c("blue", "red", "green", "orange", "purple"), pch=nmds$patient)

#lets do a test to compare 3 patients patients, or within one patient
#or subset just to test three patients
fullnmds <-merge(nmds, metadata)
pt6ss <- subset(fullnmds, patient == 6 | patient == 9 | patient == 15)
plot(pt6ss$axis1, pt6ss$axis2)

#test to get colors in points 
pt6 <- pt6ss[grep('6', pt6ss$patient), c(2,3)]
pt9 <- pt6ss[grep('9', pt6ss$patient), c(2,3)]
pt15 <- pt6ss[grep('15', pt6ss$patient), c(2,3)]

points(pt6, pch=16, col = "blue")
points(pt9, pch=16, col = "red")
points(pt15, pch=16, col = "green")

#or plot looking at organs again-- can't do all of this at once
left_mucosa <- pt6ss[grep('LB', pt6ss$group), c(2,3)]
right_mucosa <- pt6ss[grep('RB', pt6ss$group), c(2,3)]
left_stool <- pt6ss[grep('LS', pt6ss$group), c(2,3)]
right_stool <- pt6ss[grep('RS', pt6ss$group), c(2,3)]
spon_stool <- pt6ss[grep('SS', pt6ss$group), c(2,3)]

points(left_mucosa, pch=16, col = "blue")
points(right_mucosa, pch=16, col = "red")
points(left_stool, pch=16, col = "green")
points(right_stool, pch=16, col = "orange")
points(spon_stool, pch=16, col = "purple")

legend <- c("left_mucosa", "right_mucosa", "left_stool", "right_stool", "spon_stool")
legend(x="topright", legend, col = c("blue", "red", "green", "orange", "purple"), pch=16)


#thetayc comparison time

#set up file to do comparisons
rows <- c(row.names(tyc))
cols <- c(colnames(tyc))


meta <- metadata[metadata$group %in% colnames(tyc),]
left <- meta[meta$side=='left', 'group']

left_mucosa <- colnames(tyc)[grep('LB', colnames(tyc))]
left2left <- as.numeric(unlist(tyc[left_mucosa, left_mucosa]))

right_mucosa <- colnames(tyc)[grep('RB', colnames(tyc))]
left2right <- as.numeric(unlist(tyc[left_mucosa, right_mucosa]))
left2right <- na.omit(left2right)

right_stool <- colnames(tyc)[grep('RS', colnames(tyc))]
left_stool <- colnames(tyc)[grep('LS', colnames(tyc))]
spon_stool <- colnames(tyc)[grep('SS', colnames(tyc))]

LS2RS <- as.numeric(unlist(tyc[left_stool, right_stool]))
LS2RS <- na.omit(LS2RS)

right2right <- as.numeric(unlist(tyc[right_mucosa, right_mucosa]))
right2right <- na.omit(right2right)

LS2LS <- as.numeric(unlist(tyc[left_stool, left_stool]))
LS2LS <- na.omit(LS2LS)

LS2SS <- as.numeric(unlist(tyc[left_stool, spon_stool]))
LS2SS <- na.omit(LS2SS)

RS2SS <- as.numeric(unlist(tyc[right_stool, spon_stool]))
RS2SS <- na.omit(RS2SS)

RS2RS <- as.numeric(unlist(tyc[right_stool, right_stool]))
RS2RS <- na.omit(RS2RS)

SS2SS <- as.numeric(unlist(tyc[spon_stool, spon_stool]))
SS2SS <- na.omit(SS2SS)

#can take mean and sd and plot store and plot, or plot boxplots
labels <- c("LM/LM", "LM/RM", "RM/RM", "LS/LS", "LS/RS", "RS/RS", "LS/SS", "RS/SS", "SS/SS")
boxplot(left2left, left2right, right2right, LS2LS, LS2RS, RS2RS, LS2SS, RS2SS, SS2SS, main="dissimilarity between sample sites, all patients", ylab="theta YC distance", names=labels)

## now do simpson diversity per patient
metasimp <- merge(metadata, simps)

#want to plot diversity side by side for each patient.
#plotted with a loop

for(i in unique(metasimp$patient)){
  temp <- subset(metasimp, patient == i)
  plot(temp$location, temp$invsimpson, type = 'p', main = i, ylab="invsimpson", xlab="site")
}

#plot diversity per each site 
plot(metasimp$location, metasimp$invsimpson, main = "invsimpson diversity per site, all patients", ylab= "invsimpson", xlab="site")
plot(fullshan$location, fullshan$shannon, main= "Shannon diversity by location", ylab= "Shannon diversity", xlab="location", xaxt = "n")
axis(1, at=1:5, labels=c("left biopsy", "left stool", "right biopsy", "right stool", "spon. stool"))

justbiop <- subset(fullshan, location == "LB" | location == "RB" | location == "SS")
plot(justbiop$location, justbiop$shannon, main= "Shannon diversity by location", ylab= "Shannon diversity", xlab="location", xaxt = "n")
axis(1, at=1:5, labels=c("left biopsy", " ", "right biopsy", " ", "spon. stool"))





#Random forest results just trying these data sets a little. used mothur's classify.rf
#to compare LB/RB and LB/LS. takes some finagling of design files
#found that for LB/LS most predictive is Ecoli/Shigella (but for which side?)
#found for LB/RB most predictive is Enterococcus. But both of these have low values (under 1)
#found for RB/RS most predictive is OTU007 Bacteriodes 
#with the caveat that this is probably wrong/off/will retry with better dataset



