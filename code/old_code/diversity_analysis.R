#LR Colon analysis
#Kaitlin 1 24

#bring in files
meta_file <- read.table(file='data/raw/kws_metadata.tsv', header = T)
#simpson <- read.table(file='kws.an.0.03.subsample.groups.summary', header = T)
nmds <- read.table(file='data/mothur/kws_final.an.thetayc.0.03.lt.ave.nmds.axes', header = T)
simps <- read.table(file='data/mothur/kws_final.an.groups.summary', header = T)
#fullshan <- merge(shannon, metadata)
shared_file <- read.table(file='data/mothur/kws_final.an.shared', sep = '\t', header=T, row.names=2)

#load niel's script for properly reading in thetayc distances
source(file = 'code/read.dist.R')
#tyc <-read.dist(file='kws.an.thetayc.0.03.lt.dist', input = "lt", make.square=F, diag=NA)

#nmds analysis
#first separate by organ
left_mucosa <- nmds[grep('LB', nmds$group), c(2,3)]
right_mucosa <- nmds[grep('RB', nmds$group), c(2,3)]
left_stool <- nmds[grep('LS', nmds$group), c(2,3)]
right_stool <- nmds[grep('RS', nmds$group), c(2,3)]
spon_stool <- nmds[grep('SS', nmds$group), c(2,3)]

nmds <- merge(nmds, meta_file)
#nmds plots for all samples
ggplot(nmds, aes(x=axis1, y=axis2)) +geom_point(aes(color=side))
ggplot(nmds, aes(x=axis1, y=axis2)) +geom_point(aes(color=site))

#subset shared and then do distance, nmds to get comparisons of just those samples 
#could name this shared subsetting thing a function also 

#fix this function 

separate_nmds <- function(m_file, colname, observation, s_file){
  temp <- m_file$group[m_file$colname=='observation']
  tempshared <- s_file[row.names(s_file) %in% temp]
  tempshared <- add_rownames(tempshared, var = 'group')
  #tempshared <- tempshared[-2] #i dont think we need this 
  tempshared <- cbind(label=0.03, tempshared[1:ncol(tempshared)])
  #fix that, how to get variable into stringname in function 
  write.table(tempshared, file='data/process/observation.shared', sep='\t', quote=F, row.names=F)
}

separate_nmds(meta_file, side, 'left', shared_file)


left <- meta_file$group[meta_file$side=='left']
leftshared <- shared_file[row.names(shared_file) %in% left,]
leftshared <- add_rownames(leftshared, var = "group")
leftshared <- leftshared[-2] #? maybe dont do this
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

all <- meta_file$group[meta_file$side=='left' | meta_file$side=='right']
allshared <- shared_file[row.names(shared_file) %in% all,]
allshared <- add_rownames(allshared, var = "group")
allshared <- allshared[-2]
allshared <- cbind(label=0.03, allshared[1:ncol(allshared)])

write.table(leftshared, file='data/process/leftshared.shared', sep='\t', quote=F, row.names=F)
write.table(rightshared, file='data/process/rightshared.shared', sep='\t', quote=F, row.names=F)
write.table(mucosashared, file='data/process/mucosashared.shared', sep='\t', quote=F, row.names=F)
write.table(stoolshared, file='data/process/stoolshared.shared', sep='\t', quote=F, row.names=F)
write.table(allshared, file='data/process/allshared.shared', sep='\t', quote=F, row.names=F)

#run all in mothur, dist.shared and then nmds. import back in. plot. dance! 

leftnmds <- read.table(file='data/process/leftshared.thetayc.0.03.lt.ave.nmds.axes', sep='\t', header=T)
leftnmds <- merge(leftnmds, meta_file)
ggplot(leftnmds, aes(x=axis1, y=axis2)) + geom_point(aes(color=as.factor(patient), shape=site, size = 1.5)) +theme_bw() + ggtitle("Left colon comparison")

rightnmds <- read.table(file='data/process/rightshared.thetayc.0.03.lt.ave.nmds.axes', sep='\t', header=T)
rightnmds <- merge(rightnmds, meta_file)
ggplot(rightnmds, aes(x=axis1, y=axis2)) + geom_point(aes(color=as.factor(patient), shape=site, size = 1.5)) +theme_bw() + ggtitle("Right colon comparison")

mucosanmds <- read.table(file='data/process/mucosashared.thetayc.0.03.lt.ave.nmds.axes', sep='\t', header=T)
mucosanmds <- merge(mucosanmds, meta_file)
ggplot(mucosanmds, aes(x=axis1, y=axis2)) + geom_point(aes(color=as.factor(patient), shape=side, size = 1.5)) +theme_bw() + ggtitle("Mucosa samples comparison")

stoolnmds <- read.table(file='data/process/stoolshared.thetayc.0.03.lt.ave.nmds.axes', sep='\t', header=T)
stoolnmds <- merge(stoolnmds, meta_file)
ggplot(stoolnmds, aes(x=axis1, y=axis2)) + geom_point(aes(color=as.factor(patient), shape=side, size = 1.5)) +theme_bw() + ggtitle("Stool samples comparison")

allnmds <- read.table(file='data/process/allshared.thetayc.0.03.lt.nmds.axes', sep = '\t', header = T)
allnmds <- merge(allnmds, meta_file)
ggplot(allnmds, aes(axis1, axis2, group=patient)) +
  geom_point(aes(color=side), size = 2.5) +
  facet_wrap(~site, labeller=labeller(site = c(mucosa = "Muocsa", stool = "Lumen"))) +theme_bw() +
  geom_line(size=0.25, color='grey') +
  theme(legend.position = c(0.95, 0.92), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(size= 12), axis.title= element_text(size=14), 
        strip.text=element_text(size=14), legend.text=element_text(size=12)) 

alljclass <- read.table(file='data/process/allshared.jclass.0.03.lt.nmds.axes', sep = '\t', header=T)
alljclass <- merge(alljclass, meta_file)
ggplot(alljclass, aes(axis1, axis2, group=patient)) +
  geom_point(aes(color=side), size = 2.5) +
  facet_wrap(~site, labeller=labeller(site = c(mucosa = "Muocsa", stool = "Lumen"))) +theme_bw() +
  geom_line(size=0.25, color='grey') +
  theme(legend.position = c(0.95, 0.92), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(size= 12), axis.title= element_text(size=14), 
        strip.text=element_text(size=14), legend.text=element_text(size=12)) 



#nmds plot of distances faceted by location
ggplot(nmds, aes(axis1, axis2)) +geom_point(aes(color=as.factor(patient), shape=side)) + facet_wrap(~location) +theme_bw()

ggplot(leftnmds, aes(axis1, axis2)) +geom_point() +theme_bw() +ggtitle("Left colon comparison")

p <- ggplot(nmds, aes(axis1, axis2)) +geom_point(aes(color=as.factor(patient), shape=side)) + facet_wrap(~site) +theme_bw()
p <- p + scale_color_discrete(name="Patient")

noexit <- subset(nmds, side == 'left' | side == 'right')
ggplot(noexit, aes(axis1, axis2)) +geom_point(aes(color=as.factor(patient), shape=as.factor(patient)), size = 2.5) +
                                                facet_wrap(~site) +theme_bw() + scale_color_discrete(name="Patient") + scale_shape_manual(values=seq(0,20))

#testing plot options - two colors for side, each patient a shape
ggplot(noexit, aes(axis1, axis2)) +geom_point(aes(color=side, shape=as.factor(patient)), size = 3, stroke = 1.5) +
  facet_wrap(~site) +theme_bw() + scale_color_discrete(name="Patient") + scale_shape_manual(values=seq(0,20))

#testing plot - each patient a letter 
letters <- c(seq(80,90), seq(65,75))
ggplot(noexit, aes(axis1, axis2)) +geom_point(aes(color=side, shape=as.factor(patient)), size = 3, stroke =3) +
  facet_wrap(~site) +theme_bw() + scale_color_discrete(name="Patient") + scale_shape_manual(values=letters)

#NMDS connected lines plot 
ggplot(noexit, aes(axis1, axis2, group=patient)) +
  geom_point(aes(color=side), size = 2.5) +
  facet_wrap(~site, labeller=labeller(site = c(mucosa = "Muocsa", stool = "Lumen"))) +theme_bw() +
  geom_line(size=0.25, color='grey') +
  theme(legend.position = c(0.95, 0.92), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(size= 12), axis.title= element_text(size=14), 
        strip.text=element_text(size=14), legend.text=element_text(size=12)) 




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
metasimp <- merge(meta_file, simps)

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


#plotting simpson again this time in 2017
#somewhere simpmeta gets defined and idk where
simpmeta <- merge(meta_file, simps)

simpmeta$location <- factor(simpmeta$location, c("LB","RB", "LS", "RS", "SS"))
ggplot(simpmeta, aes(x=location, y=invsimpson, group =1)) +geom_point() +geom_jitter(width=0.2) +
  theme_bw() + ylab("Inverse Simpson Diversity") +
  scale_x_discrete(labels=c("L Mucosa", "R Mucosa", "L Lumen", "R Lumen", "Stool")) +
  theme(legend.position='none', axis.title.x=element_blank()) +
  stat_summary(aes(x=location, y=invsimpson), data = simpmeta, fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", width=0.4) +
  theme(axis.text = element_text(size= 16), axis.title= element_text(size=18))

ggplot(simpmeta, aes(x=site, y=invsimpson)) +geom_boxplot()






