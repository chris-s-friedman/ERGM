#Command 1: Install or update statnet
install.packages("statnet")
install.packages("ergmharris")

#Command 2: Load statnet
#this will need to be done every time R is opened to use statnet
library(statnet)
library(ergmharris)

#Command 3: Accessing data available in R
data()

#Command 4: Load the local health department R data
data("lhds")

#Command 5: Check that the data loaded properly
lhds

#Command 6: See a summary of the network and attributes
summary(lhds)

#Command 7: Visualizing the network with nodes colored by state and HIV
#program attribute (Figure 3.1)
#dev.off( ) results in a refreshed graphic window with default settings for
#the new graphic
#palette assigns the default colors using in subsequent graphics

dev.off()
par(mfrow = c(2, 1), mar = c(0,0,1,0))
palette(gray.colors(3, 0, 1))
plot(lhds, vertex.col = "state", main = "State")
plot(lhds, vertex.col = "hivscreen", main = "HIV Screening Programs")

#Command 8: Visualize the largest component with nodes colored by HIV 
#screening programming 
#First syntax block identifies largest component, assigns to object, and 
#saves object as network 
#Second syntax block creates a subset of the HIV attribute vector and assigns 
#to network object 
#Third syntax block plots the new network with the HIV attribute vector, 
#adding a legend (Figure 3.2)

lhdscomp <- component.largest(lhds)
lhdssmall <- lhds[lhdscomp, lhdscomp]
smallnet <- as.network(lhdssmall, directed = FALSE)

hivscreen <- get.vertex.attribute(lhds, "hivscreen")
subhiv <- as.vector(subset(hivscreen, lhdscomp != "FALSE"))
smallnet %v% "hivscreen" <- subhiv

dev.off()
par(mar = c(.2, .2, 1, .2))
palette(gray.colors(2, 0, 1))
plot(smallnet, vertex.col = "hivscreen", main = "HIV Screening Programs")
legend("bottomleft", c("Y", "N"), pch = c(21, 19))

#Command 9: Using node size to represent a continuous network measure 
#R-Tip: Use the up and down arrows to scroll through syntax you've already 
#used to reuse it

numties <- degree(lhds, gmode = "graph")
lhds %v% "numties" <-  numties

dev.off()
palette(gray.colors(2, 0, 1))

par(mar = c(0,0,0,0))
plot(lhds, vertex.col = "hivscreen", vertex.cex = "numties")

# Commaqnd 10: Recoding degree and plotting the network (Figure 3.3)
# par allows editing of plot window features; mar is used to specifiy border
# size.

par(mar = c(0,0,0,0))
plot(lhds, vertex.col = "hivscreen", vertex.cex = numties/6)

# Command 11 Degree and triangle distributions for the LHD network
mean(degree(lhds, gmode = "graph"))
sd(degree(lhds, gmode = "graph"))
table(degree(lhds, gmode = "graph"))
triad.census(lhds, mode = "graph")

# Command 12: Create a random network the same size and density as LHD network
# set.seed allows simulations and statistical processes to be replicated
# exactly. 
# Use goodness-of-fit procedures to get ESP and DSP values (figure 3.4)

set.seed(2)
randomg <- rgraph(1283, 1, tprob = .0033, mode = "graph")
randomnet <- as.network(randomg, directed = FALSE)

summary(randomnet)

nullrandom <- ergm(randomnet ~ edges)
gof_nullrandom <- gof(nullrandom, GOF = 
                      ~ degree + distance + espartners + dspartners, 
                      control = control.gof.ergm(seed = 567))

lhdforesp <- ergm(lhds ~ edges)
gof_lhdforesp <- gof(lhdforesp, GOF = 
                     ~ degree + distance + espartners + dspartners, 
                     control = control.gof.ergm(seed = 567))

dev.off()

par(mfrow = c(3,2), mai = c(.8, .8, .2, .2))

hist(degree(lhds, gmode = "graph"), 
     breaks = max(degree(lhds)) / 2, 
     main = "", 
     xlab = "Degree ( LHD )", 
     xlim = range(0:15),
     ylim = range(0:500))

hist(degree(randomnet, gmode = "graph"), 
     breaks = max(degree(randomnet)) / 2, 
     main = "", 
     xlab = "Degree ( Random)", 
     xlim = range(0:15), 
     ylim = range(0:500))

barplot(gof_lhdforesp$obs.dspart[2:5], 
        main = "", 
        ylab = "Frequency", 
        xlab = "DSP ( LHD )", 
        xlim = range(0:5), 
        ylim = range(0:15000))

barplot(gof_nullrandom$obs.dspart[2:5], 
        main = "",
        xlab = "DSP (random network)", 
        xlim = range(0:5), 
        ylim = range(0:15000))

barplot(gof_lhdforesp$obs.espart[2:5], 
        main = "", 
        ylab = "Frequency", 
        xlab = "ESP ( LHD )", 
        xlim = range(0:5), 
        ylim = range(0:800))

barplot(gof_nullrandom$obs.espart[2:5], 
        main = "", 
        xlab = "ESP ( random network )", 
        xlim = range(0:5), 
        ylim = range(0:800))

# Command 13: Mixing matrices for HIV screening, nutrition programs, and years
# (table 3.2)

mixingmatrix(lhds, "state")
mixingm