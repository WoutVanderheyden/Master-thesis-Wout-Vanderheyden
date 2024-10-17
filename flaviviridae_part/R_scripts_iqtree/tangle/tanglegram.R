rm(list=ls())

######################################################
##### Stuff you should change ########################
######################################################

# Stick the trees in a separate folder and put the path here
# Names on the plot will be the names of the files, minus the .tree extension
# (So best you use .tree)
treefileswd <- "/data/leuven/358/vsc35887/master_thesis/Master-thesis-Wout-Vanderheyden/flaviviridae_part/R_scripts_iqtree/tangle/data/"

# Where do you want the plot to be printed? And whats the name? (I use .png for now)
outpath <- "/data/leuven/358/vsc35887/master_thesis/Master-thesis-Wout-Vanderheyden/flaviviridae_part/R_scripts_iqtree/tangle/tanglegramE1.png"

# Dimensions in pixels. Depends on how long your trees are and how many you have.
xdim <- 2000
ydim <- 2000

# Depending on the size of your trees you might want to change the size of the x labels
xlabsize <- 4
# Same thing for the thickness of the lines connecting the tips.
# If your trees are invisible between the lines, make them thinner.
# If the lines are barely visible, make them thicker
linethickness <- 0.8 #values between 0.1 and 1.5 or something make for decent plots

# Remark: Check the appropriateness of these graphical parameters in the outputted png, 
# not in the R plot viewer because that one automatically rescales things idk.

######################################################
##### Setup ##########################################
######################################################

# Libs
library(ape)
library(ggplot2)
library(ggtree)
library(phytools)
library(dplyr)
library(khroma)
# Read tree files
treefiles <- list.files(treefileswd)
trees <- vector("list", length=length(treefiles))

#hier eventueel automatisch de bomen 'midpoint rooten'

for (t in 1:length(treefiles)){
  #trees[[t]] <- read.nexus(paste0(treefileswd, treefiles[t]))
  trees[[t]] <- read.tree(paste0(treefileswd, treefiles[t]))
  trees[[t]] <- midpoint.root(trees[[t]])
  trees[[t]][["node.label"]] <- NULL
  names(trees)[t] <- gsub(".treefile", "", treefiles[t])
}

trees <- trees[c(2,1)]

######################################################
##### My brilliant untangling code ###################
######################################################

# Just leave it as is. Can take a while for big trees though.

for (t in 1:(length(trees)-1)){
  
  print(paste0("Aligning trees ", t, " and ", (t+1)))
  
  # Take trees
  t1 <- trees[[t]]
  t2 <- trees[[(t+1)]]
  
  # Tips
  tips <- 1:length(t1$tip.label)
  
  # Nodes in both trees
  nodes <- (length(t1$tip.label)+1):(t1$Nnode+length(t1$tip.label))
  nodes1 <- data.frame(node=nodes,tips=NA)
  nodes2 <- data.frame(node=nodes,tips=NA)
  
  for (n in 1:length(nodes)){
    f1 <- sort(getDescendants(t1, nodes[n])[which(getDescendants(t1, nodes[n]) %in% tips)])
    f2 <- sort(getDescendants(t2, nodes[n])[which(getDescendants(t2, nodes[n]) %in% tips)])
    nodes1$tips[n] <- paste(f1, collapse = "-")
    nodes2$tips[n] <- paste(f2, collapse = "-")
  }
  
  # Go through each node
  for (n in nodes){
    # Node to check
    node <- n
    # Babies of tree 1 node
    babies <- nodes1$tips[which(nodes1$node==n)]
    # Check if it exists in tree 2
    correspondingnode <- nodes2$node[which(nodes2$tips == babies)]
    # If it does, check which tips are on which side of the split in both trees
    if (length(correspondingnode)>0){
      # Tree 1
      # Whole family of node
      fam1 <- as.numeric(unlist(strsplit(nodes1$tips[which(nodes1$node==n)], "-")))
      # Direct kids of node
      kids1 <- t1$edge[t1$edge[,1] == node, 2]
      # Descendants of direct kids
      upfam1 <- c()
      lowfam1 <- c()
      for (desc in fam1[which(!fam1 %in% kids1)]) {
        if (desc %in% getDescendants(t1, kids1[1])) {
          upfam1 <- c(upfam1, desc)
        } else {
          lowfam1 <- c(lowfam1, desc)
        }
      }
      # Tree 2
      # Whole family of node
      fam2 <- as.numeric(unlist(strsplit(nodes2$tips[which(nodes2$node==correspondingnode)], "-")))
      # Direct kids
      kids2 <- t2$edge[t2$edge[,1] == correspondingnode, 2]
      # Descendants of direct kids
      upfam2 <- c()
      lowfam2 <- c()
      for (desc in fam2[which(!fam2 %in% kids2)]) {
        if (desc %in% getDescendants(t2, kids2[1])) {
          upfam2 <- c(upfam2, desc)
        } else {
          lowfam2 <- c(lowfam2, desc)
        }
      }
      # If the fams are empty, they are tips, so we stick them back in
      if (is.null(upfam1)){upfam1 <- c(upfam1, kids1[1])}
      if (is.null(upfam2)){upfam2 <- c(upfam2, kids2[1])}
      if (is.null(lowfam1)){lowfam1 <- c(lowfam1, kids1[2])}
      if (is.null(lowfam2)){lowfam2 <- c(lowfam2, kids2[2])}
      # If they're flipped, rotate the tree
      if (identical(upfam1, lowfam2) & identical(lowfam1, upfam2)){
        t2 <- ape::rotate(t2, node=correspondingnode)
        print(paste0("Node ", correspondingnode, " flipped"))
      }
    }
  }
  
  # Stick em back in
  trees[[t+1]] <- t2
  
}


######################################################
##### Extract plot data ##############################
######################################################

# Make basic plot for each tree and extract the plot and the ggplot data
pdata <- vector("list", length = length(trees))
plots <- vector("list", length = length(trees))
for (t in 1:length(trees)){
  p <- ggplot(trees[[t]]) + geom_tree() + theme_tree()
  plots[[t]] <- p
  pd <- p$data
  pd$tree <- paste0("t",t)
  pdata[[t]] <- pd
}

######################################################
##### Build plotting data ############################
######################################################

# Empty label locations dataframe
lablocs <- data.frame(
  val = names(trees),
  x = NA, 
  y = 1
)

# Combine plot data of the trees and shift labels according to tree lengths
pp <- plots[[1]]
for (i in 1:length(pdata)){
  if (i == 1){
    lablocs$x[i] <- max(pdata[[i]]$x, na.rm = T)-
      (max(pdata[[i]]$x, na.rm = T)-min(pdata[[i]]$x, na.rm = T))/2
    next
  }
  pdata[[i]]$x <- pdata[[i]]$x + max(pdata[[i-1]]$x) + pdata[[i]]$x*.3
  pp <- pp + geom_tree(data=pdata[[i]])
  lablocs$x[i] <- max(pdata[[i]]$x, na.rm = T)-
    (max(pdata[[i]]$x, na.rm = T)-min(pdata[[i]]$x, na.rm = T))/2
}
dd <- bind_rows(pdata) %>% 
  filter(isTip == TRUE)
dd1 <- as.data.frame(dd)
plottree <- dd1[, c('label', 'x', 'y', 'tree')]


######################################################
##### Build plot #####################################
######################################################

# Colours (you can change this if you want)
cols <- rev(as.character(color("batlow")(length(unique(pdata[[1]]$label)))))
names(cols) <- pdata[[1]]$label[order(pdata[[1]]$y[!is.na(pdata[[1]]$label)])]

pp + geom_line(
  aes(x, y, group=label, col=label), 
  data=plottree,
  alpha = .7,
  size = linethickness
) +
  scale_color_manual(values = cols) +
  guides(col="none") + 
  annotate(
    label = lablocs$val,
    x = lablocs$x,
    y = lablocs$y,
    geom = "text",
    size= xlabsize,
    hjust = .5,
    vjust = 2,
    angle = 0,
    col="black"
  )

ggsave(outpath, 
       width = xdim,
       height = ydim,
       dpi = "retina",
       unit = "px",
       device = png
)


