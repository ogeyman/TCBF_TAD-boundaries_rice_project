#dna_te
dna_te_cons <- read.table("dna_te_cons.txt",header=F)
dna_te_non <- read.table("dna_te_non.txt",header=F)

colnames(dna_te_cons) <- c("chr","start","end","val")
dna_te_cons$conservation <- "Domesticated"
colnames(dna_te_non) <- c("chr","start","end","val")
dna_te_non$conservation <- "Wild"

dna_te <- rbind(dna_te_cons,dna_te_non)
dna_te <- dna_te[,c("val","conservation")]

write.csv(dna_te,file="dna_te_cat1_combined.csv")

#gc_content

gc_cons <- read.table("gc_cons.txt",header=F)
gc_non <- read.table("gc_non.txt",header=F)

colnames(gc_cons) <- c("chr","start","end","val")
gc_cons$conservation <- "Domesticated"
colnames(gc_non) <- c("chr","start","end","val")
gc_non$conservation <- "Wild"

gc <- rbind(gc_cons,gc_non)
gc <- dna_te[,c("val","conservation")]

write.csv(gc,file="gc-content_cat1_combined.csv")