
library(FuncBlocks)

gene_ids = c("Ank3","Cdkn1b","Cog8","Esd","Mtmr1","Pcf11","Tceb2","Yme1l1","Actn1","Cul5","Eif1a","Eif4h","Gramd1b","Kbtbd2","Ndufs2","Park2","Psmd3","Rab6a","Slc35e1","Sos1","Srsf2","Tmeff1","Zfp385b","Eif5a","Nus1","Pbrm1","Rap2c","Smarca5","Usp16")
genes = rep(1, length(gene_ids))
names(genes) = gene_ids

bg_gene_ids = c('Arsi', 'Mapk4', 'Papola', 'Tfrc', 'Bak1', 'Fopnl', 'Mus81', 'Opa3', 'Npcd')
bg_genes = rep(0, length(bg_gene_ids))
names(bg_genes) = bg_gene_ids

genes = c(genes,bg_genes)

go_mouse = go_enrich(genes, ref_genome="grcm38", n_randset=10)

##########
library(FuncBlocks)

genes.a = readRDS('test_go_func.genes.a.RDS')
genes.b = readRDS('test_go_func.genes.b.RDS')

table(genes.a)
table(genes.b)

setdiff(genes.a[genes.a == 1], genes.b[genes.b == 1])

res1 = go_enrich(genes.a, ref_genome="grcm38", n_randsets=10)
ran_out = read.table("tmp/randset_out", skip=2, header=T)
res2 = go_enrich(genes.b, ref_genome="grcm38", n_randsets=10)
ran_out2 = read.table("tmp/randset_out", skip=2, header=T)
ran_out = ran_out[,colnames(ran_out) %in% colnames(ran_out2)]
all(colnames(ran_out) == colnames(ran_out2))
all(ran_out[2,]==ran_out2[2,])

# from res2
Rap2c
Slc35e1
Park2
Psmd3
Tmeff1
Esd
Smarca5
Srsf2
Ank3
Ndufs2
Yme1l1
Mtmr1
Gramd1b
Rab6a
Cdkn1b
Pbrm1
Usp16
Cog8
Kbtbd2
Pcf11
Eif4h
Sos1
Nus1
Tceb2


# from res1
	Zfp385b
Rap2c
	Eif5a
	Cul5
Slc35e1
Park2
Psmd3
Tmeff1
Esd
Smarca5
Srsf2
Ank3
Ndufs2
Yme1l1
Mtmr1
Gramd1b
Rab6a
Cdkn1b
Pbrm1
Usp16
Cog8
Kbtbd2
Pcf11
Eif4h
Sos1
Nus1
Tceb2
	Actn1

