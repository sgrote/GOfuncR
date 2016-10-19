
setwd("/r1/people/steffi_grote/David")

go = "GO:0098800"

# read ontology graph
term = read.table("/r1/people/steffi_grote/R_packages/term_tables/term.txt" ,sep="\t", quote="", comment.char="", as.is=TRUE)
graph_path = read.table("/r1/people/steffi_grote/R_packages/term_tables/graph_path.txt" ,sep="\t", quote="", comment.char="", as.is=TRUE)

# go_id in ontology
go_id = term[term[,4]==go, 1]

# linked go-categories
go_related = graph_path[graph_path[,2]==go_id | graph_path[,3]==go_id,]
# what is the relation of term2 and term1?
go_related$relation = term[match(go_related[,4],term[,1]),4]
# exctract child_nodes (either '"is_a" our_GO' or '"part of" our_GO')
go_children = go_related[go_related[,2]==go_id & (go_related$relation %in% c("part_of", "is_a")), ]
# GO-IDs of child_nodes
child_gos = term[term[,1] %in% go_children[,3],4]

### check if any of the input_genes from low.genes gets annotated to any of the children
load("low.genes.RData")  
genes = low.genes

# this is what is done in the package:
	message("get GOs for genes")
	go = read.table("/mnt/expressions/miguel/neandertal_int/mart_export.txt",as.is=T,sep="\t",head=T)
	go = go[go[,2] %in% names(genes) & go[,1]!="",2:1]  
	go$value = genes[match(go[,1], names(genes))]
	# subset to GOs present in term.txt
	term = read.table("/r1/people/steffi_grote/R_packages/term_tables/term.txt" ,sep="\t", quote="", comment.char="", as.is=TRUE)
	# remove obsolete terms
	term = term[term[,5]==0,]
	go = go[go[,2] %in% term[,4],]

# and?
anno = go[go[,2] %in% child_gos,]
# aha!
# and the names
anno$go_name = term[match(anno[,2], term[,4]),2]
table(anno$go_name)
