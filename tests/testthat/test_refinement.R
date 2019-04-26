

test_that("hyper() - refinement hyper category test",{
    expect_true(all.equal(hyper(c(3,4,6), c(12,12,12), c(0,2,12), c(19,19,19)), c(0.0489433, 0.1366560, 0.863453), tolerance=1.5e-6))
    expect_true(all.equal(hyper(c(3,4,6), c(12,12,12), c(0,2,12), c(19,19,19), under=TRUE), c(1, 0.978307, 0.362281), tolerance=1.5e-6))

})


test_that("hyper_nodes() - refinement hyper category test for all leaves",{
	
	scores_root = structure(c(12, 17641), .Dim = 1:2, .Dimnames = list("GO:0008150", 
    NULL))
	
	# empty and non-empty leaves
	anno_nodes = structure(list(go_id = c("GO:0072205", "GO:0072205", "GO:0072205", 
	"GO:0072205", "GO:0072205", "GO:0072205", "GO:0072205", "GO:0072205", 
	"GO:0072221", "GO:0072221", "GO:0072221", "GO:0072221", "GO:0072221"
	), gene = c("AQP2", "BMP4", "CALB1", "DLG5", "PAX2", "PAX8", 
	"PKD1", "WNT7B", "CALB1", "PAX2", "PAX8", "POU3F3", "UMOD"), 
    scores = c(0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0)), row.names = 17684:17696, class = "data.frame")
	empty_nodes = c("nodeA", "nodeB")	
	pvals = hyper_nodes(anno_nodes, empty_nodes, scores_root)	
	expected = structure(list(go_id = structure(1:4, .Label = c("GO:0072205", 
	"GO:0072221", "nodeA", "nodeB"), class = "factor"), new_p = c(1.18340981614424e-05, 
	4.23125608013295e-06, 1, 1)), class = "data.frame", row.names = c(NA,-4L))	
	expect_true(all.equal(pvals, expected))
	
	# only non-empty leaves
	pvals = hyper_nodes(anno_nodes, character(0), scores_root)
	expected = structure(list(go_id = structure(1:2, .Label = c("GO:0072205", 
	"GO:0072221"), class = "factor"), new_p = c(1.18340981614424e-05, 
	4.23125608013295e-06)), class = "data.frame", row.names = c(NA, -2L))
	expect_true(all.equal(pvals, expected[1:2,]))
	
	# only empty leaves
	anno_nodes = structure(list(go_id = character(0), gene = character(0), scores = numeric(0)),
	row.names = integer(0), class = "data.frame")
	pvals = hyper_nodes(anno_nodes, empty_nodes, scores_root)
	expected = structure(list(go_id = structure(1:2, .Label = c("nodeA", "nodeB"), class = "factor"), 
	new_p = c(1, 1)), class = "data.frame", row.names = c(NA, -2L))
	expect_true(all.equal(pvals, expected))
	
})





test_that("refinement wilcoxon category test",{
	
	anno_genes_node = structure(list(gene = c("AGTR1", "BTBD3", "CACNG2", "CALB1", 
	"ENGASE", "G6PD", "GCK", "GYG1", "GYS1", "HK2", "PAX2", "PYGL", 
	"SLC2A8", "UGP2"), score = c(14L, 15L, 10L, 5L, 21L, 23L, 27L, 
	11L, 30L, 29L, 12L, 26L, 20L, 22L)), class = "data.frame", row.names = c(NA, 
	14L))

	anno_genes_root = structure(list(gene = c("AGTR1", "ANO1", "BTBD3", "CACNG2", "CALB1", 
	"ENGASE", "G6PD", "GCK", "GYG1", "GYS1", "HK2", "MTUS1", "PAX2", 
	"PYGL", "SLC2A8", "UGP2", "ZWINT"), score = c(14L, 9L, 15L, 10L, 
	5L, 21L, 23L, 27L, 11L, 30L, 29L, 13L, 12L, 26L, 20L, 22L, 28L
	)), class = "data.frame", row.names = c(NA, 17L))

	p_high_rank = wilcox(anno_genes_node, anno_genes_root)
	p_low_rank = wilcox(anno_genes_node, anno_genes_root, low=TRUE)
	
	expect_true(all.equal(c(p_high_rank, p_low_rank), c(0.329622, 0.714625), tolerance=1.5e-6))

})
