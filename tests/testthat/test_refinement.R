

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
	expected = structure(list(go_id = c("GO:0072205", "GO:0072221", "nodeA", 
	"nodeB"), new_p = c(1.18340981614424e-05, 4.23125608013295e-06, 
	1, 1)), class = "data.frame", row.names = c(NA, -4L))
	expect_true(all.equal(pvals, expected))
	
	# only non-empty leaves
	pvals = hyper_nodes(anno_nodes, character(0), scores_root)
	expected = structure(list(go_id = c("GO:0072205", "GO:0072221"), new_p = c(1.18340981614424e-05, 
	4.23125608013295e-06)), class = "data.frame", row.names = c(NA, -2L))
	expect_true(all.equal(pvals, expected[1:2,]))
	
	# only empty leaves
	anno_nodes = structure(list(go_id = character(0), gene = character(0), scores = numeric(0)),
	row.names = integer(0), class = "data.frame")
	pvals = hyper_nodes(anno_nodes, empty_nodes, scores_root)
	expected = structure(list(go_id = c("nodeA", "nodeB"), new_p = c(1, 1)), 
	class = "data.frame", row.names = c(NA, -2L))
	expect_true(all.equal(pvals, expected))
	
})




test_that("wilcox() - refinement wilcoxon category test",{
	
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



test_that("wilcox_nodes() - refinement hyper category test for all leaves",{

	scores_root = structure(list(gene = c("GYG1", "MTUS1", "NCAPG", "NGFR", "NXPH4", 
	"PAX2"), scores = c(30L, 10L, 9L, 12L, 24L, 22L)), row.names = 21:26, class = "data.frame")
	
	# empty and non-empty leaves
	anno_nodes = structure(list(go_id = c("GO:0003824", "GO:0003824", "GO:0008194"
	), gene = c("GYG1", "PAX2", "GYG1"), scores = c(30L, 22L, 30L
	), term_id = c(3010L, 3010L, 6655L), root_id = c("GO:0003674", 
	"GO:0003674", "GO:0003674")), row.names = c(27L, 28L, 67L), class = "data.frame")
	empty_nodes = c("nodeA", "nodeB")
	pvals = wilcox_nodes(anno_nodes, empty_nodes, scores_root)
	expected = structure(list(go_id = c("GO:0003824", "GO:0008194", "nodeA", 
	"nodeB"), new_p = c(0.123579986605954, 0.120783293484486, 1, 
	1)), row.names = c(NA, -4L), class = "data.frame")
	expect_true(all.equal(pvals, expected))
	
	# only non-empty leaves
	pvals = wilcox_nodes(anno_nodes, character(0), scores_root)
	expected = structure(list(go_id = c("GO:0003824", "GO:0008194"), new_p = c(0.123579986605954, 
	0.120783293484486)), class = "data.frame", row.names = c(NA, -2L))
	expect_true(all.equal(pvals, expected))
	
	# only empty leaves
	anno_nodes = structure(list(go_id = character(0), gene = character(0), scores = numeric(0)),
	row.names = integer(0), class = "data.frame")
	pvals = wilcox_nodes(anno_nodes, empty_nodes, scores_root)
	expected = structure(list(go_id = c("nodeA", "nodeB"), new_p = c(1, 1)),
	class = "data.frame", row.names = c(NA, -2L))
	expect_true(all.equal(pvals, expected))
	
})




test_that("binom() - refinement binomial category test",{

	a_node = 148
	b_node = 176
	a_root = 289
	b_root = 315
	
	p_high_a = binom(a_node, b_node, a_root, b_root)
	p_high_b = binom(a_node, b_node, a_root, b_root, low=TRUE)
	
	expect_true(all.equal(c(p_high_a, p_high_b), c(0.798633, 0.234112), tolerance=1.5e-6))
	
})


test_that("conti() - contingency table category test",{

	# high A/B
	abcd = c(55, 17, 95, 31) 
	p_high_ab = conti(abcd[1], abcd[2], abcd[3], abcd[4])
	p_high_cd = conti(abcd[1], abcd[2], abcd[3], abcd[4], low=TRUE)
	expect_true(all.equal(c(p_high_ab, p_high_cd), c(0.8754846, 1), tolerance=1.5e-6))
	
	# high C/D
	abcd = c(166, 63, 320, 75)
	p_high_ab = conti(abcd[1], abcd[2], abcd[3], abcd[4])
	p_high_cd = conti(abcd[1], abcd[2], abcd[3], abcd[4], low=TRUE)
	expect_true(all.equal(c(p_high_ab, p_high_cd), c(1, 0.01340939), tolerance=1.5e-6))
	
	# Fisher
	abcd = c(29, 1, 54, 22)
	p_high_ab = conti(abcd[1], abcd[2], abcd[3], abcd[4])
	p_high_cd = conti(abcd[1], abcd[2], abcd[3], abcd[4], low=TRUE)
	expect_true(all.equal(c(p_high_ab, p_high_cd), c(0.003322334, 1), tolerance=1.5e-6))
	

})






