import cpalgorithm as cpa
import networkx as nx


def test_BE():
	# BE =========================
	G=nx.karate_club_graph()
	print("Running BE algorithm ...")
	be = cpa.BE()
	be.detect(G)
	pair_id = be.get_pair_id()
	coreness = be.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, be)

def test_MINRES():
	# MINRES =========================
	G=nx.karate_club_graph()
	print("Running MINRES algorithm ...")
	mrs = cpa.MINRES()
	mrs.detect(G)
	pair_id = mrs.get_pair_id()
	coreness = mrs.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, mrs)
	
def test_SBM():
	# SBM =========================
	G=nx.karate_club_graph()
	print("Running SBM algorithm ...")
	sbm = cpa.SBM()
	sbm.detect(G)
	pair_id = sbm.get_pair_id()
	coreness = sbm.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, sbm)
	
def test_LowRankCore():
	# LowRankCore =========================
	G=nx.karate_club_graph()
	print("Running LowRankCore algorithm ...")
	lrc = cpa.LowRankCore()
	lrc.detect(G)
	pair_id = lrc.get_pair_id()
	coreness = lrc.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, lrc)
	
def test_LapCore():
	# LapCore =========================
	G=nx.karate_club_graph()
	print("Running LapCore algorithm ...")
	lc = cpa.LapCore()
	lc.detect(G)
	pair_id = lc.get_pair_id()
	coreness = lc.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, lc)
	
def test_LapSgnCore():
	# LapSgnCore =========================
	G=nx.karate_club_graph()
	print("Running LapSgnCore algorithm ...")
	lsc = cpa.LapSgnCore()
	lsc.detect(G)
	pair_id = lsc.get_pair_id()
	coreness = lsc.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, lsc)
	
def test_Rombach():
	# Rombach =========================
	G=nx.karate_club_graph()
	print("Running Rombach's algorithm ...")
	rb = cpa.Rombach()
	rb.detect(G)
	pair_id = rb.get_pair_id()
	coreness = rb.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, rb)
	
def test_Rossa():
	# Rossa =========================
	G=nx.karate_club_graph()
	print("Running Rossa's algorithm ...")
	rs = cpa.Rossa()
	rs.detect(G)
	pair_id = rs.get_pair_id()
	coreness = rs.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, rs)
	
def test_KM_ER():
	# KM--ER =========================
	G=nx.karate_club_graph()
	print("Running KM_ER ...")
	km = cpa.KM_ER()
	km.detect(G)
	pair_id = km.get_pair_id()
	coreness = km.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, km)
	
def test_KM_config():
	# KM--config =========================
	G=nx.karate_club_graph()
	print("Running KM_config ...")
	km = cpa.KM_config()
	km.detect(G)
	pair_id = km.get_pair_id()
	coreness = km.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, km)

def test_Divisive():
	# KM--config =========================
	G=nx.karate_club_graph()
	print("Running Divisive ...")
	dv = cpa.Divisive()
	dv.detect(G)
	pair_id = dv.get_pair_id()
	coreness = dv.get_coreness()
	sig_pair_id, sig_coreness, significance, p_values = cpa.qstest(pair_id, coreness, G, dv)
