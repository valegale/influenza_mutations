import pandas as pd
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth

def frequent_itemsets_fpgrowth(transactions, min_supp = 0.2, max_len = 8):
	""" 
	Using the function fpgrowth from the mlxtend package to compute the set of 
	frequent itemsets. 
	"""
	te = TransactionEncoder()
	te_ary = te.fit(transactions).transform(transactions)
	df = pd.DataFrame(te_ary, columns=te.columns_)
	return fpgrowth(df, min_support=min_supp, use_colnames=True, max_len = max_len)
	