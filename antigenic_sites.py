'''
Additional function for filtering only the frequent sets with at least one mutations in the 
antigenic site of HA/NA.

References:
Min-Shi Lee et al, 2004	
Colman et al, 1983

====
HA 

epitope A            122, 124, 126, 130, 131, 132, 133, 135, 137, 140, 142, 143, 144, 145, 146, 150, 152, 168
epitope B            128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198
epitope C            44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312
epitope D            96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248
epitope E            57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265

(total = 131 aa)
===
NA

epitope I            328 - 336
epitope II           339 - 347
epitope III          367 - 370
epitope IV           400 - 403
epitope V            431 - 434
epitope VI           197 - 199
epitope VII          153

'''

ep_A = [122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 142, 143, 144, 145, 146, 150, 152, 168]
ep_B = [128, 129, 155, 156, 157, 158, 159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198]
ep_C = [44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312]
ep_D = [96, 102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230, 238, 240, 242, 244, 246, 247, 248]
ep_E = [57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265]

ep_I = [*range(328,337)]
ep_II = [*range(339,348)]
ep_III = [*range(367,371)]
ep_IV = [*range(400,404)]
ep_V = [*range(431,435)]
ep_VI = [*range(197,200)]
ep_VII = [153]

all_sites_ha = ep_A + ep_B + ep_C + ep_D + ep_E
all_sites_na = ep_I + ep_II + ep_III + ep_IV + ep_V + ep_VI + ep_VII

antigenic_sites = frozenset(["ha_" + str(aa_position) for aa_position in all_sites_ha] + ["na_" + str(aa_position) for aa_position in all_sites_na])

transition_sites = frozenset(["ha_" + aa_position for aa_position in ("145", "155", "156", "158", "159", "189", "193")])

def contain_antigenic_site(itemset):
	#return True if at least one antigenic site is in the frequent itemset
	clean_itemset = frozenset([mutation[:3] + mutation[4:-1] for mutation in itemset])
	z=frozenset.intersection(clean_itemset, antigenic_sites)
	return True if len(z) != 0 else False

def contain_transition_site(itemset):
	#return True if at least one transition site is in the frequent itemset
	clean_itemset = frozenset([mutation[:3] + mutation[4:-1] for mutation in itemset])
	z=frozenset.intersection(clean_itemset, transition_sites)
	return True if len(z) != 0 else False
	


