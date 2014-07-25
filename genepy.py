import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

try :
	import seaborn
	sns = True
except :
	sns = False

try : 
	from sklearn.neighbors import KernelDensity
	skl = True
except :
	skl = False














# RNA / DNA alphabet
RNA = { "-" : 0,
		"A" : 1,
		"C" : 2,
		"G" : 3,
		"T" : 4,
		"U" : 4
}



# Interpretable filetypes for ClustalO and BioPython
extensions = { "aln" : 	"clustal",
			   "clu" : 	"clustal",
			   "fst" : 	"fasta",
			   "fasta" : "fasta",
			   "phy" : 	"phylip",
			   "gb"  : 	"genbank",
			   "gbk" : 	"genbank",
			   "st"  : 	"stockholm"
}



# Colour map and plot colour
colourmap = colors.ListedColormap([[0, 0, 0],
								   [0, 0, 1],
	 							   [1, 0, 0],
	 							   [0, 1, 0],
	 							   [1, 1, 1]])


c = (0.7686274509803922, 0.3058823529411765, 0.3215686274509804)














# Displays the alignment as an image,
# with C and G as hot colours,
# with A and T as cold colours,
# and with "-" or unknown as black
def showalignment(obj, colourmap = colourmap) :

	plt.imshow(obj, aspect = "auto", cmap = colourmap, interpolation = "nearest")
	plt.grid(False)
	plt.show()

















# Return a numpy array representing the alignment
def alignmentarray(alignment, length = None, RNA = RNA) :

	if length is None :
		length = len(alignment[0].seq)


	X = np.zeros( (len(alignment), length), dtype = "int8" )

	# Convert
	for i, record in enumerate(alignment) :
		X[i, :len(record.seq)] = [RNA.setdefault(nuc, 0) for nuc in record.seq] 

	return np.array(X)



















# Read an alignment
def readalignment(filename, extensions = extensions) :

	# Check the file exists
	if not os.path.isfile(filename) :
		print "%s not found." % filename
		return



	# Check the file is of an interpretable filetype
	if filename.split(".")[1] not in extensions.keys() :
		print "GenePy currently supports the following extensions :"
		print "\n- ".join(extensions.keys())
		return

	

	# Read in the records
	X = []

	f = open(filename, "rU")
	for record in SeqIO.parse(f, extensions[filename.split(".")[1]]) :
		X.append(record.upper())

	f.close()



	# Return the alignment as a list of sequences
	return X

















def align(filename, force, threads, full, full_iter, iter, auto) :		

	# If the data isn't on disk already
	if filename == "genepy.fasta" :
		print "genepy !"# TODO : Write to disk !

	else :
		# Generate flags
		command = "clustalo -i %s -v -o %s_aligned_genepy.phy --outfmt=phy --wrap=60" % (filename, filename.split(".")[0])
		command += " --force" if force else ""
		command += " --threads %d" % threads if threads else ""
		command += " --full" if full else ""
		command += " --full-iter" if full_iter else ""
		command += " --iter %d" % iter if iter else ""
		command += " --auto" if not (iter or full or full_iter) else ""


		# Call ClustalO
		os.system(command)


















def calcstats(seq) :

	stats = {}
	stats["A"] = []
	stats["C"] = []
	stats["G"] = []
	stats["T"] = []
	stats["transition"] = np.zeros((4, 4))

	for a in seq :
		stats["A"].append(a.seq.count("A") / float(len(a.seq)))
		stats["C"].append(a.seq.count("C") / float(len(a.seq)))
		stats["G"].append(a.seq.count("G") / float(len(a.seq)))
		stats["T"].append(a.seq.count("T") / float(len(a.seq)))

		for i, base1 in enumerate("ACGT") :
			for j, base2 in enumerate("ACGT") :
				stats["transition"][i,j] += a.seq.count(base1 + base2) / float(len(a.seq)-1)


	stats["A"] = np.array(stats["A"])
	stats["C"] = np.array(stats["C"])
	stats["G"] = np.array(stats["G"])
	stats["T"] = np.array(stats["T"])

	stats["transition"] /= len(seq)

	return stats




















def stats(s) :

	frequencies = [s["A"].mean(),
				   s["C"].mean(),
				   s["G"].mean(),
				   s["T"].mean()]

	# Nucleotide frequencies
	fig = plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=2)
	plt.bar(range(4), frequencies, width=0.9)
	fig.set_xticks(np.arange(0.45, 4.45))
	fig.set_xticklabels(("A", "C", "G", "T"))
	plt.title("Nucleotide Frequencies")


	# Nucleotide frequency distributions
	if skl :
		x = np.linspace(0, 1, 100)[:, np.newaxis]
		yA = KernelDensity(bandwidth=0.005).fit(s["A"][:, np.newaxis])
		yC = KernelDensity(bandwidth=0.005).fit(s["C"][:, np.newaxis])
		yG = KernelDensity(bandwidth=0.005).fit(s["G"][:, np.newaxis])
		yT = KernelDensity(bandwidth=0.005).fit(s["T"][:, np.newaxis])

	plt.subplot2grid((3, 4), (2, 0))
	if skl :
		plt.plot(x, np.exp(yA.score_samples(x)), lw=3, c=c)
		plt.fill_between(x.squeeze(), np.exp(yA.score_samples(x)), color=c, alpha=0.5)
	else :
		plt.hist(s["A"], normed = True, alpha = 0.7)
	plt.title("Freq. A")
	plt.xlim([0, 1])

	plt.subplot2grid((3, 4), (2, 1))
	if skl :
		plt.plot(x, np.exp(yC.score_samples(x)), lw=3, c=c)
		plt.fill_between(x.squeeze(), np.exp(yC.score_samples(x)), color=c, alpha=0.5)
	else :
		plt.hist(s["C"], normed = True, alpha = 0.7)
	plt.title("Freq. C")
	plt.xlim([0, 1])

	plt.subplot2grid((3, 4), (2, 2))
	if skl :
		plt.plot(x, np.exp(yG.score_samples(x)), lw=3, c=c)
		plt.fill_between(x.squeeze(), np.exp(yG.score_samples(x)), color=c, alpha=0.5)
	else :
		plt.hist(s["G"], normed = True, alpha = 0.7)
	plt.title("Freq. G")
	plt.xlim([0, 1])

	plt.subplot2grid((3, 4), (2, 3))
	if skl :
		plt.plot(x, np.exp(yT.score_samples(x)), lw=3, c=c)
		plt.fill_between(x.squeeze(), np.exp(yT.score_samples(x)), color=c, alpha=0.5)
	else :
		plt.hist(s["T"], normed = True, alpha = 0.7)
	plt.title("Freq. T")
	plt.xlim([0, 1])

	# Transition Matrix
	plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=2)
	plt.imshow(s["transition"], interpolation="nearest", cmap="hot")
	plt.colorbar()
	plt.title("Transition Matrix")
	plt.xticks([0, 1, 2, 3], ["A", "C", "G", "T"])
	plt.yticks([0, 1, 2, 3], ["A", "C", "G", "T"])
	plt.grid(False)

	plt.tight_layout()
	plt.show()


















def phylotree(filename, nucleotide_frequency, bootstrap, search_algorithm) :

	command = "phyml -d nt -m GTR -v e -i %s" % (filename.split(".")[0] + "_aligned_genepy.phy")
	command += " -f e" if nucleotide_frequency == "empirical" else " -f m"
	command += " -b %d" % bootstrap


	if search_algorithm == "SPR" :
		command += " -s SPR"
	elif search_algorithm == "NNI" :
		command += " -s NNI"
	else :
		command += " -s BEST"


	if bootstrap > 1 :
		try :
			os.system("mpirun -np 8 %s" % command)
		except :
			os.system(command)

	else :
		os.system(command)























def trimalignment(alignment, array = None, left = None, right = None) :

	if array is not None :
		array[np.where(array > 0)] = 1
		density = np.sum(array, axis=0)
		# Currently, no auto-guessing. Soon !


	else :
		X = []

		for seq in alignment :
			X.append(seq[left:right])

		return X



