import os
from subprocess import call
from shutil import copyfileobj

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
colourmap = colors.ListedColormap([[0.0,  0.0,  0.0],
								   [0.1,  0.6,  0.25],
	 							   [0.8,  0.1,  0.1],
	 							   [1.0,  0.7,  0.4],
	 							   [0.65, 0.85, 0.4]])


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
		X[i, :len(record.seq)] = [RNA.get(nuc, 0) for nuc in record.seq] 

	return np.array(X)



















# Read an alignment
def readalignment(filename, extensions = extensions) :

	# Check the file exists
	if not os.path.isfile(filename) :
		print "%s not found." % filename
		return



	# Check the file is of an interpretable filetype
	if os.path.splitext(filename)[1][1:] not in extensions :
		print "GenePy currently supports the following extensions :"
		print "\n- ".join(extensions.keys())
		return

	

	# Read in the records
	X = []

	with open(filename, "rU") as f :
		for record in SeqIO.parse(f, extensions[os.path.splitext(filename)[1][1:]]) :
			X.append(record.upper())




	# Return the alignment as a list of sequences
	return X

















def align(filename, force, threads, full, full_iter, it, auto) :		

	# If the data isn't on disk already
	if filename == "genepy.fasta" :
		print "genepy !"# TODO : Write to disk !

	else :
		# Generate flags :
		command = ["clustalo", "-v", "--outfmt=phy"]

		# Input file
		command.append("-i")
		command.append(filename)

		# Output file
		command.append("-o")
		command.append("temp_genepy.phy")

		# Force overwrite
		if force :
			command.append("--force")

		# Limit threads
		if threads :
			command.append("--threads")
			command.append(threads)

		# Full distance matrix
		if full :
			command.append("--full")

		# Full distance matrix during iteration only
		if full_iter :
			command.append("--full-iter")

		# Iteration
		if it :
			command.append("--iter")
			command.append(it)

		if not (it or full or full_iter) :
			command.append("--auto")


		# Call ClustalO
		print " ".join(command)
		call(command)




		# Determine number of lines in file
		with open("temp_genepy.phy", "r") as infile :
			for linecount, temp in enumerate(infile) :
				pass


		with open(os.path.splitext(filename)[0] + "_aligned_genepy.phy", "w") as outfile, open("temp_genepy.phy", "r") as infile :

			# The number of lines to change ( sequence number )
			l1 = infile.readline()
			N = int(l1.split(" ")[1])

			# Drop header in out-file
			outfile.write(l1)

			# Now write the next N lines, adding a space after the sequence name
			for i in range(N) :
				line = infile.readline()
				outfile.write(line[:10] + " " + line[10:])
			
			# Copy the rest of the file as-is
			copyfileobj(infile, outfile)

		
		os.remove("temp_genepy.phy")
		print "File rewritten as PhyML-useable input to %s" % (os.path.splitext(filename)[0] + "_aligned_genepy.phy")




		# Rewrite the file, as ClustalO output fails in PhyML
	"""
	s = readalignment(filename.split(".")[0] + "_aligned_genepy.phy")

	f = open(filename.split(".")[0] + "_aligned_genepy.phy", "w")
	SeqIO.write(s, f, "phylip")
	f.close()
	"""

















def calcstats(seq) :

	stats = {}
	stats["A"] = []
	stats["C"] = []
	stats["G"] = []
	stats["T"] = []
	stats["transition"] = np.zeros((4, 4))

	#stats["lengths"] = 

	for a in seq :
		A = a.seq.count("A")
		C = a.seq.count("C")
		G = a.seq.count("G")
		T = a.seq.count("T")
		ACGT = float(A + C + G + T)

		stats["A"].append(A / ACGT)
		stats["C"].append(C / ACGT)
		stats["G"].append(G / ACGT)
		stats["T"].append(T / ACGT)

		for i, base1 in enumerate("ACGT") :
			for j, base2 in enumerate("ACGT") :
				stats["transition"][i,j] += a.seq.count(base1 + base2) / float(len(a.seq)-1)


	stats["A"] = np.array(stats["A"])
	stats["C"] = np.array(stats["C"])
	stats["G"] = np.array(stats["G"])
	stats["T"] = np.array(stats["T"])

	stats["transition"] /= np.sum(stats["transition"])

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

	command = ["phyml", "-d", "nt", "-m", "GTR", "-v", "e"]

	# Input file
	command.append("-i")
	command.append(os.path.splitext(filename)[0] + "_aligned_genepy.phy")

	# Nucleotide frequencies
	command.append("-f")
	if nucleotide_frequency == "empirical" :
		command.append("e")
	elif nucleotide_frequency == "max_likelihood" :
		command.append("m")
	else :
		print "WARNING : Unrecognised option for nucleotide_frequency; setting to empirical."
		command.append("e")

	# Bootstrapping
	command.append("-b")
	command.append(str(bootstrap))

	# Search algorithm
	command.append("-s")
	if search_algorithm == "SPR" :
		command.append("SPR")
	elif search_algorithm == "NNI" :
		command.append("NNI")
	elif search_algorithm == "BEST" :
		command.append("BEST")
	else :
		print "WARNING : Unrecognised option for search_algorithm; setting to BEST."
		command.append("BEST")


	print " ".join(command)


	if bootstrap > 1 :
		try :
			command.insert(0, "8")
			command.insert(0, "-np")
			command.insert(0, "mpirun")
			call(command)
		except OSError :
			print "MPI not detected; running non-parallelised reconstruction."
			call(command[3:])

	else :
		call(command)























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























def dropempties(alignment, fraction) :

	def density(seq) :
		known = seq.seq.count("A") + \
				seq.seq.count("C") + \
				seq.seq.count("G") + \
				seq.seq.count("T")

		return known / float(len(seq))


	# Count ACGT
	return [seq for seq in alignment if density(seq) > fraction]


























