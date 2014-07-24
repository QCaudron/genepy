import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors




# RNA / DNA alphabet
RNA = { "-" : 0,
		"A" : 1,
		"C" : 2,
		"G" : 3,
		"T" : 4,
		"U" : 4
}



# Interpretable filetypes for ClustalO and BioPython
extensions = { "aln" : "clustal",
			   "clu" : "clustal",
			   "fst" : "fasta",
			   "phy" : "phylip",
			   "gb"  : "genbank",
			   "gbk" : "genbank",
			   "st"  : "stockholm"
}





# Colour map
colourmap = colors.ListedColormap([[0, 0, 0],
								   [0, 0, 1],
	 							   [1, 0, 0],
	 							   [0, 1, 0],
	 							   [1, 1, 1]])














# Displays the alignment as an image,
# with C and G as hot colours,
# with A and T as cold colours,
# and with "-" or unknown as black
def showalignment(obj, colourmap = colourmap) :

	# If it's an alignment array
	if type(obj) is np.ndarray :
		plt.imshow(obj, aspect = "auto", cmap = colourmap, interpolation = "nearest")
		plt.grid(False)
		plt.show()

















# Return a numpy array representing the alignment
def alignmentarray(alignment, length = None, RNA = RNA) :

	print "This might be slow !"

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
		X.append(record)

	f.close()



	# Return the alignment as a list of sequences
	return X

















def align(filename, force, threads, full, full_iter, iter, auto) :		

	# If the data isn't on disk already
	if filename == "genepy.fasta" :
		print "genepy !"# TODO : Write to disk !

	else :
		# Generate flags
		command = "clustalo -i %s -v -o %s_aligned_genepy.phy --outfmt=phy" % (filename, filename.split(".")[0])
		command += " --force" if force else ""
		command += " --threads %d" % threads if threads else ""
		command += " --full" if full else ""
		command += " --full-iter" if full_iter else ""
		command += " --iter %d" % iter if iter else ""
		command += " --auto" if not (iter or full or full_iter) else ""


		# Call ClustalO
		return(command)
