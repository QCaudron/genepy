import genepy

import numpy as np
import matplotlib.pyplot as plt
import os








# Base sequence list class
class seqset :
	"""
	Base class for a set of RNA / DNA sequences. 
	Constructor takes a filename or a list of strings
	Members :
		.len 			Number of sequences in the set
		.seq_len		np.array of lengths of sequences
		.seq			List of BioPython Seq objects
		.show()			Image representation of the sequences

	"""














	def __init__(self, source) :
		"""
		Help here ?
		"""

		# If we're reading in a sequence set from a file
		if type(source) is str :
			if os.path.isfile(source) :
				self.seq = genepy.readalignment(source)
				self.filename = source
			else :
				print "%s not found, aborting." % source

		# If we're fed a list
		elif type(source) is list :
			self.seq = [Seq(s, alphabet = generic_dna) for s in source]
			self.filename = "genepy.fasta"

		else :
			print "Expected a filename or a list of strings. Aborting."
			return


		# Generate static members
		self.update()
		

















	def __str__(self) :
		summary = \
"GenePy sequence set :\n-- %d sequences\n-- Mean length : %.01f (min %d, max %d)\n \
-- C+G content : %.03f\n-- From file : %s\n" % \
		(self.len, 
		np.array(self.seq_len).mean(), np.min(self.seq_len), np.max(self.seq_len),  
		(self.statistics["C"].mean() + self.statistics["G"].mean()),
		self.filename.split("/")[-1])

		return summary















	def __repr__(self) :
		return "GenePy seqset : %f" % self.filename

















	def update(self) :

		# Number of sequences
		self.len = len(self.seq)

		# Sequence lengths
		self.seq_len = np.array([len(s.seq) for s in self.seq])

		# Alignment numerical array
		l = self.seq_len.max() if type(self.seq_len) == np.ndarray else self.seq_len
		self.array = genepy.alignmentarray(self.seq, length = l)

		# Statistics
		self.statistics = genepy.calcstats(self.seq)














	# Show sequences
	def show(self) :
		genepy.showalignment(self.array)













	# Align sequences 
	def align(self, force = True, iter = False, full = False, full_iter = False, auto = True, threads = False) :

		# System call to ClustalO
		genepy.align(self.filename, force, threads, full, full_iter, iter, auto)

		# Read alignment back in
		self.seq = genepy.readalignment(self.filename.split(".")[0] + "_aligned_genepy.phy")
		
		# Update static members
		self.update()

















	def phylotree(self, nucleotide_frequency = "empirical", bootstrap = -4, search_algorithm = "BEST") :

		if not os.path.isfile(self.filename.split(".")[0] + "_aligned_genepy.phy") :
			print "GenePy can't find an aligned sequence file for %s.\nTry calling .align()." % \
			self.filename.split("/")[-1]

			return


		genepy.phylotree(self.filename, nucleotide_frequency, bootstrap, search_algorithm)
















	def stats(self) :

		# Display statistics
		genepy.stats(self.statistics)



