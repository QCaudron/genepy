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


		




		# Number of sequences
		self.len = len(self.seq)



		# Sequence lengths
		self.seq_len = np.array([len(s.seq) for s in self.seq])









	# Show sequences
	def show(self) :
		self.array = genepy.alignmentarray(self.seq, length = self.seq_len.max())
		genepy.showalignment(self.array)







	# Align sequences 
	def align(self, force = True, iter = False, full = False, full_iter = False, auto = True, threads = False) :

		# Generate ClustalO command
		command = genepy.align(self.filename, force, threads, full, full_iter, iter, auto)
		print "Calling :    %s" % command

		# System call
		os.system(command)












