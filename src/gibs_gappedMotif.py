import argparse
import numpy as np
from numpy.random import choice

def profileWithPseudoCounts(Dna):
    '''returns positional weight matrix from inputted list of dna'''
    As,Cs,Gs,Ts=[],[],[],[]
    
    for j in range(len(Dna[0])): # for each nt position
        frequencies={'A':1,'C':1,'G':1,'T':1}
        for i in range(len(Dna)): # for each dna sequence
            nt= Dna[i][j]
            frequencies[nt]+=1
        As+=[frequencies['A']/sum(frequencies.values())]
        Cs+=[frequencies['C']/sum(frequencies.values())]
        Gs+=[frequencies['G']/sum(frequencies.values())]
        Ts+=[frequencies['T']/sum(frequencies.values())]
    
    return [As,Cs,Gs,Ts]

# # Test profileWithPseudoCounts(Dna)
# print('profile','AAAA','AAAT','AAGA','ACAA =')
# for i in profileWithPseudoCounts(['AAAA','AAAT','AAGA','ACAA']): print(i)

def scoreMotifs(profile,motifs):
    # get consensus sequence
    consensus=[]
    for pA,pC,pG,pT in zip(profile[0],profile[1],profile[2],profile[3]):
        pMax=max([pA,pT,pG,pC])
        if pA==pMax:
            consensus+=['A']
        elif pC==pMax:
            consensus+=['C']
        elif pG==pMax:
            consensus+=['G']
        elif pT==pMax:
            consensus+=['T']
    
    # get score
    score=0
    for m in motifs:
        for ntIdx in range(len(motifs[0])):
            if m[ntIdx] != consensus[ntIdx]:
                score+=1
    
    return score

# # Test scoreMotifs()
# print(scoreMotifs(profileWithPseudoCounts(['AAA','AAA','AAA','AAA']),['AAT','AAA','TAA']))

def kmers(text,k):
    '''returns list of kmers, a sliding window of length k across a string'''
    kmers=[]
    for i in range(len(text)-k+1):
        kmers+=[text[i:i+k]]
    return kmers

def kmerProbabilities(text,k,profile):
    '''returns probability (1 per index from 0 to |text|-k)''' 
    kmerProbs=[]
    nt2idx={'A':0,'C':1,'G':2,'T':3}
    
    for kmer in kmers(text,k):
        # calculate probability of each nt in kmer
        ntProbs=[]
        for i in range(k):
            nt=kmer[i]
            ntProbs+=[profile[nt2idx[nt]][i]]

        # calculate total probability (product of nt probs)
        totalProb=1
        for prob in ntProbs:
            totalProb*=prob

        # add kmer to dict
        kmerProbs+=[totalProb]
    
    # return first kmer with max kmer probability
    return kmerProbs

def profileBiasedGappedK1K2mer(profile,dna,k1,k2,G):
	'''
	Calculates probability each k1mer and k2mer in dna given profile.
	Finds maximum probability of k1-g-k2mer at each position by taking max probaility at all possible g's separating k1mer and k2mer.
	
	Inputs:
		profile - motif profile
		dna - string
		k1 - integer length of motif 1
		k2 - integer length of motif 2
		G - collection of gaps g to consider

	Returns:
		k1-g-k2-mer
	'''
	
	# Calculate profile of k1 and k2.
	profile_k1=[p[:k1] for p in profile]
	profile_k2=[p[k1:] for p in profile]
	
	# Calculate kmer probability at each nt of dna.
	probs_k1 = kmerProbabilities(dna,k1,profile_k1)
	probs_k2 = kmerProbabilities(dna,k2,profile_k2)

	# Calculate the maximum k1-g-k2mer probability at each position out of all G.
	maxProb_ofAllG     =[] # max probability of k1-g-k2mer corresponding to index of dna
	gWithMaxProb_ofAllG=[] # g corresponding to max probability of k1-g-k2mer corresponding to index of dna
	for pos in range(len(dna)-k1-k2-min(G)+1):
		# Calculate all summed probabilities of k1/k2 for each value of g IF the value of g does not cause indexing error.
		maxProb_atPos = max([ (probs_k1[pos]+probs_k2[pos+k1+g] , g) for g in G if pos+k1+g+k2 <= len(dna) ])
		
		maxProb_ofAllG.append(maxProb_atPos[0])
		gWithMaxProb_ofAllG.append(maxProb_atPos[1])
	
	# Coerce probabilities to sum to one.
	sum_maxProb_ofAllG = sum(maxProb_ofAllG)
	maxProb_ofAllG=[i/sum_maxProb_ofAllG for i in maxProb_ofAllG]

	# print(maxProb_ofAllG)
	# print(gWithMaxProb_ofAllG)

	# Return random index of string based on probabilities in maxProb_ofAllG
	dice = choice(range(len(dna)-k1-k2-min(G)+1),p=maxProb_ofAllG)
	g = gWithMaxProb_ofAllG[dice]

	# print('dna',dna)
	# print('dice',dice)
	# print('g',g)

	# print('returned motif',dna[dice:dice+k1]+dna[dice+k1+g:dice+k1+g+k2])

	return dice, dna[dice:dice+k1]+dna[dice+k1+g:dice+k1+g+k2] , g

# # Test profileBiasedGappedK1K2mer
# dna='AGCTCTTTGTTATAG'
# k1=3
# k2=3
# G=[0,1]
# profile=profileWithPseudoCounts(['TGTTAT','TGTTAT','TGTTAT','TGTTAT'])

# print('Profile of TGTTAT')
# for i in profile: print(i)

# profileBiasedGappedK1K2mer(profile,dna,k1,k2,G)

def initializeRandomKmers(Dna,k1,k2,g_min,g_max):
	motifs         = []
	motifs_indices = []
	gaps           = []
	for d in Dna:
		randStart=np.random.choice(range(0,len(d)-(k1+k2+g_max)))
		randGap  =np.random.choice(range(g_min,g_max+1))

		k1Start=randStart
		k1End=k1+k1Start

		k2Start=randStart+k1+randGap
		k2End=randStart+k1+randGap+k2

		motifs.append(d[k1Start:k1End]+d[k2Start:k2End])
		motifs_indices.append(randStart)
		gaps.append(randGap)

	# print('motifs')
	# for i in motifs: print(i)

	# print('\nmotif_indicers')
	# for i in motifs_indices: print(i)

	# print('\ngaps')
	# for i in gaps: print(i)

	return motifs, motifs_indices, gaps

def gappedMotifGibbs(Dna,G,k1,k2,runs,iters_perRun):
	'''
	Gibbs-sampler for gapped composite motif finding from t strings dna.
	
	Inputs:
		Dna - a collection of dna strings of length t.
		G - a collection of integers corresponding to gap length.
		k1 - integer length of k1 motif.
		k2 - integer length of k2 motif.
		N - an integer specifying number of iterations of sampling to perform. 

	Returns (in order):
		bestMotifs - t x k1+k2 matrix. Contains unspaced motif of length k1+k2. Rows correspond to each dna in Dna.
		bestMotifs_indices - t x 2 matrix. col[0] k1 start in dna, col[1] k2 start in dna. Rows correspond to each dna in Dna. 
		bestGaps - t x 1 matrix. Cointains a gap in each row corresponding to gap observed in gap-free motif of bestMotifs matrix.
		gapDistribution - 2 x G matrix. row[0] contains gap lengths tested. row[1] contains frequency of each gap in G.
	'''

	t = len(Dna)

	# # Choose arbitrary motifs to instantiate motifs.
	# # Choose first k1-g-k2mer where k1 starts at 0 and gap = min(G). 
	# # Append unspaced motif to motifs.
	# motifs         = [] x
	# motifs_indices = [] x
	# gaps           = [] x
	# g_min          = min(G)
	# for d in Dna:
	# 	motifs.append(d[0:k1]+d[k1+g_min:k1+g_min+k2])
	# 	motifs_indices.append(0)
	# 	gaps.append(g_min)

	# Initialize with random motifs
	g_min=min(G)
	g_max=max(G)
	motifs, motifs_indices, gaps = initializeRandomKmers(Dna,k1,k2,g_min,g_max)

	# Initialize best motifs.
	bestMotifs = motifs
	bestMotif_indices = motifs_indices
	bestGaps   = gaps
	score_bestMotifs = scoreMotifs(profileWithPseudoCounts(bestMotifs),bestMotifs)

	for run in range(runs):
		motifs, motifs_indices, gaps = initializeRandomKmers(Dna,k1,k2,g_min,g_max)

		# Sample n times.
		for j in range(iters_perRun):
			
			# Choose random d in Dna to sample.
			i = choice(range(t))
			
			# print('\n-------------\ni',i)
			# print('\nmotifs')
			# for z in motifs: print('',z)
			# print('\ndna[i]',Dna[i])

			# Create profile of all motifs except ith motif. 
			motifs_exceptI = [m for m_index,m in enumerate(motifs) if m_index!=i ]
			profile_exceptI = profileWithPseudoCounts(motifs_exceptI)

			# print('\nmotifs except i')
			# for z in motifs_exceptI: print('',z)


			# Roll a biased die to randomly choose k1/g/k2 motif.
			motifs_indices[i], motifs[i], gaps[i] = profileBiasedGappedK1K2mer(profile_exceptI, Dna[i], k1, k2, G)

			# print('\n>> New Motifs')
			# for z in motifs: print('',z)

			# Check if new motifs are better than old motifs.
			if scoreMotifs(profileWithPseudoCounts(motifs),motifs) < score_bestMotifs:
				bestMotifs        = motifs
				bestMotif_indices = motifs_indices
				bestGaps          = gaps
				score_bestMotifs  = scoreMotifs(profileWithPseudoCounts(bestMotifs),motifs)

	# Calculate frequency of each gap length observed
	gapDist=[]
	for g in G:
		gapDist.append(bestGaps.count(g))


	return bestMotifs, bestMotif_indices, bestGaps, [G,gapDist]

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Find gapped composite motifs in a collection of dna strings.')
	parser.add_argument('--k1', metavar = 'k1', type=int, help='Integer length of k1.',required=True)
	parser.add_argument('--k2', metavar = 'k2', type=int, help='Integer length of k2.',required=True)
	parser.add_argument('--gaps', metavar = 'gaps', type=int, nargs='+', help='Valid gaps to search for. Collection of integers (min length 1) delimited by spaces.',required=True)
	parser.add_argument('--runs', metavar = 'runs', type=int, help='Number of times to randomly initialize motifs.',required=True)
	parser.add_argument('--iters', metavar = 'iters', type=int, help='Number of Gibs sampling iterations per run.',required=True)
	parser.add_argument('--input', metavar = 'input', type=str, help='A .txt file. Inputted collection Dna of string dna. Each string dna should be on its own line.',required=True)
	parser.add_argument('--outputdir', metavar = 'outputdir', type=str, help='Output /dir/ (do not append with .txt or other suffix).',required=True)
	args=parser.parse_args()

	# Read input File.
	Dna=[line.strip() for line in open(args.input,'r').readlines() if line != '\n']

	# Determine directory and basename of outputted files.
	if '.' in args.input:
		outBasename=args.input[:args.input.index('.')] 
	outBasename=args.outputdir+'/'+outBasename

	# Define filenames to write to.
	logFN=outBasename+'_log.txt'
	motifLengthsFN=outBasename+'_motif_lengths.txt'
	bestMotifsFN=outBasename+'_motif_strings.txt'
	bestMotif_indicesFN=outBasename+'_motif_indices.txt'
	bestGapsFN=outBasename+'_gap_counts.txt'
	gapDistributionFN=outBasename+'_gap_distribution.txt'

	with open(logFN,'w') as f:
		f.write('# Inputs and Parameters'+'\n')
		f.write('Dna input: '+args.input+'\n')
		f.write('gaps: '+str(args.gaps)+'\n')
		f.write('k1 length: '+str(args.k1)+'\n')
		f.write('k2 length: '+str(args.k2)+'\n')
		f.write('Number of initial random motifs (runs): '+str(args.runs)+'\n')
		f.write('Number of iterations of single run (iters): '+str(args.iters)+'\n')
	
	bestMotifs, bestMotif_indices, bestGaps, gapDistribution=gappedMotifGibbs(Dna,args.gaps,args.k1,args.k2,args.runs,args.iters)

	# Write Motif Lengths
	with open(motifLengthsFN,'w') as f:
		f.write('k1_length\tk2_length\n')
		f.write(str(args.k1)+'\t'+str(args.k2)+'\n')

	# Write bestMotifs
	with open(bestMotifsFN,'w') as f: 
		f.write('k1_motif\tk2_motif\n')
		f.write('\n'.join([unspaced[:args.k1]+'\t'+unspaced[args.k1:] for unspaced in bestMotifs]))
	# print('\nbestMotifs',bestMotifs)

	# Write bestMotif Indices
	with open(bestMotif_indicesFN,'w') as f:
		f.write('k1_start'+'\t'+'k2_start'+'\n')
		for k1gk2_index , g in zip(bestMotif_indices,bestGaps):
			f.write(str(k1gk2_index)+'\t'+str(k1gk2_index+args.k1+g)+'\n')
	# print('\nbestMotif_indices',bestMotif_indices)

	# Write bestGaps
	with open(bestGapsFN,'w') as f: 
		f.write('gap'+'\n')
		f.write('\n'.join([str(i) for i in bestGaps]))
	# print('\nbestGaps',bestGaps)

	with open(gapDistributionFN,'w') as f: 
		f.write('gap\tcounts\n')
		for g,g_freq in zip(gapDistribution[0],gapDistribution[1]):
			f.write(str(g)+'\t'+str(g_freq)+'\n')
	# print('\ngapDistribution',gapDistribution)


























