import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
sns.set_context('poster')


if __name__=='__main__':
	parser=argparse.ArgumentParser(description='Produce plots comparing truth vs discovery in simulated motif data.')
	parser.add_argument('--mlf', metavar = 'mlf', type=str, help='Truth or Discovered Motif Length File.',required=True)
	parser.add_argument('--dmif', metavar = 'dmif', type=str, help='Discovered Motif Indices File.',required=True)
	parser.add_argument('--tmif', metavar = 'tmif', type=str, help='Truth Motif Indices File.',required=True)
	parser.add_argument('--dgf', metavar = 'dgf', type=str, help='Discovered Motif Gap File.',required=True)
	parser.add_argument('--tgf', metavar = 'tgf', type=str, help='Truth Motif Gap File.',required=True)
	args=parser.parse_args()

	# 	-	-	-	-	-	-	-	-	-	-	-	-	-	
	# LOAD DATA 

	# Define k1 and k2 length.
	with open(args.mlf,'r') as f:
		f.readline()
		k1,k2=[int(i) for i in f.readline().strip().split('\t')]
	
	# Load Discovery Data.
	FNdiscMotifs=args.dmif
	discMotifsDF=pd.read_csv(FNdiscMotifs,sep='\t')
	discMotifsDF.columns=['k1_start_discovered','k2_start_discovered']
	print('dmif discovered motif indices\n',discMotifsDF)

	FNdiscGaps=args.dgf
	discGapsDF=pd.read_csv(FNdiscGaps,sep='\t')
	discGapsDF.columns=['gap_discovered']
	print('\n\ndgf discovered motif gaps\n',discGapsDF)

	# Load Test Data.
	FNtruthMotifs=args.tmif
	truthMotifsDF=pd.read_csv(FNtruthMotifs,sep='\t')
	truthMotifsDF.columns=['k1_start_truth','k2_start_truth']
	print('\n\ntgf true motif indices\n',truthMotifsDF)

	FNgaps=args.tgf
	gapsDF=pd.read_csv(FNgaps,sep='\t')
	gapsDF.columns=['gap_truth']
	print('\n\ntgf true motif gaps\n',gapsDF)

	# 	-	-	-	-	-	-	-	-	-	-	-	-	-	
	# Figure 1A: Comparing true index of motif implanted with index of discovered motifs.
	
	fig1AFileName='fig1A.png'

	# # Calculate difference between discovered motif and truth motif
	f1aDF=pd.concat([discMotifsDF,truthMotifsDF],axis='columns')
	f1aDF['k1_delta']=f1aDF['k1_start_discovered']-f1aDF['k1_start_truth']
	f1aDF['k2_delta']=f1aDF['k2_start_discovered']-f1aDF['k2_start_truth']
	print('\n\nfig 1A DF\n',f1aDF)

	xlim_bool=True
	xlim=[-50,50]

	fig,axs=plt.subplots(1,2,figsize=(24,12))
	sns.distplot(f1aDF.k1_delta,ax=axs[0],kde=False,bins=20)
	sns.distplot(f1aDF.k2_delta,ax=axs[1],kde=False,bins=20)

	if xlim_bool:
	    axs[0].set_xlim(xlim)
	    axs[1].set_xlim(xlim)

	fig.suptitle('Distribution of  Nucleotide Distance Between Discovered Motif and Truth Motif',size=35)

	axs[0].set_xlabel('\nK1 Discovered - True [nt]')
	axs[1].set_xlabel('\nK2 Discovered - True [nt]')

	axs[0].set_ylabel('Count\n')

	plt.savefig(fig1AFileName,dpi=300)

	# 	-	-	-	-	-	-	-	-	-	-	-	-	-	
	# Figure 1B: Plotting distribution of length of discovered_motif overlapped with truth_motif
	fig1BFileName='fig1B.png'

	# A A A A A
	# . . A A A A A
	# => distance = 2
	# => overlap = k1-distance = 5-2 = 3

	# Convert distance between discovered_motif and truth_motif 
	# from distance to overlap.
	# eg overlap = 0,1,2,...,len(motif)
	f1bDF=f1aDF.copy().loc[:,['k1_delta','k2_delta']]

	f1bDF['k1_overlap']=[abs(abs(dist)-k1) if abs(dist)<k1 else 0 for dist in f1bDF.k1_delta]
	f1bDF['k2_overlap']=[abs(abs(dist)-k2) if abs(dist)<k2 else 0 for dist in f1bDF.k2_delta]

	print('\n\nfig 1B DF\n',f1bDF)

	fig,axs=plt.subplots(1,2,figsize=(24,12))

	sns.distplot(f1bDF.k1_overlap,ax=axs[0],kde=False)
	sns.distplot(f1bDF.k2_overlap,ax=axs[1],kde=False)

	xlim_k1=[0,k1+1]
	xlim_k2=[0,k2+1]

	axs[0].set_xlim(xlim_k1)
	axs[1].set_xlim(xlim_k2)

	fig.suptitle('Distribution of Nucleotide Overlaps Between Discovered and Truth Motifs',size=35)

	axs[0].set_xlabel('\nK1 Discovered/True Motif Overlap')
	axs[1].set_xlabel('\nK2 Discovered/True Motif Overlap')

	axs[0].set_ylabel('Count\n')

	plt.savefig(fig1BFileName,dpi=300)

















