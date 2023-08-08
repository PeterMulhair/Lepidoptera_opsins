#Label gene tree
~/bin/hyphy/hyphy ~/bin/hyphy-analyses/LabelTrees/label-tree.bf --tree lep_LWS_nucl_trim_aln.nwk --list ../gene_labels_diurnal.txt --output lep_LWS_nucl_trim_aln_labelled_diurnal.nwk

#Run BUSTED-PH
~/bin/hyphy/hyphy ~/bin/hyphy-analyses/BUSTED-PH/BUSTED-PH.bf --alignment codon_align/lep_LWS_nucl_trim_aln.fasta --tree gene-trees/lep_LWS_nucl_trim_aln_labelled_diurnal.nwk --branches Foreground --output codon_align/lep_LWS_nucl_trim_aln.fasta.dirunal.BUSTEDPH.json
