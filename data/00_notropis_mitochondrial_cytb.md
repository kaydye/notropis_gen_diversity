# Mitochondrial analysis and phylogenetic reconstruction using cytochrome *b*


##### Obtain data
```{bash}
unzip raw_seq_kd.zip
```
###### Should have a list of the ab1 files for the forward and reverse reads of each individual

#### Geneious was used to align forward and reverse sequences (ab1 files), and create a consensus fasta file for each individual.
  * Some individuals were removed due to poor sequence quality
    ###### NM-116, NM-230, NM-259, and NM-320 were removed
  * Consensus sequences are in /data/Notorpis_consensus_seqs
  * Outgroup fasta sequences obtained from GenBank, found in /data/Notropis_outgroup_seqs

#### Remove header from outgroup fasta files
```{bash}
sed -i '1s/^>.*/>Agosia chrysogaster/' agosia_chrysogaster.fasta
sed -i '1s/^>.*/>Notropis atherinoides/' notropis_atherinoides.fasta
sed -i '1s/^>.*/>Notropis hudsonius HUD1/' notropis_hud1.fasta
sed -i '1s/^>.*/>Notropis jemezanus JEME2/' notropis_jeme2.fasta
sed -i '1s/^>.*/>Notropis jemezanus/' notropis_jemezanus.fasta
sed -i '1s/^>.*/>Notropis atherinoides NAN1/' notropis_nan1.fasta
sed -i '1s/^>.*/>Notropis percobromus NPB1/' notropis_pb1.fasta
sed -i '1s/^>.*/>Notropis stilbius/' notropis_stilbius.fasta
sed -i '1s/^>.*/>Pimephales vigilax/' pimephales_vigilax.fasta
```

#### Create multi.fasta with all consensus and outgroup sequences

###### Linearizing outgroup fasta files
```{bash}
for file in *.fasta; do
    sed -e 's/\(^>.*$\)/#\1#/' "$file" | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > "${file%.fasta}_linear.fasta"
done
```

###### Renaming of fasta files
```{bash}
cut -d " " -f1 Notropis_cytb_consensus.fasta | sed 's/Nucleotide_alignment_//g' > Notropis_cytb_consensus_rename.fasta
```

###### Create multi-fasta to be aligned
```{bash}
cat Notropis_cytb_consensus_rename.fasta ../outgroups/*linear.fasta > all.fasta
```

#### Align sequences using ClustalOmega
```{bash}
../programs/bin/bin/clustalo -i all.fasta -o test.fasta --outfmt=fasta
```

###### Check alignment manually using UGene
   * Reimport alignment as ugene_alignment.fasta

#### Trim
###### Linearize multi-fasta for trimming
```{bash}
awk '/^>/ {if (seq != "") {print seq; seq="";} print; next} {seq = seq $0} END {if (seq != "") print seq}' ugene_alignment.fasta > linear_alignment.fasta
```

###### Trim first 26 bases and last 29 bases off each sequence
```{bash}
awk '/^>/ {if (seq != "") {print substr(seq,27,length(seq)-55); seq="";} print; next} {seq = seq $0} END {if (seq != "") print substr(seq,27,length(seq)-55)}' linear_alignment.fasta > trimmed_alignment.fasta
```

-![#1589F0] Trimmed alignment contains 1,121 bp and 56 individuals + outgroups `#1589F0`

#### Run model test
```{bash}
modeltest−ng −i alignment_mega.fa 
```      

##### Phylogenetic Reconstruction (BEAST)

###### Beauti  
```{GUI}
Use Beauti to make xml
    Import alignment alignment_mega.fa
    Use log file from modeltest
        
Best model according to BIC
---------------------------
Model:              TrN+I+G4
lnL:                -4841.1534
Frequencies:        0.2590 0.3093 0.1430 0.2888
Subst. Rates:       1.0000 54.6190 1.0000 1.0000 9.8690 1.0000 
Inv. sites prop:    0.6265
Gamma shape:        2.7361
Score:              10623.2517
Weight:             0.7341

    Gamma category count : 6
    Shape : 2.7361 (estimate)
    Proportion invariant : 0.6265 (estimate)
    Subst Model : TN93 (everything estimated)
    
  MCMC 
    
    Store every 1000
    tracelog : notropis_1
    treelog : notropis_1
    
  File > save as > notropis_1.xml
```
###### run five times  

#### BEAST
```{GUI}
open xml file : notropis_1
Run
```

#### Tracer
```{GUI}
#look at log file from each BEAST run 
```

#### Tree Annotator
```{GUI}
Maximum Clade Credibility Tree
Node heights : mean heights
```    

