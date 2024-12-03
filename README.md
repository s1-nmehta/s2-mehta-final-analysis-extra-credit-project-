# s2-mehta-final-analysis-extra-credit-project-


# Lab 3: Finding homologs with BLAST KEY

Step 1: make a new directory called GFM1
```
mkdir ~/lab03-$MYGIT/GFM1
```
Step 2: Enter the new directory using the cd command and ensure you are in the right directory
```
cd GFM1
pwd
```
Step 3: Downlaod the GFM1 protein seuqence in fasta format from NCBI
```
ncbi-acc-download -F fasta -m protein "NP_079272.4"
```
Step 4: Verify the download was successful and view the file
```
ls
less NP_079272.4.fa
```
Step 5: Perform a BLAST search of GFM1 against the other animal proteomes. The command specifies the output folder and it name.
```
blastp -db ../allprotein.fas -query NP_079272.4.fa -outfmt 0 -max_hsps 1 -out GFM1.blastp.typical.out
```
Step 6: View the output file containing BLAST results
```
less GFM1.blastp.typical.out
```
Step 7: Perform a BLAST search in tabular format and save the output
```
blastp -db ../allprotein.fas -query NP_079272.4.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out GFM1.blastp.detail.out

```
Step 8: Open the detailed output file with BLAST results in tabular format
```
less -S GFM1.blastp.detail.out
```
Step 9: Filter the BLAST hits with e-value < 1e-30 and save the filtered hits
```
awk '{if ($6 < 1e-30) print $1}' GFM1.blastp.detail.out > GFM1.blastp.filtered.out
wc -l GFM1.blastp.filtered.out
```
> I had 29 homologs with an e-value less than 1e-30

Step 11: Count how many paralogs were found per species
```
grep -o -E "^[A-Z]\.[a-z]+" GFM1.blastp.filtered.out | sort | uniq -c
```
> Homo sapiens had 9 paralogs, while Mus musculus had 8. Other species had fewer hits.

# Lab 4: Gene family sequence alignment

Step 1: Create a new directory called GFM1 for lab 4
```
mkdir ~/lab04-$MYGIT/GFM1

```
Step 2: Navigate into the new directory and confirm you’re in the correct location
```
cd ~/lab04-$MYGIT/GFM1
pwd

```
Step 3: Extract and isolate the homologous sequences of interest based on the BLAST results using a specific filtering criterion
```
seqkit grep --pattern-file ~/lab03-$MYGIT/GFM1/GFM1.blastp.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/GFM1/GFM1.homologs.fas
     
```
> this loaded 29 patterns from the file

Step 4: Perform a multiple sequence alignment (MSA) using the extracted homolog sequences and save the results in the specified output file
```
muscle -align ~/lab04-$MYGIT/GFM1/GFM1.homologs.fas -output ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas

```

Step 5: Format the alignment output into a scrollable view for easier examination
```
alv -kli ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas | less -RS
```
Modify the alignment settings to highlight columns where the most frequent amino acid appears in more than 50% of sequences
```
alv -kli --majority ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas | less -RS
```
Generate a graphical representation of the multiple sequence alignment using an R script
```
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas
```
Step 8: Calculate the total width of the alignment (number of positions in the alignment)
```
alignbuddy -al ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas    
```
> length of alignment: 712

Step 9: Determine the number of columns removed after filtering out gaps from the alignment
```
alignbuddy -trm all ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas | alignbuddy -al       
```
> After removing gaps, the alignment was reduced to 231 positions, meaning 481 gaps were present.

Step 10: Measure the alignment length after removing invariant columns (positions that don’t vary across sequences)
```
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas | alignbuddy -al
```
> After filtering invariant positions, the alignment was reduced to 683 positions, leaving 29 invariant columns.

Step 11: Use T-Coffee to compute the percentage identity of the alignment excluding gaps
```
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas -output sim 
```
> percent indentity using t_coffee was 42.3

Step 12: Use AlignBuddy to calculate the percentage identity, ensuring gaps are excluded
```
alignbuddy -pi ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas | awk ' (NR>2) { for (i=2;i<=NF ;i++){ sum+=$i;num++ } } END{ print(100*sum/num) } '

```
> The percent identity calculated using AlignBuddy was 38.6%

# Lab 5: Gene Family Phylogeny using IQ-TREE

Step 1: Create a new directory for GFM1 in Lab 5 and navigate into it
```
mkdir ~/lab05-$MYGIT/GFM1
cd ~/lab05-$MYGIT/GFM1
```
Step 2: Remove any sequences with duplicate labels from the alignment file and copy the filtered alignment file from Lab 4 to this directory
```
sed 's/ /_/g' ~/lab04-$MYGIT/GFM1/GFM1.homologs.al.fas | seqkit grep -v -r -p "dupelabel" > ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.fas
```
Step 3: Generate a maximum likelihood phylogenetic tree for GFM1 using IQ-TREE. The tool will determine the best-fit amino acid substitution model, conduct a tree search, estimate branch lengths, and calculate bootstrap support values.
```
iqtree -s ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.fas -bb 1000 -nt 2
```
Step 4: View the generated tree file in Newick format using the nw_display command
```
nw_display ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.fas.treefile
```
Step 5: Use R to visualize the tree graphically and export it as a PDF file.
```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.fas.treefile ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.fas.treefile.pdf 0.4 15
```
Step 6: Reroot the tree using the midpoint rooting method to interpret it more meaningfully
```
gotree reroot midpoint -i ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile
```
Step 7: View the midpoint-rooted tree first on the command line and then save it as an SVG file for further visualization
```
nw_order -c n ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile | nw_display -
nw_order -c n ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.svg

```
Step 8: Convert the SVG file of the midpoint-rooted tree into a PDF format for documentation
```
convert ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.svg ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.pdf

```
Step 9: Change the visualization from a phylogram (branch lengths proportional to substitutions) to a cladogram (focuses on topology), and save the output in both SVG and PDF formats
```
nw_order -c n ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.midCl.treefile.svg
convert ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.midCl.treefile.pdf

```
Step 10: Reroot the tree using an outgroup to provide evolutionary context. Replace the placeholder outgroup names with specific taxa from your dataset if applicable.
```
nw_reroot ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.fas.treefile H.sapiens_GFM1_homolog1 H.sapiens_GFM1_homolog2 > ~/lab05-$MYGIT/GFM1/GFM1.homologsf.outgroup.treefile
```
Step 11: Visualize the outgroup-rooted tree and export it as an SVG file, then convert it to a PDF.
```
nw_order -c n ~/lab05-$MYGIT/GFM1/GFM1.homologsf.outgroup.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/GFM1/GFM1.homologsf.outgroup.treefile.svg
convert ~/lab05-$MYGIT/GFM1/GFM1.homologsf.outgroup.treefile.svg ~/lab05-$MYGIT/GFM1/GFM1.homologsf.outgroup.treefile.pdf

```
# Lab 6: Reconciling a Gene and Species Tree

Step 1: Create a new directory for Lab 6 and navigate into it.
```
mkdir ~/lab06-$MYGIT/GFM1
cd ~/lab06-$MYGIT/GFM1

```
Step 2: Copy the midpoint-rooted tree file for GFM1 from Lab 5 into the new directory for reconciliation
```
cp ~/lab05-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile

```

Step 3: Use NOTUNG to reconcile the gene tree with the species tree. This step identifies events like gene duplications, losses, and transfers. Save the reconciled tree output as a PNG file.
```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/GFM1/

```
Step 4: View the reconciliation events from the NOTUNG output file.
```
less ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.rec.events.txt

```
Step 5: Convert the reconciled NOTUNG file into RecPhyloXML format using a Python script for compatibility with further visualization tools.
```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.rec.ntg --include.species
```
Step 6: Use the thirdkind tool to generate an SVG file that visually integrates the reconciled gene tree and species tree.
```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.rec.ntg.xml -o ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.rec.svg

```
Step 7: Convert the SVG reconciliation tree file into a PDF format for easy access and presentation.
```
convert -density 150 ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/GFM1/GFM1.homologsf.al.mid.treefile.rec.pdf

```

# Lab 8: Protein Domain Prediction (GFM1)

Step 1: Step 1: Create a directory for Lab 8 and move into it.
```
mkdir ~/lab08-$MYGIT/GFM1
cd ~/lab08-$MYGIT/GFM1

```
Step 2: Remove stop codons from the GFM1 sequence file and save the output into the new directory.
```
sed 's/*//' ~/lab04-$MYGIT/GFM1/GFM1.homologs.fas > ~/lab08-$MYGIT/GFM1/GFM1.homologs.fas

```
Step 3: Run rpsblast to compare GFM1 protein sequences against the Pfam database to identify potential protein domains. The output file contains information like domain length, start and end positions, and e-values to assess domain significance.
```
rpsblast -query ~/lab08-$MYGIT/GFM1/GFM1.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001

```
Step 4: Copy the midpoint-rooted tree from Lab 5 into the Lab 8 directory for integrating domain annotations.
```
cp ~/lab05-$MYGIT/GFM1/GFM1.homologsf.outgroup.treefile ~/lab08-$MYGIT/GFM1

```
Step 5: Use R to integrate the midpoint-rooted tree with predicted protein domains from the rpsblast results. Save the annotated tree with domains as a PDF.
```
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/GFM1/GFM1.homologsf.outgroup.treefile ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out ~/lab08-$MYGIT/GFM1/GFM1.tree.rps.pdf

```
Step 6: Use mlr to view the rpsblast domain annotations as a tabular format in a spreadsheet or terminal for further analysis.
```
mlr --inidx --ifs "\t" --opprint cat ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out | tail -n +2 | less -S

```
Step 7: Analyze the rpsblast output to identify key domain-related statistics, such as the most commonly found domains, proteins with multiple domains, and those with the longest annotated domains.
```
cut -f 1 ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out | sort | uniq -c
cut -f 6 ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out | sort | uniq -c
awk '{a=$4-$3; print $1, "\t", a;}' ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out | sort -k2nr

```
Step 8: Continue by identifying proteins with the shortest domains and the best e-values for domain matches.
```
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out | sort -t $'\t' -k2,2g
cut -f 1,6 ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out | tr ',' '\n' | sort | uniq -c | awk '$1 > 1'
cut -f 1,6 ~/lab08-$MYGIT/GFM1/GFM1.rps-blast.out | awk '{print length($2), $0}' | sort -nr | head -n 1

```


