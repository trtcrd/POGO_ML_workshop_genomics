# Tutorial Part 2

## Prediction of functional traits using Traitar

Traitar uses __[support vector machines (SVMs)](https://en.wikipedia.org/wiki/Support_vector_machine)__, which are a class of algorithms using linear multidimensional hyperplanes based on supervised training. Though Traitar does not, SVMs can also use a [kernel function](https://en.wikipedia.org/wiki/Kernel_method), simply put data transformation, to be turned into non-linear methods. SVMs were trained using phenotype annotations from the [GIDEON database](https://doi.org/10.1186/1476-072X-4-10) based on the content of predicted [Pfam protein families](https://pfam.xfam.org/) of their genomes.


### 1. Select 2 or 3 MAGs and prepare metadata

Traitar is implemented in Python and has already been installed in the AMI used for the tutorial, along with its dependencies. Since it takes about 10 minutes per genome to run the Traitar prediction workflow, we will select only two or three of the nitrogen-fixing Tara Oceans MAGs for annotation. First, prepare a new directory for the tutorial, then download and expand the MAG sequence files into it:

```mkdir Trait_Pred_Tutorial
cd Trait_Pred_Tutorial
wget https://www.dropbox.com/s/4ufmyg3r5n0gxw6/15_NITROGEN_FIXING_MAGs.tar.gz
tar xvzf 15_NITROGEN_FIXING_MAGs.tar.gz
```

The file ``Metadata-15-nitrogen-fixing-MAGs.txt`` in ``15_NITROGEN_FIXING_MAGs`` has information about the different MAGs. Have a look at it, and select two or three MAGs that you want to analyse using Traitar. Copy these (or simply delete the others) to a new directory, and make a metadata file in tab-separated text format for Traitar [as specified in its GitHub repo (see below also)](https://github.com/hzi-bifo/traitar), with the first column giving the filename (must match the FASTA file) and the second giving a sample name. The headers must be ``sample_file_name`` and ``sample_name``.


sample_file_name|sample_name|category
--- | --- | --- 
sample1_file_name|sample1_name|sample_category1
sample2_file_name|sample2_name|sample_category2
sample3_file_name|sample3_name|sample_category3
...|...|...
sampleN_file_name|sampleN_name|sample_categoryN

The first step of Traitar is to predict protein-coding reading frames of the MAG sequences using [Prodigal](https://github.com/hyattpd/Prodigal). These are then aligned to the profile [Hidden Markov Models (HMMs)](https://en.wikipedia.org/wiki/Hidden_Markov_model), of Pfam, statistical models trained from multiple sequences of protein familieis. This is done using the program [HMMER](hmmer.org). Though HMMER has been pre-installed, we need to tell Traitar to download the Pfam database and to where, before it can be used:

```mkdir ~/Pfam
sudo traitar pfam ~/Pfam
```

### 2. Run Traitar

Now that you have selected some genomes to analyse, start the Traitar workflow with the metadata file you specified. We will run the job on four processes (or more, depending on what instance you spawned):

```cd ~/Trait_Pred_Tutorial
traitar phenotype <in_dir> <sample_file> from_nucleotides TraitarPred -c 4
```

..where ``sample_file`` is the metadata file you prepared and ``in_dir`` is the file containing your MAG sequences. This step will take about 10 minutes per genome, so it is a good time to grab a coffee, or read more about Traitar, MAGs, SVMs or something else.


### 3. Look at Traitar predictions

When Traitar has completed, copy the output directory, that we called ``TraitarPred``, to your laptop using scp so that you can look at the graphics and output files generated:

`scp  -r -i <my_aws_key.pem> ubuntu@<Your_AWS_Public_DNS>:Trait_Pred_Tutorial/TraitarPred .` 

The files in directory ``gene_prediction`` contain the predicted coding sequences by Prodigal, in FASTA and [General Feature Format (.gff)](https://en.wikipedia.org/wiki/General_feature_format) and those in ``annotation`` contain the identified PFAM hits based on the former. These were then used by Traitar's trained SVM model to predict traits using the default "pyhpat" model alone and with an additional model they call "PGL", taking into the evolutionary history of each strain in the training data, to model the loss and gain of new protein families / traits. These results are in ``phenotype_prediction``.

The file ``predictions_majority-vote_combined.txt`` lists the predicted distribution of all modelled traits across the MAGs analysed. The results are also presented as a heatmap in ``heatmap_combined.pdf`` and a detailed list can be found in the .dat files for each MAG in the directory ``feat_gffs``.

Can you find any traits that are present in only some of them? Have a look at the explanation of the treat name in [Table S1 of the paper](https://msystems.asm.org/content/1/6/e00101-16#DC1). You can compare any interesting predictions to [Supplemntary Table 6 in Delmont et al](https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-018-0176-9/MediaObjects/41564_2018_176_MOESM8_ESM.xlsx).

For example, the MAG HBD-06, appears to be capable of denitrification in addition to N fixation, which is also mentioned in the paper, consistent with the Traitar prediction "Nitrite to gas" for this MAG. To take a closer look at what protein families that were used to predict this feature, you can use:

`traitar show 'Nitrite to gas'`

## Prediction using a more general trait predictor (Farrell et al)

Unfortunately, key environmental traits such as Nitrogen fixation are missing from the traits targeted by Traitar, which used a pathogen-focused database (GIDEON) for selecting and traniing models. Luckily, there is a more general program described in [this preprint by Farrell et al](https://www.biorxiv.org/content/10.1101/307157v1) where they used a bigger set of prokaryotic genomes and the more extensive trait database [FAPROTAX](https://www.nature.com/articles/s41559-016-0015) for training. This classifier, called GenePhene is similar to Traitar, except using [LASSO regression](https://en.wikipedia.org/wiki/Lasso_(statistics)) instead of SVMs. Though it is still experimental, the implemented model has been installed on the AMI you are using so that we can test it.

We will run GenePhene using the list of identified PFAMs generated by the workflow script of Traitar (and used as input to Traitar). However, we need to do a small hack because GenePhene has one small difference in the format it expects for PFAMs. Open the file ``TraitarPred/annotation/pfam/summary.dat`` using a command line word editor, for example emacs (or one with a GUI, like gedit, if you want to try to log in with an X11 connection):

`emacs TraitarPred/annotation/pfam/summary.dat`

...then write "genome_ID" before the first tab and then replace all tabs with commas and save as "traitarPFAMs.csv" or similar. 

If you want to make life easier you can use the script we made for this instead, in the genephene directory:

`~/genephene/convertTraitarPFAMs.py TraitarPred/annotation/pfam/summary.dat > ./traitarPFAMs.csv`

Now, we can run GenePhene using the PFAMs identified earlier as input, writing predictions to e.g.  `someNFix_PFAM_GenePhene.csv` (GenePhene, by the way, is implemented in python3 and is installed under the home directory, not made available as an executable in the PATH):

`python3 ~/genephene/genephene_genome_predict.py -i traitarPFAMs.csv -g Pfam -o someNFix_PFAM_GenePhene.csv`

Copy the output file to your laptop and have a look at it. Nitrogen fixation is now a trait. Are all your MAGs predicted to have this trait? What other interesting nitrogen related or other traits can you find? 

Though GenePhene works with PFAM as input, it works better with the [KOs (KEGG Orthology terms)](https://www.genome.jp/kegg/ko.html) of the  database KEGG (Kyoto Encyclopedia of Genes and Genomes), which is very useful for functional annotation. Read through the information on the KEGG website about KEGGs briefly. However, alignment to KEGG would take too long, so we have made a file with ready formatted KOs for most of the 15 nitrogen fixing MAGs from Tara Oceans available. To run the program using this as input, use the same script as before, but without option "-g Pfam" since KOs is the default input:

`python3 ~/genephene/genephene_genome_predict.py -i ~/genephene/precomputed/TaraMAGs_Nfix_KEGG.csv -o mostNFix_KO_GenePhene.csv`

Looking in the output (``mostNFix_KO_GenePhene.csv``), can you find any traits that were predicted differently for the 2-3 MAGs you selected? How about nitrogen fixation?

If you want to learn more about the selected traits you can look them up in the [FAPROTAX database](https://www.nature.com/articles/s41559-016-0015). Othwerise, now is a good time to explore your own dataset if you bring any, or perhaps even to call it a day.

Thank you!

Anders and Tristan




