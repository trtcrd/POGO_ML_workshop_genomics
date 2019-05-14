# Tutorial Part 2

## Prediction of functional traits using Traitar

Traitar uses __[support vector machines (SVMs)](https://en.wikipedia.org/wiki/Support_vector_machine)__, which are a class of algorithms using linear multidimensional hyperplanes based on supervised training. Though Traitar does not, SVMs can also use a [kernel function](https://en.wikipedia.org/wiki/Kernel_method), simply put data transformation, to be turned into non-linear methods. SVMs were trained using phenotype annotations from the [GIDEON database](https://doi.org/10.1186/1476-072X-4-10) based on the content of predicted [Pfam protein families](https://pfam.xfam.org/) of their genomes.


### 1. Select 2 or 3 MAGs and prepare metadata

Traitar is implemented in Python and has already been installed in the AMI used for the tutorial, along with its dependencies. Since it takes about 10 minutes per genome to run the Traitar prediction workflow, we will select only two or three of the nitrogen-fixing Tara Oceans MAGs for annotation. First, prepare a new directory for the tutorial, then download and expand the MAG sequence files into it:

`mkdir Trait_Pred_Tutorial`
`cd Trait_Pred_Tutorial`
`wget https://www.dropbox.com/s/4ufmyg3r5n0gxw6/15_NITROGEN_FIXING_MAGs.tar.gz`
`tar xvzf 15_NITROGEN_FIXING_MAGs.tar.gz`

The file Metadata-15-nitrogen-fixing-MAGs.txt in ``15_NITROGEN_FIXING_MAGs`` has information about the different MAGs. Have a look at it, and select two or three MAGs that you want to analyse using Traitar. Copy these to a new directory or simply delete the others, and make a metadata file in tab-separated text format for Traitar [as specified in its GitHub repo](https://github.com/hzi-bifo/traitar), with the first column giving the filename (must match the FASTA file) and the second giving a sample name. The headers must be ``"sample_file_name" and "sample_name"``.

The first step of Traitar is to predict protein-coding reading frames of the MAG sequences using [Prodigal](https://github.com/hyattpd/Prodigal). These are then aligned to the profile [Hidden Markov Models (HMMs)](https://en.wikipedia.org/wiki/Hidden_Markov_model), of Pfam, statistical models trained from multiple sequences of protein familieis. This is done using the program [HMMER](hmmer.org). Though HMMER has been pre-installed, we need to tell Traitar to download the Pfam database and to where, before it can be used:

`mkdir ~/Pfam`
`sudo traitar pfam ~/Pfam`

### 2. Run Traitar

Now that you have selected some genomes to analyse, start the Traitar workflow with the metadata file you specified. We will run the job on four processes (or more, depending on what instance you spawned):

`cd ~/Trait_Prediction_Tutorial`
`traitar phenotype <in_dir> <sample_file> from_nucleotides <out_dir> -c 4`

..where sample_file is the metadata file you prepared and in_dir is the file containing your MAG sequences. This step will take about 10 minutes per genome, so it is a good time to grab a coffee, or read more about Traitar, MAGs, SVMs or something else. When the workflow is completed, the output will be in out_dir.


### 3. Look at Traitar predictions


