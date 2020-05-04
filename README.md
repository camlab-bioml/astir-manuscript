
# Summer 2020 automated mass cytometry project

## Setup

### 1. Clone this repo

Navigate to where you would like to work, e.g. your home directory:

```bash
cd ~
```

then

```bash
git clone https://github.com/camlab-bioml/imc-2020.git
```

### 2. Branch this repo

cd into the repo
```bash
cd imc-2020
```

```bash
git checkout -b [name of your branch]
```

This is the branch you'll work on and merge back into master. To push your branch to github, call

```bash
git push -u origin [name of your branch]
```

If you want to checkout master again then 

```bash
git checkout master
```

The usual git rules apply: do some work, add any new or modified files with `git add`, then commit your work with `git commit -m [commit message here]`. Having clear commit messages is important.

### 3. Setup conda environment

You can create the necessary conda environment via

```bash
conda create -f envs/imc.yml
```

which will install all necessary python packages. This can be then loaded via

```bash
conda activate imc
```

**Note this needs performed at every startup** (or modify your `.bashrc` file).

If you update your conda environment you can update it for all users via

```bash
conda env export > envs/imc.yml
```

and git commit.

### 4. Sync raw data

The raw data the pipeline starts at is kept in a `data-raw/` directory. This is populated with symbolic links to the raw files that are stored in `/home/home/ltri/campbell/share/datasets/`. These should **never** be modified without speaking to Kieran first. To create symbolic links to these and setup your `data-raw/` file, run

```bash
./setup-galen.sh
```

### 5. Check you can run snakemake

You can do a dry run of snakemake via

```bash
snakemake -n
```



