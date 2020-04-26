
# Summer 2020 automated mass cytometry project

## Setup

### 1. Clone this repo

Navigate to your home directory e.g.

```bash
cd ~
```

then

```bash
git clone https://github.com/camlab-bioml/imc-2020.git
```

### 2. Branch this repo

```bash
git checkout -b [name of your branch]
```

This is the branch you'll work on and merge back into master. To push your branch to github, call

```bash
git push -u origin [name of your branch]
```

If you want to checkout master again then 

```git checkout master```

The usual git rules apply: do some work, add any new or modified files with `git add`, then commit your work with `git commit -m [commit message here]`. Having clear commit messages is important.

### 3. Setup conda environment

Note if you wi

### 4. Sync raw data

The raw data the pipeline starts at is kept in a `data-raw/` directory. This is populated with symbolic links to the raw files that are stored in `/home/home/ltri/campbell/share/datasets/`. These should **never** be modified without speaking to Kieran first. To create symbolic links to these and setup your `data-raw/` file, run

```sh
./setup-galen.sh
```



