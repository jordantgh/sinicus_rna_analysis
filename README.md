# Rhinolophus sinicus RNA sequencing analysis pipelines

## Basics

- `reference_pipeline.sh` gets in the RNA seq and WGS data, and does genome indexing followed by alignment of the reads.

- `denovo.sh` currently assumes the RNA seq data are available, and generates *de novo* transcriptomes for all of the specified samples, which is followed by a quantification and tabulation step (main output is in `denovo_results` folder).

- The `tabulate_tx.py` script is used in `denovo.sh` to create a `.csv` with 
transcript sequences matching our target sequence (currently using the Sinicus "IFITM1-like" variant (X2))

- The `R_files` directory has a (very) basic plotting function, this can be improved significantly. 

- I use poetry and pyenv for specifying the local python environment (hence the `pyproject.toml` and `.python-version` files). If this is something you want to try I can provide instructions about poetry and pyenv (they are quite simple).

## Contributing 

### Setup

First you need to install Git. This is the tool to use for working with project repositories locally. You can find it at e.g. https://git-scm.com/download/mac on mac and follow the installation instructions there.

#### Configure Git

Open Terminal and configure your Git username and email using the following commands:

```{bash}
git config --global user.name "Your Name"
git config --global user.email "your.github.email@example.com"
```

#### Set up SSH Authentication: 

Using SSH keys allows you to access Github via Git in a secure way without needing to manually enter a password each time. The SSH key acts like a password stored on your machine that Git knows where to look for. Whenever you generate an SSH key, two files are created - the public key and the private key. They are linked. The public key (file ending in `.pub`) is what you will eventually store on Github. 

To generate an SSH key, running the following command in Terminal:

```{bash}
ssh-keygen -t ed25519 -C "your.github.email@example.com"
```

Follow the prompts and accept the default file location, or provide a custom location if needed (usually it will be something like `~/.ssh/id_rsa`). If you decide to use a passphrase for the key, you may end up adding some other complexities to the setup (it probably isn't necessary, unless you're worried about someone stealing your laptop and messing with your Github).

Create an SSH config file to manage your SSH connections. In Terminal, run the following command to do this:

```{bash}
touch ~/.ssh/config
```

Open the config file with your preferred text editor and add the following configuration lines:

```{bash}
# GitHub
Host github.com
    HostName github.com
    User git
    IdentityFile ~/.ssh/your_private_key # Use the actual name
```

(Note, the identityfile used here should be the one that does *not* end in `.pub` -- that's for Github)

#### Add the SSH key to your GitHub account: 

Go to account settings (click on your profile icon on the top right > `Settings`)

Select `SSH and GPG keys` on the left hand sidebar. 

Click the `New SSH Key` button (top right). Give it a title (e.g. home_laptop, work_laptop, something like that). To populate the `Key` field, you'll want the contents of the public key file.

In your Terminal, enter

```{bash}
cat ~/.ssh/your_public_key.pub # Use the actual name
```

which should print out the contents of the public key. Copy and paste this into the `Key` field back on GitHub, then click `Add SSH key` to complete the process. 

Now run 

```{bash}
ssh -T git@github.com
```

to test the connection. You should see an output like the following:

```{bash}
ssh -T git@github.com
The authenticity of host 'github.com (140.82.121.4)' can't be established.
ED25519 key fingerprint is SHA256:[ ... ].
This key is not known by any other names.
Are you sure you want to continue connecting (yes/no/[fingerprint])? 
```
You should type "yes" and will then see this message:

```{bash}
Warning: Permanently added 'github.com' (ED25519) to the list of known hosts.
```
You should now be able to begin working with this repository.

As an aside, if you are working on mac as I assume, you should follow the instructions here to make sure you don't end up flooding the repo with temporary mac files (`.DS_store` files):

https://www.theodinproject.com/lessons/foundations-setting-up-git

### Start working on the repo

#### Clone the Repository: 

Open Terminal, navigate to the desired folder where you want to keep the project files (I usually use a `repos` or `projects` folder), and run the following command:

```
git clone git@github.com:jordantgh/sinicus_rna_analysis.git
```

#### Create a New Branch: 

Before making any changes, create a new branch for your work. This keeps your contributions organized and separated from the `main` branch. Use the following command to create a new branch and switch to it:

```{bash}
git checkout -b your-new-branch-name
```

Replace your-new-branch-name with a descriptive name for your branch (usually a feature name, e.g., fix-issue-123 or add-new-feature).

#### Make Changes: 

Open the folder in your preferred code/text editor. Edit the files in the repository or make your own new files. 

After making changes, stage the modified files and create a commit with a descriptive message. Run the following commands:

```{bash}
git add .
git commit -m "Your commit message here"
```

Replace the commit message with a brief description of the changes you've made. Note the `.` after `git add` is a shorthand for "the current directory". In other words, if you made changes on several files this would stage all of them together to be committed together. It is often better to individually add files e.g. `git add /path/to/script.py` so that your commit messages are focused and relevant.

#### Sync with Main Branch: 

Before pushing your changes, it's a good idea fetch the latest changes from the main branch and merge them with your branch:

```{bash}
git fetch origin
git merge origin/main
```

#### Push Your Changes:

Push your branch with your changes to the repository on GitHub:

```{bash}
git push origin your-feature-branch-name
```

#### Create a Pull Request: 

Go to the repo on GitHub and click the `Pull Requests` tab along the top (next to `Issues`) then click the `New pull request` button. Select your branch from the dropdown menu, fill in a descriptive title and description for your pull request, and click "Create pull request".

I will then be able to review your changes and provide feedback. If it looks all good, I'll merge your pull request into the main branch.

(BTW, the nomenclature of "Pull request" might seem a bit strange; a `pull` in Github is a shorthand way of doing both `fetch` and `merge` as you did above in one operation, so "pull request" is another way of saying you are "requesting" me as the maintainer to merge your branch into `main`)

Additional Resources
To learn more about Git and GitHub, check out these resources:

**TOP Git Setup Guide** - *From a course aimed at web developers, but universally applicable*:
https://www.theodinproject.com/lessons/foundations-setting-up-git

**Pro Git Book** - *A comprehensive guide to Git*: 
https://git-scm.com/book/en/v2

**GitHub Docs** - *Official GitHub documentation*:
https://docs.github.com/en

**Git Cheatsheet** - *A handy reference for common Git commands*:
https://education.github.com/git-cheat-sheet-education.pdf 