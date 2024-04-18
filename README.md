# Team Project Template

This repository serves as the default template for analysis projects within the team.

## Main Features

Buildable Package: This template can be built into a package and easily imported into other projects.

Modules Included: Comes with several built-in modules such as datatracker, teamconnector, and findhere to streamline your project development.

Use this template to jump-start your analysis projects and ensure you have a consistent structure across the team.

## Installation

To get started with this template, follow these steps:

Create a new repository on GitHub by clicking "New" on the dashboard.

Choose a pre-defined template from the list under "Repository template."

Name your repository and set its visibility.

Click "Create repository" to initialize the repo with the chosen template.

Clone the repository to your local machine using the command `git clone <your-new-repo-url>`.

Make changes to the code as needed.

**Push** your changes to the repository using the command `git push`.

To customize this template for your project, you'll need to replace some placeholders:

Search for `{{ cookiecutter.full_name }}` and replace it with the appropriate full name associated with this repository.

Look for `{{ cookiecutter.project_slug }}` and substitute it with the repository name or slug you're using.

Replace `{{ cookiecutter.project_short_description }}` with a concise description of this repository's purpose or functionality.

## Usage

### Update pip

Update your pip package manager with the following command:

```bash
make update-pip
```

### Create Environment

To create a new environment:

```bash
make create-env
```

After that, enter the created Conda environment and only run commands in that environment.

```bash
make act
```

### Set Environment Variables

Make sure your environment variables are set appropriately according to `teamconnector` in `~/.env`.
To set up your environment variables:

```bash
make set-env
tcinit ~
```

### Install Versioneer

To install Versioneer:

```bash
make install-version
```

### Install Current Package

For an interactive installation of the current package:

```bash
make install
```

### Update packages

Running `make install` will not install the most up-to-date packages from PyPI. To update all packages in `requirements.txt` to the newest version, you can use the following command:

```bash
make install-reqs
```


### Manage Dependencies

To install or update the required packages:

```bash
make install-reqs
make install-dev
```

## Version Information

### Get Version

To get the version of this library:

```bash
make version
```

### Get Specific Library Version

To get the version of a specific Python library (e.g., hail):

```bash
make lib-version NAME=hail
make hail-version
```

## Install Java 11 using Conda Forge

Hail requires Java 11 for running on local machines, whether it's a SLURM cluster or a local Mac. You can use the Conda Forge package manager.

1. **Open Terminal**

2. **Add Conda-Forge Channel**
   ```bash
   conda config --add channels conda-forge
   ```

3. **Install Java 11**

   ```bash
   conda install -c conda-forge openjdk=11
   ```

4. **Verify Installation**

   ```bash
   java -version
   ```

## Cite

## Maintainer

TJ Singh lab @ singhlab@nygenome.org

## Acknowledgements

## Release Notes