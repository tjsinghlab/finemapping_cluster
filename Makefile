.PHONY: clean install data requirements

# init

SHELL=/bin/bash
PROJECT_NAME={{ cookiecutter.project_slug }}
CONDA_BASE:=$(shell conda info --base)

# install

update-pip:
	pip install -U pip setuptools wheel

install: clean update-pip
	pip install -e .

install-setup: clean update-pip
	python setup.py install

install-reqs:
	pip install --upgrade --no-deps --force-reinstall -r requirements.txt

install-dev:
	pip install --upgrade --no-deps --force-reinstall -r requirements-dev.txt

freeze-reqs: install-requirements
	pip freeze > requirements.txt

version:
	@python -c 'import {{ cookiecutter.project_slug }}; print({{ cookiecutter.project_slug }}.__version__)'

lib-version:
	@python -c 'import $(NAME); print($(NAME).__version__)'

hail-version:
	@python -c 'import hail; print(hail.version())'

# clean

clean: clean-py

clean-py: clean-build clean-pyc clean-test

clean-pyc:
	find . -name '*.py[co]' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

# env

act:
	@(. $(CONDA_BASE)/etc/profile.d/conda.sh && conda activate {{ cookiecutter.project_slug }})

create-env:
	conda create -y --name $(PROJECT_NAME) python=3.11.4 jupyterlab=3.6.3 pylint pep8 versioneer>=0.29

set-env:
	conda env config vars set CLOUD_ROOT=gs://$(PROJECT_NAME) PROJECT_ROOT=`pwd`

remove-env: clean
	( source $(shell conda info --base)/etc/profile.d/conda.sh && \
		conda deactivate && \
		conda remove --name $(PROJECT_NAME) --all )

install-version:
	versioneer install