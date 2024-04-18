#!/usr/bin/env python3

from setuptools import setup, find_packages
import versioneer

with open('README.md') as f:
    readme = f.read()

dependencies = []
with open('requirements.txt', 'r') as f:
    for line in f:
        dependencies.append(line.strip())

setup(
    name='{{ cookiecutter.project_slug }}',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='{{ cookiecutter.project_short_description }}',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Singh Lab',
    author_email='singhlab@nygenome.org',
    url='https://github.com/tjsinghlab/{{ cookiecutter.project_slug }}',
    license='MIT',
    python_requires='>=3.7',
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=dependencies
)
