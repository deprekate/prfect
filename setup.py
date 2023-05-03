#!/usr/bin/env python3
# encoding: utf-8

import os
#from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils import sysconfig


#os.environ["CC"] = "g++"

sources = ""

module = Extension('prfect',
			language='c++',
			extra_compile_args=['-g'],
			extra_link_args=['-lm'],
			include_dirs=[
						 '.',
						 '...',
						 os.path.join(os.getcwd(), 'include'),
			],
			sources = ["./src/"+item for item in sources.split()]
			)

def readme():
	with open("README.md", "r") as fh:
		long_desc = fh.read()
	return long_desc

def get_version():
	with open("VERSION.md", 'r') as f:
		v = f.readline().strip()
		return v

def main():
	setup (
		name = 'prfect',
		version = get_version(),
		author = "Katelyn McNair",
		author_email = "deprekate@gmail.com",
		description = 'A tool to predict programmed ribosomal frameshifts',
		long_description = readme(),
		long_description_content_type="text/markdown",
		url =  "https://github.com/deprekate/prfect",
		scripts=['prfect.py'],
		classifiers=[
			"Programming Language :: Python :: 3",
			"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
			"Operating System :: OS Independent",
		],
		python_requires='>3.5.2',
		packages=find_packages(),
		package_data = {'': ['*.pkl'],},
		install_requires=['scikit-learn>=0.24.0','sklearn','numpy','pandas','genbank','LinearFold','hotknots','score_rbs', 'packaging'],
		#ext_modules = [module],
	)


if __name__ == "__main__":
	main()
