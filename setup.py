#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
import sys

# with open('README.rst') as readme_file:
#     readme = readme_file.read()

# with open('HISTORY.rst') as history_file:
#     history = history_file.read()

with open('requirements.txt') as req_file:
    reqs = req_file.read()

requirements = reqs.splitlines()
# print(requirements, file=sys.stdout)

# setup_requirements = ['pytest-runner', ]

# test_requirements = ['pytest', ]

setup(
    author="M. Williams",
    author_email='mwilliams@protoxrd.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
    ],
    description="PROTO Rietveld Refinement",
    entry_points={
        'console_scripts': [
            'Rietveld_CCTBX=rietveld_server:main',
        ],
    },
    install_requires=requirements,
    # long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords=['PROTO', 'Rietveld Refinement'],
    name='Rietveld_CCTBX',
    # packages=find_packages(include=['./twisted_server']),
    # setup_requires=setup_requirements,
    # test_suite='tests',
    # tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
