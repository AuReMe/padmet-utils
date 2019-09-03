#!/usr/bin/env python
# Copyright (c) 2019, Dyliss team Inria <gem-aureme@inria.fr>
#
# This file is part of padmet_utils.
#
# padmet_utils is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# padmet_utils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with padmet_utils.  If not, see <http://www.gnu.org/licenses/>.
# -*- coding: utf-8 -*-import setuptools

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="padmet_utils",
    version="0.0.8",
    author="Dyliss Inria",
    author_email="gem-aureme@inria.fr",
    description="Modules to manage, explore and connect metabolic networks using Padmet and SBML formats",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AuReMe/padmet-utils",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)