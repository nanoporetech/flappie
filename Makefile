#  Copyright 2018 Oxford Nanopore Technologies, Ltd

#  This Source Code Form is subject to the terms of the Oxford Nanopore
#  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
#  was not  distributed with this file, You can obtain one at
#  http://nanoporetech.com

buildDir ?= build
hdf5Root ?= ''
openblasRedHat ?= ''
openblasRoot ?= ''
releaseType ?= Release

.PHONY: all
all: flappie runnie

flappie: ${buildDir}/flappie
	cp ${buildDir}/flappie flappie

runnie: ${buildDir}/runnie
	cp ${buildDir}/runnie runnie

${buildDir}:
	mkdir ${buildDir}

.PHONY: test
test: ${buildDir}/flappie ${buildDir}/runnie ${buildDir}/flappie_unittest
	cd ${buildDir} && \
	make test

.PHONY: clean
clean:
	rm -rf ${buildDir} flappie runnie

${buildDir}/%: ${buildDir}
	cd ${buildDir} && \
	cmake .. -DCMAKE_BUILD_TYPE=${releaseType} \
	         -DHDF5_ROOT=${hdf5Root} \
	         -DOPENBLAS_REDHAT=${openblasRedHat} \
	         -DOPENBLAS_ROOT=${openblasRoot} && \
	make $*
