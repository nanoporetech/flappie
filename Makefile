#  Copyright 2018 Oxford Nanopore Technologies, Ltd

#  This Source Code Form is subject to the terms of the Oxford Nanopore
#  Technologies, Ltd. Public License, v. 1.0. If a copy of the License 
#  was not  distributed with this file, You can obtain one at
#  http://nanoporetech.com

buildDir ?= build
hdf5Root ?= ''
openblasRoot ?=''
releaseType ?= Release

.PHONY: all
all: flappie
flappie: ${buildDir}/flappie
	cp ${buildDir}/flappie flappie

${buildDir}:
	mkdir ${buildDir}

.PHONY: test
test: ${buildDir}/flappie
	cd ${buildDir} && \
	make test

.PHONY: clean
clean:
	rm -rf ${buildDir} flappie

${buildDir}/flappie: ${buildDir}
	cd ${buildDir} && \
	cmake .. -DCMAKE_BUILD_TYPE=${releaseType} -DHDF5_ROOT=${hdf5Root} -DOPENBLAS_ROOT=${openblasRoot} && \
	make flappie
    
