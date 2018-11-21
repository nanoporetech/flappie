buildDir ?= build
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
	cmake .. -DCMAKE_BUILD_TYPE=${releaseType} && \
	make
    
