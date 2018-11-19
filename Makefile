buildDir ?= build
releaseType ?= Release

.PHONY: all flappie
all: flappie
flappie: ${buildDir}/flappie

${buildDir}:
	mkdir ${buildDir}

.PHONY: test
test: ${buildDir}/flappie
	cd ${buildDir} && \
	make test

.PHONY: clean
clean:
	rm -rf ${buildDir}

${buildDir}/flappie: ${buildDir}
	cd ${buildDir} && \
	cmake .. -DCMAKE_BUILD_TYPE=${releaseType} && \
	make
    
