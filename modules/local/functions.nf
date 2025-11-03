#!/usr/bin/env nextflow

def getNullFile() {
    return file(projectDir.resolve('assets/NO_FILE'), checkIfExists: true)
}

def isNullFile(fileObject) {
    return fileObject.getSimpleName() == 'NO_FILE'
}

def parseOptionalPath (path) {
    if (path == null) return getNullFile()
    return file(path, checkIfExists: true)
}

def getSampleLibraryPaths (sampleList) {
    return sampleList.map { sample -> sample.libraries.collect { library -> library.path } }
}

def getSampleNames (sampleList) {
    return sampleList.map { sample -> sample.name }
}

def getLibraryTypes (libraryList) {
    boolean gex = false
    boolean vdj_b = false
    boolean vdj_t = false
    boolean feat = false

    libraryList.collect { library -> library.type }.unique().each() { type ->
        if ('Gene Expression'.equals(type)) {
            gex = true
        } else if ('VDJ-B'.equals(type)) {
            vdj_b = true
        } else if ('VDJ-T'.equals(type)) {
            vdj_t = true
        } else if ('Antibody Capture'.equals(type)) {
            feat = true
        } else {
            error "Invalid feature type: ${type}"
        }
    }

    return [
        'gex': gex,
        'vdj_b': vdj_b,
        'vdj_t': vdj_t,
        'feat':feat
    ]
}
