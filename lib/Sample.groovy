import groovy.transform.Immutable
import nextflow.io.ValueObject
import nextflow.util.KryoHelper

@ValueObject
@Immutable(copyWith=true, knownImmutables = ['name', 'libraries'])
public class Sample {
    static { KryoHelper.register(Sample) }
    String name
    List libraries

    static Sample create (sampleMap) {
        return new Sample (
            sampleMap.name,
            sampleMap.libraries.collect { libraryMap -> Library.create(libraryMap) }
        )
    }

    def getFeatureTypes() {
        return libraries.collect { library -> library.type }
    }

    def hasGeneExpressionLibrary() {
        return getFeatureTypes().contains('Gene Expression')
    }

    def hasVdjBLibrary() {
        return getFeatureTypes().contains('VDJ-B')
    }

    def hasVdjTLibrary() {
        return getFeatureTypes().contains('VDJ-T')
    }

    def hasFeatureLibrary() {
        return getFeatureTypes().contains('Antibody Capture')
    }
}