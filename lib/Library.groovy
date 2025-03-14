import java.nio.file.Path
import groovy.transform.Immutable
import nextflow.io.ValueObject
import nextflow.util.KryoHelper

@ValueObject
@Immutable(copyWith=true, knownImmutables = ['id', 'type', 'path'])
public class Library {
    static { KryoHelper.register(Library) }

    String id
    String type
    String path

    static Library create (libraryId, libraryType, libraryPath) {
        // verify libraryType
        if (libraryType == 'VDJ')
            error ('Feature type "VDJ" not supported. Specify VDJ-T or VDJ-B')
        if (!['Gene Expression', 'VDJ-T', 'VDJ-B', 'Antibody Capture'].contains(libraryType))
            error ("Feature type \"${libraryType}\" not supported.")
        
        // verify libraryPath
        if (!libraryPath)
            error ("Library path cannot be falsy ($libraryPath).")
        def pathObject = Path.of(libraryPath)
        if (!pathObject.exists())
            error ("Library path $libraryPath does not exist.")
        if (!pathObject.isDirectory())
            error ("Library path $libraryPath is not a directory.")
        return new Library (libraryId, libraryType, libraryPath)
    }

    static Library create (libraryMap) {
        return Library.create (
            libraryMap.fastq_id,
            libraryMap.feature_types,
            libraryMap.fastqs
        )
    }
}