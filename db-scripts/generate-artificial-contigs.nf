#!/usr/bin/env nextflow

import java.lang.String
import java.util.Random
import java.nio.file.*


// Parameters
plasmidsPath    = Paths.get( params.plasmids ).toRealPath()
chromosomesPath = Paths.get( params.chromosomes ).toRealPath()


// Constants
final List RANDOM_CONTIG_LENGTHS = [1_000, 5_000, 10_000, 20_000, 50_000]


chInputPlasmids = Channel.create()
Channel.fromPath( plasmidsPath )
    .splitFasta( by: 1, record: [id: true, seqString: true ] )
    .map( { ['p', it] } )
    .set( {chInputPlasmids} )


chInputChromosomes = Channel.create()
Channel.fromPath( chromosomesPath )
    .splitFasta( by: 1, record: [id: true, seqString: true ] )
    .map( { ['c', it] } )
    .set( {chInputChromosomes} )


chInput = chInputPlasmids.concat( chInputChromosomes )
    .flatMap( {
        def type   = it[0]
        def record = it[1]

        def subsequences = []
        int seqLength = record.seqString.length()
        for(randomContigLength in RANDOM_CONTIG_LENGTHS) {
            if( seqLength > randomContigLength ) {
                Random rand = new Random()
                for( int i=1; i<11; i++ ) {
                    int startPos = rand.nextInt( seqLength - randomContigLength )
                    String subSequence = record.seqString.substring( startPos, startPos + randomContigLength )
                    subsequences << ">${type}_${record.id}_${randomContigLength}_${i}\n${subSequence}"
                }
            } else {
                subsequences << ">${type}_${record.id}_${randomContigLength}_1\n${record.seqString}"
            }
        }
        return subsequences
    } )
    .collectFile( sort: false, name: 'artificial-contigs.fna', storeDir: '.' , newLine: true, cache: false, tempDir: './tmp' )
