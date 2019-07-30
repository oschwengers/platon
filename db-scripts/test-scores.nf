#!/usr/bin/env nextflow

import java.lang.String
import java.util.Random
import java.nio.file.*


// Parameters
proteinScoresPath = Paths.get( params.protein_scores ).toRealPath()
plasmidsPath      = Paths.get( params.plasmids ).toRealPath()
chromosomesPath   = Paths.get( params.chromosomes ).toRealPath()
nrpcDB            = Paths.get(params.nrpcDB).toRealPath().toString() - '.inf' // workaround as ghostz expects data base name which is not a valid file name


// Constants
final int PROTEIN_PENALTY = -1.187401
final int RANDOM_CONTIG_LENGTH = 10_000


// parse protein marker information
final def PROTEIN_MARKERS = [:]
proteinScoresPath.eachLine( {
    def cols = it.split( '\t' )
    PROTEIN_MARKERS[ cols[0] ] = [
        length: cols[2] as int,
        ratio: cols[3] as float
    ]
} )


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


chRandomContigs = Channel.create()
chInput = chInputPlasmids.concat( chInputChromosomes )
    .flatMap( {
        def type   = it[0]
        def record = it[1]

        def subsequences = []
        int seqLength = record.seqString.length()
        if( seqLength > RANDOM_CONTIG_LENGTH ) {
            Random rand = new Random()
            for( int i=1; i<11; i++ ) {
                int startPos = rand.nextInt( seqLength - RANDOM_CONTIG_LENGTH )
                String subSequence = record.seqString.substring( startPos, startPos + RANDOM_CONTIG_LENGTH )
                String id = "${record.id}_${i}"
                subsequences << [ type, id, subSequence ]
            }
        } else {
            subsequences << [ type, "${record.id}_1", record.seqString ]
        }
        return subsequences
    } )
    .set( { chRandomContigs } )


process analyzeSubSeq {

    cache false
    cpus 1
    memory '2 GB'

    input:
    set val(type), val(id), val(subSequence) from chRandomContigs

    output:
    set val(type), file('orfs.faa'), file('ghostz.out') into chBlast

    script:
    """
    echo '>${id}\n${subSequence}' > seq.fasta
    prodigal -i seq.fasta -a orfs.faa -p meta
    ghostz aln -i orfs.faa -o ghostz.out -d ${nrpcDB} -b 1 -a ${task.cpus}
    """
}


chBlast.map( {
    def type       = it[0]
    def orfsPath   = it[1]
    def ghostzPath = it[2]

    def orfLengths = [:]
    orfsPath.text.split( '>' ).each( {
        def lines = it.split('\n')
        int length = lines.collect( {it.length()} ).sum() - lines[0].length()
        orfLengths[ lines[0].split(' ')[0] ] = length
    } )
    if( orfLengths.isEmpty() ) {
        return [ type, -1 ]
    } else {
        int noHits = 0
        float scoreSum = 0
        ghostzPath.eachLine( {
            noHits++
            def cols = it.split( '\t' )
            def proteinMarker = PROTEIN_MARKERS[ cols[1] ]
            if( proteinMarker != null ) {
                float alignmentLength = cols[3] as float
                int orfLength = orfLengths[ cols[0].split(' ')[0] ]
                if( (Float.parseFloat(cols[2]) >= 90)  &&  (alignmentLength >= 0.9*proteinMarker['length'])  &&  (alignmentLength >= 0.9*orfLength) )
                    scoreSum += proteinMarker['ratio']
                else
                    scoreSum += PROTEIN_PENALTY
            } else {
                scoreSum += PROTEIN_PENALTY
            }
        } )
        int noOrfs = orfLengths.size()
        if( noHits < noOrfs ) { // add a proteinPenalty for each CDS without hit
            scoreSum += ( noOrfs - noHits ) * PROTEIN_PENALTY
        }
        return [ type, (scoreSum / noOrfs) ]
    }
} )
.flatMap( {
    def type  = it[0]
    def score = it[1]

    def testedCutoffs = []
    for( def cutoff=-50.0; cutoff<=10.0; cutoff+=0.1 ) {
        if( score >= cutoff ) {
            if( type == 'p' )
                testedCutoffs << [ cutoff, [1,0,0,0] ] //tp++
            else
                testedCutoffs << [ cutoff, [0,0,1,0] ] //fp++
        } else {
            if( type == 'p' )
                testedCutoffs << [ cutoff, [0,0,0,1] ] //fn++
            else
                testedCutoffs << [ cutoff, [0,1,0,0] ] //tn++
        }
    }
    return testedCutoffs
} )
.groupTuple()
.map( {
    def cutoff = it[0]
    def list = it[1]

    int tp = 0
    int tn = 0
    int fp = 0
    int fn = 0
    list.each( {
        tp += it[0]
        tn += it[1]
        fp += it[2]
        fn += it[3]
    } )
    float sn  = (float)tp / ( tp + fn ) // sensitivity
    float sp  = (float)tn / ( tn + fp ) // specificity
    float acc = (float)( tp + tn ) / ( tp + fp + tn + fn ) // accuracy
    float ppv = (float)tp / ( tp + fp ) // positive predictive value
    float npv = (float)tn / ( tn + fn ) // negative predictive value
    return "${cutoff}\t${sn}\t${sp}\t${acc}\t${ppv}\t${npv}\t${tp}\t${tn}\t${fp}\t${fn}"
} )
.toSortedList( { a, b -> a.split('\t')[0] as float <=> b.split('\t')[0] as float } )
.flatten()
.collectFile( sort: false, name: 'rds-metrics.tsv', storeDir: '.' , newLine: true )
