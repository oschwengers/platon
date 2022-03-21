#!/usr/bin/env nextflow

import java.lang.String
import java.util.Random
import java.nio.file.*


// Parameters
mpsPath     = Paths.get( params.mps ).toRealPath()
contigsPath = Paths.get( params.contigs ).toRealPath()
pcDb        = Paths.get( params.pcDb ).toRealPath()


// Constants
final int PROTEIN_PENALTY = 0.0


// parse protein marker information
final def PROTEIN_MARKERS = [:]
mpsPath.eachLine( {
    def cols = it.split( '\t' )
    PROTEIN_MARKERS[ cols[0] ] = [
        length: cols[2] as int,
        rds: cols[3] as float
    ]
} )


chInputPlasmids = Channel.create()
Channel.fromPath( contigsPath )
    .splitFasta( by: 1, record: [id: true, seqString: true ] )
    .map( { [it.id.split('_')[0], it.id, it.seqString] } )
    .set( {chContigs} )


process searchProts {

    errorStrategy 'ignore'
    maxRetries 3
    cpus 4
    memory '4 GB'
    conda 'prodigal=2.6.3 diamond=2.0.14'

    input:
    set val(type), val(id), val(subSequence) from chContigs

    output:
    set val(type), file('cdss.faa'), file('output.tsv') into chBlast

    script:
    """
    echo '>${id}\n${subSequence}' > seq.fna
    prodigal -i seq.fna -a cdss.faa -p meta
    if [ -s cdss.faa ]
    then
        diamond blastp --query cdss.faa --db ${pcDb} --fast --threads ${task.cpus} --out output.tsv --max-target-seqs 1 --id 90 --query-cover 80 --subject-cover 80
    else
        touch output.tsv
    fi
    """
}


chBlast.map( {
    def type       = it[0]
    def orfsPath   = it[1]
    def resultPath = it[2]

    int noOrfs = orfsPath.text.findAll( { it == '>'} ).size()
    if( noOrfs == 0 ) {
        return [ type, PROTEIN_PENALTY ]
    } else {
        int noHits = 0
        float scoreSum = 0
        resultPath.eachLine( {
            noHits++
            def cols = it.split( '\t' )
            def proteinMarker = PROTEIN_MARKERS[ cols[1] ]
            if( proteinMarker != null ) {
                scoreSum += proteinMarker['rds']
            } else {
                scoreSum += PROTEIN_PENALTY
            }
        } )
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
