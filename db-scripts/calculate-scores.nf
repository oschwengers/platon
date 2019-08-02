#!/usr/bin/env nextflow


import java.nio.file.*


// Parameters
protClusterMapping = Paths.get( params.protClusterMapping ).toRealPath()
plasmidsPath       = Paths.get( params.plasmids ).toRealPath()
chromosomesPath    = Paths.get( params.chromosomes ).toRealPath()
nrpcDB             = Paths.get(params.nrpcDB).toRealPath().toString() - '.inf' // workaround as ghostz expects data base name which is not a valid file name


final def proteinInfo = [:]
protClusterMapping.eachLine( {
    if( it[0] != '#') {
        def row = it.split('\t')
        proteinInfo[ row[1] ] = [
            'definition': row[2],
            'length': (row[5] as int)
        ]
    }
} )


Channel.fromPath( plasmidsPath )
    .splitFasta( by: 1, file: true )
    .map( { ['p', it] } )
    .into( { chInputPlasmids; chInputPlasmidCounts } )


Channel.fromPath( chromosomesPath )
    .splitFasta( by: 1, file: true )
    .map( { ['c', it] } )
    .into( { chInputChromosomes; chInputChromosomeCounts } )


chInput = chInputPlasmids.concat( chInputChromosomes )


process callORFs {

    cache false

    input:
    set val(type), file('replicon.fna') from chInput

    output:
    set val(type), file('orf.faa') into chAA

    script:
    if( type == 'p' )
        """
        prodigal -i replicon.fna -a orf.faa -p meta
        """
    else
        """
        prodigal -i replicon.fna -a orf.faa
        """
}


process blastProts {

    cache false
    cpus 2
    memory 2.GB

    input:
    set val(type), file('orf.faa') from chAA

    output:
    set val(type), file('output.tsv') into chBlastResults

    """
    ghostz aln -i orf.faa -o output.tsv -d ${nrpcDB} -b 1 -F F -a ${task.cpus}
    """
}


chProteinPlasmidHits = Channel.create()
chProteinChromosomeHits = Channel.create()


chBlastResults
    .subscribe( onNext: {
        def type = it[0]
        def resultFile = it[1]
        Files.readAllLines( resultFile ).each( {
            def cols = it.split( '\t' )
            def wpId = cols[1]
            if( proteinInfo.containsKey( wpId ) ) {
                // calculate protein length according to prodigal output
                def orfLabel = cols[0].split(' # ')
                int proteinLength = ( Integer.parseInt( orfLabel[2] ) - Integer.parseInt( orfLabel[1] + 1 ) ) / 3
                // check if identity is >= 90 % and alignment length is >= 90 % of the database entry and the called ORF
                int alignmentLength = cols[3] as int
                if( (Float.parseFloat(cols[2]) >= 90)  &&  (alignmentLength >= 0.9*proteinInfo[wpId]['length'])  &&  (alignmentLength >= 0.9*proteinLength) ) {
                    if( type == 'p' )
                        chProteinPlasmidHits << wpId
                    else
                        chProteinChromosomeHits << wpId
                }
            } else {
                println( "no protein info for ${wpId}" )
            }
        } )
    },
    onComplete: {
        chProteinPlasmidHits.close()
        chProteinChromosomeHits.close()
    } )


chResults = Channel.create()


process mergeBlastResults {

    cache false

    input:
    val valBlastPlasmidCounts    from chProteinPlasmidHits.countBy()
    val valBlastChromosomeCounts from chProteinChromosomeHits.countBy()
    val noPlasmids    from chInputPlasmidCounts.count()
    val noChromosomes from chInputChromosomeCounts.count()

    exec:
    def blastPlasmidCounts = valBlastPlasmidCounts // 'variable already defined' bug workaround
    def blastChromosomeCounts = valBlastChromosomeCounts // 'variable already defined' bug workaround

    def wpIds = [].toSet()
    wpIds.addAll( blastPlasmidCounts.keySet() )
    wpIds.addAll( blastChromosomeCounts.keySet() )
    println( "# plasmid proteins: ${blastPlasmidCounts.size()}" )
    println( "# chromosome proteins: ${blastChromosomeCounts.size()}" )
    println( "# distinct proteins: ${wpIds.size()}" )
    wpIds.each( { wpId ->
        int plasmidHits = blastPlasmidCounts[ wpId ] ?: 0
        int chromosomeHits = blastChromosomeCounts[ wpId ] ?: 0
        double pHitFrequency = (double)plasmidHits / noPlasmids
        double cHitFrequency = (double)chromosomeHits / noChromosomes
        double hitFrequencyRatio = ( (pHitFrequency / ( pHitFrequency + cHitFrequency )) - 0.5 ) * 2
        double absDiff = Math.abs( pHitFrequency - cHitFrequency )
        double rds = hitFrequencyRatio * absDiff * 1000
        def prot = proteinInfo[ wpId ]
        chResults << "${wpId}\t${prot.definition}\t${prot.length}\t${plasmidHits}\t${chromosomeHits}\t${pHitFrequency}\t${cHitFrequency}\t${hitFrequencyRatio}\t${absDiff}\t${rds}"
    } )
    chResults.close()
}


chResults.collectFile( sort: false, name: 'rds.full.tsv', storeDir: '.', newLine: true )
