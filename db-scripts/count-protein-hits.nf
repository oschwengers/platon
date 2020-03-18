#!/usr/bin/env nextflow


import java.nio.file.*


// Parameters
plasmidsPath    = Paths.get( params.plasmids ).toRealPath()
chromosomesPath = Paths.get( params.chromosomes ).toRealPath()
pcDb            = Paths.get( params.pcDb ).toRealPath()


Channel.fromPath( plasmidsPath )
    .splitFasta( by: 1, file: true )
    .map( { ['p', it] } )
    .set( { chInputPlasmids } )


Channel.fromPath( chromosomesPath )
    .splitFasta( by: 1, file: true )
    .map( { ['c', it] } )
    .set( { chInputChromosomes } )


chInput = chInputPlasmids.concat( chInputChromosomes )


process callORFs {

    errorStrategy 'retry'
    maxRetries 3

    input:
    set val(type), file('replicon.fna') from chInput

    output:
    set val(type), file('cdss.faa') into chAA

    script:
    if( type == 'p' )
        """
        prodigal -i replicon.fna -a cdss.faa -p meta
        """
    else
        """
        prodigal -i replicon.fna -a cdss.faa
        """
}


process searchProts {

    errorStrategy 'retry'
    maxRetries 3
    cpus 2
    memory 10.GB

    input:
    set val(type), file(cdss) from chAA

    output:
    set val(type), file('output.tsv') into chSearchResults

    when:
    cdss.size() > 0

    """
    diamond blastp --query ${cdss} --db ${pcDb} --threads ${task.cpus} --out output.tsv --max-target-seqs 1 --id 90 --query-cover 80 --subject-cover 80 --tmpdir /var/scratch/
    """
}


chProteinPlasmidHits = Channel.create()
chProteinChromosomeHits = Channel.create()


chSearchResults
    .subscribe( onNext: {
        def type = it[0]
        def resultFile = it[1]
        Files.readAllLines( resultFile ).each( {
            def cols = it.split( '\t' )
            def protId = cols[1]
            if( type == 'p' ) {
                chProteinPlasmidHits << protId
            } else {
                chProteinChromosomeHits << protId
            }
        } )
    },
    onComplete: {
        chProteinPlasmidHits.close()
        chProteinChromosomeHits.close()
    } )


chResults = Channel.create()


process mergeSearchResults {

    cache false

    input:
    val valPlasmidHitCounts    from chProteinPlasmidHits.countBy()
    val valChromosomeHitCounts from chProteinChromosomeHits.countBy()

    exec:
    def plasmidHitCounts = valPlasmidHitCounts // 'variable already defined' bug workaround
    def chromosomeHitCounts = valChromosomeHitCounts // 'variable already defined' bug workaround

    def protIds = [].toSet()
    protIds.addAll( plasmidHitCounts.keySet() )
    protIds.addAll( chromosomeHitCounts.keySet() )
    println( "# plasmid protein hits: ${plasmidHitCounts.size()}" )
    println( "# chromosome protein hits: ${chromosomeHitCounts.size()}" )
    println( "# distinct proteins: ${protIds.size()}" )
    protIds.each( { protId ->
        int plasmidHits = plasmidHitCounts[ protId ] ?: 0
        int chromosomeHits = chromosomeHitCounts[ protId ] ?: 0
        // clusterId \t plasmidHits \t chromosomeHits
        chResults << "${protId}\t${plasmidHits}\t${chromosomeHits}"
    } )
    chResults.close()
}


chResults.collectFile( sort: false, name: 'protein-hit-counts.tsv', storeDir: '.', newLine: true )
