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

    errorStrategy 'ignore'
    maxRetries 3
    conda 'prodigal=2.6.3'

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

    errorStrategy 'ignore'
    maxRetries 3
    cpus 4
    memory '12 GB'
    conda 'diamond=2.0.14'

    input:
    set val(type), file(cdss) from chAA

    output:
    set val(type), file('output.tsv') into chSearchResults

    when:
    cdss.size() > 0

    """
    diamond blastp --query ${cdss} --db ${pcDb} --threads ${task.cpus} --fast --out output.tsv --max-target-seqs 1 --id 90 --query-cover 80 --subject-cover 80
    """
}


chProteinHits = Channel.create()
chResults = Channel.create()

chSearchResults
    .map( {
        def type = it[0]
        def resultFile = it[1]
        def protIds = []
        Files.readAllLines( resultFile ).each( {
            protIds << it.split( '\t' )[1]
        } )
        return [type, protIds]
    } )
    .branch {
        plasmids: it[0] == 'p'
        chromosomes: it[0] == 'c'
    }
    .set { chProteinHits }


chProteinHits.plasmids
    .map( { it[1] } )
    .flatten()
    .countBy()
    .set( { chProteinPlasmidHits } )


chProteinHits.chromosomes
    .map( { it[1] } )
    .flatten()
    .countBy()
    .set( { chProteinChromosomeHits } )


process mergeSearchResults {

    cache false
    executor 'local'

    input:
    val valPlasmidHitCounts    from chProteinPlasmidHits
    val valChromosomeHitCounts from chProteinChromosomeHits

    exec:
    def protIds = [].toSet()
    protIds.addAll( valPlasmidHitCounts.keySet() )
    protIds.addAll( valChromosomeHitCounts.keySet() )
    println( "# plasmid protein hits: ${valPlasmidHitCounts.size()}" )
    println( "# chromosome protein hits: ${valChromosomeHitCounts.size()}" )
    println( "# distinct proteins: ${protIds.size()}" )
    protIds.each( { protId ->
        int plasmidHits = valPlasmidHitCounts[ protId ] ?: 0
        int chromosomeHits = valChromosomeHitCounts[ protId ] ?: 0
        chResults << "${protId}\t${plasmidHits}\t${chromosomeHits}"
    } )
    chResults.close()
}


chResults.collectFile( sort: false, name: 'protein-hit-counts.tsv', storeDir: '.', newLine: true )
