
import java.nio.file.*


final int MIN_CHROMOSOME_LENGTH = 100_000


Channel.fromPath( params.assembly )
    .splitCsv( skip: 2, sep: '\t'  )
    .filter( { (it[11].toLowerCase() == 'complete genome') } )
    .map( { return [ it[0], it[7], it[19] ] } )
    .set { validGenomes }


process download {
    tag { "${acc} - ${orgName}" }

    cache false
    maxForks 4
    errorStrategy 'retry'
    maxRetries 3

    input:
    set val(acc), val(orgName), val(path) from validGenomes

    output:
    set val(acc), val(orgName), file("${acc}.fna.gz") into chDeflateFasta

    script:
    """
    wget -O ${acc}.fna.gz ${path}/${path.split('/').last()}_genomic.fna.gz
    """
}


process deflate {
    tag { "${acc} - ${orgName}" }

    cache false
    errorStrategy 'retry'
    maxRetries 3

    input:
    set val(acc), val(orgName), file('genome.fna.gz') from chDeflateFasta

    output:
    file("${acc}.fna") into chFasta

    script:
    """
    gzip -cd genome.fna.gz > ${acc}.fna
    """
}


chFasta.splitFasta( by: 1, record: [id: true, desc: true, seqString: true ] )
    .filter( { rec -> (rec.seqString.length() > MIN_CHROMOSOME_LENGTH)  &&  !rec.desc.toLowerCase().contains( 'plasmid' ) } )
    .map( { rec -> ">${rec.id} ${rec.desc}\n${rec.seqString}\n" } )
    .collectFile( sort: false, name: 'refseq-chromosomes.fna', storeDir: '.', tempDir: './nf-tmp' )
