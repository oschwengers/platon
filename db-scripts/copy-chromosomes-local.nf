
import java.nio.file.*


final int MAX_PLASMID_LENGTH = 1_000_000
final int MIN_CHROMOSOME_LENGTH = 1_000_000


db       = params.db
ncbiPath = '/nfs/biodb/ncbi_genomes'


Channel.fromPath( "${ncbiPath}/${db}/bacteria/assembly_summary.txt" )
    .splitCsv( skip: 2, sep: '\t'  )
    .filter( { (it[11].toLowerCase() == 'complete genome') } )
    .map( { return [ it[0], it[7], it[19] - 'ftp://ftp.ncbi.nlm.nih.gov/genomes/' ] } )
    .set { validGenomes }


process deflate {
    tag { "${acc} - ${orgName}" }

    errorStrategy 'retry'
    maxRetries 3

    input:
    set val(acc), val(orgName), val(path) from validGenomes

    output:
    file("${acc}.fna") into chFasta

    script:
    """
    gzip -cd ${ncbiPath}/${path}/${path.split('/').last()}_genomic.fna.gz > ${acc}.fna
    """
}

chPlasmids      = Channel.create()
chChromosomes   = Channel.create()
chFasta.splitFasta( by: 1, record: [id: true, desc: true, seqString: true ] )
    .choice( chPlasmids, chChromosomes ) { rec -> rec.desc.toLowerCase().contains( 'plasmid' ) ? 0 : 1 }

chPlasmidsSize = Channel.create()
chPlasmidsSeqs = Channel.create()
chPlasmids.filter( { rec -> rec.seqString.length() < MAX_PLASMID_LENGTH } )
    .separate( chPlasmidsSize, chPlasmidsSeqs ) { record -> ["${record.id}\t${record.desc}\t${record.seqString.length()}\n", ">${record.id} ${record.desc}\n${record.seqString}\n" ] }
chPlasmidsSize.collectFile( sort: false, name: "${db}.plasmids.tsv", storeDir: '.' )
chPlasmidsSeqs.collectFile( sort: false, name: "${db}.plasmids.fna", storeDir: '.' )

chChromosomesSize = Channel.create()
chChromosomesSeqs = Channel.create()
chChromosomes.filter( { rec -> rec.seqString.length() > MIN_CHROMOSOME_LENGTH } )
    .separate( chChromosomesSize, chChromosomesSeqs ) { record -> ["${record.id}\t${record.desc}\t${record.seqString.length()}\n", ">${record.id} ${record.desc}\n${record.seqString}\n" ] }
chChromosomesSize.collectFile( sort: false, name: "${db}.chromosomes.tsv", storeDir: '.' )
chChromosomesSeqs.collectFile( sort: false, name: "${db}.chromosomes.fna", storeDir: '.' )
