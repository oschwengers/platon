
import java.nio.file.*

// Parameters
geneClass = params.geneClass
nrpcPath  = Paths.get( params.nrpc ).toRealPath()
nrpcMappingPath = Paths.get( params.nrpcMapping ).toRealPath()
nrpPath         = Paths.get( params.nrp ).toRealPath()


Channel.fromPath( nrpcPath )
    .splitCsv( sep: '\t')
    .map( { it[0] } )
    .set( { chGeneCluster } )


process extractNRP {
    tag { "hmm - ${cluster}" }

    cache false
    errorStrategy 'ignore'

    input:
    val(cluster) from chGeneCluster

    output:
    file("${geneClass}_${cluster}.hmm") into chHMMs

    script:
    """
    grep ${cluster} ${nrpcMappingPath} | cut -f2 > cluster-proteins.txt
    seqtk subseq ${nrpPath} cluster-proteins.txt > cluster-proteins.faa
    muscle -diags -quiet -in cluster-proteins.faa -out msa.fasta
    hmmbuild -n ${geneClass}_${cluster} --amino --cpu ${task.cpus} ${geneClass}_${cluster}.hmm msa.fasta
    """
}


chHMMs
    .collectFile( sort: false, name: "${geneClass}.hmm", storeDir: '.' )