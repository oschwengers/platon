
import argparse
import gzip
import xml.etree.ElementTree as et
from pathlib import Path


parser = argparse.ArgumentParser(
    description='Filter Uniprot\'s UniRef90 XML files to bacterial subsequences and init pc db.'
)
parser.add_argument('--taxonomy', action='store', help='Path to NCBI taxonomy node.dmp file.')
parser.add_argument('--xml', action='store', help='Path to UniRef xml file.')
parser.add_argument('--fasta', action='store', help='Path to MPS fasta file.')
parser.add_argument('--tsv', action='store', help='Path to MPS tsv file.')
args = parser.parse_args()

taxonomy_path = Path(args.taxonomy).resolve()
xml_path = Path(args.xml).resolve()
fasta_path = Path(args.fasta)
tsv_path = Path(args.tsv)


def is_taxon_child(child, LCA, taxonomy):
    parent = taxonomy.get(child, None)
    while(parent is not None and parent != '1'):
        if(parent == LCA):
            return True
        else:
            parent = taxonomy.get(parent, None)
    return False


taxonomy = {}
with taxonomy_path.open() as fh:
    for line in fh:
        cols = line.split('\t|\t', maxsplit=2)
        taxonomy[cols[0]] = cols[1]


ns = 'http://uniprot.org/uniref'
with gzip.open(str(xml_path), mode='rt') as fh_xml, fasta_path.open(mode='wt') as fh_fasta, tsv_path.open(mode='wt') as fh_tsv:
    i = 0
    for event, elem in et.iterparse(fh_xml):
        # print(elem)
        if(elem.tag == "{%s}entry" % ns):
            # print(elem)
            if('Fragment' not in elem.find("./{%s}name" % ns).text):  # skip protein fragments
                tax_property = elem.find("./{%s}property[@type='common taxon ID']" % ns)
                if(tax_property is not None):
                    tax_id = tax_property.attrib['value']
                    if(is_taxon_child(tax_id, '2', taxonomy)):
                        # print("tax-id=%s" % tax_id)
                        uniref90_id = elem.attrib['id'][9:]  # remove 'UniRef90_' prefix
                        rep_member = elem.find("./{%s}representativeMember/{%s}dbReference" % (ns, ns))
                        try:
                            prot_name = rep_member.find("./{%s}property[@type='protein name']" % ns).attrib['value']
                        except:
                            prot_name = ''
                        seq = elem.find("./{%s}representativeMember/{%s}sequence" % (ns, ns)).text.upper()
                        fh_fasta.write(">%s\n%s\n" % (uniref90_id, seq))
                        fh_tsv.write("%s\t%s\t%i\n" % (uniref90_id, prot_name, len(seq)))
                        # print(seq)
            i += 1
            if((i % 1000000) == 0):
                print("%i..." % i)
            elem.clear()  # forstall out of memory errors
