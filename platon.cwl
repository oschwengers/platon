#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [platon]
id: platon
label: "Identification and characterization of bacterial plasmid contigs from short-read draft assemblies."

doc: |
      The software and documentation can be found here:
      https://github.com/oschwengers/platon

      Necessary database files can be found here:
      https://doi.org/10.5281/zenodo.3349651

hints:
  SoftwareRequirement:
    packages:
      platon:
      version: [ "1.7.0" ]

requirements:
  ResourceRequirement:
    ramMin: 4096
    coresMin: 1

inputs:
  - doc: Draft genome in fasta format
    id: fasta_file
    inputBinding: {position: 0}
    type: File
    format: edam:format_1929
  - doc: Database path (default = <platon_path>/db)
    id: db
    inputBinding: {prefix: --db}
    type: Directory
  - doc: Mode
    id: mode
    inputBinding: {prefix: --mode}
    type: string
    default: 'accuracy'
  - doc: Run in metagenome mode
    id: metagenome
    inputBinding: {prefix: --meta}
    type: boolean
  - doc: Threads
    id: threads
    inputBinding: {prefix: --threads}
    type: int

outputs:
  - doc: Plasmid contig information
    id: plasmid_info
    type: stdout
  - doc: Plasmid contigs
    id: plasmid_contigs
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.plasmid.fasta'}
  - doc: Chromosome contigs
    id: chromosome_contigs
    type: File
    format: edam:format_1929
    outputBinding: {glob: '*.chromosome.fasta'}

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0003-4216-2721
    s:email: mailto:oliver.schwengers@computational.bio.uni-giessen.de
    s:name: Oliver Schwengers

s:citation: https://doi.org/10.1099/mgen.0.000398
s:codeRepository: https://github.com/oschwengers/platon
s:license: https://spdx.org/licenses/GNU GPL3
s:programmingLanguage: Python

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
 - https://schema.org/version/latest/schema.rdf
 - http://edamontology.org/EDAM_1.18.owl
