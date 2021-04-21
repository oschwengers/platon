import pytest

from pathlib import Path
from subprocess import run

from .conftest import FILES


@pytest.mark.slow
def test_platon_w_plasmids(tmpdir):
    # full test on draft assembly containing plasmid contigs
    proc = run(['bin/platon', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', 'test/data/draft-w-plasmids.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
        assert output_path.stat().st_size > 0


@pytest.mark.slow
def test_platon_wo_plasmids(tmpdir):
    # full test on draft assembly containing no plasmid contigs
    proc = run(['bin/platon', '--db', 'test/db', '--output', tmpdir, '--prefix', 'test', 'test/data/draft-wo-plasmids.fna'])
    assert proc.returncode == 0

    tmpdir_path = Path(tmpdir)
    for file in FILES:
        output_path = tmpdir_path.joinpath(file)
        assert Path.exists(output_path)
    
    chromosome_path = tmpdir_path.joinpath('test.chromosome.fasta')
    assert chromosome_path.stat().st_size > 0

    plasmids_path = tmpdir_path.joinpath('test.plasmid.fasta')  # test if plasmid fasta file is empty
    assert plasmids_path.stat().st_size == 0