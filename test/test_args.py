import os
import pytest

from subprocess import run


@pytest.mark.parametrize(
    'parameters',
    [
        (['--db', 'test/db']),  # no parameter
        (['--db', 'test/db', '']),  # empty argument
        (['--db', 'test/db', 'foo.fasta'])  # argument not existing
    ]
)
def test_genome_failing(parameters, tmpdir):
    # test genome arguments
    cmd_line = ['bin/platon', '--output', tmpdir] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        (['test/data/draft-w-plasmids.fna']),  # no parameter
        (['--db', 'test/data/draft-w-plasmids.fna']),  # missing argument
        (['--db', '', 'test/data/draft-w-plasmids.fna']),  # empty argument
        (['--db', 'test/foo', 'test/data/draft-w-plasmids.fna']),  # argument not existing
    ]
)
def test_database_failing(parameters, tmpdir):
    # test database arguments
    cmd_line = ['bin/platon', '--output', tmpdir] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0