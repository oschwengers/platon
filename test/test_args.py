import os
import pytest

from pathlib import Path
from subprocess import run


@pytest.mark.parametrize(
    "cmd_line",
    [
        (["bin/platon", '--db', 'test/db']),  # no parameter
        (["bin/platon", '--db', 'test/db', '']),  # empty argument
        (["bin/platon", '--db', 'test/db', 'foo.fasta'])  # argument not existing
    ]
)
def test_genome_failing(cmd_line):
    # test genome arguments
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    "cmd_line",
    [
        (["bin/platon", 'test/data/draft-w-plasmids.fna']),  # no parameter
        (["bin/platon", '--db', 'test/data/draft-w-plasmids.fna']),  # missing argument
        (["bin/platon", '--db', '', 'test/data/draft-w-plasmids.fna']),  # empty argument
        (["bin/platon", '--db', 'test/foo', 'test/data/draft-w-plasmids.fna']),  # argument not existing
    ]
)
def test_database_failing(cmd_line):
    # test database arguments

    proc = run(cmd_line)
    assert proc.returncode != 0