import os
import pytest

from subprocess import run


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # no parameter
        (['']),  # empty argument
        (['foo.fasta'])  # argument not existing
    ]
)
def test_genome_failing(parameters, tmpdir):
    # test genome arguments
    cmd_line = ['bin/platon', '--db', 'test/db', '--output', tmpdir] + parameters
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'parameters',
    [
        ([]),  # not provided
        (['--db']),  # missing path
        (['--db', '', ]),  # empty
        (['--db', 'test/foo']),  # not existing
    ]
)
def test_database_failing_parameter(parameters, tmpdir):
    # test database arguments

    cmd_line = ['bin/platon', '--output', tmpdir] + parameters + ['test/data/mock-sample.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0


@pytest.mark.parametrize(
    'env_key,env_value',
    [
        ('foo', ''),  # not provided
        ('PLATON_DB', ''),  # missing path
        ('PLATON_DB', 'test/foo')  # not existing path
    ]
)
def test_database_failing_environment(env_key, env_value, tmpdir):
    # test database arguments

    env = os.environ
    env[env_key] = env_value
    cmd_line = ['bin/platon', '--output', tmpdir, 'test/data/mock-sample.fna']
    proc = run(cmd_line, env=env)
    assert proc.returncode != 0


def test_output_failing():
    # test database arguments
    cmd_line = ['bin/platon', '--output', '/', 'test/data/mock-sample.fna']
    proc = run(cmd_line)
    assert proc.returncode != 0