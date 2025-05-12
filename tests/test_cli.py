import os
import tempfile
from promb.cli import main

RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')


def test_cli_exact_antibody():
    with tempfile.TemporaryDirectory() as temp_dir:
        main([
            'exact',
            '-d', 'human-oas',
            os.path.join(RESOURCES_DIR, 'test_antibody.fa'),
            '-o', os.path.join(temp_dir, 'test_antibody.csv')
        ])
        with open(os.path.join(temp_dir, 'test_antibody.csv'), 'r') as f:
            result = f.read()
        with open(os.path.join(RESOURCES_DIR, 'test_antibody_exact.csv'), 'r') as f:
            expected = f.read()
        print(result)
        assert result == expected


def test_cli_nearest_antibody():
    with tempfile.TemporaryDirectory() as temp_dir:
        main([
            'nearest',
            '-d', 'human-oas',
            os.path.join(RESOURCES_DIR, 'test_antibody.fa'),
            '-o', os.path.join(temp_dir, 'test_antibody.csv')
        ])
        with open(os.path.join(temp_dir, 'test_antibody.csv'), 'r') as f:
            result = f.read()
        with open(os.path.join(RESOURCES_DIR, 'test_antibody_nearest.csv'), 'r') as f:
            expected = f.read()
        print(result)
        assert result == expected

