from promb import init_db
import pytest

from promb.db import PrombDB


@pytest.fixture
def human_db() -> PrombDB:
    return init_db('human-swissprot', peptide_length=9)


def test_is_human(human_db: PrombDB):
    assert human_db.contains('MEPLRLLIL') == True
    assert human_db.contains('ELVISLIVE') == False


def test_compute_peptide_content(human_db: PrombDB):
    score = human_db.compute_peptide_content('MEPLRLLILLFVTELSGAHNTTVFQGVAGQSLQVSCPYDSMKHWGRRKAWCRQLGEKGPC')
    assert score == 1.0


def test_find_peptides_matching_template(human_db: PrombDB):
    peptides = human_db.find_peptides_matching_template('TSP*YSPTS')
    assert peptides == [
        'TSPGYSPTS',
        'TSPKYSPTS',
        'TSPNYSPTS',
        'TSPSYSPTS',
        'TSPTYSPTS',
    ]
    peptides = human_db.find_peptides_matching_template('TSP*YSP*S')
    assert peptides == [
        'TSPGYSPTS',
        'TSPKYSPTS',
        'TSPNYSPTS',
        'TSPSYSPSS',
        'TSPSYSPTS',
        'TSPTYSPTS',
    ]


def test_find_point_mutant_peptides(human_db: PrombDB):
    peptides = human_db.find_point_mutant_peptides('TSPSYSPTS')
    assert peptides == [
        'MSPSYSPTS',
        'QSPSYSPTS',
        'TSPGYSPTS',
        'TSPKYSPTS',
        'TSPNYSPTS',
        'TSPTYSPTS',
        'TSPSYSPSS',
    ]


