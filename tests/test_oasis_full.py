from promb import init_db
import pytest

from promb.db import PrombDB


@pytest.fixture
def human_db() -> PrombDB:
    return init_db('human-oas')


def test_compute_pmb_with_point_mutants(human_db: PrombDB):
    pmb = human_db.compute_pmb('QSGAELARPGASVKMSCKAS')
    assert pmb == 7 / 12 # even point mutants out of 17 peptides

def test_compute_pmb_with_double_mutants(human_db: PrombDB):
    pmb = human_db.compute_pmb('QVQLQQSGAELARPGASVKMSCKAS')
    assert pmb == (8*1 + 1*2) / 17 # eight point mutants + one double mutant out of 17 peptides

def test_compute_pmb_multiple_mutants(human_db: PrombDB):
    pmb = human_db.compute_pmb('ELVISISALIVE')
    assert pmb == (2 + 3 + 3 + 3) / 4

def test_compute_peptide_wise_pmb(human_db: PrombDB):
    peptides = human_db.chop_seq_peptides('QVQLQQSGAELARPGASVKMSCKAS')
    peptide_wise_pmb = human_db.compute_peptide_wise_mutations(peptides)
    # human QVQLQQSGA
    # human VQLQQSGAE
    # human QLQQSGAEL
    # point LQQSGAELA
    # double QQSGAELAR
    # point QSGAELARP
    # point SGAELARPG
    # point GAELARPGA
    # point AELARPGAS
    # point ELARPGASV
    # point LARPGASVK
    # point ARPGASVKM
    # human RPGASVKMS
    # human PGASVKMSC
    # human GASVKMSCK
    # human ASVKMSCKA
    # human SVKMSCKAS
    assert peptide_wise_pmb == [0, 0, 0, 1, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]

def test_compute_pmb_with_double_mutants_and_max_1(human_db: PrombDB):
    pmb = human_db.compute_pmb('QVQLQQSGAELARPGASVKMSCKAS', max=1)
    assert pmb == (8*1 + 1*1) / 17 # eight point mutants + one double mutant out of 17 peptides

def test_compute_pmb_with_double_mutants_and_max_2(human_db: PrombDB):
    pmb = human_db.compute_pmb('QVQLQQSGAELARPGASVKMSCKAS', max=2)
    assert pmb == (8*1 + 1*2) / 17 # eight point mutants + one double mutant out of 17 peptides

def test_compute_pmb_with_random_and_max_3(human_db: PrombDB):
    pmb = human_db.compute_pmb('IAFHSGOIAHOGIHAOSIGHAOSIGH', max=3)
    assert pmb == 3.0
