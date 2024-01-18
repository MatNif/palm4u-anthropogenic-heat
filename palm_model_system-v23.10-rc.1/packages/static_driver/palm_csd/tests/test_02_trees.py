"""Run tests for the DomainTree class."""
import copy
import lzma
import pickle
from typing import Generator, List, Optional, Tuple
import numpy.ma as ma
import pytest

from palm_csd.csd_config import (
    CSDConfigDomain,
    CSDConfigSettings,
)
from palm_csd.vegetation import DomainTree


@pytest.fixture(scope="module")
def config() -> Generator[CSDConfigDomain, None, None]:
    """Create a domain configuration with remove_low_lai_tree=False.
    Reset counter in cleanup to not interfere with other tests.
    """
    config = CSDConfigDomain(
        name="test",
        pixel_size=5,
        input_lower_left_x=2,
        input_lower_left_y=2,
        nx=5,
        ny=5,
        dz=5,
        remove_low_lai_tree=False,
    )
    yield config
    CSDConfigDomain._reset_counter()


@pytest.fixture(scope="module")
def config_remove_low_lai(config: CSDConfigDomain) -> CSDConfigDomain:
    """Copy a domain configuration and set remove_low_lai_tree=True."""
    config_remove_low_lai = copy.copy(config)
    config_remove_low_lai.remove_low_lai_tree = True
    return config_remove_low_lai


@pytest.fixture(scope="module")
def settings() -> Generator[CSDConfigSettings, None, None]:
    """Create a settings configuration.
    Reset counter in cleanup to not interfere with other tests.
    """
    settings = CSDConfigSettings(lai_tree_lower_threshold=0.5)
    yield settings
    CSDConfigSettings._reset_counter()


@pytest.fixture(scope="module")
def domaintree_defaults() -> Generator[None, None, None]:
    """Read in default values for trees and reset in cleanup."""
    DomainTree.populate_defaults()
    yield
    DomainTree.defaults = []


@pytest.fixture(scope="module")
def trees_valid(
    domaintree_defaults: None, config: CSDConfigDomain, settings: CSDConfigSettings
) -> List[Optional[DomainTree]]:
    """Create different valid trees."""
    trees_valid = []

    # trees with different shapes
    trees_valid.append(
        DomainTree.generate_tree(
            i=2,
            j=2,
            type=4,
            shape=1,
            height=14.0,
            lai=4.0,
            crown_diameter=15.0,
            trunk_diameter=ma.masked,
            settings=settings,
            config=config,
        )
    )
    trees_valid.append(
        DomainTree.generate_tree(
            i=3,
            j=1,
            type=32,
            shape=2,
            height=12.0,
            lai=2.0,
            crown_diameter=ma.masked,
            trunk_diameter=0.5,
            settings=settings,
            config=config,
        )
    )
    trees_valid.append(
        DomainTree.generate_tree(
            i=5,
            j=4,
            type=28,
            shape=3,
            height=6.0,
            lai=ma.masked,
            crown_diameter=4.5,
            trunk_diameter=0.7,
            settings=settings,
            config=config,
        )
    )
    trees_valid.append(
        DomainTree.generate_tree(
            i=3,
            j=3,
            type=73,
            shape=4,
            height=ma.masked,
            lai=5.0,
            crown_diameter=10.5,
            trunk_diameter=1.0,
            settings=settings,
            config=config,
        )
    )
    trees_valid.append(
        DomainTree.generate_tree(
            i=3,
            j=3,
            type=74,
            shape=ma.masked,
            height=14.7,
            lai=7.3,
            crown_diameter=14.5,
            trunk_diameter=1.5,
            settings=settings,
            config=config,
        )
    )
    trees_valid.append(
        DomainTree.generate_tree(
            i=3,
            j=3,
            type=82,
            shape=6,
            height=ma.masked,
            lai=ma.masked,
            crown_diameter=ma.masked,
            trunk_diameter=ma.masked,
            settings=settings,
            config=config,
        )
    )

    # tree with only masked input values
    trees_valid.append(
        DomainTree.generate_tree(
            i=2,
            j=3,
            type=ma.masked,
            shape=ma.masked,
            height=ma.masked,
            lai=ma.masked,
            crown_diameter=ma.masked,
            trunk_diameter=ma.masked,
            settings=settings,
            config=config,
        )
    )

    # tree with too low LAI
    trees_valid.append(
        DomainTree.generate_tree(
            i=2,
            j=3,
            type=86,
            shape=ma.masked,
            height=ma.masked,
            lai=0.4,
            crown_diameter=ma.masked,
            trunk_diameter=ma.masked,
            settings=settings,
            config=config,
        )
    )

    # uncomment to update reference file
    # write_trees_fields_to_file(trees_valid, config)

    return trees_valid


@pytest.fixture(scope="module")
def trees_none(
    domaintree_defaults: None, config_remove_low_lai: CSDConfigDomain, settings: CSDConfigSettings
) -> List[Optional[DomainTree]]:
    """Try to create different invalid trees resulting in None."""

    trees_none = []

    # tree_height too low
    trees_none.append(
        DomainTree.generate_tree(
            i=2,
            j=2,
            type=4,
            shape=1,
            height=2.0,
            lai=4.0,
            crown_diameter=15.0,
            trunk_diameter=ma.masked,
            settings=settings,
            config=config_remove_low_lai,
        )
    )

    # tree_lai too low and remove_low_lai_tree=True
    trees_none.append(
        DomainTree.generate_tree(
            i=2,
            j=2,
            type=4,
            shape=1,
            height=2.0,
            lai=4.0,
            crown_diameter=15.0,
            trunk_diameter=ma.masked,
            settings=settings,
            config=config_remove_low_lai,
        )
    )

    return trees_none


@pytest.fixture(scope="module")
def trees_fields_reference() -> (
    Tuple[List[ma.MaskedArray], List[ma.MaskedArray], List[ma.MaskedArray], List[ma.MaskedArray]]
):
    """Read in trees fields reference from file."""
    with lzma.open("tests/02_trees/output/trees.pkl.xz", "rb") as f:
        lads = pickle.load(f)
        bads = pickle.load(f)
        tree_ids = pickle.load(f)
        tree_types = pickle.load(f)
    return lads, bads, tree_ids, tree_types


def calculate_fields_for_trees(
    trees: List[DomainTree], config: CSDConfigDomain
) -> Tuple[List[ma.MaskedArray], List[ma.MaskedArray], List[ma.MaskedArray], List[ma.MaskedArray]]:
    """Generate fields for trees."""

    # lists for the 3d fields
    lads: List[ma.MaskedArray] = []
    bads: List[ma.MaskedArray] = []
    tree_ids: List[ma.MaskedArray] = []
    tree_types: List[ma.MaskedArray] = []

    # single tree 3d fields
    dimensions = (7, config.ny, config.nx)
    lad = ma.masked_all(dimensions)
    bad = ma.masked_all(dimensions)
    tree_id = ma.masked_all(dimensions)
    tree_type = ma.masked_all(dimensions)

    for tree in trees:
        # single tree 3d fields
        tree.generate_store_3d_fields(lad, bad, tree_id, tree_type, config)

        # append to lists
        lads.append(lad.copy())
        bads.append(bad.copy())
        tree_ids.append(tree_id.copy())
        tree_types.append(tree_type.copy())

        # reset single tree 3d fields
        lad.mask = True
        bad.mask = True
        tree_id.mask = True
        tree_type.mask = True

    return lads, bads, tree_ids, tree_types


def write_trees_fields_to_file(trees: List[DomainTree], config: CSDConfigDomain) -> None:
    """Write the trees fields to a file."""
    lads, bads, tree_ids, tree_types = calculate_fields_for_trees(trees, config)

    with lzma.open("tests/02_trees/output/trees.pkl.xz", "wb", preset=9) as f:
        pickle.dump(lads, f)
        pickle.dump(bads, f)
        pickle.dump(tree_ids, f)
        pickle.dump(tree_types, f)


@pytest.mark.usefixtures("domaintree_defaults")
def test_defaults() -> None:
    """Test default values for trees."""
    # 86 species + 1 default
    assert len(DomainTree.defaults) == 87


def test_valid_tree_objects(trees_valid: List[Optional[DomainTree]]) -> None:
    """Test valid tree objects."""

    # check tree values with some values from the defaults
    tree = trees_valid[0]
    assert tree is not None
    assert tree.type == 4
    assert tree.shape == 1
    assert tree.height == 14.0
    assert tree.lai == 4.0
    assert tree.crown_diameter == 15.0
    assert tree.trunk_diameter == 1.30

    tree = trees_valid[1]
    assert tree is not None
    assert tree.type == 32
    assert tree.shape == 2
    assert tree.height == 12.0
    assert tree.lai == 2.0
    assert tree.crown_diameter == 4.5
    assert tree.trunk_diameter == 0.5

    tree = trees_valid[2]
    assert tree is not None
    assert tree.type == 28
    assert tree.shape == 3
    assert tree.height == 6.0
    assert tree.lai == 3.0
    assert tree.crown_diameter == 4.5
    assert tree.trunk_diameter == 0.7

    tree = trees_valid[3]
    assert tree is not None
    assert tree.type == 73
    assert tree.shape == 4
    assert tree.height == 25.0
    assert tree.lai == 5.0
    assert tree.crown_diameter == 10.5
    assert tree.trunk_diameter == 1.0

    tree = trees_valid[4]
    assert tree is not None
    assert tree.type == 74
    assert tree.shape == 5
    assert tree.height == 14.7
    assert tree.lai == 7.3
    assert tree.crown_diameter == 14.5
    assert tree.trunk_diameter == 1.5

    tree = trees_valid[5]
    assert tree is not None
    assert tree.type == 82
    assert tree.shape == 6
    assert tree.height == 12.5
    assert tree.lai == 3.0
    assert tree.crown_diameter == 7.0
    assert tree.trunk_diameter == 0.7

    tree = trees_valid[6]
    assert tree is not None
    assert tree.type == 0
    assert tree.shape == 1
    assert tree.height == 12.0
    assert tree.lai == 3.0
    assert tree.crown_diameter == 4.0
    assert tree.trunk_diameter == 0.35

    # LAI should be that of type 86's summer LAI
    tree = trees_valid[7]
    assert tree is not None
    assert tree.type == 86
    assert tree.shape == 1
    assert tree.height == 5.0
    assert tree.lai == 3.0
    assert tree.crown_diameter == 5.0
    assert tree.trunk_diameter == 0.4


def test_none_tree_objects(
    trees_none: List[Optional[DomainTree]], config: CSDConfigDomain, settings: CSDConfigSettings
) -> None:
    """Test that invalid trees are set to None. Try to create trees that raise errors."""

    for tree in trees_none:
        assert tree is None

    # test shape >= 7
    with pytest.raises(ValueError):
        DomainTree.generate_tree(
            i=2,
            j=2,
            type=4,
            shape=7,
            height=2.0,
            lai=4.0,
            crown_diameter=15.0,
            trunk_diameter=1.0,
            settings=settings,
            config=config,
        )


def test_trees_fields_generator(
    trees_valid: List[DomainTree],
    config: CSDConfigDomain,
    trees_fields_reference: Tuple[
        List[ma.MaskedArray], List[ma.MaskedArray], List[ma.MaskedArray], List[ma.MaskedArray]
    ],
):
    """Test the fields generator."""

    lads, bads, tree_ids, tree_types = calculate_fields_for_trees(trees_valid, config)

    (
        lads_reference,
        bads_reference,
        tree_ids_reference,
        tree_types_reference,
    ) = trees_fields_reference

    # all lists should have the same length
    length = len(lads_reference)
    assert len(lads) == length
    assert len(bads) == length
    assert len(tree_ids) == length
    assert len(tree_types) == length

    for i in range(length):
        assert ma.allclose(lads[i], lads_reference[i], rtol=1e-14, atol=1e-14)
        assert ma.allclose(bads[i], bads_reference[i], rtol=1e-14, atol=1e-14)
        assert ma.allclose(tree_ids[i], tree_ids_reference[i], rtol=1e-14, atol=1e-14)
        assert ma.allclose(tree_types[i], tree_types_reference[i], rtol=1e-14, atol=1e-14)
