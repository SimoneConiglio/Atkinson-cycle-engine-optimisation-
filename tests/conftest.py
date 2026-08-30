"""Shared fixtures."""

from __future__ import annotations

import matplotlib
import pytest

matplotlib.use("Agg")

from exlink import PUBLISHED_DESIGN, Design, analyse
from exlink.model import Analysis
from exlink.reference import REFINED_DESIGN


@pytest.fixture(scope="session")
def published() -> Design:
    """The design vector tabulated in the 2015 report."""
    return PUBLISHED_DESIGN


@pytest.fixture(scope="session")
def refined() -> Design:
    """The feasible design this framework produces from the published one."""
    return REFINED_DESIGN


@pytest.fixture(scope="session")
def published_analysis() -> Analysis:
    """Analysis of the published design at a good resolution."""
    return analyse(PUBLISHED_DESIGN, samples=1440)


@pytest.fixture(scope="session")
def refined_analysis() -> Analysis:
    """Analysis of the refined design at a good resolution."""
    return analyse(REFINED_DESIGN, samples=1440)
