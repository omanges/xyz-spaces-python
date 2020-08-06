# Copyright (C) 2019-2020 HERE Europe B.V.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# SPDX-License-Identifier: Apache-2.0
# License-Filename: LICENSE

"""Module for providing test fixtures for the Hub API tests."""

import json
from pathlib import Path

import pytest

from xyzspaces.apis import HubApi
from xyzspaces.datasets import get_chicago_parks_data, get_countries_data
from xyzspaces.spaces import Space
from xyzspaces.utils import get_xyz_token

XYZ_TOKEN = get_xyz_token()


@pytest.fixture()
def api():
    """Create shared XYZ Hub Api instance as a pytest fixture."""
    api = HubApi(credentials=XYZ_TOKEN)
    return api


@pytest.fixture()
def space_id():
    """Create shared XYZ space with countries data as a pytest fixture."""
    api = HubApi(credentials=XYZ_TOKEN)

    # setup, create temporary space
    res = api.post_space(
        data={
            "title": "Testing xyzspaces",
            "description": "Temporary space containing countries data.",
        }
    )
    space_id = res["id"]

    # add features to space
    gj_countries = get_countries_data()
    api.put_space_features(space_id=space_id, data=gj_countries)

    yield space_id

    # now teardown (delete temporary space)
    api.delete_space(space_id=space_id)


@pytest.fixture()
def empty_space(api):
    """Create shared empty XYZ space as a pytest fixture."""
    # setup, create temporary space
    space = Space(api=api).new(
        title="Testing xyzspaces",
        description="Temporary empty space containing no features.",
    )

    yield space

    # now teardown (delete temporary space)
    space.delete()


@pytest.fixture()
def space_object(api, space_id):
    """Create from an existing space ID."""
    space = Space.from_id(space_id)
    return space


@pytest.fixture()
def upstream_spaces(api, space_id):
    """Create a list of space_ids, to test virtual-spaces."""
    res = api.post_space(
        data={
            "title": "Testing xyzspaces with Chicago Parks data",
            "description": "Temporary space containing Chicago parks data",
        }
    )
    space_id2 = res["id"]

    # add features to space
    gj_countries = get_chicago_parks_data()
    api.put_space_features(space_id=space_id2, data=gj_countries)

    yield [space_id, space_id2]

    # now teardown (delete temporary spaces)
    api.delete_space(space_id=space_id2)


@pytest.fixture()
def shared_space():
    """Create a new shared space."""
    space = Space.new(
        title="test shared space", description="test shared space", shared=True
    )

    yield space

    # now teardown (delete temporary space)
    space.delete()


@pytest.fixture()
def large_data_space():
    """Create a new large data space."""
    space = Space.new(
        title="test large data space",
        description="test large data space",
        shared=True,
    )
    path = Path(__file__).parent.parent / "data" / "large_data.json"

    with open(path) as json_file:
        data = json.load(json_file)

    space.add_features(data)

    yield space

    # now teardown (delete temporary space)
    space.delete()
