from mne.bem import FIFF, read_surface, _ico_downsample, _check_complete_surface
import os
from pathlib import Path

def _surfaces_to_bem_no_size_check(
    surfs, ids, sigmas, ico=None, rescale=True, incomplete="raise", extra=""
):
    """Convert surfaces to a BEM."""
    # equivalent of mne_surf2bem
    # surfs can be strings (filenames) or surface dicts
    if len(surfs) not in (1, 3) or not (len(surfs) == len(ids) == len(sigmas)):
        raise ValueError(
            "surfs, ids, and sigmas must all have the same "
            "number of elements (1 or 3)"
        )
    for si, surf in enumerate(surfs):
        if isinstance(surf, (str, Path, os.PathLike)):
            surfs[si] = surf = read_surface(surf, return_dict=True)[-1]
    # Downsampling if the surface is isomorphic with a subdivided icosahedron
    if ico is not None:
        for si, surf in enumerate(surfs):
            surfs[si] = _ico_downsample(surf, ico)
    for surf, id_ in zip(surfs, ids):
        # Do topology checks (but don't save data) to fail early
        surf["id"] = id_
        _check_complete_surface(surf, copy=True, incomplete=incomplete, extra=extra)
        surf["coord_frame"] = surf.get("coord_frame", FIFF.FIFFV_COORD_MRI)
        surf.update(np=len(surf["rr"]), ntri=len(surf["tris"]))
        if rescale:
            surf["rr"] /= 1000.0  # convert to meters

    # Shifting surfaces is not implemented here...

    # Order the surfaces for the benefit of the topology checks
    for surf, sigma in zip(surfs, sigmas):
        surf["sigma"] = sigma
    surfs = _order_surfaces(surfs)

    # Check topology as best we can
    _check_surfaces(surfs, incomplete=incomplete)
    # for surf in surfs:  # REMOVED FOR BONES OF ARM!
    #     _check_surface_size(surf)
    _check_thicknesses(surfs)
    logger.info("Surfaces passed the basic topology checks.")
    return surfs