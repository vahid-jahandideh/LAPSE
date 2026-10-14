"""Validation and provenance shared by the command-line tools."""

import hashlib
import importlib.metadata
import math
import platform
import re
from dataclasses import asdict


def parse_params(text):
    params = {}
    for item in text.split(",") if text else []:
        key, sep, value = item.strip().partition("=")
        if not sep or not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", key):
            raise ValueError("parameters must be comma-separated name=integer pairs")
        if key in params:
            raise ValueError(f"duplicate parameter: {key}")
        params[key] = int(value, 0)
    return params


def parse_roles(text):
    if text is None:
        return None
    roles = set(text.split(","))
    roles = {role.strip() for role in roles}
    if not roles <= {"input", "output", "rand", "temp", "const"}:
        raise ValueError("roles must be input, output, rand, temp, or const")
    return roles


def finite_probability(value):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("probabilities must be finite and in [0, 1]")
    return value


def clean_json(value):
    """Use JSON null for unavailable/nonfinite statistics, never NaN/Infinity."""
    if isinstance(value, dict):
        return {key: clean_json(val) for key, val in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean_json(val) for val in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def result_dict(result):
    return clean_json(asdict(result))


def provenance(circuit_text, params, options):
    from lapse import __version__
    versions = {}
    for name in ("lark", "numpy", "pandas", "matplotlib"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return {
        "lapse_version": __version__,
        "python": platform.python_version(),
        "dependencies": versions,
        "circuit_sha256": hashlib.sha256(circuit_text.encode("utf-8")).hexdigest(),
        "parameters": params,
        "options": options,
    }
