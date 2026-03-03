"""Pipeline wrapper for a pure-Python m-CBAC workflow."""

from __future__ import annotations

import math
import re
import shutil
import tempfile
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path


@dataclass(frozen=True)
class AtomRecord:
    label: str
    symbol: str
    frac_x_str: str
    frac_y_str: str
    frac_z_str: str
    frac_x: float
    frac_y: float
    frac_z: float


@dataclass(frozen=True)
class CifModel:
    path: Path
    header_lines: tuple[str, ...]
    atoms: tuple[AtomRecord, ...]
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    volume: float


@dataclass(frozen=True)
class PipelineResult:
    input_count: int
    final_files: tuple[Path, ...]
    work_dir: Path


# Legacy covalent radii map from secondlayer.f90. Falls back to 1.5 for unknown symbols.
RADII = {
    "Zn": 1.193, "C": 0.757, "Cl": 1.044, "Br": 1.192, "F": 0.668, "H": 0.354,
    "I": 1.382, "N": 0.700, "O": 0.634, "P": 1.101, "S": 1.064, "Ag": 1.386,
    "Al": 1.244, "As": 1.211, "Au": 1.262, "B": 0.838, "Ba": 2.277, "Be": 1.074,
    "Bi": 1.512, "Ca": 1.761, "Cd": 1.403, "Ce": 1.841, "Co": 1.241, "Cr": 1.345,
    "Cu": 1.302, "Dy": 1.710, "Er": 1.673, "Fe": 1.335, "Ga": 1.260, "Gd": 1.735,
    "Ge": 1.197, "Hf": 1.611, "Hg": 1.340, "Ho": 1.696, "In": 1.459, "Ir": 1.371,
    "K": 1.953, "La": 1.943, "Li": 1.336, "Lu": 1.671, "Mg": 1.421, "Mn": 1.382,
    "Mo": 1.470, "Na": 1.539, "Nb": 1.473, "Nd": 1.816, "Ni": 1.164, "Np": 1.666,
    "Pb": 1.459, "Pd": 1.338, "Pr": 1.823, "Pt": 1.364, "Pu": 1.657, "Rb": 2.260,
    "Re": 1.343, "Rh": 1.332, "Ru": 1.478, "Sb": 1.407, "Sc": 1.513, "Se": 1.190,
    "Si": 1.117, "Sm": 1.780, "Sn": 1.398, "Sr": 2.052, "Te": 1.386, "Th": 1.721,
    "Ti": 1.412, "Tm": 1.660, "U": 1.684, "V": 1.402, "W": 1.392, "Y": 1.698,
    "Yb": 1.637, "He": 0.849, "Ne": 0.920, "Zr": 1.564, "Tc": 1.322, "Xe": 1.267,
    "Cs": 2.570, "Pm": 1.801, "Eu": 1.771, "Tb": 1.732, "Ta": 1.511, "Os": 1.372,
    "Po": 1.500, "At": 1.545, "Rn": 1.420, "Fr": 2.880, "Ra": 2.512, "Ac": 1.983,
    "Pa": 1.711, "Am": 1.660, "Cm": 1.801, "Bk": 1.761, "Cf": 1.750, "Es": 1.724,
    "Fm": 1.712, "Md": 1.689, "No": 1.679, "Lw": 1.698,
}


def _asset_dir() -> Path:
    return Path(str(files("mcbac").joinpath("assets")))


def _parse_float_token(raw: str) -> float:
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", raw)
    if not match:
        raise ValueError(f"Cannot parse float from token: {raw}")
    return float(match.group(0))


def _load_database_0(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 2:
            out[parts[0]] = float(parts[1])
    return out


def _load_database_1(path: Path) -> dict[tuple[str, str], float]:
    out: dict[tuple[str, str], float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 3:
            out[(parts[0], parts[1])] = float(parts[2])
    return out


def _load_database_2(path: Path) -> dict[tuple[str, str, str], float]:
    out: dict[tuple[str, str, str], float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 4:
            out[(parts[0], parts[1], parts[2])] = float(parts[3])
    return out



def _computed_volume_from_cell(params: dict[str, float]) -> float:
    rad = math.pi / 180.0
    a = params["_cell_length_a"]
    b = params["_cell_length_b"]
    c = params["_cell_length_c"]
    alpha = params["_cell_angle_alpha"]
    beta = params["_cell_angle_beta"]
    gamma = params["_cell_angle_gamma"]

    ca = math.cos(alpha * rad)
    cb = math.cos(beta * rad)
    cg = math.cos(gamma * rad)
    return a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)


def _find_atom_loop(lines: list[str]) -> tuple[int, int, list[str], int, int]:
    for idx, line in enumerate(lines):
        if line.strip() != "loop_":
            continue

        header_start = idx + 1
        headers: list[str] = []
        cursor = header_start
        while cursor < len(lines) and lines[cursor].lstrip().startswith("_"):
            headers.append(lines[cursor].strip())
            cursor += 1

        if not headers:
            continue

        needed = {"_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"}
        if not needed.issubset(set(headers)):
            continue

        if "_atom_site_label" not in headers and "_atom_site_type_symbol" not in headers:
            continue

        data_start = cursor
        data_end = data_start
        while data_end < len(lines):
            s = lines[data_end].strip()
            if not s:
                break
            if s == "loop_" or s.startswith("_") or s.startswith("data_") or s.startswith("#"):
                break
            data_end += 1

        return header_start, data_start, headers, data_start, data_end

    raise ValueError("No atom loop found in CIF")


def _parse_cif(path: Path) -> CifModel:
    lines = path.read_text(encoding="utf-8").splitlines()

    params: dict[str, float] = {}
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        key = parts[0]
        if key in {
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
            }:
            params[key] = _parse_float_token(parts[1])

    required = [
        "_cell_length_a",
        "_cell_length_b",
        "_cell_length_c",
        "_cell_angle_alpha",
        "_cell_angle_beta",
        "_cell_angle_gamma",
    ]
    missing = [k for k in required if k not in params]
    if missing:
        raise ValueError(f"Missing CIF parameters in {path.name}: {', '.join(missing)}")

    header_start, _data_start_unused, headers, data_start, data_end = _find_atom_loop(lines)
    col_index = {name: i for i, name in enumerate(headers)}

    label_idx = col_index.get("_atom_site_label")
    symbol_idx = col_index.get("_atom_site_type_symbol", label_idx)
    x_idx = col_index["_atom_site_fract_x"]
    y_idx = col_index["_atom_site_fract_y"]
    z_idx = col_index["_atom_site_fract_z"]

    atoms: list[AtomRecord] = []
    for row in lines[data_start:data_end]:
        parts = row.split()
        needed_idx = max(label_idx or 0, symbol_idx or 0, x_idx, y_idx, z_idx)
        if len(parts) <= needed_idx:
            continue
        label = parts[label_idx] if label_idx is not None else parts[symbol_idx]
        symbol = parts[symbol_idx] if symbol_idx is not None else label
        x_raw = parts[x_idx]
        y_raw = parts[y_idx]
        z_raw = parts[z_idx]
        atoms.append(
            AtomRecord(
                label=label,
                symbol=symbol,
                frac_x_str=x_raw,
                frac_y_str=y_raw,
                frac_z_str=z_raw,
                frac_x=_parse_float_token(x_raw),
                frac_y=_parse_float_token(y_raw),
                frac_z=_parse_float_token(z_raw),
            )
        )

    if not atoms:
        raise ValueError(f"No atom rows found in {path.name}")

    return CifModel(
        path=path,
        header_lines=tuple(lines[:header_start]),
        atoms=tuple(atoms),
        a=params["_cell_length_a"],
        b=params["_cell_length_b"],
        c=params["_cell_length_c"],
        alpha=params["_cell_angle_alpha"],
        beta=params["_cell_angle_beta"],
        gamma=params["_cell_angle_gamma"],
        volume=params.get("_cell_volume", _computed_volume_from_cell(params)),
    )


def _lattice_vectors(model: CifModel) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    rad = math.pi / 180.0
    a1 = model.a
    b1 = model.b * math.cos(model.gamma * rad)
    c1 = model.c * math.cos(model.beta * rad)
    a2 = 0.0
    b2 = model.b * math.sin(model.gamma * rad)
    c2 = model.c * ((math.cos(model.alpha * rad) - math.cos(model.beta * rad) * math.cos(model.gamma * rad)) / math.sin(model.gamma * rad))
    a3 = 0.0
    b3 = 0.0
    c3 = model.volume / (model.a * model.b * math.sin(model.gamma * rad))
    return (a1, b1, c1), (a2, b2, c2), (a3, b3, c3)


def _frac_to_cart(frac: tuple[float, float, float], vecs: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]) -> tuple[float, float, float]:
    (a1, b1, c1), (a2, b2, c2), (a3, b3, c3) = vecs
    x, y, z = frac
    return (
        x * a1 + y * b1 + z * c1,
        x * a2 + y * b2 + z * c2,
        x * a3 + y * b3 + z * c3,
    )


def _compute_types(model: CifModel) -> tuple[list[list[str]], list[list[int]]]:
    vecs = _lattice_vectors(model)
    n_atoms = len(model.atoms)
    frac_coords = [(a.frac_x, a.frac_y, a.frac_z) for a in model.atoms]
    symbols = [a.symbol for a in model.atoms]

    cart = [_frac_to_cart(fr, vecs) for fr in frac_coords]

    expanded_cart: list[tuple[float, float, float]] = []
    expanded_symbol: list[str] = []
    for ii in (-1, 0, 1):
        for jj in (-1, 0, 1):
            for kk in (-1, 0, 1):
                for idx in range(n_atoms):
                    fx, fy, fz = frac_coords[idx]
                    expanded_cart.append(_frac_to_cart((fx + ii, fy + jj, fz + kk), vecs))
                    expanded_symbol.append(symbols[idx])

    all_ranks: list[list[str]] = []
    all_nums: list[list[int]] = []

    for j in range(n_atoms):
        distances = [100.0] * 8
        ranks = ["null"] * 8
        nums = [0] * 8
        r1 = RADII.get(symbols[j], 1.5)

        cx, cy, cz = cart[j]
        for z, (ex, ey, ez) in enumerate(expanded_cart):
            dx = abs(cx - ex)
            dy = abs(cy - ey)
            dz = abs(cz - ez)
            dis = math.sqrt(dx * dx + dy * dy + dz * dz)

            r2 = RADII.get(expanded_symbol[z], 1.5)
            r_cut = (r1 + r2) * 1.25
            if dis >= r_cut:
                continue

            slot = -1
            for k in range(8):
                if dis < distances[k]:
                    slot = k
                    break
            if slot < 0:
                continue

            for m in range(7, slot, -1):
                distances[m] = distances[m - 1]
                ranks[m] = ranks[m - 1]
                nums[m] = nums[m - 1]

            distances[slot] = dis
            ranks[slot] = expanded_symbol[z]
            nums[slot] = (z % n_atoms) + 1

        all_ranks.append(ranks)
        all_nums.append(nums)

    return all_ranks, all_nums


def _sorted_concat(tokens: list[str]) -> str:
    if not tokens:
        return ""
    return "".join(sorted(tokens))


def _is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def _assign_charges(model: CifModel, db0: dict[str, float], db1: dict[tuple[str, str], float], db2: dict[tuple[str, str, str], float]) -> list[float]:
    ranks, nums = _compute_types(model)
    n_atoms = len(model.atoms)

    zero_stage: list[str] = []

    for i in range(n_atoms):
        m = sum(1 for x in nums[i] if x != 0)

        key0 = ranks[i][0]
        first_neighbors = ranks[i][1:m]
        key1 = _sorted_concat(first_neighbors)

        second_tokens: list[str] = []
        for idx in range(1, m):
            neighbor_id = nums[i][idx] - 1
            second_tokens.extend(ranks[neighbor_id][1:8])
        key2 = _sorted_concat(second_tokens)

        sec = db2.get((key0, key1, key2))
        sec_text = str(sec) if sec is not None else f"{key0} {key1}"

        sec_parts = sec_text.split()
        first = db1.get((sec_parts[0], sec_parts[1])) if len(sec_parts) >= 2 else None
        first_text = str(first) if first is not None else sec_parts[0]

        zero = db0.get(first_text)
        zero_text = str(zero) if zero is not None else first_text
        zero_stage.append(zero_text)

    known_sum = sum(float(v) for v in zero_stage if _is_number(v))
    unknown_count = sum(1 for v in zero_stage if not _is_number(v))

    if unknown_count > 0:
        unknown_charge = (0.0 - known_sum) / unknown_count
        assigned = [float(v) if _is_number(v) else unknown_charge for v in zero_stage]
    else:
        assigned = [float(v) for v in zero_stage]

    abs_sum = sum(abs(v) for v in assigned)
    real_sum = sum(assigned)
    if abs_sum == 0.0:
        return assigned

    return [v - (real_sum * abs(v) / abs_sum) for v in assigned]


def _write_final_cif(model: CifModel, charges: list[float], output_path: Path) -> None:
    header = list(model.header_lines)
    header.extend(
        [
            " _atom_site_label",
            " _atom_site_type_symbol",
            " _atom_site_fract_x",
            " _atom_site_fract_y",
            " _atom_site_fract_z",
            " _atom_site_charge",
        ]
    )

    lines = header[:]
    for atom, charge in zip(model.atoms, charges):
        lines.append(
            f"{atom.label} {atom.symbol} {atom.frac_x_str} {atom.frac_y_str} {atom.frac_z_str} {charge:.12g}"
        )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(
    data_dir: Path,
    output_dir: Path,
    log_path: Path | None = None,
    work_dir: Path | None = None,
    keep_work_dir: bool = False,
    gfortran: str = "gfortran",
) -> PipelineResult:
    """Run the m-CBAC pipeline over CIF files in ``data_dir`` using pure Python."""
    _ = gfortran  # retained for CLI compatibility; no longer used.

    data_dir = data_dir.resolve()
    output_dir = output_dir.resolve()

    if not data_dir.exists() or not data_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {data_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    created_temp_workdir = work_dir is None
    run_dir = Path(tempfile.mkdtemp(prefix="mcbac-")) if created_temp_workdir else work_dir.resolve()
    if not created_temp_workdir:
        run_dir.mkdir(parents=True, exist_ok=True)

    asset_dir = _asset_dir()
    db0 = _load_database_0(asset_dir / "Database_0th.txt")
    db1 = _load_database_1(asset_dir / "Database_1st.txt")
    db2 = _load_database_2(asset_dir / "Database_2nd.txt")

    log_lines: list[str] = []
    produced: list[Path] = []

    cif_files = sorted(data_dir.glob("*.cif"))
    try:
        for cif_path in cif_files:
            model = _parse_cif(cif_path)
            charges = _assign_charges(model, db0=db0, db1=db1, db2=db2)

            run_final = run_dir / f"FINAL_{cif_path.name}"
            _write_final_cif(model, charges, run_final)

            final_out = output_dir / run_final.name
            shutil.copy2(run_final, final_out)
            produced.append(final_out)
            log_lines.append(f"Processed {cif_path.name} -> {final_out.name} ({len(model.atoms)} atoms)")

        if log_path is not None:
            log_path = log_path.resolve()
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text("\n".join(log_lines) + ("\n" if log_lines else ""), encoding="utf-8")

        return PipelineResult(input_count=len(cif_files), final_files=tuple(produced), work_dir=run_dir)
    finally:
        if created_temp_workdir and not keep_work_dir:
            shutil.rmtree(run_dir, ignore_errors=True)