"""
NANOGrav DR15 Data Loader

Load public NANOGrav DR15 timing residuals and pulsar parameters.
Provides unified interface for residual-first audit.

Reference: 
  NANOGrav Collaboration (2023)
  "The NANOGrav 15 yr Data Set: Searching for Signals from Supermassive 
   Black Hole Binaries in Pulsar Timing Data"
  ApJ (The Astrophysical Journal)

Public data: https://data.nanoGrav.org/

Author: Dimitar Kretski
License: MIT
"""

import re
import warnings
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np

PSR_RE = re.compile(r"([JB]\d{4}[+-]\d{2,4})")
OBLIQUITY_DEG = 84381.406 / 3600.0      # IERS2010 obliquity (PINT default for ecliptic coords)


def normalize_psr_name(name: str) -> str:
    """Unicode minus → ASCII hyphen (NANOGrav file names use '-')."""
    return name.replace("\u2212", "-").strip()


def _sexagesimal_to_deg(value: str, hours: bool) -> float:
    parts = [float(p) for p in value.replace(":", " ").split()]
    sign = -1.0 if value.strip().startswith("-") else 1.0
    parts[0] = abs(parts[0])
    deg = parts[0] + (parts[1] if len(parts) > 1 else 0) / 60 + (parts[2] if len(parts) > 2 else 0) / 3600
    return sign * deg * (15.0 if hours else 1.0)


def ecliptic_to_equatorial(elong_deg: float, elat_deg: float) -> Tuple[float, float]:
    lam, beta, eps = np.radians([elong_deg, elat_deg, OBLIQUITY_DEG])
    sin_dec = np.sin(beta) * np.cos(eps) + np.cos(beta) * np.sin(eps) * np.sin(lam)
    dec = np.arcsin(sin_dec)
    ra = np.arctan2(np.sin(lam) * np.cos(eps) - np.tan(beta) * np.sin(eps), np.cos(lam))
    return float(np.degrees(ra) % 360.0), float(np.degrees(dec))


def read_par_file(path) -> Dict[str, float]:
    """Minimal .par reader: position (RAJ/DECJ or ELONG/ELAT), DM, F0, F1, PX."""
    raw = {}
    with open(path, errors="replace") as f:
        for line in f:
            tok = line.split()
            if len(tok) >= 2 and not tok[0].startswith("#"):
                raw.setdefault(tok[0].upper(), tok[1])

    def num(key, default=np.nan):
        try:
            return float(raw[key].replace("D", "E"))
        except (KeyError, ValueError):
            return default

    if "RAJ" in raw and "DECJ" in raw:
        ra, dec = _sexagesimal_to_deg(raw["RAJ"], True), _sexagesimal_to_deg(raw["DECJ"], False)
    elif "ELONG" in raw and "ELAT" in raw:
        ra, dec = ecliptic_to_equatorial(num("ELONG"), num("ELAT"))
    elif "LAMBDA" in raw and "BETA" in raw:
        ra, dec = ecliptic_to_equatorial(num("LAMBDA"), num("BETA"))
    else:
        ra = dec = np.nan
    f0, f1 = num("F0"), num("F1", 0.0)
    return {"ra": ra, "dec": dec, "dm": num("DM"), "parallax": num("PX", 0.0),
            "period": 1.0 / f0 if f0 == f0 else np.nan,
            "period_derivative": -f1 / f0 ** 2 if f0 == f0 else np.nan}


@dataclass
class PulsarTimingData:
    """Container for a single pulsar's timing data."""
    
    name: str                    # PSR name (e.g., "J0023+0923")
    toas: np.ndarray            # Times of arrival (MJD)
    residuals: np.ndarray       # Post-fit timing residuals (SECONDS)
    errors: np.ndarray          # Measurement uncertainties (SECONDS)
    frequencies: np.ndarray     # Observation frequencies (MHz)
    
    # Pulsar parameters (from best-fit timing model)
    ra: float                   # Right ascension (degrees)
    dec: float                  # Declination (degrees)
    parallax: float             # Parallax (mas)
    dm: float                   # Dispersion measure (pc/cm³)
    period: float               # Spin period (seconds)
    period_derivative: float    # Period derivative (s/s)
    
    # Metadata
    n_observations: int         # Number of TOAs
    time_span: float            # Observation span (years)
    cadence: float              # Mean TOA cadence (days)
    
    def __post_init__(self):
        """Validate data consistency."""
        if not (len(self.toas) == len(self.residuals) == len(self.errors)):
            raise ValueError("Mismatched array lengths")
        if len(self.toas) == 0:
            raise ValueError("Empty pulsar data")
    
    def residual_rms(self) -> float:
        """RMS of timing residuals."""
        return np.sqrt(np.mean(self.residuals**2))
    
    def weighted_rms(self) -> float:
        """Uncertainty-weighted RMS."""
        weights = 1.0 / (self.errors**2)
        return np.sqrt(np.sum(weights * self.residuals**2) / np.sum(weights))


class NANOGravDR15Loader:
    """
    Mock loader for NANOGrav DR15 data.
    
    In production, this would:
      1. Download from data.nanoGrav.org
      2. Parse TEMPO2 .par and .tim files
      3. Extract residuals from published posterior chains
    
    For now, demonstrates the interface and provides synthetic data
    consistent with published parameters.
    """
    
    # Published NANOGrav DR15 pulsar list and key parameters
    PULSARS_DR15 = {
        'J0023+0923': {
            'ra': 5.832,
            'dec': 9.383,
            'parallax': 0.57,
            'dm': 14.396,
            'period': 0.00344,
            'n_obs': 281,
            'time_span': 15.0,
        },
        'J0030+0451': {
            'ra': 7.639,
            'dec': 4.850,
            'parallax': 3.61,
            'dm': 4.326,
            'period': 0.00749,
            'n_obs': 302,
            'time_span': 15.0,
        },
        'J0340+4130': {
            'ra': 55.034,
            'dec': 41.506,
            'parallax': 0.33,
            'dm': 49.291,
            'period': 0.00346,
            'n_obs': 276,
            'time_span': 15.0,
        },
        'J0437−4715': {
            'ra': 69.314,
            'dec': -47.252,
            'parallax': 5.76,
            'dm': 2.6437,
            'period': 0.00575,
            'n_obs': 451,
            'time_span': 15.0,
        },
        'J0613−0200': {
            'ra': 93.356,
            'dec': -2.019,
            'parallax': 1.52,
            'dm': 38.881,
            'period': 0.00324,
            'n_obs': 322,
            'time_span': 15.0,
        },
        # ... add more as needed
    }
    
    SYNTH_START_MJD = 57000.0

    def __init__(self, use_synthetic: bool = True, seed: Optional[int] = None,
                 red_amplitude: float = 50e-9, red_index: float = 2.5,
                 common_amplitude: float = 3e-9, common_index: float = 13 / 3,
                 data_dir: str = "NANOGrav15yr", epoch_averaged: bool = True,
                 whitened: bool = False, column_map: Optional[Dict[str, int]] = None,
                 residual_unit: Optional[str] = None):
        """
        Initialize loader.

        Args:
            use_synthetic: If True, generate synthetic data. If False, attempt
                           to load real data from disk (not implemented yet).
            seed: RNG seed for reproducible synthetic data.
            red_amplitude / red_index: intrinsic red noise (RMS in s, β),
                           independent for every pulsar.
            common_amplitude / common_index: common red process (RMS in s, β).
                           ONE realisation shared by all pulsars (monopole,
                           i.e. clock-error-like — NOT Hellings–Downs correlated).
                           Set 0 for a pure-null ensemble.
            data_dir: folder produced by fetch_nanograv15.py (real data).
            epoch_averaged: use the epoch-averaged residual tables (recommended:
                           the full narrowband tables have ~10^4 TOAs per pulsar).
            whitened: use whitened residuals. Keep False — whitening removes the
                           red noise that the audit is supposed to find.
            column_map: override column detection, e.g.
                           {'mjd': 0, 'residual': 1, 'error': 2, 'frequency': 3}
            residual_unit: 's' or 'us'; None = read from header, else guess.
        """
        self.use_synthetic = use_synthetic
        self.data = {}
        self.rng = np.random.default_rng(seed)
        self.red_amplitude, self.red_index = red_amplitude, red_index
        self.common_amplitude, self.common_index = common_amplitude, common_index
        self._common_grid = None
        self._common_series = None
        self.data_dir = Path(data_dir)
        self.epoch_averaged = epoch_averaged
        self.whitened = whitened
        self.column_map = column_map
        self.residual_unit = residual_unit
    
    def load_pulsar(self, psr_name: str) -> PulsarTimingData:
        """
        Load timing data for a single pulsar.
        
        Args:
            psr_name: Pulsar name (e.g., "J0023+0923")
        
        Returns:
            PulsarTimingData object
        """
        if psr_name in self.data:
            return self.data[psr_name]

        if self.use_synthetic:
            if psr_name not in self.PULSARS_DR15:
                raise ValueError(f"Pulsar {psr_name} not in synthetic list")
            data = self._generate_synthetic_residuals(psr_name, self.PULSARS_DR15[psr_name])
        else:
            data = self._load_from_disk(normalize_psr_name(psr_name))
        
        self.data[psr_name] = data
        return data
    
    def load_all_pulsars(self) -> Dict[str, PulsarTimingData]:
        """Load all DR15 pulsars."""
        all_data = {}
        names = self.PULSARS_DR15.keys() if self.use_synthetic else self.available_pulsars()
        for psr_name in names:
            try:
                all_data[psr_name] = self.load_pulsar(psr_name)
            except Exception as e:
                warnings.warn(f"Failed to load {psr_name}: {e}")
        return all_data
    
    def _generate_synthetic_residuals(self, psr_name: str, params: Dict) -> PulsarTimingData:
        """
        Generate synthetic timing residuals consistent with NANOGrav parameters.
        
        Incorporates:
          1. White noise (measurement error)
          2. Red noise (spin-down, interstellar propagation)
          3. Common red noise (GW background, if present)
        
        Args:
            psr_name: Pulsar name
            params: Dictionary of pulsar parameters
        
        Returns:
            PulsarTimingData with synthetic residuals
        """
        n_obs = params['n_obs']
        time_span = params['time_span']
        rng = self.rng

        # Uniform cadence (real TOAs are irregular — kept simple here)
        toas = np.linspace(self.SYNTH_START_MJD, self.SYNTH_START_MJD + time_span * 365.25, n_obs)
        dt = np.diff(toas).mean()

        # White noise (measurement errors), seconds
        white_noise_level = 0.3e-6
        white_noise = rng.normal(0, white_noise_level, n_obs)
        errors = np.full(n_obs, white_noise_level)

        # Intrinsic red noise — independent per pulsar
        red_noise = self._generate_red_noise(n_obs, self.red_index, self.red_amplitude)

        # Common red process — the SAME realisation for every pulsar
        common_noise = self._common_signal(toas)

        residuals = white_noise + red_noise + common_noise
        frequencies = rng.uniform(400, 1400, n_obs)

        # Construct dataclass
        data = PulsarTimingData(
            name=psr_name,
            toas=toas,
            residuals=residuals,
            errors=errors,
            frequencies=frequencies,
            ra=params['ra'],
            dec=params['dec'],
            parallax=params['parallax'],
            dm=params['dm'],
            period=params['period'],
            period_derivative=1e-18,  # Typical for millisecond pulsars
            n_observations=n_obs,
            time_span=time_span,
            cadence=dt
        )
        
        return data
    
    def _generate_red_noise(self, n_obs: int, spectral_index: float, amplitude: float) -> np.ndarray:
        """
        Gaussian noise with PSD ∝ f^(-β), scaled to the given RMS [s].
        (Previous version set f[0] = 1e-10, which put ~100% of the power in the
        DC bin — the "red noise" was effectively a constant offset.)
        """
        if amplitude <= 0:
            return np.zeros(n_obs)
        f = np.fft.rfftfreq(n_obs)
        amp = np.zeros_like(f)
        amp[1:] = f[1:] ** (-spectral_index / 2.0)          # DC bin = 0
        spec = amp * (self.rng.standard_normal(len(f)) + 1j * self.rng.standard_normal(len(f)))
        x = np.fft.irfft(spec, n=n_obs)
        return amplitude * x / np.std(x)

    def _common_signal(self, toas: np.ndarray) -> np.ndarray:
        """One common realisation on a fine time grid, interpolated to each pulsar's TOAs."""
        if self.common_amplitude <= 0:
            return np.zeros(len(toas))
        if self._common_series is None:
            span = max(p['time_span'] for p in self.PULSARS_DR15.values()) * 365.25
            self._common_grid = np.linspace(self.SYNTH_START_MJD, self.SYNTH_START_MJD + span, 4096)
            self._common_series = self._generate_red_noise(4096, self.common_index,
                                                           self.common_amplitude)
        return np.interp(toas, self._common_grid, self._common_series)

    # ------------------------------------------------------------------ real data
    def _residual_files(self) -> List[Path]:
        """NANOGrav 15-yr v2.1.0: residuals/<PSR>_..._nb.avg.res and ..._nb.full.res"""
        if not self.data_dir.exists():
            raise FileNotFoundError(
                f"{self.data_dir} not found — run  python fetch_nanograv15.py  first")
        return sorted(p for p in self.data_dir.rglob("*.res") if p.is_file())

    def available_pulsars(self) -> List[str]:
        names = set()
        for p in self._residual_files():
            m = PSR_RE.search(p.name)
            if m:
                names.add(m.group(1))
        return sorted(names)

    def _pick_residual_file(self, psr: str) -> Path:
        tag = ".avg." if self.epoch_averaged else ".full."
        cands = [p for p in self._residual_files() if psr in p.name]
        good = [p for p in cands if tag in p.name]
        if not good:
            listing = "\n  ".join(str(p.relative_to(self.data_dir)) for p in cands) or "(none)"
            raise FileNotFoundError(f"No '{tag}res' file for {psr} in {self.data_dir}. "
                                    f"Candidates:\n  {listing}")
        return good[0]

    def _find_par(self, psr: str) -> Optional[Path]:
        pars = [p for p in self.data_dir.rglob("*.par") if psr in p.name]
        nb = [p for p in pars if "narrowband" in str(p).lower()]
        pars = nb or pars
        pref = [p for p in pars if not any(k in p.name.lower() for k in ("norednoise", "predictive"))]
        return (pref or pars or [None])[0]

    def _read_table(self, path: Path):
        header = None
        with open(path, errors="replace") as f:
            for line in f:
                s_ = line.strip()
                if not s_:
                    continue
                if s_.startswith("#"):
                    header = s_.lstrip("#").strip()
                    continue
                try:
                    [float(x) for x in s_.split()[:2]]
                    break                            # first data line reached
                except ValueError:
                    header = s_                      # un-commented header line
        data = np.genfromtxt(path, comments="#", invalid_raise=False)
        if data.ndim == 1:
            data = data[None, :]
        data = data[np.all(np.isfinite(data[:, :min(3, data.shape[1])]), axis=1)]
        cols = header.replace(",", " ").split() if header else []
        return data, cols, header or ""

    def _map_columns(self, cols: List[str], path: Path, header: str) -> Dict[str, int]:
        if self.column_map:
            return self.column_map
        low = [c.lower() for c in cols]

        def is_white(c):
            return "white" in c and not any(k in c for k in ("unwhite", "un-white", "un_white",
                                                             "nowhite", "not_white", "raw"))

        want_white = self.whitened
        rules = {
            "mjd": lambda c: "mjd" in c or c in ("toa", "toas", "time", "day"),
            "error": lambda c: any(k in c for k in ("err", "unc", "sigma", "uncert")),
            "residual": lambda c: (("res" in c) and not any(k in c for k in ("err", "unc", "sigma"))
                                   and is_white(c) == want_white),
            "frequency": lambda c: "freq" in c or "mhz" in c,
        }
        found = {}
        for key, rule in rules.items():
            hits = [i for i, c in enumerate(low) if rule(c) and i not in found.values()]
            if hits:
                found[key] = hits[0]
        missing = {"mjd", "residual", "error"} - set(found)
        if missing:
            raise ValueError(
                f"Cannot identify columns {sorted(missing)} in {path.name}.\n"
                f"Header seen: '{header}'.\n"
                f"Pass column_map={{'mjd':i,'residual':j,'error':k,'frequency':l}} to the loader, "
                f"or send the first lines of the file.")
        print(f"[{path.name}] columns used: " +
              ", ".join(f"{k}='{cols[i]}'(#{i})" for k, i in found.items()))
        return found

    def _unit_scale(self, header: str, err: np.ndarray) -> Tuple[float, str]:
        unit = self.residual_unit
        h = header.lower()
        if unit is None:
            if any(k in h for k in ("(us)", "[us]", "_us", " us", "µs", "usec", "micro")):
                unit = "us"
            elif any(k in h for k in ("(s)", "[s]", "sec)")):
                unit = "s"
        if unit is None:                       # heuristic: TOA errors are ~1e-7…1e-5 s
            unit = "us" if np.median(err) > 1e-3 else "s"
            warnings.warn(f"Residual unit not stated in header — assuming '{unit}' "
                          f"(median error = {np.median(err):.3g}). Override with residual_unit=.")
        return (1e-6 if unit == "us" else 1.0), unit

    def _load_from_disk(self, psr_name: str) -> PulsarTimingData:
        path = self._pick_residual_file(psr_name)
        data, cols, header = self._read_table(path)
        cmap = self._map_columns(cols, path, header)
        scale, unit = self._unit_scale(header, data[:, cmap["error"]])

        order = np.argsort(data[:, cmap["mjd"]])
        data = data[order]
        toas = data[:, cmap["mjd"]]
        residuals = data[:, cmap["residual"]] * scale
        errors = data[:, cmap["error"]] * scale
        freqs = data[:, cmap["frequency"]] if "frequency" in cmap else np.full(len(toas), np.nan)

        par_path = self._find_par(psr_name)
        par = read_par_file(par_path) if par_path else {}
        if not par:
            warnings.warn(f"No .par file for {psr_name}: sky position / DM unknown")

        span_days = toas[-1] - toas[0]
        return PulsarTimingData(
            name=psr_name, toas=toas, residuals=residuals, errors=errors, frequencies=freqs,
            ra=par.get("ra", np.nan), dec=par.get("dec", np.nan),
            parallax=par.get("parallax", np.nan), dm=par.get("dm", np.nan),
            period=par.get("period", np.nan), period_derivative=par.get("period_derivative", np.nan),
            n_observations=len(toas), time_span=span_days / 365.25,
            cadence=float(np.median(np.diff(np.unique(np.round(toas))))) if len(toas) > 1 else np.nan,
        )


class PulsarTimingEnsemble:
    """
    Ensemble of pulsars for common residual searches.
    """
    
    def __init__(self, pulsars: Dict[str, PulsarTimingData]):
        """
        Initialize ensemble.
        
        Args:
            pulsars: Dictionary mapping PSR names to PulsarTimingData
        """
        self.pulsars = pulsars
        self.n_pulsars = len(pulsars)
        
        # Find common time range
        all_toas = np.concatenate([p.toas for p in pulsars.values()])
        self.toa_min = np.min(all_toas)
        self.toa_max = np.max(all_toas)
    
    def pulsar_names(self) -> List[str]:
        """Pulsar names in the SAME order as the rows of common_residual_matrix()."""
        return sorted(self.pulsars.keys())

    def common_observing_epochs(self, n_bins: int = 100) -> np.ndarray:
        """
        Find common observing epochs across all pulsars.
        
        Args:
            n_bins: Number of time bins
        
        Returns:
            Array of common observation times (MJD)
        """
        epochs = np.linspace(self.toa_min, self.toa_max, n_bins)
        return epochs
    
    def align_to_epochs(self, epochs: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Interpolate all pulsars' residuals to common epochs.
        
        Args:
            epochs: Common observation times
        
        Returns:
            Dictionary mapping PSR name → residuals at common epochs
        """
        aligned = {}
        for name, pulsar_data in self.pulsars.items():
            # Linear interpolation
            residuals_interp = np.interp(
                epochs,
                pulsar_data.toas,
                pulsar_data.residuals,
                left=np.nan,
                right=np.nan
            )
            aligned[name] = residuals_interp
        
        return aligned
    
    def common_residual_matrix(self, epochs: Optional[np.ndarray] = None, n_bins: int = 100) -> np.ndarray:
        """
        Construct matrix: rows = pulsars, columns = time bins.
        
        This is the input for PCA, cross-correlation, etc.
        
        Args:
            epochs: Common epochs (default: auto-generate)
            n_bins: Number of time bins if epochs is None (default: 100)
        
        Returns:
            Matrix of shape (n_pulsars, n_epochs)
        """
        if epochs is None:
            epochs = self.common_observing_epochs(n_bins=n_bins)
        
        aligned = self.align_to_epochs(epochs)
        
        # Stack into matrix
        residual_matrix = np.array([aligned[name] for name in self.pulsar_names()])
        
        return residual_matrix


# ============================================================================
# Example Usage
# ============================================================================

if __name__ == "__main__":
    print("NANOGrav DR15 Loader Demo\n")
    
    import sys
    if "--real" in sys.argv:
        loader = NANOGravDR15Loader(use_synthetic=False)
        names = loader.available_pulsars()
        print(f"Real NANOGrav 15-yr pulsars found: {len(names)}")
        for n in names:
            try:
                d = loader.load_pulsar(n)
                print(f"  {n}: {d.n_observations:5d} rows, {d.time_span:5.1f} yr, "
                      f"RMS {d.residual_rms()*1e6:7.3f} µs, RA {d.ra:7.2f}, Dec {d.dec:+6.2f}")
            except Exception as e:
                print(f"  {n}: FAILED — {e}")
        sys.exit(0)

    # Load 5 example pulsars
    loader = NANOGravDR15Loader(use_synthetic=True, seed=1)
    
    psr_names = ['J0023+0923', 'J0030+0451', 'J0340+4130', 'J0437−4715', 'J0613−0200']
    pulsars = {}
    
    for psr_name in psr_names:
        data = loader.load_pulsar(psr_name)
        pulsars[psr_name] = data
        
        print(f"{psr_name}:")
        print(f"  Observations: {data.n_observations}")
        print(f"  Time span: {data.time_span:.1f} years")
        print(f"  Residual RMS: {data.residual_rms()*1e6:.3f} µs")
        print(f"  Weighted RMS: {data.weighted_rms()*1e6:.3f} µs")
        print()
    
    # Create ensemble
    ensemble = PulsarTimingEnsemble(pulsars)
    print(f"Ensemble: {ensemble.n_pulsars} pulsars")
    print(f"Time span: {ensemble.toa_min:.0f} – {ensemble.toa_max:.0f} MJD")
    
    # Construct common residual matrix
    matrix = ensemble.common_residual_matrix(n_bins=100)
    print(f"Residual matrix shape: {matrix.shape}")
    print(f"Matrix RMS: {np.nanstd(matrix)*1e9:.2f} ns")
