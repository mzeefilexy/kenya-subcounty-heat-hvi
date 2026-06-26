"""Download the remaining 2025 ERA5 inputs for a focused Heat Index comparison."""
from calendar import monthrange
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent / "pydeps"))
import cdsapi

out = Path(__file__).resolve().parents[1] / "data" / "era5_sensitivity"
out.mkdir(exist_ok=True)
client = cdsapi.Client(quiet=True)

for month in range(2, 13):
    target = out / f"era5_kenya_hourly_2025_{month:02d}.nc"
    if target.exists():
        continue
    client.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": ["reanalysis"],
            "variable": ["2m_temperature", "2m_dewpoint_temperature"],
            "year": ["2025"],
            "month": [f"{month:02d}"],
            "day": [f"{day:02d}" for day in range(1, monthrange(2025, month)[1] + 1)],
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "area": [5, 34, -5, 42],
            "data_format": "netcdf",
            "download_format": "unarchived",
        },
        str(target),
    )
    print(f"complete {month:02d}", flush=True)
