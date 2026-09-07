"""Local FastAPI application for MagneticTB."""

from __future__ import annotations


def main() -> None:
    """Run the single-process local web application."""

    try:
        import uvicorn
    except ImportError as error:  # pragma: no cover - exercised without the extra
        raise SystemExit(
            "Web dependencies are missing; reinstall the local MagneticTB wheel "
            "with pip without --no-deps"
        ) from error
    uvicorn.run("magnetictb.web.app:app", host="127.0.0.1", port=8000, workers=1)


__all__ = ["main"]
