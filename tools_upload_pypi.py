# -*- coding: utf-8 -*-
"""
Build and upload synthbench to PyPI. Run this on your own machine, not in the sandbox:
upload.pypi.org is on the sandbox network denylist and cannot be reached from there.

Before running, set the distribution name in pyproject.toml. The name "synthbench" is
already taken on PyPI by an unrelated package (synthetic datasets for ML benchmarking),
so the `name = ` line must be changed; the import name stays `synthbench` either way.

    python -m pip install --upgrade build twine
    python tools_upload_pypi.py --check        # build and validate, upload nothing
    python tools_upload_pypi.py --test         # upload to test.pypi.org first
    python tools_upload_pypi.py --release      # upload to pypi.org

Credentials come from a PyPI API token, read from the TWINE_PASSWORD environment
variable or ~/.pypirc. Never paste a token into this file.
"""
import argparse
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).parent
DIST = HERE / "dist"


def run(*cmd):
    print("$", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=HERE)
    if r.returncode:
        sys.exit(r.returncode)


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true", help="build and validate only")
    g.add_argument("--test", action="store_true", help="upload to test.pypi.org")
    g.add_argument("--release", action="store_true", help="upload to pypi.org")
    a = ap.parse_args()

    name = next((l.split("=", 1)[1].strip().strip('"')
                 for l in (HERE / "pyproject.toml").read_text(encoding="utf-8").splitlines()
                 if l.strip().startswith("name")), None)
    if name == "synthbench":
        sys.exit("pyproject.toml still says name = \"synthbench\", which is taken on PyPI. "
                 "Change it before uploading.")
    print(f"distribution name: {name}")

    shutil.rmtree(DIST, ignore_errors=True)
    run(sys.executable, "-m", "build")
    run(sys.executable, "-m", "twine", "check", "dist/*")
    if a.check:
        print("built and validated; nothing uploaded")
        return
    repo = "testpypi" if a.test else "pypi"
    run(sys.executable, "-m", "twine", "upload", "-r", repo, "dist/*")


if __name__ == "__main__":
    main()
