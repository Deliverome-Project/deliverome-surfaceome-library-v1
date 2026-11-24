# Digest-pinned offline replay image for the L1 surface-protein ORF library review repository.
#
# Both base images are pinned by immutable OCI manifest digest (verified 2026-11-22):
#   - python:3.14.10-slim-bookworm (matches .python-version)
#   - ghcr.io/astral-sh/uv:0.13.6  (matches >=0.13.4 in pyproject.toml)
#
# Build once while PyPI is reachable, then save the image tarball alongside the frozen S3 bundle:
#   docker build -t deliverome-surfaceome-library-v1:frozen .
#   docker save deliverome-surfaceome-library-v1:frozen | gzip -n > deliverome-surfaceome-library-v1-image.tar.gz
#
# Replay years later with zero network access:
#   docker load -i deliverome-surfaceome-library-v1-image.tar.gz
#   docker run --rm --network=none \
#     -v /path/to/l1-reviewer-bundle:/bundle:ro \
#     deliverome-surfaceome-library-v1:frozen

FROM python:3.14.10-slim-bookworm@sha256:3fe1a0b6c2d6b2ebb572906cdabeca75aeb558c0c17db5ffefbcf8094871c32b

COPY --from=ghcr.io/astral-sh/uv:0.13.6@sha256:af22a52f5eafef3d78e551599d405c02df61919d18028f0cea8eda93bc9986d0 /uv /uvx /usr/local/bin/

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=0 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=never

WORKDIR /workspace

# Install locked dependencies first for layer caching.
# Flow/modal groups are excluded: L1 does not use agentflow-cytometry or Modal at replay time.
COPY pyproject.toml uv.lock .python-version ./
COPY src/ ./src/
RUN uv sync --locked --no-install-project --group dev --no-group flow --no-group flow-gui --no-group modal

COPY . .
RUN uv sync --locked --group dev --no-group flow --no-group flow-gui --no-group modal

ENTRYPOINT ["uv", "run", "--locked", "--offline", "python", "library/l1/orf-design/reviewer-package/run_review.py"]
CMD ["--bundle", "/bundle"]
