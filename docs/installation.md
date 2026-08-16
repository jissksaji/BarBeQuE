# Installation

## Requirements

- Nextflow DSL2, version 24.10.5 or newer
- Java compatible with your Nextflow version
- One software provisioning backend:
  - Singularity
  - Apptainer
  - Docker
  - Podman
  - Conda

Container backends are recommended for production. Conda is useful for development and small local tests, but it is less reproducible across platforms and time.

## Install Nextflow

Follow the official Nextflow installation guide:

```text
https://www.nextflow.io/docs/latest/getstarted.html#installation
```

Check the installed version:

```bash
nextflow -version
```

## Choose A Profile

Common local profiles:

```bash
-profile singularity
-profile apptainer
-profile docker
-profile podman
-profile conda
```

For cluster or cloud execution, provide a site-specific Nextflow config with `-c` or use a shared profile from the bio-raum config repository.

## Install References

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --build_references \
  --reference_base /path/to/references
```

Do not include `barbeque/<reference_version>` in `--reference_base`; the pipeline adds that structure itself.

## Verify The Installation

List installed/known databases:

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --list_dbs \
  --reference_base /path/to/references
```

List named primer sets from the upstream catalog:

```bash
nextflow run bio-raum/BarBeQuE \
  -profile singularity \
  --list_primers
```

Run a small local test when the required software backend is available:

```bash
bash run_all_tests.sh
```

The module tests require the selected backend and external bioinformatics packages. The Python unit tests at the start of `run_all_tests.sh` only require `python3`.
