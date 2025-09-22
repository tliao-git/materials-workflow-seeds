# materials-workflow-seeds

Two minimal seeds for AI/DFT materials workflows:

- `chgnet_relax_cli/`: CHGNet + ASE pre-relaxation and optional short MD “shake test” with sanity checks and provenance. Seed for atomate2 integration to reduce wasted DFT jobs.
- `hull_tools/`: Lightweight pymatgen convex-hull tools from a CSV of compositions and energies; computes formation energies and distance-to-hull, optional plot.

See each subfolder for setup and usage.

## Authorship & tools
Developed by me. AI code-completion was used for boilerplate.
