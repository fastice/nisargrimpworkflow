# makeMaster — moved to mosaicworkflow

`makeMaster` now lives in the **mosaicworkflow** package, which owns the
`makeMaster` console script. See
[mosaicworkflow/Documents/makeMaster.md](../../mosaicworkflow/Documents/makeMaster.md).

Nothing changes for callers here. The command name, its flags
(`--projectDir`, `--inputPaths`, `--outputPath`, `--verbose`) and its output
are identical — verified byte-for-byte against the version this replaced, for
both the Antarctica40 + Antarctica80 invocation and the `newGreenlandProject`
defaults. `autoupdate.py` invokes it by name and needed no change, so the
`masterInputPaths` / `masterOutputPath` keys of `autoupdate.yaml` work exactly
as documented in [autoupdate.md](autoupdate.md).

The move generalised the tool so one program builds both the NISAR master and
the Sentinel-1 / TSX masters, the latter from a YAML config listing the source
regions. It replaced the `createMaster.ipynb` notebook at the same time.

Two NISAR-relevant notes now live in the new doc:

- **Offsets-only (`nophase`) products** are still excluded from the NISAR
  master, but now through an explicit `nonFileTokens` setting rather than a
  misparse. What that exclusion costs is quantified there — worth a look.
- **Bare-name vs `.vrt` resolution** and the **baseline/VRT ionosphere
  consistency check** are unchanged and are described in full there.
