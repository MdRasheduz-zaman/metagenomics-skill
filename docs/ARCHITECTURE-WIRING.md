# metagx wiring — the mental map (and how the validation reference stays in scope)

This is the visual companion to `metagx/wiring.py` (run `metagx wiring`). It shows **every moving
part a tool or module touches**, so adding one and forgetting another is caught — by the audit, not
by memory. Part 2 answers a specific question: *how does blastn validate against the same references
as kraken2, and what does kraken2 actually expose?*

---

## Part 1 — The wiring DAG

A tool/module is wired across many independently-defined parts. The **single sources of truth** are
the per-tool registries (`parameters/*.yaml`), `DEFAULT_MODULES`, and `tool_advisor.MODULE_TOOLS`;
everything else consumes them. `metagx wiring` cross-checks each edge and fails on any gap.

```mermaid
flowchart TD
    subgraph SOT["Sources of truth"]
        REG["parameters/{tool}.yaml<br/>(registry: flags, managed, version_probe)"]
        MODS["config_builder.DEFAULT_MODULES<br/>(the module toggles)"]
        MT["tool_advisor.MODULE_TOOLS<br/>(canonical module → tools)"]
    end

    subgraph CFG["Config layer"]
        CLI["config_builder.build_config<br/>(CLI tool sections)"]
        MCP["mcp_server.build_config<br/>(MCP tool sections)"]
        DBK["config_builder.DB_EXTRA_KEYS<br/>(accepted db.{key})"]
    end

    subgraph RUN["Workflow + DB layer"]
        SNK["workflow/Snakefile<br/>(include guard + targets)"]
        RULE["workflow/rules/{module}.smk<br/>(render_args → command)"]
        SPECS["dbprovision.SPECS<br/>(module-DB fetchers)"]
        DOC["doctor.needed_dbs / checks<br/>(fail-fast guards)"]
    end

    subgraph OUT["Provenance + surfaces"]
        ACT["report.active_tools<br/>(version capture)"]
        CIT["report.CITATIONS<br/>(methods/paper)"]
        IV["prompts/INTERVIEW.md"]
        SK["SKILL.md"]
    end

    REG -->|"A: user tool ⇒ CLI section"| CLI
    REG -->|"B: user tool ⇒ MCP section"| MCP
    REG --> RULE
    SPECS -->|"C: SPEC ⊆ DB_EXTRA_KEYS"| DBK
    SPECS --> DOC
    MODS -->|"D: module ⇒ MODULE_TOOLS"| MT
    MODS -->|"E: module ⇒ include guard"| SNK
    MODS -->|"F: module ⇒ docs"| IV
    MODS -->|"F: module ⇒ docs"| SK
    MT -->|"G: tools ⇒ version capture"| ACT
    ACT -->|"H: tool ⇒ citation"| CIT
    SNK --> RULE
    DBK --> DOC

    classDef sot fill:#1f6feb,stroke:#0d419d,color:#fff;
    classDef chk fill:#f6f8fa,stroke:#8b949e,color:#24292f;
    class REG,MODS,MT sot;
```

**The labelled edges are the audit invariants** (A–H in `wiring.audit()`):

| Edge | Invariant | What a gap means |
|------|-----------|------------------|
| **A** | every registry *user* tool is a `build_config` kwarg (CLI) | a tool you can't actually configure |
| **B** | …and a `mcp_server.build_config` kwarg (MCP) | CLI and web/MCP surfaces drifted |
| **C** | every `dbprovision.SPECS` key ∈ `config_builder.DB_EXTRA_KEYS` | a provisioner whose `db.<key>` path the config rejects |
| **D** | every module ∈ `tool_advisor.MODULE_TOOLS` | `recommend`/`advise` blind to the module's tools |
| **E** | every module referenced in `workflow/Snakefile` | a toggle that includes no rule |
| **F** | every module documented in INTERVIEW.md **and** SKILL.md | the interview can't ask for it |
| **G** | each enabled module's tools appear in `report.active_tools` | tool version missing from the provenance manifest |
| **H** | every captured tool has a `report.CITATIONS` entry | methods/paper can't cite it |

> `kraken2-build`/`bracken-build` are DB-construction tools, not user sections (excluded from A/B).
> Routing modules (qc/assembly/phylogenetics/…) resolve tools per-platform, so G is checked over a
> kitchen-sink config. The MCP edge (B) is parsed from source with `ast`, so the audit runs even
> without the optional `mcp` extra installed.

### Adding a TOOL — the checklist the DAG enforces

```mermaid
flowchart LR
    A["parameters/{tool}.yaml<br/>(scaffold → curate)"] --> B["build_config kwarg<br/>(CLI + MCP)"]
    B --> C{"needs a DB?"}
    C -->|yes| D["dbprovision.SPECS<br/>+ DB_EXTRA_KEYS<br/>+ needed_dbs"]
    C -->|no| E
    D --> E["report.active_tools<br/>+ CITATIONS"]
    E --> F["metagx wiring → 0 gaps"]
```

### Adding a MODULE — the checklist the DAG enforces

```mermaid
flowchart LR
    A["DEFAULT_MODULES toggle"] --> B["workflow/rules/{m}.smk<br/>+ Snakefile include + targets"]
    B --> C["tool_advisor.MODULE_TOOLS"]
    C --> D["config block validator<br/>(config_builder)"]
    D --> E["report.active_tools<br/>+ doctor guard"]
    E --> F["INTERVIEW.md + SKILL.md"]
    F --> G["metagx wiring → 0 gaps"]
```

---

## Part 2 — How blastn validates against the *same references* as kraken2

**Short answer:** they don't share a database file — they share the **source genomes**. kraken2's
built DB is an opaque k-mer hash you cannot BLAST. But kraken2-build leaves the genomes it ingested
on disk as plain FASTA, and we run `makeblastdb` on **that same FASTA**. Same input genomes ⇒ both
tools cover the same organism set ⇒ the BLAST cross-check is a fair, in-scope benchmark.

### What kraken2 actually writes (verified on `local_databases/viral_custom`)

```
viral_custom/
├── hash.k2d              ← OPAQUE binary: minimizer → LCA taxon hash  (you CANNOT BLAST this)
├── opts.k2d  taxo.k2d    ← build options + taxonomy tree (binary)
├── custom_library.fasta  ← the GENOMES kraken2 ingested  (>acc|kraken:taxid|N …)  ← shareable
├── library/added/*.fna   ← (same sequences, as added)                              ← shareable
├── seqid2taxid.map       ← accession → taxid
└── taxonomy/{names,nodes}.dmp
```

`file hash.k2d` → `data` (binary). It stores *k-mers → taxon*, not retrievable sequences — so blastn
has nothing to align against in `.k2d`. The **sequences** live only in `custom_library.fasta` /
`library/`. That FASTA is exactly what `makeblastdb` needs.

```mermaid
flowchart TD
    SRC["Source genomes (FASTA)<br/>e.g. genomes.fasta / db.build source"]
    SRC -->|"kraken2-build --add-to-library → --build"| K2LIB["custom_library.fasta + library/<br/>(kept on disk)"]
    K2LIB -->|"kraken2-build hashes k-mers"| K2["hash.k2d / taxo.k2d<br/>(opaque — classify only)"]
    K2 -->|"kraken2 classify reads"| CALLS["per-read taxon calls (.kraken)"]

    K2LIB -->|"makeblastdb -dbtype nucl -parse_seqids<br/>(metagx build-blast-db / validate.build_from)"| BDB["BLAST nucl DB<br/>(*.nin/.nhr/.nsq)"]
    CALLS -->|"top taxa → read subsample"| Q["query reads"]
    Q -->|"blastn vs SAME genomes"| AGREE["per-taxon agreement + verdict<br/>(in scope ✓)"]
    BDB --> AGREE

    NT["full NCBI nt (~200 GB)"] -.->|"different organism set"| WARN["⚠ different benchmark<br/>(false 'disagreements')"]

    classDef opaque fill:#6e7681,stroke:#30363d,color:#fff;
    classDef share fill:#1a7f37,stroke:#116329,color:#fff;
    class K2 opaque;
    class K2LIB,SRC,BDB share;
```

### Empirical validation (what I ran — you can re-run it)

```bash
# kraken2's own ingested library:
grep -c '^>' local_databases/viral_custom/custom_library.fasta          # → 30 genomes

# build the BLAST validation DB from that SAME file:
metagx build-blast-db --from local_databases/viral_custom/custom_library.fasta --out /tmp/insync
#   → makeblastdb "added 30 sequences"

# both DBs hold the same 30 accessions:
grep '^>' local_databases/viral_custom/custom_library.fasta | sed 's/|.*//;s/>//' | sort      # 30
blastdbcmd -db /tmp/insync -entry all -outfmt '%a' | sed 's/\..*//' | sort -u                 # 30  (identical set)
```

Result: **30 ⇄ 30, same accessions** — the BLAST reference is exactly the classifier's organism set.

### How metagx wires this for you

- `validate.build_from: <FASTA|folder>` — the genomes you used for the classifier; the
  `build_validate_blast_db` rule runs `makeblastdb` on them before validating (no separate `db.blast`).
- `validate.build_from: classifier` — reuse the `db.build` **source** automatically, for
  `custom-fasta` / `custom-folder` / `spike-in` builds (where a local source FASTA exists).
- `db.blast: <path>` / `validate.remote: true` — only when you *deliberately* want a broader
  reference (e.g. nt). That is a different benchmark, by design.

```mermaid
flowchart LR
    Q{"how is the classifier DB built?"}
    Q -->|"db.build custom-fasta/folder/spike-in"| A["validate.build_from: classifier<br/>(auto-reuses db.build source)"]
    Q -->|"your own genomes FASTA"| B["validate.build_from: refs.fasta<br/>(or metagx build-blast-db)"]
    Q -->|"standard / prebuilt index (fetch-db)"| C["point build_from at the genome FASTAs you used<br/>⏳ auto-derive from index library/ = future pass"]
    Q -->|"want a broader benchmark on purpose"| D["db.blast: nt  /  validate.remote: true"]
```

> **Limit (next pass):** for a *standard* or *prebuilt-downloaded* kraken2 index, metagx does not yet
> auto-derive the BLAST DB from the index's `library/*/library.fna`. Until then, set `build_from` to
> the genome FASTA(s) you classified against. Tracked in the ROADMAP.

---

*Keep this in sync with `metagx/wiring.py`. If you change an invariant there, update Part 1 here;*
*`metagx wiring` is the executable version of this picture.*
