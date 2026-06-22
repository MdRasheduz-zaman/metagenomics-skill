# Auto-provision the per-tool module DBs listed in config db.provision (genomad/checkv/
# checkm2/gtdbtk/bakta/amrfinderplus). One wildcard rule writes a sentinel under
# {OUT}/dbprovision/<tool>.done; the actual download is idempotent (dbprovision.provision skips
# when the DB is already present), so this is cheap to re-evaluate. Domain/functional rules
# depend on the sentinel via provision_ready(<tool>), so the DB is fetched before it's used.


rule provision_module_db:
    output:
        sentinel=f"{OUT}/dbprovision/{{tool}}.done",
    params:
        tool=lambda wc: wc.tool,
        db_dir=lambda wc: DB[wc.tool],
    wildcard_constraints:
        tool="|".join(DB_PROVISION) if DB_PROVISION else "NONE",
    script:
        "../scripts/provision_db.py"
