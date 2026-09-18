//
// Enrich, gene set descriptions and Export for AI
// MIT Open source
// -
// MSigDB's one-line description of every gene set (web_data/msigdb_descriptions.json,
// built from the MSigDB release JSON) is loaded once, on first use, and shown
// wherever a set is named. Export for AI writes a JSON file with the loaded
// data, the settings and the results, plus instructions for a language model,
// in the same spirit as Correlate's Export for AI.
//

const DESC = { map: null, loading: null };

function DESC_load() {
    if (DESC.map) return Promise.resolve(DESC.map);
    if (!DESC.loading) {
        DESC.loading = fetch('web_data/msigdb_descriptions.json?v=2023.2')
            .then(r => r.ok ? r.json() : {})
            .then(m => { DESC.map = m; return m; })
            .catch(() => { DESC.map = {}; return DESC.map; });
    }
    return DESC.loading;
}

// { d: brief description, p: PubMed id, s: exact source } or null
function DESC_get(name) {
    return DESC.map ? (DESC.map[name] || null) : null;
}

Object.assign(GSEAApp.prototype, {

    // Plain-text explanation of a gene set for a title attribute
    describeSet(name) {
        const e = DESC_get(name);
        const coll = this._getSetCollection(name);
        if (!e) return `${name} (${coll}). Description not loaded yet.`;
        return `${name} (${coll}): ${e.d}${e.p ? ' PubMed ' + e.p + '.' : ''}`;
    },

    // ---------------- Export for AI ----------------

    openAIExportDialog() {
        if (!this.results) { alert('Run an analysis first.'); return; }
        const dlg = document.getElementById('aiAnalysisDialog');
        document.getElementById('aiQuestion').value = '';
        document.getElementById('aiExportStatus').textContent = '';
        const nSig = this.results.filter(r => r.fdr < 0.25).length;
        document.getElementById('aiDialogSource').textContent = `${this._dataDescription()} ${this.results.length.toLocaleString()} gene sets tested, ${nSig} with FDR below 0.25.`;
        document.getElementById('aiDataTierInfo').innerHTML =
            `The file will contain: what data was analysed and how (source, ranking metric, settings, collections); the 25 top and bottom genes of the ranked list; every tested gene set with its collection, one-line description, size, ES, NES, p-value, FDR and leading-edge genes (all significant sets in full, the rest in brief); and instructions that tell the model how to read it.`;
        dlg.style.display = 'flex';
        document.getElementById('aiAnalysisDialogClose').onclick = () => { dlg.style.display = 'none'; };
        dlg.onclick = (e) => { if (e.target === dlg) dlg.style.display = 'none'; };
        document.getElementById('aiExportBtn').onclick = () => this.exportForAI();
    },

    _dataDescription() {
        const st = document.getElementById('uploadStatus');
        const txt = st && !st.classList.contains('hidden') ? st.textContent.trim() : '';
        if (txt) return txt.replace(/\s+/g, ' ');
        const f = document.getElementById('fileInput');
        return f && f.files && f.files[0] ? `Uploaded file ${f.files[0].name}.` : 'Uploaded data.';
    },

    async exportForAI() {
        const status = document.getElementById('aiExportStatus');
        status.textContent = 'Building...';
        try { await DESC_load(); } catch (e) { /* descriptions are optional */ }
        const question = document.getElementById('aiQuestion').value.trim() || null;
        const dt = this.settings.dataType || document.getElementById('dataType')?.value || 'expression';
        const metric = document.getElementById('metricColumn').value || 'metric';
        const genes = this.rankedList.genes, metrics = this.rankedList.metrics;
        const N = genes.length;
        const topGenes = genes.slice(0, 25).map((g, i) => ({ gene: g, value: +metrics[i].toFixed(4) }));
        const bottomGenes = genes.slice(-25).map((g, i) => ({ gene: g, value: +metrics[N - 25 + i].toFixed(4) })).reverse();
        const collections = [];
        ['checkHallmark', 'checkC2kegg', 'checkC2reactome', 'checkC2wp', 'checkC2biocarta', 'checkC2pid', 'checkC2cgp', 'checkC3', 'checkC5bp', 'checkC5cc', 'checkC5mf', 'checkC5hpo', 'checkC6', 'checkC7', 'checkC8'].forEach(id => {
            const cb = document.getElementById(id);
            if (cb && cb.checked) { const lab = cb.closest('label')?.querySelector('strong'); collections.push(lab ? lab.textContent.trim() : id.replace('check', '')); }
        });
        if (this.customGeneSets && Object.keys(this.customGeneSets).length) collections.push(`custom GMT (${Object.keys(this.customGeneSets).length} sets)`);
        if (this.useCustomSelection && this.selectedGeneSets && this.selectedGeneSets.size) collections.push(`custom selection of ${this.selectedGeneSets.size} sets`);
        const sorted = this.results.slice().sort((a, b) => a.fdr - b.fdr || Math.abs(b.nes) - Math.abs(a.nes));
        const sets = sorted.map(r => {
            const d = DESC_get(r.name);
            const sig = r.fdr < 0.25;
            const row = {
                name: r.name, collection: r.collection || this._getSetCollection(r.name),
                description: d ? d.d : null, pubmed: d && d.p ? d.p : null,
                size: r.size, ES: +r.es.toFixed(4), NES: +r.nes.toFixed(3), pValue: +r.pvalue.toExponential(3), FDR: +r.fdr.toExponential(3)
            };
            if (sig) row.leadingEdge = r.leadingEdge || [];
            else row.leadingEdgeCount = (r.leadingEdge || []).length;
            return row;
        });
        const hidden = this._hiddenSets ? [...this._hiddenSets] : [];
        const payload = {
            tool: 'Enrich, a gene set enrichment analysis (GSEA) tool from the Wermeling Lab, Karolinska Institutet',
            version: (document.getElementById('versionBadge')?.textContent || '').trim(),
            exportedAt: new Date().toISOString(),
            question,
            aiInstructions:
                'WHO YOU ARE TALKING TO: a biomedical researcher who ran gene set enrichment analysis (GSEA) in a tool called Enrich and attached this file. They have not read the file and do not know the names of its fields. Never mention a field name or key in your reply; say what it means in ordinary words.\n\n'
                + 'HOW TO READ IT: `context` says what was analysed. The ranking metric is what genes were sorted by: for expression data a positive value means higher in the sample or group A; for CRISPR screen data (Chronos score) a NEGATIVE value means the cells depend on the gene (knocking it out is harmful), so gene sets with negative NES are the ones the cells need. `results.sets` lists the gene sets, most significant first, with an FDR (false discovery rate) and a normalised enrichment score (NES). FDR below 0.25 is the usual GSEA threshold, below 0.05 is strong. A positive NES means the set sits at the top of the ranking, negative at the bottom. `leadingEdge` are the genes that drive that set\'s signal. Each set has a one-line description from MSigDB; use it to explain what a set is, since names like RB_P107_DN.V1_UP are not self-explanatory.\n\n'
                + 'OPEN YOUR REPLY with one short paragraph, no heading, saying in plain words what was analysed, how many sets were tested and how many were significant, and quoting their question back (or saying that no question came with the file).\n\n'
                + 'THEN answer their question if there is one, in a paragraph a colleague would understand, before any numbers. If there is no question, give one headline finding and one caveat, then ask what they want to look at. Group related gene sets together (many MSigDB sets overlap; say when several significant sets describe the same biology) and name the leading-edge genes that recur. Numbers support the answer, they are not the answer: give the few that matter, rounded, with what they mean. When the data is a single cell line compared with the DepMap median, remember that the result describes what makes this cell line different from an average cancer cell line, not a treatment effect. Do not invent gene sets or genes that are not in the file.',
            context: {
                description: this._dataDescription(),
                dataType: dt === 'crispr' ? 'CRISPR knockout screen (Chronos gene effect; negative = essential)' : 'gene expression (higher metric = higher expression in the sample or group A)',
                rankingMetric: metric,
                genesRanked: N,
                settings: {
                    permutations: this.settings.permutations, minSetSize: this.settings.minSize, maxSetSize: this.settings.maxSize,
                    weightP: this.settings.weightP, collections
                },
                comparison: this._loadedCellLine && this._loadedCellLine.compare ? this._loadedCellLine : null,
                topGenes, bottomGenes
            },
            results: {
                setsTested: this.results.length,
                setsSignificantFDR025: this.results.filter(r => r.fdr < 0.25).length,
                setsSignificantFDR005: this.results.filter(r => r.fdr < 0.05).length,
                hiddenByUser: hidden.length ? hidden : undefined,
                sets
            },
            whatIsInThisFile: ['context: the data and settings', 'context.topGenes and bottomGenes: the 25 genes at each end of the ranking', 'results.sets: every tested gene set, sorted by FDR, with description and leading-edge genes for significant sets'],
            notIncluded: ['the full ranked gene list (only the 25 genes at each end)', 'the gene membership of gene sets beyond the leading edge', 'any image of the plots']
        };
        const fmt = document.querySelector('input[name="aiExportFormat"]:checked')?.value || 'plain';
        const json = JSON.stringify(payload, null, 1);
        const stem = `enrich_for_ai_${exportStamp()}`;
        if (fmt === 'gz' && typeof CompressionStream !== 'undefined') {
            const gz = await new Response(new Blob([json]).stream().pipeThrough(new CompressionStream('gzip'))).blob();
            this._saveBlob(gz, `${stem}.json.gz`);
            status.textContent = `Saved ${stem}.json.gz (${Math.round(gz.size / 1024)} KB).`;
        } else {
            const blob = new Blob([json], { type: 'application/json' });
            this._saveBlob(blob, `${stem}.json`);
            status.textContent = `Saved ${stem}.json (${Math.round(blob.size / 1024)} KB).`;
        }
    }
});
