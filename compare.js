//
// Enrich, compare cell lines
// MIT Open source
// -
// Builds a ranked gene list from two sets of DepMap cell lines, A and B:
// the difference of their means per gene (expression: log2 fold change of A
// over B; CRISPR: difference in Chronos score, negative = more essential in
// A), or Welch's t-statistic when both groups have at least two lines. The
// list then goes through GSEA like an uploaded file.
//

Object.assign(GSEAApp.prototype, {

    openCompareDialog() {
        if (!this._cmp) this._cmp = { A: [], B: [], type: 'expression' };
        document.getElementById('comparePopup').classList.add('open');
        document.getElementById('howToUseBackdrop').classList.add('open');
        DEPMAP_index().then(() => this._cmpRender()).catch(e => {
            document.getElementById('cmpStatus').textContent = 'Could not load the DepMap index: ' + e.message;
        });
        this._cmpRender();
    },

    closeCompareDialog() {
        document.getElementById('comparePopup').classList.remove('open');
        document.getElementById('howToUseBackdrop').classList.remove('open');
    },

    _cmpRender() {
        const c = this._cmp;
        const idx = DEPMAP.index;
        for (const g of ['A', 'B']) {
            const el = document.getElementById('cmpGroup' + g);
            el.innerHTML = c[g].length
                ? c[g].map(id => {
                    const line = idx ? idx.all.find(x => x.id === id) : null;
                    const name = line ? line.name : id;
                    return `<span class="coll-chip" style="margin: 2px; cursor: pointer;" title="Remove ${this._escText(name)} from group ${g}" onclick="app._cmpRemove('${g}','${id}')">${this._escText(name)} &times;</span>`;
                }).join('')
                : '<span style="color: var(--gray-400); font-size: 0.85em;">No cell lines yet. Search below and press +A.</span>';
            const n = document.getElementById('cmpCount' + g);
            if (n) n.textContent = c[g].length ? `${c[g].length}` : '';
        }
        const both = c.A.length && c.B.length;
        const tOk = c.A.length >= 2 && c.B.length >= 2;
        document.getElementById('cmpBuild').disabled = !both;
        const tOpt = document.querySelector('#cmpMetric option[value="t"]');
        if (tOpt) tOpt.disabled = !tOk;
        if (!tOk && document.getElementById('cmpMetric').value === 't') document.getElementById('cmpMetric').value = 'diff';
        document.querySelectorAll('input[name="cmpType"]').forEach(r => { r.checked = r.value === c.type; });
    },

    _cmpAdd(group, id) {
        const c = this._cmp;
        const other = group === 'A' ? 'B' : 'A';
        if (!c[group].includes(id)) c[group].push(id);
        c[other] = c[other].filter(x => x !== id);
        this._cmpRender();
        this._cmpSearch();
    },

    _cmpAddAll(group) {
        const res = this._cmpLastResults || [];
        for (const line of res) if (line.rows[this._cmp.type] !== undefined) this._cmpAdd(group, line.id);
    },

    _cmpRemove(group, id) {
        this._cmp[group] = this._cmp[group].filter(x => x !== id);
        this._cmpRender();
        this._cmpSearch();
    },

    _cmpSetType(type) {
        this._cmp.type = type;
        // Lines without data of this type cannot stay in a group
        for (const g of ['A', 'B']) this._cmp[g] = this._cmp[g].filter(id => (DEPMAP.index?.all.find(x => x.id === id) || { rows: {} }).rows[type] !== undefined);
        this._cmpRender();
        this._cmpSearch();
    },

    _cmpSearch() {
        const q = document.getElementById('cmpSearch').value.trim();
        const out = document.getElementById('cmpResults');
        if (!q || !DEPMAP.index) { out.innerHTML = ''; return; }
        const res = DEPMAP_search(q, 200);
        const type = this._cmp.type;
        const rows = res.rows.filter(l => l.rows[type] !== undefined);
        this._cmpLastResults = rows;
        if (!rows.length) { out.innerHTML = `<div style="color: var(--gray-500); padding: 2px 4px;">No cell line with ${type} data matches "${this._escText(q)}".</div>`; return; }
        let html = '';
        if (rows.length > 1) html += `<div style="padding: 2px 4px; font-size: 0.85em; color: var(--gray-600);">${rows.length} matches${res.total > rows.length ? ` (${res.total - rows.length} more without ${type} data)` : ''}: <a href="#" onclick="app._cmpAddAll('A'); return false;">add all to A</a> &middot; <a href="#" onclick="app._cmpAddAll('B'); return false;">add all to B</a></div>`;
        for (const l of rows.slice(0, 60)) {
            const where = [l.lineage, l.disease].filter(Boolean).join(', ');
            const inA = this._cmp.A.includes(l.id), inB = this._cmp.B.includes(l.id);
            html += `<div class="example-row"><div class="example-cell-info"><span class="example-cell-name">${this._escText(l.name)}</span><span class="example-cancer-type">${this._escText(where)}</span></div><div class="example-buttons">`
                + `<button class="btn btn-outline btn-xs" ${inA ? 'style="background: var(--green-100);"' : ''} onclick="app._cmpAdd('A','${l.id}')" title="Put this cell line in group A">${inA ? '&#10003; ' : '+'}A</button>`
                + `<button class="btn btn-outline btn-xs" ${inB ? 'style="background: var(--green-100);"' : ''} onclick="app._cmpAdd('B','${l.id}')" title="Put this cell line in group B">${inB ? '&#10003; ' : '+'}B</button></div></div>`;
        }
        if (rows.length > 60) html += `<div style="padding: 2px 4px; font-size: 0.85em; color: var(--gray-500);">Showing 60 of ${rows.length}. Type more to narrow, or add all.</div>`;
        out.innerHTML = html;
    },

    async buildComparison() {
        const c = this._cmp;
        const type = c.type;
        const metricKind = document.getElementById('cmpMetric').value;
        const status = document.getElementById('cmpStatus');
        const idx = await DEPMAP_index();
        const t = idx[type];
        const rowOf = (id) => t.cellLines.findIndex(x => x.id === id);
        const nameOf = (id) => (idx.all.find(x => x.id === id) || { name: id }).name;
        const fetchGroup = async (ids, label) => {
            const rows = [];
            for (let i = 0; i < ids.length; i++) {
                status.textContent = `Loading group ${label}: ${i + 1} of ${ids.length}...`;
                rows.push(await DEPMAP_rowValues(type, rowOf(ids[i])));
            }
            return rows;
        };
        try {
            const A = await fetchGroup(c.A, 'A');
            const B = await fetchGroup(c.B, 'B');
            status.textContent = 'Computing...';
            const nG = t.genes.length;
            const stat = (rows, g) => {
                let n = 0, s = 0, ss = 0;
                for (const r of rows) { const v = r[g]; if (!isNaN(v)) { n++; s += v; ss += v * v; } }
                const m = n ? s / n : NaN;
                const varc = n > 1 ? Math.max(0, (ss - n * m * m) / (n - 1)) : NaN;
                return { n, m, v: varc };
            };
            const out = [];
            const col = metricKind === 't' ? (type === 'crispr' ? 'Chronos_t_A_vs_B' : 'Expression_t_A_vs_B') : (type === 'crispr' ? 'Chronos_A_minus_B' : 'log2FC_A_vs_B');
            for (let g = 0; g < nG; g++) {
                const a = stat(A, g), b = stat(B, g);
                if (!a.n || !b.n) continue;
                let val;
                if (metricKind === 't') {
                    if (a.n < 2 || b.n < 2) continue;
                    const se = Math.sqrt(a.v / a.n + b.v / b.n);
                    if (!(se > 0)) continue;
                    val = (a.m - b.m) / se;
                } else {
                    val = a.m - b.m;
                }
                out.push({ Gene: t.genes[g], [col]: Math.round(val * 10000) / 10000 });
            }
            out.sort((x, y) => y[col] - x[col]);
            this.rawData = out;
            this.populateColumnDropdowns(['Gene', col]);
            document.getElementById('geneColumn').value = 'Gene';
            document.getElementById('metricColumn').value = col;
            const dtVal = type === 'crispr' ? 'crispr' : 'expression';
            document.getElementById('dataType').value = dtVal;
            const dtInline = document.getElementById('dataTypeInline');
            if (dtInline) dtInline.value = dtVal;
            this.settings.dataType = dtVal;
            const h = document.getElementById('checkHallmark');
            if (h && !h.checked) { h.checked = true; await this.onCollectionChange(); }
            const list = (ids) => ids.length <= 4 ? ids.map(nameOf).join(', ') : `${ids.slice(0, 3).map(nameOf).join(', ')} and ${ids.length - 3} more`;
            const how = metricKind === 't'
                ? `Welch t-statistic per gene (positive = higher in A)`
                : (type === 'crispr' ? 'difference of mean Chronos score, A minus B (negative = more essential in A)' : 'log2 fold change of A over B (mean of A minus mean of B)');
            this.showStatus('uploadStatus', 'success', `Comparison built: A = ${list(c.A)} (${c.A.length}) vs B = ${list(c.B)} (${c.B.length}), ${type} data, ${how}. ${out.length.toLocaleString()} genes.`);
            this._loadedCellLine = { compare: true, type, A: c.A.slice(), B: c.B.slice() };
            this.checkReady();
            status.textContent = '';
            this.closeCompareDialog();
        } catch (e) {
            status.textContent = 'Failed: ' + e.message;
        }
    }
});
