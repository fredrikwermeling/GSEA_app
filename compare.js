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
        Promise.all([DEPMAP_index(), DEPMAP_mutations().catch(() => null)]).then(() => { this._cmpFillFilters(); this._cmpRender(); this._cmpSearch(); }).catch(e => {
            document.getElementById('cmpStatus').textContent = 'Could not load the DepMap index: ' + e.message;
        });
        this._cmpRender();
    },

    // Lineage and disease menus from the index; disease follows the chosen lineage
    _cmpFillFilters() {
        const idx = DEPMAP.index; if (!idx) return;
        const lin = document.getElementById('cmpLineage'), dis = document.getElementById('cmpDisease');
        const count = (key, filter) => {
            const m = new Map();
            for (const c of idx.all) { if (filter && !filter(c)) continue; const k = c[key] || ''; if (k) m.set(k, (m.get(k) || 0) + 1); }
            return [...m.entries()].sort((a, b) => a[0].localeCompare(b[0]));
        };
        const curLin = lin.value;
        lin.innerHTML = '<option value="">Any tissue</option>' + count('lineage').map(([k, n]) => `<option value="${k}">${k} (${n})</option>`).join('');
        lin.value = curLin;
        const curDis = dis.value;
        dis.innerHTML = '<option value="">Any disease</option>' + count('disease', lin.value ? (c => c.lineage === lin.value) : null).map(([k, n]) => `<option value="${k}">${k} (${n})</option>`).join('');
        dis.value = [...dis.options].some(o => o.value === curDis) ? curDis : '';
        const sub = document.getElementById('cmpSubtype');
        const curSub = sub.value;
        sub.innerHTML = '<option value="">Any subtype</option>' + count('subtype', c => (!lin.value || c.lineage === lin.value) && (!dis.value || c.disease === dis.value)).map(([k, n]) => `<option value="${k}">${k} (${n})</option>`).join('');
        sub.value = [...sub.options].some(o => o.value === curSub) ? curSub : '';
        const mut = DEPMAP.mut;
        const note = document.getElementById('cmpMutNote');
        if (note) note.textContent = mut ? `Hotspot calls for ${Object.keys(mut.hotspot).length} cancer genes, damaging calls for ${Object.keys(mut.damaging).length} genes; ${mut.profiled.length} lines have mutation data.` : 'Mutation data not available.';
    },

    // The cell lines that pass the tissue, disease and mutation filters (and the text box)
    _cmpCandidates() {
        const idx = DEPMAP.index; if (!idx) return [];
        const type = this._cmp.type;
        const lin = document.getElementById('cmpLineage').value, dis = document.getElementById('cmpDisease').value, sub = document.getElementById('cmpSubtype').value;
        const gene = document.getElementById('cmpMutGene').value.trim().toUpperCase();
        const status = document.getElementById('cmpMutStatus').value;
        const q = document.getElementById('cmpSearch').value.trim();
        let rows = q ? DEPMAP_search(q, 5000).rows : idx.all.slice();
        rows = rows.filter(c => c.rows[type] !== undefined && (!lin || c.lineage === lin) && (!dis || c.disease === dis) && (!sub || c.subtype === sub));
        const mut = DEPMAP.mut;
        if (gene && status !== 'any' && mut) {
            const hs = mut.hotspotSets[gene] || new Set(), dm = mut.damagingSets[gene] || new Set();
            rows = rows.filter(c => {
                if (!mut.profiledSet.has(c.id)) return false;           // unknown mutation status: never counted as either
                const isHot = hs.has(c.id), isDam = dm.has(c.id);
                if (status === 'hotspot') return isHot;
                if (status === 'mutated') return isHot || isDam;
                if (status === 'wt') return !isHot && !isDam;
                return true;
            });
        }
        return rows;
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
        const parts = [];
        const lin = document.getElementById('cmpLineage').value, dis = document.getElementById('cmpDisease').value, sub = document.getElementById('cmpSubtype').value;
        const gene = document.getElementById('cmpMutGene').value.trim().toUpperCase(), st = document.getElementById('cmpMutStatus').value;
        const q = document.getElementById('cmpSearch').value.trim();
        if (q) parts.push(`"${q}"`); if (lin) parts.push(lin); if (dis) parts.push(dis); if (sub) parts.push(sub);
        if (gene && st !== 'any') parts.push(`${gene} ${st === 'wt' ? 'wild type' : st === 'hotspot' ? 'hotspot mutated' : 'mutated'}`);
        if (!this._cmp.labels) this._cmp.labels = {};
        this._cmp.labels[group] = parts.join(', ');
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

    // Typing a gene is meant as a filter, so the status menu follows: hotspot
    // if the gene has hotspot calls, otherwise any mutation. The note under the
    // box says what the data holds for that gene, or that it holds nothing.
    _cmpGeneChanged() {
        const gene = document.getElementById('cmpMutGene').value.trim().toUpperCase();
        const st = document.getElementById('cmpMutStatus');
        const note = document.getElementById('cmpMutGeneNote');
        const mut = DEPMAP.mut;
        if (!gene) { if (note) note.textContent = ''; st.value = 'any'; this._cmpSearch(); return; }
        if (!mut) { if (note) note.textContent = 'Mutation data is still loading...'; DEPMAP_mutations().then(() => this._cmpGeneChanged()).catch(() => { if (note) note.textContent = 'Mutation data could not be loaded.'; }); return; }
        const nHot = (mut.hotspot[gene] || []).length, nDam = (mut.damaging[gene] || []).length;
        if (!nHot && !nDam) {
            if (note) note.innerHTML = `<span style="color:#b45309;">No mutation calls for ${this._escText(gene)} in DepMap. Check the symbol; hotspot calls exist for ${Object.keys(mut.hotspot).length} cancer genes, damaging calls for most genes.</span>`;
            st.value = 'any';
        } else {
            if (note) note.textContent = `${gene}: ${nHot} lines with a hotspot mutation, ${nDam} with a damaging mutation, ${mut.profiled.length} lines have mutation data.`;
            if (st.value === 'any') st.value = nHot ? 'hotspot' : 'mutated';
        }
        this._cmpSearch();
    },

    _cmpSearch() {
        const q = document.getElementById('cmpSearch').value.trim();
        const out = document.getElementById('cmpResults');
        if (!DEPMAP.index) { out.innerHTML = ''; return; }
        const anyFilter = q || document.getElementById('cmpLineage').value || document.getElementById('cmpDisease').value || document.getElementById('cmpSubtype').value
            || (document.getElementById('cmpMutGene').value.trim() && document.getElementById('cmpMutStatus').value !== 'any');
        if (!anyFilter) { out.innerHTML = '<div style="color: var(--gray-400); padding: 2px 4px; font-size: 0.85em;">Type a name, or choose a tissue, disease or mutation, to list cell lines.</div>'; this._cmpLastResults = []; return; }
        const type = this._cmp.type;
        const rows = this._cmpCandidates();
        this._cmpLastResults = rows;
        if (!rows.length) { out.innerHTML = `<div style="color: var(--gray-500); padding: 2px 4px;">No cell line with ${type} data matches these filters.</div>`; return; }
        let html = '';
        if (rows.length > 1) html += `<div style="padding: 2px 4px; font-size: 0.85em; color: var(--gray-600);">${rows.length} matching cell lines: <a href="#" onclick="app._cmpAddAll('A'); return false;">add all to A</a> &middot; <a href="#" onclick="app._cmpAddAll('B'); return false;">add all to B</a></div>`;
        for (const l of rows.slice(0, 60)) {
            const where = [l.lineage, l.disease, l.subtype && l.subtype !== l.disease ? l.subtype : ''].filter(Boolean).join(', ');
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
            const list = (ids, g) => {
                const lab = this._cmp.labels && this._cmp.labels[g];
                if (ids.length > 4 && lab) return `${lab} (${ids.length} lines: ${ids.slice(0, 3).map(nameOf).join(', ')}, ...)`;
                return ids.length <= 4 ? ids.map(nameOf).join(', ') : `${ids.slice(0, 3).map(nameOf).join(', ')} and ${ids.length - 3} more`;
            };
            const how = metricKind === 't'
                ? `Welch t-statistic per gene (positive = higher in A)`
                : (type === 'crispr' ? 'difference of mean Chronos score, A minus B (negative = more essential in A)' : 'log2 fold change of A over B (mean of A minus mean of B)');
            this.showStatus('uploadStatus', 'success', `Comparison built: A = ${list(c.A, 'A')} vs B = ${list(c.B, 'B')}, ${type} data, ${how}. ${out.length.toLocaleString()} genes.`);
            this._loadedCellLine = { compare: true, type, A: c.A.slice(), B: c.B.slice() };
            this.checkReady();
            status.textContent = '';
            this.closeCompareDialog();
        } catch (e) {
            status.textContent = 'Failed: ' + e.message;
        }
    }
});
