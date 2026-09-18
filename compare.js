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
        // Counts are for the chosen data type only, so the menus say how many
        // lines can actually go into a group.
        const type = this._cmp.type;
        const count = (key, filter) => {
            const m = new Map();
            for (const c of idx.all) { if (c.rows[type] === undefined) continue; if (filter && !filter(c)) continue; const k = c[key] || ''; if (k) m.set(k, (m.get(k) || 0) + 1); }
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
        this._cmpFillFilters();
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
        const type = this._cmp.type, has = (id) => (DEPMAP.index?.all.find(x => x.id === id) || { rows: {} }).rows[type] !== undefined;
        const nHot = (mut.hotspot[gene] || []).filter(has).length, nDam = (mut.damaging[gene] || []).filter(has).length;
        if (!nHot && !nDam) {
            if (note) note.innerHTML = `<span style="color:#b45309;">No mutation calls for ${this._escText(gene)} in DepMap. Check the symbol; hotspot calls exist for ${Object.keys(mut.hotspot).length} cancer genes, damaging calls for most genes.</span>`;
            st.value = 'any';
        } else {
            if (note) note.textContent = `${gene}: ${nHot} lines with a hotspot mutation, ${nDam} with a damaging mutation, among the ${DEPMAP.index.all.filter(c => c.rows[type] !== undefined && mut.profiledSet.has(c.id)).length} ${type} lines with mutation data.`;
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
            this._cmpData = { type, genes: t.genes, geneIndex: new Map(t.genes.map((g, i) => [g, i])), mutGene: (document.getElementById('cmpMutGene').value || '').trim().toUpperCase(),
                A: { ids: c.A.slice(), names: c.A.map(nameOf), rows: A, label: (this._cmp.labels && this._cmp.labels.A) || '' },
                B: { ids: c.B.slice(), names: c.B.map(nameOf), rows: B, label: (this._cmp.labels && this._cmp.labels.B) || '' } };
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


// ---------------- Gene set heatmap: genes x cell lines of the two groups ----------------
Object.assign(GSEAApp.prototype, {

    hasComparisonData() { return !!(this._cmpData && this._cmpData.A.rows.length && this._cmpData.B.rows.length); },

    openGeneSetHeatmap(geneSetName) {
        if (!this.hasComparisonData()) { alert('The heatmap needs a comparison built from DepMap cell lines (Compare cell lines in the sidebar).'); return; }
        const sel = document.getElementById('geneSetSelector');
        const name = geneSetName || (sel && sel.value);
        if (!name) return;
        this._hmSet = name;
        document.getElementById('gsHeatmapPopup').classList.add('open');
        document.getElementById('howToUseBackdrop').classList.add('open');
        this.renderGeneSetHeatmap();
    },

    closeGeneSetHeatmap() {
        document.getElementById('gsHeatmapPopup').classList.remove('open');
        document.getElementById('howToUseBackdrop').classList.remove('open');
    },

    renderGeneSetHeatmap() {
        const d = this._cmpData, name = this._hmSet;
        const result = this.results && this.results.find(r => r.name === name);
        if (!d || !result) return;
        const leOnly = document.getElementById('hmLeadingOnly').checked;
        const zscore = document.getElementById('hmZscore').checked;
        const sortBy = document.getElementById('hmSort').value;
        const sortGene = (document.getElementById('hmSortGene').value || '').trim().toUpperCase();
        const stripGeneEl = document.getElementById('hmMutGene');
        if (!stripGeneEl.value && d.mutGene) stripGeneEl.value = d.mutGene;
        const stripGene = (stripGeneEl.value || '').trim().toUpperCase();
        const le = new Set((result.leadingEdge || []).map(g => g.toUpperCase()));
        let genes = (result.hits || []).map(i => this.rankedList.genes[i]).filter(g => d.geneIndex.has(g));
        if (leOnly) genes = genes.filter(g => le.has(g));
        const isCrispr = d.type === 'crispr';
        const idx = DEPMAP.index, mut = DEPMAP.mut;
        const lineOf = (id) => (idx && idx.all.find(x => x.id === id)) || {};
        const mutStatus = (id) => {
            if (!stripGene || !mut) return '';
            if (!mut.profiledSet.has(id)) return 'no data';
            const h = (mut.hotspotSets[stripGene] || new Set()).has(id), dm = (mut.damagingSets[stripGene] || new Set()).has(id);
            return h ? 'hotspot' : dm ? 'damaging' : 'wild type';
        };
        // one descriptor per column (cell line)
        let cols = [];
        for (const g of ['A', 'B']) d[g].ids.forEach((id, k) => {
            const l = lineOf(id);
            cols.push({ id, name: d[g].names[k], group: g, row: d[g].rows[k], lineage: l.lineage || '', disease: l.disease || '', subtype: l.subtype || '', mut: mutStatus(id), order: cols.length });
        });
        // values per gene, optionally z-scored across all shown cell lines
        const geneVals = new Map();
        for (const g of genes) {
            const gi = d.geneIndex.get(g);
            let vals = cols.map(c => c.row[gi]);
            if (zscore) {
                const ok = vals.filter(v => !isNaN(v)); const m = ok.reduce((a, b) => a + b, 0) / (ok.length || 1);
                const sd = Math.sqrt(ok.reduce((a, b) => a + (b - m) * (b - m), 0) / (ok.length > 1 ? ok.length - 1 : 1)) || 1;
                vals = vals.map(v => isNaN(v) ? null : (v - m) / sd);
            } else vals = vals.map(v => isNaN(v) ? null : v);
            geneVals.set(g, vals);
        }
        const mean = (arr) => { const ok = arr.filter(v => v !== null && v !== undefined); return ok.length ? ok.reduce((a, b) => a + b, 0) / ok.length : null; };
        // column score for sorting: one typed gene, else the mean of the
        // leading-edge genes (the signal GSEA found), else the mean of all shown genes
        const desc = document.getElementById('hmSortDir').value !== 'asc';
        const leGenes = genes.filter(g => le.has(g));
        const scoreGenes = sortGene && geneVals.has(sortGene) ? [sortGene] : (sortBy.startsWith('le') && leGenes.length ? leGenes : genes);
        cols.forEach((c, k) => { c.score = mean(scoreGenes.map(g => geneVals.get(g)[k])); });
        const byGroup = (a, b) => (a.group < b.group ? -1 : a.group > b.group ? 1 : 0);
        const num = (a, b) => desc ? ((b.score ?? -Infinity) - (a.score ?? -Infinity)) : ((a.score ?? Infinity) - (b.score ?? Infinity));
        const str = (key) => (a, b) => String(a[key]).localeCompare(String(b[key])) || byGroup(a, b) || a.order - b.order;
        const sorters = {
            group: (a, b) => byGroup(a, b) || a.order - b.order,
            leWithin: (a, b) => byGroup(a, b) || num(a, b) || a.order - b.order,
            le: (a, b) => num(a, b) || a.order - b.order,
            valueWithin: (a, b) => byGroup(a, b) || num(a, b) || a.order - b.order,
            value: (a, b) => num(a, b) || a.order - b.order,
            mutation: (a, b) => byGroup(a, b) || String(a.mut).localeCompare(String(b.mut)) || a.order - b.order,
            disease: str('disease'), subtype: str('subtype'), lineage: str('lineage')
        };
        cols.sort(sorters[sortBy] || sorters.group);
        const scoreLabel = sortGene && geneVals.has(sortGene) ? sortGene : (scoreGenes === leGenes ? 'leading-edge score' : 'mean of shown genes');
        const order = cols.map(c => c.order);
        const colNames = cols.map(c => c.name);
        const nA = cols.filter(c => c.group === 'A').length, nB = cols.length - nA;
        // rows: genes sorted by mean A minus mean B, leading edge starred
        const rows = genes.map(g => {
            const vals = order.map(k => geneVals.get(g)[k]);
            const mA = mean(cols.map((c, k) => c.group === 'A' ? vals[k] : null)), mB = mean(cols.map((c, k) => c.group === 'B' ? vals[k] : null));
            return { g, vals, mA, mB, diff: (mA ?? 0) - (mB ?? 0) };
        }).sort((a, b) => b.diff - a.diff);
        const metricName = zscore ? 'z-score per gene' : (isCrispr ? 'Chronos gene effect' : 'log2 expression vs DepMap median');
        const yLabels = rows.map(r => (le.has(r.g) ? '\u2605 ' : '') + r.g);
        const z = rows.map(r => r.vals);
        const absMax = Math.max(0.5, ...z.flat().filter(v => v !== null).map(Math.abs));
        const zmax = zscore ? Math.min(absMax, 3) : Math.min(absMax, isCrispr ? 2 : 6);
        const showX = cols.length <= 80;
        // annotation strips below the map: group, mutation, disease, subtype, lineage
        const strips = [{ key: 'group', label: 'Group', vals: cols.map(c => c.group) }];
        if (stripGene) strips.push({ key: 'mut', label: `${stripGene} mutation`, vals: cols.map(c => c.mut) });
        for (const [key, label] of [['disease', 'Disease'], ['subtype', 'Subtype'], ['lineage', 'Tissue']]) {
            const vals = cols.map(c => c[key]);
            if (new Set(vals).size > 1) strips.push({ key, label, vals });
        }
        const palette = ['#7ab950', '#4472c4', '#ed7d31', '#a5a5a5', '#ffc000', '#5b9bd5', '#70ad47', '#9e480e', '#636363', '#997300', '#264478', '#43682b', '#c9c9c9', '#f4b183', '#8faadc', '#c5e0b4'];
        const fixed = { A: '#dc2626', B: '#2563eb', hotspot: '#b91c1c', damaging: '#f59e0b', 'wild type': '#cbd5e1', 'no data': '#f3f4f6' };
        const stripZ = [], stripText = [], stripColors = []; let catIndex = 0; const catColor = new Map();
        for (const st of strips) {
            const cats = [...new Set(st.vals)];
            const rowZ = [], rowT = [];
            for (const v of st.vals) {
                const key = st.key + ':' + v;
                if (!catColor.has(key)) catColor.set(key, { i: catIndex++, color: fixed[v] || palette[catColor.size % palette.length] });
                rowZ.push(catColor.get(key).i); rowT.push(`${st.label}: ${v || 'n/a'}`);
            }
            stripZ.push(rowZ); stripText.push(rowT);
        }
        const nCat = Math.max(1, catIndex);
        const catScale = [...catColor.values()].sort((a, b) => a.i - b.i).flatMap(c => [[c.i / nCat, c.color], [(c.i + 1) / nCat, c.color]]);
        const stripH = strips.length * 16;
        const mapH = Math.max(320, Math.min(900, 40 + rows.length * 14));
        const total = mapH + stripH + (showX ? 120 : 40) + 120;
        const stripFrac = stripH / total, xlabFrac = (showX ? 110 : 30) / total;
        const yMain = [stripFrac + xlabFrac + 0.02, 1], yStrip = [xlabFrac, xlabFrac + stripFrac];
        const lab = (grp, n, l) => `<b>Group ${grp}</b> (${n} lines)${l ? ': ' + l : ''}`;
        const layout = {
            title: { text: `${this.cleanName(name)}<br><span style="font-size:11px;color:#6b7280">${metricName}; \u2605 = leading edge; rows sorted by mean A minus mean B</span>`, font: { size: 14 } },
            xaxis: { domain: [0, 0.86], tickangle: -60, showticklabels: showX, tickfont: { size: 9 }, side: 'bottom', anchor: 'y2' },
            xaxis2: { domain: [0.875, 0.875 + Math.max(0.03, Math.min(0.07, 2 * 0.86 / cols.length * 1.3))], tickfont: { size: 9 }, side: 'top', anchor: 'y' },
            yaxis: { domain: yMain, autorange: 'reversed', tickfont: { size: 10 }, automargin: true },
            yaxis2: { domain: yStrip, autorange: 'reversed', tickfont: { size: 9 }, automargin: true },
            margin: { l: 120, r: 70, t: 110, b: showX ? 110 : 30 },
            height: total,
            shapes: !['value', 'le', 'disease', 'subtype', 'lineage'].includes(sortBy)
                ? [{ type: 'line', x0: nA - 0.5, x1: nA - 0.5, y0: 0, y1: 1, xref: 'x', yref: 'paper', line: { color: '#111', width: 2 } }] : [],
            annotations: [
                { text: lab('A', nA, this._escText(d.A.label)), x: 0.0, y: 1.0, xref: 'paper', yref: 'paper', xanchor: 'left', yanchor: 'bottom', showarrow: false, font: { size: 11, color: '#dc2626' } },
                { text: lab('B', nB, this._escText(d.B.label)), x: 0.86, y: 1.0, xref: 'paper', yref: 'paper', xanchor: 'right', yanchor: 'bottom', showarrow: false, font: { size: 11, color: '#2563eb' } }
            ],
            paper_bgcolor: '#fff', plot_bgcolor: '#fff', font: { family: this.settings.fontFamily + ', sans-serif' }
        };
        const main = { type: 'heatmap', z, x: colNames, y: yLabels, zmin: -zmax, zmax, zmid: 0, colorscale: 'RdBu', reversescale: true,
            colorbar: { title: { text: metricName.length > 24 ? 'value' : metricName, side: 'right' }, thickness: 12, x: 1.0, len: yMain[1] - yMain[0], y: (yMain[0] + yMain[1]) / 2 },
            hovertemplate: '%{y}<br>%{x}<br>%{z:.2f}<extra></extra>', hoverongaps: false };
        // Means of many z-scores are small numbers; on the map's scale they all
        // looked pale, so the two mean columns use their own range (mean A vs B).
        const mAbs = Math.max(0.05, ...rows.flatMap(r => [r.mA, r.mB]).filter(v => v !== null).map(Math.abs));
        const means = { type: 'heatmap', z: rows.map(r => [r.mA, r.mB]), x: ['mean A', 'mean B'], y: yLabels, xaxis: 'x2', yaxis: 'y', zmin: -mAbs, zmax: mAbs, zmid: 0,
            colorscale: 'RdBu', reversescale: true, showscale: false, hovertemplate: '%{y}<br>%{x}: %{z:.2f} (own colour range \u00b1' + mAbs.toFixed(2) + ')<extra></extra>', hoverongaps: false };
        const stripTrace = { type: 'heatmap', z: stripZ, x: colNames, y: strips.map(st => st.label), text: stripText, xaxis: 'x', yaxis: 'y2',
            zmin: 0, zmax: nCat, colorscale: catScale, showscale: false, hovertemplate: '%{x}<br>%{text}<extra></extra>', xgap: 0.5, ygap: 2 };
        Plotly.newPlot('gsHeatmap', [main, means, stripTrace], layout, { responsive: true, displayModeBar: false, displaylogo: false });
        // legend for the strips
        const legend = document.getElementById('hmLegend');
        legend.innerHTML = strips.map(st => {
            const cats = [...new Set(st.vals)];
            return `<span style="margin-right: 12px;"><b>${st.label}:</b> ` + cats.map(v => `<span style="display:inline-block; width:10px; height:10px; background:${catColor.get(st.key + ':' + v).color}; border:1px solid #ccc; vertical-align:middle; margin: 0 3px 0 6px;"></span>${this._escText(v || 'n/a')}`).join('') + '</span>';
        }).join('');
        const sortNote = ['group', 'mutation', 'disease', 'subtype', 'lineage'].includes(sortBy) ? '' : `, sorted by ${scoreLabel} ${desc ? 'high to low' : 'low to high'}${sortBy.endsWith('Within') ? ' within each group' : ''}`;
        document.getElementById('hmInfo').textContent = `${rows.length} genes of ${result.size} in the set${leOnly ? ' (leading edge only)' : ''}; ${cols.length} cell lines${sortNote}. The two columns on the right are the mean of each group, on their own colour range (\u00b1${mAbs.toFixed(2)}) so small differences stay visible.`;
    }
});
