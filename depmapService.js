//
// Enrich, DepMap cell line service
// MIT Open source
// -
// Reads one cell line at a time from two binary matrices built by
// prepare_depmap_matrix.py (int16, row-major, one row per cell line):
//   web_data/depmap_crispr.bin      Chronos gene effect x 1000
//   web_data/depmap_expression.bin  log2FC vs the DepMap median x 1000
// with web_data/depmap_index.json naming the genes and cell lines.
// A row is fetched with an HTTP Range request (about 37 KB). A server that
// ignores Range (python's http.server does) sends the whole file; that is
// kept in memory so later rows are free.
//

const DEPMAP = { index: null, loading: null, full: {}, v: "26Q1b" };

async function DEPMAP_index() {
    if (DEPMAP.index) return DEPMAP.index;
    if (!DEPMAP.loading) {
        DEPMAP.loading = fetch(`web_data/depmap_index.json?v=${DEPMAP.v}`)
            .then(r => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.json(); })
            .then(j => {
                // One list of cell lines across both data types, for the search box.
                const byId = new Map();
                for (const type of ['crispr', 'expression']) {
                    (j[type]?.cellLines || []).forEach((c, i) => {
                        const e = byId.get(c.id) || { ...c, rows: {} };
                        e.rows[type] = i;
                        byId.set(c.id, e);
                    });
                }
                j.all = [...byId.values()].sort((a, b) => a.name.localeCompare(b.name, undefined, { numeric: true }));
                DEPMAP.index = j;
                return j;
            });
    }
    return DEPMAP.loading;
}

// Case-insensitive match on name, lineage or disease; names that start with
// the query come first.
function DEPMAP_search(query, limit = 25) {
    const idx = DEPMAP.index;
    if (!idx) return [];
    const q = String(query || '').trim().toLowerCase();
    if (!q) return { total: 0, rows: [] };
    // Names are compared without punctuation or spaces, so a375, A-375 and
    // "a 375" all find A-375; tissue and disease are matched word by word.
    const norm = (t) => String(t).toLowerCase().replace(/[^a-z0-9]/g, '');
    const qn = norm(q);
    const words = q.split(/\s+/);
    const hit = (c) => {
        if (qn && norm(c.name).includes(qn)) return true;
        const hay = `${c.name} ${c.lineage} ${c.disease} ${c.subtype || ''}`.toLowerCase();
        return words.every(w => hay.includes(w));
    };
    const out = idx.all.filter(hit);
    out.sort((a, b) => {
        const sa = norm(a.name).startsWith(qn) ? 0 : 1;
        const sb = norm(b.name).startsWith(qn) ? 0 : 1;
        return sa - sb || a.name.localeCompare(b.name, undefined, { numeric: true });
    });
    return { total: out.length, rows: out.slice(0, limit) };
}

// Rows of { Gene, <metric column> } for one cell line, like an uploaded file.
async function DEPMAP_row(type, rowIndex) {
    const idx = await DEPMAP_index();
    const t = idx[type];
    if (!t) throw new Error(`unknown data type ${type}`);
    const url = `web_data/depmap_${type}.bin?v=${DEPMAP.v}`;
    const start = rowIndex * t.rowBytes, end = start + t.rowBytes - 1;
    let bytes;
    if (DEPMAP.full[type]) {
        bytes = DEPMAP.full[type].slice(start, end + 1);
    } else {
        const r = await fetch(url, { headers: { Range: `bytes=${start}-${end}` } });
        if (!r.ok && r.status !== 206) throw new Error(`HTTP ${r.status}`);
        const buf = await r.arrayBuffer();
        if (r.status === 206 && buf.byteLength === t.rowBytes) {
            bytes = buf;
        } else {
            // Whole file came back: keep it, slice the row out.
            DEPMAP.full[type] = buf;
            bytes = buf.slice(start, end + 1);
        }
    }
    const vals = new Int16Array(bytes);
    const rows = [];
    for (let i = 0; i < t.genes.length; i++) {
        const v = vals[i];
        if (v === idx.na) continue;
        rows.push({ Gene: t.genes[i], [t.metric]: v / idx.scale });
    }
    rows.sort((a, b) => b[t.metric] - a[t.metric]);
    return rows;
}

// One cell line as a Float32Array aligned to the gene list (NaN where DepMap
// has no value), for comparisons between cell lines.
async function DEPMAP_rowValues(type, rowIndex) {
    const idx = await DEPMAP_index();
    const t = idx[type];
    const url = `web_data/depmap_${type}.bin?v=${DEPMAP.v}`;
    const start = rowIndex * t.rowBytes, end = start + t.rowBytes - 1;
    let bytes;
    if (DEPMAP.full[type]) bytes = DEPMAP.full[type].slice(start, end + 1);
    else {
        const r = await fetch(url, { headers: { Range: `bytes=${start}-${end}` } });
        if (!r.ok && r.status !== 206) throw new Error(`HTTP ${r.status}`);
        const buf = await r.arrayBuffer();
        if (r.status === 206 && buf.byteLength === t.rowBytes) bytes = buf;
        else { DEPMAP.full[type] = buf; bytes = buf.slice(start, end + 1); }
    }
    const vals = new Int16Array(bytes);
    const out = new Float32Array(vals.length);
    for (let i = 0; i < vals.length; i++) out[i] = vals[i] === idx.na ? NaN : vals[i] / idx.scale;
    return out;
}


// Hotspot and damaging mutation calls per gene (DepMap somatic mutation
// matrices), loaded when the compare dialog needs them.
async function DEPMAP_mutations() {
    if (DEPMAP.mut) return DEPMAP.mut;
    if (!DEPMAP.mutLoading) {
        DEPMAP.mutLoading = fetch(`web_data/depmap_mutations.json?v=${DEPMAP.v}`)
            .then(r => { if (!r.ok) throw new Error(`HTTP ${r.status}`); return r.json(); })
            .then(m => {
                m.profiledSet = new Set(m.profiled);
                m.hotspotSets = {}; for (const g in m.hotspot) m.hotspotSets[g] = new Set(m.hotspot[g]);
                m.damagingSets = {}; for (const g in m.damaging) m.damagingSets[g] = new Set(m.damaging[g]);
                DEPMAP.mut = m; return m;
            });
    }
    return DEPMAP.mutLoading;
}
