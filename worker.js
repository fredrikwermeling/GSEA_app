// ============================================================
// GSEA Web Worker — Preranked Gene Set Enrichment Analysis
// Implements fgsea-style weighted enrichment score computation
// with gene-label permutation testing.
// ============================================================

self.onmessage = function(e) {
    if (e.data.type === 'run') {
        try {
            runGSEA(e.data);
        } catch (err) {
            self.postMessage({ type: 'error', message: err.message });
        }
    }
};

// ------------------------------------------------------------
// Main GSEA pipeline
// ------------------------------------------------------------
function runGSEA({ rankedGenes, rankedMetrics, geneSets, settings }) {
    const { permutations, minSize, maxSize, weightP } = settings;
    const N = rankedGenes.length;

    // Build gene -> rank index lookup (all uppercase)
    const geneRankMap = new Map();
    for (let i = 0; i < N; i++) {
        geneRankMap.set(rankedGenes[i], i);
    }

    // Pre-compute |metric|^p for weighted ES
    const absMetricP = new Float64Array(N);
    for (let i = 0; i < N; i++) {
        absMetricP[i] = Math.pow(Math.abs(rankedMetrics[i]), weightP);
    }

    // Filter gene sets by size, map gene names to ranked indices
    const validSets = [];
    const totalInput = Object.keys(geneSets).length;
    let excludedTooSmall = 0;
    let excludedTooLarge = 0;
    let excludedNoOverlap = 0;
    for (const [name, genes] of Object.entries(geneSets)) {
        const hitIndices = [];
        for (const g of genes) {
            const idx = geneRankMap.get(g.toUpperCase());
            if (idx !== undefined) hitIndices.push(idx);
        }
        if (hitIndices.length === 0) {
            excludedNoOverlap++;
        } else if (hitIndices.length < minSize) {
            excludedTooSmall++;
        } else if (hitIndices.length > maxSize) {
            excludedTooLarge++;
        } else {
            hitIndices.sort((a, b) => a - b);
            validSets.push({ name, hitIndices });
        }
    }

    const totalSets = validSets.length;
    if (totalSets === 0) {
        self.postMessage({
            type: 'error',
            message: 'No gene sets passed the size filter. Try adjusting min/max size or check gene name format.'
        });
        return;
    }

    // Build exclusion summary
    const exclusions = [];
    if (excludedTooSmall > 0) exclusions.push(`${excludedTooSmall} too few genes (<${minSize})`);
    if (excludedTooLarge > 0) exclusions.push(`${excludedTooLarge} too many genes (>${maxSize})`);
    if (excludedNoOverlap > 0) exclusions.push(`${excludedNoOverlap} no matching genes`);
    const excludedTotal = totalInput - totalSets;

    self.postMessage({
        type: 'progress',
        percent: 2,
        text: `${totalSets} of ${totalInput} gene sets passed size filter${excludedTotal > 0 ? ' (' + exclusions.join(', ') + ')' : ''}. Computing enrichment scores...`,
        sizeFilterInfo: { totalInput, totalPassed: totalSets, excludedTooSmall, excludedTooLarge, excludedNoOverlap }
    });

    // ---- STEP 1: Observed enrichment scores ----
    const observedResults = [];
    for (let si = 0; si < totalSets; si++) {
        const gs = validSets[si];
        const { es, runningES, leadingEdgeIndices } = computeES(
            gs.hitIndices, absMetricP, N
        );

        // Map leading edge indices back to gene names
        const leadingEdge = leadingEdgeIndices.map(idx => rankedGenes[idx]);

        observedResults.push({
            name: gs.name,
            es,
            runningES,
            hits: new Int32Array(gs.hitIndices),
            size: gs.hitIndices.length,
            leadingEdge
        });

        if (si % 100 === 0) {
            self.postMessage({
                type: 'progress',
                percent: 2 + (si / totalSets) * 8,
                text: `Computing ES: ${si + 1}/${totalSets} gene sets...`
            });
        }
    }

    // ---- STEP 2: Adaptive permutation testing ----
    // fgsea-style: stop testing a gene set early once we have enough evidence
    // it's non-significant. Significant hits run full permutations for max precision.
    const permStats = observedResults.map(r => ({
        nPosPerms: 0,
        sumPosPerms: 0,
        nNegPerms: 0,
        sumNegPerms: 0,
        nMoreExtreme: 0
    }));

    // Adaptive stopping: gene sets stop when nMoreExtreme >= earlyStopCount
    // This means clearly non-significant sets (p > ~0.05) stop early,
    // while significant/borderline sets run full permutations.
    const active = new Uint8Array(totalSets).fill(1);
    let activeSets = totalSets;
    const earlyStopCount = 50;
    const checkInterval = 100;

    for (let perm = 0; perm < permutations; perm++) {
        for (let si = 0; si < totalSets; si++) {
            if (!active[si]) continue;

            const setSize = observedResults[si].size;
            const randomHits = sampleWithoutReplacement(N, setSize);
            randomHits.sort((a, b) => a - b);

            const permES = computeESFast(randomHits, absMetricP, N);
            const stats = permStats[si];
            const obsES = observedResults[si].es;

            if (permES >= 0) {
                stats.nPosPerms++;
                stats.sumPosPerms += permES;
            } else {
                stats.nNegPerms++;
                stats.sumNegPerms += permES;
            }

            // Count permuted ES more extreme than observed
            if (obsES >= 0) {
                if (permES >= obsES) stats.nMoreExtreme++;
            } else {
                if (permES <= obsES) stats.nMoreExtreme++;
            }
        }

        // Check adaptive stopping periodically (not on last perm — let it finish)
        if ((perm + 1) % checkInterval === 0 && perm + 1 < permutations) {
            for (let si = 0; si < totalSets; si++) {
                if (!active[si]) continue;
                if (permStats[si].nMoreExtreme >= earlyStopCount) {
                    active[si] = 0;
                    activeSets--;
                }
            }
        }

        if (perm % 10 === 0 || perm === permutations - 1) {
            const refining = activeSets < totalSets
                ? ` (${activeSets} of ${totalSets} still refining)`
                : '';
            self.postMessage({
                type: 'progress',
                percent: 10 + (perm / permutations) * 80,
                text: `Permutation ${perm + 1}/${permutations}${refining}`
            });
        }
    }

    // ---- STEP 3: NES and p-values ----
    self.postMessage({ type: 'progress', percent: 92, text: 'Computing NES and p-values...' });

    for (let si = 0; si < totalSets; si++) {
        const obs = observedResults[si];
        const stats = permStats[si];

        // NES normalization
        if (obs.es >= 0) {
            const meanPos = stats.nPosPerms > 0
                ? stats.sumPosPerms / stats.nPosPerms
                : 1;
            obs.nes = meanPos > 0 ? obs.es / meanPos : 0;
        } else {
            const meanNeg = stats.nNegPerms > 0
                ? stats.sumNegPerms / stats.nNegPerms
                : -1;
            obs.nes = meanNeg < 0 ? obs.es / Math.abs(meanNeg) : 0;
        }

        // Empirical p-value with Laplace correction
        if (obs.es >= 0) {
            obs.pvalue = stats.nPosPerms > 0
                ? (stats.nMoreExtreme + 1) / (stats.nPosPerms + 1)
                : 1;
        } else {
            obs.pvalue = stats.nNegPerms > 0
                ? (stats.nMoreExtreme + 1) / (stats.nNegPerms + 1)
                : 1;
        }

        // Cap p-value at 1
        obs.pvalue = Math.min(obs.pvalue, 1);
    }

    // ---- STEP 4: FDR (Benjamini-Hochberg) ----
    self.postMessage({ type: 'progress', percent: 96, text: 'Computing FDR...' });

    // Sort by p-value ascending
    const sorted = observedResults.slice().sort((a, b) => a.pvalue - b.pvalue);
    const m = sorted.length;

    // BH procedure with monotonicity enforcement
    for (let i = m - 1; i >= 0; i--) {
        const bh = (sorted[i].pvalue * m) / (i + 1);
        sorted[i].fdr = i < m - 1
            ? Math.min(bh, sorted[i + 1].fdr)
            : Math.min(bh, 1);
        sorted[i].fdr = Math.min(sorted[i].fdr, 1);
    }

    // Copy FDR back to observedResults
    for (const r of sorted) {
        const obs = observedResults.find(o => o.name === r.name);
        if (obs) obs.fdr = r.fdr;
    }

    // Sort by |NES| descending for final output
    observedResults.sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes));

    self.postMessage({ type: 'progress', percent: 100, text: 'Done!' });

    // Convert runningES from Float64Array to regular array for transfer
    self.postMessage({
        type: 'complete',
        results: observedResults.map(r => ({
            name: r.name,
            es: r.es,
            nes: r.nes,
            pvalue: r.pvalue,
            fdr: r.fdr,
            size: r.size,
            leadingEdge: r.leadingEdge,
            runningES: Array.from(r.runningES),
            hits: Array.from(r.hits)
        }))
    });
}

// ------------------------------------------------------------
// Compute enrichment score with full tracking
// Returns: es, runningES array, leading edge indices
// ------------------------------------------------------------
function computeES(hitIndices, absMetricP, N) {
    const NH = hitIndices.length;
    const hitSet = new Set(hitIndices);

    // NR = sum of |metric|^p for genes in the set
    let NR = 0;
    for (const idx of hitIndices) {
        NR += absMetricP[idx];
    }
    if (NR === 0) NR = 1;

    const missPenalty = 1 / (N - NH);

    const runningES = new Float64Array(N);
    let maxES = -Infinity;
    let minES = Infinity;
    let maxIdx = 0;
    let minIdx = 0;
    let cumul = 0;

    for (let i = 0; i < N; i++) {
        if (hitSet.has(i)) {
            cumul += absMetricP[i] / NR;
        } else {
            cumul -= missPenalty;
        }
        runningES[i] = cumul;

        if (cumul > maxES) { maxES = cumul; maxIdx = i; }
        if (cumul < minES) { minES = cumul; minIdx = i; }
    }

    // ES = maximum deviation from zero
    const es = Math.abs(maxES) >= Math.abs(minES) ? maxES : minES;
    const peakIdx = Math.abs(maxES) >= Math.abs(minES) ? maxIdx : minIdx;

    // Leading edge: hits contributing to ES
    const leadingEdgeIndices = [];
    if (es >= 0) {
        for (const idx of hitIndices) {
            if (idx <= peakIdx) leadingEdgeIndices.push(idx);
        }
    } else {
        for (const idx of hitIndices) {
            if (idx >= peakIdx) leadingEdgeIndices.push(idx);
        }
    }

    return { es, runningES, leadingEdgeIndices };
}

// ------------------------------------------------------------
// Fast ES computation for permutations (no running array)
// Uses sorted-hit-pointer for cache-friendly inner loop
// ------------------------------------------------------------
function computeESFast(sortedHitIndices, absMetricP, N) {
    const NH = sortedHitIndices.length;
    const missPenalty = 1 / (N - NH);

    let NR = 0;
    for (let i = 0; i < NH; i++) {
        NR += absMetricP[sortedHitIndices[i]];
    }
    if (NR === 0) NR = 1;

    let cumul = 0;
    let maxES = -Infinity;
    let minES = Infinity;
    let hitPtr = 0;

    for (let i = 0; i < N; i++) {
        if (hitPtr < NH && sortedHitIndices[hitPtr] === i) {
            cumul += absMetricP[i] / NR;
            hitPtr++;
        } else {
            cumul -= missPenalty;
        }
        if (cumul > maxES) maxES = cumul;
        if (cumul < minES) minES = cumul;
    }

    return Math.abs(maxES) >= Math.abs(minES) ? maxES : minES;
}

// ------------------------------------------------------------
// Sample k indices from [0, n) without replacement
// Partial Fisher-Yates with sparse map to avoid O(n) allocation
// ------------------------------------------------------------
function sampleWithoutReplacement(n, k) {
    const swapped = new Map();
    const result = new Int32Array(k);

    for (let i = 0; i < k; i++) {
        const j = i + Math.floor(Math.random() * (n - i));
        const valJ = swapped.has(j) ? swapped.get(j) : j;
        const valI = swapped.has(i) ? swapped.get(i) : i;
        result[i] = valJ;
        swapped.set(j, valI);
        swapped.delete(i);
    }
    return result;
}
