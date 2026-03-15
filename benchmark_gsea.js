#!/usr/bin/env node
/**
 * GSEA Benchmark: Run GSEA vs Smart Run
 *
 * Tests speed and accuracy at different permutation counts and gene set sizes.
 * Uses the actual worker.js algorithm to ensure results are representative.
 *
 * Usage: node benchmark_gsea.js
 */

const fs = require('fs');
const path = require('path');

// ============================================================
// Import GSEA functions from worker.js (copy core functions)
// ============================================================

function computeES(hitIndices, absMetricP, N) {
    const NH = hitIndices.length;
    const hitSet = new Set(hitIndices);
    let NR = 0;
    for (const idx of hitIndices) NR += absMetricP[idx];
    if (NR === 0) NR = 1;
    const missPenalty = 1 / (N - NH);
    const runningES = new Float64Array(N);
    let maxES = -Infinity, minES = Infinity, maxIdx = 0, minIdx = 0, cumul = 0;
    for (let i = 0; i < N; i++) {
        if (hitSet.has(i)) cumul += absMetricP[i] / NR;
        else cumul -= missPenalty;
        runningES[i] = cumul;
        if (cumul > maxES) { maxES = cumul; maxIdx = i; }
        if (cumul < minES) { minES = cumul; minIdx = i; }
    }
    const es = Math.abs(maxES) >= Math.abs(minES) ? maxES : minES;
    const peakIdx = Math.abs(maxES) >= Math.abs(minES) ? maxIdx : minIdx;
    const leadingEdgeIndices = [];
    if (es >= 0) { for (const idx of hitIndices) if (idx <= peakIdx) leadingEdgeIndices.push(idx); }
    else { for (const idx of hitIndices) if (idx >= peakIdx) leadingEdgeIndices.push(idx); }
    return { es, runningES, leadingEdgeIndices };
}

function computeESFast(sortedHitIndices, absMetricP, N) {
    const NH = sortedHitIndices.length;
    const missPenalty = 1 / (N - NH);
    let NR = 0;
    for (let i = 0; i < NH; i++) NR += absMetricP[sortedHitIndices[i]];
    if (NR === 0) NR = 1;
    let cumul = 0, maxES = -Infinity, minES = Infinity, hitPtr = 0;
    for (let i = 0; i < N; i++) {
        if (hitPtr < NH && sortedHitIndices[hitPtr] === i) { cumul += absMetricP[i] / NR; hitPtr++; }
        else cumul -= missPenalty;
        if (cumul > maxES) maxES = cumul;
        if (cumul < minES) minES = cumul;
    }
    return Math.abs(maxES) >= Math.abs(minES) ? maxES : minES;
}

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

function runGSEADirect(rankedGenes, rankedMetrics, geneSets, permutations, minSize = 15, maxSize = 500) {
    const N = rankedGenes.length;
    const geneRankMap = new Map();
    for (let i = 0; i < N; i++) geneRankMap.set(rankedGenes[i], i);

    const absMetricP = new Float64Array(N);
    for (let i = 0; i < N; i++) absMetricP[i] = Math.abs(rankedMetrics[i]);

    const validSets = [];
    for (const [name, genes] of Object.entries(geneSets)) {
        const hitIndices = [];
        for (const g of genes) {
            const idx = geneRankMap.get(g.toUpperCase());
            if (idx !== undefined) hitIndices.push(idx);
        }
        if (hitIndices.length >= minSize && hitIndices.length <= maxSize) {
            hitIndices.sort((a, b) => a - b);
            validSets.push({ name, hitIndices });
        }
    }

    // Step 1: Observed ES
    const observedResults = [];
    for (const gs of validSets) {
        const { es, leadingEdgeIndices } = computeES(gs.hitIndices, absMetricP, N);
        observedResults.push({
            name: gs.name, es, size: gs.hitIndices.length,
            leadingEdge: leadingEdgeIndices.map(idx => rankedGenes[idx])
        });
    }

    // Step 2: Permutation testing (with adaptive stopping)
    const totalSets = observedResults.length;
    const permStats = observedResults.map(() => ({
        nPosPerms: 0, sumPosPerms: 0, nNegPerms: 0, sumNegPerms: 0, nMoreExtreme: 0
    }));
    const active = new Uint8Array(totalSets).fill(1);
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
            if (permES >= 0) { stats.nPosPerms++; stats.sumPosPerms += permES; }
            else { stats.nNegPerms++; stats.sumNegPerms += permES; }
            if (obsES >= 0) { if (permES >= obsES) stats.nMoreExtreme++; }
            else { if (permES <= obsES) stats.nMoreExtreme++; }
        }
        if ((perm + 1) % checkInterval === 0 && perm + 1 < permutations) {
            for (let si = 0; si < totalSets; si++) {
                if (!active[si]) continue;
                if (permStats[si].nMoreExtreme >= earlyStopCount) active[si] = 0;
            }
        }
    }

    // Step 3: NES and p-values
    for (let si = 0; si < totalSets; si++) {
        const obs = observedResults[si];
        const stats = permStats[si];
        if (obs.es >= 0) {
            const meanPos = stats.nPosPerms > 0 ? stats.sumPosPerms / stats.nPosPerms : 1;
            obs.nes = meanPos > 0 ? obs.es / meanPos : 0;
        } else {
            const meanNeg = stats.nNegPerms > 0 ? stats.sumNegPerms / stats.nNegPerms : -1;
            obs.nes = meanNeg < 0 ? obs.es / Math.abs(meanNeg) : 0;
        }
        if (obs.es >= 0) {
            obs.pvalue = stats.nPosPerms > 0 ? (stats.nMoreExtreme + 1) / (stats.nPosPerms + 1) : 1;
        } else {
            obs.pvalue = stats.nNegPerms > 0 ? (stats.nMoreExtreme + 1) / (stats.nNegPerms + 1) : 1;
        }
        obs.pvalue = Math.min(obs.pvalue, 1);
    }

    // Step 4: FDR (BH)
    const sorted = observedResults.slice().sort((a, b) => a.pvalue - b.pvalue);
    const m = sorted.length;
    for (let i = m - 1; i >= 0; i--) {
        const bh = (sorted[i].pvalue * m) / (i + 1);
        sorted[i].fdr = i < m - 1 ? Math.min(bh, sorted[i + 1].fdr) : Math.min(bh, 1);
        sorted[i].fdr = Math.min(sorted[i].fdr, 1);
    }
    for (const r of sorted) {
        const obs = observedResults.find(o => o.name === r.name);
        if (obs) obs.fdr = r.fdr;
    }

    return observedResults;
}

// ============================================================
// Smart Run simulation (3-phase)
// ============================================================
function computeJaccard(genesA, genesB) {
    const setB = new Set(genesB);
    let inter = 0;
    for (const g of genesA) if (setB.has(g)) inter++;
    return inter / (genesA.length + genesB.length - inter);
}

function runSmartGSEA(rankedGenes, rankedMetrics, allGeneSets, userPerms, minSize = 15, maxSize = 500) {
    const nTotal = Object.keys(allGeneSets).length;

    // Phase 1: Diverse subset (Jaccard < 0.1 from already-selected)
    const diverseSets = {};
    const diverseGenes = [];
    const entries = Object.entries(allGeneSets);
    for (const [name, genes] of entries) {
        let isDiverse = true;
        for (const prev of diverseGenes) {
            if (computeJaccard(genes, prev) > 0.1) { isDiverse = false; break; }
        }
        if (isDiverse) {
            diverseSets[name] = genes;
            diverseGenes.push(genes);
        }
    }

    const phase1Results = runGSEADirect(rankedGenes, rankedMetrics, diverseSets, 100, minSize, maxSize);
    const phase1Hits = phase1Results.filter(r => r.fdr < 0.25);

    // Phase 2: Expand around hits
    let allResults = [...phase1Results];
    if (phase1Hits.length > 0) {
        const hitNames = new Set(phase1Hits.map(r => r.name));
        const phase1Names = new Set(phase1Results.map(r => r.name));
        const expansionSets = {};
        for (const [name, genes] of entries) {
            if (phase1Names.has(name)) continue;
            for (const hit of phase1Hits) {
                if (allGeneSets[hit.name] && computeJaccard(genes, allGeneSets[hit.name]) > 0.1) {
                    expansionSets[name] = genes;
                    break;
                }
            }
        }
        if (Object.keys(expansionSets).length > 0) {
            const phase2Results = runGSEADirect(rankedGenes, rankedMetrics, expansionSets, 100, minSize, maxSize);
            const existingNames = new Set(allResults.map(r => r.name));
            for (const r of phase2Results) {
                if (!existingNames.has(r.name)) allResults.push(r);
            }
        }
    }

    // Phase 3: Refine hits with high permutations
    const allHits = allResults.filter(r => r.fdr < 0.25);
    if (allHits.length > 0) {
        const refinePerm = Math.max(1000, userPerms);
        const refineSets = {};
        for (const r of allHits) {
            if (allGeneSets[r.name]) refineSets[r.name] = allGeneSets[r.name];
        }
        const refined = runGSEADirect(rankedGenes, rankedMetrics, refineSets, refinePerm, minSize, maxSize);
        const refinedByName = {};
        for (const r of refined) refinedByName[r.name] = r;
        allResults = allResults.map(r => refinedByName[r.name] || r);
    }

    return allResults;
}

// ============================================================
// Benchmark runner
// ============================================================
function loadJSON(filepath) {
    return JSON.parse(fs.readFileSync(filepath, 'utf-8'));
}

function run() {
    console.log('='.repeat(70));
    console.log('GSEA BENCHMARK: Run GSEA vs Smart Run');
    console.log('='.repeat(70));
    console.log();

    // Load example data
    const dataDir = path.join(__dirname, 'web_data');
    const exprData = loadJSON(path.join(dataDir, 'depmap_expression_A375.json'));
    const rankedGenes = exprData.map(d => d.Gene.toUpperCase());
    const rankedMetrics = exprData.map(d => d.Expression_zscore);

    console.log(`Ranked list: ${rankedGenes.length} genes`);
    console.log(`Metric range: [${Math.min(...rankedMetrics).toFixed(2)}, ${Math.max(...rankedMetrics).toFixed(2)}]`);
    console.log();

    // Load gene set collections
    const hallmark = loadJSON(path.join(dataDir, 'h.all.v2023.2.Hs.json'));
    const c2kegg = loadJSON(path.join(dataDir, 'c2.cp.kegg.v2023.2.Hs.json'));
    const c2reactome = loadJSON(path.join(dataDir, 'c2.cp.reactome.v2023.2.Hs.json'));

    const collections = {
        'Hallmark (50 sets)': hallmark,
        'Hallmark + KEGG (~855 sets)': { ...hallmark, ...c2kegg },
        'Hallmark + KEGG + Reactome (~2547 sets)': { ...hallmark, ...c2kegg, ...c2reactome }
    };

    // ---- Ground truth: 10,000 permutations on Hallmark ----
    console.log('─'.repeat(70));
    console.log('GROUND TRUTH: Hallmark with 10,000 permutations');
    console.log('─'.repeat(70));
    const t0 = Date.now();
    const groundTruth = runGSEADirect(rankedGenes, rankedMetrics, hallmark, 10000);
    const gtTime = Date.now() - t0;
    console.log(`  Time: ${(gtTime / 1000).toFixed(1)}s`);
    console.log(`  Sets tested: ${groundTruth.length}`);
    const gtSig = groundTruth.filter(r => r.fdr < 0.25);
    console.log(`  Significant (FDR < 0.25): ${gtSig.length}`);
    console.log(`  Significant (FDR < 0.05): ${groundTruth.filter(r => r.fdr < 0.05).length}`);
    console.log();

    // Build ground truth lookup
    const gtByName = {};
    for (const r of groundTruth) gtByName[r.name] = r;

    // ---- Test different permutation counts on Hallmark ----
    console.log('─'.repeat(70));
    console.log('ACCURACY vs PERMUTATION COUNT (Hallmark, 50 sets)');
    console.log('Compared against 10,000-permutation ground truth');
    console.log('─'.repeat(70));
    console.log();

    const permCounts = [100, 200, 500, 1000, 2000, 5000];
    const nTrials = 5;

    console.log(`${'Perms'.padStart(7)} | ${'Time(s)'.padStart(8)} | ${'Sig025'.padStart(6)} | ${'Sig005'.padStart(6)} | ${'NES_r'.padStart(7)} | ${'FDR RMSE'.padStart(9)} | ${'Missed'.padStart(7)} | ${'False+'.padStart(7)}`);
    console.log('-'.repeat(78));

    for (const nPerm of permCounts) {
        const times = [];
        const nesCorrs = [];
        const fdrRMSEs = [];
        const sigCounts025 = [];
        const sigCounts005 = [];
        const missedCounts = [];
        const falsePCounts = [];

        for (let trial = 0; trial < nTrials; trial++) {
            const t1 = Date.now();
            const results = runGSEADirect(rankedGenes, rankedMetrics, hallmark, nPerm);
            times.push(Date.now() - t1);

            // Match results to ground truth
            const byName = {};
            for (const r of results) byName[r.name] = r;

            let sumSqFDR = 0, nCompared = 0;
            const nesObs = [], nesGT = [];
            let missed = 0, falsePos = 0;

            for (const gt of groundTruth) {
                const test = byName[gt.name];
                if (!test) continue;
                nesObs.push(test.nes);
                nesGT.push(gt.nes);
                sumSqFDR += (test.fdr - gt.fdr) ** 2;
                nCompared++;

                // Missed: GT says significant, test says not
                if (gt.fdr < 0.25 && test.fdr >= 0.25) missed++;
                // False positive: test says significant, GT says not
                if (test.fdr < 0.25 && gt.fdr >= 0.25) falsePos++;
            }

            // Pearson correlation of NES
            const meanObs = nesObs.reduce((a, b) => a + b, 0) / nesObs.length;
            const meanGT = nesGT.reduce((a, b) => a + b, 0) / nesGT.length;
            let num = 0, denObs = 0, denGT = 0;
            for (let i = 0; i < nesObs.length; i++) {
                num += (nesObs[i] - meanObs) * (nesGT[i] - meanGT);
                denObs += (nesObs[i] - meanObs) ** 2;
                denGT += (nesGT[i] - meanGT) ** 2;
            }
            nesCorrs.push(num / Math.sqrt(denObs * denGT));
            fdrRMSEs.push(Math.sqrt(sumSqFDR / nCompared));
            sigCounts025.push(results.filter(r => r.fdr < 0.25).length);
            sigCounts005.push(results.filter(r => r.fdr < 0.05).length);
            missedCounts.push(missed);
            falsePCounts.push(falsePos);
        }

        const avg = arr => arr.reduce((a, b) => a + b, 0) / arr.length;
        console.log(
            `${String(nPerm).padStart(7)} | ` +
            `${(avg(times) / 1000).toFixed(2).padStart(8)} | ` +
            `${avg(sigCounts025).toFixed(1).padStart(6)} | ` +
            `${avg(sigCounts005).toFixed(1).padStart(6)} | ` +
            `${avg(nesCorrs).toFixed(4).padStart(7)} | ` +
            `${avg(fdrRMSEs).toFixed(5).padStart(9)} | ` +
            `${avg(missedCounts).toFixed(1).padStart(7)} | ` +
            `${avg(falsePCounts).toFixed(1).padStart(7)}`
        );
    }

    console.log();
    console.log(`Ground truth: Sig025=${gtSig.length}, Sig005=${groundTruth.filter(r => r.fdr < 0.05).length}`);

    // ---- Speed benchmark: Run GSEA vs Smart Run on larger collections ----
    console.log();
    console.log('─'.repeat(70));
    console.log('SPEED BENCHMARK: Run GSEA (1000 perms) vs Smart Run');
    console.log('─'.repeat(70));
    console.log();

    console.log(`${'Collection'.padEnd(40)} | ${'Sets'.padStart(5)} | ${'Mode'.padEnd(10)} | ${'Time(s)'.padStart(8)} | ${'Tested'.padStart(7)} | ${'Sig025'.padStart(6)} | ${'Sig005'.padStart(6)}`);
    console.log('-'.repeat(95));

    for (const [collName, geneSets] of Object.entries(collections)) {
        const nSets = Object.keys(geneSets).length;

        // Standard Run GSEA (1000 perms)
        const t1 = Date.now();
        const stdResults = runGSEADirect(rankedGenes, rankedMetrics, geneSets, 1000);
        const stdTime = Date.now() - t1;
        const stdSig025 = stdResults.filter(r => r.fdr < 0.25).length;
        const stdSig005 = stdResults.filter(r => r.fdr < 0.05).length;

        console.log(
            `${collName.padEnd(40)} | ${String(nSets).padStart(5)} | ${'Run GSEA'.padEnd(10)} | ` +
            `${(stdTime / 1000).toFixed(2).padStart(8)} | ${String(stdResults.length).padStart(7)} | ` +
            `${String(stdSig025).padStart(6)} | ${String(stdSig005).padStart(6)}`
        );

        // Smart Run
        const t2 = Date.now();
        const smartResults = runSmartGSEA(rankedGenes, rankedMetrics, geneSets, 1000);
        const smartTime = Date.now() - t2;
        const smartSig025 = smartResults.filter(r => r.fdr < 0.25).length;
        const smartSig005 = smartResults.filter(r => r.fdr < 0.05).length;

        console.log(
            `${''.padEnd(40)} | ${String(nSets).padStart(5)} | ${'Smart Run'.padEnd(10)} | ` +
            `${(smartTime / 1000).toFixed(2).padStart(8)} | ${String(smartResults.length).padStart(7)} | ` +
            `${String(smartSig025).padStart(6)} | ${String(smartSig005).padStart(6)}`
        );

        // Compare: which significant sets did Smart Run miss?
        const stdSigNames = new Set(stdResults.filter(r => r.fdr < 0.25).map(r => r.name));
        const smartSigNames = new Set(smartResults.filter(r => r.fdr < 0.25).map(r => r.name));
        const smartTestedNames = new Set(smartResults.map(r => r.name));

        const missed = [...stdSigNames].filter(n => !smartSigNames.has(n));
        const notTested = [...stdSigNames].filter(n => !smartTestedNames.has(n));
        const extraFound = [...smartSigNames].filter(n => !stdSigNames.has(n));

        const speedup = stdTime / smartTime;
        console.log(
            `${''.padEnd(40)} |       | ${'Comparison'.padEnd(10)} | ` +
            `${speedup.toFixed(1).padStart(5)}x   | ` +
            `miss:${missed.length} ntst:${notTested.length} extra:${extraFound.length}`
        );
        console.log('-'.repeat(95));
    }

    // ---- Reproducibility test ----
    console.log();
    console.log('─'.repeat(70));
    console.log('REPRODUCIBILITY: Run-to-run variance (Hallmark, 1000 perms, 10 runs)');
    console.log('─'.repeat(70));
    console.log();

    const reproResults = [];
    for (let i = 0; i < 10; i++) {
        reproResults.push(runGSEADirect(rankedGenes, rankedMetrics, hallmark, 1000));
    }

    // For each gene set, compute NES variance across runs
    const setNames = reproResults[0].map(r => r.name);
    console.log(`${'Gene Set'.padEnd(50)} | ${'NES mean'.padStart(9)} | ${'NES SD'.padStart(8)} | ${'FDR mean'.padStart(9)} | ${'FDR SD'.padStart(8)}`);
    console.log('-'.repeat(95));

    // Show top 10 by |NES|
    const avgNES = {};
    for (const name of setNames) {
        const nesValues = reproResults.map(run => {
            const r = run.find(x => x.name === name);
            return r ? r.nes : 0;
        });
        avgNES[name] = nesValues.reduce((a, b) => a + b, 0) / nesValues.length;
    }
    const topSets = setNames.sort((a, b) => Math.abs(avgNES[b]) - Math.abs(avgNES[a])).slice(0, 15);

    for (const name of topSets) {
        const nesValues = reproResults.map(run => run.find(x => x.name === name)?.nes || 0);
        const fdrValues = reproResults.map(run => run.find(x => x.name === name)?.fdr || 1);
        const nesMean = nesValues.reduce((a, b) => a + b, 0) / nesValues.length;
        const nesSD = Math.sqrt(nesValues.reduce((a, b) => a + (b - nesMean) ** 2, 0) / nesValues.length);
        const fdrMean = fdrValues.reduce((a, b) => a + b, 0) / fdrValues.length;
        const fdrSD = Math.sqrt(fdrValues.reduce((a, b) => a + (b - fdrMean) ** 2, 0) / fdrValues.length);

        const displayName = name.replace(/^HALLMARK_/, '').replace(/_/g, ' ').substring(0, 48);
        console.log(
            `${displayName.padEnd(50)} | ${nesMean.toFixed(3).padStart(9)} | ${nesSD.toFixed(4).padStart(8)} | ` +
            `${fdrMean.toFixed(4).padStart(9)} | ${fdrSD.toFixed(4).padStart(8)}`
        );
    }

    console.log();
    console.log('='.repeat(70));
    console.log('BENCHMARK COMPLETE');
    console.log('='.repeat(70));
}

run();
