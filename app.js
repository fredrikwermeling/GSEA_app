// ============================================================
// GSEA Web Tool — Main Application
// ============================================================

class GSEAApp {
    constructor() {
        // State
        this.rawData = null;
        this.rankedList = null;
        this.geneSets = {};          // { collectionId: { setName: [genes] } }
        this.customGeneSets = {};    // uploaded GMT sets
        this.results = null;
        this.worker = null;
        this.analysisDate = null;

        // Settings
        this.settings = {
            permutations: 1000,
            minSize: 15,
            maxSize: 500,
            weightP: 1,
            topN: 20,
            fdrDisplayThreshold: 0.25,
            colorScale: 'RdBu',
            exportScale: 4,
            transparentBg: false
        };

        // Table sort state
        this.sortCol = 'nes';
        this.sortAsc = false;

        this.init();
    }

    // --------------------------------------------------------
    // Initialization
    // --------------------------------------------------------
    init() {
        this.createWorker();
        this.bindEvents();
        this.loadDefaultGeneSets();
    }

    createWorker() {
        if (this.worker) this.worker.terminate();
        this.worker = new Worker('worker.js');
        this.worker.onmessage = (e) => this.handleWorkerMessage(e);
        this.worker.onerror = (e) => this.handleWorkerError(e);
    }

    bindEvents() {
        // File upload
        document.getElementById('fileInput').addEventListener('change', (e) => {
            if (e.target.files[0]) this.handleFileUpload(e.target.files[0]);
        });

        // Example data
        document.getElementById('loadExampleBtn').addEventListener('click', () => {
            this.loadExampleData();
        });

        // Column selection
        document.getElementById('geneColumn').addEventListener('change', () => this.checkReady());
        document.getElementById('metricColumn').addEventListener('change', () => this.checkReady());

        // Gene set collection checkboxes
        ['checkHallmark', 'checkC2', 'checkC5'].forEach(id => {
            document.getElementById(id).addEventListener('change', () => this.onCollectionChange());
        });

        // Custom GMT toggle
        document.getElementById('customGmtToggle').addEventListener('click', () => {
            const header = document.getElementById('customGmtToggle');
            const body = document.getElementById('customGmtBody');
            header.classList.toggle('open');
            body.classList.toggle('open');
        });

        // Custom GMT upload
        document.getElementById('gmtInput').addEventListener('change', (e) => {
            if (e.target.files[0]) this.handleGMTUpload(e.target.files[0]);
        });

        // Run / Cancel
        document.getElementById('runBtn').addEventListener('click', () => this.runGSEA());
        document.getElementById('cancelBtn').addEventListener('click', () => this.cancelAnalysis());

        // Tabs
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.addEventListener('click', () => this.showTab(btn.dataset.tab));
        });

        // Gene set selector for ES plot
        document.getElementById('geneSetSelector').addEventListener('change', (e) => {
            if (e.target.value) this.renderESPlot(e.target.value);
        });

        // Table filter
        document.getElementById('tableFilter').addEventListener('input', () => this.filterAndRenderTable());
        document.getElementById('fdrFilter').addEventListener('change', () => this.filterAndRenderTable());
        document.getElementById('directionFilter').addEventListener('change', () => this.filterAndRenderTable());

        // Table column sort
        document.querySelectorAll('#resultsTable thead th').forEach(th => {
            th.addEventListener('click', () => this.sortTable(th.dataset.col));
        });

        // Download CSV
        document.getElementById('downloadCSVBtn').addEventListener('click', () => this.downloadCSV());

        // Copy methods
        document.getElementById('copyMethodsBtn').addEventListener('click', () => this.copyMethods());

        // Figure settings
        document.getElementById('updatePlotsBtn').addEventListener('click', () => this.updateSettings());

        // Settings inputs
        ['topN', 'fdrDisplayThreshold', 'colorScale', 'exportScale', 'transparentBg'].forEach(id => {
            const el = document.getElementById(id);
            if (el) el.addEventListener('change', () => this.readSettings());
        });
    }

    // --------------------------------------------------------
    // Gene Set Loading
    // --------------------------------------------------------
    async loadDefaultGeneSets() {
        try {
            const resp = await fetch('web_data/h.all.v2023.2.Hs.json');
            if (!resp.ok) throw new Error('Failed to load Hallmark gene sets');
            this.geneSets['hallmark'] = await resp.json();
            this.updateGeneSetStatus();
        } catch (err) {
            this.showStatus('geneSetStatus', 'error', 'Failed to load Hallmark gene sets: ' + err.message);
        }
        document.getElementById('loadingOverlay').classList.add('hidden');
    }

    async onCollectionChange() {
        const collections = {
            checkHallmark: { id: 'hallmark', file: 'h.all.v2023.2.Hs.json' },
            checkC2: { id: 'c2', file: 'c2.all.v2023.2.Hs.json' },
            checkC5: { id: 'c5', file: 'c5.all.v2023.2.Hs.json' }
        };

        for (const [checkId, info] of Object.entries(collections)) {
            const checked = document.getElementById(checkId).checked;
            if (checked && !this.geneSets[info.id]) {
                this.showStatus('geneSetStatus', 'info', `Loading ${info.id.toUpperCase()} gene sets...`);
                try {
                    const resp = await fetch('web_data/' + info.file);
                    if (!resp.ok) throw new Error('File not found');
                    this.geneSets[info.id] = await resp.json();
                } catch (err) {
                    document.getElementById(checkId).checked = false;
                    this.showStatus('geneSetStatus', 'error',
                        `Failed to load ${info.id.toUpperCase()}: ${err.message}. ` +
                        'You can upload a GMT file instead.');
                    return;
                }
            }
        }
        this.updateGeneSetStatus();
        this.checkReady();
    }

    updateGeneSetStatus() {
        const parts = [];
        let total = 0;
        const collections = {
            checkHallmark: { id: 'hallmark', label: 'Hallmark' },
            checkC2: { id: 'c2', label: 'C2' },
            checkC5: { id: 'c5', label: 'C5' }
        };
        for (const [checkId, info] of Object.entries(collections)) {
            if (document.getElementById(checkId).checked && this.geneSets[info.id]) {
                const n = Object.keys(this.geneSets[info.id]).length;
                parts.push(`${info.label}: ${n}`);
                total += n;
            }
        }
        if (Object.keys(this.customGeneSets).length > 0) {
            const n = Object.keys(this.customGeneSets).length;
            parts.push(`Custom: ${n}`);
            total += n;
        }
        if (parts.length > 0) {
            this.showStatus('geneSetStatus', 'info',
                `${parts.join(' | ')} — ${total} total gene sets`);
        } else {
            this.showStatus('geneSetStatus', 'warning', 'No gene set collections selected');
        }
    }

    // --------------------------------------------------------
    // GMT Upload
    // --------------------------------------------------------
    handleGMTUpload(file) {
        const reader = new FileReader();
        reader.onload = (e) => {
            try {
                this.customGeneSets = this.parseGMT(e.target.result);
                const n = Object.keys(this.customGeneSets).length;
                this.showStatus('gmtStatus', 'success', `Loaded ${n} custom gene sets`);
                this.updateGeneSetStatus();
                this.checkReady();
            } catch (err) {
                this.showStatus('gmtStatus', 'error', 'Failed to parse GMT: ' + err.message);
            }
        };
        reader.readAsText(file);
    }

    parseGMT(text) {
        const sets = {};
        const lines = text.split('\n');
        for (const line of lines) {
            if (!line.trim()) continue;
            const parts = line.split('\t');
            if (parts.length < 3) continue;
            const name = parts[0].trim();
            const genes = parts.slice(2)
                .map(g => g.trim().toUpperCase())
                .filter(g => g !== '');
            if (genes.length > 0) {
                sets[name] = genes;
            }
        }
        return sets;
    }

    // --------------------------------------------------------
    // File Upload & Column Selection
    // --------------------------------------------------------
    handleFileUpload(file) {
        this.showStatus('uploadStatus', 'info', 'Parsing file...');

        Papa.parse(file, {
            header: true,
            dynamicTyping: true,
            skipEmptyLines: true,
            delimiter: '',
            complete: (result) => {
                if (result.errors.length > 0 && result.data.length === 0) {
                    this.showStatus('uploadStatus', 'error',
                        'Parse error: ' + result.errors[0].message);
                    return;
                }
                this.rawData = result.data;
                this.populateColumnDropdowns(result.meta.fields);
                this.showStatus('uploadStatus', 'success',
                    `Loaded ${result.data.length} rows, ${result.meta.fields.length} columns`);
            },
            error: (err) => {
                this.showStatus('uploadStatus', 'error', 'Parse error: ' + err.message);
            }
        });
    }

    populateColumnDropdowns(fields) {
        const geneCol = document.getElementById('geneColumn');
        const metricCol = document.getElementById('metricColumn');

        geneCol.innerHTML = '<option value="">Select gene name column...</option>';
        metricCol.innerHTML = '<option value="">Select ranking metric column...</option>';

        // Guess which columns are gene names vs numeric
        const numericCols = [];
        const stringCols = [];

        for (const field of fields) {
            // Check first non-null value
            let isNumeric = false;
            for (const row of this.rawData) {
                const val = row[field];
                if (val !== null && val !== undefined && val !== '') {
                    isNumeric = typeof val === 'number' || !isNaN(parseFloat(val));
                    break;
                }
            }
            if (isNumeric) {
                numericCols.push(field);
            } else {
                stringCols.push(field);
            }
        }

        // Gene column: show all, but put string columns first
        for (const f of stringCols) {
            geneCol.appendChild(new Option(f, f));
        }
        for (const f of numericCols) {
            geneCol.appendChild(new Option(f + ' (numeric)', f));
        }

        // Metric column: show only numeric
        for (const f of numericCols) {
            metricCol.appendChild(new Option(f, f));
        }

        geneCol.disabled = false;
        metricCol.disabled = false;

        // Auto-select if obvious
        if (stringCols.length === 1) geneCol.value = stringCols[0];
        if (numericCols.length === 1) metricCol.value = numericCols[0];

        // Try common gene column names
        const geneNames = ['gene', 'gene_name', 'genename', 'symbol', 'gene_symbol',
                           'genesymbol', 'id', 'gene_id', 'name', 'genes'];
        for (const name of geneNames) {
            const match = fields.find(f => f.toLowerCase() === name);
            if (match) { geneCol.value = match; break; }
        }

        // Try common metric column names
        const metricNames = ['log2foldchange', 'log2fc', 'logfc', 'lfc',
                            'stat', 't', 'score', 'neg.log10.fdr', 'rank'];
        for (const name of metricNames) {
            const match = fields.find(f => f.toLowerCase().replace(/[._]/g, '') === name.replace(/[._]/g, ''));
            if (match) { metricCol.value = match; break; }
        }

        this.checkReady();
    }

    checkReady() {
        const geneCol = document.getElementById('geneColumn').value;
        const metricCol = document.getElementById('metricColumn').value;
        const hasData = this.rawData && geneCol && metricCol;
        const hasGeneSets = this.getActiveGeneSets() !== null;
        document.getElementById('runBtn').disabled = !(hasData && hasGeneSets);
    }

    getActiveGeneSets() {
        const all = {};
        const collections = {
            checkHallmark: 'hallmark',
            checkC2: 'c2',
            checkC5: 'c5'
        };
        for (const [checkId, id] of Object.entries(collections)) {
            if (document.getElementById(checkId).checked && this.geneSets[id]) {
                Object.assign(all, this.geneSets[id]);
            }
        }
        Object.assign(all, this.customGeneSets);
        return Object.keys(all).length > 0 ? all : null;
    }

    // --------------------------------------------------------
    // Build Ranked List
    // --------------------------------------------------------
    buildRankedList() {
        const geneCol = document.getElementById('geneColumn').value;
        const metricCol = document.getElementById('metricColumn').value;

        let pairs = [];
        for (const row of this.rawData) {
            const gene = String(row[geneCol] || '').trim().toUpperCase();
            const metric = parseFloat(row[metricCol]);
            if (gene && gene !== 'NAN' && gene !== 'NA' && gene !== '' && !isNaN(metric)) {
                pairs.push({ gene, metric });
            }
        }

        // Remove duplicates (keep first)
        const seen = new Set();
        pairs = pairs.filter(p => {
            if (seen.has(p.gene)) return false;
            seen.add(p.gene);
            return true;
        });

        // Sort descending by metric, ties broken randomly
        pairs.sort((a, b) => {
            if (b.metric !== a.metric) return b.metric - a.metric;
            return Math.random() - 0.5;
        });

        this.rankedList = {
            genes: pairs.map(p => p.gene),
            metrics: pairs.map(p => p.metric)
        };
    }

    // --------------------------------------------------------
    // GSEA Execution
    // --------------------------------------------------------
    runGSEA() {
        this.buildRankedList();

        if (this.rankedList.genes.length === 0) {
            this.showStatus('runStatus', 'error', 'No valid gene-metric pairs found. Check column selection.');
            return;
        }

        const geneSets = this.getActiveGeneSets();
        if (!geneSets) {
            this.showStatus('runStatus', 'error', 'No gene sets loaded.');
            return;
        }

        // Read settings
        this.readSettings();

        // Show progress UI
        document.getElementById('runBtn').style.display = 'none';
        document.getElementById('cancelBtn').style.display = '';
        const container = document.getElementById('progressContainer');
        container.classList.add('active');
        document.getElementById('progressBar').style.width = '0%';
        document.getElementById('progressText').textContent = 'Starting...';
        this.hideStatus('runStatus');

        // Warn about large jobs
        const nSets = Object.keys(geneSets).length;
        if (nSets > 5000) {
            document.getElementById('progressText').textContent =
                `Processing ${nSets} gene sets — this may take several minutes...`;
        }

        this.analysisDate = new Date();

        this.worker.postMessage({
            type: 'run',
            rankedGenes: this.rankedList.genes,
            rankedMetrics: this.rankedList.metrics,
            geneSets: geneSets,
            settings: {
                permutations: this.settings.permutations,
                minSize: this.settings.minSize,
                maxSize: this.settings.maxSize,
                weightP: this.settings.weightP
            }
        });
    }

    cancelAnalysis() {
        this.createWorker();
        document.getElementById('runBtn').style.display = '';
        document.getElementById('cancelBtn').style.display = 'none';
        document.getElementById('progressContainer').classList.remove('active');
        this.showStatus('runStatus', 'warning', 'Analysis cancelled.');
    }

    handleWorkerMessage(e) {
        const data = e.data;

        if (data.type === 'progress') {
            document.getElementById('progressBar').style.width = data.percent + '%';
            document.getElementById('progressText').textContent = data.text;
        }

        if (data.type === 'complete') {
            this.results = data.results;
            document.getElementById('runBtn').style.display = '';
            document.getElementById('cancelBtn').style.display = 'none';
            document.getElementById('progressContainer').classList.remove('active');

            const nSig = this.results.filter(r => r.fdr < 0.25).length;
            this.showStatus('runStatus', 'success',
                `Done! ${this.results.length} gene sets tested, ${nSig} significant (FDR < 0.25)`);

            this.displayResults();
        }

        if (data.type === 'error') {
            document.getElementById('runBtn').style.display = '';
            document.getElementById('cancelBtn').style.display = 'none';
            document.getElementById('progressContainer').classList.remove('active');
            this.showStatus('runStatus', 'error', data.message);
        }
    }

    handleWorkerError(e) {
        document.getElementById('runBtn').style.display = '';
        document.getElementById('cancelBtn').style.display = 'none';
        document.getElementById('progressContainer').classList.remove('active');
        this.showStatus('runStatus', 'error', 'Worker error: ' + (e.message || 'Unknown error'));
    }

    // --------------------------------------------------------
    // Display Results
    // --------------------------------------------------------
    displayResults() {
        // Show results panels, hide empty states
        document.getElementById('overviewEmpty').style.display = 'none';
        document.getElementById('overviewResults').style.display = '';
        document.getElementById('enrichmentEmpty').style.display = 'none';
        document.getElementById('enrichmentResults').style.display = '';
        document.getElementById('tableEmpty').style.display = 'none';
        document.getElementById('tableResults').style.display = '';
        document.getElementById('methodsEmpty').style.display = 'none';
        document.getElementById('methodsResults').style.display = '';
        document.getElementById('settingsCard').style.display = '';

        // Populate gene set selector
        this.populateGeneSetSelector();

        // Render all
        this.renderBubblePlot();
        this.renderRankedPlot();
        this.filterAndRenderTable();
        this.generateMethods();

        // Auto-show first gene set in ES plot
        if (this.results.length > 0) {
            const topSet = this.results[0].name;
            document.getElementById('geneSetSelector').value = topSet;
            this.renderESPlot(topSet);
        }
    }

    populateGeneSetSelector() {
        const sel = document.getElementById('geneSetSelector');
        sel.innerHTML = '<option value="">Select a gene set...</option>';

        // Group by positive/negative NES
        const pos = this.results.filter(r => r.nes > 0).sort((a, b) => b.nes - a.nes);
        const neg = this.results.filter(r => r.nes <= 0).sort((a, b) => a.nes - b.nes);

        if (pos.length > 0) {
            const group = document.createElement('optgroup');
            group.label = `Upregulated (${pos.length})`;
            for (const r of pos) {
                const opt = new Option(
                    `${this.cleanName(r.name)} (NES: ${r.nes.toFixed(2)}, FDR: ${this.formatPval(r.fdr)})`,
                    r.name
                );
                group.appendChild(opt);
            }
            sel.appendChild(group);
        }

        if (neg.length > 0) {
            const group = document.createElement('optgroup');
            group.label = `Downregulated (${neg.length})`;
            for (const r of neg) {
                const opt = new Option(
                    `${this.cleanName(r.name)} (NES: ${r.nes.toFixed(2)}, FDR: ${this.formatPval(r.fdr)})`,
                    r.name
                );
                group.appendChild(opt);
            }
            sel.appendChild(group);
        }
    }

    // --------------------------------------------------------
    // Bubble Plot
    // --------------------------------------------------------
    renderBubblePlot() {
        this.readSettings();
        const fdrThresh = parseFloat(this.settings.fdrDisplayThreshold);
        const topN = this.settings.topN;

        let filtered = this.results;
        if (fdrThresh < 1) {
            filtered = filtered.filter(r => r.fdr < fdrThresh);
        }
        // Sort by |NES| and take top N
        const top = filtered
            .sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes))
            .slice(0, topN);

        if (top.length === 0) {
            Plotly.newPlot('bubblePlot', [], {
                annotations: [{
                    text: 'No gene sets pass the current FDR threshold',
                    xref: 'paper', yref: 'paper', x: 0.5, y: 0.5,
                    showarrow: false, font: { size: 14, color: '#666' }
                }],
                height: 200
            });
            return;
        }

        // Sort by NES for display (most negative at bottom, most positive at top)
        top.sort((a, b) => a.nes - b.nes);

        // FDR color scale: use a proper diverging scale with good contrast
        const fdrColors = top.map(r => {
            const fdr = Math.max(r.fdr, 1e-10);
            const logFdr = -Math.log10(fdr);
            return logFdr;
        });
        const maxLogFdr = Math.max(...fdrColors, 1.5);

        // Dot + lollipop style: horizontal lines from 0 to NES
        const shapes = top.map((r, i) => ({
            type: 'line',
            x0: 0, x1: r.nes,
            y0: i, y1: i,
            line: { color: '#d1d5db', width: 1.5 }
        }));
        // Add vertical zero line
        shapes.push({
            type: 'line', x0: 0, x1: 0,
            yref: 'paper', y0: 0, y1: 1,
            line: { color: '#9ca3af', width: 1.5 }
        });

        const trace = {
            x: top.map(r => r.nes),
            y: top.map((_, i) => i),
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: top.map(r => Math.max(10, Math.min(30, Math.sqrt(r.size) * 2))),
                color: fdrColors,
                colorscale: [
                    [0, '#d4d4d8'],
                    [0.2, '#93c5fd'],
                    [0.4, '#3b82f6'],
                    [0.6, '#f97316'],
                    [0.8, '#ef4444'],
                    [1.0, '#991b1b']
                ],
                cmin: 0,
                cmax: maxLogFdr,
                showscale: true,
                colorbar: {
                    title: { text: 'FDR q-value', font: { size: 11, family: 'Open Sans' } },
                    thickness: 12,
                    len: 0.5,
                    y: 0.5,
                    tickvals: [0, -Math.log10(0.25), -Math.log10(0.1), -Math.log10(0.05), -Math.log10(0.01)].filter(v => v <= maxLogFdr),
                    ticktext: ['1', '0.25', '0.1', '0.05', '0.01'].slice(0, [0, -Math.log10(0.25), -Math.log10(0.1), -Math.log10(0.05), -Math.log10(0.01)].filter(v => v <= maxLogFdr).length),
                    tickfont: { size: 10 },
                    outlinewidth: 0
                },
                line: { width: 1.5, color: 'white' }
            },
            text: top.map(r =>
                `<b>${this.cleanName(r.name)}</b><br>` +
                `NES: ${r.nes.toFixed(3)}<br>` +
                `FDR: ${this.formatPval(r.fdr)}<br>` +
                `p-value: ${this.formatPval(r.pvalue)}<br>` +
                `Size: ${r.size} genes`
            ),
            hoverinfo: 'text'
        };

        const layout = {
            xaxis: {
                title: { text: 'Normalized Enrichment Score (NES)', font: { size: 12 } },
                zeroline: false,
                gridcolor: '#f0f0f0',
                side: 'bottom'
            },
            yaxis: {
                tickvals: top.map((_, i) => i),
                ticktext: top.map(r => this.cleanName(r.name)),
                tickfont: { size: 11 },
                automargin: true,
                gridwidth: 0,
                showgrid: false
            },
            height: Math.max(420, top.length * 26 + 100),
            margin: { l: 10, r: 30, t: 20, b: 55 },
            font: { family: 'Open Sans, sans-serif' },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            shapes: shapes
        };

        Plotly.newPlot('bubblePlot', [trace], layout, {
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d']
        });

        // Click to show ES plot
        document.getElementById('bubblePlot').on('plotly_click', (data) => {
            if (data.points && data.points.length > 0) {
                const idx = data.points[0].pointIndex;
                const name = top[idx].name;
                document.getElementById('geneSetSelector').value = name;
                this.renderESPlot(name);
                this.showTab('enrichment');
            }
        });
    }

    // --------------------------------------------------------
    // Ranked List Plot
    // --------------------------------------------------------
    renderRankedPlot() {
        const N = this.rankedList.genes.length;
        const metrics = this.rankedList.metrics;

        // Subsample for performance
        const maxPts = 3000;
        const step = N > maxPts ? Math.ceil(N / maxPts) : 1;
        const xVals = [];
        const yVals = [];
        for (let i = 0; i < N; i += step) {
            xVals.push(i);
            yVals.push(metrics[i]);
        }

        // Filled area with red/blue coloring
        // Split into positive and negative traces for dual color
        const posX = [], posY = [], negX = [], negY = [];
        for (let i = 0; i < xVals.length; i++) {
            if (yVals[i] >= 0) {
                posX.push(xVals[i]);
                posY.push(yVals[i]);
                negX.push(xVals[i]);
                negY.push(0);
            } else {
                posX.push(xVals[i]);
                posY.push(0);
                negX.push(xVals[i]);
                negY.push(yVals[i]);
            }
        }

        const posTrace = {
            x: posX, y: posY,
            type: 'scatter', mode: 'lines',
            fill: 'tozeroy',
            fillcolor: 'rgba(220, 38, 38, 0.35)',
            line: { color: 'rgba(220, 38, 38, 0.6)', width: 0.5 },
            showlegend: false, hoverinfo: 'skip'
        };

        const negTrace = {
            x: negX, y: negY,
            type: 'scatter', mode: 'lines',
            fill: 'tozeroy',
            fillcolor: 'rgba(37, 99, 235, 0.35)',
            line: { color: 'rgba(37, 99, 235, 0.6)', width: 0.5 },
            showlegend: false, hoverinfo: 'skip'
        };

        const metricLabel = document.getElementById('metricColumn').value || 'Ranking Metric';

        const layout = {
            xaxis: {
                title: { text: 'Gene Rank', font: { size: 11 } },
                showgrid: false
            },
            yaxis: {
                title: { text: metricLabel, font: { size: 11 } },
                zeroline: true,
                zerolinewidth: 1.5,
                zerolinecolor: '#333',
                gridcolor: '#e5e5e5'
            },
            height: 220,
            margin: { l: 65, r: 20, t: 15, b: 45 },
            font: { family: 'Open Sans, sans-serif' },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            shapes: [{
                type: 'rect', xref: 'paper', yref: 'paper',
                x0: 0, x1: 1, y0: 0, y1: 1,
                line: { color: '#333', width: 1.5 },
                fillcolor: 'rgba(0,0,0,0)'
            }]
        };

        Plotly.newPlot('rankedPlot', [posTrace, negTrace], layout, {
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d']
        });
    }

    // --------------------------------------------------------
    // Enrichment Score Plot — Classic GSEA style
    // 3 panels: ES curve (top), hit markers (middle), ranked metric (bottom)
    // --------------------------------------------------------
    renderESPlot(geneSetName) {
        const result = this.results.find(r => r.name === geneSetName);
        if (!result) return;

        const N = this.rankedList.genes.length;
        const metrics = this.rankedList.metrics;

        // Subsample for performance if very large
        const maxPts = 4000;
        const step = N > maxPts ? Math.ceil(N / maxPts) : 1;
        const xSampled = [];
        const esSampled = [];
        const metSampled = [];
        for (let i = 0; i < N; i += step) {
            xSampled.push(i);
            esSampled.push(result.runningES[i]);
            metSampled.push(metrics[i]);
        }
        // Always include the last point
        if (xSampled[xSampled.length - 1] !== N - 1) {
            xSampled.push(N - 1);
            esSampled.push(result.runningES[N - 1]);
            metSampled.push(metrics[N - 1]);
        }

        // Find zero-crossing index in ranked list
        let zeroCross = -1;
        for (let i = 0; i < N - 1; i++) {
            if ((metrics[i] >= 0 && metrics[i + 1] < 0) || (metrics[i] > 0 && metrics[i + 1] <= 0)) {
                zeroCross = i;
                break;
            }
        }

        // ---- Panel 1: Running Enrichment Score ----
        // Green line with fill toward zero
        const esLine = {
            x: xSampled,
            y: esSampled,
            type: 'scatter',
            mode: 'lines',
            line: { color: '#15a04a', width: 2.5 },
            fill: 'tozeroy',
            fillcolor: 'rgba(21, 160, 74, 0.12)',
            showlegend: false,
            xaxis: 'x',
            yaxis: 'y',
            hovertemplate: 'Rank: %{x}<br>ES: %{y:.4f}<extra></extra>'
        };

        const esZero = {
            x: [0, N - 1], y: [0, 0],
            type: 'scatter', mode: 'lines',
            line: { color: '#999', width: 1 },
            showlegend: false, xaxis: 'x', yaxis: 'y', hoverinfo: 'skip'
        };

        // ---- Panel 2: Hit markers ----
        // Each hit is a vertical line spanning the full height of the rug band
        const hitShapes = result.hits.map(idx => ({
            type: 'line',
            x0: idx, x1: idx,
            y0: 0, y1: 1,
            xref: 'x2', yref: 'y2',
            line: { color: '#111', width: 1.2 }
        }));

        // Invisible trace to anchor the subplot
        const rugAnchor = {
            x: [0, N - 1], y: [0.5, 0.5],
            type: 'scatter', mode: 'lines',
            line: { color: 'rgba(0,0,0,0)', width: 0 },
            showlegend: false, xaxis: 'x2', yaxis: 'y2', hoverinfo: 'skip'
        };

        // ---- Panel 3: Ranked list metric with red-to-blue gradient ----
        // Create a filled area plot with the metric values
        const metricLine = {
            x: xSampled,
            y: metSampled,
            type: 'scatter',
            mode: 'lines',
            line: { color: '#999', width: 0.8 },
            fill: 'tozeroy',
            fillcolor: 'rgba(180,180,180,0.35)',
            showlegend: false,
            xaxis: 'x3',
            yaxis: 'y3',
            hoverinfo: 'skip'
        };

        const metZero = {
            x: [0, N - 1], y: [0, 0],
            type: 'scatter', mode: 'lines',
            line: { color: '#333', width: 0.8 },
            showlegend: false, xaxis: 'x3', yaxis: 'y3', hoverinfo: 'skip'
        };

        // ---- Layout ----
        // Panel borders via shapes
        const panelShapes = [
            // ES panel border
            { type: 'rect', xref: 'paper', yref: 'paper', x0: 0, x1: 1, y0: 0.44, y1: 1,
              line: { color: '#333', width: 1.5 }, fillcolor: 'rgba(0,0,0,0)' },
            // Hit panel border
            { type: 'rect', xref: 'paper', yref: 'paper', x0: 0, x1: 1, y0: 0.32, y1: 0.44,
              line: { color: '#333', width: 1.5 }, fillcolor: 'rgba(0,0,0,0)' },
            // Metric panel border
            { type: 'rect', xref: 'paper', yref: 'paper', x0: 0, x1: 1, y0: 0, y1: 0.32,
              line: { color: '#333', width: 1.5 }, fillcolor: 'rgba(0,0,0,0)' },
            ...hitShapes
        ];

        // Red-blue gradient bar for positively/negatively correlated region
        if (zeroCross >= 0) {
            const crossFrac = zeroCross / (N - 1);
            // Red band (positive)
            panelShapes.push({
                type: 'rect', xref: 'x3', yref: 'paper',
                x0: 0, x1: zeroCross,
                y0: 0.25, y1: 0.32,
                fillcolor: 'rgba(239, 68, 68, 0.2)', line: { width: 0 }
            });
            // Blue band (negative)
            panelShapes.push({
                type: 'rect', xref: 'x3', yref: 'paper',
                x0: zeroCross, x1: N - 1,
                y0: 0.25, y1: 0.32,
                fillcolor: 'rgba(59, 130, 246, 0.2)', line: { width: 0 }
            });
        }

        const metricLabel = document.getElementById('metricColumn').value || 'Ranking metric';

        // Annotations
        const annotations = [
            // Title
            {
                text: `<b>Enrichment plot: ${geneSetName.replace(/_/g, ' ')}</b>`,
                xref: 'paper', yref: 'paper', x: 0.5, y: 1.08,
                showarrow: false, font: { size: 13, family: 'Open Sans' },
                xanchor: 'center'
            },
            // ES panel label
            {
                text: 'Enrichment score (ES)',
                xref: 'paper', yref: 'paper', x: -0.01, y: 0.72,
                showarrow: false, font: { size: 11, color: '#333' },
                textangle: -90, xanchor: 'right'
            },
            // Hits label
            {
                text: 'Hits',
                xref: 'paper', yref: 'paper', x: -0.01, y: 0.38,
                showarrow: false, font: { size: 10, color: '#333' },
                textangle: -90, xanchor: 'right'
            },
            // Ranked metric label
            {
                text: `Ranked list metric<br>(${metricLabel})`,
                xref: 'paper', yref: 'paper', x: -0.01, y: 0.16,
                showarrow: false, font: { size: 10, color: '#333' },
                textangle: -90, xanchor: 'right'
            },
            // Stats annotation (top right of ES panel)
            {
                text: `NES = ${result.nes.toFixed(2)}<br>FDR = ${this.formatPval(result.fdr)}<br>p = ${this.formatPval(result.pvalue)}<br>Size = ${result.size}`,
                xref: 'paper', yref: 'paper',
                x: result.nes >= 0 ? 0.98 : 0.02,
                y: 0.97,
                showarrow: false,
                font: { size: 10, family: 'Roboto Mono, monospace', color: '#333' },
                align: result.nes >= 0 ? 'right' : 'left',
                xanchor: result.nes >= 0 ? 'right' : 'left',
                yanchor: 'top',
                bgcolor: 'rgba(255,255,255,0.85)',
                bordercolor: '#ccc',
                borderwidth: 1,
                borderpad: 4
            }
        ];

        // Zero cross annotation
        if (zeroCross >= 0) {
            annotations.push({
                text: `Zero cross at ${zeroCross}`,
                x: zeroCross, y: 0,
                xref: 'x3', yref: 'y3',
                showarrow: true, arrowhead: 0, arrowsize: 0.8, arrowwidth: 1,
                ax: 0, ay: -25,
                font: { size: 9, color: '#666' }
            });
        }

        const layout = {
            xaxis: {
                range: [0, N - 1], showticklabels: false, showgrid: false, zeroline: false,
                domain: [0.08, 1]
            },
            yaxis: {
                domain: [0.44, 0.98],
                gridcolor: '#e5e5e5', gridwidth: 1,
                zeroline: false,
                tickfont: { size: 10 }
            },
            xaxis2: {
                range: [0, N - 1], showticklabels: false, showgrid: false, zeroline: false,
                domain: [0.08, 1], matches: 'x'
            },
            yaxis2: {
                domain: [0.32, 0.44],
                range: [0, 1],
                showticklabels: false, showgrid: false, zeroline: false, fixedrange: true
            },
            xaxis3: {
                range: [0, N - 1],
                title: { text: 'Rank in Ordered Dataset', font: { size: 11 }, standoff: 5 },
                domain: [0.08, 1], matches: 'x',
                showgrid: false,
                tickfont: { size: 10 }
            },
            yaxis3: {
                domain: [0, 0.32],
                gridcolor: '#e5e5e5', gridwidth: 1,
                zeroline: false,
                tickfont: { size: 10 }
            },
            height: 600,
            margin: { l: 80, r: 20, t: 55, b: 50 },
            font: { family: 'Open Sans, sans-serif' },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            shapes: panelShapes,
            annotations: annotations
        };

        Plotly.newPlot('esPlot',
            [esLine, esZero, rugAnchor, metricLine, metZero],
            layout,
            { responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d'] }
        );
    }

    // --------------------------------------------------------
    // Results Table
    // --------------------------------------------------------
    filterAndRenderTable() {
        if (!this.results) return;

        const query = document.getElementById('tableFilter').value.toLowerCase();
        const fdrVal = document.getElementById('fdrFilter').value;
        const dirVal = document.getElementById('directionFilter').value;

        let filtered = this.results.slice();

        if (query) {
            filtered = filtered.filter(r => r.name.toLowerCase().includes(query));
        }
        if (fdrVal !== 'all') {
            const thresh = parseFloat(fdrVal);
            filtered = filtered.filter(r => r.fdr < thresh);
        }
        if (dirVal === 'up') {
            filtered = filtered.filter(r => r.nes > 0);
        } else if (dirVal === 'down') {
            filtered = filtered.filter(r => r.nes < 0);
        }

        // Sort
        filtered.sort((a, b) => {
            let valA = a[this.sortCol];
            let valB = b[this.sortCol];
            if (this.sortCol === 'name') {
                valA = valA.toLowerCase();
                valB = valB.toLowerCase();
            }
            if (this.sortCol === 'leadingEdge') {
                valA = (valA || []).length;
                valB = (valB || []).length;
            }
            if (valA < valB) return this.sortAsc ? -1 : 1;
            if (valA > valB) return this.sortAsc ? 1 : -1;
            return 0;
        });

        // Update count
        document.getElementById('resultCount').textContent =
            `${filtered.length} of ${this.results.length} gene sets`;

        // Render rows
        const tbody = document.getElementById('resultsBody');
        tbody.innerHTML = '';

        for (const r of filtered) {
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td title="${r.name}">${this.cleanName(r.name)}</td>
                <td>${r.size}</td>
                <td class="${r.es >= 0 ? 'positive' : 'negative'}">${r.es.toFixed(4)}</td>
                <td class="${r.nes >= 0 ? 'positive' : 'negative'}"><strong>${r.nes.toFixed(3)}</strong></td>
                <td>${this.formatPval(r.pvalue)}</td>
                <td>${this.formatPval(r.fdr)}</td>
                <td class="leading-edge-cell" title="${(r.leadingEdge || []).join(', ')}">${(r.leadingEdge || []).length} genes</td>
            `;
            tr.addEventListener('click', () => {
                document.getElementById('geneSetSelector').value = r.name;
                this.renderESPlot(r.name);
                this.showTab('enrichment');
            });
            tbody.appendChild(tr);
        }

        // Update sort indicators
        document.querySelectorAll('#resultsTable thead th').forEach(th => {
            const arrow = th.querySelector('.sort-arrow');
            if (th.dataset.col === this.sortCol) {
                th.classList.add('sorted');
                arrow.textContent = this.sortAsc ? '\u25B2' : '\u25BC';
            } else {
                th.classList.remove('sorted');
                arrow.textContent = '\u25B2';
            }
        });
    }

    sortTable(col) {
        if (this.sortCol === col) {
            this.sortAsc = !this.sortAsc;
        } else {
            this.sortCol = col;
            this.sortAsc = col === 'name'; // ascending for name, descending for numbers
        }
        this.filterAndRenderTable();
    }

    // --------------------------------------------------------
    // Methods Section
    // --------------------------------------------------------
    generateMethods() {
        const metricCol = document.getElementById('metricColumn').value;
        const collections = [];
        if (document.getElementById('checkHallmark').checked) collections.push('Hallmark (H)');
        if (document.getElementById('checkC2').checked) collections.push('Curated (C2)');
        if (document.getElementById('checkC5').checked) collections.push('GO (C5)');
        if (Object.keys(this.customGeneSets).length > 0) collections.push('Custom');

        const nSets = this.results.length;
        const nSig = this.results.filter(r => r.fdr < 0.25).length;
        const date = this.analysisDate
            ? this.analysisDate.toISOString().split('T')[0]
            : new Date().toISOString().split('T')[0];

        const text = `Preranked Gene Set Enrichment Analysis (GSEA) was performed using a ` +
            `client-side JavaScript implementation of the GSEA algorithm ` +
            `(Subramanian et al., PNAS 2005; Mootha et al., Nature Genetics 2003). ` +
            `Genes were ranked by ${metricCol} in descending order. ` +
            `Enrichment scores were computed using a weighted running-sum statistic ` +
            `(weight parameter p = ${this.settings.weightP}). ` +
            `Gene set collections from MSigDB v2023.2 ` +
            `(${collections.join(', ')}) were used. ` +
            `${nSets} gene sets passed the size filter ` +
            `(minimum: ${this.settings.minSize}, maximum: ${this.settings.maxSize} genes). ` +
            `Statistical significance was assessed using ${this.settings.permutations} ` +
            `permutations of gene set membership labels. ` +
            `Normalized enrichment scores (NES) were computed by dividing the observed ` +
            `enrichment score by the mean of permuted enrichment scores of the same sign. ` +
            `P-values were adjusted for multiple testing using the Benjamini-Hochberg ` +
            `procedure. ${nSig} gene sets were significant at FDR < 0.25. ` +
            `Analysis was performed on ${date} using the GSEA Web Tool ` +
            `(Wermeling Lab, Karolinska Institutet).`;

        document.getElementById('methodsText').textContent = text;
    }

    copyMethods() {
        const text = document.getElementById('methodsText').textContent;
        navigator.clipboard.writeText(text).then(() => {
            const btn = document.getElementById('copyMethodsBtn');
            btn.textContent = 'Copied!';
            setTimeout(() => { btn.textContent = 'Copy'; }, 2000);
        });
    }

    // --------------------------------------------------------
    // Export
    // --------------------------------------------------------
    exportPlot(plotId, format) {
        const plotEl = document.getElementById(plotId);
        if (!plotEl || !plotEl.data) return;

        const scale = this.settings.exportScale;
        const bgColor = this.settings.transparentBg ? 'rgba(0,0,0,0)' : 'white';

        if (format === 'svg') {
            Plotly.toImage(plotEl, { format: 'svg', width: plotEl.offsetWidth, height: plotEl.offsetHeight })
                .then(url => {
                    // Convert data URL to blob
                    const svgData = atob(url.split(',')[1]);
                    const blob = new Blob([svgData], { type: 'image/svg+xml' });
                    this.downloadBlob(blob, `${plotId}.svg`);
                });
        } else {
            Plotly.toImage(plotEl, {
                format: 'png',
                width: plotEl.offsetWidth * scale,
                height: plotEl.offsetHeight * scale,
                scale: 1
            }).then(url => {
                const link = document.createElement('a');
                link.download = `${plotId}.png`;
                link.href = url;
                link.click();
            });
        }
    }

    downloadCSV() {
        if (!this.results) return;

        const headers = ['Gene Set', 'Size', 'ES', 'NES', 'p-value', 'FDR', 'Leading Edge'];
        const rows = this.results.map(r => [
            r.name,
            r.size,
            r.es.toFixed(6),
            r.nes.toFixed(6),
            r.pvalue.toExponential(4),
            r.fdr.toExponential(4),
            '"' + (r.leadingEdge || []).join(', ') + '"'
        ]);

        const csv = [headers.join(','), ...rows.map(r => r.join(','))].join('\n');
        const blob = new Blob([csv], { type: 'text/csv' });
        this.downloadBlob(blob, 'gsea_results.csv');
    }

    downloadBlob(blob, filename) {
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
    }

    // --------------------------------------------------------
    // Example Data
    // --------------------------------------------------------
    loadExampleData() {
        // Generate a simple example: gene names with simulated log2FC values
        this.showStatus('uploadStatus', 'info', 'Loading example data...');

        // Create example data with some known Hallmark pathway genes upregulated
        const exampleGenes = this.generateExampleData();
        this.rawData = exampleGenes;
        this.populateColumnDropdowns(['Gene', 'log2FoldChange', 'pvalue']);
        document.getElementById('geneColumn').value = 'Gene';
        document.getElementById('metricColumn').value = 'log2FoldChange';
        this.showStatus('uploadStatus', 'success',
            `Loaded example data: ${exampleGenes.length} genes with simulated log2FC values`);
        this.checkReady();
    }

    generateExampleData() {
        // Use Hallmark gene sets to create biologically meaningful example data
        const hallmark = this.geneSets['hallmark'];
        if (!hallmark) return [];

        // Collect all unique genes from Hallmark sets
        const allGenes = new Set();
        for (const genes of Object.values(hallmark)) {
            for (const g of genes) allGenes.add(g.toUpperCase());
        }

        // Pick sets to be "enriched" (upregulated) and "depleted" (downregulated)
        const setNames = Object.keys(hallmark);
        const upSets = ['HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_INFLAMMATORY_RESPONSE',
                        'HALLMARK_INTERFERON_GAMMA_RESPONSE'];
        const downSets = ['HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'HALLMARK_FATTY_ACID_METABOLISM'];

        const upGenes = new Set();
        const downGenes = new Set();
        for (const name of upSets) {
            if (hallmark[name]) {
                for (const g of hallmark[name]) upGenes.add(g.toUpperCase());
            }
        }
        for (const name of downSets) {
            if (hallmark[name]) {
                for (const g of hallmark[name]) downGenes.add(g.toUpperCase());
            }
        }

        // Generate data
        const data = [];
        for (const gene of allGenes) {
            let log2fc;
            if (upGenes.has(gene)) {
                log2fc = 0.5 + Math.random() * 3;
            } else if (downGenes.has(gene)) {
                log2fc = -(0.5 + Math.random() * 3);
            } else {
                log2fc = (Math.random() - 0.5) * 1.5;
            }
            const pval = Math.pow(10, -Math.abs(log2fc) * (1 + Math.random() * 2));
            data.push({
                Gene: gene,
                log2FoldChange: Math.round(log2fc * 10000) / 10000,
                pvalue: pval
            });
        }

        return data;
    }

    // --------------------------------------------------------
    // Settings
    // --------------------------------------------------------
    readSettings() {
        this.settings.permutations = parseInt(document.getElementById('permutations').value) || 1000;
        this.settings.minSize = parseInt(document.getElementById('minSize').value) || 15;
        this.settings.maxSize = parseInt(document.getElementById('maxSize').value) || 500;
        this.settings.weightP = parseFloat(document.getElementById('weightP').value) || 1;
        this.settings.topN = parseInt(document.getElementById('topN').value) || 20;
        this.settings.fdrDisplayThreshold = parseFloat(document.getElementById('fdrDisplayThreshold').value);
        this.settings.colorScale = document.getElementById('colorScale').value;
        this.settings.exportScale = parseInt(document.getElementById('exportScale').value) || 4;
        this.settings.transparentBg = document.getElementById('transparentBg').checked;
    }

    updateSettings() {
        this.readSettings();
        if (this.results) {
            this.renderBubblePlot();
            this.renderRankedPlot();
            const selected = document.getElementById('geneSetSelector').value;
            if (selected) this.renderESPlot(selected);
        }
    }

    // --------------------------------------------------------
    // Tab Navigation
    // --------------------------------------------------------
    showTab(tabName) {
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.tab === tabName);
        });
        document.querySelectorAll('.tab-panel').forEach(panel => {
            panel.classList.toggle('active', panel.id === 'tab-' + tabName);
        });
    }

    // --------------------------------------------------------
    // Utility
    // --------------------------------------------------------
    showStatus(elementId, type, message) {
        const el = document.getElementById(elementId);
        el.className = `status-box status-${type}`;
        el.textContent = message;
        el.classList.remove('hidden');
    }

    hideStatus(elementId) {
        document.getElementById(elementId).classList.add('hidden');
    }

    cleanName(name) {
        return name
            .replace(/^HALLMARK_/, '')
            .replace(/^GO_/, '')
            .replace(/^KEGG_/, '')
            .replace(/^REACTOME_/, '')
            .replace(/^WP_/, '')
            .replace(/^BIOCARTA_/, '')
            .replace(/_/g, ' ');
    }

    formatPval(val) {
        if (val === undefined || val === null) return 'N/A';
        if (val < 0.001) return val.toExponential(2);
        if (val < 0.01) return val.toFixed(4);
        return val.toFixed(3);
    }
}

// Initialize app
const app = new GSEAApp();
