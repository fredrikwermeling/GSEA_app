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

        // Sort by NES for display (positive at top)
        top.sort((a, b) => a.nes - b.nes);

        const trace = {
            x: top.map(r => r.nes),
            y: top.map(r => this.cleanName(r.name)),
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: top.map(r => Math.max(8, Math.sqrt(r.size) * 2.5)),
                color: top.map(r => -Math.log10(Math.max(r.fdr, 1e-10))),
                colorscale: this.settings.colorScale,
                reversescale: false,
                showscale: true,
                colorbar: {
                    title: { text: '-log10(FDR)', font: { size: 12 } },
                    thickness: 14,
                    len: 0.6
                },
                line: { width: 1, color: '#555' }
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
            title: { text: 'GSEA Enrichment Overview', font: { size: 16 } },
            xaxis: {
                title: 'Normalized Enrichment Score (NES)',
                zeroline: true,
                zerolinewidth: 1,
                zerolinecolor: '#aaa'
            },
            yaxis: {
                automargin: true,
                tickfont: { size: 11 }
            },
            height: Math.max(400, top.length * 28 + 100),
            margin: { l: 10, r: 80, t: 50, b: 60 },
            font: { family: 'Open Sans, sans-serif' },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fafafa',
            shapes: [{
                type: 'line', x0: 0, x1: 0,
                yref: 'paper', y0: 0, y1: 1,
                line: { color: '#999', width: 1, dash: 'dash' }
            }]
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

        // For large lists, subsample for performance
        let xVals, yVals, colors;
        if (N > 5000) {
            const step = Math.ceil(N / 3000);
            xVals = [];
            yVals = [];
            colors = [];
            for (let i = 0; i < N; i += step) {
                xVals.push(i);
                yVals.push(metrics[i]);
                colors.push(metrics[i] >= 0 ? 'rgba(220, 38, 38, 0.6)' : 'rgba(37, 99, 235, 0.6)');
            }
        } else {
            xVals = Array.from({ length: N }, (_, i) => i);
            yVals = metrics;
            colors = metrics.map(v => v >= 0 ? 'rgba(220, 38, 38, 0.6)' : 'rgba(37, 99, 235, 0.6)');
        }

        const trace = {
            x: xVals,
            y: yVals,
            type: 'bar',
            marker: { color: colors },
            hoverinfo: 'skip',
            showlegend: false
        };

        const layout = {
            title: { text: 'Ranked Gene List', font: { size: 14 } },
            xaxis: {
                title: 'Gene Rank',
                showgrid: false
            },
            yaxis: {
                title: document.getElementById('metricColumn').value || 'Ranking Metric',
                zeroline: true,
                zerolinewidth: 1,
                zerolinecolor: '#333'
            },
            height: 250,
            margin: { l: 70, r: 20, t: 40, b: 50 },
            font: { family: 'Open Sans, sans-serif' },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fafafa',
            bargap: 0
        };

        Plotly.newPlot('rankedPlot', [trace], layout, {
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d']
        });
    }

    // --------------------------------------------------------
    // Enrichment Score Plot (3 subplots)
    // --------------------------------------------------------
    renderESPlot(geneSetName) {
        const result = this.results.find(r => r.name === geneSetName);
        if (!result) return;

        const N = this.rankedList.genes.length;
        const xAxis = Array.from({ length: N }, (_, i) => i);
        const metrics = this.rankedList.metrics;

        // Subplot 1: Running ES curve
        const esTrace = {
            x: xAxis,
            y: result.runningES,
            type: 'scatter',
            mode: 'lines',
            line: { color: '#2ca02c', width: 2.5 },
            name: 'Running ES',
            xaxis: 'x',
            yaxis: 'y',
            showlegend: false
        };

        const zeroTrace = {
            x: [0, N - 1],
            y: [0, 0],
            type: 'scatter',
            mode: 'lines',
            line: { color: '#aaa', width: 1, dash: 'dash' },
            showlegend: false,
            xaxis: 'x',
            yaxis: 'y',
            hoverinfo: 'skip'
        };

        // Subplot 2: Hit markers (rug plot)
        const rugX = result.hits;
        const rugY = rugX.map(() => 1);
        const rugTrace = {
            x: rugX,
            y: rugY,
            type: 'scatter',
            mode: 'markers',
            marker: {
                symbol: 'line-ns',
                size: 10,
                color: 'black',
                line: { width: 0.8 }
            },
            showlegend: false,
            xaxis: 'x',
            yaxis: 'y2',
            hoverinfo: 'skip'
        };

        // Subplot 3: Ranked metric values (as filled area for performance)
        const barColors = metrics.map(v =>
            v >= 0 ? 'rgba(220, 38, 38, 0.5)' : 'rgba(37, 99, 235, 0.5)'
        );

        // For large N, use area fill instead of bars for performance
        const metricTrace = {
            x: xAxis,
            y: metrics,
            type: 'scatter',
            mode: 'lines',
            fill: 'tozeroy',
            fillcolor: 'rgba(150,150,150,0.3)',
            line: { width: 0.5, color: '#888' },
            showlegend: false,
            xaxis: 'x',
            yaxis: 'y3',
            hoverinfo: 'skip'
        };

        const zeroTrace3 = {
            x: [0, N - 1],
            y: [0, 0],
            type: 'scatter',
            mode: 'lines',
            line: { color: '#333', width: 0.5 },
            showlegend: false,
            xaxis: 'x',
            yaxis: 'y3',
            hoverinfo: 'skip'
        };

        const layout = {
            grid: {
                rows: 3,
                columns: 1,
                subplots: [['xy'], ['xy2'], ['xy3']],
                roworder: 'top to bottom',
                pattern: 'independent'
            },
            xaxis: { showticklabels: false, showgrid: false, zeroline: false, domain: [0, 1] },
            yaxis: {
                title: { text: 'Enrichment Score', font: { size: 12 } },
                domain: [0.42, 1],
                zeroline: true
            },
            xaxis2: { showticklabels: false, showgrid: false, zeroline: false, domain: [0, 1], matches: 'x' },
            yaxis2: {
                domain: [0.32, 0.40],
                showticklabels: false,
                showgrid: false,
                zeroline: false,
                range: [0, 2]
            },
            xaxis3: {
                title: { text: 'Gene Rank', font: { size: 12 } },
                domain: [0, 1],
                matches: 'x'
            },
            yaxis3: {
                title: { text: document.getElementById('metricColumn').value || 'Metric', font: { size: 11 } },
                domain: [0, 0.30],
                zeroline: true
            },
            height: 550,
            margin: { l: 70, r: 30, t: 60, b: 50 },
            font: { family: 'Open Sans, sans-serif' },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fafafa',
            title: {
                text: this.cleanName(geneSetName) +
                      `<br><span style="font-size:12px; color:#666;">` +
                      `NES = ${result.nes.toFixed(3)} &nbsp; FDR = ${this.formatPval(result.fdr)} &nbsp; ` +
                      `p = ${this.formatPval(result.pvalue)} &nbsp; Size = ${result.size}</span>`,
                font: { size: 14 }
            }
        };

        Plotly.newPlot('esPlot',
            [esTrace, zeroTrace, rugTrace, metricTrace, zeroTrace3],
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
