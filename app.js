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
            transparentBg: false,
            // Figure customization
            fontFamily: 'Open Sans',
            fontSize: 12,
            esLineColor: '#15a04a',
            esLineWidth: 2.5,
            showStatsBox: true,
            showZeroCross: true,
            showCorrelationLabels: true,
            showPanelBorders: true,
            showHitMarkers: true,
            showMetricFill: true
        };

        // Table sort state
        this.sortCol = 'nes';
        this.sortAsc = false;

        // Gene detail table sort state
        this.detailSortCol = 'rank';
        this.detailSortAsc = true;

        // Gene info cache (MyGene.info)
        this.geneInfoCache = {};

        // Active tab for settings panel
        this.activeTab = 'overview';

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
            if (e.target.value) {
                this.renderESPlot(e.target.value);
                this.renderGeneSetInfo(e.target.value);
                this.renderGeneDetailTable(e.target.value);
            }
        });

        // Prev/Next gene set buttons
        document.getElementById('prevGeneSetBtn').addEventListener('click', () => this.navigateGeneSet(-1));
        document.getElementById('nextGeneSetBtn').addEventListener('click', () => this.navigateGeneSet(1));

        // Gene search
        document.getElementById('searchGenesBtn').addEventListener('click', () => this.searchGenes());
        document.getElementById('clearSearchBtn').addEventListener('click', () => this.clearGeneSearch());

        // Top/Bottom N selector
        document.getElementById('topBottomN').addEventListener('change', () => this.renderTopBottomGenes());

        // Gene detail CSV download
        document.getElementById('downloadGeneDetailCSV').addEventListener('click', () => this.downloadGeneDetailCSV());

        // Overlap heatmap
        document.getElementById('updateOverlapBtn').addEventListener('click', () => this.renderOverlapHeatmap());

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

        // Methods card buttons
        document.getElementById('copyMethodsBtn').addEventListener('click', () => this.copyMethods());
        document.getElementById('downloadMethodsBtn').addEventListener('click', () => this.downloadMethods());

        // How to Use popup
        const howToUseLink = document.getElementById('howToUseLink');
        const howToUsePopup = document.getElementById('howToUsePopup');
        const howToUseBackdrop = document.getElementById('howToUseBackdrop');
        const howToUseClose = document.getElementById('howToUseClose');
        howToUseLink.addEventListener('click', () => {
            howToUsePopup.classList.add('open');
            howToUseBackdrop.classList.add('open');
        });
        const closeHowToUse = () => {
            howToUsePopup.classList.remove('open');
            howToUseBackdrop.classList.remove('open');
        };
        howToUseClose.addEventListener('click', closeHowToUse);
        howToUseBackdrop.addEventListener('click', closeHowToUse);

        // Info tooltips — position with JS (fixed) to avoid clipping
        // Use event delegation so dynamically added info-icons also work
        document.addEventListener('mouseenter', (e) => {
            const icon = e.target.closest('.info-icon');
            if (!icon) return;
            const tooltip = icon.querySelector('.info-tooltip');
            if (!tooltip) return;
            const rect = icon.getBoundingClientRect();
            tooltip.style.display = 'block';
            let left = rect.right + 10;
            let top = rect.top - 4;
            if (left + 330 > window.innerWidth) {
                left = rect.left - 340;
            }
            if (top + tooltip.offsetHeight > window.innerHeight - 10) {
                top = window.innerHeight - tooltip.offsetHeight - 10;
            }
            if (top < 10) top = 10;
            tooltip.style.left = left + 'px';
            tooltip.style.top = top + 'px';
        }, true);
        document.addEventListener('mouseleave', (e) => {
            const icon = e.target.closest('.info-icon');
            if (!icon) return;
            const tooltip = icon.querySelector('.info-tooltip');
            if (tooltip) tooltip.style.display = 'none';
        }, true);

        // Real-time settings: update plots on every change
        document.querySelectorAll('.settings-live').forEach(el => {
            const evType = (el.type === 'checkbox' || el.tagName === 'SELECT') ? 'change' : 'input';
            el.addEventListener(evType, () => this.updateSettings());
        });

        // Methods toggle (collapsible)
        document.getElementById('methodsToggle').addEventListener('click', () => {
            const body = document.getElementById('methodsBody');
            const chevron = document.getElementById('methodsChevron');
            if (body.style.display === 'none') {
                body.style.display = '';
                chevron.style.transform = 'rotate(90deg)';
            } else {
                body.style.display = 'none';
                chevron.style.transform = '';
            }
        });

        // Reset app
        document.getElementById('resetBtn').addEventListener('click', () => this.resetApp());

        // Floating settings panel open/close
        document.getElementById('openSettingsBtn').addEventListener('click', () => {
            const panel = document.getElementById('settingsPanel');
            panel.classList.toggle('open');
            this.updateSettingsTabVisibility();
        });
        document.getElementById('closeSettingsBtn').addEventListener('click', () => {
            document.getElementById('settingsPanel').classList.remove('open');
        });

        // Make settings panel draggable
        this.initSettingsDrag();

        // Gene detail table sorting
        document.querySelectorAll('#geneDetailTable thead th').forEach(th => {
            th.addEventListener('click', () => {
                const col = th.dataset.detailCol;
                if (!col) return;
                if (this.detailSortCol === col) {
                    this.detailSortAsc = !this.detailSortAsc;
                } else {
                    this.detailSortCol = col;
                    this.detailSortAsc = col === 'gene'; // ascending for gene name
                }
                const sel = document.getElementById('geneSetSelector').value;
                if (sel) this.renderGeneDetailTable(sel);
            });
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
        document.getElementById('overlapEmpty').style.display = 'none';
        document.getElementById('overlapResults').style.display = '';
        document.getElementById('openSettingsBtn').style.display = '';
        document.getElementById('methodsCard').style.display = '';
        this.updateSettingsTabVisibility();

        // Populate gene set selector
        this.populateGeneSetSelector();

        // Render all
        this.renderBubblePlot();
        this.renderRankedPlot();
        this.filterAndRenderTable();
        this.generateMethods();
        this.renderTopBottomGenes();

        // Auto-show first gene set in ES plot
        if (this.results.length > 0) {
            const topSet = this.results[0].name;
            document.getElementById('geneSetSelector').value = topSet;
            this.renderESPlot(topSet);
            this.renderGeneSetInfo(topSet);
            this.renderGeneDetailTable(topSet);
        }

        // Render overlap heatmap
        this.renderOverlapHeatmap();
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
        const fontFam = this.settings.fontFamily + ', sans-serif';
        const baseFontSize = this.settings.fontSize;

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
                    thickness: 18,
                    len: 0.6,
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
                title: { text: 'Normalized Enrichment Score (NES)', font: { size: baseFontSize, family: fontFam } },
                zeroline: false,
                gridcolor: '#f0f0f0',
                side: 'bottom',
                tickfont: { size: baseFontSize - 2, family: fontFam }
            },
            yaxis: {
                tickvals: top.map((_, i) => i),
                ticktext: top.map(r => this.cleanName(r.name)),
                tickfont: { size: baseFontSize - 1, family: fontFam },
                automargin: true,
                gridwidth: 0,
                showgrid: false,
                range: [-0.8, top.length - 0.2]
            },
            height: Math.max(440, top.length * 26 + 120),
            margin: { l: 10, r: 40, t: 20, b: 65 },
            font: { family: fontFam },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            shapes: shapes
        };

        // Apply custom dimensions
        const dims = this.getPlotDimensions(undefined, layout.height);
        if (dims.width) layout.width = dims.width;
        if (dims.height) layout.height = dims.height;

        Plotly.newPlot('bubblePlot', [trace], layout, {
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            edits: { annotationPosition: true, annotationText: true, axisTitleText: false, titleText: false, legendPosition: true }
        });

        // Click to show ES plot (dot or y-axis label)
        const bubblePlotEl = document.getElementById('bubblePlot');
        bubblePlotEl.on('plotly_click', (data) => {
            if (data.points && data.points.length > 0) {
                const idx = data.points[0].pointIndex;
                const name = top[idx].name;
                document.getElementById('geneSetSelector').value = name;
                this.renderESPlot(name);
                this.renderGeneSetInfo(name);
                this.renderGeneDetailTable(name);
                this.showTab('enrichment');
            }
        });
        // Also handle clicking y-axis labels (tick text)
        bubblePlotEl.addEventListener('click', (e) => {
            const yTick = e.target.closest('.ytick text, .yaxislayer-above text');
            if (yTick) {
                const label = yTick.textContent.trim();
                const match = top.find(r => this.cleanName(r.name) === label);
                if (match) {
                    document.getElementById('geneSetSelector').value = match.name;
                    this.renderESPlot(match.name);
                    this.renderGeneSetInfo(match.name);
                    this.renderGeneDetailTable(match.name);
                    this.showTab('enrichment');
                }
            }
        });
    }

    // --------------------------------------------------------
    // Ranked List Plot — Classic GSEA waterfall style
    // --------------------------------------------------------
    renderRankedPlot() {
        this.readSettings();
        const fontFam = this.settings.fontFamily + ', sans-serif';
        const baseFontSize = this.settings.fontSize;
        const N = this.rankedList.genes.length;
        const metrics = this.rankedList.metrics;

        // Subsample to ~600 bars for the waterfall effect
        const maxBars = 600;
        const step = Math.max(1, Math.ceil(N / maxBars));
        const xVals = [];
        const yVals = [];
        const barColors = [];

        for (let i = 0; i < N; i += step) {
            xVals.push(i);
            yVals.push(metrics[i]);
            // Gradient: red for positive, blue for negative, intensity by magnitude
            if (metrics[i] >= 0) {
                const intensity = Math.min(1, metrics[i] / (Math.abs(metrics[0]) || 1));
                barColors.push(`rgba(${Math.round(180 + 60 * intensity)}, ${Math.round(50 * (1 - intensity))}, ${Math.round(50 * (1 - intensity))}, 0.85)`);
            } else {
                const intensity = Math.min(1, Math.abs(metrics[i]) / (Math.abs(metrics[N - 1]) || 1));
                barColors.push(`rgba(${Math.round(50 * (1 - intensity))}, ${Math.round(60 + 50 * (1 - intensity))}, ${Math.round(180 + 60 * intensity)}, 0.85)`);
            }
        }

        const trace = {
            x: xVals,
            y: yVals,
            type: 'bar',
            marker: {
                color: barColors,
                line: { width: 0 }
            },
            showlegend: false,
            hovertemplate: 'Rank: %{x}<br>Metric: %{y:.3f}<extra></extra>'
        };

        const metricLabel = document.getElementById('metricColumn').value || 'Ranking Metric';

        // Find zero-crossing for annotation
        let zeroCross = -1;
        for (let i = 0; i < N - 1; i++) {
            if ((metrics[i] >= 0 && metrics[i + 1] < 0)) {
                zeroCross = i;
                break;
            }
        }

        const annotations = [];
        if (zeroCross >= 0) {
            annotations.push({
                text: 'Zero cross',
                x: zeroCross, y: 0,
                xref: 'x', yref: 'y',
                showarrow: true, arrowhead: 0, arrowsize: 0.8, arrowwidth: 1,
                ax: 0, ay: -20,
                font: { size: 9, color: '#666' }
            });
        }
        // Positive / Negative labels
        annotations.push({
            text: '<b>Positively correlated</b>',
            xref: 'paper', yref: 'paper', x: 0.02, y: 1.02,
            showarrow: false, font: { size: 9, color: '#dc2626' },
            xanchor: 'left', yanchor: 'bottom'
        });
        annotations.push({
            text: '<b>Negatively correlated</b>',
            xref: 'paper', yref: 'paper', x: 0.98, y: 1.02,
            showarrow: false, font: { size: 9, color: '#2563eb' },
            xanchor: 'right', yanchor: 'bottom'
        });

        // Add range padding so extreme bars aren't clipped at edges
        const yPad = Math.max(Math.abs(metrics[0]), Math.abs(metrics[N - 1])) * 0.06;
        const xPad = N * 0.015;

        const layout = {
            xaxis: {
                title: { text: 'Rank in Ordered Dataset', font: { size: baseFontSize - 1, family: fontFam } },
                showgrid: false,
                tickfont: { size: baseFontSize - 2, family: fontFam },
                range: [-xPad, N - 1 + xPad]
            },
            yaxis: {
                title: { text: 'Ranked list metric (' + metricLabel + ')', font: { size: baseFontSize - 1, family: fontFam } },
                zeroline: true,
                zerolinewidth: 1.5,
                zerolinecolor: '#333',
                gridcolor: '#e5e5e5',
                tickfont: { size: baseFontSize - 2, family: fontFam },
                range: [Math.min(0, metrics[N - 1]) - yPad, Math.max(0, metrics[0]) + yPad]
            },
            height: 250,
            margin: { l: 65, r: 20, t: 25, b: 50 },
            font: { family: fontFam },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            bargap: 0,
            shapes: [{
                type: 'rect', xref: 'paper', yref: 'paper',
                x0: 0, x1: 1, y0: 0, y1: 1,
                line: { color: '#333', width: 1.5 },
                fillcolor: 'rgba(0,0,0,0)'
            }],
            annotations: annotations
        };

        // Apply custom dimensions
        const dims = this.getPlotDimensions(undefined, layout.height);
        if (dims.width) layout.width = dims.width;
        if (dims.height) layout.height = dims.height;

        Plotly.newPlot('rankedPlot', [trace], layout, {
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d'],
            edits: { annotationPosition: true, annotationText: true, axisTitleText: true, titleText: false, legendPosition: true }
        });
    }

    // --------------------------------------------------------
    // Enrichment Score Plot — Classic GSEA style
    // 3 panels: ES curve (top), hit markers (middle), ranked metric (bottom)
    // --------------------------------------------------------
    renderESPlot(geneSetName) {
        const result = this.results.find(r => r.name === geneSetName);
        if (!result) return;
        this.readSettings();

        const N = this.rankedList.genes.length;
        const metrics = this.rankedList.metrics;
        const s = this.settings;
        const fontFam = s.fontFamily + ', sans-serif';
        const baseFontSize = s.fontSize;

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

        // ---- Panel domain proportions (classic GSEA style) ----
        // Top: ES curve (55%), Middle: Hit markers (8%), Bottom: Ranked metric (28%)
        // Gaps between panels to prevent overlap
        const esTop = 0.97;
        const esBot = 0.42;
        const hitTop = 0.39;
        const hitBot = 0.31;
        const metTop = 0.28;
        const metBot = 0.0;
        const xLeft = 0.0;  // Use yaxis titles instead of annotations
        const xRight = 1.0;

        // Range padding so extreme ends aren't clipped at borders
        const xPad = N * 0.015;

        // ---- Panel 1: Running Enrichment Score ----
        const esLine = {
            x: xSampled,
            y: esSampled,
            type: 'scatter',
            mode: 'lines',
            line: { color: s.esLineColor, width: s.esLineWidth },
            fill: 'tozeroy',
            fillcolor: s.esLineColor + '1F',  // 12% opacity hex
            showlegend: false,
            xaxis: 'x',
            yaxis: 'y',
            hovertemplate: 'Rank: %{x}<br>ES: %{y:.4f}<extra></extra>'
        };

        const esZero = {
            x: [0, N - 1], y: [0, 0],
            type: 'scatter', mode: 'lines',
            line: { color: '#aaa', width: 1, dash: 'dot' },
            showlegend: false, xaxis: 'x', yaxis: 'y', hoverinfo: 'skip'
        };

        // ---- Panel 2: Hit markers ----
        const hitShapes = [];
        if (s.showHitMarkers) {
            for (const idx of result.hits) {
                hitShapes.push({
                    type: 'line',
                    x0: idx, x1: idx,
                    y0: 0, y1: 1,
                    xref: 'x2', yref: 'y2',
                    line: { color: '#111', width: 1 }
                });
            }
        }

        // Invisible trace to anchor the subplot
        const rugAnchor = {
            x: [0, N - 1], y: [0.5, 0.5],
            type: 'scatter', mode: 'lines',
            line: { color: 'rgba(0,0,0,0)', width: 0 },
            showlegend: false, xaxis: 'x2', yaxis: 'y2', hoverinfo: 'skip'
        };

        // ---- Panel 3: Ranked list metric ----
        const traces3 = [];
        if (s.showMetricFill) {
            // Separate positive and negative for coloring
            const posX = [], posY = [], negX = [], negY = [];
            for (let i = 0; i < xSampled.length; i++) {
                if (metSampled[i] >= 0) {
                    posX.push(xSampled[i]);
                    posY.push(metSampled[i]);
                } else {
                    negX.push(xSampled[i]);
                    negY.push(metSampled[i]);
                }
            }
            if (posX.length > 0) {
                traces3.push({
                    x: posX, y: posY,
                    type: 'bar',
                    marker: { color: 'rgba(220, 38, 38, 0.5)', line: { width: 0 } },
                    showlegend: false, xaxis: 'x3', yaxis: 'y3', hoverinfo: 'skip'
                });
            }
            if (negX.length > 0) {
                traces3.push({
                    x: negX, y: negY,
                    type: 'bar',
                    marker: { color: 'rgba(37, 99, 235, 0.5)', line: { width: 0 } },
                    showlegend: false, xaxis: 'x3', yaxis: 'y3', hoverinfo: 'skip'
                });
            }
        } else {
            traces3.push({
                x: xSampled, y: metSampled,
                type: 'scatter', mode: 'lines',
                line: { color: '#999', width: 0.8 },
                fill: 'tozeroy',
                fillcolor: 'rgba(180,180,180,0.3)',
                showlegend: false, xaxis: 'x3', yaxis: 'y3', hoverinfo: 'skip'
            });
        }

        const metZero = {
            x: [0, N - 1], y: [0, 0],
            type: 'scatter', mode: 'lines',
            line: { color: '#333', width: 0.8 },
            showlegend: false, xaxis: 'x3', yaxis: 'y3', hoverinfo: 'skip'
        };

        // ---- Shapes ----
        const panelShapes = [...hitShapes];

        // Panel borders
        if (s.showPanelBorders) {
            panelShapes.push(
                { type: 'rect', xref: 'paper', yref: 'paper', x0: xLeft, x1: xRight, y0: esBot, y1: esTop,
                  line: { color: '#333', width: 1 }, fillcolor: 'rgba(0,0,0,0)' },
                { type: 'rect', xref: 'paper', yref: 'paper', x0: xLeft, x1: xRight, y0: hitBot, y1: hitTop,
                  line: { color: '#333', width: 1 }, fillcolor: 'rgba(0,0,0,0)' },
                // Metric panel: top + left + right lines only (no bottom — avoids overlap with x-axis)
                { type: 'line', xref: 'paper', yref: 'paper', x0: xLeft, x1: xRight, y0: metTop, y1: metTop,
                  line: { color: '#333', width: 1 } },
                { type: 'line', xref: 'paper', yref: 'paper', x0: xLeft, x1: xLeft, y0: metBot, y1: metTop,
                  line: { color: '#333', width: 1 } },
                { type: 'line', xref: 'paper', yref: 'paper', x0: xRight, x1: xRight, y0: metBot, y1: metTop,
                  line: { color: '#333', width: 1 } }
            );
        }

        const metricLabel = document.getElementById('metricColumn').value || 'Ranking metric';

        // Annotations
        const annotations = [
            // Title
            {
                text: `<b>${geneSetName.replace(/_/g, ' ')}</b>`,
                xref: 'paper', yref: 'paper', x: 0.5, y: 1.06,
                showarrow: false,
                font: { size: baseFontSize + 1, family: fontFam },
                xanchor: 'center'
            }
        ];

        // Stats annotation — position in the quadrant least occupied by the ES curve
        if (s.showStatsBox) {
            const isPositive = result.nes >= 0;
            // Sample the ES curve at 4 quadrants to find the emptiest corner
            const half = Math.floor(N / 2);
            // Average ES in left half vs right half
            let leftSum = 0, rightSum = 0;
            const sampleStep = Math.max(1, Math.floor(N / 200));
            let leftCount = 0, rightCount = 0;
            for (let i = 0; i < N; i += sampleStep) {
                if (i < half) { leftSum += Math.abs(result.runningES[i]); leftCount++; }
                else { rightSum += Math.abs(result.runningES[i]); rightCount++; }
            }
            const leftAvg = leftCount > 0 ? leftSum / leftCount : 0;
            const rightAvg = rightCount > 0 ? rightSum / rightCount : 0;

            // Find peak for vertical placement
            let peakVal = result.runningES[0];
            for (let i = 1; i < result.runningES.length; i++) {
                if (Math.abs(result.runningES[i]) > Math.abs(peakVal)) {
                    peakVal = result.runningES[i];
                }
            }

            // Place horizontally opposite to where the curve is densest
            const curveIsLeft = leftAvg > rightAvg;
            const boxX = curveIsLeft ? 0.97 : 0.03;
            const boxXanchor = curveIsLeft ? 'right' : 'left';
            // Place vertically opposite to the peak direction
            let boxY, boxYanchor;
            if (peakVal >= 0) {
                // Peak is positive (top), so place box at bottom of ES panel
                boxY = esBot + (esTop - esBot) * 0.04;
                boxYanchor = 'bottom';
            } else {
                // Peak is negative (bottom), so place box at top of ES panel
                boxY = esTop - (esTop - esBot) * 0.04;
                boxYanchor = 'top';
            }

            annotations.push({
                text: `NES = ${result.nes.toFixed(2)}<br>FDR = ${this.formatPval(result.fdr)}<br>p = ${this.formatPval(result.pvalue)}<br>Size = ${result.size}`,
                xref: 'paper', yref: 'paper',
                x: boxX, y: boxY,
                showarrow: false,
                font: { size: Math.max(9, baseFontSize - 2), family: 'Roboto Mono, monospace', color: '#333' },
                align: boxXanchor === 'right' ? 'right' : 'left',
                xanchor: boxXanchor,
                yanchor: boxYanchor,
                bgcolor: 'rgba(255,255,255,0.9)',
                bordercolor: '#ccc',
                borderwidth: 1,
                borderpad: 5
            });
        }

        // Zero cross annotation
        if (zeroCross >= 0 && s.showZeroCross) {
            annotations.push({
                text: `Zero cross at ${zeroCross}`,
                x: zeroCross, y: 0,
                xref: 'x3', yref: 'y3',
                showarrow: true, arrowhead: 0, arrowsize: 0.8, arrowwidth: 1,
                ax: 0, ay: -22,
                font: { size: Math.max(8, baseFontSize - 4), color: '#666', family: fontFam }
            });
        }

        // Correlation labels — position below the x-axis title to avoid overlap
        if (s.showCorrelationLabels && zeroCross >= 0) {
            annotations.push(
                {
                    text: '<i>Positively correlated</i>',
                    xref: 'paper', yref: 'paper', x: 0.0, y: -0.07,
                    showarrow: false,
                    font: { size: Math.max(8, baseFontSize - 3), color: '#dc2626', family: fontFam },
                    xanchor: 'left', yanchor: 'top'
                },
                {
                    text: '<i>Negatively correlated</i>',
                    xref: 'paper', yref: 'paper', x: 1.0, y: -0.07,
                    showarrow: false,
                    font: { size: Math.max(8, baseFontSize - 3), color: '#2563eb', family: fontFam },
                    xanchor: 'right', yanchor: 'top'
                }
            );
        }

        const tickFontSize = Math.max(8, baseFontSize - 3);

        const layout = {
            // ES panel
            xaxis: {
                range: [-xPad, N - 1 + xPad], showticklabels: false, showgrid: false, zeroline: false,
                domain: [xLeft, xRight], anchor: 'y', showline: false
            },
            yaxis: {
                title: { text: 'Enrichment score (ES)', font: { size: baseFontSize - 1, family: fontFam }, standoff: 8 },
                domain: [esBot, esTop], anchor: 'x',
                gridcolor: '#eee', gridwidth: 1,
                zeroline: false,
                tickfont: { size: tickFontSize, family: fontFam }
            },
            // Hit marker panel
            xaxis2: {
                range: [-xPad, N - 1 + xPad], showticklabels: false, showgrid: false, zeroline: false,
                domain: [xLeft, xRight], anchor: 'y2', showline: false
            },
            yaxis2: {
                domain: [hitBot, hitTop], anchor: 'x2',
                range: [0, 1],
                showticklabels: false, showgrid: false, zeroline: false, fixedrange: true,
                title: ''
            },
            // Ranked metric panel
            xaxis3: {
                range: [-xPad, N - 1 + xPad],
                title: { text: 'Rank in Ordered Dataset', font: { size: baseFontSize - 1, family: fontFam }, standoff: 4 },
                domain: [xLeft, xRight], anchor: 'y3',
                showgrid: false,
                tickfont: { size: tickFontSize, family: fontFam },
                side: 'bottom'
            },
            yaxis3: {
                title: { text: metricLabel, font: { size: baseFontSize - 2, family: fontFam }, standoff: 5 },
                domain: [metBot, metTop], anchor: 'x3',
                gridcolor: '#eee', gridwidth: 1,
                zeroline: true, zerolinecolor: '#333', zerolinewidth: 0.8,
                tickfont: { size: tickFontSize, family: fontFam },
                range: [Math.min(0, ...metSampled) * 1.08, Math.max(0, ...metSampled) * 1.08]
            },
            height: 580,
            margin: { l: 65, r: 15, t: 45, b: 60 },
            font: { family: fontFam },
            paper_bgcolor: s.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            shapes: panelShapes,
            annotations: annotations,
            bargap: 0
        };

        // Apply custom dimensions
        const dims = this.getPlotDimensions(undefined, layout.height);
        if (dims.width) layout.width = dims.width;
        if (dims.height) layout.height = dims.height;

        Plotly.newPlot('esPlot',
            [esLine, esZero, rugAnchor, ...traces3, metZero],
            layout,
            {
                responsive: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d'],
                edits: { annotationPosition: true, annotationTail: true, annotationText: true, axisTitleText: false, titleText: false, legendPosition: true, colorbarPosition: true }
            }
        );
    }

    // --------------------------------------------------------
    // Gene Set Navigation
    // --------------------------------------------------------
    navigateGeneSet(direction) {
        const sel = document.getElementById('geneSetSelector');
        const options = Array.from(sel.options).filter(o => o.value);
        const currentIdx = options.findIndex(o => o.value === sel.value);
        const newIdx = Math.max(0, Math.min(options.length - 1, currentIdx + direction));
        if (options[newIdx]) {
            sel.value = options[newIdx].value;
            this.renderESPlot(sel.value);
            this.renderGeneSetInfo(sel.value);
            this.renderGeneDetailTable(sel.value);
        }
    }

    // --------------------------------------------------------
    // Gene Set Info Panel
    // --------------------------------------------------------
    renderGeneSetInfo(geneSetName) {
        const result = this.results.find(r => r.name === geneSetName);
        if (!result) return;
        const el = document.getElementById('geneSetInfoContent');

        // Format the name nicely
        const displayName = this.cleanName(geneSetName);
        const direction = result.nes > 0 ? 'Upregulated' : 'Downregulated';
        const dirColor = result.nes > 0 ? '#dc2626' : '#2563eb';
        const sigLabel = result.fdr < 0.05 ? '★ Significant' : result.fdr < 0.25 ? '◆ Suggestive' : '○ Not significant';
        const sigColor = result.fdr < 0.05 ? '#15803d' : result.fdr < 0.25 ? '#d97706' : '#6b7280';

        el.innerHTML = `
            <div style="font-weight: 600; margin-bottom: 6px; line-height: 1.3; word-break: break-word;">${displayName}</div>
            <div style="display: grid; grid-template-columns: auto 1fr; gap: 2px 10px; font-size: 0.95em;">
                <span style="color: var(--gray-500);">NES:</span>
                <span style="font-weight: 600; color: ${dirColor};">${result.nes.toFixed(3)}</span>
                <span style="color: var(--gray-500);">FDR:</span>
                <span style="font-family: 'Roboto Mono', monospace;">${this.formatPval(result.fdr)}</span>
                <span style="color: var(--gray-500);">p-value:</span>
                <span style="font-family: 'Roboto Mono', monospace;">${this.formatPval(result.pvalue)}</span>
                <span style="color: var(--gray-500);">ES:</span>
                <span>${result.es.toFixed(4)}</span>
                <span style="color: var(--gray-500);">Size:</span>
                <span>${result.size} genes</span>
                <span style="color: var(--gray-500);">Leading Edge:</span>
                <span>${(result.leadingEdge || []).length} genes</span>
                <span style="color: var(--gray-500);">Direction:</span>
                <span style="color: ${dirColor}; font-weight: 500;">${direction}</span>
                <span style="color: var(--gray-500);">Status:</span>
                <span style="color: ${sigColor}; font-weight: 500;">${sigLabel}</span>
            </div>
        `;
    }

    // --------------------------------------------------------
    // Gene Detail Table (hits in gene set with metrics)
    // --------------------------------------------------------
    renderGeneDetailTable(geneSetName) {
        const result = this.results.find(r => r.name === geneSetName);
        if (!result || !this.rankedList) return;

        const tbody = document.getElementById('geneDetailBody');
        tbody.innerHTML = '';

        // Update column headers with context
        const metricColName = document.getElementById('metricColumn').value || 'Metric';
        const metricHeader = document.getElementById('geneDetailMetricHeader');
        if (metricHeader) {
            metricHeader.innerHTML = `${metricColName} <span class="detail-col-context">(your data)</span> <span class="sort-arrow">&#9650;</span>`;
        }
        const rankHeader = document.getElementById('geneDetailRankHeader');
        if (rankHeader) {
            rankHeader.innerHTML = `Rank <span class="detail-col-context">(gene set)</span> <span class="sort-arrow">&#9650;</span>`;
        }

        const leadingEdgeSet = new Set((result.leadingEdge || []).map(g => g.toUpperCase()));
        const genes = this.rankedList.genes;
        const metrics = this.rankedList.metrics;

        // Build rows for hits
        const rows = [];
        for (const hitIdx of result.hits) {
            const gene = genes[hitIdx];
            const metric = metrics[hitIdx];
            const es = result.runningES[hitIdx];
            const isLE = leadingEdgeSet.has(gene.toUpperCase());
            rows.push({ rank: hitIdx, gene, metric, es, isLE });
        }

        // Sort
        const col = this.detailSortCol;
        const asc = this.detailSortAsc;
        rows.sort((a, b) => {
            let va = a[col], vb = b[col];
            if (col === 'gene') { va = va.toLowerCase(); vb = vb.toLowerCase(); }
            if (col === 'isLE') { va = va ? 1 : 0; vb = vb ? 1 : 0; }
            if (va < vb) return asc ? -1 : 1;
            if (va > vb) return asc ? 1 : -1;
            return 0;
        });

        for (const row of rows) {
            const tr = document.createElement('tr');
            const leStyle = row.isLE ? 'font-weight:600; color: var(--green-700);' : 'color: var(--gray-500);';
            tr.innerHTML = `
                <td style="text-align: center;">${row.rank + 1}</td>
                <td class="gene-hover" data-gene="${row.gene}" style="font-family: 'Roboto Mono', monospace; font-size: 0.95em; cursor: help; text-align: center;">${row.gene}</td>
                <td style="font-family: 'Roboto Mono', monospace; text-align: center;">${row.metric.toFixed(4)}</td>
                <td style="font-family: 'Roboto Mono', monospace; text-align: center;">${row.es.toFixed(4)}</td>
                <td style="${leStyle} text-align: center;">${row.isLE ? 'Yes' : 'No'}</td>
            `;
            tbody.appendChild(tr);
        }

        // Update sort indicators
        document.querySelectorAll('#geneDetailTable thead th').forEach(th => {
            const arrow = th.querySelector('.sort-arrow');
            if (!arrow) return;
            if (th.dataset.detailCol === this.detailSortCol) {
                th.classList.add('sorted');
                arrow.textContent = this.detailSortAsc ? '\u25B2' : '\u25BC';
            } else {
                th.classList.remove('sorted');
                arrow.textContent = '\u25B2';
            }
        });

        // Attach gene hover tooltips
        this.attachGeneTooltips(tbody);

        this._currentGeneDetailRows = rows;
    }

    downloadGeneDetailCSV() {
        const sel = document.getElementById('geneSetSelector').value;
        if (!sel || !this._currentGeneDetailRows) return;

        const headers = ['Rank', 'Gene', 'Ranking Metric', 'Running ES', 'Leading Edge'];
        const csvRows = this._currentGeneDetailRows.map(r =>
            [r.rank + 1, r.gene, r.metric.toFixed(6), r.es.toFixed(6), r.isLE ? 'Yes' : 'No']
        );
        const csv = [headers.join(','), ...csvRows.map(r => r.join(','))].join('\n');
        const setName = sel.replace(/[^a-zA-Z0-9_]/g, '_').substring(0, 60);
        this.downloadBlob(new Blob([csv], { type: 'text/csv' }), `${setName}_genes.csv`);
    }

    // --------------------------------------------------------
    // Gene Search — highlight genes in the ES plot
    // --------------------------------------------------------
    searchGenes() {
        const input = document.getElementById('geneSearchInput').value.trim();
        if (!input || !this.rankedList) return;

        // Parse gene list (comma, newline, space, tab separated)
        const searchGenes = input.split(/[\n,\t\s]+/)
            .map(g => g.trim().toUpperCase())
            .filter(g => g.length > 0);

        const genes = this.rankedList.genes;
        const metrics = this.rankedList.metrics;
        const resultsEl = document.getElementById('geneSearchResults');

        // Get current gene set to check membership
        const currentSetName = document.getElementById('geneSetSelector').value;
        const currentResult = currentSetName ? this.results.find(r => r.name === currentSetName) : null;
        const hitsSet = currentResult ? new Set(currentResult.hits.map(i => genes[i])) : new Set();

        let html = '';
        const foundInSet = [];
        const foundInData = [];
        const notFound = [];

        for (const sg of searchGenes) {
            const idx = genes.indexOf(sg);
            if (idx >= 0) {
                const inGeneSet = hitsSet.has(sg);
                const entry = { gene: sg, rank: idx, metric: metrics[idx], inGeneSet };
                if (inGeneSet) {
                    foundInSet.push(entry);
                } else {
                    foundInData.push(entry);
                }
            } else {
                notFound.push(sg);
            }
        }

        const allFound = [...foundInSet, ...foundInData].sort((a, b) => a.rank - b.rank);

        if (foundInSet.length > 0) {
            html += `<div style="color: var(--green-700); font-weight: 600; margin-bottom: 3px;">In gene set (${foundInSet.length}):</div>`;
            html += '<table style="width:100%; border-collapse: collapse;">';
            for (const f of foundInSet.sort((a, b) => a.rank - b.rank)) {
                const color = f.metric >= 0 ? '#dc2626' : '#2563eb';
                html += `<tr style="border-bottom: 1px solid #eee;">
                    <td style="padding: 1px 4px; font-family: Roboto Mono, monospace;">${f.gene}</td>
                    <td style="padding: 1px 4px; text-align: right;">Rank ${f.rank + 1}</td>
                    <td style="padding: 1px 4px; text-align: right; color: ${color};">${f.metric.toFixed(3)}</td>
                </tr>`;
            }
            html += '</table>';
        }
        if (foundInData.length > 0) {
            html += `<div style="color: #6b7280; font-weight: 600; margin-top: 4px; margin-bottom: 3px;">In data only (${foundInData.length}):</div>`;
            html += '<table style="width:100%; border-collapse: collapse;">';
            for (const f of foundInData.sort((a, b) => a.rank - b.rank)) {
                const color = f.metric >= 0 ? '#dc2626' : '#2563eb';
                html += `<tr style="border-bottom: 1px solid #eee;">
                    <td style="padding: 1px 4px; font-family: Roboto Mono, monospace;">${f.gene}</td>
                    <td style="padding: 1px 4px; text-align: right;">Rank ${f.rank + 1}</td>
                    <td style="padding: 1px 4px; text-align: right; color: ${color};">${f.metric.toFixed(3)}</td>
                </tr>`;
            }
            html += '</table>';
        }
        if (notFound.length > 0) {
            html += `<div style="color: #dc2626; margin-top: 4px; font-size: 0.95em;">Not found: ${notFound.join(', ')}</div>`;
        }

        resultsEl.innerHTML = html;

        // Highlight on the ES plot if it's showing
        this._highlightGenesOnPlot(allFound.map(f => ({ gene: f.gene, rank: f.rank, inGeneSet: f.inGeneSet })));
    }

    _highlightGenesOnPlot(geneData) {
        // geneData: array of { gene, rank, inGeneSet }
        const plotEl = document.getElementById('esPlot');
        if (!plotEl || !plotEl.data) return;

        const currentLayout = plotEl.layout || {};
        const currentShapes = (currentLayout.shapes || []).filter(s => !s._isHighlight);
        const currentAnnotations = (currentLayout.annotations || []).filter(a => !a._isHighlight);

        for (const g of geneData) {
            const lineColor = g.inGeneSet ? 'rgba(255, 140, 0, 0.7)' : 'rgba(100, 100, 255, 0.5)';
            const lineStyle = g.inGeneSet ? 'dot' : 'dashdot';
            const labelColor = g.inGeneSet ? '#e65100' : '#4444cc';

            currentShapes.push({
                type: 'line',
                x0: g.rank, x1: g.rank,
                y0: 0, y1: 1,
                xref: 'x', yref: 'paper',
                line: { color: lineColor, width: 2, dash: lineStyle },
                _isHighlight: true
            });
            // Gene name label BELOW the graph (in the margin area)
            const suffix = g.inGeneSet === false ? ' *' : '';
            currentAnnotations.push({
                text: `<b>${g.gene}${suffix}</b>`,
                x: g.rank,
                y: -0.02,
                xref: 'x3', yref: 'paper',
                showarrow: true,
                arrowhead: 0,
                arrowwidth: 1,
                arrowcolor: lineColor,
                ax: 0, ay: 18,
                font: { size: 9, color: labelColor, family: this.settings.fontFamily + ', sans-serif' },
                bgcolor: 'rgba(255,255,255,0.85)',
                borderpad: 2,
                xanchor: 'center',
                _isHighlight: true
            });
        }

        Plotly.relayout('esPlot', { shapes: currentShapes, annotations: currentAnnotations });
    }

    clearGeneSearch() {
        document.getElementById('geneSearchInput').value = '';
        document.getElementById('geneSearchResults').innerHTML = '';
        // Remove highlight shapes and annotations
        const plotEl = document.getElementById('esPlot');
        if (plotEl && plotEl.layout) {
            const shapes = (plotEl.layout.shapes || []).filter(s => !s._isHighlight);
            const annotations = (plotEl.layout.annotations || []).filter(a => !a._isHighlight);
            Plotly.relayout('esPlot', { shapes, annotations });
        }
    }

    // --------------------------------------------------------
    // Top & Bottom Genes
    // --------------------------------------------------------
    renderTopBottomGenes() {
        if (!this.rankedList) return;
        const n = parseInt(document.getElementById('topBottomN').value) || 10;
        const genes = this.rankedList.genes;
        const metrics = this.rankedList.metrics;
        const metricLabel = document.getElementById('metricColumn').value || 'metric';
        const el = document.getElementById('topBottomGenes');

        let html = '<div style="font-weight: 600; color: #dc2626; margin-bottom: 3px;">Top ' + n + ' in your dataset <span style="font-weight:400;color:#888;font-size:0.9em;">(highest ' + metricLabel + ')</span></div>';
        html += '<table style="width: 100%; border-collapse: collapse; margin-bottom: 8px;">';
        for (let i = 0; i < Math.min(n, genes.length); i++) {
            html += `<tr style="border-bottom: 1px solid #f0f0f0;">
                <td class="gene-hover" data-gene="${genes[i]}" style="padding: 1px 2px; font-family: Roboto Mono, monospace; font-size: 0.95em; cursor: help;">${genes[i]}</td>
                <td style="padding: 1px 2px; text-align: right; color: #dc2626;">${metrics[i].toFixed(3)}</td>
            </tr>`;
        }
        html += '</table>';

        html += '<div style="font-weight: 600; color: #2563eb; margin-bottom: 3px;">Bottom ' + n + ' in your dataset <span style="font-weight:400;color:#888;font-size:0.9em;">(lowest ' + metricLabel + ')</span></div>';
        html += '<table style="width: 100%; border-collapse: collapse;">';
        for (let i = genes.length - 1; i >= Math.max(0, genes.length - n); i--) {
            html += `<tr style="border-bottom: 1px solid #f0f0f0;">
                <td class="gene-hover" data-gene="${genes[i]}" style="padding: 1px 2px; font-family: Roboto Mono, monospace; font-size: 0.95em; cursor: help;">${genes[i]}</td>
                <td style="padding: 1px 2px; text-align: right; color: #2563eb;">${metrics[i].toFixed(3)}</td>
            </tr>`;
        }
        html += '</table>';

        el.innerHTML = html;
        // Attach gene info tooltips
        this.attachGeneTooltips(el);
    }

    // --------------------------------------------------------
    // Overlap Heatmap
    // --------------------------------------------------------
    renderOverlapHeatmap() {
        if (!this.results) return;
        this.readSettings();

        const maxSets = parseInt(document.getElementById('overlapMaxSets').value) || 20;
        const collapseThresh = parseInt(document.getElementById('overlapCollapseThresh').value) || 50;
        const condensed = document.getElementById('overlapCondensed').checked;
        const fontFam = this.settings.fontFamily + ', sans-serif';

        // Get significant gene sets sorted by |NES|
        let sigSets = this.results
            .filter(r => r.fdr < 0.25)
            .sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes));

        if (sigSets.length === 0) {
            Plotly.newPlot('overlapHeatmap', [], {
                annotations: [{
                    text: 'No significant gene sets (FDR < 0.25) to show',
                    xref: 'paper', yref: 'paper', x: 0.5, y: 0.5,
                    showarrow: false, font: { size: 14, color: '#666' }
                }],
                height: 200
            });
            return;
        }

        // Get the active gene sets for gene lists
        const activeGeneSets = this.getActiveGeneSets();

        // Build gene lists for significant results
        const setGenes = {};
        for (const r of sigSets) {
            if (activeGeneSets[r.name]) {
                const rankedGenesUpper = new Set(this.rankedList.genes.map(g => g.toUpperCase()));
                setGenes[r.name] = activeGeneSets[r.name]
                    .map(g => g.toUpperCase())
                    .filter(g => rankedGenesUpper.has(g));
            }
        }

        // Condensed view: collapse highly overlapping sets
        let collapsedInfo = {};
        if (condensed) {
            const kept = [];
            const collapsed = {};

            for (const r of sigSets) {
                if (!setGenes[r.name]) continue;
                let shouldCollapse = false;
                for (const keptSet of kept) {
                    const overlap = this._computeOverlap(setGenes[r.name], setGenes[keptSet.name]);
                    const smaller = Math.min(setGenes[r.name].length, setGenes[keptSet.name].length);
                    if (smaller > 0 && (overlap / smaller) * 100 >= collapseThresh) {
                        if (!collapsed[keptSet.name]) collapsed[keptSet.name] = [];
                        collapsed[keptSet.name].push(r.name);
                        shouldCollapse = true;
                        break;
                    }
                }
                if (!shouldCollapse) {
                    kept.push(r);
                }
            }
            sigSets = kept;
            collapsedInfo = collapsed;
        }

        // Cap at maxSets
        sigSets = sigSets.slice(0, maxSets);

        if (sigSets.length < 2) {
            Plotly.newPlot('overlapHeatmap', [], {
                annotations: [{
                    text: 'Need at least 2 gene sets for overlap analysis',
                    xref: 'paper', yref: 'paper', x: 0.5, y: 0.5,
                    showarrow: false, font: { size: 14, color: '#666' }
                }],
                height: 200
            });
            return;
        }

        // Compute overlap matrix (Jaccard similarity)
        const n = sigSets.length;
        const labels = sigSets.map(r => this.cleanName(r.name));
        const matrix = [];
        const textMatrix = [];

        for (let i = 0; i < n; i++) {
            const row = [];
            const textRow = [];
            const genesI = setGenes[sigSets[i].name] || [];
            for (let j = 0; j < n; j++) {
                const genesJ = setGenes[sigSets[j].name] || [];
                if (i === j) {
                    row.push(1);
                    textRow.push(`${genesI.length} genes`);
                } else {
                    const overlap = this._computeOverlap(genesI, genesJ);
                    const union = new Set([...genesI, ...genesJ]).size;
                    const jaccard = union > 0 ? overlap / union : 0;
                    row.push(jaccard);
                    textRow.push(`${overlap} shared<br>J=${jaccard.toFixed(2)}`);
                }
            }
            matrix.push(row);
            textMatrix.push(textRow);
        }

        const trace = {
            z: matrix,
            x: labels,
            y: labels,
            type: 'heatmap',
            colorscale: [
                [0, '#ffffff'],
                [0.2, '#e0f2e9'],
                [0.4, '#a8d89a'],
                [0.6, '#5a9f4a'],
                [0.8, '#3a7333'],
                [1.0, '#2d5a27']
            ],
            text: textMatrix,
            hovertemplate: '%{y} vs %{x}<br>%{text}<extra></extra>',
            showscale: true,
            colorbar: {
                title: { text: 'Jaccard Index', font: { size: 11, family: fontFam } },
                thickness: 12, len: 0.5
            }
        };

        const size = Math.max(450, n * 28 + 150);
        const layout = {
            height: size,
            width: size + 50,
            margin: { l: 10, r: 50, t: 20, b: 10 },
            xaxis: { tickfont: { size: 10, family: fontFam }, tickangle: -45, automargin: true, showgrid: false },
            yaxis: { tickfont: { size: 10, family: fontFam }, automargin: true, showgrid: false, autorange: 'reversed' },
            font: { family: fontFam },
            paper_bgcolor: '#fff',
            plot_bgcolor: '#fff'
        };

        Plotly.newPlot('overlapHeatmap', [trace], layout, {
            responsive: true, displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d']
        });

        // Show collapsed info
        const condensedEl = document.getElementById('condensedInfo');
        const condensedList = document.getElementById('condensedList');
        if (condensed && Object.keys(collapsedInfo).length > 0) {
            condensedEl.style.display = '';
            let html = '';
            let idx = 1;
            for (const [parent, children] of Object.entries(collapsedInfo)) {
                html += `<div style="margin-bottom: 6px;">
                    <span style="font-weight: 600; color: var(--green-700);">${idx}. ${this.cleanName(parent)}</span>
                    <span style="color: var(--gray-500);"> also includes:</span>
                    <div style="padding-left: 12px; color: var(--gray-600);">
                        ${children.map(c => this.cleanName(c)).join('<br>')}
                    </div>
                </div>`;
                idx++;
            }
            condensedList.innerHTML = html;
        } else {
            condensedEl.style.display = 'none';
        }
    }

    _computeOverlap(genesA, genesB) {
        const setB = new Set(genesB);
        let count = 0;
        for (const g of genesA) {
            if (setB.has(g)) count++;
        }
        return count;
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
        const nGenes = this.rankedList ? this.rankedList.genes.length : 0;
        const date = this.analysisDate
            ? this.analysisDate.toISOString().split('T')[0]
            : new Date().toISOString().split('T')[0];

        // Get the file name if available
        const fileInput = document.getElementById('fileInput');
        const fileName = fileInput.files && fileInput.files[0]
            ? fileInput.files[0].name
            : 'example data';

        const text = `Preranked Gene Set Enrichment Analysis (GSEA) was performed using a ` +
            `client-side JavaScript implementation of the GSEA algorithm ` +
            `(Subramanian et al., PNAS 2005; Mootha et al., Nature Genetics 2003). ` +
            `Input data (${fileName}, ${nGenes} genes) was ranked by ${metricCol} in descending order. ` +
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
            `(https://fredrikwermeling.github.io/GSEA_app/), ` +
            `built with Plotly.js v2.27.0, PapaParse v5.4.1, and Web Workers ` +
            `(Wermeling Lab, Karolinska Institutet).`;

        this.methodsText = text;
        document.getElementById('methodsPopupText').textContent = text;
    }

    copyMethods() {
        const text = this.methodsText || document.getElementById('methodsPopupText').textContent;
        navigator.clipboard.writeText(text).then(() => {
            const btn = document.getElementById('copyMethodsBtn');
            const orig = btn.textContent;
            btn.textContent = 'Copied!';
            setTimeout(() => { btn.textContent = orig; }, 2000);
        });
    }

    downloadMethods() {
        const text = this.methodsText || document.getElementById('methodsPopupText').textContent;
        const blob = new Blob([text], { type: 'text/plain' });
        this.downloadBlob(blob, 'gsea_methods.txt');
    }

    // --------------------------------------------------------
    // Export
    // --------------------------------------------------------
    exportPlot(plotId, format) {
        const plotEl = document.getElementById(plotId);
        if (!plotEl || !plotEl.data) {
            console.warn('No plot data found for', plotId);
            return;
        }

        this.readSettings();
        const scale = this.settings.exportScale;

        const opts = {
            format: format,
            width: plotEl.offsetWidth || 800,
            height: plotEl.offsetHeight || 500
        };

        if (format === 'png') {
            opts.scale = scale;
        }

        Plotly.toImage(plotEl, opts)
            .then(dataUrl => {
                if (format === 'svg') {
                    // SVG data URL: data:image/svg+xml,... (URL-encoded) or data:image/svg+xml;base64,...
                    let svgContent;
                    if (dataUrl.includes(';base64,')) {
                        svgContent = atob(dataUrl.split(';base64,')[1]);
                    } else {
                        svgContent = decodeURIComponent(dataUrl.split(',').slice(1).join(','));
                    }
                    const blob = new Blob([svgContent], { type: 'image/svg+xml' });
                    this.downloadBlob(blob, `${plotId}.svg`);
                } else {
                    // PNG: data URL approach
                    this.downloadBlob(this.dataURLtoBlob(dataUrl), `${plotId}.png`);
                }
            })
            .catch(err => {
                console.error('Export failed:', err);
            });
    }

    dataURLtoBlob(dataUrl) {
        const parts = dataUrl.split(';base64,');
        const contentType = parts[0].split(':')[1];
        const raw = atob(parts[1]);
        const arr = new Uint8Array(raw.length);
        for (let i = 0; i < raw.length; i++) {
            arr[i] = raw.charCodeAt(i);
        }
        return new Blob([arr], { type: contentType });
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
    // Reset App
    // --------------------------------------------------------
    resetApp() {
        // Terminate worker
        if (this.worker) this.worker.terminate();

        // Reset state
        this.rawData = null;
        this.rankedList = null;
        this.customGeneSets = {};
        this.results = null;
        this.analysisDate = null;
        this._currentGeneDetailRows = null;

        // Reset UI
        document.getElementById('fileInput').value = '';
        const geneCol = document.getElementById('geneColumn');
        geneCol.innerHTML = '<option value="">Upload a file first...</option>';
        geneCol.disabled = true;
        const metricCol = document.getElementById('metricColumn');
        metricCol.innerHTML = '<option value="">Upload a file first...</option>';
        metricCol.disabled = true;
        document.getElementById('runBtn').disabled = true;
        document.getElementById('runBtn').style.display = '';
        document.getElementById('cancelBtn').style.display = 'none';
        document.getElementById('progressContainer').classList.remove('active');
        this.hideStatus('uploadStatus');
        this.hideStatus('runStatus');

        // Hide results
        document.getElementById('overviewEmpty').style.display = '';
        document.getElementById('overviewResults').style.display = 'none';
        document.getElementById('enrichmentEmpty').style.display = '';
        document.getElementById('enrichmentResults').style.display = 'none';
        document.getElementById('tableEmpty').style.display = '';
        document.getElementById('tableResults').style.display = 'none';
        document.getElementById('overlapEmpty').style.display = '';
        document.getElementById('overlapResults').style.display = 'none';
        document.getElementById('openSettingsBtn').style.display = 'none';
        document.getElementById('settingsPanel').classList.remove('open');
        document.getElementById('methodsCard').style.display = 'none';

        // Reset checkboxes
        document.getElementById('checkHallmark').checked = true;
        document.getElementById('checkC2').checked = false;
        document.getElementById('checkC5').checked = false;

        // Clear plots
        ['bubblePlot', 'rankedPlot', 'esPlot', 'overlapHeatmap'].forEach(id => {
            const el = document.getElementById(id);
            if (el) Plotly.purge(el);
        });

        // Show default tab
        this.showTab('overview');

        // Re-create worker
        this.createWorker();
        this.updateGeneSetStatus();
    }

    // --------------------------------------------------------
    // Gene Info Tooltips (MyGene.info)
    // --------------------------------------------------------
    async fetchGeneInfo(gene) {
        if (this.geneInfoCache[gene]) return this.geneInfoCache[gene];
        try {
            const res = await fetch(`https://mygene.info/v3/query?q=symbol:${gene}&species=human&fields=symbol,name,summary&size=1`);
            const data = await res.json();
            if (data.hits && data.hits.length > 0) {
                const hit = data.hits[0];
                const info = {
                    symbol: hit.symbol || gene,
                    name: hit.name || '',
                    summary: hit.summary || ''
                };
                this.geneInfoCache[gene] = info;
                return info;
            }
        } catch (e) {
            console.warn('Gene info fetch failed:', e);
        }
        return null;
    }

    showGeneTooltip(event, gene) {
        this.hideGeneTooltip();
        const tooltip = document.createElement('div');
        tooltip.id = 'geneTooltip';
        tooltip.style.cssText = 'position:fixed; z-index:10001; background:white; border:1px solid #d1d5db; border-radius:8px; padding:10px 14px; max-width:350px; box-shadow:0 4px 12px rgba(0,0,0,0.15); font-size:12px; line-height:1.5; pointer-events:none;';
        tooltip.innerHTML = `<div style="color:#6b7280;">Loading ${gene} info...</div>`;

        const x = Math.min(event.clientX + 12, window.innerWidth - 370);
        const y = Math.min(event.clientY + 12, window.innerHeight - 200);
        tooltip.style.left = x + 'px';
        tooltip.style.top = y + 'px';
        document.body.appendChild(tooltip);

        this.fetchGeneInfo(gene).then(info => {
            const el = document.getElementById('geneTooltip');
            if (!el) return;
            if (!info) {
                el.innerHTML = `<b>${gene}</b><br><span style="color:#999;">No info available</span>`;
                return;
            }
            let html = `<div style="margin-bottom:4px;"><b style="color:#5a9f4a; font-size:13px;">${info.symbol}</b> <span style="color:#374151;">${info.name}</span></div>`;
            if (info.summary) {
                const short = info.summary.length > 200 ? info.summary.substring(0, 200) + '...' : info.summary;
                html += `<div style="color:#4b5563; font-size:11px;">${short}</div>`;
            }
            el.innerHTML = html;
            // Reposition if needed
            const rect = el.getBoundingClientRect();
            if (rect.bottom > window.innerHeight) {
                el.style.top = (window.innerHeight - rect.height - 10) + 'px';
            }
        });
    }

    hideGeneTooltip() {
        const existing = document.getElementById('geneTooltip');
        if (existing) existing.remove();
    }

    attachGeneTooltips(container) {
        container.querySelectorAll('.gene-hover').forEach(el => {
            el.addEventListener('mouseenter', (e) => {
                this._tooltipTimer = setTimeout(() => {
                    this.showGeneTooltip(e, el.dataset.gene);
                }, 400);
            });
            el.addEventListener('mouseleave', () => {
                clearTimeout(this._tooltipTimer);
                this.hideGeneTooltip();
            });
        });
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
        // Figure customization
        const fontEl = document.getElementById('figFontFamily');
        if (fontEl) this.settings.fontFamily = fontEl.value;
        const fontSizeEl = document.getElementById('figFontSize');
        if (fontSizeEl) this.settings.fontSize = parseInt(fontSizeEl.value) || 12;
        const esColorEl = document.getElementById('esLineColor');
        if (esColorEl) this.settings.esLineColor = esColorEl.value;
        const esWidthEl = document.getElementById('esLineWidth');
        if (esWidthEl) this.settings.esLineWidth = parseFloat(esWidthEl.value) || 2.5;
        // Toggle elements
        const cb = (id) => { const el = document.getElementById(id); return el ? el.checked : true; };
        this.settings.showStatsBox = cb('showStatsBox');
        this.settings.showZeroCross = cb('showZeroCross');
        this.settings.showCorrelationLabels = cb('showCorrelationLabels');
        this.settings.showPanelBorders = cb('showPanelBorders');
        this.settings.showHitMarkers = cb('showHitMarkers');
        this.settings.showMetricFill = cb('showMetricFill');
    }

    updateSettings() {
        this.readSettings();
        if (this.results) {
            // Re-render the plot for the active tab
            if (this.activeTab === 'overview') {
                this.renderBubblePlot();
                this.renderRankedPlot();
            } else if (this.activeTab === 'enrichment') {
                const selected = document.getElementById('geneSetSelector').value;
                if (selected) this.renderESPlot(selected);
            } else if (this.activeTab === 'overlap') {
                this.renderOverlapHeatmap();
            }
        }
    }

    updateSettingsTabVisibility() {
        // Show/hide per-tab settings sections
        const overviewSection = document.getElementById('settingsOverview');
        const enrichmentSection = document.getElementById('settingsEnrichment');
        if (overviewSection) overviewSection.style.display = (this.activeTab === 'overview') ? '' : 'none';
        if (enrichmentSection) enrichmentSection.style.display = (this.activeTab === 'enrichment') ? '' : 'none';
    }

    initSettingsDrag() {
        const panel = document.getElementById('settingsPanel');
        const header = document.getElementById('settingsPanelHeader');
        let isDragging = false, startX, startY, startLeft, startTop;

        header.addEventListener('mousedown', (e) => {
            if (e.target.tagName === 'BUTTON') return;
            isDragging = true;
            const rect = panel.getBoundingClientRect();
            startX = e.clientX;
            startY = e.clientY;
            startLeft = rect.left;
            startTop = rect.top;
            e.preventDefault();
        });
        document.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            const dx = e.clientX - startX;
            const dy = e.clientY - startY;
            panel.style.left = (startLeft + dx) + 'px';
            panel.style.top = (startTop + dy) + 'px';
            panel.style.right = 'auto';
        });
        document.addEventListener('mouseup', () => { isDragging = false; });
    }

    getPlotDimensions(defaultW, defaultH) {
        const w = parseInt(document.getElementById('plotWidth').value) || 0;
        const h = parseInt(document.getElementById('plotHeight').value) || 0;
        return {
            width: w > 0 ? w : undefined,
            height: h > 0 ? h : defaultH
        };
    }

    // --------------------------------------------------------
    // Tab Navigation
    // --------------------------------------------------------
    showTab(tabName) {
        this.activeTab = tabName;
        document.querySelectorAll('.tab-btn').forEach(btn => {
            btn.classList.toggle('active', btn.dataset.tab === tabName);
        });
        document.querySelectorAll('.tab-panel').forEach(panel => {
            panel.classList.toggle('active', panel.id === 'tab-' + tabName);
        });
        // Update settings panel to show relevant controls
        this.updateSettingsTabVisibility();
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
